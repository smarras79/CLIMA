using MPI
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.DGBalanceLawDiscretizations
using CLIMA.DGBalanceLawDiscretizations.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.ODESolvers
using CLIMA.GenericCallbacks
using CLIMA.ParametersType
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates
using CLIMA.Vtk
using DelimitedFiles
using Dierckx
using Random

using TimerOutputs

const to = TimerOutput()

if haspkg("CuArrays")
    using CUDAdrv
    using CUDAnative
    using CuArrays
    CuArrays.allowscalar(false)
    const ArrayType = CuArray
else
    const ArrayType = Array
end

# Prognostic equations: ρ, (ρu), (ρv), (ρw), (ρe_tot), (ρq_tot)
# For the dry example shown here, we load the moist thermodynamics module 
# and consider the dry equation set to be the same as the moist equations but
# with total specific humidity = 0. 
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters

# State labels 
const _nstate = 6
const _ρ, _U, _V, _W, _E, _QT = 1:_nstate
const stateid = (ρid = _ρ, Uid = _U, Vid = _V, Wid = _W, Eid = _E, QTid = _QT)
const statenames = ("RHO", "U", "V", "W", "E", "QT")

# Viscous state labels
const _nviscstates = 23
const _τ11, _τ22, _τ33, _τ12, _τ13, _τ23, _qx, _qy, _qz, _Tx, _Ty, _Tz, _Trefx, _Trefy, _Trefz,  _Tpx, _Tpy, _Tpz, _θx, _θy, _θz, _SijSij, _ν_e = 1:_nviscstates

# Gradient state labels
# Gradient state labels
const _states_for_gradient_transform = (_ρ, _U, _V, _W, _E, _QT)

const _nauxstate = 16
const _a_x, _a_y, _a_z, _a_sponge, _a_02z, _a_z2inf, _a_rad, _a_ν_e, _a_LWP_02z, _a_LWP_z2inf,_a_q_liq,_a_θ, _a_P,_a_T, _a_Tref, _a_soundspeed_air = 1:_nauxstate

if !@isdefined integration_testing
    const integration_testing =
        parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end

# Problem constants (TODO: parameters module (?))
@parameter Prandtl 71 // 100 "Prandtl (molecular)"
@parameter  μ_sgs 100.0 "Constant dynamic viscosity"
@parameter Prandtl_t 1//3 "Prandtl_t"
@parameter cp_over_prandtl cp_d / Prandtl_t "cp_over_prandtl"

# Random number seed
const seed = MersenneTwister(0)


function global_max(A::MPIStateArray, states=1:size(A, 2))
    host_array = Array ∈ typeof(A).parameters
    h_A = host_array ? A : Array(A)
    locmax = maximum(view(h_A, :, states, A.realelems)) 
    MPI.Allreduce([locmax], MPI.MAX, A.mpicomm)[1]
end

function global_mean(A::MPIStateArray, states=1:size(A,2))
    host_array = Array ∈ typeof(A).parameters
    h_A = host_array ? A : Array(A) 
    (Np, nstate, nelem) = size(A) 
    numpts = (nelem * Np) + 1
    localsum = sum(view(h_A, :, states, A.realelems)) 
    MPI.Allreduce([localsum], MPI.SUM, A.mpicomm)[1] / numpts 
end

# User Input
const numdims = 2
const Npoly = 4

# Define grid size 
const Δx    = 125
const Δy    = 125
const Δz    = 125

const stretch_coe = 2.25

# Physical domain extents 
const (xmin, xmax) = (0,  5000)
const (ymin, ymax) = (0,  6000)
const (zmin, zmax) = (0,  6000)

#Get Nex, Ney from resolution
const Lx = xmax - xmin
const Ly = ymax - ymin
const Lz = zmax - ymin

#
# User defines the grid size:
#
Nex = ceil(Int64, Lx / (Δx * Npoly))
Ney = ceil(Int64, Ly / (Δy * Npoly))
Nez = ceil(Int64, Lz / (Δz * Npoly))


DoF = (Nex*Ney*Nez)*(Npoly+1)^numdims*(_nstate)
DoFstorage = (Nex*Ney*Nez)*(Npoly+1)^numdims*(_nstate + _nviscstates + _nauxstate + CLIMA.Mesh.Grids._nvgeo) +
    (Nex*Ney*Nez)*(Npoly+1)^(numdims-1)*2^numdims*(CLIMA.Mesh.Grids._nsgeo)


# Smagorinsky model requirements : TODO move to SubgridScaleTurbulence module 
@parameter C_smag 0.23 "C_smag"
# Equivalent grid-scale
#Δ = (Δx * Δy * Δz)^(1/3)
#Δ = max(Δx, Δy)
Δ = sqrt(Δx^2 + Δy^2)
const Δsqr = Δ * Δ

# Surface values to calculate surface fluxes:
const SST         = 292.5
const p_sfc       = 1017.8e2      # Pa
const q_tot_sfc   = 13.84e-3      # qs(sst) using Teten's formula
const ρ_sfc       = 1.22          #kg/m^3
const ft          =  15.0
const fq          = 115.0
const Cd          = 0.0011        #Drag coefficient
const first_node_level   = 0.5*Δy
# -------------------------------------------------------------------------
# Preflux calculation: This function computes parameters required for the 
# DG RHS (but not explicitly solved for as a prognostic variable)
# In the case of the rising_thermal_bubble example: the saturation
# adjusted temperature and pressure are such examples. Since we define
# the equation and its arguments here the user is afforded a lot of freedom
# around its behaviour. 
# The preflux function interacts with the following  
# Modules: NumericalFluxes.jl 
# functions: wavespeed, cns_flux!, bcstate!
# -------------------------------------------------------------------------
@inline function preflux(Q, aux)
    @inbounds begin
        ρ, U, V, W = Q[_ρ], Q[_U], Q[_V], Q[_W]
        ρinv = 1 / ρ
        u, v, w = ρinv * U, ρinv * V, ρinv * W
    end
end

#-------------------------------------------------------------------------
#md # Soundspeed computed using the thermodynamic state TS
# max eigenvalue
@inline function wavespeed(n, Q, aux, t)
    ρ, U, V, W = Q[_ρ], Q[_U], Q[_V], Q[_W]
    u, v, w = U/ρ, V/ρ, W/ρ
    @inbounds begin
        abs(n[1]*u + n[2]*v + n[3]*w) + aux[_a_soundspeed_air]
    end
end


# -------------------------------------------------------------------------
# ### read sounding
#md # 
#md # The sounding file contains the following quantities along a 1D column.
#md # It needs to have the following structure:
#md #
#md # z[m]   theta[K]  q[g/kg]   u[m/s]   v[m/s]   p[Pa]
#md # ...      ...       ...      ...      ...      ...
#md #
#md #
# -------------------------------------------------------------------------
function read_sounding()
    #read in the original squal sounding
    fsounding  = open(joinpath(@__DIR__, "../soundings/SOUNDING_PYCLES_Z_T_P.dat"))
    sounding = readdlm(fsounding)
    close(fsounding)
    (nzmax, ncols) = size(sounding)
    if nzmax == 0
        error("SOUNDING ERROR: The Sounding file is empty!")
    end
    return (sounding, nzmax, ncols)
end

# -------------------------------------------------------------------------
# ### Physical Flux (Required)
#md # Here, we define the physical flux function, i.e. the conservative form
#md # of the equations of motion for the prognostic variables ρ, U, V, W, E, QT
#md # $\frac{\partial Q}{\partial t} + \nabla \cdot \boldsymbol{F} = \boldsymbol {S}$
#md # $\boldsymbol{F}$ contains both the viscous and inviscid flux components
#md # and $\boldsymbol{S}$ contains source terms.
#md # Note that the preflux calculation is splatted at the end of the function call
#md # to cns_flux!
# -------------------------------------------------------------------------
@inline function cns_flux!(F, Q, VF, aux, t)
    @inbounds begin
        u, v, w = preflux(Q,aux)
        DFloat = eltype(F)
        D_subsidence = 3.75e-6
        
        ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
        
        P = aux[_a_P]
        
        xvert = aux[_a_y]
        
        # Inviscid contributions
        F[1, _ρ],  F[2, _ρ],  F[3, _ρ]  = U          , V          , W
        F[1, _U],  F[2, _U],  F[3, _U]  = u * U  + P , v * U      , w * U
        F[1, _V],  F[2, _V],  F[3, _V]  = u * V      , v * V + P  , w * V
        F[1, _W],  F[2, _W],  F[3, _W]  = u * W      , v * W      , w * W + P
        F[1, _E],  F[2, _E],  F[3, _E]  = u * (E + P), v * (E + P), w * (E + P)
        F[1, _QT], F[2, _QT], F[3, _QT] = u * QT     , v * QT     , w * QT

        #Derivative of T and Q:
        vqx, vqy, vqz = VF[_qx], VF[_qy], VF[_qz]
        vTx, vTy, vTz = VF[_Tx], VF[_Ty], VF[_Tz]
        vTrefx, vTrefy, vTrefz = VF[_Trefx], VF[_Trefy], VF[_Trefz]
        vTpx, vTpy, vTpz = VF[_Tpx], VF[_Tpy], VF[_Tpz]

        # Radiation contribution
        SijSij = VF[_SijSij]

        #Dynamic eddy viscosity
        #μ_e = ρ*VF[_ν_e] #Vreman
        μ_e = ρ*sqrt(2SijSij) * C_smag^2 * Δsqr  # Smagorinsky   
        D_e = μ_e / Prandtl_t

        # Multiply stress tensor by viscosity coefficient:
        τ11, τ22, τ33 = VF[_τ11] * μ_e, VF[_τ22]* μ_e, VF[_τ33] * μ_e
        τ12 = τ21 = VF[_τ12] * μ_e
        τ13 = τ31 = VF[_τ13] * μ_e
        τ23 = τ32 = VF[_τ23] * μ_e

        # Viscous velocity flux (i.e. F^visc_u in Giraldo Restelli 2008)
        F[1, _U] += τ11; F[2, _U] += τ12; F[3, _U] += τ13
        F[1, _V] += τ21; F[2, _V] += τ22; F[3, _V] += τ23
        F[1, _W] += τ31; F[2, _W] += τ32; F[3, _W] += τ33

        # Viscous Energy flux (i.e. F^visc_e in Giraldo Restelli 2008)
        F[1, _E] += u * τ11 + v * τ12 + w * τ13 + cp_over_prandtl * vTx * μ_e
        F[2, _E] += u * τ21 + v * τ22 + w * τ23 + cp_over_prandtl * vTy * μ_e
        F[3, _E] += u * τ31 + v * τ32 + w * τ33 + cp_over_prandtl * vTz * μ_e

        #if (xvert < 0.0001 || xvert > ymax - 0.0001)
        #    F[1, _E] = 0.0
        #    F[2, _E] = cp_v*vTrefy/Prandtl
        #end

        
        # Viscous contributions to mass flux terms
        #F[1, _ρ]  -=  vqx * D_e
        #F[2, _ρ]  -=  vqy * D_e
        #F[3, _ρ]  -=  vqz * D_e
        F[1, _QT] -=  vqx * D_e
        F[2, _QT] -=  vqy * D_e
        F[3, _QT] -=  vqz * D_e
        
    end
end

# -------------------------------------------------------------------------
#md # Here we define a function to extract the velocity components from the 
#md # prognostic equations (i.e. the momentum and density variables). This 
#md # function is not required in general, but provides useful functionality 
#md # in some cases. 
# -------------------------------------------------------------------------
# Compute the velocity from the state
const _ngradstates = 8
@inline function gradient_vars!(gradient_list, Q, aux, t)
    @inbounds begin
        u, v, w = preflux(Q,aux)
        T = aux[_a_T]
        θ = aux[_a_θ]
        xvert = aux[_a_y]
        
        ρ, QT =Q[_ρ], Q[_QT]
        # ordering should match states_for_gradient_transform
        gradient_list[1], gradient_list[2], gradient_list[3] = u, v, w
        gradient_list[4], gradient_list[5], gradient_list[6] = θ, QT/ρ, T

        #Get T_ref:
        θ_ref            = 300;
        π_exner_ref      = 1.0 - grav / (cp_d * θ_ref) * xvert # exner pressure
        ρ_ref            = MSLP / (R_d * θ_ref) * (π_exner_ref)^(cv_d / R_d) # density
        P_ref            = MSLP * (R_d * (ρ_ref * θ_ref) / MSLP)^(cp_d/cv_d) # pressure (absolute)
        T_ref            = P_ref / (ρ_ref * R_d) # T reference    
        gradient_list[7] = T_ref
        gradient_list[8] = T - T_ref #T perturbation
    end
end

# -------------------------------------------------------------------------
#md ### Viscous fluxes. 
#md # The viscous flux function compute_stresses computes the components of 
#md # the velocity gradient tensor, and the corresponding strain rates to
#md # populate the viscous flux array VF. SijSij is calculated in addition
#md # to facilitate implementation of the constant coefficient Smagorinsky model
#md # (pending)
@inline function compute_stresses!(VF, grad_mat, _...)
    @inbounds begin
        dudx, dudy, dudz = grad_mat[1, 1], grad_mat[2, 1], grad_mat[3, 1]
        dvdx, dvdy, dvdz = grad_mat[1, 2], grad_mat[2, 2], grad_mat[3, 2]
        dwdx, dwdy, dwdz = grad_mat[1, 3], grad_mat[2, 3], grad_mat[3, 3]
        # compute gradients of moist vars and temperature
        dθdx, dθdy, dθdz = grad_mat[1, 4], grad_mat[2, 4], grad_mat[3, 4]
        dqdx, dqdy, dqdz = grad_mat[1, 5], grad_mat[2, 5], grad_mat[3, 5]
        dTdx, dTdy, dTdz = grad_mat[1, 6], grad_mat[2, 6], grad_mat[3, 6]
        dTrefdx, dTrefdy, dTrefdz = grad_mat[1, 7], grad_mat[2, 7], grad_mat[3, 7]
        dTpdx,   dTpdy,   dTpdz   = grad_mat[1, 8], grad_mat[2, 8], grad_mat[3, 8]
        
        # virtual potential temperature gradient: for richardson calculation
        # strains
        # --------------------------------------------
        # SMAGORINSKY COEFFICIENT COMPONENTS
        # --------------------------------------------
        S11 = dudx
        S22 = dvdy
        S33 = dwdz
        S12 = (dudy + dvdx) / 2
        S13 = (dudz + dwdx) / 2
        S23 = (dvdz + dwdy) / 2
        # --------------------------------------------
        # SMAGORINSKY COEFFICIENT COMPONENTS
        # --------------------------------------------
        # FIXME: Grab functions from module SubgridScaleTurbulence 
        SijSij = S11^2 + S22^2 + S33^2 + 2S12^2 + 2S13^2 + 2S23^2

        Dij = SMatrix{3,3}(dudx, dudy, dudz, 
                           dvdx, dvdy, dvdz,
                           dwdx, dwdy, dwdz)
        
        DF = eltype(SijSij)
        αij = MMatrix{3,3}(DF(0) ,DF(0), DF(0),
                           DF(0), DF(0), DF(0),
                           DF(0), DF(0), DF(0))
        
        DISS = sum(Dij .^ 2)
        βij = similar(αij)
        βij = Δsqr * (Dij' * Dij)
        
        Bβ = βij[1,1]*βij[2,2] - βij[1,2]^2 + βij[1,1]*βij[3,3] - βij[1,3]^2 + βij[2,2]*βij[3,3] - βij[2,3]^2 
        ν_e = max(0,(C_smag^2 * 2.5) * sqrt(abs(Bβ/(DISS+1e-16)))) 
        VF[_ν_e] = ν_e
        #--------------------------------------------
        # deviatoric stresses
        # Fix up index magic numbers
        VF[_τ11] = -2 * (S11 - (S11 + S22 + S33) / 3)
        VF[_τ22] = -2 * (S22 - (S11 + S22 + S33) / 3)
        VF[_τ33] = -2 * (S33 - (S11 + S22 + S33) / 3)
        VF[_τ12] = -2 * S12
        VF[_τ13] = -2 * S13
        VF[_τ23] = -2 * S23

        # TODO: Viscous stresse come from SubgridScaleTurbulence module
        VF[_qx],    VF[_qy],    VF[_qz]    = dqdx,    dqdy,    dqdz
        VF[_Tx],    VF[_Ty],    VF[_Tz]    = dTdx,    dTdy,    dTdz
        VF[_Trefx], VF[_Trefy], VF[_Trefz] = dTrefdx, dTrefdy, dTrefdz
        VF[_Tpx],   VF[_Tpy],   VF[_Tpz]   = dTpdx,   dTpdy,   dTpdz
        VF[_θx],    VF[_θy],    VF[_θz]    = dθdx,    dθdy,    dθdz
        VF[_SijSij] = SijSij
    end
end

# -------------------------------------------------------------------------
#md ### Auxiliary Function (Not required)
#md # In this example the auxiliary function is used to store the spatial
#md # coordinates. This may also be used to store variables for which gradients
#md # are needed, but are not available through teh prognostic variable 
#md # calculations. (An example of this will follow - in the Smagorinsky model, 
#md # where a local Richardson number via potential temperature gradient is required)
# -------------------------------------------------------------------------
@inline function auxiliary_state_initialization!(aux, x, y, z)
    @inbounds begin
        DFloat = eltype(aux)
        xvert = y
        aux[_a_y] = xvert
       
    end
end

# -------------------------------------------------------------------------
# generic bc for 2d , 3d
@inline function bcstate!(QP, VFP, auxP, nM, QM, VFM, auxM, bctype, t)
    DFloat = eltype(QP)
    @inbounds begin
        ρM, UM, VM, WM, EM, QTM = QM[_ρ], QM[_U], QM[_V], QM[_W], QM[_E], QM[_QT]
        uM, vM, wM = UM/ρM, VM/ρM, WM/ρM
        
        UnM = nM[1] * UM + nM[2] * VM + nM[3] * WM
        QP[_U] = UM - 2 * nM[1] * UnM
        QP[_V] = VM - 2 * nM[2] * UnM
        QP[_W] = WM - 2 * nM[3] * UnM
        q1 = (uM*nM[2] - vM*nM[1])*nM[2]
        q2 = (vM*nM[1] - uM*nM[2])*nM[1]
        #QP[_ρ] = ρM
        #QP[_U] = ρM*q1
        #QP[_V] = ρM*q2
        QP[_E] = QM[_E] + 0.5*ρM*(q1^2 + q2^2)

        #if bctype == 3
        #    #VFP .= 0
        #    QP[_U] = 0.0
        #    QP[_V] = 0.0
        #    QP[_W] = 0.0
        #
        #    T = 300 + 5;
        #    QP[_E] = ρ*internal_energy(T, PhasePartition(DFloat(0)))
        #end
        
        #if bctype == 3 ||  bctype == 4
        #    x = auxM[_a_x]
        #    y = auxM[_a_y]
        #    dTrefdy = -grav/cp_d
        #    Pr = 1/3
        #    #VFP[_Trefy] = cp_d * dTrefdy / Pr
        #    VFP[_Trefy] = VFM[_Trefy]
        #end
            
        #QP[_ρ] = ρM
        #QP[_QT] = QTM
        #=
        # noslip boundary 
        QP[_U] = 0
        QP[_V] = 0
        QP[_W] = 0 
        QP[_ρ] = ρM 
        QP[_QT] = QTM 
        =#

        #Dirichlet on \theta at bottom bvoundary 
        #=if bctype == 3
            y       = auxM[_a_y]
            θ_ref   = 300.0
            θ_c     = 10.0
            Δθ      = θ_c + θ_ref #* sin(pi*x.*5/Lx).*cos(x/pi);
            p0      = MSLP
            θ       = Δθ # potential temperature
            π_exner = 1.0 - grav / (cp_d * θ) * y # exner pressure
            ρ       = p0 / (R_d * θ) * (π_exner)^ (cv_d / R_d) # density
            P       = p0 * (R_d * (ρ * θ) / p0) ^(cp_d/cv_d) # pressure (absolute)
            T       = P / (ρ * R_d) # temperature
            
            QP[_E]  = ρ*internal_energy(T, PhasePartition(DFloat(0))) + 0.5*ρ*(u^2 + v^2 + w^2)
            #elseif bctype == 4
            #    y = auxM[_a_y]
            #    T = Tc
            #    QP[_E] = internal_energy(T, PhasePartition(DFloat(0))) + grav * y 
        end=#
        nothing
    end
end

# -------------------------------------------------------------------------
@inline function stresses_boundary_penalty!(VF, _...) 
    VF .= 0
end

@inline function stresses_penalty!(VF, nM, gradient_listM, QM, aM, gradient_listP, QP, aP, t)
    @inbounds begin
        n_Δgradient_list = similar(VF, Size(3, _ngradstates))
        for j = 1:_ngradstates, i = 1:3
            n_Δgradient_list[i, j] = nM[i] * (gradient_listP[j] - gradient_listM[j]) / 2
        end
        compute_stresses!(VF, n_Δgradient_list)
    end
end
# -------------------------------------------------------------------------

@inline function source!(S,Q,aux,t)
    # Initialise the final block source term 
    S .= 0

    # Typically these sources are imported from modules
    @inbounds begin
        source_geopot!(S, Q, aux, t)
        #source_sponge!(S, Q, aux, t)
    end
end

"""
        Geostrophic wind forcing
    """

@inline function source_sponge!(S,Q,aux,t)
    @inbounds begin
        U, V, W  = Q[_U], Q[_V], Q[_W]
        beta     = aux[_a_sponge]
        #S[_U] -= beta * U
        S[_V] -= beta * V             
    end
end

@inline function source_geopot!(S,Q,aux,t)
    @inbounds S[_V] += - Q[_ρ] * grav
end

# Test integral exactly according to the isentropic vortex example
@inline function integrand(val, Q, aux)
    κ = 85.0
    @inbounds begin
        @inbounds ρ, QT = Q[_ρ], Q[_QT]
        ρinv = 1 / ρ
        q_tot = QT * ρinv * 0.0
        # Establish the current thermodynamic state using the prognostic variables
        q_liq = aux[_a_q_liq]
        val[1] = ρ * κ * q_liq 
        val[2] = ρ * q_liq # Liquid Water Path Integrand
    end
end

function preodefun!(disc, Q, t)
    DGBalanceLawDiscretizations.dof_iteration!(disc.auxstate, disc, Q) do R, Q, QV, aux
        @inbounds let
            ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
            xvert = aux[_a_y]
            e_int = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * grav * xvert) / ρ
            q_tot = QT / ρ * 0.0 
            TS = PhaseEquil(e_int, q_tot, ρ)
            T = air_temperature(TS)
            P = air_pressure(TS) # Test with dry atmosphere
            q_liq = PhasePartition(TS).liq

            R[_a_T] = T
            R[_a_P] = P
            R[_a_q_liq] = q_liq
            R[_a_soundspeed_air] = soundspeed_air(TS)
            R[_a_θ] = virtual_pottemp(TS)
        end
    end

    integral_computation(disc, Q, t)
end

function integral_computation(disc, Q, t)
    DGBalanceLawDiscretizations.indefinite_stack_integral!(disc, integrand, Q,
                                                           (_a_02z, _a_LWP_02z))
    
    DGBalanceLawDiscretizations.reverse_indefinite_stack_integral!(disc,
                                                                   _a_z2inf,
                                                                   _a_02z)

    DGBalanceLawDiscretizations.reverse_indefinite_stack_integral!(disc,
                                                                   _a_LWP_z2inf,
                                                                   _a_LWP_02z)
end

# initial condition
"""
        User-specified. Required.
        This function specifies the initial conditions
        for the dry_benchmark driver.
    """
const seed = MersenneTwister(0)

# NEW FUNCTION
function dry_benchmark!(dim, Q, t, x, y, z, _...)
    DFloat                = eltype(Q)
    
    randnum1   = rand(seed, DFloat) / 100
    randnum2   = rand(seed, DFloat) / 100
    
    # can override default gas constants 
    # to moist values later in the driver 
    R_gas::DFloat         = R_d
    c_p::DFloat           = cp_d
    c_v::DFloat           = cv_d
    p0::DFloat            = MSLP
    gravity::DFloat       = grav
    # initialise with dry domain 
    q_tot::DFloat         = 0 
    q_liq::DFloat         = 0
    q_ice::DFloat         = 0

    
    # perturbation parameters for rising bubble
    xc                    = 0.0
    yc                    = 2100.0
    r                     = sqrt((x - xc)^2 + (y - yc)^2)
    rc::DFloat            = 2000.0
    θ_ref::DFloat         =  300.0
    θ_c::DFloat           =    -2.0
    Δθ::DFloat            =    0.0
    if r <= rc 
        Δθ = θ_c * (1 - (r/rc)) + randnum1*Δθ
    end

    #Reference state:
    π_exner_ref               = 1.0 - gravity / (c_p * θ_ref) * y # exner pressure
    ρ_ref                     = p0 / (R_gas * θ_ref) * (π_exner_ref)^ (c_v / R_gas) # density
    P_ref                     = p0 * (R_gas * (ρ_ref * θ_ref) / p0) ^(c_p/c_v) # pressure (absolute)
    T_ref                     = P_ref / (ρ_ref * R_gas) # temperature    
    
    #Perturbed state:
    θ                     = θ_ref + Δθ # potential temperature
    π_exner               = 1.0 - gravity / (c_p * θ) * y # exner pressure
    ρ                     = p0 / (R_gas * θ) * (π_exner)^ (c_v / R_gas) # density
    #P                     = T * ρ * R_gas
    P                     = p0 * (R_gas * (ρ * θ) / p0) ^(c_p/c_v) # pressure (absolute)
    T                     = P / (ρ * R_gas) # temperature
    
    U, V, W               = 0.0 , 0.0 , 0.0  # momentum components
    # energy definitions
    e_kin                 = (U^2 + V^2 + W^2) / (2*ρ)/ ρ
    e_pot                 = gravity * y
    e_int                 = c_v*(T - T_0)
    E                     = ρ * (e_int + e_kin + e_pot)  #* total_energy(e_kin, e_pot, T, q_tot, q_liq, q_ice)
    @inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]= ρ, U, V, W, E, ρ * q_tot
end

function run(mpicomm, dim, Ne, N, timeend, DFloat, dt)

    #current stretching options:
    # "boundary_stretching" --> node clustering by the bottom boundary
    # "top_stretching" --> node clustering by the top wall
    # "interior_stretching" --> node clustering around the coordinate of the `attractor_value`
    
    y_range     = grid_stretching_1d(ymin, ymax, Ne[end], stretch_coe, "no_boundary_stretching")    
    brickrange  = (range(DFloat(xmin), length=Ne[1]+1, DFloat(xmax)),
                   y_range)
    
    # User defined periodicity in the topl assignment
    # brickrange defines the domain extents
    @timeit to "Topo init" topl = StackedBrickTopology(mpicomm, brickrange, periodicity=(false,false))

    @timeit to "Grid init" grid = DiscontinuousSpectralElementGrid(topl,
                                                                   FloatType = DFloat,
                                                                   DeviceArray = ArrayType,
                                                                   polynomialorder = N)

    numflux!(x...) = NumericalFluxes.rusanov!(x..., cns_flux!, wavespeed)
    numbcflux!(x...) = NumericalFluxes.rusanov_boundary_flux!(x..., cns_flux!, bcstate!, wavespeed)

    # spacedisc = data needed for evaluating the right-hand side function
    @timeit to "Space Disc init" spacedisc = DGBalanceLaw(grid = grid,
                                                          length_state_vector = _nstate,
                                                          flux! = cns_flux!,
                                                          numerical_flux! = numflux!,
                                                          numerical_boundary_flux! = numbcflux!, 
                                                          number_gradient_states = _ngradstates,
                                                          number_viscous_states = _nviscstates,
                                                          gradient_transform! = gradient_vars!,
                                                          viscous_transform! = compute_stresses!,
                                                          viscous_penalty! = stresses_penalty!,
                                                          viscous_boundary_penalty! = stresses_boundary_penalty!,
                                                          auxiliary_state_length = _nauxstate,
                                                          auxiliary_state_initialization! = (x...) ->
                                                          auxiliary_state_initialization!(x...),
                                                          source! = source!,
                                                          preodefun! = preodefun!)

    # This is a actual state/function that lives on the grid
    @timeit to "IC init" begin
        # ----------------------------------------------------
        # GET DATA FROM INTERPOLATED ARRAY ONTO VECTORS
        # This driver accepts data in 6 column format
        # ----------------------------------------------------
        (sounding, _, ncols) = read_sounding()

        # WARNING: Not all sounding data is formatted/scaled
        # the same. Care required in assigning array values
        # height theta qv    u     v     pressure
        zinit, tinit, pinit = sounding[:, 1], sounding[:, 2], sounding[:, 3]
        thetainit, qinit = sounding[:, 4], sounding[:, 5]
        #------------------------------------------------------
        # GET SPLINE FUNCTION
        #------------------------------------------------------
        spl_tinit    = Spline1D(zinit, tinit; k=1) #sensible T (K)
        spl_pinit    = Spline1D(zinit, pinit; k=1) #pressure P (Pa)
        
        spl_thetainit= Spline1D(zinit, thetainit; k=1) #sensible T (K)
        spl_qinit    = Spline1D(zinit, qinit; k=1) #sensible T (K)
        
        # Set type of filter 
        filter_dry_benchmark = CLIMA.Mesh.Grids.CutoffFilter(spacedisc.grid)

        initialcondition(Q, x...) = dry_benchmark!(Val(dim), Q, DFloat(0), x...)
        Q = MPIStateArray(spacedisc, initialcondition)
    end

    @timeit to "Time stepping init" begin
        lsrk = LSRK54CarpenterKennedy(spacedisc, Q; dt = dt, t0 = 0)
        
        # Set up the information callback
        starttime = Ref(now())
        cbinfo = GenericCallbacks.EveryXWallTimeSeconds(10, mpicomm) do (s=false)
            if s
                starttime[] = now()
            else
                maxVertVelo = global_max(Q, _V)
                @info @sprintf("""Update
                               simtime = %.16e
                               runtime = %s
                               max(vert) = %.16e
                               """,
                               ODESolvers.gettime(lsrk),
                               Dates.format(convert(Dates.DateTime,
                                                    Dates.now()-starttime[]),
                                            Dates.dateformat"HH:MM:SS"),
                               maxVertVelo
                               )
            end
        end
        
        npoststates = 5
        _o_u, _o_v, _o_w, _o_T, _o_θ= 1:npoststates
        postnames = ("u", "v", "w", "T", "THETA")
        postprocessarray = MPIStateArray(spacedisc; nstate=npoststates)

        cbfilter = GenericCallbacks.EveryXSimulationSteps(1) do
            DGBalanceLawDiscretizations.apply!(Q, 1:_nstate, spacedisc,
                                               filter_dry_benchmark;
                                               horizontal=true,
                                               vertical=true)
            nothing
        end
        
        step = [0]
        cbvtk = GenericCallbacks.EveryXSimulationSteps(100) do (init=false)
            DGBalanceLawDiscretizations.dof_iteration!(postprocessarray, spacedisc, Q) do R, Q, QV, aux
                @inbounds let
                    u, v, w     = preflux(Q, aux)
                    R[_o_u]     = u
                    R[_o_v]     = v
                    R[_o_w]     = w
                    R[_o_T]     = aux[_a_T]
                    R[_o_θ]     = aux[_a_θ]
                end
            end
            
            mkpath("./CLIMA-output-scratch/ahmad/")
            outprefix = @sprintf("./CLIMA-output-scratch/ahmad/ahmad_%dD_mpirank%04d_step%04d", dim,
                                 MPI.Comm_rank(mpicomm), step[1])
            @debug "doing VTK output" outprefix
            writevtk(outprefix, Q, spacedisc, statenames,
                     postprocessarray, postnames)
            
            step[1] += 1
            nothing
        end
    end

@info @sprintf """Starting...
        norm(Q) = %25.16e""" norm(Q)

# Initialise the integration computation. Kernels calculate this at every timestep?? 
@timeit to "initial integral" integral_computation(spacedisc, Q, 0) 
@timeit to "solve" solve!(Q, lsrk; timeend=timeend, callbacks=(cbinfo, cbvtk))


@info @sprintf """Finished...
        norm(Q) = %25.16e""" norm(Q)

end

using Test
let
    MPI.Initialized() || MPI.Init()
    Sys.iswindows() || (isinteractive() && MPI.finalize_atexit())
    mpicomm = MPI.COMM_WORLD
    ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
    loglevel = ll == "DEBUG" ? Logging.Debug :
        ll == "WARN"  ? Logging.Warn  :
        ll == "ERROR" ? Logging.Error : Logging.Info
    logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
    global_logger(ConsoleLogger(logger_stream, loglevel))
    @static if haspkg("CUDAnative")
        device!(MPI.Comm_rank(mpicomm) % length(devices()))
    end
    # User defined number of elements
    # User defined timestep estimate
    # User defined simulation end time
    # User defined polynomial order 
    numelem = (Nex, Ney)
    dt = 0.005
    timeend = 1020
    polynomialorder = Npoly
    DFloat = Float64
    dim = numdims

    if MPI.Comm_rank(mpicomm) == 0
        @info @sprintf """ ------------------------------------------------------"""
        @info @sprintf """   ______ _      _____ __  ________                    """     
        @info @sprintf """  |  ____| |    |_   _|  ...  |  __  |                 """  
        @info @sprintf """  | |    | |      | | |   .   | |  | |                 """ 
        @info @sprintf """  | |    | |      | | | |   | | |__| |                 """
        @info @sprintf """  | |____| |____ _| |_| |   | | |  | |                 """
        @info @sprintf """  | _____|______|_____|_|   |_|_|  |_|                 """
        @info @sprintf """                                                       """
        @info @sprintf """ ------------------------------------------------------"""
        @info @sprintf """ Rayleigh-Benard                                       """
        @info @sprintf """   Resolution:                                         """ 
        @info @sprintf """     (Δx, Δy, Δz)   = (%.2e, %.2e, %.2e)               """ Δx Δy Δz
        @info @sprintf """     (Nex, Ney, Nez) = (%d, %d, %d)                    """ Nex Ney Nez
        @info @sprintf """     DoF = %d                                          """ DoF
        @info @sprintf """     Minimum necessary memory to run this test: %g GBs """ (DoFstorage * sizeof(DFloat))/1000^3
        @info @sprintf """     Time step dt: %.2e                                """ dt
        @info @sprintf """     End time  t : %.2e                                """ timeend
        @info @sprintf """ ------------------------------------------------------"""
    end

    engf_eng0 = run(mpicomm, dim, numelem[1:dim], polynomialorder, timeend,
                    DFloat, dt)

    show(to)
end

isinteractive() || MPI.Finalize()

nothing
