
# Load Modules 
using MPI
using CLIMA
using CLIMA.Topologies
using CLIMA.Grids
using CLIMA.DGBalanceLawDiscretizations
using CLIMA.DGBalanceLawDiscretizations.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.ODESolvers
using CLIMA.GenericCallbacks
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates
using CLIMA.Vtk
using DelimitedFiles
using Dierckx

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
using CLIMA.PlanetParameters: R_d, cp_d, grav, cv_d, MSLP, T_0, Omega

# State labels 
const _nstate = 6
const _ρ, _U, _V, _W, _E, _QT = 1:_nstate
const stateid = (ρid = _ρ, Uid = _U, Vid = _V, Wid = _W, Eid = _E, QTid = _QT)
const statenames = ("RHO", "U", "V", "W", "E", "QT")

# Viscous state labels
const _nviscstates = 16
const _τ11, _τ22, _τ33, _τ12, _τ13, _τ23, _qx, _qy, _qz, _Tx, _Ty, _Tz, _θx, _θy, _θz, _SijSij = 1:_nviscstates

# Gradient state labels
const _states_for_gradient_transform = (_ρ, _U, _V, _W, _E, _QT)


const _nauxstate = 10
const _a_x, _a_y, _a_z, _a_sponge, _a_02z, _a_z2inf, _a_rad, _a_ν_e, _a_LWP_02z, _a_LWP_z2inf = 1:_nauxstate


if !@isdefined integration_testing
    const integration_testing =
        parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
    using Random
end

# Problem constants (TODO: parameters module (?))
const μ_sgs           = 100.0
const Prandtl         = 71 // 100
const Prandtl_t       = 1 // 3
const Ri_c            = 1 // 3
const cp_over_prandtl = cp_d / Prandtl_t

# Problem description 
# --------------------
# 2D thermal perturbation (cold bubble) in a neutrally stratified atmosphere
# No wall-shear, lateral periodic boundaries with no-flux walls at the domain
# top and bottom. 
# Inviscid, Constant viscosity, StandardSmagorinsky, MinimumDissipation
# filters are tested against this benchmark problem
# TODO: link to module SubGridScaleTurbulence

# Random number seed
const seed = MersenneTwister(0)

#
# User Input
#
const numdims = 3
const Npoly = 4

#
# Define grid size 
#
Δx    = 35
Δy    = 35
Δz    = 10
#
# OR:
#
# Set Δx < 0 and define  Nex, Ney, Nez:
#
(Nex, Ney, Nez) = (10, 10, 1)

# Physical domain extents 
const (xmin, xmax) = (0, 840)
const (ymin, ymax) = (0, 840)
const (zmin, zmax) = (0, 1500)

#Get Nex, Ney from resolution
const Lx = xmax - xmin
const Ly = ymax - ymin
const Lz = zmax - ymin

if ( Δx > 0)
    #
    # User defines the grid size:
    #
    ratiox = (Lx/Δx - 1)/Npoly
    ratioy = (Ly/Δy - 1)/Npoly
    ratioz = (Lz/Δz - 1)/Npoly
    Nex = ceil(Int64, ratiox)
    Ney = ceil(Int64, ratioy)
    Nez = ceil(Int64, ratioz)
    
else
    #
    # User defines the number of elements:
    #
    Δx = Lx / ((Nex * Npoly) + 1)
    Δy = Ly / ((Ney * Npoly) + 1)
    Δz = Lz / ((Nez * Npoly) + 1)
end

DoF = (Nex*Ney*Nez)*(Npoly+1)^numdims*(_nstate)
DoFstorage = (Nex*Ney*Nez)*(Npoly+1)^numdims*(_nstate + _nviscstates + _nauxstate + CLIMA.Grids._nvgeo) +
    (Nex*Ney*Nez)*(Npoly+1)^(numdims-1)*2^numdims*(CLIMA.Grids._nsgeo)


# Smagorinsky model requirements : TODO move to SubgridScaleTurbulence module 
const C_smag = 0.15
# Equivalent grid-scale
Δ = (Δx * Δy * Δz)^(1/3)
const Δsqr = Δ * Δ


function deardorff_zero_equation_model(modSij, θv, dθvdz, Δ)

    #
    # 0-equation parameterzied Deardorff model by Kirkpatrick et al JAS 2006
    #
    # Brunt-Vaisala frequency
    N2 = grav / θv * dθvdz
    
    # Ri number
    Ri = N2 / (2 * normSij + eps(normSij))     
    # Buoyancy correction factor
     
     ca  = 0.10
     ce1 = 0.19
     ce2 = 0.51
     c1  = ca*0.75^2
     c2  = ce2 + 2*c1

     CB = 3*Ri >= 1 ? 0 :  sqrt(1 - 3*Ri)
     L  = N2 <=0 ? Δ : Δ*(c1*(1/Ri - 1) - ce1)/c2
     
     cϵ = ce1 + ce2*buoyancy_factor
     c3 = sqrt(ca^3)

    Km = L <= 0 ? 0 : c3*L^2*(2*modSik)*CB / sqrt(cϵ)
     
    return Km
  end

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
@inline function preflux(Q,VF, aux, _...)
    R_gas::eltype(Q) = R_d
    @inbounds ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
    ρinv = 1 / ρ
    q_tot = ρinv * QT
    
    x, y, z = aux[_a_x], aux[_a_y], aux[_a_z]
    u, v, w = ρinv * U, ρinv * V, ρinv * W

    #energies:
    e_tot = ρinv * E
    e_pot = grav * z
    e_kin = 0.5 * (u^2 + v^2 + w^2)
    e_int = e_tot - e_kin - e_pot
    
    
    # Establish the current thermodynamic state using the prognostic variables
    TS = PhaseEquil(e_int, q_tot, ρ)
   
    #Update E and obtain q_liq at current time step:
    E  = ρ * total_energy(e_kin, e_pot, TS)
   
    #Diagnose q_liq and θv for posprocessing:
    P     = air_pressure(TS) # Test with dry atmosphere
    T     = air_temperature(TS)
    q_liq = PhasePartition(TS).liq
    θ     = dry_pottemp(TS)
    θv    = virtual_pottemp(TS)

    #Return:
    (P, u, v, w, ρinv, q_liq, T, θ)
end

#-------------------------------------------------------------------------
#md # Soundspeed computed using the thermodynamic state TS
# max eigenvalue
@inline function wavespeed(n, Q, aux, t, P, u, v, w, ρinv, q_liq, T, θ)
    @inbounds begin 
        ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
        x,y,z = aux[_a_x], aux[_a_y], aux[_a_z]
        u, v, w = ρinv * U, ρinv * V, ρinv * W
        e_int = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * grav * z) / ρ
        q_tot = QT / ρ
        TS = PhaseEquil(e_int, q_tot, ρ)
        (n[1] * u + n[2] * v + n[3] * w) + soundspeed_air(TS)
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

    #fsounding  = open(joinpath(@__DIR__, "../soundings/sounding_DYCOMS_TEST1.dat"))
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
cns_flux!(F, Q, VF, aux, t) = cns_flux!(F, Q, VF, aux, t, preflux(Q,VF, aux)...)
@inline function cns_flux!(F, Q, VF, aux, t, P, u, v, w, ρinv, q_liq, T, θ)
    @inbounds begin
        ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
        # Inviscid contributions
        F[1, _ρ], F[2, _ρ], F[3, _ρ] = U          , V          , W
        F[1, _U], F[2, _U], F[3, _U] = u * U  + P , v * U      , w * U
        F[1, _V], F[2, _V], F[3, _V] = u * V      , v * V + P  , w * V
        F[1, _W], F[2, _W], F[3, _W] = u * W      , v * W      , w * W + P
        F[1, _E], F[2, _E], F[3, _E] = u * (E + P), v * (E + P), w * (E + P)
        F[1, _QT], F[2, _QT], F[3, _QT] = u * QT  , v * QT     , w * QT 

        #Derivative of T and Q:
        vqx, vqy, vqz = VF[_qx], VF[_qy], VF[_qz]
        vTx, vTy, vTz = VF[_Tx], VF[_Ty], VF[_Tz]
        vθz = VF[_θz]
      
        # Radiation contribution 
        F_rad = ρ * radiation(aux)  
        aux[_a_rad] = F_rad
       
        SijSij = VF[_SijSij]
        f_R = 1.0# buoyancy_correction_smag(SijSij, θ, dθdy)

        #Dynamic eddy viscosity from Smagorinsky:
        ν_e = sqrt(2.0 * SijSij) * C_smag^2 * Δsqr
        
        #N2 = grav * vθz / 289
        #Ri = N2 / (2 * SijSij + eps(SijSij))
        #buoyancy_factor = N2 <=0 ? 1 : max(0.0, 1 - Ri/Ri_c)
        #ν_e = C_smag^2 * Δsqr * sqrt(2*SijSij * buoyancy_factor)
        D_e = ν_e / Prandtl_t
        aux[_a_ν_e] = ν_e
        
        # Multiply stress tensor by viscosity coefficient:
        τ11, τ22, τ33 = VF[_τ11] * ν_e, VF[_τ22]* ν_e, VF[_τ33] * ν_e
        τ12 = τ21 = VF[_τ12] * ν_e 
        τ13 = τ31 = VF[_τ13] * ν_e               
        τ23 = τ32 = VF[_τ23] * ν_e
        
        # Viscous velocity flux (i.e. F^visc_u in Giraldo Restelli 2008)
        F[1, _U] -= τ11 * f_R ; F[2, _U] -= τ12 * f_R ; F[3, _U] -= τ13 * f_R
        F[1, _V] -= τ21 * f_R ; F[2, _V] -= τ22 * f_R ; F[3, _V] -= τ23 * f_R
        F[1, _W] -= τ31 * f_R ; F[2, _W] -= τ32 * f_R ; F[3, _W] -= τ33 * f_R
        
        # Viscous Energy flux (i.e. F^visc_e in Giraldo Restelli 2008)
        F[1, _E] -= u * τ11 + v * τ12 + w * τ13 + cp_over_prandtl * vTx * ν_e
        F[2, _E] -= u * τ21 + v * τ22 + w * τ23 + cp_over_prandtl * vTy * ν_e
        F[3, _E] -= u * τ31 + v * τ32 + w * τ33 + cp_over_prandtl * vTz * ν_e
        
        F[1, _E] -= 0
        F[2, _E] -= 0
        F[3, _E] -= F_rad
        # Viscous contributions to mass flux terms
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
gradient_vars!(vel, Q, aux, t, _...) = gradient_vars!(vel, Q, aux, t, preflux(Q,~,aux)...)

const _ngradstates = 7
@inline function gradient_vars!(vel, Q, aux, t, P, u, v, w, ρinv, q_liq, T, θ)
  @inbounds begin
      # ordering should match states_for_gradient_transform
    ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
    vel[1], vel[2], vel[3] = u, v, w
    vel[4], vel[5], vel[6] = E, QT/ρ, T
    vel[7] = θ
  end
end

@inline function radiation(aux)
  zero_to_z = aux[_a_02z]
  z_to_inf = aux[_a_z2inf]
  z = aux[_a_z]
  z_i = 840  # Start with constant inversion height of 840 meters then build in check based on q_tot
  (z - z_i) >=0 ? Δz_i = (z - z_i) : Δz_i = 0 
  # Constants 
  F_0 = 70 
  F_1 = 22
  α_z = 1
  ρ_i = 1.22
  D_subsidence = 3.75e-6
  term1 = F_0 * exp(-z_to_inf) 
  term2 = F_1 * exp(-zero_to_z)
  term3 = ρ_i * cp_d * D_subsidence * α_z * (0.25 * (cbrt(Δz_i))^4 + z_i * cbrt(Δz_i))
  F_rad = term1 + term2 + term3  
  return F_rad
end

# -------------------------------------------------------------------------
#md ### Viscous fluxes. 
#md # The viscous flux function compute_stresses computes the components of 
#md # the velocity gradient tensor, and the corresponding strain rates to
#md # populate the viscous flux array VF. SijSij is calculated in addition
#md # to facilitate implementation of the constant coefficient Smagorinsky model
#md # (pending)
@inline function compute_stresses!(VF, grad_vel, _...)
  @inbounds begin
    dudx, dudy, dudz = grad_vel[1, 1], grad_vel[2, 1], grad_vel[3, 1]
    dvdx, dvdy, dvdz = grad_vel[1, 2], grad_vel[2, 2], grad_vel[3, 2]
    dwdx, dwdy, dwdz = grad_vel[1, 3], grad_vel[2, 3], grad_vel[3, 3]
    # compute gradients of moist vars and temperature
    dqdx, dqdy, dqdz = grad_vel[1, 5], grad_vel[2, 5], grad_vel[3, 5]
    dTdx, dTdy, dTdz = grad_vel[1, 6], grad_vel[2, 6], grad_vel[3, 6]
    dθdx, dθdy, dθdz = grad_vel[1, 7], grad_vel[2, 7], grad_vel[3, 7]
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
    SijSij = (S11^2 + S22^2 + S33^2
              + 2.0 * S12^2
              + 2.0 * S13^2 
              + 2.0 * S23^2) 
    modSij = sqrt(2.0 * SijSij)
    
    #--------------------------------------------
    # deviatoric stresses
    # Fix up index magic numbers
    VF[_τ11] = 2 * (S11 - (S11 + S22 + S33) / 3)
    VF[_τ22] = 2 * (S22 - (S11 + S22 + S33) / 3)
    VF[_τ33] = 2 * (S33 - (S11 + S22 + S33) / 3)
    VF[_τ12] = 2 * S12
    VF[_τ13] = 2 * S13
    VF[_τ23] = 2 * S23
    
    # TODO: Viscous stresse come from SubgridScaleTurbulence module
    VF[_qx], VF[_qy], VF[_qz] = dqdx, dqdy, dqdz
    VF[_Tx], VF[_Ty], VF[_Tz] = dTdx, dTdy, dTdz
    VF[_θx], VF[_θy], VF[_θz] = dθdx, dθdy, dθdz
    VF[_SijSij] = SijSij
  end
end
# -------------------------------------------------------------------------
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
        aux[_a_x] = x
        aux[_a_y] = y
        aux[_a_z] = z
        
        #Vertical sponge
        ctop        = 0.0
        ct          = 0.5
        zd          = 300.0       
        sponge_type = 1
        
        top_sponge  = zmax - zd
        if sponge_type == 1
            if z >= top_sponge
                ctop = ct * (sinpi(0.5 * (z - top_sponge)/(zmax - top_sponge)))^4
            end
        elseif sponge_type == 2
            ct          = 0.25
            if z >= top_sponge
                ctop = ct * sinpi(0.5 * (1.0 - (zmax - z)/zd))^2.0
            end
        end
        
        beta  = 1.0 - (1.0 - ctop) #*(1.0 - csleft)*(1.0 - csright)*(1.0 - csfront)*(1.0 - csback)
        beta  = min(beta, 1.0)
        aux[_a_sponge] = beta
        
    end
end

# -------------------------------------------------------------------------
# generic bc for 2d , 3d

@inline function bcstate!(QP, VFP, auxP, nM, QM, VFM, auxM, bctype, t, PM, uM, vM, wM, ρinvM, q_liqM, TM, θM)
    @inbounds begin
        x, y, z = auxM[_a_x], auxM[_a_y], auxM[_a_z]
        ρM, UM, VM, WM, EM, QTM = QM[_ρ], QM[_U], QM[_V], QM[_W], QM[_E], QM[_QT]
        # No flux boundary conditions
        # No shear on walls (free-slip condition)
        UnM = nM[1] * UM + nM[2] * VM + nM[3] * WM
        QP[_U] = UM - 2 * nM[1] * UnM
        QP[_V] = VM - 2 * nM[2] * UnM
        QP[_W] = WM - 2 * nM[3] * UnM
        #QP[_ρ] = ρM
        #QP[_QT] = QTM
        VFP .= 0

        
        
        nothing
    end
end

# -------------------------------------------------------------------------
@inline function stresses_boundary_penalty!(VF, _...) 
    VF .= 0
end

@inline function stresses_penalty!(VF, nM, velM, QM, aM, velP, QP, aP, t)
    @inbounds begin
        n_Δvel = similar(VF, Size(3, _ngradstates))
        for j = 1:_ngradstates, i = 1:3
            n_Δvel[i, j] = nM[i] * (velP[j] - velM[j]) / 2
        end
        compute_stresses!(VF, n_Δvel)
    end
end
# -------------------------------------------------------------------------

@inline function source!(S,Q,aux,t)
    # Initialise the final block source term 
    S .= 0

    # Typically these sources are imported from modules
    @inbounds begin
        source_geopot!(S, Q, aux, t)
        source_sponge!(S, Q, aux, t)
        source_coriolis!(S, Q, aux, t)
        source_geostrophic!(S, Q, aux, t)
    end
end

"""
Coriolis force
"""
const f_coriolis = 7.62e-5
const u_geostrophic = 7.0
const v_geostrophic = -5.5 
const Ω = Omega
@inline function source_coriolis!(S,Q,aux,t)
  @inbounds begin
    U, V, W = Q[_U], Q[_V], Q[_W]
    S[_U] -= 0
    S[_V] -= 0
    S[_W] -= 0
  end
end

"""
    Geostrophic wind forcing
"""
@inline function source_geostrophic!(S,Q,aux,t)
    @inbounds begin
        ρ = Q[_ρ]
        U = Q[_U]
        V = Q[_V]        
        W = Q[_W]
        
        S[_U] -= f_coriolis * (U - ρ*u_geostrophic)
        S[_V] -= f_coriolis * (V - ρ*v_geostrophic)
    end
end

@inline function source_sponge!(S,Q,aux,t)
    @inbounds begin
        U, V, W  = Q[_U], Q[_V], Q[_W]        
        beta     = aux[_a_sponge]
        S[_U] -= beta * U
        S[_V] -= beta * V
        S[_W] -= beta * W
        
    end
end
#---END SPONGE

@inline function source_geopot!(S,Q,aux,t)
    @inbounds begin
        ρ, U, V, W, E  = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E]
        S[_W] += - ρ * grav
    end
end

# Test integral exactly according to the isentropic vortex example
@inline function integral_knl(val, Q, aux)
  κ = 85.0
  @inbounds begin
    @inbounds ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
    ρinv = 1 / ρ
    x,y,z = aux[_a_x], aux[_a_y], aux[_a_z]
    u, v, w = ρinv * U, ρinv * V, ρinv * W
    e_int = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * grav * z) / ρ
    q_tot = QT / ρ
    # Establish the current thermodynamic state using the prognostic variables
    TS = PhaseEquil(e_int, q_tot, ρ)
    q_liq = PhasePartition(TS).liq
      
    val[1] = ρ * κ * (q_liq / (1.0 - q_tot))
    val[2] = ρ * (q_liq / (1.0 - q_tot))     #LWP
  end
end

function integral_computation(disc, Q, t)
  DGBalanceLawDiscretizations.indefinite_stack_integral!(disc, integral_knl, Q,
                                                         (_a_02z, _a_LWP_02z))
    
  DGBalanceLawDiscretizations.reverse_indefinite_stack_integral!(disc,
                                                                 _a_z2inf,
                                                                 _a_02z)

  DGBalanceLawDiscretizations.reverse_indefinite_stack_integral!(disc,
                                                                 _a_LWP_z2inf,
                                                                 _a_LWP_02z)
end

# ------------------------------------------------------------------
# -------------END DEF SOURCES-------------------------------------# 

# initial condition
"""
    User-specified. Required. 
    This function specifies the initial conditions
    for the dycoms driver. 
  """
#=
function dycoms!(dim, Q, t, spl_tinit, spl_qinit, spl_uinit, spl_vinit,
                 spl_pinit, x, y, z, _...)
  DFloat         = eltype(Q)
  # --------------------------------------------------
  # INITIALISE ARRAYS FOR INTERPOLATED VALUES
  # --------------------------------------------------
  xvert          = z

  datat          = DFloat(spl_tinit(xvert))
  dataq          = DFloat(spl_qinit(xvert))
  datau          = DFloat(spl_uinit(xvert))
  datav          = DFloat(spl_vinit(xvert))
  datap          = DFloat(spl_pinit(xvert))
  dataq          = dataq / 1000

  randnum1   = rand(seed, DFloat) / 100
  randnum2   = rand(seed, DFloat) / 100

  θ_liq = datat + randnum1 * datat
  q_tot = dataq + randnum2 * dataq
  P     = datap
  T     = air_temperature_from_liquid_ice_pottemp(θ_liq, P, PhasePartition(q_tot))
  ρ     = air_density(T, P)

  # energy definitions
  u, v, w     = datau, datav, zero(DFloat) #geostrophic. TO BE BUILT PROPERLY if Coriolis is considered
  U           = ρ * u
  V           = ρ * v
  W           = ρ * w
  e_kin       = (u^2 + v^2 + w^2) / 2
  e_pot       = grav * xvert
  e_int       = internal_energy(T, PhasePartition(q_tot))
  E           = ρ * total_energy(e_kin, e_pot, T, PhasePartition(q_tot))

  @inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT] = ρ, U, V, W, E, ρ * q_tot
end
=#
function dycoms!(dim, Q, t, spl_tinit, spl_pinit, spl_thetainit, spl_qinit, x, y, z, _...)
    
    DFloat     = eltype(Q)
    p0::DFloat = MSLP
    
    randnum1   = rand(seed, DFloat) / 150
    randnum2   = rand(seed, DFloat) / 150
    
    xvert  = z
    P      = spl_pinit(xvert)     #P
    θ_l    = spl_thetainit(xvert) #θ_l
    q_tot  = spl_qinit(xvert)     #qtot
    T      = spl_tinit(xvert)    #T
        
    zi = 840.0
    if ( z <= zi)
	θ_lx   = 289.0;
	q_totx = 9.0e-3; #specific humidity
    else
	θ_lx   = 297.5 + (z - zi)^(1/3);
	q_totx = 1.5e-3; #kg/kg  specific humidity --> approx. to mixing ratio is ok
    end  
   
    q_liq = 0.0
    if z >= 600.0 && z <= 840.0
        q_liq = (z - 600)*0.00045/200.0 
    end
    θ_l   = θ_l   + randnum1 * θ_l
    #q_tot = q_tot + randnum2 * q_tot
    q_liq = q_liq + randnum2 * q_liq
    
    q_partition = PhasePartition(q_tot, q_liq, 0.0)
    e_int  = internal_energy(T, q_partition)
    
    #TS = LiquidIcePotTempSHumEquil_no_ρ(θ_l, q_tot, P, T)
    #ρ  = air_density(TS)
    ρ  = air_density(T, P, q_partition)

    #u, v
    u, v, w = 7.0, -5.5, 0.0 #geostrophic. TO BE BUILT PROPERLY if Coriolis is considered
    
    e_kin = (u^2 + v^2 + w^2) / 2
    e_pot = grav * xvert
    E     = ρ * (e_int + e_kin + e_pot)

    U, V, W = ρ * u, ρ * v, ρ * w
    
    @inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]= ρ, U, V, W, E, ρ * q_tot
    #try the filter
end


# ------------------------------------------------------------------
# -------------END DEF SOURCES-------------------------------------# 

function run(mpicomm, dim, Ne, N, timeend, DFloat, dt)

    brickrange = (range(DFloat(xmin), length=Ne[1]+1, DFloat(xmax)),
                  range(DFloat(ymin), length=Ne[2]+1, DFloat(ymax)),
                  range(DFloat(zmin), length=Ne[3]+1, DFloat(zmax)))
    
    
    # User defined periodicity in the topl assignment
    # brickrange defines the domain extents
    topl = StackedBrickTopology(mpicomm, brickrange, periodicity=(true,true,false))

    grid = DiscontinuousSpectralElementGrid(topl,
                                            FloatType = DFloat,
                                            DeviceArray = ArrayType,
                                            polynomialorder = N)
    
    numflux!(x...) = NumericalFluxes.rusanov!(x..., cns_flux!, wavespeed, preflux)
    numbcflux!(x...) = NumericalFluxes.rusanov_boundary_flux!(x..., cns_flux!, bcstate!, wavespeed, preflux)

    # spacedisc = data needed for evaluating the right-hand side function
    spacedisc = DGBalanceLaw(grid = grid,
                             length_state_vector = _nstate,
                             flux! = cns_flux!,
                             numerical_flux! = numflux!,
                             numerical_boundary_flux! = numbcflux!, 
                             number_gradient_states = _ngradstates,
                             states_for_gradient_transform =
                             _states_for_gradient_transform,
                             number_viscous_states = _nviscstates,
                             gradient_transform! = gradient_vars!,
                             viscous_transform! = compute_stresses!,
                             viscous_penalty! = stresses_penalty!,
                             viscous_boundary_penalty! = stresses_boundary_penalty!,
                             auxiliary_state_length = _nauxstate,
                             auxiliary_state_initialization! = (x...) ->
                             auxiliary_state_initialization!(x...),
                             source! = source!,
                             preodefun! = integral_computation)

     # This is a actual state/function that lives on the grid
    # ----------------------------------------------------
    # GET DATA FROM INTERPOLATED ARRAY ONTO VECTORS
    # This driver accepts data in 6 column format
    # ----------------------------------------------------
    (sounding, _, ncols) = read_sounding()

    # WARNING: Not all sounding data is formatted/scaled
    # the same. Care required in assigning array values
    # height theta qv    u     v     pressure
    #zinit, tinit, qinit, uinit, vinit, pinit  =
    #    sounding[:, 1], sounding[:, 2], sounding[:, 3], sounding[:, 4], sounding[:, 5], sounding[:, 6]

    zinit, tinit, pinit = sounding[:, 1], sounding[:, 2], sounding[:, 3]
    thetainit, qinit = sounding[:, 4], sounding[:, 5]
    
    #------------------------------------------------------
    # GET SPLINE FUNCTION
    #------------------------------------------------------
    spl_tinit    = Spline1D(zinit, tinit; k=1) #sensible T (K)
    spl_pinit    = Spline1D(zinit, pinit; k=1) #pressure P (Pa)
    
    spl_thetainit= Spline1D(zinit, thetainit; k=1) #sensible T (K)
    spl_qinit    = Spline1D(zinit, qinit; k=1) #sensible T (K)
    
    #initialcondition(Q, x...) = dycoms!(Val(dim), Q, DFloat(0), spl_tinit,
    #                                    spl_qinit, spl_uinit, spl_vinit,
    #                                    spl_pinit, x...)

    initialcondition(Q, x...) = dycoms!(Val(dim), Q, DFloat(0), spl_tinit, spl_pinit, spl_thetainit, spl_qinit, x...)
    
    Q = MPIStateArray(spacedisc, initialcondition)     
    # This is a actual state/function that lives on the grid
    #initialcondition(Q, x...) = dycoms!(Val(dim), Q, DFloat(0), x...)
    #Q = MPIStateArray(spacedisc, initialcondition)
    
    lsrk = LSRK54CarpenterKennedy(spacedisc, Q; dt = dt, t0 = 0)

    #=eng0 = norm(Q)
    @info @sprintf """Starting
      norm(Q₀) = %.16e""" eng0
    =#
    # Set up the information callback
    starttime = Ref(now())
    cbinfo = GenericCallbacks.EveryXWallTimeSeconds(10, mpicomm) do (s=false)
        if s
            starttime[] = now()
        else
            #energy = norm(Q)
            #globmean = global_mean(Q, _ρ)
            @info @sprintf("""Update
                         simtime = %.16e
                         runtime = %s""",
                           ODESolvers.gettime(lsrk),
                           Dates.format(convert(Dates.DateTime,
                                                Dates.now()-starttime[]),
                                        Dates.dateformat"HH:MM:SS")) #, energy )#, globmean)
        end
    end

    npoststates = 9
    _LWP_out, _P, _u, _v, _w, _ρinv, _q_liq, _T, _θ = 1:npoststates
    postnames = ("LWP", "P", "u", "v", "w", "ρinv", "_q_liq", "T", "THETA")
    postprocessarray = MPIStateArray(spacedisc; nstate=npoststates)

    step = [0]
    mkpath("./CLIMA-output-scratch/dycoms-today-Ri/")
    cbvtk = GenericCallbacks.EveryXSimulationSteps(2000) do (init=false)
        DGBalanceLawDiscretizations.dof_iteration!(postprocessarray, spacedisc,
                                                   Q) do R, Q, QV, aux
                                                       @inbounds let
                                                           F_rad_out = radiation(aux)
                                                           (R[_LWP_out], R[_P], R[_u], R[_v], R[_w], R[_ρinv], R[_q_liq], R[_T], R[_θ]) = ( aux[_a_LWP_02z] + aux[_a_LWP_z2inf], preflux(Q, QV, aux)...)
                                                       end
                                                   end

        outprefix = @sprintf("./CLIMA-output-scratch/dycoms-today-Ri/dyc_%dD_mpirank%04d_step%04d", dim,
                             MPI.Comm_rank(mpicomm), step[1])
        @debug "doing VTK output" outprefix
        writevtk(outprefix, Q, spacedisc, statenames,
                 postprocessarray, postnames)
        
        step[1] += 1
        nothing
    end

    # Initialise the integration computation. Kernels calculate this at every timestep?? 
    integral_computation(spacedisc, Q, 0) 
    solve!(Q, lsrk; timeend=timeend, callbacks=(cbinfo, cbvtk))



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
    numelem = (Nex,Ney,Nez)
    dt = 0.005
    timeend = 14400
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
        @info @sprintf """ Dycoms                                                """
        @info @sprintf """   Resolution:                                         """ 
        @info @sprintf """     (Δx, Δy, Δz)   = (%.2e, %.2e, %.2e)               """ Δx Δy Δz
        @info @sprintf """     (Nex, Ney, Nez) = (%d, %d, %d)                    """ Nex Ney Nez
        @info @sprintf """     DoF = %d                                          """ DoF
        @info @sprintf """     Minimum necessary memory to run this test: %g GBs """ (DoFstorage * sizeof(DFloat))/1000^3
        @info @sprintf """     Time step dt: %.2e                                """ dt
        @info @sprintf """     End time  t : %d                                  """ timeend
        @info @sprintf """ ------------------------------------------------------"""
    end
    
    engf_eng0 = run(mpicomm, dim, numelem[1:dim], polynomialorder, timeend,
                    DFloat, dt)
end

isinteractive() || MPI.Finalize()

nothing
