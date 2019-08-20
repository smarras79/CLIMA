# Load Modules 
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
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates
using CLIMA.Vtk
using DelimitedFiles
using Dierckx
using Random
using CLIMA.RootSolvers

if haspkg("CuArrays")
    using CUDAdrv
    using CUDAnative
    using CuArrays
    CuArrays.allowscalar(false)
    const ArrayType = CuArray
else
    const ArrayType = Array
end

using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters

# State labels 
const _nstate = 6
const _ρ, _U, _V, _W, _E, _QT = 1:_nstate
const stateid = (ρid = _ρ, Uid = _U, Vid = _V, Wid = _W, Eid = _E, QTid = _QT)
const statenames = ("RHO", "U", "V", "W", "E", "QT")

# Viscous state labels
const _nviscstates = 30
const _τ11, _τ22, _τ33, _τ12, _τ13, _τ23, _qtx, _qty, _qtz, _JplusDx, _JplusDy, _JplusDz, _θx, _θy, _θz, _uz, _vz, _SijSij, _ν_e, _qvx, _qvy, _qvz, _qlx, _qly, _qlz, _ν_smago, _ν_smago_fR, _ν_vreman, _μ_e, _f_R = 1:_nviscstates

const _nauxstate = 28
const _a_x, _a_y, _a_z, _a_sponge, _a_02z, _a_z2inf, _a_rad, _a_LWP_02z, _a_LWP_z2inf,_a_q_liq, _a_θ, _a_θ_l, _a_P,_a_T, _a_soundspeed_air, _a_y_FN, _a_z_FN, _a_ρ_FN, _a_U_FN, _a_V_FN, _a_W_FN, _a_E_FN, _a_QT_FN, _a_Rm, _a_f_R, _a_ν_smago, _a_ν_smago_fR, _a_ν_vreman = 1:_nauxstate

if !@isdefined integration_testing
    const integration_testing =
        parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
    using Random
end

const μ_sgs           = 100.0
const Prandtl         = 71 // 100
const Prandtl_t       = 1 // 3
const cp_over_prandtl = cp_d / Prandtl_t
#const seed = MersenneTwister(0)
# User Input
const numdims = 2
const Npoly = 4

# Define grid size 
Δx = 20
Δy = 10
Δz = 10

const h_first_layer = Δy

#B.C.
const bc_fix_T   = 0
const bc_no_slip = 0
const bc_fix_bott_flux = 1

# OR:
# Set Δx < 0 and define  Nex, Ney, Nez:
(Nex, Ney) = (2, 1)

# Physical domain extents 
const (xmin, xmax) = (0, 2000) #820)
const (ymin, ymax) = (0, 1500) #820)
const (zmin, zmax) = (0, 1500) #820)

const zi = 840
#Get Nex, Ney from resolution
const Lx = xmax - xmin
const Ly = ymax - ymin
const Lz = zmax - zmin
if ( Δx > 0)
    # User defines the grid size:
    ratiox = (Lx/Δx - 1)/Npoly
    ratioy = (Ly/Δy - 1)/Npoly
    ratioz = (Lz/Δz - 1)/Npoly
    Nex = ceil(Int64, ratiox)
    Ney = ceil(Int64, ratioy)
    Nez = ceil(Int64, ratioz)
else
    # User defines the number of elements:
    Δx = Lx / ((Nex * Npoly) + 1)
    Δy = Ly / ((Ney * Npoly) + 1)
    Δz = Lz / ((Nez * Npoly) + 1)
end

DoF = (Nex*Ney)*(Npoly+1)^numdims*(_nstate)
DoFstorage = (Nex*Ney)*(Npoly+1)^numdims*(_nstate + _nviscstates + _nauxstate + CLIMA.Mesh.Grids._nvgeo) +
    (Nex*Ney)*(Npoly+1)^(numdims-1)*2^numdims*(CLIMA.Mesh.Grids._nsgeo)

const C_smag = 0.15
# Equivalent grid-scale
Δ = (Δx * Δy * Δz)^(1/3)
const Δsqr = Δ * Δ


# Surface values to calculate surface fluxes:
const SST        = 292.5
const psfc       = 1017.8e2      # Pa
const qtot_sfc   = 13.84e-3      # qs(sst) using Teten's formula
const ρsfc       = 1.22          #kg/m^3
const Cd         = 0.0011        #Drag coefficient
const first_node_level   = 0.0001

const D_subsidence = 3.75e-6

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

function strainrate_tensor_components(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    # Assemble components of the strain-rate tensor 
    S11  = dudx
    S12  = (dudy + dvdx) / 2
    S13  = (dudz + dwdx) / 2
    S22  = dvdy
    S23  = (dvdz + dwdy) / 2
    S33  = dwdz
    normSij = S11^2 + S22^2 + S33^2 + 2 * (S12^2 + S13^2 + S23^2)  
    return (S11, S22, S33, S12, S13, S23, normSij)
  end


# -------------------------------------------------------------------------
# Diagnostics: e.g. thermodynamics properties, preflux no longer used in 
# `bcstate`
# -------------------------------------------------------------------------
@inline function preflux(Q,aux)
    R_gas::eltype(Q) = R_d
    @inbounds ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
    ρinv = 1 / ρ
    x,y,z = aux[_a_x], aux[_a_y], aux[_a_z]
    xvert = aux[_a_y]
    u, v, w = ρinv * U, ρinv * V, ρinv * W
    e_int = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * grav * xvert) / ρ
    q_tot = QT / ρ
    # Establish the current thermodynamic state using the prognostic variables
    TS = PhaseEquil(e_int, q_tot, ρ)
    T = air_temperature(TS)
    Rm = gas_constant_air(TS)
    P = air_pressure(TS) # Test with dry atmosphere
    q_liq = PhasePartition(TS).liq
    θ = virtual_pottemp(TS)
    (u, v, w, T, θ, Rm, P)
end

#-------------------------------------------------------------------------
#md # Soundspeed computed using the thermodynamic state TS
# max eigenvalue
@inline function wavespeed(n, Q, aux, t)
    @inbounds begin 
        P = aux[_a_P]
        T = aux[_a_T]
        θ = aux[_a_θ]
        (u, v, w, _, _, _, _) = preflux(Q,aux)
        ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
        ρinv = 1 / ρ
        x,y,z = aux[_a_x], aux[_a_y], aux[_a_z]
        xvert = aux[_a_y]
        u, v, w = ρinv * U, ρinv * V, ρinv * W
        e_int = E/ρ - (u^2 + v^2+ w^2)/2 - grav * xvert
        q_tot = QT / ρ
        TS = PhaseEquil(e_int, q_tot, ρ)
        abs(n[1] * u + n[2] * v + n[3] * w) + soundspeed_air(TS)
    end
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
# ------------------------------------------------------------------------
function buoyancy_correction(modSij, θv, dθvdz)
    # Brunt-Vaisala frequency
    N2 = grav * dθvdz / θv 
    # Richardson number
    Richardson = N2 / (modSij^2 + 1e-12)
    # Buoyancy correction factor
    buoyancy_factor = N2 <=0 ? 1 : sqrt(max(0.0, 1 - Richardson/Prandtl_t))
    return buoyancy_factor
end


@inline function cns_flux!(F, Q, VF, aux, t)
    @inbounds begin
        P = aux[_a_P]
        T = aux[_a_T]
        θ = aux[_a_θ]
        q_liq = aux[_a_q_liq]
        ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
        u,v,w = U/ρ, V/ρ, W/ρ
        xvert = aux[_a_y]
        θ     = aux[_a_θ]
        w -= D_subsidence * xvert
        W  = v*ρ
        # Inviscid contributions
        F[1, _ρ], F[2, _ρ], F[3, _ρ] = U          , V          , W
        F[1, _U], F[2, _U], F[3, _U] = u * U  + P , v * U      , w * U
        F[1, _V], F[2, _V], F[3, _V] = u * V      , v * V + P  , w * V
        F[1, _W], F[2, _W], F[3, _W] = u * W      , v * W      , w * W + P
        F[1, _E], F[2, _E], F[3, _E] = u * (E + P), v * (E + P), w * (E + P)
        F[1, _QT], F[2, _QT], F[3, _QT] = u * QT  , v * QT     , w * QT 
        #Derivative of T and Q:
        vqx, vqy, vqz                     = VF[_qtx],  VF[_qty],  VF[_qtz]
        vqvx, vqvy, vqvz                  = VF[_qvx], VF[_qvy], VF[_qvz]
        vqlx, vqly, vqlz                  = VF[_qlx], VF[_qly], VF[_qlz]    
        vJplusDx, vJplusDy, vJplusDz      = VF[_JplusDx], VF[_JplusDy], VF[_JplusDz]
        vθz                               = VF[_θz]
        
        # Radiation contribution
        F_rad                     = ρ * radiation(aux)  
        aux[_a_rad]               = F_rad

        #Buoyancy correction
        modSij                    = VF[_SijSij]
        dudz                      = VF[_uz]
        dvdz                      = VF[_vz]
        f_R                      = buoyancy_correction(modSij, θ, vθz)
        
        if (xvert > 820)
            f_R = 0.0
        end
        
        #Dynamic eddy viscosity from Smagorinsky:
        ν_e                       = VF[_ν_smago]
        μ_e                       = ρ * ν_e * f_R
        D_e                       = μ_e / Prandtl_t

        #Dynamic eddy viscosity from Vreman: 
        #ν_vreman                  = VF[_ν_vreman]
        #μ_e                       = ν_vreman * ρ * f_R
        #D_e                       = 3ν_vreman
        
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
        F[1, _E] += u * τ11 + v * τ12 + w * τ13 + vJplusDx * D_e  #dTd should not be diffused.
        F[2, _E] += u * τ21 + v * τ22 + w * τ23 + vJplusDy * D_e
        F[3, _E] += u * τ31 + v * τ32 + w * τ33 + vJplusDz * D_e
        F[3, _E] += F_rad
        # Viscous contributions to mass flux terms
        F[1, _ρ]  +=  vqx * D_e
        F[2, _ρ]  +=  vqy * D_e
        F[3, _ρ]  +=  vqz * D_e
        F[1, _QT] +=  vqx * D_e
        F[2, _QT] +=  vqy * D_e
        F[3, _QT] +=  vqz * D_e
    end
end

# -------------------------------------------------------------------------
#md # Here we define a function to extract the velocity components from the 
#md # prognostic equations (i.e. the momentum and density variables). This 
#md # function is not required in general, but provides useful functionality 
#md # in some cases. 
# -------------------------------------------------------------------------
# Compute the velocity from the state
# Gradient state labels
const _ngradstates = 7
@inline function gradient_vars!(grad_vars, Q, aux, t)
    @inbounds begin
        Rm = aux[_a_Rm]
        T = aux[_a_T]
        ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
        grad_vars[1], grad_vars[2], grad_vars[3] = U/ρ, V/ρ, W/ρ
        grad_vars[4], grad_vars[5], grad_vars[6] = E, QT/ρ, (E/ρ + Rm*T)
        grad_vars[7] = aux[_a_θ]
    end
end

@inline function radiation(aux)
    zero_to_z = aux[_a_02z]
    z_to_inf = aux[_a_z2inf]
    xvert = aux[_a_y]
    z_i = 840  # Start with constant inversion height of 840 meters then build in check based on q_tot
    (xvert - z_i) >=0 ? Δz_i = (xvert - z_i) : Δz_i = 0 
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
@inline function compute_stresses!(VF, grad_list, _...)
    @inbounds begin
        # get gradients of velocities 
        dudx, dudy, dudz = grad_list[1, 1], grad_list[2, 1], grad_list[3, 1]
        dvdx, dvdy, dvdz = grad_list[1, 2], grad_list[2, 2], grad_list[3, 2]
        dwdx, dwdy, dwdz = grad_list[1, 3], grad_list[2, 3], grad_list[3, 3]
        # get gradients of moist vars and potential temperature
        dqtdx, dqtdy, dqtdz             = grad_list[1, 5], grad_list[2, 5], grad_list[3, 5]
        dJplusDdx, dJplusDdy, dJplusDdz = grad_list[1, 6], grad_list[2, 6], grad_list[3, 6]
        dθdx, dθdy, dθdz                = grad_list[1, 7], grad_list[2, 7], grad_list[3, 7]

        # --------------------------------------------
        # STRAINRATE TENSOR COMPONENTS
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
        SijSij = (S11^2 + S22^2 + S33^2
                  + 2.0 * S12^2
                  + 2.0 * S13^2 
                  + 2.0 * S23^2)
        modSij = sqrt(2.0 * SijSij)
        
        ν_smago = modSij * C_smag^2 * Δsqr
        
        # -------------------------------------------
        # VREMAN COEFFICIENT
        # -------------------------------------------
        Dij = SMatrix{3,3}(dudx, dudy, dudz, 
                           dvdx, dvdy, dvdz,
                           dwdx, dwdy, dwdz)
        DF = eltype(SijSij)
        DISS = sum(Dij .^ 2)
        βij = similar(Dij)
        βij = Δsqr * (Dij' * Dij)
        Bβ = βij[1,1]*βij[2,2] - βij[1,2]^2 + βij[1,1]*βij[3,3] - βij[1,3]^2 + βij[2,2]*βij[3,3] - βij[2,3]^2 
        ν_vreman = max(0,(C_smag^2 * 2.5) * sqrt(abs(Bβ/(DISS+1e-16)))) 
        
        #--------------------------------------------
        # STRESS COMPONENTS
        #--------------------------------------------
        VF[_τ11] = -2 * (S11 - (S11 + S22 + S33) / 3)
        VF[_τ22] = -2 * (S22 - (S11 + S22 + S33) / 3)
        VF[_τ33] = -2 * (S33 - (S11 + S22 + S33) / 3)
        VF[_τ12] = -2 * S12
        VF[_τ13] = -2 * S13
        VF[_τ23] = -2 * S23
        
        VF[_qtx], VF[_qty], VF[_qtz]             = -dqtdx,  -dqtdy,  -dqtdz        
        VF[_JplusDx], VF[_JplusDy], VF[_JplusDz] = -dJplusDdx, -dJplusDdy, -dJplusDdz
        VF[_θx], VF[_θy], VF[_θz]                = dθdx, dθdy, dθdz
        VF[_SijSij]                              = modSij
        VF[_uz]                                  = dudz
        VF[_vz]                                  = dvdz
                
        VF[_ν_vreman]                            = ν_vreman
        VF[_ν_smago]                             = ν_smago

    end
end
# -------------------------------------------------------------------------
#md ### Auxiliary Function 
# -------------------------------------------------------------------------
@inline function auxiliary_state_initialization!(aux, x, y, z)
    @inbounds begin
        DFloat = eltype(aux)
        xvert = y
        aux[_a_y] = xvert
        #Sponge 
        ctop    = zero(DFloat)
        cs_left_right = zero(DFloat)
        cs_front_back = zero(DFloat)
        ct            = DFloat(0.75)
        domain_bott  = ymin
        domain_top   = ymax
        
        #Vertical sponge:
        sponge_type = 1
        if sponge_type == 1
            ct = 0.9
            bc_zscale  = 450.0
            zd = domain_top - bc_zscale       
            if xvert >= zd
                ctop = ct * (sinpi(0.5*(xvert - zd)/(domain_top - zd)))^4
            end
        else
            aux[_a_x] = x
            aux[_a_y] = y
            #Sponge
            csleft  = 0.0
            csright = 0.0
            csfront = 0.0
            csback  = 0.0
            ctop    = 0.0
            cs_left_right = 0.0
            cs_front_back = 0.0
            ct            = 0.75
            #BEGIN  User modification on domain parameters.
            domain_left  = xmin 
            domain_right = xmax
            domain_front = ymin 
            domain_back  = ymax 
            domain_bott  = zmin 
            domain_top   = zmax 
            #END User modification on domain parameters.
            # Define Sponge Boundaries      
            xc       = 0.5 * (domain_right + domain_left)
            yc       = 0.5 * (domain_back  + domain_front)
            zc       = 0.5 * (domain_top   + domain_bott)
            top_sponge  = 0.85 * domain_top
            xsponger    = domain_right - 0.15 * (domain_right - xc)
            xspongel    = domain_left  + 0.15 * (xc - domain_left)
            ysponger    = domain_back  - 0.15 * (domain_back - yc)
            yspongel    = domain_front + 0.15 * (yc - domain_front)
            #x left and right
            #xsl
            if x <= xspongel
                csleft = cs_left_right * (sinpi(1/2 * (x - xspongel)/(domain_left - xspongel)))^4
            end
            #xsr
            if x >= xsponger
                csright = cs_left_right * (sinpi(1/2 * (x - xsponger)/(domain_right - xsponger)))^4
            end        
            #y left and right
            #ysl
            if y <= yspongel
                csfront = cs_front_back * (sinpi(1/2 * (y - yspongel)/(domain_front - yspongel)))^4
            end
            #ysr
            if y >= ysponger
                csback = cs_front_back * (sinpi(1/2 * (y - ysponger)/(domain_back - ysponger)))^4
            end
            #Vertical sponge:         
            if z >= top_sponge
                ctop = ct * (sinpi(0.5 * (z - top_sponge)/(domain_top - top_sponge)))^4
            end
        end
        beta  = 1.0 - (1.0 - ctop) #*(1.0 - csleft)*(1.0 - csright)*(1.0 - csfront)*(1.0 - csback)
        beta  = min(beta, 1.0)
        aux[_a_sponge] = beta
    end
end

# -------------------------------------------------------------------------
@inline function bcstate!(QP, VFP, auxP, nM, QM, VFM, auxM, bctype, t)
    @inbounds begin
        ρM, UM, VM, WM, EM, QTM = QM[_ρ], QM[_U], QM[_V], QM[_W], QM[_E], QM[_QT]
        uM, vM, wM  = UM/ρM, VM/ρM, WM/ρM
        q_totM = QTM/ρM
        q_liqM = auxM[_a_q_liq]
        UnM = nM[1] * UM + nM[2] * VM + nM[3] * WM
        QP[_U] = UM - 2 * nM[1] * UnM
        QP[_V] = VM - 2 * nM[2] * UnM
        QP[_W] = WM - 2 * nM[3] * UnM
        QP[_ρ] = ρM
        QP[_QT] = QTM
        VFP .= 0
        xvert = auxM[_a_y]
        if xvert < 0.00001
            if bc_fix_bott_flux == 1
                # TODO specify boundary keyword and get correct bctype for general topography
                # ------------------------------------------------------------------------
                # First node quantities (first-model level here represents the first node)
                # ------------------------------------------------------------------------
                z_FN             = auxM[_a_y_FN]
                ρ_FN             = auxM[_a_ρ_FN]
                U_FN             = auxM[_a_U_FN]
                V_FN             = auxM[_a_V_FN]
                W_FN             = auxM[_a_W_FN]
                E_FN             = auxM[_a_E_FN]
                u_FN, v_FN, w_FN = U_FN/ρ_FN, V_FN/ρ_FN, W_FN/ρ_FN
                windspeed_FN     = sqrt(u_FN^2 + v_FN^2 + w_FN^2)
                q_tot_FN         = auxM[_a_QT_FN] / ρ_FN
                e_int_FN         = E_FN/ρ_FN - 0.5*windspeed_FN^2 - grav*z_FN
                TS_FN            = PhaseEquil(e_int_FN, q_tot_FN, ρ_FN) 
                T_FN             = air_temperature(TS_FN)
                PhPart           = PhasePartition(TS_FN)
                q_liq_FN         = PhasePartition(TS_FN).liq
                q_vap_FN         = q_tot_FN - q_liq_FN
                
                # -----------------------------------
                # Bottom boundary quantities 
                # -----------------------------------
                zM          = xvert
                q_totM      = QTM/ρM
                q_liqM      = auxM[_a_q_liq]
                windspeed   = sqrt(uM^2 + vM^2)
                e_intM      = EM/ρM - 0.5*windspeed^2 - grav*zM
                TSM         = PhaseEquil(e_intM, q_totM, ρM)
                TM          = air_temperature(TSM)
                
                #q_v_sfc     = q_totM - PhasePartition(TSM).liq
                q_v_sfc     = q_vap_saturation(SST, ρM, PhPart)   #Dedf. between eq. (59) and (eq 60) in the CLIMA-doc
                
                
                # ----------------------------------------------
                # Assigning calculated values to boundary states
                # ----------------------------------------------
                VFP[_τ33] = 0  
                # Case specific for flat bottom topography, normal vector is n⃗ = k⃗ = [0, 0, 1]ᵀ
                
                # A more general implementation requires (n⃗ ⋅ ∇A) to be defined where A is replaced by the appropriate flux terms
                VFP[_τ13]     = -ρM * Cd * windspeed_FN * u_FN 
                VFP[_τ23]     = -ρM * Cd * windspeed_FN * v_FN 
                VFP[_qtz]     = +115 /(ρM * LH_v0)
                #VFP[_qtz]     = -ρM * Cd * windspeed_FN * (q_vap_FN - q_v_sfc)
                VFP[_JplusDz] = +130 / ρM
                #=
                h_FN  = internal_energy(T_FN, PhPart) + Rm * T_FN
                h_sfc = internal_energy(SST, PhPart) + Rm * SST
                Φ_FN = grav * z_FN
                Φ_sfc = grav * xvert
                VFP[_JplusDz] = -ρM * Cd * windspeed_FN * (h_FN - h_sfc + Φ_FN - Φ_sfc)
                =#
                
            end #bc_fix_bott_flux

            if bc_fix_T == 1
                UnM = nM[1] * UM + nM[2] * VM + nM[3] * WM
                QP[_W] = WM - 2 * nM[3] * UnM
                QP[_U] = QM[_U]
                QP[_V] = QM[_V]
                u, v, w = QP[_U]/ρM, QP[_V]ρM, QP[_W]/ρM
                
                if bc_no_slip == 1
                    QP[_U], QP[_V], QP[_W] = 0.0, 0.0, 0.0
                    u, v, w = QP[_U]/ρM, QP[_V]ρM, QP[_W]/ρM
                end
                #VFP .= VFM
                #VFP[_Tz] = VFM[_Tz]
                
                #Dirichlet on T: SST
                T     = SST
                ρP    = ρM
                e_kin = 0.5 * (u^2 + v^2)
                e_pot = grav * xvert
                e_int = internal_energy(T, PhasePartition(q_totM, q_liqM, 0.0))
                E     = ρM * total_energy(e_kin, e_pot, T, PhasePartition(q_totM, q_liqM, 0.0))
                QP[_E]  = E
                #QP[_QT] = qtot_sfc
                QP[_QT] = QM[_QT]
            end
            
        end
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
        W = Q[_W]
        U = Q[_U]
        V = Q[_V]
        S[_U] -= f_coriolis * (U/ρ - u_geostrophic)
        S[_V] -= f_coriolis * (V/ρ - v_geostrophic)
    end
end

"""
    Sponge source (coefficient from aux)
    """
@inline function source_sponge!(S,Q,aux,t)
    @inbounds begin
        U, V   = Q[_U], Q[_V]       
        beta   = aux[_a_sponge]
        S[_U] -= beta * U
        S[_V] -= beta * V
    end
end

"""
    Gravity source 
    """
@inline function source_geopot!(S,Q,aux,t)
    @inbounds begin
        ρ, U, V, E  = Q[_ρ], Q[_U], Q[_V], Q[_E]
        S[_V] -= ρ * grav
    end
end

"""
    Defines integrand for DYCOMS radiation
    """
@inline function integral_knl(val, Q, aux)
    κ = 85.0
    @inbounds begin
        @inbounds ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
        ρinv = 1 / ρ
        x,y,z = aux[_a_x], aux[_a_y], aux[_a_z]        
        xvert = aux[_a_y]
        u, v, w = ρinv * U, ρinv * V, ρinv * W
        e_int = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * grav * xvert) / ρ
        q_tot = QT / ρ
        # Establish the current thermodynamic state using the prognostic variables
        TS     = PhaseEquil(e_int, q_tot, ρ)
        q_liq  = PhasePartition(TS).liq
        val[1] = ρ * κ * q_liq 
        val[2] = ρ * q_liq       # LWP Integrand
    end
end

"""
    Stores thermodynamic properties and computes other quantities required prior to time update
    """
function preodefun!(disc, Q, t)
    DGBalanceLawDiscretizations.dof_iteration!(disc.auxstate, disc, Q) do R, Q, QV, aux
        @inbounds let
            ρ, U, V, W, E, QT = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]
            xvert = aux[_a_y]
            e_int = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * grav * xvert) / ρ
            q_tot                = QT / ρ
            TS                   = PhaseEquil(e_int, q_tot, ρ)
            T                    = air_temperature(TS)
            P                    = air_pressure(TS) # Test with dry atmosphere
            q_liq                = PhasePartition(TS).liq
            PhPart               = PhasePartition(TS)
            R[_a_T]              = T
            R[_a_P]              = P
            R[_a_q_liq]          = q_liq
            R[_a_soundspeed_air] = soundspeed_air(TS)
            R[_a_θ]              = virtual_pottemp(TS)
            R[_a_θ_l]            = liquid_ice_pottemp(TS)
            R[_a_Rm]             = gas_constant_air(TS)

            #Viscosity:
            modSij               = QV[_SijSij]

            dudz                 = QV[_uz]
            dvdz                 = QV[_vz]
            vθz                  = QV[_θz]           
            θ                    = aux[_a_θ]
            f_R                  = buoyancy_correction(modSij, θ, vθz)
            
            R[_a_ν_smago]        = QV[_ν_smago]*f_R           
            R[_a_ν_vreman]       = QV[_ν_vreman]
            
        end
    end
    firstnode_info(disc,Q,t)
    integral_computation(disc, Q, t)
end

function firstnode_info(disc,Q,t)
    # User specified kernel to allow access to first-interior points of specific auxiliary and state variables
    DGBalanceLawDiscretizations.aux_firstnode_values!(disc, Q,
                                                      (_a_y_FN), (_a_y))
    DGBalanceLawDiscretizations.state_firstnode_values!(disc, Q,
                                                        (_a_ρ_FN, _a_U_FN, _a_V_FN, _a_W_FN, _a_E_FN, _a_QT_FN), (_ρ, _U, _V, _W, _E, _QT))
end
function integral_computation(disc, Q, t)
    # Kernel to compute vertical integrals
    DGBalanceLawDiscretizations.indefinite_stack_integral!(disc, integral_knl, Q,
                                                           (_a_02z, _a_LWP_02z))
    
    DGBalanceLawDiscretizations.reverse_indefinite_stack_integral!(disc,
                                                                   _a_z2inf,
                                                                   _a_02z)
    
    DGBalanceLawDiscretizations.reverse_indefinite_stack_integral!(disc,
                                                                   _a_LWP_z2inf,
                                                                   _a_LWP_02z)
end

function dycoms!(dim, Q, t, x, y, z, _...)
    
    DFloat         = eltype(Q)
    xvert::DFloat  = y

    #These constants are those used by Stevens et al. (2005)
    R_d::DFloat     = 287.0
    cp_d::DFloat    = 1015.0
    cp_v::DFloat    = 1859.0
    cp_l::DFloat    = 4181.0
    Lv::DFloat      = 2.47e6
    epsdv::DFloat   = 1.61
    g::DFloat       = grav
    p0::DFloat      = 1.0178e5
    ρ0::DFloat      = 1.22
    r_tot_sfc::DFloat=8.1e-3
    Rm_sfc          = R_d * (1.0 + (epsdv - 1.0)*r_tot_sfc)
    ρ_sfc::DFloat   = 1.22
    P_sfc           = 1.0178e5
    T_0::DFloat     = 285.0
    T_sfc           = P_sfc/(ρ_sfc * Rm_sfc);
    
    # --------------------------------------------------
    # INITIALISE ARRAYS FOR INTERPOLATED VALUES
    # --------------------------------------------------
    randnum1   = rand(1)[1] / 100
    randnum2   = rand(1)[1] / 100

    
    q_liq      = 0.0
    q_ice      = 0.0
    zb         = 600.0    #initial cloud bottom
    zi         = 840.0    #initial cloud top
    dz_cloud   = zi - zb
    q_liq_peak = 0.00045 #cloud mixing ratio at z_i    
    if xvert > zb && xvert <= zi	
	q_liq = (xvert - zb)*q_liq_peak/dz_cloud
    end

    if ( xvert <= zi)
	θ_liq  = 289.0
	r_tot      = 8.1e-3                  #kg/kg  specific humidity --> approx. to mixing ratio is ok
	q_tot      = r_tot #/(1.0 - r_tot)     #total water mixing ratio
    else
	θ_liq = 297.5 + (xvert - zi)^(1/3)
	r_tot     = 1.5e-3                    #kg/kg  specific humidity --> approx. to mixing ratio is ok
	q_tot     = r_tot #/(1.0 - r_tot)      #total water mixing ratio
    end

    if xvert <= 200.0
        θ_liq += randnum1 * θ_liq 
        q_tot += randnum2 * q_tot
    end
    
    Rm       = R_d * (1 + (epsdv - 1)*q_tot - epsdv*q_liq);
    cpm     = cp_d + (cp_v - cp_d)*q_tot + (cp_l - cp_v)*q_liq;

    #Pressure
    H = Rm_sfc * T_0 / g;
    P = P_sfc * exp(-xvert/H);
    
    #Exner
    exner = (P/P_sfc)^(R_d/cp_d);
    
    #T, Tv 
    T     = exner*θ_liq + Lv*q_liq/(cpm*exner);
    Tv    = T*(1 + (epsdv - 1)*q_tot - epsdv*q_liq);
    
    #Density
    ρ  = P/(Rm*T);
    
    #θ, θv
    θ      = T/exner;
    θv     = θ*(1 + (epsdv - 1)*q_tot - epsdv*q_liq);
    PhPart = PhasePartition(q_tot, q_liq, q_ice)
    
    # energy definitions   
    u, v, w     = 7, 0.0, 0.0 #geostrophic. TO BE BUILT PROPERLY if Coriolis is considered
    U           = ρ * u
    V           = ρ * v
    W           = ρ * w
    e_kin       = 0.5 * (u^2 + v^2)
    e_pot       = grav * xvert
    E           = ρ * total_energy(e_kin, e_pot, T, PhPart)

    @inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]= ρ, U, V, W, E, ρ * q_tot
    #@inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]=  T, U, V, W, P, q_liq   #for initial state plottin only:
    
end


# ------------------------------------------------------------------
# -------------END DEF SOURCES-------------------------------------# 

function run(mpicomm, dim, Ne, N, timeend, DFloat, dt)

    brickrange = (range(DFloat(xmin), length=Ne[1]+1, DFloat(xmax)),
                  range(DFloat(ymin), length=Ne[2]+1, DFloat(ymax)))
    
    
    # User defined periodicity in the topl assignment
    # brickrange defines the domain extents
    topl = StackedBrickTopology(mpicomm, brickrange, periodicity=(true,false))

    grid = DiscontinuousSpectralElementGrid(topl,
                                            FloatType = DFloat,
                                            DeviceArray = ArrayType,
                                            polynomialorder = N)
    
    numflux!(x...) = NumericalFluxes.rusanov!(x..., cns_flux!, wavespeed)
    numbcflux!(x...) = NumericalFluxes.rusanov_boundary_flux!(x..., cns_flux!, bcstate!, wavespeed)

    # spacedisc = data needed for evaluating the right-hand side function
    spacedisc = DGBalanceLaw(grid = grid,
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
    
 
    initialcondition(Q, x...) = dycoms!(Val(dim), Q, 0, x...)

    Q = MPIStateArray(spacedisc, initialcondition)     
    
    lsrk = LSRK54CarpenterKennedy(spacedisc, Q; dt = dt, t0 = 0)

    # Set up the information callback
    starttime = Ref(now())
    cbinfo = GenericCallbacks.EveryXWallTimeSeconds(10, mpicomm) do (s=false)
        if s
            starttime[] = now()
        else
            ql_max = global_max(spacedisc.auxstate, _a_q_liq)
            @info @sprintf("""Update
                             simtime = %.16e
                             runtime = %s
                             max(ql) = %.16e""",
                           ODESolvers.gettime(lsrk),
                           Dates.format(convert(Dates.DateTime,
                                                Dates.now()-starttime[]),
                                        Dates.dateformat"HH:MM:SS"), ql_max)
        end
    end

    npoststates = 6
    _o_LWP, _o_q_liq, _o_θ_l, _o_P, _o_smago, _o_vreman = 1:npoststates
    postnames = ("LWP", "q_liq", "THETA_L", "P", "MU_SMAGO", "MU_VREMAN")
    postprocessarray = MPIStateArray(spacedisc; nstate=npoststates)

    step = [0]
    cbvtk = GenericCallbacks.EveryXSimulationSteps(1500) do (init=false)
        DGBalanceLawDiscretizations.dof_iteration!(postprocessarray, spacedisc, Q) do R, Q, QV, aux
            @inbounds let
                R[_o_LWP]   = aux[_a_LWP_02z] + aux[_a_LWP_z2inf]
                R[_o_q_liq] = aux[_a_q_liq]
                R[_o_θ_l]   = aux[_a_θ_l]
                R[_o_P]     = aux[_a_P]
                
                #R[_o_f_R]   = aux[_a_f_R]
                R[_o_smago] = aux[_a_ν_smago]
                R[_o_vreman]= aux[_a_ν_vreman] #This is actually smago*f_R
                
            end
        end
        
        mkpath("./CLIMA-output-scratch/dycoms-bc-SMAGO-FR-analytic-TUNED-dx35m-dy35m-dz10m-2D-FLUXBC/")
        outprefix = @sprintf("./CLIMA-output-scratch/dycoms-bc-SMAGO-FR-analytic-TUNED-dx35m-dy35m-dz10m-2D-FLUXBC/dy_%dD_mpirank%04d_step%04d", dim,
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

    numelem = (Nex,Ney)
    dt = 0.001
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
        @info @sprintf """     (xmin, xmax)    = (%.2e, %.2e)                    """ xmin xmax
        @info @sprintf """     (ymin, ymax)    = (%.2e, %.2e)                    """ ymin ymax
        @info @sprintf """     (Δx, Δy)        = (%.2e, %.2e)                    """ Δx Δy
        @info @sprintf """     (Nex, Ney, Nez) = (%d, %d)                        """ Nex Ney
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
