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
using Random

#using CLIMA.SubgridScaleTurbulence
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters: R_d, cp_d, grav, cv_d, MSLP, T_0
using CLIMA.Microphysics

if haspkg("CuArrays")
    using CUDAdrv
    using CUDAnative
    using CuArrays
    CuArrays.allowscalar(false)
    const ArrayType = CuArray
else
    const ArrayType = Array
end

"""
State labels
"""
const _nstate = 8
const _ρ, _U, _V, _W, _E, _QT, _QL, _QR = 1:_nstate
const stateid = (ρid = _ρ, Uid = _U, Vid = _V, Wid = _W, Eid = _E, QTid = _QT, QLid = _QL, QRid = _QR,)
const statenames = ("RHO", "U", "V", "W", "E", "QT", "QL", "QR")


"""
Viscous state labels
"""
const _nviscstates = 16
const _τ11, _τ22, _τ33, _τ12, _τ13, _τ23, _qx, _qy, _qz, _Tx, _Ty, _Tz, _ρx, _ρy, _ρz, _SijSij = 1:_nviscstates

"""
Number of variables of which gradients are required 
"""
const _ngradstates = 8

"""
Number of states being loaded for gradient computation
"""
const _states_for_gradient_transform = (_ρ, _U, _V, _W, _E, _QT, _QL, _QR)


if !@isdefined integration_testing
    const integration_testing =
        parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
    using Random
end

"""
Problem constants 
"""
const Prandtl   = 71 // 100
const k_μ       = cp_d / Prandtl

"""
Problem Description
-------------------
2 Dimensional falling thermal bubble (cold perturbation in a warm neutral atmosphere)
"""



#
# Define grid size
#
const numdims = 3
Δx    =  25
Δy    =  25
Δz    =  25
Npoly = 4

#
# OR:
#
# Set Δx < 0 and define  Nex, Ney, Nez:
#
(Nex, Ney, Nez) = (10, 10, 10)

# Physical domain extents
const (xmin, xmax) = (-5000, 5000)
const (ymin, ymax) = (-5000, 5000)
const (zmin, zmax) = (   0,  5000)


#Get Nex, Ney from resolution
const Lx = xmax - xmin
const Ly = ymax - ymin
const Lz = zmax - zmin

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


# Equivalent grid-scale


# -------------------------------------------------------------------------
#md ### Auxiliary Function (Not required)
#md # In this example the auxiliary function is used to store the spatial
#md # coordinates and the equivalent grid lengthscale coefficient. 
# -------------------------------------------------------------------------
const _nauxstate = 4
const _a_x, _a_y, _a_z, _a_sponge = 1:_nauxstate
@inline function auxiliary_state_initialization!(aux, x, y, z)
    @inbounds begin
        aux[_a_x] = x
        aux[_a_y] = y
        aux[_a_z] = z
        
        #Sponge
        csleft  = 0.0
        csright = 0.0
        csfront = 0.0
        csback  = 0.0
        ctop    = 0.0

        cs_left_right = 0.0
        cs_front_back = 0.0
        ct            = 0.9

        #BEGIN  User modification on domain parameters.
        #Only change the first index of brickrange if your axis are
        #oriented differently:
        #x, y, z = aux[_a_x], aux[_a_y], aux[_a_z]
        #TODO z is the vertical coordinate
        #
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

        sponge_type = 2
        if sponge_type == 1

            bc_zscale   = 7000.0
            top_sponge  = 0.85 * domain_top
            zd          = domain_top - bc_zscale
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

        elseif sponge_type == 2


            alpha_coe = 0.5
            bc_zscale = 7500.0
            zd        = domain_top - bc_zscale
            xsponger  = domain_right - 0.15 * (domain_right - xc)
            xspongel  = domain_left  + 0.15 * (xc - domain_left)
            ysponger  = domain_back  - 0.15 * (domain_back - yc)
            yspongel  = domain_front + 0.15 * (yc - domain_front)

            #
            # top damping
            # first layer: damp lee waves
            #
            ctop = 0.0
            ct   = 0.5
            if z >= zd
                zid = (z - zd)/(domain_top - zd) # normalized coordinate
                if zid >= 0.0 && zid <= 0.5
                    abstaud = alpha_coe*(1.0 - cos(zid*pi))

                else
                    abstaud = alpha_coe*( 1.0 + cos((zid - 0.5)*pi) )

                end
                ctop = ct*abstaud
            end

        end #sponge_type

        beta  = 1.0 - (1.0 - ctop) #*(1.0 - csleft)*(1.0 - csright)*(1.0 - csfront)*(1.0 - csback)
        beta  = min(beta, 1.0)
        aux[_a_sponge] = beta
    end
end

# -------------------------------------------------------------------------
# Preflux calculation: This function computes parameters required for the 
# DG RHS (but not explicitly solved for as a prognostic variable)
# In the case of the rising_thermal_bubble example: the saturation
# adjusted temperature and pressure are such examples. Since we define
# the equation and its arguments here the user is afforded a lot of freedom
# around its behaviour. Future drivers won't use the preflux function.  
# The preflux function interacts with the following  
# Modules: NumericalFluxes.jl 
# functions: wavespeed, cns_flux!, bcstate!
# -------------------------------------------------------------------------
@inline function preflux(Q,VF, aux, _...)
    gravity::eltype(Q) = grav
    R_gas::eltype(Q) = R_d
    @inbounds ρ, U, V, W, E, QT, QL, QR = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT], Q[_QL], Q[_QR]

    DF = eltype(ρ)
    
    ρinv    = 1 / ρ
    x,y,z   = aux[_a_x], aux[_a_y], aux[_a_z]
    u, v, w = ρinv * U, ρinv * V, ρinv * W
    e_int   = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * gravity * z) / ρ
    
    q_tot = QT / ρ
    q_liq = QL / ρ
    q_rai = QR / ρ
    q_ice = DF(0)

    # compute rain fall speed    
    #if(q_rai >= DF(0)) #TODO - need a way to prevent negative values
    #    w_rain = terminal_velocity(q_rai, ρ)
    #else
        w_rain = DF(0)
    #end
    
    # Establish the current thermodynamic state using the prognostic variables
    q     = PhasePartition(q_tot, q_liq, q_ice)
    T     = air_temperature(e_int, q)
    P     = air_pressure(T, ρ, q)

    #TS = PhaseEquil(e_int, q_tot, ρ)
    #T = air_temperature(TS)
    #P = air_pressure(TS) # Test with dry atmosphere
    
    return (P, u, v, w, w_rain, ρinv, q_tot, q_liq, q_rai)
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
    #fsounding  = open(joinpath(@__DIR__, "../soundings/sounding_JCP2013_with_pressure.dat"))
    fsounding  = open(joinpath(@__DIR__, "../soundings/sounding_gabersek.dat"))
    sounding = readdlm(fsounding)
    close(fsounding)
    (nzmax, ncols) = size(sounding)
    if nzmax == 0
        error("SOUNDING ERROR: The Sounding file is empty!")
    end
    return (sounding, nzmax, ncols)
end


#-------------------------------------------------------------------------
#md # Soundspeed computed using the thermodynamic state TS
# max eigenvalue
@inline function wavespeed(n, Q, aux, t, P, u, v, w, w_rain, ρinv, q_tot, q_liq,
                           q_rai)
    
    gravity::eltype(Q) = grav
    DF = eltype(Q)
    
  @inbounds begin 
      ρ, U, V, W, E = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E]
      
      x,y,z   = aux[_a_x], aux[_a_y], aux[_a_z]
      u, v, w = ρinv * U, ρinv * V, ρinv * W
      e_int   = (E - (U^2 + V^2+ W^2)/(2*ρ) - ρ * gravity * z) / ρ
      TS      = PhaseEquil(e_int, q_tot, ρ)
      
      abs(n[1] * u + n[2] * v + n[3] * max(abs(w), abs(w_rain), abs(w - w_rain))) + soundspeed_air(TS)
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
# -------------------------------------------------------------------------
cns_flux!(F, Q, VF, aux, t) = cns_flux!(F, Q, VF, aux, t, preflux(Q,VF, aux)...)
@inline function cns_flux!(F, Q, VF, aux, t, P, u, v, w, w_rain, ρinv, q_tot, q_liq, q_rai)
  gravity::eltype(Q) = grav
  @inbounds begin
    ρ, U, V, W, E, QT, QL, QR = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT], Q[_QL], Q[_QR]

    # Inviscid contributions
    F[1, _ρ], F[2, _ρ], F[3, _ρ] = U          , V          , W
    F[1, _U], F[2, _U], F[3, _U] = u * U  + P , v * U      , w * U
    F[1, _V], F[2, _V], F[3, _V] = u * V      , v * V + P  , w * V
    F[1, _W], F[2, _W], F[3, _W] = u * W      , v * W      , w * W + P
    F[1, _E], F[2, _E], F[3, _E] = u * (E + P), v * (E + P), w * (E + P)
    F[1, _QT], F[2, _QT], F[3, _QT] = u * ρ * q_tot, v * ρ * q_tot, w * ρ * q_tot
    F[1, _QL], F[2, _QL], F[3, _QL] = u * ρ * q_liq, v * ρ * q_liq, w * ρ * q_liq
    F[1, _QR], F[2, _QR], F[3, _QR] = u * ρ * q_rai, v * ρ * q_rai, (w - w_rain) * ρ * q_rai

    #Derivative of T and Q:
    vqx, vqy, vqz = VF[_qx], VF[_qy], VF[_qz]        
    vTx, vTy, vTz = VF[_Tx], VF[_Ty], VF[_Tz]
    vρy = VF[_ρy]
    SijSij = VF[_SijSij]
    
    (ν_e, D_e) = 200, 200 #SubgridScaleTurbulence.standard_smagorinsky(SijSij, Δsqr)
    
    #Richardson contribution:
    f_R = 1 #SubgridScaleTurbulence.buoyancy_correction(SijSij, ρ, vρy)
    
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
    F[1, _E] -= u * τ11 + v * τ12 + w * τ13 + ν_e * k_μ * vTx 
    F[2, _E] -= u * τ21 + v * τ22 + w * τ23 + ν_e * k_μ * vTy
    F[3, _E] -= u * τ31 + v * τ32 + w * τ33 + ν_e * k_μ * vTz 
  end
end

# -------------------------------------------------------------------------
#md # Here we define a function to extract the velocity components from the 
#md # prognostic equations (i.e. the momentum and density variables). This 
#md # function is not required in general, but provides useful functionality 
#md # in some cases. 
# -------------------------------------------------------------------------
# Compute the velocity from the state
gradient_vars!(gradient_list, Q, aux, t, _...) = gradient_vars!(gradient_list, Q, aux, t, preflux(Q,~,aux)...)
@inline function gradient_vars!(gradient_list, Q, aux, t, P, 
u, v, w, w_rain, ρinv, q_tot, q_liq, q_rai)
    @inbounds begin
         z= aux[_a_z]
        # ordering should match states_for_gradient_transform
        ρ, U, V, W, E = Q[_ρ],  Q[_U],  Q[_V], Q[_W], Q[_E]
        QT, QL, QR    = Q[_QT], Q[_QL], Q[_QR]
        ρinv = 1 / ρ
        gradient_list[1], gradient_list[2], gradient_list[3] = u, v, w
        gradient_list[4], gradient_list[5], gradient_list[6], gradient_list[7] = E, QT, QL, QR
        gradient_list[8] = ρ
    end
end

# -------------------------------------------------------------------------
#md ### Viscous fluxes. 
#md # The viscous flux function compute_stresses computes the components of 
#md # the velocity gradient tensor, and the corresponding strain rates to
#md # populate the viscous flux array VF. SijSij is calculated in addition
#md # to facilitate implementation of the constant coefficient Smagorinsky model
#md # (pending)
@inline function compute_stresses!(VF, grad_vars, _...)
    gravity::eltype(VF) = grav
    @inbounds begin
        dudx, dudy, dudz = grad_vars[1, 1], grad_vars[2, 1], grad_vars[3, 1]
        dvdx, dvdy, dvdz = grad_vars[1, 2], grad_vars[2, 2], grad_vars[3, 2]
        dwdx, dwdy, dwdz = grad_vars[1, 3], grad_vars[2, 3], grad_vars[3, 3]
        # compute gradients of moist vars and temperature
        dqdx, dqdy, dqdz = grad_vars[1, 5], grad_vars[2, 5], grad_vars[3, 5]
        dTdx, dTdy, dTdz = grad_vars[1, 6], grad_vars[2, 6], grad_vars[3, 6]
        dρdx, dρdy, dρdz = grad_vars[1, 7], grad_vars[2, 7], grad_vars[3, 7]
        # virtual potential temperature gradient: for richardson calculation
        # strains
        # --------------------------------------------
        S11 = dudx
        S22 = dvdy
        S33 = dwdz
        S12 = (dudy + dvdx) / 2
        S13 = (dudz + dwdx) / 2
        S23 = (dvdz + dwdy) / 2
        SijSij = (S11^2 + S22^2 + S33^2
                  + 2.0 * S12^2
                  + 2.0 * S13^2 
                  + 2.0 * S23^2) 
        modSij = sqrt(2.0 * SijSij)
        
        #--------------------------------------------
        # deviatoric stresses
        VF[_τ11] = 2 * (S11 - (S11 + S22 + S33) / 3)
        VF[_τ22] = 2 * (S22 - (S11 + S22 + S33) / 3)
        VF[_τ33] = 2 * (S33 - (S11 + S22 + S33) / 3)
        VF[_τ12] = 2 * S12
        VF[_τ13] = 2 * S13
        VF[_τ23] = 2 * S23
        
        VF[_Tx], VF[_Ty], VF[_Tz] = dTdx, dTdy, dTdz
        VF[_ρx], VF[_ρy], VF[_ρz] = dρdx, dρdy, dρdz
        VF[_SijSij] = SijSij
    end
end

# -------------------------------------------------------------------------
# generic bc for 2d , 3d

@inline function bcstate!(QP, VFP, auxP, nM, QM, VFM, auxM, bctype, t,
                          PM, uM, vM, wM, w_rainM, ρinvM, q_totM, q_liqM, q_raiM)
    @inbounds begin
        x, y, z = auxM[_a_x], auxM[_a_y], auxM[_a_z]
        ρM, UM, VM, WM, EM, QTM = QM[_ρ], QM[_U], QM[_V], QM[_W], QM[_E], QM[_QT]
        # No flux boundary conditions
        # No shear on walls (free-slip condition)
        UnM = nM[1] * UM + nM[2] * VM + nM[3] * WM
        QP[_U] = UM - 2 * nM[1] * UnM
        QP[_V] = VM - 2 * nM[2] * UnM
        QP[_W] = WM - 2 * nM[3] * UnM
        VFP .= 0 
        nothing
    end
end

"""
Boundary correction for Neumann boundaries
"""
@inline function stresses_boundary_penalty!(VF, _...) 
  compute_stresses!(VF, 0) 
end


"""
Gradient term flux correction 
"""
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

@inline function source_geopot!(S,Q,aux,t)
    gravity::eltype(Q) = grav
    @inbounds begin
        ρ, U, V, W, E  = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E]
        S[_V] += - ρ * gravity
    end
end


@inline function source_sponge!(S,Q,aux,t)
    @inbounds begin
        ρu, ρv, ρw  = Q[_U], Q[_V], Q[_W]
        
        beta   = aux[_a_sponge]
        S[_U] -= beta * ρu
        S[_V] -= beta * ρv
        S[_W] -= beta * ρw
    end
end


# ------------------------------------------------------------------
# -------------END DEF SOURCES-------------------------------------# 

# initial condition
#function rising_bubble!(dim, Q, t, x, y, z, _...)
function squall_line!(dim, Q, t, spl_tinit, spl_qinit, spl_uinit, spl_vinit,
                 spl_pinit, x, y, z, _...)
    DFloat         = eltype(Q)
 
    DFloat                = eltype(Q)
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
    θ_ref::DFloat         = 303.0
    θ_c::DFloat           =   3.0
    Δθ::DFloat            =   0.0
    a::DFloat             =  50.0
    s::DFloat             = 100.0

    # INITIALISE ARRAYS FOR INTERPOLATED VALUES
    # --------------------------------------------------
    xvert          = z

    datat          = DFloat(spl_tinit(xvert))
    dataq          = DFloat(spl_qinit(xvert))
    datau          = DFloat(spl_uinit(xvert))
    datav          = DFloat(spl_vinit(xvert))
    datap          = DFloat(spl_pinit(xvert))
    dataq          = dataq / 1000

    if xvert >= 14000
        dataq = 0.0
    end

    rx  = 350  #10000.0
    ry  = 350  # 1200.0
    rz  = 250  # 1500.0
    xc  = 0.5*(xmax + xmin)
    yc  = 0.5*(ymax + ymin)
    zc  = 400 #2000.0
    
    r   = sqrt( (x - xc)^2/rx^2 + 0*(y - yc)^2/ry^2 + (z - zc)^2/rz^2)
    Δθ  = 0.0
    if r <= 1.0
        Δθ = θ_c * (cospi(0.5*r))^2
    end

    θ_liq = datat + Δθ
    q_tot = dataq
    p     = datap
    T     = air_temperature_from_liquid_ice_pottemp(θ_liq, p, PhasePartition(q_tot))
    ρ     = air_density(T, p)

    # energy definitions
    u, v, w     = datau, datav, zero(DFloat) #geostrophic. TO BE BUILT PROPERLY if Coriolis is considered
    
    U      = ρ * u
    V      = ρ * v
    W      = ρ * w
    e_kin  = (u^2 + v^2 + w^2) / 2
    e_pot  = grav * xvert
    e_int  = internal_energy(T, PhasePartition(q_tot))
    E      = ρ * total_energy(e_kin, e_pot, T, PhasePartition(q_tot))
    QT     = ρ * q_tot

     @inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT], Q[_QL], Q[_QR] = ρ, U, V, W, E, QT, DFloat(0), DFloat(0)
    
    #=
    qvar                  = PhasePartition(q_tot)
    θ                     = θ_ref + Δθ # potential temperature
    π_exner               = 1.0 - gravity / (c_p * θ) * z # exner pressure
    ρ                     = p0 / (R_gas * θ) * (π_exner)^ (c_v / R_gas) # density

    P                     = p0 * (R_gas * (ρ * θ) / p0) ^(c_p/c_v) # pressure (absolute)
    T                     = P / (ρ * R_gas) # temperature
    U, V, W               = 0.0 , 0.0 , 0.0  # momentum components
    # energy definitions
    e_kin                 = (U^2 + V^2 + W^2) / (2*ρ)/ ρ
    e_pot                 = gravity * z
    e_int                 = internal_energy(T, qvar)
    E                     = ρ * (e_int + e_kin + e_pot)  #* total_energy(e_kin, e_pot, T, q_tot, q_liq, q_ice)
    @inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E], Q[_QT]= ρ, U, V, W, E, ρ * q_tot
    =#
end


function run(mpicomm, dim, Ne, N, timeend, DFloat, dt)

    #-----------------------------------------------------------------
    # build physical range to be stratched
    #-----------------------------------------------------------------
    x_range = range(DFloat(xmin), length=Ne[1]   + 1, DFloat(xmax))
    y_range = range(DFloat(ymin), length=Ne[2]   + 1, DFloat(ymax))
    z_range = range(DFloat(zmin), length=Ne[end] + 1, DFloat(zmax))

    
    #-----------------------------------------------------------------
    # Build grid stretching along each direction
    # (ONLY Z for now. We need to decide what function we want to use for x and y)
    #-----------------------------------------------------------------
    xstretch_flg, ystretch_flg, zstretch_flg = 0, 0, 1
    xstretch_coe, ystretch_coe, zstretch_coe = 1.5, 1.5, 2.6
    (x_range_stretched, y_range_stretched, z_range_stretched) = grid_stretching_cube(xmin, xmax,
                                                                                     ymin, ymax,
                                                                                     zmin, zmax,
                                                                                     Ne,
                                                                                     xstretch_flg, ystretch_flg, zstretch_flg,
                                                                                     xstretch_coe, ystretch_coe, zstretch_coe,
                                                                                     dim)
    
    #-----------------------------------------------------------------
    x_range, y_range, z_range = x_range_stretched, y_range_stretched, z_range_stretched
    #
    # END grid stretching 
    #-----------------------------------------------------------------


    
    #-----------------------------------------------------------------
    #Build grid:
    #-----------------------------------------------------------------
    brickrange = (x_range, y_range, z_range)
    
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
                             auxiliary_state_initialization! =
                             auxiliary_state_initialization!,
                             source! = source!)


    # This is a actual state/function that lives on the grid
    #@timeit to "IC init" begin
    # ----------------------------------------------------
    # GET DATA FROM INTERPOLATED ARRAY ONTO VECTORS
    # This driver accepts data in 6 column format
    # ----------------------------------------------------
    (sounding, _, ncols) = read_sounding()

    # WARNING: Not all sounding data is formatted/scaled
    # the same. Care required in assigning array values
    # height theta qv    u     v     pressure
    zinit, tinit, qinit, uinit, vinit, pinit  =
        sounding[:, 1], sounding[:, 2], sounding[:, 3], sounding[:, 4], sounding[:, 5], sounding[:, 6]
    #------------------------------------------------------
    # GET SPLINE FUNCTION
    #------------------------------------------------------
    spl_tinit    = Spline1D(zinit, tinit; k=1)
    spl_qinit    = Spline1D(zinit, qinit; k=1)
    spl_uinit    = Spline1D(zinit, uinit; k=1)
    spl_vinit    = Spline1D(zinit, vinit; k=1)
    spl_pinit    = Spline1D(zinit, pinit; k=1)
    
    # This is a actual state/function that lives on the grid
    #initialcondition(Q, x...) = rising_bubble!(Val(dim), Q, DFloat(0), x...)
    initialcondition(Q, x...) = squall_line!(Val(dim), Q, DFloat(0), spl_tinit,
                                            spl_qinit, spl_uinit, spl_vinit,
                                            spl_pinit, x...)
    Q = MPIStateArray(spacedisc, initialcondition)

    lsrk = LSRK54CarpenterKennedy(spacedisc, Q; dt = dt, t0 = 0)

    eng0 = norm(Q)
    @info @sprintf """Starting
      norm(Q₀) = %.16e""" eng0

    # Set up the information callback
    starttime = Ref(now())
    cbinfo = GenericCallbacks.EveryXWallTimeSeconds(10, mpicomm) do (s=false)
        if s
            starttime[] = now()
        else
            energy = norm(Q)
            #globmean = global_mean(Q, _ρ)
            @info @sprintf("""Update
                         simtime = %.16e
                         runtime = %s
                         norm(Q) = %.16e""", 
                           ODESolvers.gettime(lsrk),
                           Dates.format(convert(Dates.DateTime,
                                                Dates.now()-starttime[]),
                                        Dates.dateformat"HH:MM:SS"),
                           energy )#, globmean)
        end
    end

    npoststates = 8
    _P, _u, _v, _w, _ρinv, _q_liq, _T, _θ = 1:npoststates
    postnames = ("P", "u", "v", "w", "rhoinv", "_q_liq", "T", "THETA")
    postprocessarray = MPIStateArray(spacedisc; nstate=npoststates)

    step = [0]
    mkpath("vtk-RTB")
    cbvtk = GenericCallbacks.EveryXSimulationSteps(1000) do (init=false)
        DGBalanceLawDiscretizations.dof_iteration!(postprocessarray, spacedisc,
                                                   Q) do R, Q, QV, aux
                                                       @inbounds let
                                                           (R[_P], R[_u], R[_v], R[_w], R[_ρinv], R[_q_liq], R[_T], R[_θ]) = (preflux(Q, QV, aux))
                                                       end
                                                   end

        outprefix = @sprintf("vtk-RTB/cns_%dD_mpirank%04d_step%04d", dim,
                             MPI.Comm_rank(mpicomm), step[1])
        @debug "doing VTK output" outprefix
        writevtk(outprefix, Q, spacedisc, statenames,
                 postprocessarray, postnames)
        
        step[1] += 1
        nothing
    end
    
    solve!(Q, lsrk; timeend=timeend, callbacks=(cbinfo, cbvtk))


    # Print some end of the simulation information
    engf = norm(Q)
    if integration_testing
        Qe = MPIStateArray(spacedisc,
                           (Q, x...) -> initialcondition!(Val(dim), Q,
                                                          DFloat(timeend), x...))
        engfe = norm(Qe)
        errf = euclidean_distance(Q, Qe)
        @info @sprintf """Finished
        norm(Q)                 = %.16e
        norm(Q) / norm(Q₀)      = %.16e
        norm(Q) - norm(Q₀)      = %.16e
        norm(Q - Qe)            = %.16e
        norm(Q - Qe) / norm(Qe) = %.16e
        """ engf engf/eng0 engf-eng0 errf errf / engfe
    else
        @info @sprintf """Finished
        norm(Q)            = %.16e
        norm(Q) / norm(Q₀) = %.16e
        norm(Q) - norm(Q₀) = %.16e""" engf engf/eng0 engf-eng0
    end
integration_testing ? errf : (engf / eng0)
end

using Test
let
    MPI.Initialized() || MPI.Init()
    Sys.iswindows() || (isinteractive() && MPI.finalize_atexit())
    mpicomm = MPI.COMM_WORLD
    if MPI.Comm_rank(mpicomm) == 0
        ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
        loglevel = ll == "DEBUG" ? Logging.Debug :
            ll == "WARN"  ? Logging.Warn  :
            ll == "ERROR" ? Logging.Error : Logging.Info
        global_logger(ConsoleLogger(stderr, loglevel))
    else
        global_logger(NullLogger())
    end
    # User defined number of elements
    # User defined timestep estimate
    # User defined simulation end time
    # User defined polynomial order 
    numelem = (Nex, Ney, Nez)
    dt = 0.005
    timeend = 900
    polynomialorder = Npoly
    DFloat = Float64
    dim = numdims

    
    DoF = (Nex*Ney*Nez)*(Npoly+1)^numdims*(_nstate)
    DoFstorage = (Nex*Ney*Nez) *
        (Npoly+1)^numdims *
        (_nstate + _nviscstates + _nauxstate + CLIMA.Grids._nvgeo) +
        (Nex*Ney*Nez) * (Npoly+1)^(numdims-1) *
        2^numdims*(CLIMA.Grids._nsgeo)
    

    @info @sprintf """ ----------------------------------------------------"""
    @info @sprintf """   ______ _      _____ __  ________                  """     
    @info @sprintf """  |  ____| |    |_   _|  ...  |  __  |               """  
    @info @sprintf """  | |    | |      | | |   .   | |  | |               """ 
    @info @sprintf """  | |    | |      | | | |   | | |__| |               """
    @info @sprintf """  | |____| |____ _| |_| |   | | |  | |               """
    @info @sprintf """  | _____|______|_____|_|   |_|_|  |_|               """
    @info @sprintf """                                                     """
    @info @sprintf """ ----------------------------------------------------"""
    @info @sprintf """ 3D Rising Bubble with GRID STRETCHING               """
    @info @sprintf """   Resolution:                                         """
    @info @sprintf """     (Δx, Δy, Δz)   = (%.2e, %.2e, %.2e)               """ Δx Δy Δz
    @info @sprintf """     (Nex, Ney, Nez) = (%d, %d, %d)                    """ Nex Ney Nez
    @info @sprintf """     DoF = %d                                          """ DoF
    @info @sprintf """     Minimum necessary memory to run this test: %g GBs """ (DoFstorage * sizeof(DFloat))/1000^3
    @info @sprintf """     Time step dt: %.2e                                """ dt
    @info @sprintf """     End time  t : %.2e                                """ timeend
    @info @sprintf """ ------------------------------------------------------"""
    
    
    engf_eng0 = run(mpicomm, dim, numelem[1:dim], polynomialorder, timeend,
                    DFloat, dt)
end

isinteractive() || MPI.Finalize()

nothing
