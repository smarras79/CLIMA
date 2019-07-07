
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

# Prognostic equations: ρ, (U), (V), (W), (E_tot), (Q_tot)
# For the dry example shown here, we load the moist thermodynamics module
# and consider the dry equation set to be the same as the moist equations but
# with total specific humidity = 0.
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters
using CLIMA.Microphysics

# State labels
const _nstate = 9
const _ρ, _U, _V, _W, _E_tot, _Q_tot, _Q_liq, _Q_ice, _Q_rai =1:_nstate
const stateid = (ρid = _ρ, U_id = _U, V_id = _V, W_id = _W,
                 E_tot_id = _E_tot, Q_tot_id = _Q_tot, Q_liq_id = _Q_liq,
                 Q_ice_id = _Q_ice, Q_rai_id = _Q_rai)
const statenames = ("ρ", "U", "V", "W", "E_tot", "Q_tot", "Q_liq",
                    "Q_ice", "Q_rai")

# Viscous state labels
const _nviscstates = 22
const _τ11, _τ22, _τ33, _τ12, _τ13, _τ23,
      _q_tot_x, _q_tot_y, _q_tot_z,
      _q_liq_x, _q_liq_y, _q_liq_z,
      _q_ice_x, _q_ice_y, _q_ice_z,
      _q_rai_x, _q_rai_y, _q_rai_z,
      _Tx, _Ty, _Tz, _SijSij = 1:_nviscstates

# Gradient state labels
const _ngradstates = 9
const _states_for_gradient_transform = (_ρ, _U, _V, _W, _E_tot, _Q_tot,
                                        _Q_liq, _Q_ice, _Q_rai)

const _nauxstate = 11
const _a_z, _a_dx, _a_dy, _a_dz, _a_sponge, _a_02z, _a_z2inf, _a_T, _a_p, _a_soundspeed_air, _a_timescale = 1:_nauxstate

if !@isdefined integration_testing
    const integration_testing =
        parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end

# Problem constants (TODO: parameters module (?))
@parameter Prandtl_t 1//3 "Prandtl_t"
@parameter cp_over_prandtl cp_d / Prandtl_t "cp_over_prandtl"

# Random number seed
const seed = MersenneTwister(0)

# Problem description
# --------------------
# 2D thermal perturbation (cold bubble) in a neutrally stratified atmosphere
# No wall-shear, lateral periodic boundaries with no-flux walls at the domain
# top and bottom.
# Inviscid, Constant viscosity, StandardSmagorinsky, MinimumDissipation
# filters are tested against this benchmark problem
# TODO: link to module SubGridScaleTurbulence

#
# User Input
#
const numdims = 3
const Npoly = 4

#
# Define grid size
#
Δx    =  1000
Δy    =  1000
Δz    =  500

#
# OR:
#
# Set Δx < 0 and define  Nex, Ney, Nez:
#
(Nex, Ney, Nez) = (40, 30, 24)

# Physical domain extents
const (xmin, xmax) = (-20000, 20000)
const (ymin, ymax) = (-20000, 20000)
const (zmin, zmax) = (     0, 24000)

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


DoF = (Nex*Ney*Nez)*(Npoly+1)^numdims*(_nstate)
DoFstorage = (Nex*Ney*Nez) *
             (Npoly+1)^numdims *
             (_nstate + _nviscstates + _nauxstate + CLIMA.Grids._nvgeo) +
             (Nex*Ney*Nez) * (Npoly+1)^(numdims-1) *
             2^numdims*(CLIMA.Grids._nsgeo)


# Smagorinsky model requirements : TODO move to SubgridScaleTurbulence module
@parameter C_smag 0.15 "C_smag"
# Equivalent grid-scale
Δ = (Δx * Δy * Δz)^(1/3)
const Δsqr = Δ * Δ

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
@inline function preflux(Q, VF, aux, _...)

    DFloat = eltype(Q)

    @inbounds begin

      # unpack model variables
        ρ, U, V, W, E_tot          = Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E_tot]
        Q_tot, Q_liq, Q_ice, Q_rai = Q[_Q_tot], Q[_Q_liq], Q[_Q_ice], Q[_Q_rai]
        u, v, w                    = U/ρ, V/ρ, W/ρ       
        e_tot                      = E_tot/ρ
        q_tot, q_liq, q_ice, q_rai = Q_tot/ρ, Q_liq/ρ, Q_ice/ρ, Q_rai/ρ


      # compute rain fall speed
      DF = eltype(ρ)
      #if(q_rai >= DF(0)) #TODO - need a way to prevent negative values
      #  rain_w = terminal_velocity(q_rai, ρ)
      #else
        rain_w = Float64(0)
      #end

      return (u, v, w, rain_w, ρ, q_tot, q_liq, q_ice, q_rai, e_tot)
    end
end

#-------------------------------------------------------------------------
#md # Soundspeed computed using the thermodynamic state TS
# max eigenvalue
@inline function wavespeed(n, Q, aux, t, u, v, w, rain_w, ρ,
                           q_tot, q_liq, q_ice, q_rai, e_tot)
    @inbounds begin
        (n[1] * u + n[2] * v + n[3] * max(abs(w), abs(rain_w), abs(w-rain_w))) +
          aux[_a_soundspeed_air]
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
    fsounding  = open(joinpath(@__DIR__, "../soundings/sounding_gabersek.dat"))
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
@inline function cns_flux!(F, Q, VF, aux, t, u, v, w, rain_w, ρ,
                           q_tot, q_liq, q_ice, q_rai, e_tot)
    @inbounds begin

        DFloat = eltype(F)
        p = aux[_a_p]

        # Inviscid contributions
        F[1, _ρ],  F[2, _ρ],  F[3, _ρ]  = ρ * u, ρ * v, ρ * w

        F[1, _U],  F[2, _U],  F[3, _U]  = u * ρ * u  + p , v * ρ * u      , w * ρ * u
        F[1, _V],  F[2, _V],  F[3, _V]  = u * ρ * v      , v * ρ * v + p  , w * ρ * v
        F[1, _W],  F[2, _W],  F[3, _W]  = u * ρ * w      , v * ρ * w      , w * ρ * w + p

        F[1, _E_tot],  F[2, _E_tot],  F[3, _E_tot]  = u * (ρ * e_tot + p), v * (ρ * e_tot + p), w * (ρ * e_tot + p)

        F[1, _Q_tot], F[2, _Q_tot], F[3, _Q_tot] = u * ρ * q_tot, v * ρ * q_tot, w * ρ * q_tot
        F[1, _Q_liq], F[2, _Q_liq], F[3, _Q_liq] = 0.0 #u * ρ * q_liq, v * ρ * q_liq, w * ρ * q_liq
        F[1, _Q_ice], F[2, _Q_ice], F[3, _Q_ice] = 0.0 #u * ρ * q_ice, v * ρ * q_ice, w * ρ * q_ice
        F[1, _Q_rai], F[2, _Q_rai], F[3, _Q_rai] = 0.0 #u * ρ * q_rai, v * ρ * q_rai, (w - rain_w) * ρ * q_rai

        #Derivative of q_tot, q_liq, q_ice, q_rai, T:
        vq_tot_x, vq_tot_y, vq_tot_z = VF[_q_tot_x], VF[_q_tot_y], VF[_q_tot_z]
        vq_liq_x, vq_liq_y, vq_liq_z = VF[_q_liq_x], VF[_q_liq_y], VF[_q_liq_z]
        vq_ice_x, vq_ice_y, vq_ice_z = VF[_q_ice_x], VF[_q_ice_y], VF[_q_ice_z]
        vq_rai_x, vq_rai_y, vq_rai_z = VF[_q_rai_x], VF[_q_rai_y], VF[_q_rai_z]
        vTx, vTy, vTz = VF[_Tx], VF[_Ty], VF[_Tz]

        
        SijSij = VF[_SijSij]

        #Dynamic eddy viscosity from Smagorinsky:
        ν_e    = sqrt(2SijSij) * C_smag^2 * DFloat(Δsqr)
        D_e    = 200 #ν_e / Prandtl_t
        D_e_ql = 200
        D_e_qi = 200
        D_e_qr = 200

        # Multiply stress tensor by viscosity coefficient:
        τ11, τ22, τ33 = VF[_τ11] * ν_e, VF[_τ22]* ν_e, VF[_τ33] * ν_e
        τ12 = τ21 = VF[_τ12] * ν_e
        τ13 = τ31 = VF[_τ13] * ν_e
        τ23 = τ32 = VF[_τ23] * ν_e

        # Viscous velocity flux (i.e. F^visc_u in Giraldo Restelli 2008)
        F[1, _U] -= τ11; F[2, _U] -= τ12; F[3, _U] -= τ13
        F[1, _V] -= τ21; F[2, _V] -= τ22; F[3, _V] -= τ23
        F[1, _W] -= τ31; F[2, _W] -= τ32; F[3, _W] -= τ33

        # Viscous Energy flux (i.e. F^visc_e in Giraldo Restelli 2008)
        # TODO should it also depend on q_tot gradients?
        F[1, _E_tot] -= u * τ11 + v * τ12 + w * τ13 + cp_over_prandtl * vTx * ν_e
        F[2, _E_tot] -= u * τ21 + v * τ22 + w * τ23 + cp_over_prandtl * vTy * ν_e
        F[3, _E_tot] -= u * τ31 + v * τ32 + w * τ33 + cp_over_prandtl * vTz * ν_e

        # Viscous contributions to mass flux terms
        F[1, _Q_tot] -= vq_tot_x * D_e; F[2, _Q_tot] -= vq_tot_y * D_e; F[3, _Q_tot] -= vq_tot_z * D_e
        F[1, _Q_liq] -= vq_liq_x * D_e; F[2, _Q_liq] -= vq_liq_y * D_e; F[3, _Q_liq] -= vq_liq_z * D_e
        F[1, _Q_ice] -= vq_ice_x * D_e; F[2, _Q_ice] -= vq_ice_y * D_e; F[3, _Q_ice] -= vq_ice_z * D_e
        F[1, _Q_rai] -= vq_rai_x * D_e; F[2, _Q_rai] -= vq_rai_y * D_e; F[3, _Q_rai] -= vq_rai_z * D_e
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
@inline function gradient_vars!(vel, Q, aux, t, u, v, w, rain_w, ρ,
                           q_tot, q_liq, q_ice, q_rai, e_tot)

    @inbounds begin
        T = aux[_a_T]

        # TODO
        # ordering should match states_for_gradient_transform
        #_states_for_gradient_transform = (_ρ, _U, _V, _W, _E_tot, _Q_tot, _Q_liq, _Q_ice, _Q_rai)

        vel[1], vel[2], vel[3] = u, v, w

        vel[4], vel[5], vel[6], vel[7], vel[8], vel[9]  = e_tot, q_tot, q_liq, q_ice, q_rai, T
    end
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
        # TODO - missing e_tot?

        # compute gradients of moist vars and temperature
        dq_tot_dx, dq_tot_dy, dq_tot_dz = grad_vel[1, 5], grad_vel[2, 5], grad_vel[3, 5]
        dq_liq_dx, dq_liq_dy, dq_liq_dz = grad_vel[1, 6], grad_vel[2, 6], grad_vel[3, 6]
        dq_ice_dx, dq_ice_dy, dq_ice_dz = grad_vel[1, 7], grad_vel[2, 7], grad_vel[3, 7]
        dq_rai_dx, dq_rai_dy, dq_rai_dz = grad_vel[1, 8], grad_vel[2, 8], grad_vel[3, 8]
        dTdx,      dTdy,      dTdz      = grad_vel[1, 9], grad_vel[2, 9], grad_vel[3, 9]

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
        VF[_q_tot_x], VF[_q_tot_y], VF[_q_tot_z] = dq_tot_dx, dq_tot_dy, dq_tot_dz
        VF[_q_liq_x], VF[_q_liq_y], VF[_q_liq_z] = dq_liq_dx, dq_liq_dy, dq_liq_dz
        VF[_q_ice_x], VF[_q_ice_y], VF[_q_ice_z] = dq_ice_dx, dq_ice_dy, dq_ice_dz
        VF[_q_rai_x], VF[_q_rai_y], VF[_q_rai_z] = dq_rai_dx, dq_rai_dy, dq_rai_dz
        VF[_Tx],      VF[_Ty],      VF[_Tz]      = dTdx,      dTdy,      dTdz
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
        DFloat = eltype(aux)
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

        ##BEGIN  User modification on domain parameters.
        ##Only change the first index of brickrange if your axis are
        ##oriented differently:
        ##x, y, z = aux[_a_x], aux[_a_y], aux[_a_z]
        
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
# generic bc for 2d , 3d
#
@inline function bcstate!(QP, VFP, auxP, nM, QM, VFM, auxM, bctype, t,
                          u, v, w, rain_wM, ρ, q_tot, q_liq, q_ice, q_rai, e_tot)
    @inbounds begin
        U_M, V_M, W_M = QM[_U], QM[_V], QM[_W]
        # No flux boundary conditions
        # No shear on walls (free-slip condition)
        U_nM = nM[1] * U_M + nM[2] * V_M + nM[3] * W_M
        QP[_U] = U_M - 2 * nM[1] * U_nM
        QP[_V] = V_M - 2 * nM[2] * U_nM
        QP[_W] = W_M - 2 * nM[3] * U_nM
        VFP .= 0
        nothing
    end
end

# -------------------------------------------------------------------------
@inline function stresses_boundary_penalty!(VF, _...)
    #VF .= 0
    compute_stresses!(VF, 0) #
end

@inline function stresses_penalty!(VF, nM, velM, QM, aM, velP, QP, aP, t)
    @inbounds begin
        n_Δvel = similar(VF, Size(3, 3))
        for j = 1:_ngradstates, i = 1:3
            n_Δvel[i, j] = nM[i] * (velP[j] - velM[j]) / 2
        end
        compute_stresses!(VF, n_Δvel)
    end
end
# -------------------------------------------------------------------------

source!(S, Q, aux, t) = source!(S, Q, aux, t, preflux(Q, ~, aux)...)
@inline function source!(S, Q, aux, t, u, v, w, rain_w, ρ,
                         q_tot, q_liq, q_ice, q_rai, e_tot)
    # Initialise the final block source term
    S .= 0

    # Typically these sources are imported from modules
    @inbounds begin
        #source_microphysics!(S, Q, aux, t, u, v, w, rain_w, ρ,
        #                     q_tot, q_liq, q_ice, q_rai, e_tot)
        source_geopot!(S, Q, aux, t)
        #source_sponge!(S, Q, aux, t)
    end
end

@inline function source_microphysics!(S, Q, aux, t, u, v, w, rain_w, ρ,
                             q_tot, q_liq, q_ice, q_rai, e_tot)

  DF = eltype(Q)

  @inbounds begin

    z = aux[_a_z]
    p = aux[_a_p]

    #TODO - tmp
    q_tot = max(DF(0), q_tot)
    q_liq = max(DF(0), q_liq)
    q_ice = max(DF(0), q_ice)
    q_rai = max(DF(0), q_rai)

    # current state
    e_int = e_tot - 1//2 * (u^2 + v^2 + w^2) - grav * z
    q     = PhasePartition(q_tot, q_liq, q_ice)
    T     = air_temperature(e_int, q)
    # equilibrium state at current T
    q_eq = PhasePartition_equil(T, ρ, q_tot)

    # cloud water condensation/evaporation
    src_q_liq = conv_q_vap_to_q_liq(q_eq, q)
    #src_q_ice = conv_q_vap_to_q_ice(q_eq, q)
    S[_Q_liq] += ρ * src_q_liq
    #S[_Q_ice] += ρ * src_q_ice

    # tendencies from rain
    # TODO - ensure positive definite
    # TODO - temporary handling ice
    #if(q_tot >= DF(0) && q_liq >= DF(0) && q_rai >= DF(0))

    src_q_rai_evap = conv_q_rai_to_q_vap(q_rai, q, T , p, ρ)

    src_q_rai_acnv_liq = conv_q_liq_to_q_rai_acnv(q.liq)
    src_q_rai_accr_liq = conv_q_liq_to_q_rai_accr(q.liq, q_rai, ρ)

    #src_q_rai_acnv_ice = conv_q_liq_to_q_rai_acnv(q.ice)
    #src_q_rai_accr_ice = conv_q_liq_to_q_rai_accr(q.ice, q_rai, ρ)

    src_q_rai_tot = src_q_rai_acnv_liq + src_q_rai_accr_liq + src_q_rai_evap# + src_q_rai_acnv_ice + src_q_rai_accr_ice

    S[_Q_liq] -= ρ * (src_q_rai_acnv_liq + src_q_rai_accr_liq)
    #S[_Q_ice] -= ρ * (src_q_rai_acnv_ice + src_q_rai_accr_ice)

    S[_Q_rai] += ρ * src_q_rai_tot
    S[_Q_tot] -= ρ * src_q_rai_tot

    S[_E_tot] -= (
                    src_q_rai_evap * (DF(cv_v) * (T - DF(T_0)) + e_int_v0) -
                    (src_q_rai_acnv_liq + src_q_rai_accr_liq) * DF(cv_l) * (T - DF(T_0))# -
                    #(src_q_rai_acnv_ice + src_q_rai_accr_ice) * DF(cv_i) * (T - DF(T_0))
                  ) * ρ
    #end
  end
end
"""
        Geostrophic wind forcing
        """
@inline function source_geostrophic!(S,Q,aux,t)
    DFloat = eltype(S)
    f_coriolis = DFloat(7.62e-5)
    u_geostrophic = DFloat(7)
    v_geostrophic = DFloat(-5.5)
    @inbounds begin
        ρ = Q[_ρ]
        U = Q[_U]
        V = Q[_V]
        S[_U] -= f_coriolis * (U - ρ * u_geostrophic)
        S[_V] -= f_coriolis * (U - ρ * v_geostrophic)
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

@inline function source_geopot!(S,Q,aux,t)
    @inbounds S[_W] += - Q[_ρ] * grav
end

# Test integral exactly according to the isentropic vortex example
@inline function integral_knl(val, Q, aux)
    κ = 85
    @inbounds begin
        ρ, Q_liq = Q[_ρ], Q[_Q_liq]
        q_liq = Q_liq / ρ
        val[1] = ρ * κ * q_liq
    end
end

function preodefun!(disc, Q, t)
    DGBalanceLawDiscretizations.dof_iteration!(disc.auxstate, disc, Q) do R, Q, QV, aux
        @inbounds let
            ρ, U, V, W, E_tot, Q_tot, Q_liq, Q_ice, Q_rai =
              Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E_tot], Q[_Q_tot], Q[_Q_liq],
              Q[_Q_ice], Q[_Q_rai]

            z = aux[_a_z]
            #dx, dy, dz = aux[_a_dx], aux[_a_dy], aux[_a_dz]

            q_tot = Q_tot / ρ; q_liq = Q_liq / ρ; q_ice = Q_ice / ρ
            u = U / ρ; v = V / ρ; w = W / ρ
            e_tot = E_tot / ρ

            e_int = e_tot - 1//2 * (u^2 + v^2 + w^2) - grav * z
            q     = PhasePartition(q_tot, q_liq, q_ice)
            T     = air_temperature(e_int, q)
            p     = air_pressure(T, ρ, q)
            
            
            R[_a_T] = T
            R[_a_p] = p
            R[_a_soundspeed_air] = 340 # soundspeed_air(T, q)
            #u_wavespeed = (abs(u) + soundspeed) / dx
            #v_wavespeed = (abs(v) + soundspeed) / dy 
            #w_wavespeed = (abs(w) + soundspeed) / dz
            #R[_a_timescale] = max(u_wavespeed,v_wavespeed,w_wavespeed)
        end
    end

    integral_computation(disc, Q, t)
end

function integral_computation(disc, Q, t)
    DGBalanceLawDiscretizations.indefinite_stack_integral!(disc, integral_knl, Q,
                                                           (_a_02z))
    DGBalanceLawDiscretizations.reverse_indefinite_stack_integral!(disc,
                                                                   _a_z2inf,
                                                                   _a_02z)
end

# initial condition
"""
            User-specified. Required.
            This function specifies the initial conditions
            for the dycoms driver.
        """
function squall_line!(dim, Q, t, spl_tinit, spl_qinit, spl_uinit, spl_vinit,
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

    if xvert >= 14000
        dataq = 0.0
    end

    θ_c =     5.0
    rx  = 10000.0
    ry  =  1200.0
    rz  =  1500.0
    xc  = 0.5*(xmax + xmin)
    yc  = 0.5*(ymax + ymin)
    zc  = 2000.0
    
    r   = sqrt( (x - xc)^2/rx^2 + (y - yc)^2/ry^2 + (z - zc)^2/rz^2)
    Δθ  = 0.0
    if r <= 1.0
        Δθ = θ_c * (cospi(0.5*r))^2
    end

    randnum1   = rand(seed, DFloat) / 100
        
    θ_liq = datat + Δθ*randnum1
    q_tot = dataq
    p     = datap
    T     = air_temperature_from_liquid_ice_pottemp(θ_liq, p, PhasePartition(q_tot))
    ρ     = air_density(T, p)

    # energy definitions
    u, v, w     = datau, datav, zero(DFloat) #geostrophic. TO BE BUILT PROPERLY if Coriolis is considered
    
    U          = ρ * u
    V          = ρ * v
    W          = ρ * w
    e_kin       = (u^2 + v^2 + w^2) / 2
    e_pot       = grav * xvert
    e_int       = internal_energy(T, PhasePartition(q_tot))
    E_tot      = ρ * total_energy(e_kin, e_pot, T, PhasePartition(q_tot))
    Q_tot      = ρ * q_tot

    @inbounds Q[_ρ], Q[_U], Q[_V], Q[_W], Q[_E_tot], Q[_Q_tot],
                Q[_Q_liq], Q[_Q_ice], Q[_Q_rai] = ρ, U, V, W, E_tot, Q_tot, DFloat(0), DFloat(0), DFloat(0)
end

function grid_stretching(DFloat,
                         xmin, xmax, ymin, ymax, zmin, zmax,
                         Ne,
                         xstretch_flg, ystretch_flg, zstretch_flg)

    #build physical range to be stratched
    x_range_stretched = (range(DFloat(xmin), length=Ne[1]+1, DFloat(xmax)))
    y_range_stretched = (range(DFloat(ymin), length=Ne[2]+1, DFloat(ymax)))
    z_range_stretched = (range(DFloat(zmin), length=Ne[3]+1, DFloat(zmax)))

    #build logical space
    ksi  = (range(DFloat(0), length=Ne[1]+1, DFloat(1)))
    eta  = (range(DFloat(0), length=Ne[2]+1, DFloat(1)))
    zeta = (range(DFloat(0), length=Ne[3]+1, DFloat(1)))

    xstretch_coe = 0.0
    if xstretch_flg == 1
        xstretch_coe = 1.5
        x_range_stretched = (xmax - xmin).*(exp.(xstretch_coe * ksi)  .- 1.0)./(exp(xstretch_coe) - 1.0)
    end

    ystretch_coe = 0.0
    if ystretch_flg == 1
        ystretch_coe = 1.5
        y_range_stretched = (ymax - ymin).*(exp.(ystretch_coe * eta)  .- 1.0)./(exp(ystretch_coe) - 1.0)
    end

    zstretch_coe = 0.0
    if zstretch_flg == 1
        zstretch_coe = 2.5
        z_range_stretched = (zmax - zmin).*(exp.(zstretch_coe * zeta) .- 1.0)./(exp(zstretch_coe) - 1.0)
    end

    return (x_range_stretched, y_range_stretched, z_range_stretched)

end

function run(mpicomm, dim, Ne, N, timeend, DFloat, dt)


    #Build stretching along each direction
    (x_range_stretched, y_range_stretched, z_range_stretched) = grid_stretching(DFloat, xmin, xmax, ymin, ymax, zmin, zmax, Ne, 0, 0, 1)

    #Build (stretched) grid:
    brickrange = (x_range_stretched, y_range_stretched, z_range_stretched)

@info " A"
    # User defined periodicity in the topl assignment
    # brickrange defines the domain extents
    #@timeit to "Topo init" topl = StackedBrickTopology(mpicomm, brickrange, periodicity=(true,true,false))
    topl = StackedBrickTopology(mpicomm, brickrange, periodicity=(true,true,false))
@info " B"
    #@timeit to "Grid init" grid = DiscontinuousSpectralElementGrid(topl,
    grid = DiscontinuousSpectralElementGrid(topl,
                                                                   FloatType = DFloat,
                                                                   DeviceArray = ArrayType,
                                                                   polynomialorder = N)
@info " C"
    numflux!(x...) = NumericalFluxes.rusanov!(x..., cns_flux!, wavespeed, preflux)
    @info " D"
    numbcflux!(x...) = NumericalFluxes.rusanov_boundary_flux!(x..., cns_flux!, bcstate!, wavespeed, preflux)
@info " E"
    # spacedisc = data needed for evaluating the right-hand side function
    #@timeit to "Space Disc init" spacedisc = DGBalanceLaw(grid = grid,
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
                                                          preodefun! = preodefun!)
@info " F"
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

        initialcondition(Q, x...) = squall_line!(Val(dim), Q, DFloat(0), spl_tinit,
                                            spl_qinit, spl_uinit, spl_vinit,
                                            spl_pinit, x...)
    Q = MPIStateArray(spacedisc, initialcondition)
    @info " G"
    #end

    #@timeit to "Time stepping init" begin
        
        lsrk = LSRK54CarpenterKennedy(spacedisc, Q; dt = dt, t0 = 0)
    @info " H"
        
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
@info " III"
npoststates = 6
_u, _v, _w, _q_tot, _q_liq, _q_rai = 1:npoststates
postnames = ("u", "v", "w", "q_tot", "q_liq", "q_rai")
postprocessarray = MPIStateArray(spacedisc; nstate=npoststates)

step = [0]
mkpath("vtk-RTB")
cbvtk = GenericCallbacks.EveryXSimulationSteps(1) do (init=false)
    @info " III1"
    #=DGBalanceLawDiscretizations.dof_iteration!(postprocessarray, spacedisc,
                                               Q) do R, Q, QV, aux
                                                   @inbounds let
                                                       (R[_u], R[_v], R[_w], R[_w_rain], R[_ρ], R[_q_tot], R[_q_liq], R[_q_ice], R[_q_rai], R[_e_tot]) = (preflux(Q, QV, aux))
                                                   end
                                               end
    @info " III2"
    outprefix = @sprintf("vtk-RTB/cns_%dD_mpirank%04d_step%04d", dim,
                         MPI.Comm_rank(mpicomm), step[1])
    @debug "doing VTK output" outprefix
    writevtk(outprefix, Q, spacedisc, statenames,
             postprocessarray, postnames)
    @info " III3"
=#
    step[1] += 1
    nothing
end
@info " L"

        #npoststates = 18
        #out_z, out_u, out_v, out_w, out_e_tot, out_e_int, out_e_kin, out_e_pot, out_p, out_beta, out_T, out_q_tot, out_q_vap, out_q_liq, out_q_ice, out_q_rai, out_rain_w, out_tht = 1:npoststates
        #postnames = ("height", "u", "v", "w", "e_tot", "e_int", "e_kin", "e_pot", "p", "beta", "T", "q_tot", "q_vap", "q_liq", "q_ice", "q_rai", "rain_w", "theta")

#=
        npoststates = 9
        out_u, out_v, out_w, out_T, out_q_tot, out_q_vap, out_q_liq, out_q_ice, out_q_rai = 1:npoststates
        postnames = ("u", "v", "w", "T", "q_tot", "q_vap", "q_liq", "q_ice", "q_rai")
        
        postprocessarray = MPIStateArray(spacedisc; nstate=npoststates)
@info " F"
        step = [0]
        mkpath("./CLIMA-output-scratch/vtk-sq-3d")
        cbvtk = GenericCallbacks.EveryXSimulationSteps(3600) do (init=false) #every 1 min = (0.025) * 40 * 60 * 1min
            DGBalanceLawDiscretizations.dof_iteration!(postprocessarray, spacedisc, Q) do R, Q, QV, aux
                @inbounds let
                    DF = eltype(Q)
              @info " F1"      
                    u, v, w, rain_w, ρ, q_tot, q_liq, q_ice, q_rai, e_tot =
                      preflux(Q, QV, aux)
@info " F2"
                    e_kin = 1//2 * (u^2 + v^2 + w^2)
                    @info " F3"
                    e_pot = grav * aux[_a_z]
                    @info " F4"
                    e_int = e_tot - e_kin - e_pot
                    @info " F5"
                    q = PhasePartition(q_tot, q_liq, q_ice)
@info " F6"
                    T = air_temperature(e_int, q)
                    @info " F7"
                    p = aux[_a_p]
                    @info " F8"
                    tht = liquid_ice_pottemp(T, p, q)
@info " F9"
                    R[out_u] = u
                    R[out_v] = v
                    R[out_w] = w
@info " F10"
                    R[out_q_tot] = q_tot
                    R[out_q_liq] = q_liq
                    R[out_q_rai] = q_rai
          @info " F11"          
                end
            end
@info " G"
            outprefix = @sprintf("./CLIMA-output-scratch/vtk-sq-3d/squall_%dD_mpirank%04d_step%04d", dim,
                                 MPI.Comm_rank(mpicomm), step[1])
            @debug "doing VTK output" outprefix
            writevtk(outprefix, Q, spacedisc, statenames,
                     postprocessarray, postnames)

            step[1] += 1
            nothing
        end
    #end
=#

#@info @sprintf """Starting...
#            norm(Q) = %25.16e""" norm(Q)
 @info " M"
# Initialise the integration computation. Kernels calculate this at every timestep??
#@timeit to "initial integral" integral_computation(spacedisc, Q, 0)
#@timeit to "solve" solve!(Q, lsrk; timeend=timeend, callbacks=(cbinfo, cbvtk))
solve!(Q, lsrk; timeend=timeend, callbacks=(cbinfo, cbvtk))


#@info @sprintf """Finished...
#            norm(Q) = %25.16e""" norm(Q)
 @info " N"
#=
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
=#

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
         @info " O"
    else
         @info " P"
        global_logger(NullLogger())
         @info " Q"
    end
    # User defined number of elements
    # User defined timestep estimate
    # User defined simulation end time
    # User defined polynomial order
    numelem = (Nex, Ney, Nez)
    dt = 0.01
    timeend = 12000 # 2h 30 min
    polynomialorder = Npoly
    DFloat = Float64
    dim = numdims

    if MPI.Comm_rank(mpicomm) == 0
        @info @sprintf """ ------------------------------------------------------"""
        @info @sprintf """   ______ _      _____ __  ________                    """
        @info @sprintf """  |  ____| |    |_   _|  ...  |  __  |      _____      """
        @info @sprintf """  | |    | |      | | |   .   | |  | |     (     )     """
        @info @sprintf """  | |    | |      | | | |   | | |__| |    (       )    """
        @info @sprintf """  | |____| |____ _| |_| |   | | |  | |   (         )   """
        @info @sprintf """  | _____|______|_____|_|   |_|_|  |_|  (___________)  """
        @info @sprintf """                                                       """
        @info @sprintf """ ------------------------------------------------------"""
        @info @sprintf """ Squall line                                           """
        @info @sprintf """   Resolution:                                         """
        @info @sprintf """     (Δx, Δy, Δz)   = (%.2e, %.2e, %.2e)               """ Δx Δy Δz
        @info @sprintf """     (Nex, Ney, Nez) = (%d, %d, %d)                    """ Nex Ney Nez
        @info @sprintf """     DoF = %d                                          """ DoF
        @info @sprintf """     Minimum necessary memory to run this test: %g GBs """ (DoFstorage * sizeof(DFloat))/1000^3
        @info @sprintf """     Time step dt: %.2e                                """ dt
        @info @sprintf """     End time  t : %.2e                                """ timeend
        @info @sprintf """ ------------------------------------------------------"""
    end
@info " HHH"
    engf_eng0 = run(mpicomm, dim, numelem[1:dim], polynomialorder, timeend,
                    DFloat, dt)
@info " GGG"
    show(to)
end

isinteractive() || MPI.Finalize()

nothing
