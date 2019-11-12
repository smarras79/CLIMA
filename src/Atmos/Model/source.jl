using CLIMA.PlanetParameters: Omega
export Gravity, RayleighSponge, Subsidence, GeostrophicForcing, Coriolis

# kept for compatibility
# can be removed if no functions are using this
function atmos_source!(f::Function, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  f(source, state, aux, t)
end
function atmos_source!(::Nothing, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
end
# sources are applied additively
function atmos_source!(stuple::Tuple, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  map(s -> atmos_source!(s, m, source, state, aux, t), stuple)
end

abstract type Source
end

struct Gravity <: Source
end
function atmos_source!(::Gravity, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  source.ρu -= state.ρ * aux.orientation.∇Φ
end

struct Subsidence <: Source
end
function atmos_source!(::Subsidence, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  source.ρu -= state.ρ * m.radiation.D_subsidence
end

struct Coriolis <: Source
end
function atmos_source!(::Coriolis, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  # note: this assumes a SphericalOrientation
  source.ρu -= SVector(0, 0, 2*Omega) × state.ρu
end

struct GeostrophicForcing{FT} <: Source
  f_coriolis::FT
  u_geostrophic::FT
  v_geostrophic::FT
end
function atmos_source!(s::GeostrophicForcing, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  u = state.ρu / state.ρ
  u_geo = SVector(s.u_geostrophic, s.v_geostrophic, 0)
  source.ρu -= state.ρ * s.f_coriolis * (u - u_geo)
end


"""
  RayleighSponge{FT} <: Source
Rayleigh Damping (Relaxation) for top wall momentum components
Assumes laterally periodic boundary conditions for LES flows. Momentum components
are relaxed to reference values (user defined) within the sponge layer whose thickness
is user defined in the driver.
The relaxation function now is a simple quartic function. Additional sponges may be added.
"""
struct RayleighSponge{FT} <: Source
  "Domain maximum height [m]"
  zmax::FT
  "Vertical extent at with sponge starts [m]"
  zsponge::FT
  "Sponge Strength 0 ⩽ c_sponge ⩽ 1"
  c_sponge::FT
  "Lateral reference relaxation velocity, u [m/s]"
  u_ref::FT
  "Lateral reference relaxation velocity, v [m/s]"
  v_ref::FT
  "Vertical reference relaxation velocity, w [m/s]"
  w_ref::FT
    
end

function atmos_source!(s::RayleighSponge, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  FT = eltype(state)
  z = aux.orientation.Φ / grav
  coeff = FT(0)
  if z >= s.zsponge
    coeff_top = s.c_sponge * (sinpi(FT(1/2)*(z - s.zsponge)/(s.zmax-s.zsponge)))^FT(4)
    coeff = min(coeff_top, 1.0)
  end

  u = state.ρu / state.ρ
  u_relaxation = SVector(s.u_ref, s.v_ref, s.w_ref)
  source.ρu -= state.ρ * coeff * (u - u_relaxation)
end
