"""
    Diagnostics

Accumulate mean fields and covariance statistics on the computational grid.

"""

# TODO: reduce MPI communication; record local averages only?

module Diagnostics

using FileIO
using JLD2
using MPI
using StaticArrays

using ..Atmos
using ..DGmethods: num_state, vars_state, num_aux, vars_aux
using ..Mesh.Topologies
using ..Mesh.Grids
using ..MoistThermodynamics
using ..MPIStateArrays
using ..PlanetParameters
using ..VariableTemplates

export gather_diagnostics

include("diagnostic_vars.jl")

function extract_state(dg, localQ, ijk, e)
    bl = dg.balancelaw
    FT = eltype(localQ)
    nstate = num_state(bl, FT)
    l_Q = MArray{Tuple{nstate},FT}(undef)
    for s in 1:nstate
        l_Q[s] = localQ[ijk,s,e]
    end
    return Vars{vars_state(bl, FT)}(l_Q)
end

function extract_aux(dg, auxstate, ijk, e)
    bl = dg.balancelaw
    FT = eltype(auxstate)
    nauxstate = num_aux(bl, FT)
    l_aux = MArray{Tuple{nauxstate},FT}(undef)
    for s in 1:naux
        l_aux[s] = auxstate[ijk,s,e]
    end
    return Vars{vars_aux(bl, FT)}(l_aux)
end

# thermodynamic variables of interest
function vars_thermo(FT)
  @vars begin
    q_liq::FT
    q_ice::FT
    q_vap::FT
    T::FT
    θ_liq_ice::FT
    θ_dry::FT
    θ_v::FT
  end
end
num_thermo(FT) = varsize(vars_thermo(FT))
thermo_vars(array) = Vars{vars_thermo(eltype(array))}(array)

function compute_thermo!(FT, state, i, j, k, ijk, ev, eh, e,
                         localvgeo, x3id, zvals, thermoQ)
    zvals[k,ev] = localvgeo[ijk,x3id,e]

    u̅ = state.ρu[1] / state.ρ
    v̅ = state.ρu[2] / state.ρ
    w̅ = state.ρu[3] / state.ρ
    e̅_tot = state.ρe / state.ρ
    q̅_tot = state.moisture.ρq_tot / state.ρ

    e_int = e̅_tot - 1//2 * (u̅^2 + v̅^2 + w̅^2) - grav * zvals[k,ev]

    ts = PhaseEquil(convert(FT, e_int), q̅_tot, state.ρ)
    Phpart = PhasePartition(ts)

    th = thermo_vars(thermoQ[ijk,e])
    th.q_liq     = Phpart.liq
    th.q_ice     = Phpart.ice
    th.q_vap     = q̅_tot-Phpart.liq-Phpart.ice
    th.T         = ts.T
    th.θ_liq_ice = liquid_ice_pottemp(ts)
    th.θ_dry     = dry_pottemp(ts)
    th.θ_v       = virtual_pottemp(ts)
end

# horizontal averages
function vars_horzavg(FT)
  @vars begin
    ρ::FT
    ρu::FT
    ρv::FT
    ρw::FT
    ρq_tot::FT
    q_liq::FT
    q_vap::FT
    θ_liq_ice::FT
    θ_dry::FT
    θ_v::FT
  end
end
num_horzavg(FT) = varsize(vars_horzavg(FT))
horzavg_vars(array) = Vars{vars_horzavg(eltype(array))}(array)

function compute_horzsums!(FT, state, i, j, k, ijk, ev, eh, e,
                           Nqk, nvertelem, localaux, κ, hdvsr, LWP,
                           thermoQ, horzsums)
    th = thermo_vars(thermoQ[ijk,e])
    hs = horzavg_vars(horzsums[k,ev])
    hs.ρ         += state.ρ
    hs.ρu        += state.ρu[1]
    hs.ρv        += state.ρu[2]
    hs.ρw        += state.ρu[3]
    hs.ρq_tot    += state.moisture.ρq_tot
    hs.q_liq     += th.q_liq
    hs.q_vap     += th.q_vap
    hs.θ_liq_ice += th.θ_liq_ice
    hs.θ_dry     += th.θ_dry
    hs.θ_v       += th.θ_v

    # liquid water path
    if ev == floor(nvertelem/2) && k == floor(Nqk/2)
        # TODO: uncomment the line below after rewriting the LWP assignment below using aux.∫dz...?
        # aux = extract_aux(dg, localaux, ijk, e)
        LWP[1] += (localaux[ijk,1,e] + localaux[ijk,2,e]) / κ / hdvsr
    end
end

function compute_diagnosticsums!(FT, state, i, j, k, ijk, ev, eh, e,
                                 zvals, thermoQ, horzavgs, dsums)
    th = thermo_vars(thermoQ[ijk,e])
    ha = horzavg_vars(horzavgs[k,ev])
    ds = diagnostic_vars(dsums[k,ev])

    u̅ = state.ρu[1] / state.ρ
    v̅ = state.ρu[2] / state.ρ
    w̅ = state.ρu[3] / state.ρ
    q̅_tot = state.moisture.ρq_tot / state.ρ
    ũ = ha.ρu / ha.ρ
    ṽ = ha.ρv / ha.ρ
    w̃ = ha.ρw / ha.ρ
    q̃_tot = ha.ρq_tot / ha.ρ

    # vertical coordinate
    ds.z         = zvals[k,ev]

    # state and functions of state
    ds.u        += ũ
    ds.v        += ṽ
    ds.w        += w̃
    ds.q_tot    += ha.ρq_tot / ha.ρ
    ds.q_liq    += ha.q_liq
    ds.θ        += ha.θ_dry
    ds.θ_liq    += ha.θ_liq_ice

    # vertical fluxes
    ds.w′ρ′     += (w̅ - w̃) * (state.ρ - ha.ρ)
    ds.w′u′     += (w̅ - w̃) * (u̅ - ha.ρu / ha.ρ)
    ds.w′v′     += (w̅ - w̃) * (v̅ - ha.ρv / ha.ρ)
    ds.w′q_tot′ += (w̅ - w̃) * (q̅_tot - q̃_tot)
    ds.w′q_liq′ += (w̅ - w̃) * (th.q_liq - ha.q_liq)
    ds.w′q_vap′ += (w̅ - w̃) * (th.q_vap - ha.q_vap)
    ds.w′θ′     += (w̅ - w̃) * (th.θ_dry - ha.θ_dry)
    ds.w′θ_v′   += (w̅ - w̃) * (th.θ_v - ha.θ_v)
    ds.w′θ_liq′ += (w̅ - w̃) * (th.θ_liq_ice - ha.θ_liq_ice)

    # variances
    ds.u′u′     += (u̅ - ũ)^2
    ds.v′v′     += (v̅ - ṽ)^2
    ds.w′w′     += (w̅ - w̃)^2

    # skewness
    ds.w′w′w′   += (w̅ - w̃)^3

    # turbulent kinetic energy
    ds.TKE      += 0.5 * (ds.u′u′ + ds.v′v′ + ds.w′w′)
end

# TODO: make this a single reduction
function horz_average_all(FT, mpicomm, num, (Nqk, nvertelem), sums, dvsr)
    mpirank = MPI.Comm_rank(mpicomm)
    nranks = MPI.Comm_size(mpicomm)
    avgs = [zeros(FT, num) for _ in 1:Nqk, _ in 1:nvertelem]
    for ev in 1:nvertelem
        for k in 1:Nqk
            for n in 1:num
                avgs[k,ev][n] = sums[k,ev][n] / dvsr
                avgs[k,ev][n] = MPI.Reduce(avgs[k,ev][n], +, 0, mpicomm)
                if mpirank == 0
                    avgs[k,ev][n] /= nranks
                end
            end
        end
    end
    return avgs
end

"""
    gather_diagnostics(mpicomm, dg, Q, current_time_string, κ, out_dir)

Compute various diagnostic variables and write them to JLD2 files in `out_dir`,
indexed by `current_time_string`.
"""
function gather_diagnostics(mpicomm, dg, Q, current_time_string, κ, out_dir)
    mpirank = MPI.Comm_rank(mpicomm)
    nranks = MPI.Comm_size(mpicomm)

    # extract grid information
    bl = dg.balancelaw
    grid = dg.grid
    topology = grid.topology
    N = polynomialorder(grid)
    Nq = N + 1
    Nqk = dimensionality(grid) == 2 ? 1 : Nq
    npoints = Nq * Nq * Nqk
    nrealelem = length(topology.realelems)
    nvertelem = topology.stacksize
    nhorzelem = div(nrealelem, nvertelem)

    # get the state, auxiliary and geo variables onto the host if needed
    if Array ∈ typeof(Q).parameters
        localQ = Q.realdata
        localaux = dg.auxstate.realdata
        localvgeo = grid.vgeo
    else
        localQ = Array(Q.realdata)
        localaux = Array(dg.auxstate.realdata)
        localvgeo = Array(grid.vgeo)
    end
    FT = eltype(localQ)

    nstate = num_state(bl, FT)
    nauxstate = num_aux(bl, FT)

    # divisor for horizontal averages
    hdvsr = Nq * Nq * nhorzelem

    # traverse the grid, running each of `funs` on each node
    function visitQ(FT, funs::Vector{Function})
        for eh in 1:nhorzelem
            for ev in 1:nvertelem
                e = ev + (eh - 1) * nvertelem

                for k in 1:Nqk
                    for j in 1:Nq
                        for i in 1:Nq
                            ijk = i + Nq * ((j-1) + Nq * (k-1)) 
                            state = extract_state(dg, localQ, ijk, e)
                            for f in funs
                                f(FT, state, i, j, k, ijk, ev, eh, e)
                            end
                        end
                    end
                end
            end
        end
    end

    # record the vertical coordinates and compute thermo variables
    zvals = zeros(Nqk, nvertelem)
    thermoQ = [zeros(FT, num_thermo(FT)) for _ in 1:npoints, _ in 1:nrealelem]
    thermo_visitor(FT, state, i, j, k, ijk, ev, eh, e) =
        compute_thermo!(FT, state, i, j, k, ijk, ev, eh, e, localvgeo, grid.x3id,
                        zvals, thermoQ)

    # compute the horizontal sums and the liquid water path
    l_LWP = zeros(FT, 1)
    horzsums = [zeros(FT, num_horzavg(FT)) for _ in 1:Nqk, _ in 1:nvertelem]
    horzsum_visitor(FT, state, i, j, k, ijk, ev, eh, e) =
        compute_horzsums!(FT, state, i, j, k, ijk, ev, eh, e,
                          Nqk, nvertelem, localaux, κ, hdvsr, l_LWP,
                          thermoQ, horzsums)

    # run both in one grid traversal
    visitQ(FT, Function[thermo_visitor, horzsum_visitor])

    # compute the horizontal and LWP averages
    horzavgs = horz_average_all(FT, mpicomm, num_horzavg(FT), (Nqk, nvertelem),
                                horzsums, hdvsr)
    LWP = zero(FT)
    LWP = MPI.Reduce(l_LWP[1], +, 0, mpicomm)
    if mpirank == 0
        LWP /= nranks
    end

    # compute the diagnostics with the previous computed variables
    dsums = [zeros(FT, num_diagnostic(FT)) for _ in 1:Nqk, _ in 1:nvertelem]
    dsum_visitor(FT, state, i, j, k, ijk, ev, eh, e) =
        compute_diagnosticsums!(FT, state, i, j, k, ijk, ev, eh, e,
                                zvals, thermoQ, horzavgs, dsums)

    # another grid traversal
    visitQ(FT, Function[dsum_visitor])

    # compute the averages
    davgs = horz_average_all(FT, mpicomm, num_diagnostic(FT), (Nqk, nvertelem),
                             dsums, hdvsr)

    if mpirank == 0
        jldopen(joinpath(out_dir, "diagnostics.jld2"), "a+") do file
            file[current_time_string] = davgs
        end
        jldopen(joinpath(out_dir, "liquid_water_path.jld2"), "a+") do file
            file[current_time_string] = LWP
        end
    end
end # function gather_diagnostics

end # module Diagnostics

