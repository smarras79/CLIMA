"""
    Diagnostics

Accumulate mean fields and covariance statistics on the computational grid

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

include("thermo_vars.jl")
include("horzavg_vars.jl")
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

function gather_diagnostics(mpicomm, dg, Q, current_time_string, κ, out_dir)
    mpirank = MPI.Comm_rank(mpicomm)
    nranks = MPI.Comm_size(mpicomm)

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

    function state_visitor(state, i, j, k, ev, eh, accum)
    end

    function visitQ(funs::Vector{Function}, accum = nothing)
        for eh in 1:nhorzelem
            for ev in 1:nvertelem
                e = ev + (eh - 1) * nvertelem

                for k in 1:Nqk
                    for j in 1:Nq
                        for i in 1:Nq
                            ijk = i + Nq * ((j-1) + Nq * (k-1)) 
                            state = extract_state(dg, localQ, ijk, e)
                            for f in funs
                                f(state, i, j, k, ev, eh, accum)
                            end
                        end
                    end
                end
            end
        end
    end

    zvals = zeros(Nqk, nvertelem)
    thermoQ = fill(zeros(num_thermo(FT)), (npoints, nrealelem))
    thermo_vars(th) = Vars{vars_thermo(FT)}(th)

    function compute_thermo!(state, i, j, k, ev, eh, accum)
        ijk = i + Nq * ((j-1) + Nq * (k-1))
        e = ev + (eh - 1) * nvertelem

        zvals[k,ev] = localvgeo[ijk,grid.x3id,e]

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

    horzsums = fill(zeros(num_horzavg(FT)), (Nqk, nvertelem))
    horzavg_vars(ha) = Vars{vars_horzavg(FT)}(ha)
    hdvsr = Nq * Nq * nhorzelem

    LWP = zero(FT)
    l_LWP = zero(FT)

    function compute_horzsums!(state, i, j, k, ev, eh, accum)
        ijk = i + Nq * ((j-1) + Nq * (k-1))
        e = ev + (eh - 1) * nvertelem

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
            # TODO: uncomment the line below after rewriting the l_LWP assignment below using aux.∫dz...?
            # aux = extract_aux(dg, localaux, ijk, e)
            l_LWP += (localaux[ijk,1,e] + localaux[ijk,2,e]) / κ / hdvsr
        end
    end

    visitQ(Function[compute_thermo!, compute_horzsums!])

    horzavgs = fill(zeros(num_horzavg(FT)), (Nqk, nvertelem))
    for ev in 1:nvertelem
        for k in 1:Nqk
            for n in 1:num_horzavg(FT)
                horzavgs[k,ev][n] = horzsums[k,ev][n] / hdvsr
                horzavgs[k,ev][n] = MPI.Reduce(horzavgs[k,ev][n], +, 0, mpicomm)
                if mpirank == 0
                    horzavgs[k,ev][n] /= nranks
                end
            end
        end
    end

    LWP = MPI.Reduce(l_LWP, +, 0, mpicomm)
    if mpirank == 0
        LWP /= nranks
    end

    dsums = fill(zeros(num_diagnostic(FT)), (Nqk, nvertelem))
    diagnostic_vars(dv) = Vars{vars_diagnostic(FT)}(dv)

    function compute_diagnosticsums!(state, i, j, k, ev, eh, accum)
        ijk = i + Nq * ((j-1) + Nq * (k-1))
        e = ev + (eh - 1) * nvertelem

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

    visitQ(Function[compute_diagnosticsums!])

    davgs = fill(zeros(num_diagnostic(FT)), (Nqk, nvertelem))
    for ev in 1:nvertelem
        for k in 1:Nqk
            for n in 1:num_diagnostic(FT)
                davgs[k,ev][n] = dsums[k,ev][n] / hdvsr
                davgs[k,ev][n] = MPI.Reduce(davgs[k,ev][n], +, 0, mpicomm)
                if mpirank == 0
                    davgs[k,ev][n] /= nranks
                end
            end
        end
    end

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

