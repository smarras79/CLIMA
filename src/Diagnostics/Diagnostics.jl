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
    e_int::FT
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
    th.e_int = e_int
end

# horizontal averages
function vars_horzavg(FT)
  @vars begin
    ρ::FT
    ρu::FT
    ρv::FT
    ρw::FT
    e_tot::FT
    ρq_tot::FT
    q_liq::FT
    q_vap::FT
    θ_liq_ice::FT
    θ_dry::FT
    θ_v::FT
    e_int::FT
  end
end
num_horzavg(FT) = varsize(vars_horzavg(FT))
horzavg_vars(array) = Vars{vars_horzavg(eltype(array))}(array)

function compute_horzsums!(FT, state, i, j, k, ijk, ev, eh, e,
                           Nqk, nvertelem, localaux, κ, LWP,
                           thermoQ, horzsums, localvgeo, xmax, ymax, repdvsr, x1id, x2id, Nq)
    #the next lines are used to properly take into account repeating nodes while also factoring in nodes on the boundary
    x = localvgeo[ijk,x1id,e]
    y = localvgeo[ijk,x2id,e]

    if ((x == 0 || abs(x - xmax)<=0.001) && (ymax == 0 || abs(y - 1500)<=0.001))
      bound = 4 #node on the corner of domain thus is not repeated and we cancel out the 4 times repetition below
    elseif (x == 0 || abs(x - xmax)<= 0.001 || ymax == 0 || abs(y - 1500)<=0.001)
      bound = 2 #node on the edge of the domain thus we halve it's repetition
    else
      bound = 1 #node not on any considered boundary and thus nothing done.
    end

    if ((i == 1 || i == Nq) && (j ==1 || j==Nq))
      rep = 1/4 * bound #corner node repeated 4 times for horizontal considerations
    elseif (i == 1 || i == Nq || j==1 || j==Nq)
      rep = 1/2 * bound #edge node repeated 2 times for horizontal considerations
    else
      rep = 1 * bound # inner node nothing done
    end
    #end of repetition consideration

    th = thermo_vars(thermoQ[ijk,e])
    hs = horzavg_vars(horzsums[k,ev])
    hs.ρ         += rep * state.ρ
    hs.ρu        += rep * state.ρu[1]
    hs.ρv        += rep * state.ρu[2]
    hs.ρw        += rep * state.ρu[3]
    hs.e_tot     += rep * state.ρe
    hs.ρq_tot    += rep * state.moisture.ρq_tot
    hs.q_liq     += rep * th.q_liq
    hs.q_vap     += rep * th.q_vap
    hs.θ_liq_ice += rep * th.θ_liq_ice
    hs.θ_dry     += rep * th.θ_dry
    hs.θ_v       += rep * th.θ_v
    hs.e_int     += rep * th.e_int

    # liquid water path
    # This condition is also going to be used to get the number of points that exist on a horizontal plane provided all planes have the same number of points
    # TODO adjust for possibility of non equivalent horizontal slabs
    if ev == floor(nvertelem/2) && k == floor(Nqk/2)
        # TODO: uncomment the line below after rewriting the LWP assignment below using aux.∫dz...?
        # aux = extract_aux(dg, localaux, ijk, e)
        LWP[1] += rep * (localaux[ijk,1,e] + localaux[ijk,2,e]) / κ 
        repdvsr[1] += rep #number of points to be divided by
    end
end

function compute_diagnosticsums!(FT, state, i, j, k, ijk, ev, eh, e,
                                 zvals, thermoQ, horzavgs, dsums, localvgeo, xmax, ymax, x1id, x2id, Nq)

    #the next lines are used to properly take into account repeating nodes while also factoring in nodes on the boundary
    x = localvgeo[ijk,x1id,e]
    y = localvgeo[ijk,x2id,e]

    if ((x == 0 || abs(x - xmax)<=0.001) && (ymax == 0 || abs(y - 1500)<=0.001))
      bound = 4 #node on the corner of domain thus is not repeated and we cancel out the 4 times repetition below
    elseif (x == 0 || abs(x - xmax)<= 0.001 || ymax == 0 || abs(y - 1500)<=0.001)
      bound = 2 #node on the edge of the domain thus we halve it's repetition
    else
      bound = 1 #node not on any considered boundary and thus nothing done.
    end

    if ((i == 1 || i == Nq) && (j ==1 || j==Nq))
      rep = 1/4 * bound #corner node repeated 4 times for horizontal considerations
    elseif (i == 1 || i == Nq || j==1 || j==Nq)
      rep = 1/2 * bound #edge node repeated 2 times for horizontal considerations
    else
      rep = 1 * bound # inner node nothing done
    end
    #end of repetition consideration

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
    ẽ = ha.e_tot / ha.ρ
    q̃_tot = ha.ρq_tot / ha.ρ

    # vertical coordinate
    ds.z         += zvals[k,ev]

    # state and functions of state
    ds.u        += rep * ũ
    ds.v        += rep * ṽ
    ds.w        += rep * w̃
    ds.e_tot    += rep * ẽ
    ds.q_tot    += rep * ha.ρq_tot / ha.ρ
    ds.q_liq    += rep * ha.q_liq
    ds.θ        += rep * ha.θ_dry
    ds.θ_liq    += rep * ha.θ_liq_ice
    ds.θ_v      += rep * ha.θ_v
    ds.e_int    += rep * ha.e_int

    # vertical fluxes
    ds.w′ρ′     += rep * (w̅ - w̃) * (state.ρ - ha.ρ)
    ds.w′u′     += rep * (w̅ - w̃) * (u̅ - ha.ρu / ha.ρ)
    ds.w′v′     += rep * (w̅ - w̃) * (v̅ - ha.ρv / ha.ρ)
    ds.w′q_tot′ += rep * (w̅ - w̃) * (q̅_tot - q̃_tot)
    ds.w′q_liq′ += rep * (w̅ - w̃) * (th.q_liq - ha.q_liq)
    ds.w′q_vap′ += rep * (w̅ - w̃) * (th.q_vap - ha.q_vap)
    ds.w′θ′     += rep * (w̅ - w̃) * (th.θ_dry - ha.θ_dry)
    ds.w′θ_v′   += rep * (w̅ - w̃) * (th.θ_v - ha.θ_v)
    ds.w′θ_liq′ += rep * (w̅ - w̃) * (th.θ_liq_ice - ha.θ_liq_ice)

    # variances
    ds.u′u′     += rep * (u̅ - ũ)^2
    ds.v′v′     += rep * (v̅ - ṽ)^2
    ds.w′w′     += rep * (w̅ - w̃)^2

    # skewness
    ds.w′w′w′   += rep * (w̅ - w̃)^3

    # turbulent kinetic energy
    ds.TKE      += rep * 0.5 * (ds.u′u′ + ds.v′v′ + ds.w′w′)
end

# TODO: make this a single reduction
function horz_average_all(FT, mpicomm, num, (Nqk, nvertelem), sums, repdvsr)
    mpirank = MPI.Comm_rank(mpicomm)
    nranks = MPI.Comm_size(mpicomm)
    avgs = [zeros(FT, num) for _ in 1:Nqk, _ in 1:nvertelem]
    for ev in 1:nvertelem
        for k in 1:Nqk
            for n in 1:num
                avgs[k,ev][n] = sums[k,ev][n]
                avgs[k,ev][n] = MPI.Reduce(avgs[k,ev][n], +, 0, mpicomm)
                if mpirank == 0
                    avgs[k,ev][n] /= repdvsr
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
function gather_diagnostics(mpicomm, dg, Q, current_time_string, κ, xmax, ymax ,out_dir)
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
    l_repdvsr = zeros(FT, 1) 

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
                          Nqk, nvertelem, localaux, κ, l_LWP,
                          thermoQ, horzsums, localvgeo, xmax, ymax, l_repdvsr, grid.x1id, grid.x2id, Nq)

    # run both in one grid traversal
    visitQ(FT, Function[thermo_visitor, horzsum_visitor])
    
    # compute the full number of points on a slab
    repdvsr = zero(FT)
    repdvsr = MPI.Reduce(l_repdvsr[1], +, 0, mpicomm)
    # compute the horizontal and LWP averages
    horzavgs = horz_average_all(FT, mpicomm, num_horzavg(FT), (Nqk, nvertelem),
                                horzsums, repdvsr)
    LWP = zero(FT)
    LWP = MPI.Reduce(l_LWP[1], +, 0, mpicomm)
    if mpirank == 0
        LWP /= repdvsr
    end

    # compute the diagnostics with the previous computed variables
    dsums = [zeros(FT, num_diagnostic(FT)) for _ in 1:Nqk, _ in 1:nvertelem]
    dsum_visitor(FT, state, i, j, k, ijk, ev, eh, e) =
        compute_diagnosticsums!(FT, state, i, j, k, ijk, ev, eh, e,
                                zvals, thermoQ, horzavgs, dsums, localvgeo, xmax, ymax, grid.x1id, grid.x2id, Nq)

    # another grid traversal
    visitQ(FT, Function[dsum_visitor])

    # compute the averages
    davgs = horz_average_all(FT, mpicomm, num_diagnostic(FT), (Nqk, nvertelem),
                             dsums, repdvsr)

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

