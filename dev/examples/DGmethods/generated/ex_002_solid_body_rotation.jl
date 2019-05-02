using MPI
using CLIMA.Topologies
using CLIMA.Grids
using CLIMA.DGBalanceLawDiscretizations
using CLIMA.MPIStateArrays
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.ODESolvers
using CLIMA.GenericCallbacks
using LinearAlgebra
using Logging
using Dates
using Printf
using StaticArrays

MPI.Initialized() || MPI.Init()

const num_aux_states = 3
function velocity_initilization!(uvec::MVector{num_aux_states, DFloat},
                                 x, y, z) where DFloat
  @inbounds begin
    r = hypot(x, y)
    θ = atan(y, x)
    uvec .= 2DFloat(π) * r .* (-sin(θ), cos(θ), 0)
  end
end

function advectionflux!(F, state, _, uvec, _)
  DFloat = eltype(state) # get the floating point type we are using
  @inbounds begin
    q = state[1]
    F[:, 1] = uvec * q
  end
end

function upwindflux!(fs, nM, stateM, viscM, uvecM, stateP, viscP, uvecP, t)
  DFloat = eltype(fs)
  @inbounds begin
    # determine the advection speed and direction
    un = dot(nM, uvecM)
    qM = stateM[1]
    qP = stateP[1]
    # Determine which state is "upwind" of the minus side
    fs[1] = un ≥ 0 ? un * qM : un * qP
  end
end

function upwindboundaryflux!(fs, nM, stateM, viscM, uvecM, stateP, viscP, uvecP,
                             bctype, t)
  DFloat = eltype(fs)
  @inbounds begin
    # determine the advection speed and direction
    un = dot(nM, uvecM)
    qM = stateM[1]
    # Determine which state is "upwind" of the minus side
    fs[1] = un ≥ 0 ? un * qM : 0
  end
end

function initialcondition!(Q, x, y, z, _)
  @inbounds Q[1] = exp(-(8 * hypot(x - 1//2, y, z))^2)
end

function exactsolution!(Q, t, x, y, z, uvec)
  @inbounds begin
    DFloat = eltype(Q)

    r = hypot(x, y)
    θ = atan(y, x) - 2DFloat(π) * t

    x, y = r * cos(θ), r * sin(θ)

    initialcondition!(Q, x, y, z, uvec)
  end
end

function setupDG(mpicomm, dim, Ne, polynomialorder, DFloat=Float64,
                 ArrayType=Array)

  @assert ArrayType === Array

  brickrange = (range(DFloat(-1); length=Ne+1, stop=1),
                range(DFloat(-1); length=Ne+1, stop=1),
                range(DFloat(-1); length=Ne+1, stop=1))

  topology = BrickTopology(mpicomm, brickrange[1:dim])

  grid = DiscontinuousSpectralElementGrid(topology; polynomialorder =
                                          polynomialorder, FloatType = DFloat,
                                          DeviceArray = ArrayType,)

  spatialdiscretization = DGBalanceLaw(grid = grid, length_state_vector = 1,
                                       flux! = advectionflux!,
                                       numerical_flux! = upwindflux!,
                                       numerical_boundary_flux! =
                                       upwindboundaryflux!,
                                       auxiliary_state_length = num_aux_states,
                                       auxiliary_state_initialization! =
                                       velocity_initilization!)

end

let
  mpicomm = MPI.COMM_WORLD
  mpi_logger = ConsoleLogger(MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull)
  dim = 2
  Ne = 20
  polynomialorder = 4
  spatialdiscretization = setupDG(mpicomm, dim, Ne, polynomialorder)
  Q = MPIStateArray(spatialdiscretization, initialcondition!)
  filename = @sprintf("initialcondition_mpirank%04d", MPI.Comm_rank(mpicomm))
  DGBalanceLawDiscretizations.writevtk(filename, Q, spatialdiscretization,
                                       ("q",))

  h = 1 / Ne

  CFL = h / (2π)
  dt = CFL / polynomialorder^2
  lsrk = LowStorageRungeKutta(spatialdiscretization, Q; dt = dt, t0 = 0)
  finaltime = 1.0

  vtk_step = 0
  mkpath("vtk")
  cb_vtk = GenericCallbacks.EveryXSimulationSteps(20) do
    vtk_step += 1
    filename = @sprintf("vtk/solid_body_rotation_mpirank%04d_step%04d",
                         MPI.Comm_rank(mpicomm), vtk_step)
    DGBalanceLawDiscretizations.writevtk(filename, Q, spatialdiscretization,
                                         ("q",))
    nothing
  end

  solve!(Q, lsrk; timeend = finaltime, callbacks = (cb_vtk, ))

  filename = @sprintf("finalsolution_mpirank%04d", MPI.Comm_rank(mpicomm))
  DGBalanceLawDiscretizations.writevtk(filename, Q, spatialdiscretization,
                                       ("q",))

  Qe = MPIStateArray(spatialdiscretization) do Qin, x, y, z, uvec
    exactsolution!(Qin, finaltime, x, y, z, uvec)
  end

  error = euclidean_distance(Q, Qe)
  with_logger(mpi_logger) do
    @info @sprintf("""Run with
                   dim              = %d
                   Ne               = %d
                   polynomial order = %d
                   error            = %e
                   """, dim, Ne, polynomialorder, error)
  end
end

let
  mpicomm = MPI.COMM_WORLD
  mpi_logger = ConsoleLogger(MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull)

  dim = 2
  polynomialorder = 4
  finaltime = 1.0

  with_logger(mpi_logger) do
    @info @sprintf("""Running with
                   dim              = %d
                   polynomial order = %d
                   """, dim, polynomialorder)
  end

  base_Ne = 5
  lvl_error = zeros(4) # number of levels to compute is length(lvl_error)
  for lvl = 1:length(lvl_error)
    # `Ne` for this mesh level
    Ne = base_Ne * 2^(lvl-1)
    spatialdiscretization = setupDG(mpicomm, dim, Ne, polynomialorder)

    Q = MPIStateArray(spatialdiscretization, initialcondition!)
    h = 1 / Ne
    CFL = h / (2π)
    dt = CFL / polynomialorder^2
    lsrk = LowStorageRungeKutta(spatialdiscretization, Q; dt = dt, t0 = 0)

    solve!(Q, lsrk; timeend = finaltime)

    Qe = MPIStateArray(spatialdiscretization) do Qin, x, y, z, uvec
      exactsolution!(Qin, finaltime, x, y, z, uvec)
    end

    lvl_error[lvl] = euclidean_distance(Q, Qe)
    msg =  @sprintf   "Level      = %d" lvl
    msg *= @sprintf "\nNe               = %d" Ne
    msg *= @sprintf "\nerror            = %.4e" lvl_error[lvl]
    if lvl > 1
      rate = log2(lvl_error[lvl-1]) - log2(lvl_error[lvl])
      msg *= @sprintf "\nconvergence rate = %.4e" rate
    end
    with_logger(mpi_logger) do
      @info msg
    end
  end
end

Sys.iswindows() || MPI.finalize_atexit()
Sys.iswindows() && !isinteractive() && MPI.Finalize()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

