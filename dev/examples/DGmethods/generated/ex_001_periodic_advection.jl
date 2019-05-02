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

const uvec = (1, 2, 3)

function advectionflux!(F, state, _...)
  DFloat = eltype(state) # get the floating point type we are using
  @inbounds begin
    q = state[1]
    F[:, 1] = SVector{3, DFloat}(uvec) * q
  end
end

function upwindflux!(fs, nM, stateM, viscM, auxM, stateP, viscP, auxP, t)
  DFloat = eltype(fs)
  @inbounds begin
    # determine the advection speed and direction
    un = dot(nM, DFloat.(uvec))
    qM = stateM[1]
    qP = stateP[1]
    # Determine which state is "upwind" of the minus side
    fs[1] = un ≥ 0 ? un * qM : un * qP
  end
end

function initialcondition!(Q, x_1, x_2, x_3)
  @inbounds Q[1] = exp(sin(2π * x_1)) * exp(sin(2π * x_2)) * exp(sin(2π * x_3))
end

function exactsolution!(dim, Q, t, x_1, x_2, x_3)
  @inbounds begin
    DFloat = eltype(Q)

    # trace back the point (x_1, x_2, x_3) in the velocity field and
    # determine where in our "original" [0, L_1] X [0, L_2] X [0, L_3] domain
    # this point is located
    y_1 = mod(x_1 - DFloat(uvec[1]) * t, 1)
    y_2 = mod(x_2 - DFloat(uvec[2]) * t, 1)

    # if we are really just 2-D we do not want to change the x_3 coordinate
    y_3 = dim == 3 ? mod(x_3 - DFloat(uvec[3]) * t, 1) : x_3

    initialcondition!(Q, y_1, y_2, y_3)
  end
end

function setupDG(mpicomm, dim, Ne, polynomialorder, DFloat=Float64,
                 ArrayType=Array)

  @assert ArrayType === Array

  brickrange = (range(DFloat(0); length=Ne+1, stop=1), # x_1 corner locations
                range(DFloat(0); length=Ne+1, stop=1), # x_2 corner locations
                range(DFloat(0); length=Ne+1, stop=1)) # x_3 corner locations

  periodicity = (true, true, true)

  topology = BrickTopology(mpicomm, brickrange[1:dim];
                           periodicity=periodicity[1:dim])

  grid = DiscontinuousSpectralElementGrid(topology; polynomialorder =
                                          polynomialorder, FloatType = DFloat,
                                          DeviceArray = ArrayType,)

  spatialdiscretization = DGBalanceLaw(grid = grid, length_state_vector = 1,
                                       flux! = advectionflux!,
                                       numerical_flux! = upwindflux!)

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

  h = 1 / Ne                           # element size
  CFL = h / maximum(abs.(uvec[1:dim])) # time to cross the element once
  dt = CFL / polynomialorder^2         # DG time step scaling (for this
                                       # particular RK scheme could go with a
                                       # factor of ~2 larger time step)
  lsrk = LowStorageRungeKutta(spatialdiscretization, Q; dt = dt, t0 = 0)

  finaltime = 1.0
  solve!(Q, lsrk; timeend = finaltime)

  filename = @sprintf("finalsolution_mpirank%04d", MPI.Comm_rank(mpicomm))
  DGBalanceLawDiscretizations.writevtk(filename, Q, spatialdiscretization,
                                       ("q",))

  Qe = MPIStateArray(spatialdiscretization) do Qin, x, y, z
    exactsolution!(dim, Qin, finaltime, x, y, z)
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
  Ne = 20
  polynomialorder = 4
  spatialdiscretization = setupDG(mpicomm, dim, Ne, polynomialorder)
  Q = MPIStateArray(spatialdiscretization, initialcondition!)
  filename = @sprintf("initialcondition_mpirank%04d", MPI.Comm_rank(mpicomm))
  DGBalanceLawDiscretizations.writevtk(filename, Q, spatialdiscretization,
                                       ("q",))
  h = 1 / Ne
  CFL = h / maximum(abs.(uvec[1:dim]))
  dt = CFL / polynomialorder^2
  lsrk = LowStorageRungeKutta(spatialdiscretization, Q; dt = dt, t0 = 0)
  finaltime = 1.0

  store_norm_index = 0
  normQ = Array{Float64}(undef, ceil(Int, finaltime / dt))
  function cb_store_norm()
    store_norm_index += 1
    normQ[store_norm_index] = norm(Q)
    nothing
  end

  vtk_step = 0
  mkpath("vtk")
  cb_vtk = GenericCallbacks.EveryXSimulationSteps(20) do
    vtk_step += 1
    filename = @sprintf("vtk/advection_mpirank%04d_step%04d",
                         MPI.Comm_rank(mpicomm), vtk_step)
    DGBalanceLawDiscretizations.writevtk(filename, Q, spatialdiscretization,
                                         ("q",))
    nothing
  end

  starttime = Ref(now())
  cb_info = GenericCallbacks.EveryXWallTimeSeconds(1, mpicomm) do (init=false)
    if init
      starttime[] = now()
    else
      with_logger(mpi_logger) do
        @info @sprintf("""Update
                       simtime = %.16e
                       runtime = %s
                       norm(Q) = %.16e""", ODESolvers.gettime(lsrk),
                       Dates.format(convert(Dates.DateTime,
                                            Dates.now()-starttime[]),
                                    Dates.dateformat"HH:MM:SS"),
                       norm(Q))
      end
    end
  end

  solve!(Q, lsrk; timeend = finaltime,
         callbacks = (cb_store_norm, cb_vtk, cb_info))

  filename = @sprintf("finalsolution_mpirank%04d", MPI.Comm_rank(mpicomm))
  DGBalanceLawDiscretizations.writevtk(filename, Q, spatialdiscretization,
                                       ("q",))

  Qe = MPIStateArray(spatialdiscretization) do Qin, x, y, z
    exactsolution!(dim, Qin, finaltime, x, y, z)
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
    CFL = h / maximum(abs.(uvec[1:dim]))
    dt = CFL / polynomialorder^2
    lsrk = LowStorageRungeKutta(spatialdiscretization, Q; dt = dt, t0 = 0)

    solve!(Q, lsrk; timeend = finaltime)

    Qe = MPIStateArray(spatialdiscretization) do Qin, x, y, z
      exactsolution!(dim, Qin, finaltime, x, y, z)
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

