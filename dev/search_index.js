var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#CLIMA-1",
    "page": "Home",
    "title": "CLIMA",
    "category": "section",
    "text": "Climate Machine"
},

{
    "location": "Utilities/RootSolvers/#",
    "page": "RootSolvers",
    "title": "RootSolvers",
    "category": "page",
    "text": ""
},

{
    "location": "Utilities/RootSolvers/#CLIMA.RootSolvers",
    "page": "RootSolvers",
    "title": "CLIMA.RootSolvers",
    "category": "module",
    "text": "RootSolvers\n\nModule containing functions for solving roots of non-linear equations. The returned result is a tuple of the root and a Bool indicating convergence.\n\nfind_zero(f::F,\n           x_0::T,\n           x_1::T,\n           args::Tuple,\n           iter_params::IterParams{R, Int},\n           method::RootSolvingMethod\n           )::Tuple{T, Bool} where {F, R, T <: Union{R, AbstractArray{R}}}\n\n\n\nInterface\n\nfind_zero compute x^* such that f(x^*) = 0\n\nArguments\n\nf equation roots function, where f is callable via f(x, args...)\nx_0, x_1 initial guesses\nIterParams struct containing absolute tolerance on f(x^*) and maximum iterations\nRootSolvingMethod Algorithm to solve roots of the equation:\nSecantMethod Secant method\nRegulaFalsiMethod Regula Falsi Method\n\nSingle example\n\njulia> using RootSolvers\nx_0 = 0.0\nx_1 = 1.0\nf(x, y) = x^2 - y\nx_star2 = 10000.0\nargs = Tuple(x_star2)\nx_star = sqrt(x_star2)\ntol_abs = 1.0e-3\niter_max = 100\n\nx_root, converged = find_zero(f,\n                              x_0,\n                              x_1,\n                              args,\n                              IterParams(tol_abs, iter_max),\n                              SecantMethod())\n\nBroadcast example\n\nTo broadcast, wrap arguments that need not be broadcasted in a Ref.\n\njulia> using RootSolvers\nx_0 = rand(5,5)\nx_1 = rand(5,5).+2.0\nf(x, y) = x^2 - y\nx_star2 = 10000.0\nargs = Tuple(x_star2)\nx_star = sqrt(x_star2)\ntol_abs = 1.0e-3\niter_max = 100\nx_root, converged = find_zero.(f,\n                               x_0,\n                               x_1,\n                               Ref(args),\n                               Ref(IterParams(tol_abs, iter_max)),\n                               Ref(SecantMethod()))\n\n\n\n\n\n"
},

{
    "location": "Utilities/RootSolvers/#CLIMA.RootSolvers.find_zero",
    "page": "RootSolvers",
    "title": "CLIMA.RootSolvers.find_zero",
    "category": "function",
    "text": "find_zero(f, x_0, x_1, args, iter_params, SecantMethod)\n\nSolves the root equation, f, using Secant method. See RootSolvers for more information.\n\n\n\n\n\nfind_zero(f, x_0, x_1, args, iter_params, RegulaFalsiMethod)\n\nSolves the root equation, f, using Regula Falsi method. See RootSolvers for more information.\n\n\n\n\n\n"
},

{
    "location": "Utilities/RootSolvers/#RootSolvers-1",
    "page": "RootSolvers",
    "title": "RootSolvers",
    "category": "section",
    "text": "CurrentModule = CLIMA.RootSolversRootSolvers\nfind_zero"
},

{
    "location": "Utilities/MoistThermodynamics/#",
    "page": "MoistThermodynamics",
    "title": "MoistThermodynamics",
    "category": "page",
    "text": ""
},

{
    "location": "Utilities/MoistThermodynamics/#MoistThermodynamics-Module-1",
    "page": "MoistThermodynamics",
    "title": "MoistThermodynamics Module",
    "category": "section",
    "text": "The MoistThermodynamics module provides all thermodynamic functions needed for the atmosphere and functions shared across model components. The functions are general for a moist atmosphere that includes suspended cloud condensate in the working fluid; the special case of a dry atmosphere is obtained for zero specific humidities (or simply by omitting the optional specific humidity arguments in the functions that are needed for a dry atmosphere). The general formulation assumes that there are tracers for the total water specific humidity q_t, the liquid specific humidity q_l, and the ice specific humidity q_i to characterize the thermodynamic state and composition of moist air.There are several types of functions:Equation of state (ideal gas law):\nair_pressure\nSpecific gas constant and isobaric and isochoric specific heats of moist air:\ngas_constant_air\ncp_m\ncv_m\nSpecific latent heats of vaporization, fusion, and sublimation:\nlatent_heat_vapor\nlatent_heat_fusion\nlatent_heat_sublim\nSaturation vapor pressure and specific humidity over liquid and ice:\nsat_vapor_press_liquid\nsat_vapor_press_ice\nsat_shum\nFunctions computing energies and inverting them to obtain temperatures\ntotal_energy\ninternal_energy\nair_temperature\nFunctions to compute temperatures and partitioning of water into phases in thermodynamic equilibrium (when Gibbs\' phase rule implies that the entire thermodynamic state of moist air, including the liquid and ice specific humidities, can be calculated from the 3 thermodynamic state variables, such as energy, pressure, and total specific humidity)\nliquid_fraction (fraction of condensate that is liquid)\nsaturation_adjustment (compute temperature from energy, density, and total specific humidity)\nAuxiliary functions for diagnostic purposes, e.g., other thermodynamic quantities\nliquid_ice_pottemp (liquid-ice potential temperature)A moist dynamical core that assumes equilibrium thermodynamics can be obtained from a dry dynamical core with total energy as a prognostic variable by including a tracer for the total specific humidity q_t, using the functions, e.g., for the energies in the module, and computing the temperature T and the liquid and ice specific humidities q_l and q_i from the internal energy E_int by saturation adjustment:    T = saturation_adjustment(E_int, ρ, q_t, T_init);\n    phase_partitioning_eq!(q_l, q_i, T, ρ, q_t);here, ρ is the density of the moist air, T_init is an initial temperature guess for the saturation adjustment iterations, and the internal energy E_int = E_tot - KE - geopotential is the total energy E_tot minus kinetic energy KE and potential energy geopotential (all energies per unit mass). No changes to the \"right-hand sides\" of the dynamical equations are needed for a moist dynamical core that supports clouds, as long as they do not precipitate. Additional source-sink terms arise from precipitation.Schematically, the workflow in such a core would look as follows:\n    # initialize\n    geopotential = grav * z\n    T_prev       = ...\n    q_t          = ...\n    ρ            = ...\n\n    (u, v, w)    = ...\n    KE           = 0.5 * (u.^2 .+ v.^2 .+ w.^2)\n\n    E_tot        = total_energy(KE, geopotential, T, q_t)\n\n    do timestep   # timestepping loop\n\n      # advance dynamical variables by a timestep (temperature typically\n      # appears in terms on the rhs, such as radiative transfer)\n      advance(u, v, w, ρ, E_tot, q_t)  \n\n      # compute internal energy from dynamic variables\n      E_int = E_tot - 0.5 * (u.^2 .+ v.^2 .+ w.^2) - geopotential\n\n      # compute temperature, pressure and condensate specific humidities,\n      # using T_prev as initial condition for iterations\n      T = saturation_adjustment(E_int, ρ, q_t, T_prev);\n      phase_partitioning_eq!(q_l, q_i, T, ρ, q_t);\n      p = air_pressure(T, ρ, q_t, q_l, q_i)\n\n      # update temperature for next timestep\n      T_prev = T;  \n    endFor a dynamical core that additionally uses the liquid and ice specific humidities q_l and q_i as prognostic variables, and thus explicitly allows the presence of non-equilibrium phases such as supercooled water, the saturation adjustment in the above workflow is replaced by a direct calculation of temperature and pressure:    T = air_temperature(E_int, q_t, q_l, q_i)\n    p = air_pressure(T, ρ, q_t, q_l, q_i)"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.gas_constant_air",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.gas_constant_air",
    "category": "function",
    "text": "gas_constant_air([q_t=0, q_l=0, q_i=0])\n\nReturn the specific gas constant of moist air, given the total specific humidity q_t, and, optionally, the liquid specific humidity q_l, and the ice specific humidity q_i. When no input argument is given, it returns the specific gas constant of dry air.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.air_pressure",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.air_pressure",
    "category": "function",
    "text": "air_pressure(T, ρ[, q_t=0, q_l=0, q_i=0])\n\nReturn the air pressure from the equation of state (ideal gas law), given the air temperature T, the (moist-)air density ρ, and, optionally, the total specific humidity q_t, the liquid specific humidity q_l, and the ice specific humidity q_i. Without the specific humidity arguments, it returns the air pressure from the equation of state of dry air.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.air_density",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.air_density",
    "category": "function",
    "text": "air_density(T, p[, q_t=0, q_l=0, q_i=0])\n\nReturn the (moist-)air density from the equation of state (ideal gas law), given the air temperature T, the pressure p, and, optionally, the total specific humidity q_t, the liquid specific humidity q_l, and the ice specific humidity q_i. Without the specific humidity arguments, it returns the (moist-)air density from the equation of state of dry air.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.cp_m",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.cp_m",
    "category": "function",
    "text": "cp_m([q_t=0, q_l=0, q_i=0])\n\nReturn the isobaric specific heat capacity of moist air, given the total water specific humidity q_t, liquid specific humidity q_l, and ice specific humidity q_i. Without the specific humidity arguments, it returns the isobaric specific heat capacity of dry air.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.cv_m",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.cv_m",
    "category": "function",
    "text": "cv_m([q_t=0, q_l=0, q_i=0])\n\nReturn the isochoric specific heat capacity of moist air, given the total water specific humidity q_t, liquid specific humidity q_l, and ice specific humidity q_i. Without the specific humidity arguments, it returns the isochoric specific heat capacity of dry air.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.air_temperature",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.air_temperature",
    "category": "function",
    "text": "air_temperature(e_int[, q_t=0, q_l=0, q_i=0])\n\nReturn the air temperature, given the internal energy e_int per unit mass, and, optionally, the total specific humidity q_t, the liquid specific humidity q_l, and the ice specific humidity q_i.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.internal_energy",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.internal_energy",
    "category": "function",
    "text": "internal_energy(T[, q_t=0, q_l=0, q_i=0])\n\nReturn the internal energy per unit mass, given the temperature T, and, optionally, the total specific humidity q_t, the liquid specific humidity q_l, and the ice specific humidity q_i.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.internal_energy_sat",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.internal_energy_sat",
    "category": "function",
    "text": "internal_energy_sat(T, ρ, q_t)\n\nReturn the internal energy per unit mass in thermodynamic equilibrium at saturation, given the temperature T, (moist-)air density ρ, and total specific humidity q_t.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.total_energy",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.total_energy",
    "category": "function",
    "text": "total_energy(e_kin, e_pot, T[, q_t=0, q_l=0, q_i=0])\n\nReturn the total energy per unit mass, given the kinetic energy per unit mass e_kin, the potential energy per unit mass e_pot, the temperature T, and, optionally, the total specific humidity q_t, the liquid specific humidity q_l, and the ice specific humidity q_i.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.latent_heat_vapor",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.latent_heat_vapor",
    "category": "function",
    "text": "latent_heat_vapor(T)\n\nReturn the specific latent heat of vaporization at temperature T.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.latent_heat_sublim",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.latent_heat_sublim",
    "category": "function",
    "text": "latent_heat_sublim(T)\n\nReturn the specific latent heat of sublimation at temperature T.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.latent_heat_fusion",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.latent_heat_fusion",
    "category": "function",
    "text": "latent_heat_fusion(T)\n\nReturn the specific latent heat of fusion at temperature T.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.latent_heat_generic",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.latent_heat_generic",
    "category": "function",
    "text": "latent_heat_generic(T, LH_0, cp_diff)\n\nReturn the specific latent heat of a generic phase transition between two phases using Kirchhoff\'s relation.\n\nThe latent heat computation assumes constant isobaric specifc heat capacities of the two phases. T is the temperature, LH_0 is the latent heat of the phase transition at T_0, and cp_diff is the difference between the isobaric specific heat capacities (heat capacity in the higher-temperature phase minus that in the lower-temperature phase).\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.saturation_vapor_pressure",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.saturation_vapor_pressure",
    "category": "function",
    "text": "`saturation_vapor_pressure(T, Liquid())`\n\nReturn the saturation vapor pressure over a plane liquid surface at temperature T.\n\n`saturation_vapor_pressure(T, Ice())`\n\nReturn the saturation vapor pressure over a plane ice surface at temperature T.\n\n`saturation_vapor_pressure(T, LH_0, cp_diff)`\n\nCompute the saturation vapor pressure over a plane surface by integration of the Clausius-Clepeyron relation.\n\nThe Clausius-Clapeyron relation\n\ndlog(p_vs)/dT = [LH_0 + cp_diff * (T-T_0)]/(R_v*T^2)\n\nis integrated from the triple point temperature T_triple, using Kirchhoff\'s relation\n\nL = LH_0 + cp_diff * (T - T_0)\n\nfor the specific latent heat L with constant isobaric specific heats of the phases. The linear dependence of the specific latent heat on temperature T allows analytic integration of the Clausius-Clapeyron relation to obtain the saturation vapor pressure p_vs as a function of the triple point pressure press_triple.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.saturation_shum_generic",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.saturation_shum_generic",
    "category": "function",
    "text": "saturation_shum_generic(T, ρ[; phase=Liquid()])\n\nCompute the saturation specific humidity over a plane surface of condensate, given the temperature T and the (moist-)air density ρ.\n\nThe optional argument phase can be Liquid() or ice and indicates the condensed phase.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.saturation_shum",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.saturation_shum",
    "category": "function",
    "text": "saturation_shum(T, ρ[, q_l=0, q_i=0])\n\nCompute the saturation specific humidity, given the temperature T and (moist-)air density ρ.\n\nIf the optional liquid, and ice specific humdities q_t and q_l are given, the saturation specific humidity is that over a mixture of liquid and ice, computed in a thermodynamically consistent way from the weighted sum of the latent heats of the respective phase transitions (Pressel et al., JAMES, 2015). That is, the saturation vapor pressure and from it the saturation specific humidity are computed from a weighted mean of the latent heats of vaporization and sublimation, with the weights given by the fractions of condensate q_l/(q_l + q_i) and q_i/(q_l + q_i) that are liquid and ice, respectively.\n\nIf the condensate specific humidities q_l and q_i are not given or are both zero, the saturation specific humidity is that over a mixture of liquid and ice, with the fraction of liquid given by temperature dependent liquid_fraction(T) and the fraction of ice by the complement 1 - liquid_fraction(T).\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.saturation_shum_from_pressure",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.saturation_shum_from_pressure",
    "category": "function",
    "text": "saturation_shum_from_pressure(T, ρ, p_vs)\n\nCompute the saturation specific humidity, given the ambient air temperature T, density ρ, and the saturation vapor pressure p_vs.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.liquid_fraction",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.liquid_fraction",
    "category": "function",
    "text": "liquid_fraction(T[, q_l=0, q_i=0])\n\nReturn the fraction of condensate that is liquid.\n\nIf the optional input arguments q_l and q_i are not given or are zero, the fraction of liquid is a function that is 1 above T_freeze and goes to zero below T_freeze. If q_l or q_i are nonzero, the liquid fraction is computed from them.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.heaviside",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.heaviside",
    "category": "function",
    "text": "heaviside(t)\n\nReturn the Heaviside step function at t.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.phase_partitioning_eq!",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.phase_partitioning_eq!",
    "category": "function",
    "text": "phase_partitioning_eq!(q_l, q_i, T, ρ, q_t)\n\nReturn the partitioning of the phases in equilibrium.\n\nGiven the temperature T and (moist-)air density ρ, phase_partitioning_eq! partitions the total specific humidity q_t into the liquid specific humidity q_l and ice specific humiditiy q_l using the liquid_fraction function. The residual q_t - q_l - q_i is the vapor specific humidity.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.saturation_adjustment",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.saturation_adjustment",
    "category": "function",
    "text": "saturation_adjustment(e_int, ρ, q_t[, T_init = T_triple])\n\nReturn the temperature that is consistent with the internal energy e_int, (moist-)air density ρ, and total specific humidity q_t.\n\nThe optional input value of the temperature T_init is taken as the initial value of the saturation adjustment iterations.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#CLIMA.MoistThermodynamics.liquid_ice_pottemp",
    "page": "MoistThermodynamics",
    "title": "CLIMA.MoistThermodynamics.liquid_ice_pottemp",
    "category": "function",
    "text": "liquid_ice_pottemp(T, p[, q_t=0, q_l=0, q_i=0])\n\nReturn the liquid-ice potential temperature, given the temperature T, pressure p, total specific humidity q_t, liquid specific humidity q_l, and ice specific humidity q_i.\n\n\n\n\n\n"
},

{
    "location": "Utilities/MoistThermodynamics/#Functions-1",
    "page": "MoistThermodynamics",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = CLIMA.MoistThermodynamicsgas_constant_air\nair_pressure\nair_density\ncp_m\ncv_m\nair_temperature\ninternal_energy\ninternal_energy_sat\ntotal_energy\nlatent_heat_vapor\nlatent_heat_sublim\nlatent_heat_fusion\nlatent_heat_generic\nsaturation_vapor_pressure\nsaturation_shum_generic\nsaturation_shum\nsaturation_shum_from_pressure\nliquid_fraction\nheaviside\nphase_partitioning_eq!\nsaturation_adjustment\nliquid_ice_pottemp"
},

{
    "location": "ODESolvers/#",
    "page": "ODESolvers",
    "title": "ODESolvers",
    "category": "page",
    "text": ""
},

{
    "location": "ODESolvers/#ODESolvers-1",
    "page": "ODESolvers",
    "title": "ODESolvers",
    "category": "section",
    "text": ""
},

{
    "location": "ODESolvers/#LowStorageRungeKutta-1",
    "page": "ODESolvers",
    "title": "LowStorageRungeKutta",
    "category": "section",
    "text": "LowStorageRungeKuttaMethod.LowStorageRungeKutta\nLowStorageRungeKuttaMethod.updatedt!"
},

{
    "location": "ODESolvers/#GenericCallbacks-1",
    "page": "ODESolvers",
    "title": "GenericCallbacks",
    "category": "section",
    "text": "GenericCallbacks.EveryXWallTimeSeconds\nGenericCallbacks.EveryXSimulationSteps"
},

{
    "location": "ODESolvers/#ODESolvers-2",
    "page": "ODESolvers",
    "title": "ODESolvers",
    "category": "section",
    "text": "ODESolvers.solve!"
},

{
    "location": "Mesh/#",
    "page": "Mesh",
    "title": "Mesh",
    "category": "page",
    "text": ""
},

{
    "location": "Mesh/#Meshing-Stuff-1",
    "page": "Mesh",
    "title": "Meshing Stuff",
    "category": "section",
    "text": ""
},

{
    "location": "Mesh/#Topologies-1",
    "page": "Mesh",
    "title": "Topologies",
    "category": "section",
    "text": "Topologies encode the connectivity of the elements, spatial domain interval and MPI communication.Topologies.BrickTopology\nTopologies.StackedBrickTopology\nTopologies.CubedShellTopology\nTopologies.cubedshellmesh\nTopologies.cubedshellwarp\nTopologies.StackedCubedSphereTopology"
},

{
    "location": "Mesh/#Grids-1",
    "page": "Mesh",
    "title": "Grids",
    "category": "section",
    "text": "Grids specify the approximation within each element, and any necessary warping.Grids.DiscontinuousSpectralElementGrid"
},

{
    "location": "AtmosDycore/#",
    "page": "AtmosDycore",
    "title": "AtmosDycore",
    "category": "page",
    "text": ""
},

{
    "location": "AtmosDycore/#CLIMA.CLIMAAtmosDycore.getrhsfunction",
    "page": "AtmosDycore",
    "title": "CLIMA.CLIMAAtmosDycore.getrhsfunction",
    "category": "function",
    "text": "getrhsfunction(disc::AbstractAtmosDiscretization)\n\nThe spatial discretizations are of the form Q = f(Q), and this function returns the handle to right-hand side function f of the disc\n\n\n\n\n\n"
},

{
    "location": "AtmosDycore/#CLIMAAtmosDycore-1",
    "page": "AtmosDycore",
    "title": "CLIMAAtmosDycore",
    "category": "section",
    "text": "CurrentModule = CLIMA.CLIMAAtmosDycoregetrhsfunction"
},

{
    "location": "AtmosDycore/#CLIMA.CLIMAAtmosDycore.VanillaAtmosDiscretizations.VanillaAtmosDiscretization",
    "page": "AtmosDycore",
    "title": "CLIMA.CLIMAAtmosDycore.VanillaAtmosDiscretizations.VanillaAtmosDiscretization",
    "category": "type",
    "text": "VanillaAtmosDiscretization{nmoist, ntrace}(grid; gravity = true,\nviscosity = 0)\n\nGiven a \'grid <: AbstractGrid\' this construct all the data necessary to run a vanilla discontinuous Galerkin discretization of the the compressible Euler equations with nmoist moisture variables and ntrace tracer variables. If the boolean keyword argument gravity is true then gravity is used otherwise it is not. Isotropic viscosity can be used if viscosity is set to a positive constant.\n\n\n\n\n\n"
},

{
    "location": "AtmosDycore/#CLIMA.CLIMAAtmosDycore.VanillaAtmosDiscretizations.estimatedt",
    "page": "AtmosDycore",
    "title": "CLIMA.CLIMAAtmosDycore.VanillaAtmosDiscretizations.estimatedt",
    "category": "function",
    "text": "estimatedt(disc::VanillaAtmosDiscretization, Q::MPIStateArray)\n\nGiven a discretization disc and a state Q compute an estimate for the time step\n\ntodo: Todo\nThis estimate is currently very conservative, needs to be revisited\n\n\n\n\n\n"
},

{
    "location": "AtmosDycore/#VanillaAtmosDiscretizations-1",
    "page": "AtmosDycore",
    "title": "VanillaAtmosDiscretizations",
    "category": "section",
    "text": "A discretization adds additional information for the atmosphere problem.VanillaAtmosDiscretizations.VanillaAtmosDiscretization\nVanillaAtmosDiscretizations.estimatedt"
},

{
    "location": "AtmosDycore/#AtmosStateArray-1",
    "page": "AtmosDycore",
    "title": "AtmosStateArray",
    "category": "section",
    "text": "Storage for the state of a discretization.AtmosStateArrays.AtmosStateArray\nAtmosStateArrays.postrecvs!\nAtmosStateArrays.startexchange!\nAtmosStateArrays.finishexchange!"
},

{
    "location": "CodingConventions/#",
    "page": "Coding Conventions",
    "title": "Coding Conventions",
    "category": "page",
    "text": ""
},

{
    "location": "CodingConventions/#Coding-Conventions-1",
    "page": "Coding Conventions",
    "title": "Coding Conventions",
    "category": "section",
    "text": "A list of recommended coding conventions.There are good recommendations in the Julia style-guide:\nhttps://docs.julialang.org/en/v1/manual/style-guide/index.html\nhttps://docs.julialang.org/en/v0.6/manual/packages/#Guidelines-for-naming-a-package-1\nPlease only use Unicode characters that are within our list of acceptable Unicode characters (in AcceptableUnicode.md).\nModules, and class names (structs), should follow TitleCase convention. Note that class names cannot coincide with module names.\nFunction names should be lowercase, with words separated by underscores as necessary to improve readability.\nVariable names follow the same convention as function names. Follow CMIP conventions (http://clipc-services.ceda.ac.uk/dreq/) where possible and practicable.\nMake names consistent, distinctive, and meaningful.\nDocument design and purpose, rather than mechanics and implementation (document interfaces and embed documentation in code).\nAvoid variable names that coincide with module and class names, as well as function/variable names that are natively supported.\nNever use the characters \'l\' (lowercase letter el), \'O\' (uppercase letter oh), or \'I\' (uppercase letter eye) as single character variable names.\nTwo white spaces are used for indent. This is not part of the standard convention, but recent development efforts have been using this indentation style (e.g., Google\'s Tensorflow), and this style is being used here also.\nKISS (keep it simple stupid).\nTry to limit all lines to a maximum of 79 characters.\nSingle access point - if a variable/constant is defined more than once, then move it into a module and import (or \"using\") to that module to access the variable in order to enforce a single access point (to avoid consistency issues). Any time a chunk of code is used more than once, or when several similar versions exist across the codebase, consider generalizing this functionality and using a new function to avoid replicating code\n\"import\"/\"using\" should be grouped in the following order:\nStandard library imports.\nRelated third party imports.\nLocal application/library specific imports.\nUse a blank line between each group of imports."
},

{
    "location": "CodingConventions/#Why-do-we-limit-our-Unicode-use?-1",
    "page": "Coding Conventions",
    "title": "Why do we limit our Unicode use?",
    "category": "section",
    "text": "Some characters are visibly indistinguishable. Capital \"a\" and capital alpha are visibly indistinguishable, but are recognized as separate characters (e.g., search distinguishable).\nSome characters are difficult to read. Sometimes, the overline/overdot/hats overlap with characters making them difficult to see.\nPortability issues. Unicode does not render in Jupyter notebook natively (on OSX).\nIf it does improve readability enough, and are not worried about portability, we may introduce a list of permissible characters that are commonly used."
},

{
    "location": "AcceptableUnicode/#",
    "page": "Acceptable Unicode characters:",
    "title": "Acceptable Unicode characters:",
    "category": "page",
    "text": ""
},

{
    "location": "AcceptableUnicode/#Acceptable-Unicode-characters:-1",
    "page": "Acceptable Unicode characters:",
    "title": "Acceptable Unicode characters:",
    "category": "section",
    "text": "Using Unicode seems to be irresistible. However, we must ensure avoiding problematic Unicode usage.Below is a list of acceptable Unicode characters. All characters not listed below are forbidden. We forbid the use of accents (dot, hat, vec, etc.), because this can lead to visually ambiguous characters."
},

{
    "location": "AcceptableUnicode/#Acceptable-lower-case-Greek-letters:-1",
    "page": "Acceptable Unicode characters:",
    "title": "Acceptable lower-case Greek letters:",
    "category": "section",
    "text": "α # (alpha)\nβ # (beta)\nδ # (delta)\nϵ # (epsilon)\nε # (varepsilon)\nγ # (gamma)\nκ # (kappa)\nλ # (lambda)\nμ # (mu)\nν # (nu)\nη # (eta)\nω # (omega)\nπ # (pi)\nρ # (rho)\nσ # (sigma)\nθ # (theta)\nχ # (chi)\nξ # (xi)\nζ # (zeta)\nϕ # (psi)\nφ # (varphi)"
},

{
    "location": "AcceptableUnicode/#Acceptable-upper-case-Greek-letters:-1",
    "page": "Acceptable Unicode characters:",
    "title": "Acceptable upper-case Greek letters:",
    "category": "section",
    "text": "Δ # (Delta)\n∑ # (Sigma)\nΓ # (Gamma)\nΩ # (Omega)\nΨ # (Psi)\n<!-- Φ # (Phi) removed in favor of lowercase psi -->"
},

{
    "location": "AcceptableUnicode/#Acceptable-mathematical-symbols:-1",
    "page": "Acceptable Unicode characters:",
    "title": "Acceptable mathematical symbols:",
    "category": "section",
    "text": "∫ # (int)\n∬ # (iint)\n∭ # (iiint)\n∞ # (infinity)\n≈ # (approx)\n∂ # (partial)\n∇ # (nabla/del), note that nabla and del are indistinguishable\n∀ # (forall)\n≥ # (greater than equal to)\n≤ # (less than equal to)\n<!-- ∈ # (in) removed in favor of epsilon -->"
},

]}
