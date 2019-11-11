"""
    Diagnostics variable template

Container for diagnostic variables of interest. Useful also for 
post-processing.

"""

using CLIMA.VariableTemplates

function vars_diagnostic(FT)
  @vars begin
    # vertical coordinate
    z::FT
    # state and functions of state
    u::FT
    v::FT
    w::FT
    q_tot::FT
    e_tot::FT
    q_liq::FT
    θ::FT
    θ_liq::FT
    θ_v::FT
    e_int::FT
    # vertical fluxes
    w′ρ′::FT
    w′u′::FT
    w′v′::FT
    w′q_tot′::FT
    w′q_liq′::FT
    w′q_vap′::FT
    w′θ′::FT
    w′θ_v′::FT
    w′θ_liq′::FT
    # variances
    u′u′::FT
    v′v′::FT
    w′w′::FT
    # skewness
    w′w′w′::FT
    # turbulent kinetic energy
    TKE::FT
  end
end
num_diagnostic(FT) = varsize(vars_diagnostic(FT))
diagnostic_vars(array) = Vars{vars_diagnostic(eltype(array))}(array)

