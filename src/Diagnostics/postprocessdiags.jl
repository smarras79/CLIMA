using VegaLite, DataFrames
using JLD2,FileIO
using CLIMA.VariableTemplates
include("diagnostic_vars.jl")
d = load("output/diagnostics.jld2")
Nqk = size(d["0.0"],1)
Nev = size(d["0.0"],2)
FT = Float32
diagnostic_vars(dv) = Vars{vars_diagnostic(FT)}(dv)
S = zeros(Nqk*Nev)
Z = zeros(Nqk*Nev)
Denom = 0
for key in keys(d)
  for ev in 1:Nev
    for k in 1:Nqk
      dv = diagnostic_vars(d[key][k,ev])
      S[k+(ev-1)*Nqk] += dv.w′w′
      
      Z[k+(ev-1)*Nqk] = dv.z
    end
  end
  Denom += 1
end
S = S /Denom
data = DataFrame(a = S, b = Z)
data |> @vlplot(:line, x=:a, y=:b)

  
