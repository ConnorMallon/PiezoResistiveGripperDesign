using Gridap
using GridapTopOpt

using Gridap, Gridap.MultiField, Gridap.Algebra, GridapSolvers

using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.TensorValues

n = 60
# Numerical Parameters
order=1
domain = (0, 0.15, 0, 0.07)
partition = (2*n, n) # MAKE SURE EITHER NUMBER IS GREATER THAN n 

model = CartesianDiscreteModel(domain, partition)

f_Γ_d1(x) = (0 <= x[1] <= 0.01 + eps()) && (x[2] ≈ 0)
f_Γ_d2(x) = (0.035 <= x[1] <= 0.045 + eps()) && (x[2] ≈ 0)
f_Γ_e1(x) = 0.075 <= x[1] <= 0.085 + eps() && (x[2] ≈ 0.0)
f_Γ_e2(x) = (0.12 <= x[1] <= 0.13 + eps()) && (x[2] ≈ 0.07)
f_Γ_f1(x) = (0.14 <= x[1] <= 0.15 + eps()) && (x[2] ≈ 0.07)
f_Γ_D1(x) = (0.0 <= x[1] <= 0.085 + eps()) && (-eps() < x[2] < 0.005)
f_Γ_D2(x) = (0.14 <= x[1] <= 0.15 + eps()) && (0.065 < x[2] < 0.07 + eps())
f_Γ_D1s(x) = (0 <= x[1] <= 0.085 + eps()) && (x[2] ≈ 0.0)
f_Γ_D2s(x) = (0.14 <= x[1] <= 0.15 + eps()) && (x[2] ≈ 0.07)
f_Γ_topright(x) = (x[1] ≈ 0.15) && (x[2] ≈ 0.07)

update_labels!(1, model, f_Γ_d1, "Gamma_d1")
update_labels!(3, model, f_Γ_e1, "Gamma_e1")
update_labels!(5, model, f_Γ_f1, "Gamma_f1")
writevtk(model,"2model")
