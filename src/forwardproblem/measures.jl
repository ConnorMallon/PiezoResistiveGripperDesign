function get_measures(θ)

# function optimise_pr(θ)
α2 = θ["α2"]
α3 = θ["α3"]
η_coeff_TO = θ["η_coeff_TO"]
η_coeff_IP = θ["η_coeff_IP"]
n = θ["n"]
order = θ["order"]

# Numerical Parameters
order=1
domain = (0, 0.15, 0, 0.07)
partition = (2*n, n) # MAKE SURE EITHER NUMBER IS GREATER THAN n 

# # model 
# model = simplexify(CartesianDiscreteModel(domain, partition))
# #el_Δ = get_el_Δ(model)
# f_Γ_d1(x) = (0 <= x[1] <= 0.01 + eps()) && (x[2] ≈ 0)
# f_Γ_d2(x) = (0.035 <= x[1] <= 0.045 + eps()) && (x[2] ≈ 0)
# f_Γ_e1(x) = 0.075 <= x[1] <= 0.085 + eps() && (x[2] ≈ 0.0)
# #f_Γ_e1(x) = (x[1] ≈ 0.085) && (x[2] ≈ 0.0)
# f_Γ_e2(x) = (0.12 <= x[1] <= 0.125 + eps()) && (x[2] ≈ 0.07)
# f_Γ_f1(x) = (0.145 - eps() <= x[1] <= 0.15 + eps()) && (x[2] ≈ 0.07)
# f_Γ_D1(x) = (0.0 <= x[1] <= 0.085 + eps()) && (-eps() < x[2] < 0.005)
# f_Γ_D2(x) = (0.12 <= x[1] <= 0.15 + eps()) && (0.065 < x[2] < 0.07 + eps())
# f_Γ_D1s(x) = (0 <= x[1] <= 0.085 + eps()) && (x[2] ≈ 0.0)
# f_Γ_D2s(x) = (0.12 <= x[1] <= 0.15 + eps()) && (x[2] ≈ 0.07)
# f_Γ_topright(x) = (x[1] ≈ 0.15) && (x[2] ≈ 0.07)
# @assert n % 15 == 0 "n must be a multiple of 15"
# f_Γ_e3(x) = (0.02-eps() <= x[1] <= 0.025 + eps()) && (x[2] ≈ 0.0046666666666667)
# f_Γ_e4(x) = (0.06-eps() <= x[1] <= 0.065 + eps()) && (x[2] ≈ 0.0046666666666667)
# f_Γ_e5(x) = (0.0 <= x[1] <= 0.005 + eps()) && (x[2] ≈ 0.0046666666666667)
# f_Γ_e6(x) = (0.04-eps() <= x[1] <= 0.045 + eps()) && (x[2] ≈ 0.0046666666666667)
# f_Ω_fk(x) = (0.14 <= x[1] <= 0.15 + eps()) && (0.065 <= x[2] <= 0.07 + eps())
# f_Ω_dk(x) = (0.075 <= x[1] <= 0.085 + eps()) && (0.0 <= x[2] <= 0.005 + eps())

# update_labels!(8, model, f_Γ_D1, "Gamma_D1")
# update_labels!(9, model, f_Γ_D2, "Gamma_D2")
# #update_labels!(14, model, f_Ω_fk, "Gamma_Ω_fk")
# #update_labels!(15, model, f_Ω_dk, "Gamma_Ω_dk")

# update_labels!(6, model, f_Γ_D1s, "Gamma_D1s")
# update_labels!(7, model, f_Γ_D2s, "Gamma_D2s")
# update_labels!(1, model, f_Γ_d1, "Gamma_d1")
# update_labels!(2, model, f_Γ_d2, "Gamma_d2")
# update_labels!(3, model, f_Γ_e1, "Gamma_e1")
# update_labels!(4, model, f_Γ_e2, "Gamma_e2")
# update_labels!(5, model, f_Γ_f1, "Gamma_f1")
# update_labels!(10, model, f_Γ_e3, "Gamma_e3")
# update_labels!(11, model, f_Γ_e4, "Gamma_e4")
# update_labels!(12, model, f_Γ_e5, "Gamma_e5")
# update_labels!(13, model, f_Γ_e6, "Gamma_e6")




model = GmshDiscreteModel(srcdir("geometry/domain_with_e.msh"))

point_coords = Dict(
  # the four corners
  "P00"    => (0.00,   0.00),
  #"P1500"  => (0.15,   0.00),
  "P15070" => (0.15,   0.07),
  #"P0070"  => (0.00,   0.07),
  # bottom‐split points
  "B001"   => (0.01,   0.00),
  "B035"   => (0.035,  0.00),
  "B045"   => (0.045,  0.00),
  "B075"   => (0.075,  0.00),
  "B085"   => (0.085,  0.00),
  # top‐split points
  "T012070"  => (0.12,  0.07),
  "T0125070" => (0.125, 0.07),
  "T0145070" => (0.145, 0.07),
  # D1 band ends
  "D1A"    => (0.00,   0.0046666667),
  #"D1B"    => (0.15,   0.0046666667),
  # D2 band ends
  #"D2A"    => (0.00,   0.065),
  "D2B"    => (0.15,   0.065),
  # vertical splits
  "V1"     => (0.085,  0.0046666667),
  "V2"     => (0.12,   0.065),
  # interior feature‐line endpoints
  "E3A"    => (0.02,   0.0046666667),
  # "E3B"    => (0.03,  0.0046666667),
  "E3B"    => (0.025,  0.0046666667),
  
  "E4A"    => (0.06,   0.0046666667),
  #"E4B"    => (0.07,  0.0046666667),
  "E4B"    => (0.065,  0.0046666667),
  
  "E5A"    => (0.00,   0.0046666667),
  # "E5B"    => (0.01,  0.0046666667),
  "E5B"    => (0.005,  0.0046666667),
  
  "E6A"    => (0.04,   0.0046666667),
  #"E6B"    => (0.05,  0.0046666667)
  "E6B"    => (0.045,  0.0046666667)
)
# 2) One test‐function that returns true if a point sits at any of those coords
f_all_pts(x) = any( x[1]≈x0 && x[2]≈y0
                    for (name, (x0,y0)) in point_coords )
# 3) One single group
update_labels!(500, model, f_all_pts, "AllMyPoints")



## Triangulations and measures
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Γ_d1 = BoundaryTriangulation(model,tags="Gamma_d1")
Γ_d2 = BoundaryTriangulation(model,tags="Gamma_d2")
Γ_e1 = BoundaryTriangulation(model,tags="Gamma_e1")
Γ_e2 = BoundaryTriangulation(model,tags="Gamma_e2")
Γ_f1 = BoundaryTriangulation(model,tags="Gamma_f1")
Γ_c = BoundaryTriangulation(model,tags=["Gamma_D1s","Gamma_d1","Gamma_d2","Gamma_e1"])
Γ_e3 = BoundaryTriangulation(model,tags="Gamma_e3")
Γ_e4 = BoundaryTriangulation(model,tags="Gamma_e4")

Ω_D2 = Triangulation(model,tags="Gamma_D2")
#Ω_fk = Triangulation(model,tags="Gamma_Ω_fk")
#Ω_dk = Triangulation(model,tags="Gamma_Ω_dk")

dΩ = Measure(Ω,2order)
dΓ = Measure(Γ_d1,2order)
dΓ_d1 = Measure(Γ_d1,2order)
dΓ_d2 = Measure(Γ_d2,2order)
dΓ_e1 = Measure(Γ_e1,2order)
dΓ_e2 = Measure(Γ_e2,2order)
dΓ_f1 = Measure(Γ_f1,2order)
dΓ_c = Measure(Γ_c,2order)
n_Γc = get_normal_vector(Γ_c)
dΓ_e3 = Measure(Γ_e3,2order)
dΓ_e4 = Measure(Γ_e4,2order)
dΓ_d1_meas = sum(∫(1)dΓ_d1)
dΓ_f1_meas = sum(∫(1)dΓ_f1)
dΩ_D2 = Measure(Ω_D2,2order)
#dΩ_fk = Measure(Ω_fk,2order)
#dΩ_dk = Measure(Ω_dk,2order)

measures = Dict(
  "dΩ"=>dΩ,
  "dΓ"=>dΓ,
  "dΓ_d1"=>dΓ_d1,
  "dΓ_d2"=>dΓ_d2,
  "dΓ_e1"=>dΓ_e1,
  "dΓ_e2"=>dΓ_e2,
  "dΓ_f1"=>dΓ_f1,
  "dΓ_c"=>dΓ_c,
  "n_Γc"=>n_Γc,
  "dΓ_e3"=>dΓ_e3,
  "dΓ_e4"=>dΓ_e4,
  "dΩ_D2"=>dΩ_D2,
#  "dΩ_fk"=>dΩ_fk,  
#  "dΩ_dk"=>dΩ_dk,
  )
end 