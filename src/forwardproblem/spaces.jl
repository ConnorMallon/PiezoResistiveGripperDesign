function attach_top_right_constraint(Vd,measures)
 
  Γin = measures[:dΓin].quad.trian

  cell_inids = (get_cell_dof_ids(V,Γin))
  cell_outids = (get_cell_dof_ids(V,Γout))

  inids = getindex.(cell_inids,1) # getting bottom left dof of linear simplex on left edge
  outids = getindex.(cell_outids,3) # getting top right dof of linear simplex on right edge
  
  pushfirst!(outids,cell_outids[1][1]) # adding bottom right element
  pop!(outids) # get rid top right
  
  fdof_to_val = collect(Float64,1:num_free_dofs(V))
  ddof_to_val = -collect(Float64,1:num_dirichlet_dofs(V))
  vh = FEFunction(V,fdof_to_val,ddof_to_val)
  get_cell_dof_values(vh)

  sDOF_to_dof = outids
  sDOF_to_dofs   = Table([[id] for id in inids])
  sDOF_to_coeffs = Table([[1 ] for id in inids])

  Vc = FESpaceWithLinearConstraints(
    sDOF_to_dof,
    sDOF_to_dofs,
    sDOF_to_coeffs,
    V)
end

function get_spaces(θ,measures)

order = θ["order"]
order_d = 1 
println("using first order elements for d")
Sϕ = θ["Sϕ"]

Ω = measures["dΩ"].quad.trian
Γ_c = measures["dΓ_c"].quad.trian
D = 2 

# Boundary Conditions
d_D(x) = VectorValue(0,0)
Vin(x) = 1.0
Vout(x) = 0.0

## Spaces
reffe_scalar = ReferenceFE(lagrangian,Float64,order)
reffed = ReferenceFE(lagrangian,VectorValue{D,Float64},order_d)
reffeϕ = reffe_scalar

V_φ = TestFESpace(Ω,reffe_scalar)
V_φ_IP = TestFESpace(Ω, reffe_scalar)

V_φ_ = TestFESpace(Ω, reffe_scalar;
  dirichlet_tags=[
      "Gamma_D2_boundary_vertical","Gamma_D1_boundary_vertical","Gamma_D1_vertical","Gamma_D2_vertical", #"Gamma_D1_gaps",#"Gamma_D2_gaps",
      "AllMyPoints", #"Gamma_D1_extra",
  "Gamma_D1","Gamma_D1s","Gamma_d1", "Gamma_d2", "Gamma_e1","Gamma_e3","Gamma_e4","Gamma_e5","Gamma_e6","Gamma_D2","Gamma_D2s" ,"Gamma_e2", "Gamma_f1"],#,"Gamma_O1"],
  )
#f_φ(x) = - x[1]*10
U_φ_ = TrialFESpace(V_φ_,-1)

#p0 =  interpolate(1,U_φ_)
#writevtk(Ω,"p0",cellfields=["p0"=>p0])


T = eltype(ReverseDiff.GradientConfig(reffeϕ[2][1][]).input)
reffeϕT = ReferenceFE(lagrangian, T, order)
V_diff = TestFESpace(Ω, reffeϕT; vector_type=Vector{T})

V_reg = V_φ_
U_reg = TrialFESpace(V_reg, 0.0)

V_reg_IP = TestFESpace(Ω, reffe_scalar;
  dirichlet_tags=["Gamma_f1"],
  )
U_reg_IP = TrialFESpace(V_reg_IP, 0.0)

Vk = TestFESpace(Γ_c,reffe_scalar)

Vd_(Ω) = TestFESpace(Ω, reffed;conformity=:H1,
  dirichlet_tags="Gamma_e1",
  )
Γ_f1 = measures["dΓ_f1"].quad.trian
# ids = get_cell_dof_ids(Vd_,Γ_f1) 

# polys = get_polytopes(get_background_model(Ω))  
# @assert length(polys) == 1 "Only one cell type is currently supported"
# poly = first(polys)
# if poly == TRI
#   u1_tr = last(ids)[3]
#   u2_tr = last(ids)[6]
# elseif poly == QUAD
#   u1_tr = last(ids)[4]
#   u2_tr = last(ids)[8]
# end

# sDOF_to_dof = [u1_tr]
# sDOF_to_dofs = Table([[u2_tr]])
# sDOF_to_coeffs = Table([[1]])
# Vd = FESpaceWithLinearConstraints(
#   sDOF_to_dof,
#   sDOF_to_dofs,
#   sDOF_to_coeffs,
#   Vd_)

Vd(Ω) = Vd_(Ω)
Ud(Ω) = TrialFESpace(Vd(Ω),VectorValue(0.0,0.0))

num_comp = D * (D + 1) ÷ 2 # Number of components of a symmetric tensor in D-dim space
sym_tensor_type = Gridap.Fields.SymTensorValue{D,Float64,num_comp}
reffes = ReferenceFE(lagrangian,sym_tensor_type,order)
Vs = TestFESpace(Ω,reffes;conformity=:H1)
Us = TrialFESpace(Vs)

matrix_type = Gridap.Fields.SymTensorValue{D,Float64}
order_dg0 = 1
reffe_matrix = ReferenceFE(lagrangian,matrix_type,order_dg0)
VIk(Ω) = TestFESpace(Ω,reffe_matrix;conformity=:H1)
UIk(Ω) = TrialFESpace(VIk(Ω))

function dirichlet_tags_ϕ(i) 
  if i==1
    return ["Gamma_e2","Gamma_f1"]
  elseif i==2
    return ["Gamma_e4","Gamma_e2"]
  elseif i==3 
    return ["Gamma_e6","Gamma_e4"]
  elseif i==4
    return ["Gamma_e3","Gamma_e6"]
  elseif i==5 
    return ["Gamma_e5","Gamma_e3"]
  elseif i==6
    return ["Gamma_f1","Gamma_e5"]
  elseif i==7
    return ["Gamma_e2","Gamma_e3"]
  end
end

# function dirichlet_tags_ϕ(i) 
#   #@assert n_ϕ == 8 #10
#   if i == 1
#     return ["Gamma_f1","Gamma_e2"]
#   elseif i == 2 
#     return ["Gamma_e2","Gamma_e5"] # changed
#   elseif i == 3
#     return ["Gamma_d1","Gamma_e5"]
#   elseif i == 4
#     return ["Gamma_d1","Gamma_d2"]
#   elseif i == 5
#     return ["Gamma_e1","Gamma_f1"]
#   elseif i == 6
#     return ["Gamma_e5","Gamma_e3"]
#   elseif i == 7
#     return ["Gamma_e3","Gamma_e6"]
#   elseif i == 8 
#     return ["Gamma_e4","Gamma_e6"]
#   # elseif i == 9 
#   #   return ["Gamma_e4","Gamma_e1"]
#   # elseif i == 10
#   #   return ["Gamma_d2","Gamma_e1"]
#   end
# end

Vϕ(Ω,i) = TestFESpace(Ω,reffeϕ;conformity=:H1,dirichlet_tags=dirichlet_tags_ϕ(i))
Uϕ(Ω,i) = TrialFESpace(Vϕ(Ω,i), [Vin,Vout])

# #mfs = BlockMultiFieldStyle()
# X =  MultiFieldFESpace([Ud,Ud,Ud])#;style=mfs)
# Y =  MultiFieldFESpace([Vd,Vd,Vd])#;style=mfs)

model = get_background_model(Ω)
Ωs = EmbeddedCollection(model,zero(V_φ)) do cutgeo,_,_
  Ωin = DifferentiableTriangulation(Triangulation(cutgeo,PHYSICAL),V_φ)
  Γ = DifferentiableTriangulation(EmbeddedBoundary(cutgeo),V_φ)
  Γg = GhostSkeleton(cutgeo)
  Ωact = Triangulation(cutgeo,ACTIVE)
  chi_pairs = [
    Symbol("χ$(i)") => GridapTopOpt.get_isolated_volumes_mask(cutgeo, dirichlet_tags_ϕ(i))
    for i ∈ Sϕ
  ]
  chi_nt = NamedTuple(chi_pairs)
  base_nt = 
  (;
    :Ωin  => Ωin,
    :dΩin => Measure(Ωin,2*order),
    :Γg   => Γg,
    :dΓg  => Measure(Γg,2*order),
    :n_Γg => get_normal_vector(Γg),
    :Γ    => Γ,
    :dΓ   => Measure(Γ,2*order),
    :n_Γ  => get_normal_vector(Γ), # Note, need to recompute inside obj/constraints to compute derivs
    :Ωact => Ωact,
    :χd => GridapTopOpt.get_isolated_volumes_mask(cutgeo, "Gamma_e1"),
  )
  return merge(base_nt, chi_nt)
end

Ωs_IP = EmbeddedCollection(model,zero(V_φ)) do cutgeo,_,_
  Γ =  DifferentiableTriangulation(EmbeddedBoundary(cutgeo),V_φ)
  (;
  :dΓ   => Measure(Γ,2*order),
  )
  end

spaces = Dict(
  "V_φ"=>V_φ,
  "V_φ_IP"=> V_φ_IP,
  "V_φ_"=>V_φ_,
  "U_φ_"=>U_φ_,
  "V_diff"=>V_diff,
  "V_reg"=>V_reg,
  "U_reg"=>U_reg,
  "V_reg_IP"=>V_reg_IP,
  "U_reg_IP"=>U_reg_IP,
  "Vd"=>Vd,
  "Ud"=>Ud,
  "Vϕ"=>Vϕ,
  "Uϕ"=>Uϕ,
  "Vs"=>Vs,
  "Us"=>Us,
  "VIk"=>VIk,
  "UIk"=>UIk,
  # "X"=>X,
  # "Y"=>Y,
  "Vk" => Vk,
  "Ωs" => Ωs,
  "Ωs_IP" => Ωs_IP,
  )

spaces 

end
