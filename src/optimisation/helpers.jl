# Degree of freedoms weights matrices computation
function get_dofs_weights_matrix(uh, points)
  Ω = get_triangulation(uh)
  Ug = get_fe_space(uh)
  cell_basis = get_data(get_fe_basis(Ug.space))
  cache = _point_to_cell_cache(KDTreeSearch(), Ω)
  x_to_cell(x) = _point_to_cell!(cache, x)
  point_to_cell = map(x_to_cell, points)
  cell_to_points, _ = make_inverse_table(point_to_cell, num_cells(Ω))
  cell_to_xs = lazy_map(Broadcasting(Reindex(points)), cell_to_points)
  cell_map = get_cell_map(Ω)
  inv_map = lazy_map(inverse_map, cell_map)
  cell_xs_at_ref_space = lazy_map(evaluate, inv_map, cell_to_xs)
  cell_to_weights = lazy_map(evaluate, cell_basis, cell_xs_at_ref_space)
  assemble_dof_weights_matrix(length(points), cell_to_weights, cell_to_points, Ug)
end

# Degree of freedoms weights matrices assembly
function assemble_dof_weights_matrix(npoints, cell_to_weights, cell_to_points, Ug)
  cell_to_dof_ids = get_cell_dof_ids(Ug)
  idxdofs, idxdofspoints, dofentries = Int64[], Int64[], eltype(eltype(cell_to_weights))[]
  idxdirs, idxdirspoints, direntries = Int64[], Int64[], eltype(eltype(cell_to_weights))[]

  function assemble_cell_wisely(weights, ipoints, dof_ids)
    if isempty(weights)
      return
    end

    npts, ndofs = size(weights)
    for ipt in 1:npts
      for idof in 1:ndofs
        if dof_ids[idof] > 0
          push!(idxdofs, dof_ids[idof])
          push!(idxdofspoints, ipoints[ipt])
          push!(dofentries, weights[ipt, idof])
        else
          push!(idxdirs, -dof_ids[idof])
          push!(idxdirspoints, ipoints[ipt])
          push!(direntries, weights[ipt, idof])
        end
      end
    end
  end
  map(assemble_cell_wisely, cell_to_weights, cell_to_points, cell_to_dof_ids)

  ndofs, ndirs = num_free_dofs(Ug), num_dirichlet_dofs(Ug)
  FW = sparse(idxdofspoints, idxdofs, dofentries, npoints, ndofs)
  DW = sparse(idxdirspoints, idxdirs, direntries, npoints, ndirs)
  return FW, DW
end

# FE model construction
function init_model(order=2, file=joinpath(@__DIR__, "square.msh"))
  reffe = ReferenceFE(lagrangian, Float64, order)
  model = GmshDiscreteModel(file)
  V0 = TestFESpace(model, reffe, dirichlet_tags=["dirichlet_points", "dirichlet_surfaces"])
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 2 * order)
  Γ = BoundaryTriangulation(Ω, tags=["neumann_surface"])
  dΓ = Measure(Γ, 2 * order)
  n_Γ = get_normal_vector(Γ)
  dv = get_fe_basis(V0)
  return model, V0, Ω, dΩ, Γ, dΓ, n_Γ, dv
end

function interpolate_param_function(func_of_p_and_x, param, Ph)
  func_of_x(x) = func_of_p_and_x(param, x)
  free_values = Vector{eltype(param)}(undef, num_free_dofs(Ph))
  diri_values = Vector{eltype(param)}(undef, num_dirichlet_dofs(Ph))
  interpolate_everywhere!(func_of_x, free_values, diri_values, Ph)
  return free_values
end

# shortcut builders
function SPIWM(j, U, V, Uφ, Ureg)
  assem_U = SparseMatrixAssembler(U, V)
  assem_deriv = SparseMatrixAssembler(Ureg, Ureg)
  objective = StateParamMap(j, U, Uφ, Ureg, assem_U, assem_deriv)
end

function AFESM(a, l, U, V, U_φ_, U_reg)
  φh0 = zero(U_φ_)
  AffineFEStateMap(a, l, U, V, U_φ_, U_reg, φh0)
end

function get_initial_φh(θ, spaces)
  ξ_ls = θ["ξ_ls"]
  a_ls = θ["a_ls"]
  b_ls = θ["b_ls"]
  φ0 = initial_lsf(ξ_ls, a_ls; b=b_ls)
  #φ0h = interpolate(interpolate(φ0,spaces["U_φ_"]),spaces["V_φ"])

  # signed “distance”‐like functions for each rectangle
  fD1 = x -> -min(x[1] - 0.0,
    0.085 - x[1],
    x[2] - 0.0,
    0.0081 - x[2])

  fD2 = x -> -min(x[1] - 0.12,
    0.15 - x[1],
    x[2] - 0.0059, #0.065-0.001,
    0.07 - x[2])


  f1(x) = min(fD1(x), fD2(x))*100
  f2(x) = min(f1(x), φ0(x))
  # combine so φ<0 on D1∪D2, φ>0 elsewhere
  #f = x -> minimum(( fD1(x), fD2(x)  ))#, φ0(x) ))

  #φ0h = interpolate(f1,spaces["V_φ"])
  φ0h = interpolate(f2, spaces["V_φ"])
  Ω = get_triangulation(spaces["V_φ"])
  #writevtk(Ω,"φ0qqq",cellfields=["φ0"=>φ0h])
  # dadad
  φ0h

end

function potato_level_set(x)

  a = 0.075 / 4
  b = 0.035 / 2
  A = 0.04 * 40
  B = 0.03 * 10
  k1 = 40
  k2 = 10
  phi1 = 0.3
  phi2 = 0.9

  term1 = ((x[1] - 0.02) / a)^2
  term2 = ((x[2] - 0.015) / b)^2
  bump_x = A * cos(k1 * x[1] + phi1)
  bump_y = B * cos(k2 * x[2] + phi2)

  return -1 * (term1 + term2 + bump_x + bump_y - 1)
end

function potato_level_set2(x)

  a = 0.075 / 3.5
  b = 0.035 / 2.5
  A = 0.04 * 20
  B = 0.03 * 10
  k1 = 20
  k2 = 20
  phi1 = 0.3
  phi2 = 0.9

  term1 = ((x[1] - 0.015) / a)^2
  term2 = ((x[2] - 0.012) / b)^2
  bump_x = A * cos(k1 * x[1] + phi1)
  bump_y = B * cos(k2 * x[2] + phi2)

  return -1 * (term1 + term2 + bump_x + bump_y - 1)
end

function get_actual_o(θ, measures, spaces, ls_evo)
  dΩ = measures["dΩ"]
  V_φ = spaces["V_φ"]
  order = θ["order"]
  γ_reinit = θ["γ_reinit"]
  shape = θ["shape"]

  δ_square = θ["δ_square"]
  δ_star8 = θ["δ_star8"]
  δ_circle_actual = θ["δ_circle"]
  δ_triangle = θ["δ_triangle"]

  if shape == "potato"
    object = potato_level_set
  elseif shape == "circle"
    δ = δ_circle_actual
    R = 0.019
    object = circle
  elseif shape == "square"
    δ = δ_square #0.007#0.02#0.015
    R = 0.013#0.025 
    object = square
  elseif shape == "triangle"
    δ = δ_triangle
    R = 0.023
    object = triangle
  elseif shape == "star6"
    δ = 0.024
    R = 0.022
    object = star6
  elseif shape == "diamond"
    δ = 0.001
    R = 0.025
    object = diamond
  elseif shape == "star8"
    δ = δ_star8#0.02
    R = 0.017
    object = star8
  else
    error("Unknown shape: $shape ")
  end

  l1 = (0.005, 0.005)
  l2 = (0.06, 0.05)
  #δ = 0 # translating along symmetry line
  o_obsh = interpolate(x -> object(R, x, δ), V_φ)

  writevtk(measures["dΩ"].quad.trian, "0", cellfields=["o_obsh" => o_obsh])
  update_collection!(spaces["Ωs_IP"], o_obsh)
  reinit!(ls_evo, o_obsh)
  writevtk(measures["dΩ"].quad.trian, "1", cellfields=["o_obsh" => o_obsh])
  o_obsh
end


function get_initial_os0h_bg(θ, measures, spaces, ls_evo)
  dΩ = measures["dΩ"]
  V_φ = spaces["V_φ"]
  order = θ["order"]
  γ_reinit = θ["γ_reinit"]
  initial_shape = θ["initial_shape"]
  @show initial_shape
  δ_circle = θ["δ_circle"]

  #o0 = initial_lsf(10,0.01) #initial_lsf(25,0.2)
  object = initial_shape

  if initial_shape == "potato"
    object = potato_level_set
  elseif initial_shape == "circle"
    δ = δ_circle #0.012 # avoid NaN at x=0 (not sure why it exists)
    R = 0.016
    object = circle
    # elseif initial_shape == "square"
    #   δ = 0.015
    #   R = 0.02 
    #   object = square
    # elseif initial_shape == "triangle"
    #   δ = 0.00
    #   R = 0.025
    #   object = triangle
    # elseif initial_shape == "star6"
    #   δ = 0.02
    #   R = 0.02
    #   object = star6
    # elseif initial_shape == "diamond"
    #   δ = 0.001
    #   R = 0.025
    #   object = diamond
  else
    error("Unknown INITIAL shape: $initial_shape")
  end

  l1 = (0.005, 0.005)
  l2 = (0.06, 0.05)
  os0h_bg = interpolate(x -> object(R, x, δ), V_φ)




  println("tmp")
  update_collection!(spaces["Ωs_IP"], os0h_bg)
  reinit!(ls_evo, os0h_bg)#,γ_reinit)


  #   reinit!(ls_evo,os0h_bg)#,γ_reinit)
  #   reinit!(ls_evo,os0h_bg)#,γ_reinit)
  #   reinit!(ls_evo,interpolate(x->x[1]-0.0001,V_φ))#,γ_reinit)
  #   println("reinitied2")

  writevtk(measures["dΩ"].quad.trian,"tmpe",cellfields=["os0h_bg"=>os0h_bg])

  #  # dada

  os0h_bg
end

function get_αf(θ, spaces)
  order = θ["order"]
  γ = θ["γ"]
  V_φ = spaces["V_φ"]
  model = get_background_model(get_triangulation(V_φ))
  el_size = model.grid.node_coords.data.partition # same as "partition"
  max_steps = floor(Int, order * minimum(el_size) / 10)
  tol = 1 / (5order^2) / minimum(el_size)
  α_coeff = 4max_steps * γ
  el_Δ = get_el_Δ(model)
  αf = α_coeff * maximum(el_Δ)
end