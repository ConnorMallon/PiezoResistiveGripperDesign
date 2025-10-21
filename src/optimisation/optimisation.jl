function optimise1(θ,measures,spaces,φh0,U_ve_reg,V_ve_reg,maxiter,functionals)
	
    order = θ["order"]
    γ = θ["γ"]
    γ_reinit = θ["γ_reinit"]
    iter_mod = θ["iter_mod"]
	reg_coeff = θ["reg_coeff"]
	reg_coeff_IP = θ["reg_coeff_IP"]
	n = θ["n"]
	ν_reinit = θ["ν_reinit"]
	max_steps_factor = θ["max_steps_factor"]
    dΩ = measures["dΩ"]
	V_φ = spaces["V_φ"]
	U_φ_ = spaces["U_φ_"]
	U_reg = spaces["U_reg"]
	V_reg = spaces["V_reg"]

	Ω = dΩ.quad.trian
	model = get_background_model(Ω)
	h = minimum(get_element_diameters(model))	
	max_steps = max_steps_factor/h # /(maximum(partition))  # Time-steps for evolution equation
    hₕ = get_element_diameter_field(model)
	evo = CutFEMEvolver(V_φ,spaces["Ωs"],measures["dΩ"],hₕ;max_steps=max_steps,γg=0.1)
	reinit = StabilisedReinitialiser(V_φ,spaces["Ωs"],measures["dΩ"],hₕ;
	stabilisation_method=ArtificialViscosity(ν_reinit))
	ls_evo = LevelSetEvolution(evo,reinit)


	evo_IP = CutFEMEvolver(V_φ,spaces["Ωs_IP"],measures["dΩ"],hₕ;max_steps,γg=0.1)
	reinit_IP = StabilisedReinitialiser(V_φ,spaces["Ωs_IP"],measures["dΩ"],hₕ;
	stabilisation_method=ArtificialViscosity(ν_reinit))
	ls_evo_IP = LevelSetEvolution(evo_IP,reinit_IP)
	
	# square(t) = t^2

	# reinit_IP = HeatReinitialiser(V_φ,model; t =  4*h , boundary_tags = [
	# 	"Gamma_d1",                # bottom, left segment
	# 	"Gamma_d2",                # bottom, middle segment
	# 	"Gamma_e1",                # bottom, right segment
	# 	"Gamma_e2",                # top, left segment
	# 	"Gamma_f1",                # top, right segment
	# 	"Gamma_D1_boundary_vertical",  # left & right vertical in the D1 band
	# 	"Gamma_D2_boundary_vertical",  # left & right vertical in the D2 band
	# 	"Gamma_remainder"          # whatever outer‑boundary pieces weren’t in any of the above
	#   ])
	# ls_evo_IP = LevelSetEvolution(evo_IP,reinit_IP)

	## Hilbertian extension-regularisation problems
	α_coeff = 4max_steps*γ * reg_coeff
	α = α_coeff * maximum(get_element_diameters(model))
	a_hilb(p,q) =∫(α^2*∇(p)⋅∇(q) + p*q)dΩ

	α_coeff_IP = 4max_steps*γ * reg_coeff_IP
	α_IP = α_coeff_IP * maximum(get_element_diameters(model))
	a_hilb_IP(p,q) =∫(α_IP^2*∇(p)⋅∇(q) + p*q)dΩ

	vel_ext = VelocityExtension(a_hilb,U_ve_reg,V_ve_reg)
	vel_ext_IP = VelocityExtension(a_hilb_IP,V_φ,V_φ)
	# update_collection!(spaces["Ωs"],φh0)
	
	# writevtk(Ω,"init0",cellfields=["φ"=>φh0])
	# writevtk(spaces["Ωs"].Ωin,"init00",cellfields=["φ"=>φh0])


	#reinit!(ls_evo,φh0)
	# writevtk(Ω,"init",cellfields=["φ"=>φh0])
	# writevtk(spaces["Ωs"].Ωin,"init1",cellfields=["φ"=>φh0])


	φh0, ls_evo, vel_ext, ls_evo_IP, vel_ext_IP
	end



	function optimise(θ,measures,spaces,pcfs,oh0,maxiter,ls_evo,vel_ext,functionals)

			
		order = θ["order"]
		γ = θ["γ"]
		γ_reinit = θ["γ_reinit"]
		iter_mod = θ["iter_mod"]
		reg_coeff = θ["reg_coeff"]
		n = θ["n"]
		
		dΩ = measures["dΩ"]
	
		V_φ = spaces["V_φ"]
		U_φ_ = spaces["U_φ_"]
		U_reg = spaces["U_reg"]
		V_reg = spaces["V_reg"]
	
		Ω = dΩ.quad.trian
		model = get_background_model(Ω)

	## Optimiser	
	get_last_iteration(h) = h.niter
	get_dof_Δ(m) = m.params.Δ
	ph0 = FEFunction(V_φ,deepcopy(oh0.free_values))





	# optimiser = HilbertianProjection(pcfs,ls_evo,vel_ext,ph0;
	# 	γ,γ_reinit,verbose=true,debug=false,
	# 	constraint_names=[Symbol("C$i") for i in 1:1],#length(pcfs.dC)],
	# 	#converged=custom_hp_converged,
	# 	maxiter = maxiter)




	@show maxiter

	function custom_al_converged(
		m::AugmentedLagrangian
	  )
		false
	  end
	
	optimiser = AugmentedLagrangian(pcfs,ls_evo,vel_ext,ph0;
		γ,verbose=true,debug=false,
		converged=custom_al_converged,

		#constraint_names=[Symbol("C$i") for i in 1:1], # length(pcfs.dC)],
		maxiter = maxiter)

	path=datadir("results")
	mkpath(path) 

	for (it,uh,φh) in optimiser
		data = ["φ"=>φh,"H(φ)"=>(functionals["H"] ∘ φh),"‖∇(φ)‖"=> norm∘(∇(φh)) ]
		#iszero(it % iter_mod) && 
		#writevtk(Ω,"outin$it",cellfields=data)	
		# iszero(it % iter_mod) && writevtk(spaces["Ωs"].dΩin.quad.trian,"out$it")
		#write_history(path*"/history.txt",optimiser.history)	
	end
	#opt_results = Dict("φh_bg"=>oh0,"state_maps"=>state_maps)
	Js = optimiser.history[:J]

	# options = Optim.Options(
	# 	iterations = 10,
	# 	store_trace = true,
	# 	show_trace = true,
	# 	show_warnings = true)

	# 	println("optim opt")
	# result = Optim.optimize(p->pcfs.φ_to_jc(p), ph0.free_values, BFGS(),options)
	# pf = Optim.minimizer(result)
	# ph0 = FEFunction(V_φ, pf)
	# Js = Optim.minimum(results)



	ph0, Js, ls_evo

end