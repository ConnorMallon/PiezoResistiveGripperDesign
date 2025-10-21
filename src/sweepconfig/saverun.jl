function saverun(θ, φh, ofh, oih, oah, dah, dih, spaces, functionals, sname, TO_pcfs,state_maps,d_to_l,j1s, j3s, c1s, ls)
	model_name = θ["model_name"]
	Sϕ = θ["Sϕ"]
	σ = functionals["σ"]
	σ_vm = functionals["σ_vm"]
	K = functionals["K"]
	K0 = functionals["K_0"]
	H = functionals["H"]
	H_IP = functionals["H_IP"]
	DH_IP = functionals["DH_IP"]
	I_ϕ = functionals["I_ϕ"]
	
	# TO results 
	Ωact = spaces["Ωs"].Ωact
	φ = φh.free_values 
	da = TO_pcfs.embedded_collection.state_map_da(φ)
	dafh  = FEFunction(spaces["Ud"](Ωact),da)
	Ikf = TO_pcfs.embedded_collection.state_map_Ik([φ;da])
	Ikh = FEFunction(spaces["UIk"](Ωact),Ikf)
	Ik0f = TO_pcfs.embedded_collection.state_map_Ik0([φ;da*0])
	Ik0h = FEFunction(spaces["UIk"](Ωact),Ik0f)
	Uϕs = [spaces["Uϕ"](Ωact,i) for i in Sϕ]
	ϕfs = [TO_pcfs.embedded_collection.state_map_ϕs[j]([φ;da]) for j in 1:length(Sϕ)]
	ϕfhs = [FEFunction(Uϕs[j],ϕfs[j]) for j in 1:length(Sϕ)]	
	ϕ0fs = [TO_pcfs.embedded_collection.state_map_ϕ0s[j]([φ;da*0]) for j in 1:length(Sϕ)]
	ϕ0fhs = [FEFunction(Uϕs[j],ϕ0fs[j]) for j in 1:length(Sϕ)]
	dϕhs_squared = [(abs2∘(ϕ - ϕ0)) for (ϕ,ϕ0) in zip(ϕfhs,ϕ0fhs)]
	#Cs = [Ikh⋅∇(ϕ) for ϕ in ϕfhs]
	#C0s = [Ik0h⋅∇(ϕ0) for ϕ0 in ϕ0fhs]
	σ_vmh = σ_vm(dafh)
	Vc = ConstantFESpace(get_background_model(get_triangulation(spaces["V_φ"])))

	data_TO = [
		"φh" => φh,"dafh" => dafh, #"ddfh" => ddfh, 
		"σ_vmh" => σ_vmh, #,"daafh" => daafh, #"safh"=>safh, #"j1h" => j1h,"j2h" => j2h,
		"H(φ)"=>(H ∘ φh),"|∇(φ)|"=>(norm ∘ ∇(φh)),"Iσ_vm"=> (1-H∘φh)*(σ_vm(dafh)), #"Is" => (1-H∘φh)*safh,
		"σ" => σ∘(ε(dafh)),
		#[ "ψfhs_$i"  => f for (i, f) in enumerate(ψfhs) ]...,
		[ "ϕfhs_$i"  => f for (i, f) in enumerate(ϕfhs) ]...,
		[ "ϕ0fhs_$i" => f for (i, f) in enumerate(ϕ0fhs) ]...,
		[ "dϕhs_squared$i"   => d for (i, d) in enumerate(dϕhs_squared) ]...,
		#[ "C$i" => C for (i, C) in enumerate(Cs) ]...,
		#[ "C0$i" => C0 for (i, C0) in enumerate(C0s) ]...,	
		"Ik0" => Ik0h,
		"Ik" => Ikh,
		]
	data_TO_in = ["dafh" => dafh, #"ddfh" => ddfh, 
				"σ_vmh" => σ_vmh,
				"σ" => σ∘(ε(dafh)),
				[ "ϕfhs_$i"  => f for (i, f) in enumerate(ϕfhs) ]...,
				[ "ϕ0fhs_$i" => f for (i, f) in enumerate(ϕ0fhs) ]...,
				[ "dϕhs_squared$i"   => d for (i, d) in enumerate(dϕhs_squared) ]...,
				#[ "C$i" => C for (i, C) in enumerate(Cs) ]...,
				#[ "C0$i" => C0 for (i, C0) in enumerate(C0s) ]...,	
				"Iσ_vm"=> (σ_vm(dafh)), #"Is" => (1-H∘φh)*safh,
				#"Ik0" => Ik0h,
				#"Ik" => Ikh,
				]
	plots_directory = plotsdir(model_name)
	mkpath(plots_directory)
	writevtk(get_triangulation(spaces["V_φ"]),joinpath(plots_directory,"$(sname)_TO"),cellfields = data_TO)
	writevtk(spaces["Ωs"].dΩin.quad.trian, joinpath(plots_directory,"$(sname)_TO_Omega"),cellfields= data_TO_in)
	println("donew writing TO")

	## Reconstruction Results

	# final state
	osf = ofh.free_values #deepcopy(get_state(state_maps["of"]).free_values)
	df = deepcopy(get_state(state_maps["d"]).free_values)
	ϕf = deepcopy([(get_state(state_maps["ϕr"][j]).free_values) for j in 1:length(Sϕ)])
	osfh = FEFunction(spaces["V_φ"],osf)
	dfh = FEFunction(spaces["Ud"](Ωact),df)# fxx(x) = x[1]≈0.0 && x[2]≈0.0
	Ik_IPf = state_maps["Ik"](df)
	Ik_IPfh = FEFunction(spaces["UIk"](Ωact),Ik_IPf)
	Uϕs = [spaces["Uϕ"](Ωact,i) for i in Sϕ]
	ϕfh = [FEFunction(Uϕs[j],ϕf[j]) for j in 1:length(Sϕ)]
	# Compute loss 
	L = sum(d_to_l(osf, df, j) for j ∈ 1:length(Sϕ))

	# unloaded state
	sh0 = interpolate(SymTensorValue((0.0, 0.0, 0.0)), spaces["Us"])
	dh0 = interpolate(VectorValue(0.0, 0.0), spaces["Ud"](Ωact))
	Ik_IP0 = deepcopy(state_maps["Ik"](dh0))
	Ik_IP0h = FEFunction(spaces["UIk"](Ωact),Ik_IP0)
	ϕ₀s = [state_maps["ϕ0"][j](Ik_IP0) for j in 1:length(Sϕ)] 
	ϕ₀hs = [FEFunction(Uϕs[j],ϕ₀s[j]) for j in 1:length(Sϕ)]

	# initial state
	osi = oih.free_values
	di = dih.free_values 
	Ik_IPi = state_maps["Ik"](di) 
	Ik_IPih = FEFunction(spaces["UIk"](Ωact),Ik_IPi)
	ϕi = [deepcopy(state_maps["ϕr"][j](Ik_IPi)) for j in 1:length(Sϕ)]
	osih = FEFunction(spaces["V_φ"],osi)
	dih = FEFunction(spaces["Ud"](Ωact),di)
	ϕih = [FEFunction(Uϕs[j],ϕi[j]) for j in 1:length(Sϕ)]

	# actual/observed state
	osa = oah.free_values 
	da = dah.free_values 
	Ik_IPa = state_maps["Ik"](da) 
	Ik_IPah = FEFunction(spaces["UIk"](Ωact),Ik_IPa)
	ϕa = [state_maps["ϕr"][j](Ik_IPa) for j in 1:length(Sϕ)]
	osah= FEFunction(spaces["V_φ"],osa)
	dah = FEFunction(spaces["Ud"](Ωact),da)
	ϕah = [FEFunction(Uϕs[j],ϕa[j]) for j in 1:length(Sϕ)]

	data_IP = [
		"df" => dfh, 
		"di" => dih,
		"da" => dah,
		"di-da" => (dih-dah),
		"df-da" => (dfh-dah),

		"σvmf" => σ_vm(dafh),
		"σvmi" => σ_vm(dih),
		"σvma" => σ_vm(dah),
		
		["ϕ0_$i" => f for (i, f) in enumerate(ϕ₀hs)]...,
		
		["ϕf_$i" => f for (i, f) in enumerate(ϕfh)]...,
		["ϕi_$i" => f for (i, f) in enumerate(ϕih)]...,
		["ϕa_$i" => f for (i, f) in enumerate(ϕah)]...,
		
		["dϕf_$i" => f for (i, f) in enumerate(ϕfh.-ϕ₀hs)]...,
		["dϕi_$i" => f for (i, f) in enumerate(ϕih.-ϕ₀hs)]...,
		["dϕa_$i" => f for (i, f) in enumerate(ϕah.-ϕ₀hs)]...,
			
		"DH(of)" => (DH_IP∘ofh),
		"DH(oi)" => (DH_IP∘oih),
		"DH(oa)" => (DH_IP∘oah),
		
		"H(osf)" => (H_IP ∘ osfh),
		"H(of)"=>(H_IP ∘ ofh), 
		"H(oi)"=>(H_IP ∘ oih), 
		"H(oa)"=>(H_IP ∘ oah), 	
	
		"Ik_IP0" => Ik_IP0h,
		"Ik_IPf" => Ik_IPfh,
		"Ik_IPi" => Ik_IPih,
		"Ik_IPa" => Ik_IPah,
		]
	# Saving geometry
	writevtk(get_triangulation(spaces["V_φ"]),joinpath(plots_directory,"$(sname)_IP"),cellfields = data_IP)
end