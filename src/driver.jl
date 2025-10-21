function run_problem(θ,sname)  
  measures = get_measures(θ)
  spaces = get_spaces(θ,measures)
  functionals = get_functionals(θ,measures,spaces)

  φh0 = get_initial_φh(θ,spaces) 
  φh0, ls_evo, vel_ext, ls_evo_IP, vel_ext_IP = optimise1(θ,measures,spaces,φh0,spaces["U_reg"],spaces["V_reg"],θ["maxiter_TO"],functionals) # updates internally the state_maps
  TO_pcfs, j1s, j3s, c1s = TOforwardpass(θ,measures,spaces,functionals,φh0)
	update_collection!(spaces["Ωs"],φh0)
	reinit!(ls_evo,φh0)
  φh_bg, Js_TO, ls_evo = optimise(θ,measures,spaces,TO_pcfs,φh0,θ["maxiter_TO"],ls_evo,vel_ext,functionals) # updates internally the state_maps

  φh = φh_bg 
  oah = get_actual_o(θ,measures,spaces, ls_evo_IP) 
  oih = get_initial_os0h_bg(θ,measures,spaces,ls_evo_IP) 
  IP_pcfs, dah, dih, d_to_l, state_maps_IP, ls = IPforwardpass(θ,spaces,functionals,oah,oih,φh)
  update_collection!(spaces["Ωs_IP"],oih)
  ofh, Js = optimise(θ,measures,spaces,IP_pcfs,oih,θ["maxiter_IP"],ls_evo_IP,vel_ext_IP,functionals) # updates internally the state_maps
  saverun(θ, φh , ofh, oih ,oah, dah,dih, spaces, functionals, sname, TO_pcfs, state_maps_IP,d_to_l,j1s, j3s, c1s, ls)
end