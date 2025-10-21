function TOforwardpass(θ,measures,spaces,functionals,φh0)

    α1 = θ["α1"]
    α2 = θ["α2"]
    α3 = θ["α3"]
    α4 = θ["α4"]
    σ_lim = θ["σ_lim"]
    Sϕ = θ["Sϕ"]
    j3val = θ["j3val"]
    hsnf = θ["hsnf"]
    Π_TO = θ["Π_TO"] 

    V_φ = spaces["V_φ"]
    U_reg = spaces["U_reg"]
    V_reg = spaces["V_reg"]
    V_φ_ = spaces["V_φ_"]
    U_φ_ = spaces["U_φ_"]
    Uϕ = spaces["Uϕ"]
    Vϕ = spaces["Vϕ"]
    Ud = spaces["Ud"]
    Vd = spaces["Vd"]
    Us = spaces["Us"]
    Vs = spaces["Vs"]
    UIk = spaces["UIk"]
    VIk = spaces["VIk"]
    Ωs = spaces["Ωs"]

    σ_vm = functionals["σ_vm"]
    σ = functionals["σ"]
    K = functionals["K"]
    Π = functionals["Π"]
    a_da = functionals["a_da"]
    l_da = functionals["l_da"]
    ada_cut = functionals["ada_cut"]  
    add_cut = functionals["add_cut"]
    a_dd = functionals["a_dd"]
    l_dd = functionals["l_dd"]
    a_daa = functionals["a_daa"]
    l_daa = functionals["l_daa"]
    a_s = functionals["a_s"]
    l_s = functionals["l_s"]
    a_Ik = functionals["a_Ik"]
    l_Ik = functionals["l_Ik"]
    a_ϕ = functionals["a_ϕ"]
    l_ϕ = functionals["l_ϕ"]
    a_ϕ_in = functionals["a_ϕ_in"]
    l_ϕ_in = functionals["l_ϕ_in"]
    a_ψ = functionals["a_ψ"]
    l_ψ = functionals["l_ψ"]
    j1f = functionals["j1f"]
    j2f = functionals["j2f"]
    j3f = functionals["j3f"]
    j4f = functionals["j4f"]
    Vol = functionals["Vol"]
    Cσ = functionals["Cσ"]

    obs_points = electrode_locations(hsnf)
    model = get_background_model(get_triangulation(V_φ))
    φh = zero(V_φ)
    
    state_collection = GridapTopOpt.EmbeddedCollection_in_φh(model,φh) do _φh 

      update_collection!(Ωs,_φh)
      #writevtk(Ωs.dΩin.quad.trian,"current_Ω1.vtu",cellfields=["χd"=>Ωs.χd])
      Udi = Ud(Ωs.Ωact)
      Vdi = Vd(Ωs.Ωact)
      U_φdi = MultiFieldFESpace([V_φ,Udi])
      UIki = UIk(Ωs.Ωact)
      VIki = VIk(Ωs.Ωact)
      U_φIki = MultiFieldFESpace([V_φ,UIki])
      Vϕs = [Vϕ(Ωs.Ωact,i) for i ∈ Sϕ]
      Uϕs = [Uϕ(Ωs.Ωact,i) for i ∈ Sϕ]
      χs = [Ωs[Symbol("χ$(i)")] for i ∈ Sϕ]
      φdh = interpolate([_φh; zero(Udi)],U_φdi)
      #φIkh = interpolate([_φh; zero(UIki)],U_φIki)
      #Ikh = zero(UIki)

      state_map_dd = AffineFEStateMap(add_cut,l_dd,Udi,Vdi,V_φ)#,_φh)
      state_map_da = AffineFEStateMap(ada_cut,l_da,Udi,Vdi,V_φ)#,_φh)
      function objective_j4(Uϕi,Vϕi) 
        assem_U_i = SparseMatrixAssembler(Uϕi,Vϕi)
        assem_deriv_i = SparseMatrixAssembler(V_φ,V_φ) 
        StateParamMap(j4f,Uϕi,Vϕi,assem_U_i,assem_deriv_i)
      end

      (;
        :state_map_da => state_map_da,
        :state_map_dd => state_map_dd,
        :state_map_Ik => AffineFEStateMap(a_Ik(Π_TO),l_Ik(Π_TO),UIki,VIki,U_φdi),
        :state_map_Ik0 => AffineFEStateMap(a_Ik(Π_TO),l_Ik(Π_TO),UIki,VIki,U_φdi),
        :state_map_ϕs => map((Uϕi,Vϕi,χ) -> AffineFEStateMap(a_ϕ_in(χ,Π_TO),l_ϕ_in,Uϕi,Vϕi,U_φdi), Uϕs,Vϕs,χs),
        :state_map_ϕ0s => map((Uϕi,Vϕi,χ) -> AffineFEStateMap(a_ϕ_in(χ,Π_TO),l_ϕ_in,Uϕi,Vϕi,U_φdi), Uϕs,Vϕs,χs),
        :objective_j1 => StateParamMap(j1f,state_map_da),
        :objective_j2 => StateParamMap(j2f,state_map_dd),
        :obs_operators_1s => map(Uϕi -> FEObservationOperator(obs_points, Uϕi), Uϕs),
        :obs_operators_2s => map(Uϕi -> FEObservationOperator(obs_points, Uϕi), Uϕs),
        :objective_j4s => objective_j4.(Uϕs,Vϕs),
        :C => map(Ci -> StateParamMap(Ci,state_map_dd),[Vol,]),
      )
    end

    function s_to_j3j4(φ,d,j)
        #Ik = state_collection.state_map_Ik([φ;d]).+0
        #Ik0 = state_collection.state_map_Ik0([φ;d*0]).+0
        ϕ = state_collection.state_map_ϕs[j]([φ;d]).+0
        ϕ0 = state_collection.state_map_ϕ0s[j]([φ;d*0]).+0
        e = state_collection.obs_operators_1s[j](ϕ0) - state_collection.obs_operators_2s[j](ϕ)
        j3 = -norm(e) 
        # j4 = state_collection.objective_j4s[j](ϕ,φ)
        # j3,j4
        # j3=0
        j4=0
        j3,j4
    end

    # state_maps = Dict(
    #     "da"=>state_collection.state_map_da,
    #     "dd"=>state_collection.state_map_dd,
    #     "ϕ0"=>state_collection.state_map_ϕ0s,
    #     "ϕ" => state_collection.state_map_ϕs,
    #     "Ik"=>state_collection.state_map_Ik,
    #     "Ik0"=>state_collection.state_map_Ik0,
    #     "obs_operators_1"=>state_collection.obs_operators_1s,
    #     "obs_operators_2"=>state_collection.obs_operators_2s,)

    i=1
    j1s = Float64[]
    j3s = Float64[]
    c1s = Float64[]

    function φ_to_jc(φ)
      da = state_collection.state_map_da(φ)
      j1 = state_collection.objective_j1(da,φ)

      #s = state_map_sa(da).+0
      #Ik0 = state_collection.state_map_Ik0([φ;da]).+0
      #Ik = state_collection.state_map_Ik([φ;da]).+0)

      #dd = state_collection.state_map_dd(φ)

      # i += 1
      # Zygote.ignore() do   
      #   writevtk(Ωs.dΩin.quad.trian,"tmpvca_$(i)",cellfields=["dd" => FEFunction(Ud(Ωs.Ωact),dd)])
      # end

      j2=0.0 #state_collection.objective_j2(dd,φ)      
      j3=0.0
      j4=0.0
      for j ∈ eachindex(Sϕ)
        j3j,j4j = s_to_j3j4(φ,da,j) 
        j3 += j3j
        j4 += j4j
      end


      #j3,j4 = s_to_j3j4(φ,da,1) 

      j = α1*j1 + α2*j2 + α3*j3 + α4*j4 

      @show j1,j2,j3,j4
      @show α1*j1,α2*j2,α3*j3,α4*j4


      c1 = state_collection.C[1](da,φ)
      #c2 = constraints[2](u,φ)^(1/p) - σ_lim # pᵗʰ root for p-norm

      Zygote.ignore() do   
        push!(j1s,j1)
        push!(j3s,j3)
        push!(c1s,c1)
      end

      [j,c1]#,-(j3-j3val)]
    end

    pcfs = CustomEmbeddedPDEConstrainedFunctionals(φ_to_jc,1,state_collection) # atm we just use the representive state map for the correct space for the input variable (the chain rule will interpolate on the state map space)

    pcfs, j1s, j3s, c1s 

  end