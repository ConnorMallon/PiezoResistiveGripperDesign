function IPforwardpass(θ, spaces, functionals, oah, oh0, φh)
    Sϕ = θ["Sϕ"] 
    hsnf = θ["hsnf"]
    αr = θ["αr"]
    k∞ = θ["k∞"]
    Π_IP = θ["Π_IP"]
    f_in = θ["f_in"]
    σ_n = θ["σ_n"]
    V_φ = spaces["V_φ"]
    V_φ_IP = spaces["V_φ_IP"]
    V_diff = spaces["V_diff"]
    Vd = spaces["Vd"]
    Ud = spaces["Ud"]
    Vs = spaces["Vs"]
    Us = spaces["Us"]
    UIk = spaces["UIk"]
    VIk = spaces["VIk"]
    Vϕ = spaces["Vϕ"]
    Uϕ = spaces["Uϕ"]
    Vk = spaces["Vk"]
    U_φ_IP = spaces["V_reg_IP"]
    U_reg_IP = spaces["U_reg_IP"]
    a_d = functionals["a_d"]
    l_d = functionals["l_d"]
    a_s = functionals["a_s"]
    l_s = functionals["l_s"]
    a_Ik = functionals["a_Ik"]
    l_Ik = functionals["l_Ik"]
    a_ϕ = functionals["a_ϕ"]
    l_ϕ = functionals["l_ϕ"]
    a_ϕ_in_I = functionals["a_ϕ_in_I"]
    l_ϕ_in_I = functionals["l_ϕ_in_I"]
    a_f = functionals["a_f"]
    l_f = functionals["l_f"]
    res_d = functionals["res_d"]
    loss = functionals["loss"]
    r = functionals["r"]
    Ωs = spaces["Ωs"]

    UIki = UIk(Ωs.Ωact)
    VIki = VIk(Ωs.Ωact)
    Udi = Ud(Ωs.Ωact)
    Vdi = Vd(Ωs.Ωact)
    Vϕs = [Vϕ(Ωs.Ωact,i) for i ∈ Sϕ]
    Uϕs = [Uϕ(Ωs.Ωact,i) for i ∈ Sϕ]
    χs = [Ωs[Symbol("χ$(i)")] for i ∈ Sϕ]
    dh0 = interpolate(VectorValue(0, 0), Udi)
    d0 = dh0.free_values

    nls = NLSolver(show_trace=true, method=:newton,iterations=20,linesearch=BackTracking())
    state_map_d = ReverseNonlinearFEStateMap(res_d(f_in,k∞,φh),Udi,Vdi,V_φ_IP,V_diff,nls=nls)
    state_map_Ik = AffineFEStateMap((Ik,q,d)->a_Ik(Π_IP)(Ik,q,(φh,d)),(q,d)->l_Ik(Π_IP)(q,(φh,d)),UIki, VIki, Udi)
    state_map_Ik0 = AffineFEStateMap((Ik,q,d)->a_Ik(Π_IP)(Ik,q,(φh,d)),(q,d)->l_Ik(Π_IP)(q,(φh,d)),UIki, VIki, Udi)
    state_map_ϕ0s = 
        map(
            (Uϕ,Vϕ,χ) -> 
        AffineFEStateMap(
        (ϕ,η,d) -> a_ϕ_in_I(ϕ, η, (φh,d), χ, Π_IP ),
        (η,d) -> l_ϕ_in_I(η, (φh,d) ),
        Uϕ, Vϕ, UIki)
        ,Uϕs,Vϕs,χs)
    state_map_ϕs = 
        map(
            (Uϕ,Vϕ,χ) -> 
        AffineFEStateMap(
        (ϕ,η,d) -> a_ϕ_in_I(ϕ, η, (φh,d), χ, Π_IP ),
        (η,d) -> l_ϕ_in_I(η, (φh,d) ),
        Uϕ, Vϕ, UIki)
        ,Uϕs,Vϕs,χs)

    #ramping d0 for initial
    ofh0 = oh0 #FEFunction(V_φ,state_map_of(oh0))
    opd(f,k∞) = FEOperator((u,v)->res_d(f,k∞,φh)(u,v,ofh0),Udi,Vdi)

    opdf(f,k∞,ofh) = FEOperator((u,v)->res_d(f,k∞,φh)(u,v,ofh),Udi,Vdi)
    
    function o_to_d(o)
        of = o

        # Zygote.ignore() do 
        #     d_0 = state_map_d.cache.fwd_cache[3]
        #     Udi = GridapTopOpt.get_trial_space(state_map_d)
        #     dh_0 = FEFunction(Udi,d_0)
        #     ofh = FEFunction(V_φ,of)

        #     dh_ramp = solve(nls,opdf(f_in,k∞/100,ofh))
        #     #dh_ramp,_ = solve!(dh_0,nls,opdf(f_in,k∞/100,ofh))
        #     dh_ramp,_ = solve!(dh_0,nls,opdf(f_in,k∞/10,ofh))
        #     state_map_d.cache.fwd_cache[3] .= dh_ramp.free_values
            
        # end
        # @show sum(state_map_d.cache.fwd_cache[3])

        d = state_map_d(of) .+ 0
    end


    ## ===================== ##
    ## ramping d0 for actual ##
    ## ===================== ##

    #osah = FEFunction(V_φ,state_map_of(oah))
    opda(f,k∞) = FEOperator((u,v)->res_d(f,k∞,φh)(u,v,oah),Udi,Vdi)
    println("ramping1")
    d0 = ramp(opda,nls,f_in,k∞)# zero(Udi) #
    #@save "da.jld2" d0
    # @load "da.jld2" d0
    # d0 = FEFunction(Udi, d0.free_values)
    GridapTopOpt.build_cache!(state_map_d,oah)
    x = state_map_d.cache.fwd_cache[3]
    x .= d0.free_values 
    # finished ramping 
    da = state_map_d(oah) .+ 0 # d0
    dah = FEFunction(Udi, da)
    Ika = state_map_Ik(da) .+ 0
    ϕa = deepcopy([state_map_ϕs[j](Ika) for j ∈ 1:length(Sϕ)])
    ϕah = [FEFunction(Uϕs[j],ϕa[j]) for j in 1:length(Sϕ)]
    obs_points = electrode_locations(hsnf)
    @show σ_n*randn(1)
    @show typeof(σ_n*randn(1))
    obs_values = [[ϕah[j](x) + σ_n*randn() for x in obs_points] for j in 1:length(Sϕ)]
    obs_operators = [FEObservationOperator(obs_points, Uϕs[j]) for j in 1:length(Sϕ)]

    function Ik_to_l(o, Ik, j)
        ϕ = state_map_ϕs[j](Ik)
        e = obs_values[j] - obs_operators[j](ϕ)
        l = norm(e)

        # c = obs_values[j].*0 - obs_operators[j](ϕ)
        # norm(c)
        # sum(ϕ)
        #sum(d)
    end

    state_maps_IP = Dict(
        #"of" => state_map_of,
        "d" => state_map_d,
        "Ik" => state_map_Ik,
        "Ik0" => state_map_Ik0,
        "ϕr" => state_map_ϕs,
        "ϕ0" => state_map_ϕ0s,
        "obs_operators" => obs_operators,
    )


    assem_U = SparseMatrixAssembler(V_φ,V_φ)
    assem_deriv = SparseMatrixAssembler(V_φ,V_φ) 
    objective_r = StateParamMap(r,V_φ,V_φ,assem_U,assem_deriv)

    ls=Float64[]
    function o_to_l(o)

        # Zygote.ignore() do 
        #     ofhi = FEFunction(V_φ,o)
        #     writevtk(get_triangulation(V_φ),"tmpdi",cellfields=["oih"=>functionals["H_IP"]∘ofhi])
        # end

        d = o_to_d(o)
        r = objective_r(o,o)
        Ik = state_map_Ik(d)

        l = sum(Ik_to_l(o, Ik, j) for j ∈ 1:length(Sϕ))

        l = sum(Ik_to_l(o, Ik, j) for j ∈ 1:6)

        #foldl(+, (d_to_l(o, d, j) for j in 1:6))

        # # @show l 

        # Zygote.ignore() do 
        #     push!(ls,l)
        # end
    
        [l+αr*r] 

        #cs = [d_to_l(o,d,j) for j ∈ 1:length(Sϕ)] 
        #c =  d_to_l(o,d,i) 

        #objective_tmp(d,o)

        
    end

    function o_to_l2(o)

        # Zygote.ignore() do 
        #     ofhi = FEFunction(V_φ,o)
        #     writevtk(get_triangulation(V_φ),"tmpdi",cellfields=["oih"=>functionals["H_IP"]∘ofhi])
        # end

        d = o_to_d(o)
        Ik = state_map_Ik(d)
        r = objective_r(o,o)

        cs = [Ik_to_l(o, Ik, j) for j ∈ 1:length(Sϕ)]
    end

    # squared(x) = x^2
    # @show sum(squared.(Zygote.gradient(o->o_to_l2(o)[1], oah.free_values)[1]))
    # @show sum(squared.(Zygote.gradient(o->o_to_l2(o)[2], oah.free_values)[1]))
    # @show sum(squared.(Zygote.gradient(o->o_to_l2(o)[3], oah.free_values)[1]))
    # @show sum(squared.(Zygote.gradient(o->o_to_l2(o)[4], oah.free_values)[1]))
    # @show sum(squared.(Zygote.gradient(o->o_to_l2(o)[5], oah.free_values)[1]))
    # @show sum(squared.(Zygote.gradient(o->o_to_l2(o)[6], oah.free_values)[1]))

    #@load "di.jld2" d0
    println("ramping2")
    d0 = ramp(opd,nls,f_in,k∞)#  zero(Udi) # 
    #@save "di.jld2" d0
    #@load "di.jld2" d0
    d0 = FEFunction(Udi, d0.free_values)
    x = state_map_d.cache.fwd_cache[3]
    x .= d0.free_values 
    dih = d0
    
    pcf = CustomPDEConstrainedFunctionals(o_to_l,0)


    
    function d_to_l(o,d,j)
        Ik = state_map_Ik(d)
        Ik_to_l(o, Ik, j)
    end

    pcf, dah, dih, d_to_l, state_maps_IP,ls
end 
