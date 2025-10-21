function get_functionals(θ, measures,spaces)

  η_coeff_TO = θ["η_coeff_TO"]
  η_coeff_IP = θ["η_coeff_IP"]
  η_coeff_ϕ = θ["η_coeff_ϕ"]
  n = θ["n"]
  σ_lim = θ["σ_lim"]
  k_ad = θ["k_ad"]
  k_af = θ["k_af"]
  k_df = θ["k_df"]
  k_aad = θ["k_aad"]
  r_σ = θ["r_σ"]
  E = θ["E"]
  ν = θ["ν"]
  vf = θ["vf"]
  ϵ_ϕ = θ["ϵ_ϕ"]
  ϵ_ψ = θ["ϵ_ψ"]
  α_ψ = θ["α_ψ"]
  #k∞ = θ["k∞"]
  Π_TO = θ["Π_TO"]
  Π_IP = θ["Π_IP"]
  α_Gd = θ["α_Gd"]
  αfd = θ["αfd"]

  dΓ_d1 = measures["dΓ_d1"]
  dΓ_d1_meas = sum(∫(1)dΓ_d1)
  dΓ_d2 = measures["dΓ_d2"]
  dΓ_d2_meas = sum(∫(1)dΓ_d2)
  dΓ_f1 = measures["dΓ_f1"]
  dΓ_f1_meas = sum(∫(1)dΓ_f1)
  dΓ_c = measures["dΓ_c"]
  dΓ_c_meas = sum(∫(1)dΓ_c)
  dΩ = measures["dΩ"]
  dΩ_D2 = measures["dΩ_D2"]
  dΩ_D2_meas = sum(∫(1)dΩ_D2)
  n_Γc = measures["n_Γc"]
  #dΩ_fk = measures["dΩ_fk"]
  #dΩ_fk_meas = sum(∫(1)dΩ_fk)
  #dΩ_dk = measures["dΩ_dk"]

  Ωs = spaces["Ωs"]

  vol_D = sum(∫(1)dΩ)

  ## Interpolation
  model = get_background_model(dΩ.quad.trian)

  polys = get_polytopes(model)  
  poly = first(polys)
  if poly == TRI
    h = minimum(get_element_diameters(model))
    hₕ = get_element_diameter_field(model)
  elseif poly == QUAD
    el_Δ = get_el_Δ(model)
    h = maximum(el_Δ)
    hₕ = h
  end

  interp = SmoothErsatzMaterialInterpolation(η=η_coeff_TO * h)
  I, H, DH, ρ = interp.I, interp.H, interp.DH, interp.ρ

  interp_ϕ = SmoothErsatzMaterialInterpolation(η=η_coeff_ϕ * h, ϵ=ϵ_ϕ)
  I_ϕ, H_ϕ, DH_ϕ, ρ_ϕ = interp_ϕ.I, interp_ϕ.H, interp_ϕ.DH, interp_ϕ.ρ

  interp_ψ = SmoothErsatzMaterialInterpolation(η=η_coeff_TO * h, ϵ=ϵ_ψ)
  I_ψ, H_ψ, DH_ψ, ρ_ψ = interp_ψ.I, interp_ψ.H, interp_ψ.DH, interp_ψ.ρ

  interp_IP = SmoothErsatzMaterialInterpolation(η=η_coeff_IP * h, ϵ=0)
  I_IP, H_IP, DH_IP, ρ_IP = interp_IP.I, interp_IP.H, interp_IP.DH, interp_IP.ρ

  # Structural params
  λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
  μ = E / (2 * (1 + ν))
  σ(ε) = λ * tr(ε) * one(ε) + 2 * μ * ε
  function σ_vm(uh)
    hydrostatic_σ = tr(σ ∘ (ε(uh))) / 3
    deviatoric_σ = σ ∘ (ε(uh)) - hydrostatic_σ
    j_2 = 1 / 2 * deviatoric_σ ⊙ deviatoric_σ        # deviatoric stress invariant = j_2 = 1/2 ∑σ_dev_ij^2
    sigma_vm = sqrt ∘ (3 / 2 * j_2)
  end
  fa(x) = VectorValue(-1, -1) # actual fore - need to scale this by the area of the region
  fd(x) = VectorValue(0, αfd) # dummy force
  fd2(x) = VectorValue(0, 1) # dummy force

  # Electric params
  k = 1 # 1e-5  # conductivity constant
  function Π(π_i)
    Π₁₁ = π_i # 1e-7 #-102.2e-11 *1e8#6.6e-11 * 1e5 # (100) silicon 
    Π₁₂ = Π₁₁ #53.5e-11 *1e8#-1.1e-11 * 1e5
    Π₄₄ = Π₁₂ #-13.6e-11 *1e8# 138.1e-11 * 1e5
    Π_i = SymFourthOrderTensorValue((Π₁₁, Π₁₂, 0, Π₁₂, Π₁₁, 0, 0, 0, Π₄₄)) # piezoresistive tensor
  end

  D = 2
  Id = one(SymTensorValue{D,Float64})
  K(Π, σ) = k * (Id - Π ⊙ (σ)) # conductivity (or +???))
  K(Π) = σ -> K(Π, σ)
  K0 = k * Id

  # Stress prjection 
  r_σ_scaled = r_σ * h

  # Misc
  p = 12
  power(x) = x^p
  downward_disp(d) = abs(d[2])
  power2(x) = x^2
  myrelu(x) = max(0, 1e13x)

  
  x1 = 0.145
  x2 = 0.15
  #@assert x2-x1 ≈ dΓ_f1_meas "dΓ_f1_meas must be equal to x2 - x1 (change x1 to match)"
  wf(x) = 1/(x2-x1) 
  #wf(x) = 6*(x[1] - x1)*(x2 - x[1])/((x2 - x1)^3)
  #wf(x) = 2*(x[1] - x1)/((x2 - x1)^2)
  #wf(x) = 3*(x[1] - x1)^2/((x2 - x1)^3)
  #wk(x) = wf(x) * (x2-x1)

  α_Gd = α_Gd
  γ_Gd(h) = α_Gd*(λ + μ)*h^3
  a_s_Ω(d,s) = ε(s) ⊙ (σ ∘ ε(d)) # Elasticity
  j_s_k(d,s) = mean(γ_Gd ∘ hₕ)*(jump(Ωs.n_Γg ⋅ ∇(s)) ⋅ jump(Ωs.n_Γg ⋅ ∇(d)))
  v_s_ψ(d,s) = (Ωs.χd)*(d⋅s) # Isolated volume term

  # Weak forms
  #a_el(d, v, φ) = ∫((I ∘ φ) * (ε(v) ⊙ (σ ∘ ε(d))))dΩ

  a_da(d, v, φ) = a_el(d, v, φ) + ∫(k_ad * d ⋅ v) * dΓ_d1  + ∫(k_af*d⋅v)dΓ_f1 
  ada_cut(d, s, φ) = ∫(a_s_Ω(d, s) + v_s_ψ(d, s))Ωs.dΩin + ∫(j_s_k(d, s))Ωs.dΓg + ∫(k_ad * d ⋅ s) * dΓ_d1  + ∫(k_af*d⋅s)dΓ_f1
  l_da(v, φ) = ∫(fa⋅v / dΓ_f1_meas)dΓ_f1 

  a_dd(d, v, φ) = a_el(d, v, φ) + ∫(k_df * d ⋅ v) * dΓ_f1
  add_cut(d,s,φ) = ∫(a_s_Ω(d,s) + v_s_ψ(d,s))Ωs.dΩin + ∫(j_s_k(d,s))Ωs.dΓg + ∫(k_df * d ⋅ s) * dΓ_f1

  l_dd(v, φ) = ∫(fd ⋅ v / dΓ_d1_meas)dΓ_d1

  a_daa(d, v, φ) = a_el(d, v, φ) + ∫(k_aad * d ⋅ v) * dΓ_d1 + ∫(k_af * d ⋅ v)dΓ_f1
  l_daa(v, φ) = ∫(  fa ⋅ v )dΓ_f1

  a_d(φ) = (d, v, o) -> a_el(d, v, φ) + ∫(1e13 * (H ∘ (o)) * d ⋅ v)dΩ  #+ ∫(k_aad*d⋅v) * dΓ_d1 #+ ∫(k*d⋅v)*dΓ_c
  l_d(φ) = (v, o) -> ∫(wf * ( fa ⋅ v) )dΓ_f1

  a_s(p, q, d) = ∫(p ⊙ q)dΩ + ∫(r_σ_scaled * ∇(p) ⊙ ∇(q))dΩ
  l_s(q, d) = ∫((σ ∘ ε(d)) ⊙ q)dΩ

  a_Ik(Π, Ik, q, (φ, d)) = ∫(Ik ⊙ q)dΩ
  l_Ik(Π, q, (φ, d)) = ∫( (K(Π) ∘ ( (σ ∘ ε(d)))) ⊙ q)dΩ

  a_Ik(Π) = (Ik, q, (φ, d)) -> a_Ik(Π, Ik, q, (φ, d))
  l_Ik(Π) = (q, (φ, d)) -> l_Ik(Π, q, (φ, d))

  a_ϕ(ϕ, η, Ik) = ∫((Ik + 1e-10) ⋅ ∇(ϕ) ⋅ ∇(η))dΩ
  l_ϕ(η, Ik) = ∫(0 * η)dΩ

  γg = 10
  a_ϕ_in(ϕ, η, (φ,Ik), χ, Π ) = ∫((K(Π) ∘ ( (σ ∘ ε(Ik)))) ⋅ ∇(ϕ) ⋅ ∇(η))Ωs.dΩin +
                     ∫((γg * mean(hₕ)) * jump(Ωs.n_Γg ⋅ ∇(η)) * jump(Ωs.n_Γg ⋅ ∇(ϕ)))Ωs.dΓg +
                     ∫(χ * η * ϕ)Ωs.dΩin
  a_ϕ_in(χ,Π) = (ϕ, η, (φ,Ik)) -> a_ϕ_in(ϕ, η, (φ,Ik), χ,Π )

  l_ϕ_in(η, (φ,Ik)) = ∫(0 * η)Ωs.dΩin

  γg = 10
  a_ϕ_in_I(ϕ, η, (φ,Ik), χ, Π ) = ∫(Ik ⋅ ∇(ϕ) ⋅ ∇(η))Ωs.dΩin +
                     ∫((γg * mean(hₕ)) * jump(Ωs.n_Γg ⋅ ∇(η)) * jump(Ωs.n_Γg ⋅ ∇(ϕ)))Ωs.dΓg +
                     ∫(χ * η * ϕ)Ωs.dΩin
  a_ϕ_in_I(χ,Π) = (ϕ, η, (φ,Ik)) -> a_ϕ_in_I(ϕ, η, (φ,Ik), χ,Π )

  l_ϕ_in_I(η, (φ,Ik)) = ∫(0 * η)Ωs.dΩin



  a_ψ(ψ, ξ, φ) = ∫((I_ψ ∘ φ) * (∇(ψ) ⋅ ∇(ξ)))dΩ + ∫((H ∘ (φ)) * (α_ψ * ψ * ξ))dΩ
  l_ψ(ξ, φ) = ∫(0 * ξ)dΩ

  a_f(α) = (of, q, o) -> ∫(α^2 * ∇(of) ⋅ ∇(q) + of * q)dΩ
  l_f(α) = (q, o) -> ∫(o * q)dΩ

  #k∞ = k∞# 1e8
  #oh = interpolate(potato_level_set,V_φ)
  function kf(d, x, oh)
    x_in = VectorValue(min(max(x[1] + d[1], 0.0), 0.06), min(max(x[2] - d[2], 0.0), 0.07))
    #H_IP(oh(x_in))
    #H_IP(oh(x_in)) #
    max(0,oh(x_in)) # only take posistive part 
  end
  kf(o) = (d, x) -> kf(d, x, o)
  x(x) = x
  m(k∞, w, d, v, o) = ∫(k∞ * (kf(o) ∘ (w, x)) * (v⋅n_Γc) )dΓ_c # + ∫(k∞*(H_IP∘(o))*d*v)Ωs.dΩin #
  a(d, s, φ) = ∫(a_s_Ω(d, s))Ωs.dΩin  + ∫(k_af*d⋅s)dΓ_f1  + ∫(j_s_k(d, s))Ωs.dΓg # + ∫(k_ad * d ⋅ s) * dΓ_d1
  l(f, s) = ∫(f * (fa ⋅ s))dΓ_f1
  res_d(f, k∞, φ) = (d, s, o) ->   m(k∞, d, d, s, o) + a(d, s, φ) - l(f, s)

  # Objectives 
  j1f(da, φ) = ∫(-(downward_disp ∘ (da)))dΓ_d1
  j2f(dd, φ) = ∫((ε(dd) ⊙ (σ ∘ ε(dd))))Ωs.dΩin
  j3f((ψ, ϕ, ϕ0), φ) = ∫(-(1 - (H ∘ φ)) * (power2 ∘ (ϕ - ϕ0)) * ψ)dΩ
  j4f(ϕ, φ) = ∫((1) * (∇(ϕ) ⋅ ∇(ϕ)))dΩ

  loss(dϕ0) = (dϕ, o) -> ∫(power2 ∘ (dϕ - dϕ0))dΩ

  r(of, o) = ∫(∇(of) ⋅ ∇(of))dΩ

  # Constraints 
  #Vol(da, φ) = ∫(((ρ ∘ φ) - vf) / vol_D)dΩ
  Vol(d,φ) = ∫(1/vol_D)Ωs.dΩin - ∫(vf/vol_D)dΩ
  Cσ(da, φ) = ∫((power ∘ ((I ∘ φ) * (σ_vm(da)))) / vol_D)dΩ

  functionals = Dict(
    "σ_vm" => σ_vm,
    "σ" => σ,
    "K" => K,
    "Π" => Π,
    "a_da" => a_da,
    "ada_cut" => ada_cut,
    "l_da" => l_da,
    "a_dd" => a_dd,
    "add_cut" => add_cut,
    "l_dd" => l_dd,
    "a_daa" => a_daa,
    "l_daa" => l_daa,
    "a_s" => a_s,
    "l_s" => l_s,
    "a_Ik" => a_Ik,
    "l_Ik" => l_Ik,
    "a_ϕ" => a_ϕ,
    "l_ϕ" => l_ϕ,
    "a_ϕ_in" => a_ϕ_in,
    "l_ϕ_in" => l_ϕ_in,
    "a_ϕ_in_I" => a_ϕ_in_I,
    "l_ϕ_in_I" => l_ϕ_in_I,
    "a_ψ" => a_ψ,
    "l_ψ" => l_ψ,
    "a_d" => a_d,
    "l_d" => l_d,
    "a_f" => a_f,
    "l_f" => l_f,
    "res_d" => res_d,
    "j1f" => j1f,
    "j2f" => j2f,
    "j3f" => j3f,
    "j4f" => j4f,
    "loss" => loss,
    "r" => r,
    "Vol" => Vol,
    "Cσ" => Cσ,
    "K_0" => K0,
    "H" => H,
    "H_IP" => H_IP,
    "DH_IP" => DH_IP,
    "I_ϕ" => I_ϕ,
  )
  functionals

end