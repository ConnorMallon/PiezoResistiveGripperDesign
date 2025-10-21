using GridapTopOpt: AbstractFEStateMap, AffineFEStateMap
using GridapTopOpt: get_trial_space, get_deriv_space, get_deriv_assembler, get_measure
using Gridap
using Gridap.Fields
using ChainRulesCore
using ChainRulesCore: InplaceableThunk
using Zygote

import GridapTopOpt: adjoint_solve!

struct PDEConstrainedLoss2{A}
    loss 
    dloss
    state_map::A

    function PDEConstrainedLoss2(
        loss :: Function,
        state_map   :: AbstractFEStateMap)

        dloss = state_map.plb_caches[1] # allocate DJ equivalent to the state map pullback cache for the derivate dRdφ 

        T = typeof(state_map)
        return new{T}(loss,dloss,state_map)
    end
end

function Fields.evaluate!(pcl::PDEConstrainedLoss2,φh)
    loss, dloss = pcl.loss,pcl.dloss
    U = get_trial_space(pcl.state_map)

    U_reg = get_deriv_space(pcl.state_map)
    deriv_assem = get_deriv_assembler(pcl.state_map)
    dΩ = get_measure(pcl.state_map)
  
    ## Foward problem
    u, u_pullback = rrule(pcl.state_map,φh)
    uh = FEFunction(U,u)
  
    function ∇!(F::Function,dF)
        j_val, (∂j_∂p, ∂j_∂u) = Zygote.withgradient(loss, φh.free_values, u) 

        @show j_val, ∂j_∂p, ∂j_∂u
        _, dφ_adj         = u_pullback(∂j_∂u) # Compute -dFdu*dudφ via adjoint
        copy!(dF,dφ_adj)
        dF .+= ∂j_∂p
        return j_val
    end

    j = ∇!(loss,dloss)
    return j,dloss
  end

abstract type ParameterisedObjective end

# struct CustomPDEConstrainedFunctionals{N,A} <: ParameterisedObjective
#   φ_to_jc :: Function
#   dj :: Vector{Float64}
#   dc :: Vector{Vector{Float64}}
#   state_map :: A

#     function CustomPDEConstrainedFunctionals(
#       φ_to_jc :: Function,
#       state_map :: AbstractFEStateMap,
#       φ0,
#     )
    
#     # Pre-allocaitng
#     grad = Zygote.jacobian(φ_to_jc, φ0.free_values)
#     dj = grad[1][1,:]
#     dc = [collect(row) for row in eachrow(grad[1][2:end,:])]
#     N = length(dc)
#     A = typeof(state_map)

#     return new{N,A}(φ_to_jc,dj,dc,state_map)
#   end
# end

# function Fields.evaluate!(pcf::CustomPDEConstrainedFunctionals,φh)
#   φ_to_jc,dj,dc = pcf.φ_to_jc,pcf.dj,pcf.dc

#   obj,grad = Zygote.withjacobian(φ_to_jc, φh.free_values)
#   j = obj[1]
#   c = obj[2:end]
#   dj = copy!(dj,grad[1][1,:])
#   dc = copy!(dc,[collect(row) for row in eachrow(grad[1][2:end,:])])

#   return j,c,dj,dc
# end

# function adjoint_solve!(state_map::AffineFEStateMap, du::InplaceableThunk)
#   return adjoint_solve!(state_map, unthunk(du))
# end