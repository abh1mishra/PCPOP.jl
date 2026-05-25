include("traceGrobner.jl")
using PolyChaos: radau, rm_jacobi

"""
    Theorem 5 [BFF24, DI bounds conditional entropy]

    H(A|x₀E) ≥ ∑wᵢ∑inf ρ[M(a|x₀)⊗(Z(a)+Z(a)'+(1−tᵢ)*Z(a)'Z(a)+tᵢZ(a)Z(a)']
        
                   s.t p(ab|xy) = ρ[M(a|x) ⊗ N(b|y)] 
                       M & N & Z commuting
                       B(p) = B₀
                       ZᵢZᵢ' ≤ αᵢ = (3/2)max(1/tᵢ, 1/(1-tᵢ))

    tᵢ, wᵢ nodes and weights Gauss-Radau quadrature

    PolyChaos implments Gauss-Radau

    ∫ g(s) (1-s)^a (1+s)^b ds  for  s ∈ [-1,1]

    ∫ f(t) (1-t)^a t^b dt 2^(1+a+b) for  t ∈ [0,1]

    Where s = 2t - 1 and f(t) = g(2t-1) = g(s)

"""

function gauss_radau(m::Int64)
    a, b = rm_jacobi(m+1, 0, 1)
    t, w = radau(m-1,a,b,1)
    return t, 4*w
end

f(z, t) = (1/t)*(one(z) + (z + z') + (1-t)*(z'*z) + t*(z*z'))

function bff(t, α, γ, k::Int; primal=true)
    # Build monoid
    @pcmonoid M Z[0, 2] a0 a1 b0 b1
    z = M.vertices[1:2]
    @comms [a0, a1] [b0, b1] z
    Projector.([a0, a1, b0, b1])
    build(M)

    # Objective function
    obj = a0*f(z[1], t) + (1-a0)*f(z[2], t)

    # Constraints CHSH violation B(p) ≥ γ
    B  = (1-2*a0)*(1-2*b0) + (1-2*a0)*(1-2*b1)
    B += (1-2*a1)*(1-2*b0) - (1-2*a1)*(1-2*b1)
    tr_ge = [(B, γ)]

    # Constraints on Z
    op_ge =[α - z[1]*z[1]',
            α - z[2]*z[2]',
            α - z[1]'*z[1],
            α - z[2]'*z[2]]

    if primal
        return npa(obj, k; tr_ge=tr_ge, op_ge=op_ge)
    else
        return pcpop(-obj, k, inequalities=op_ge, moments=tr_ge)
    end
end

function entropy_bound(m::Int, γ, k::Int)
    t, w = gauss_radau(m)
    α = [3/2*max(1/t[i], 1/(1-t[i])) for i in 1:m] 

    H = sum(w[i]/(t[i]*log(2))*(1 + bff(t[i], α[i], γ, k)[1]) for i in 1:m)
end


# Parameters
m = 5
γ = 2*sqrt(2)
k = 1

#H = entropy_bund(m, γ, k)

t, w = gauss_radau(m)
α = [3/2*max(1/t[i], 1/(1-t[i])) for i in 1:m] 

model = bff(t[2], α[2], γ, 2, primal=true)[2]
println("Termination Status :", termination_status(model))
println("Objective Value    :", objective_value(model))

sos_model = bff(t[2], α[2], γ, 1, primal=false)
set_optimizer(sos_model, Mosek.Optimizer)
optimize!(sos_model)
println("Termination Status :", termination_status(sos_model))
println("Objective Value    :", objective_value(sos_model))


""""
Primal SDP level 2:

Optimizer terminated. Time: 209.97 
OPTIMAL -15.172788542276617, 
A JuMP Model
├ solver: Mosek
├ objective_sense: MIN_SENSE
│ └ objective_function_type: AffExpr
├ num_variables: 7659
├ num_constraints: 7
│ ├ AffExpr in MOI.EqualTo{Float64}: 1
│ ├ AffExpr in MOI.GreaterThan{Float64}: 1
│ └ Vector{AffExpr} in MOI.PositiveSemidefiniteConeSquare: 5

Dual SDP level 2:

Optimizer terminated. Time: 4.72 
OPTIMAL 15.172788544618173
A JuMP Model
├ solver: Mosek
├ objective_sense: MIN_SENSE
│ └ objective_function_type: AffExpr
├ num_variables: 6127
├ num_constraints: 6642
│ ├ AffExpr in MOI.EqualTo{Float64}: 6637
│ └ Vector{VariableRef} in MOI.PositiveSemidefiniteConeTriangle: 5

(!) Primal is bigger, dual is faster
"""