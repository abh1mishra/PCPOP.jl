include("../../traceGrobner.jl")
import FastGaussQuadrature

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

function gaussradau(m)
    x, v = FastGaussQuadrature.gaussradau(m);
    t = 0.5*(1 .- x);
    w = 0.5*v;

    return (t, w)
end

# f(z, t) = (1/t)*(one(z) + (z + z') + (1-t)*(z'*z) + t*(z*z'))

f(A,z,t) = A * (z + z' + (1-t)*z'*z) + t*z*z'
# F[a] * (Z[a] + Dagger(Z[a]) + (1-ti)*Dagger(Z[a])*Z[a]) + ti*Z[a]*Dagger(Z[a])

function bff(t,w, γ, k::Int; primal=true,canonical=true)
    # Build monoid
    setup_start = time()
    @pcmonoid M Z[0, 2] a0 a1 b0 b1
    z = M.vertices[1:2]
    @comms [a0, a1] [b0, b1] z
    Projector.([a0, a1, b0, b1])
    build(M)
    Id = one(M)
    # Constraints CHSH violation B(p) ≥ γ
    B  = (1-2*a0)*(1-2*b0) + (1-2*a0)*(1-2*b1)
    B += (1-2*a1)*(1-2*b0) - (1-2*a1)*(1-2*b1)
    tr_ge = [[B, γ]]
    basis_principal = mons_at_level(M, k)
    basis = basis_principal
    if !primal
        H = 0.0
        obj = f(a0, z[1], t[1]) + f(1-a0, z[2], t[1])

        model,S,V,mons,LMI = npa_dual(obj,basis,basis_principal; tr_ge=tr_ge, min=true,change_objective=true)
        setup_stop = time()
        elapsed_setup = setup_stop - setup_start
        solve_start = time()
        set_optimizer(model, Mosek.Optimizer)
        set_silent(model)
        optimize!(model)
        ov1 = objective_value(model)
        H += w[1]/(t[1]*log(2))*(1 + ov1)
        old_obj = obj
        for i in 2:length(t)
            obj = f(a0, z[1], t[i]) + f(1-a0, z[2], t[i])
            S = S+old_obj-obj
            model,V = model_new_obj(model,S,V,mons,LMI,-1)
            optimize!(model);ovi = objective_value(model)
            old_obj = obj
            H += w[i]/(t[i]*log(2))*(1 + ovi)
        end
        solve_stop = time()
        elapsed_solve = solve_stop - solve_start
        return H, elapsed_setup, elapsed_solve
    else
        if canonical
            model,D,_ = npa(0, basis,basis_principal; tr_ge=tr_ge,min=true)
        else
            model,D,_ = npa_nc(0, basis,basis_principal; tr_ge=tr_ge,min=true)
        end
        setup_stop = time()
        elapsed_setup = setup_stop - setup_start
        set_silent(model)
        solve_start = time()
        H = 0.0
        for i in 1:length(t)
            obj = f(a0, z[1], t[i]) + f(1-a0, z[2], t[i])
            obj_p = 0.0
            obj_poly=real_rep(Polynomial(obj))
            for (m,c) in obj_poly
                obj_p += c*D[m]
            end
            @objective(model, Min, obj_p)
            set_optimizer(model, Mosek.Optimizer)
            optimize!(model)
            ovi = objective_value(model)
            H += w[i]/(t[i]*log(2))*(1 + ovi)
        end
        solve_stop = time()
        elapsed_solve = solve_stop - solve_start
        return H, elapsed_setup, elapsed_solve
    end

end

function entropy_bound(m::Int, γ, k::Int;t=[],w=[],primal=true,canonical=true)
    if isempty(t) || isempty(w)
        t, w = gaussradau(m)
        t = t[2:end]
        w = w[2:end]
    end
    return bff(t,w, γ, k; primal=primal, canonical=canonical)
end


# Parameters
m = 8
γ = 2*sqrt(2)
k = 2

t=[0.02247939, 0.11467905, 0.26578982, 0.45284637, 0.64737528,
              0.81975931, 0.94373744, 1.        ]
w=[0.05725441, 0.12482395, 0.1735074 , 0.19578608, 0.18825877,
              0.15206531, 0.09267908, 0.015625  ]