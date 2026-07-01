using QuantumNPA, MosekTools, JuMP
import FastGaussQuadrature

"""
    Theorem 5 [BFF24, DI bounds conditional entropy] — QuantumNPA implementation.

    H(A|x₀E) ≥ ∑wᵢ∑inf ρ[M(a|x₀)⊗(Z(a)+Z(a)'+(1−tᵢ)*Z(a)'Z(a)+tᵢZ(a)Z(a)']

                   s.t p(ab|xy) = ρ[M(a|x) ⊗ N(b|y)]
                       M & N & Z commuting
                       B(p) = B₀

    tᵢ, wᵢ : nodes and weights of the Gauss-Radau quadrature.

    Mirrors Benchmark/BFF/pcpop_imp.jl. The Eve operators Z are non-Hermitian
    `generic` operators placed on a separate party (3), so they automatically
    commute with Alice (party 1) and Bob (party 2) but Z and Z' do not commute.

    The NPA relaxation (moment matrix + CHSH inequality) is built once. We keep
    the monomial→variable map `v`, then iterate over the quadrature nodes,
    re-setting only the objective and summing the optimal values.

    NB: this QuantumNPA version's `npa2jump`/`sdp2jump` chokes on a *scalar* `ge`
    inequality (its `add_constraint!` calls `size` on an `AffExpr`), so we build
    the JuMP model directly from `npa2sdp` and add the scalar CHSH constraint
    ourselves.
"""

# Same Gauss-Radau quadrature as pcpop_imp.jl
function gaussradau(m)
    x, v = FastGaussQuadrature.gaussradau(m)
    t = 0.5*(1 .- x)
    w = 0.5*v

    return (t, w)
end

# Same objective building block as pcpop_imp.jl:
# f(A,z,t) = A*(z + z' + (1-t)*z'*z) + t*z*z'
f(A, z, t) = A*(z + conj(z) + (1-t)*conj(z)*z) + t*z*conj(z)

bff_objective(a0, z, t) = f(a0, z[1], t) + f(Id - a0, z[2], t)

function bff(t, w, γ, k::Int;optimize=true)
    # Build operators
    t_start = time()
    a0 = projector(1, 1, 1)
    a1 = projector(1, 1, 2)
    b0 = projector(2, 1, 1)
    b1 = projector(2, 1, 2)
    z = generic(3, 1:2)                       # Eve's non-Hermitian operators

    # CHSH violation constraint B(p) ≥ γ
    B  = (Id - 2*a0)*(Id - 2*b0) + (Id - 2*a0)*(Id - 2*b1)
    B += (Id - 2*a1)*(Id - 2*b0) - (Id - 2*a1)*(Id - 2*b1)

    # Build the relaxation once. Passing the objective and the CHSH polynomial as
    # `ge` ensures Bob's operators enter the moment matrix, so every monomial of
    # both the objective and the constraint appears in `mons`.
    obj0 = bff_objective(a0, z, t[1])
    model = npa2jump_d(obj0, k;sense=:minimise, ge=[B - γ*Id],solver=Mosek.Optimizer)
    if !optimize
        stop_setup = time()
        return nothing, stop_setup - t_start, 0.0
    end
    unset_silent(model)  # don't solve yet, we will iterate over quadrature nodes
    optimize!(model)
    println("Termination status ", termination_status(model))
    H = w[1]/(t[1]*log(2))*(1 + objective_value(model))
    for i in 2:length(t)
        obji = bff_objective(a0, z, t[i])
        model = npa2jump_d(obji, k;sense=:minimise, ge=[B - γ*Id],solver=Mosek.Optimizer)
        optimize!(model)
        println("Termination status ", termination_status(model))
        ovi = objective_value(model)
        H += w[i]/(t[i]*log(2))*(1 + ovi)
    end
    t_stop = time()
    elapsed_total = t_stop - t_start

    return H,elapsed_total
end

function avg_time(total_runs,m::Int, γ, k::Int;t=[],w=[],optimize=true)
    if isempty(t) || isempty(w)
        t, w = gaussradau(m)
        t = t[2:end]
        w = w[2:end]
    end
    bff(t,w, γ, 2)
    total_time = 0.0
    for i in 1:total_runs
        H,t_time = bff(t,w, γ, k;optimize=optimize)
        total_time += t_time
        println("H:", H)
    end
    return total_time / total_runs

end

# Parameters
m = 8
γ = 2*sqrt(2)
k = 2

t = [0.02247939, 0.11467905, 0.26578982, 0.45284637, 0.64737528,
     0.81975931, 0.94373744, 1.0]
w = [0.05725441, 0.12482395, 0.1735074 , 0.19578608, 0.18825877,
     0.15206531, 0.09267908, 0.015625]
