include("traceGrobner.jl")

# Build the monoid
@pcmonoid M a[2,0] b[2,0] c[2,0]
Unipotent.(union([a, b, c]...))
@comms a b c
build(M)
k = 2
TM = make_trace_monoid(M, 2*k, tracial=false) 
# Objective function.
α = sum(state(a[i]*b[1]*c[j], TM) for i in 1:2 for j in 1:2)
γ = sum(state(a[i]*b[1]*c[j], TM)*(-1)^(i+j) for i in 1:2 for j in 1:2)
p = (α - γ)^2/8 - (α + γ)
# Equality constraints
basis = trace_monomials(TM, 0:k)
wα = mons_at_level(a, k)
wγ = mons_at_level(c, k)
R = [state(u*v, TM) - state(u, TM)*state(v, TM) for u in wα for wb in wγ
model = tpop(p,basis,equalities=R)
set_optimizer(model, Mosek.Optimizer)
optimize!(model)
println("Termination status ", termination_status(model))
println("Optimal value is   ", val)