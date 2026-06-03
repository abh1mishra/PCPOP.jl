include("traceGrobner.jl")

println("Building monoid...")
# Build the monoid
@pcmonoid M a[2,0] b[2,0] c[2,0]
Unipotent.(M.vertices)
@comms a b c
build(M)
k = 2
TM = make_trace_monoid(M, 2*k+2, tracial=false) 
# Objective function.
α = sum(state(a[i]*b[1]*c[j], TM) for i in 1:2 for j in 1:2)
γ = sum(state(a[i]*b[1]*c[j], TM)*(-1)^(i+j) for i in 1:2 for j in 1:2)
p = (α - γ)^2/8 - (α + γ)
# Equality constraints
basis = union(trace_monomials(TM, 0:k), [state(a[i]*b[1]*c[j], TM) for i in 1:2 for j in 1:2])
wα = mons_at_level(a, k-1)
wγ = mons_at_level(c, k-1)
R = [state(u*v, TM) - state(u, TM)*state(v, TM) for u in wα for v in wγ]
R = unique([r for r in R if !(r==0)])
println("Building model...")
model = tpop(p, TM, basis, equalities=R)
set_optimizer(model, Mosek.Optimizer)
println("Optimizing...")
optimize!(model)
println("Termination status ", termination_status(model))
println("Optimal value is   ", objective_value(model))