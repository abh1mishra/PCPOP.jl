using PCPOP, JuMP, Mosek, MosekTools

println("Building monoid...")
# Build the monoid
@pcmonoid M a[2, 0] b[2, 0] c[2, 0]
Unipotent.(M.vertices)
@comms a b c
build(M)
k = 2
TM = make_trace_monoid(M, 2*k, tracial = false)
# Objective function.
p = a[2]*b[1]*c[1] + a[1]*b[2]*c[1] + a[1]*b[1]*c[2] - a[2]*b[2]*c[2]
p = state(p, TM)
# Equality constraints
basis = trace_monomials(TM, 0:k)
wα = mons_at_level(a, k)
wγ = mons_at_level(c, k)
R = [state(u*v, TM) - state(u, TM)*state(v, TM) for u in wα for v in wγ]
R = unique([r for r in R if !(r==0)])
println("Building model...")
model = tpop(p, TM, basis, op_eq = R)
set_optimizer(model, Mosek.Optimizer)
println("Optimizing...")
optimize!(model)
println("Termination status ", termination_status(model))
println("Optimal value is   ", objective_value(model))
