include("traceGrobner.jl")

# Build the monoid
@pcmonoid M a b c
Projector.([a, b, c])
@comms a b c
build(M)
# Objective function.
p = (1-a)*b*c + a*(1-b)*c + a*b*(1-c) - (1-a)*(1-b)*(1-c)

k = 3
TM = make_trace_monoid(M, 2*k, tracial=false)
# Constraints
S = []
T = []

mA = mons_at_level(a, k)
mB = mons_at_level(b, k)
mC = mons_at_level(c, k)

R = [state(u*v*w, TM) - state(u, TM)*state(v, TM)*state(w, TM) for u in mA for v in mB for w in mC]

R = [state(a*b*c,TM) - state(a,TM)*state(b,TM)*state(c,TM)]

basis = trace_monomials(TM, 0:k)

# Optimize semidefinite relaxation
model = tpop(p, TM, basis, equalities=R)
set_optimizer(model, Mosek.Optimizer)
optimize!(model)