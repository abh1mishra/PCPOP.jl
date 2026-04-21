include("traceGrobner.jl")

@pcmonoid M a b c d
@comms [a, b] [c, d]
Projector.([a,b,c,d])
build(M)

f = a*c + a*d + b*c - b*d

TM = make_trace_monoid(M, 2)
basis = trace_monomials(TM, 0:1)
p = state_embedding(f, TM)

model_sos = tpop_sos(p, TM, basis)
set_optimizer(model_sos, Mosek.Optimizer)
set_silent(model_sos)
optimize!(model_sos)
println(termination_status(model_sos))
println(objective_value(model_sos))

model = tpop(p, TM, basis)
set_optimizer(model, Mosek.Optimizer)
set_silent(model)
optimize!(model)
println(termination_status(model))
println(objective_value(model))