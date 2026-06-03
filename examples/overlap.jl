include("traceGrobner.jl")
# Initialize local monoids
@ncmonoid A a1 a2
@ncmonoid B b1 b2
@ncmonoid C c1 c2
@ncmonoid BC x1 x2
Unipotent.([a1, a2, b1, b2, c1, c2, x1, x2])
# Build global monoid
M = GraphProductMonoid("M",[A, B, C, BC])
@comms A B C
@comms A BC
build(M)
# Objective function.
p  = a1*(b1 + b2) + a2*(b1 - b2)
p += a1*(c1 + c2) + a2*(c1 - c2)
p += a1*(x1 + x2) + a2*(x1 - x2) 
# Optimize semidefinite relaxation
val,model,_ = pcpop!(p,2;min=false) 
println("Termination status ", termination_status(model))
println("Optimal value is   ", val)
sos_model= pcpop(p, 2)
set_optimizer(sos_model, Mosek.Optimizer)
optimize!(sos_model)
println("Termination status ", termination_status(sos_model))
println("Optimal value is   ", objective_value(sos_model))