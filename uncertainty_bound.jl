# Build base monoid in variables x, y, z
@pcmonoid M x y z
Unipotent.([x, y, z])
build(M)
# Build state monoid over M
TM = make_trace_monoid(M, 2, tracial=false)
# Objective function
ρx = state(x, TM)
ρy = state(y, TM)
ρz = state(z, TM)
p = ρx^2 + ρy^2 + ρz^2
# Anti-commutation relations
R = [μx*μy + μy*μx,
	 μy*μz + μz*μy,
	 μz*μx + μx*μz]
add_relations!(R)
# Optimize semidefinite relaxation
basis = [one(TM.state_monoid), μx*ρx, μy*ρy, μz*ρz]
sos_model = tpop(p, TM, basis, tracial=false)
set_optimizer(sos_model, Mosek.Optimizer)
optimize!(sos_model)
println("Termination status ", termination_status(sos_model))
println("Optimal value is   ", objective_value(sos_model))