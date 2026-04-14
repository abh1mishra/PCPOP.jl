# Initialize local monoids
@ncmonoid A a1 a2
@ncmonoid B b1 b2
@ncmonoid AB x1 x2
vars = [a1,a2,b1,b2,x1,x2]
Projector.(vars)
# Build global monoid
M = GraphProductMonoid("M",[A, B, AB])
@comms A B
build(M)
# Objective function.
obj  = a1*a2 + a2*a1 + b1*b2 + b1*b1 + x1*x2 + x2*x1 
obj += a1*b2 + a1*x2 + b1*a2 + b1*x2 + x1*a2 + x1*b2 
# Inequality constraints
op_ge=[v-v^2+0.5 for v in vars]
# Optimize semidefinite relaxation
val,_ = npa(obj,1;op_ge=op_ge) 
println("Optimal value is ", val)