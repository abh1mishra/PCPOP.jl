# Build the monoid
@pcmonoid M ρ[3,0] a b c
Projector.([a, b, c])
@comms a b c
@comms ρ
@comms ρ[1] c
@comms ρ[2] a
@comms ρ[3] b
#@comms M.vertices...
build(M)
# Objective function.
r = ρ[1]*ρ[2]*ρ[3]
p = r*((1-a)*b*c + a*(1-b)*c + a*b*(1-c) - (1-a)*(1-b)*(1-c))
# Constraints
S = [ρ[1] - ρ[1]^2,
	 ρ[2] - ρ[2]^2,
	 ρ[3] - ρ[3]^2]
T = [[ρ[1]-ρ[2], 0],
	 [ρ[2]-ρ[3], 0],
	 [ρ[3]-ρ[1], 0],]
T = [[ρ[1], 1],
	 [ρ[2], 1],
	 [ρ[3], 1]]
# Optimize semidefinite relaxation
val,model,_ = npa(p,2, min=false,
					   op_ge = S,
					   tr_eq = T,
					   normalize=false,
					   cyclic=true) 
println("Termination status ", termination_status(model))
println("Optimal value is   ", val)

println("Maximal classical value: 2")
println("Maximal quantum value: [2√2, 3.085]")
