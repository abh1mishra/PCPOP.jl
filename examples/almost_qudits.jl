include("../traceGrobner.jl")

#Parameters
d = 2
ϵ = 0.01
c = 2*(2+sqrt(2))

# Build monoid
@pcmonoid M ρ[4,0] B[2,0] P[1,0]
Projector.(ρ)
Projector.(B)
Projector.(P)
build(M)

# Objective function
obj = ρ[1]*B[1]
# Random access code
rac = ρ[1]*B[1] + ρ[2]*B[1] + ρ[3]*(1-B[1]) + ρ[4]*(1-B[1])
rac+= ρ[1]*B[2] + ρ[2]*(1-B[2]) + ρ[3]*B[2] + ρ[4]*(1-B[2])
# Linear equalities on the moments
tr_eq = [[ρ[1], 1],
         [ρ[2], 1],
         [ρ[3], 1],
         [ρ[4], 1],
         [P[1], d],
         [rac, c]]
# Linear inequalities on the moments
tr_ge = [[ρ[1]*P[1], 1 - ϵ],
         [ρ[2]*P[1], 1 - ϵ],
         [ρ[3]*P[1], 1 - ϵ],
         [ρ[4]*P[1], 1 - ϵ]]

# Level of semidefinite relaxation
k = 3
# Optimization of the semidefinite relaxation
val,model,_= pcpop!(obj,k;tr_eq=tr_eq, tr_ge=tr_ge, min=false, tracial=true, normalize=false)
println("Termination status ", termination_status(model))
println("Optimal value      ", val)
