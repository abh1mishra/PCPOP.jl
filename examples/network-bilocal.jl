using PCPOP, JuMP

# Build the monoid
@pcmonoid M ρ[2, 0] a[2, 0] b[2, 0] c[2, 0]
Unipotent.([a b c])
Projector.(ρ)
@comms a b c
@comms ρ
@comms ρ[1] c
@comms ρ[2] a
#@comms M.vertices...
build(M)
# Objective function.
σ = ρ[1]*ρ[2]
p = σ*(a[2]*b[1]*c[1] + a[1]*b[2]*c[1] + a[1]*b[1]*c[2] - a[2]*b[2]*c[2])
# Constraints
S = []

d = 8

T = [[σ, 1], [ρ[1], 1], [ρ[2], 1]]

U = [[-ρ[1], -d], [-ρ[2], -d]]
U = []

# Optimize semidefinite relaxation
val, model, _ = pcpop!(
    p,
    3,
    min = false,
    op_ge = S,
    tr_eq = T,
    tr_ge = U,
    normalize = false,
    tracial = true,
    lvl_lm = 1,
)
println("Termination status ", termination_status(model))
println("Optimal value is   ", val)
