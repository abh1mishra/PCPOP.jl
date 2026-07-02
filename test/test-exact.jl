using PCPOP
include("src/exact.jl")

@pcmonoid M a[2, 0] b[2, 0]
# Set variables to unitaries
Unipotent.([a; b])
# Set commutation relations
@comms a b
# Build the monoid
build(M)
# Objective function
p = a[1]*b[1] + a[1]*b[2] + a[2]*b[1] - a[2]*b[2]
# Level of semidefinite relaxation
level=1
# Optimization of the semidefinite relaxation
val, model, vars, moments = npa_exact(p, level)

exact_solve(model, "save-model")
