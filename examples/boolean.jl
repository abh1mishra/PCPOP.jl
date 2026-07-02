using PCPOP, JuMP

# Build the monoid
@pcmonoid M a b
Projector.([a, b])
@comms a b
build(M)
# Formula (a ∨ b) ∧ (a ∨ ¬b) ∧ (¬a ∨ b) ∧ (¬a ∨ ¬b)
p = (a+b)*(1+a-b)*(1-a+b)*(2-a-b)
# Optimize semidefinite relaxation
val, model, _ = pcpop(p, 1; min = false, list_vars = M.vertices)
println("Termination status ", termination_status(model))
println("Optimal value      ", val)
