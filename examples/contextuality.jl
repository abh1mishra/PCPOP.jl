using PCPOP, JuMP
# Build the monoid
@pcmonoid M X[9, 0]
Unipotent.(X)
x = reshape(X, (3, 3))
for i in 1:3
    @comms x[i, 1] x[i, 2] x[i, 3]
    @comms x[1, i] x[2, i] x[3, i]
end
build(M)
# Conditions on the game
R = [
    x[1, 1]*x[1, 2]*x[1, 3] - 1,
    x[2, 1]*x[2, 2]*x[2, 3] - 1,
    x[3, 1]*x[3, 2]*x[3, 3] - 1,
    x[1, 1]*x[2, 1]*x[3, 1] + 1,
    x[1, 2]*x[2, 2]*x[3, 2] + 1,
    x[1, 3]*x[2, 3]*x[3, 3] + 1,
]
# Optimize semidefinite relaxation
val, model, _ = pcpop(0, 2; op_eq = R, lvl_lm = 0)
println("Termination status ", termination_status(model))
