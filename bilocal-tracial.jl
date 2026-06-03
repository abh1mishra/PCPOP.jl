include("traceGrobner.jl")

println("Building monoid...")
# Build the monoid
@pcmonoid M ρ[2,0] a[2,0] b[2,0] c[2,0]
Unipotent.(M.vertices)
@comms a b c
@comms a ρ[2]
@comms c ρ[1]
@comms ρ
build(M)

# Objective function
r = ρ[1]*ρ[2]
α = r*sum(a[i]*b[1]*c[j] for i in 1:2 for j in 1:2)
γ = r*sum(a[i]*b[2]*c[j]*(-1)^(i+j) for i in 1:2 for j in 1:2)
p = (α - γ)^2/8 - (α + γ)
# Inequality constraints
S = [ρ[1] - ρ[1]^2,
     ρ[2] - ρ[2]^2]
# Moment constraints
T = [(ρ[1], 1),
     (ρ[2], 1)]
println("Building model...")
k = 5
model = npa(p, k, 
op_ge=S, 
tr_eq=T, 
normalize=false, 
cyclic=true,
lvl_lm = 1)
set_optimizer(model, Mosek.Optimizer)
println("Optimizing...")
optimize!(model)
println("Termination status ", termination_status(model))
println("Optimal value is   ", objective_value(model))
