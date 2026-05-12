include("traceGrobner.jl")

println("Building monoid...")
# Build the monoid
@pcmonoid M a[16,0]
Projector.(a)
A = reshape(a, 4, 4)
for i in 1:4
    @comms A[i, :]
    @comms A[:, i]
end
build(M)
# Objective function
p = A[2,2] + A[1,3] + A[3,1] + A[1,1]
p+= A[2,4] + A[4,2] + A[4,3] + A[3,4]
p+=-A[1,2] - A[1,4] - A[2,1] - A[2,3]
p+=-A[3,2] - A[3,3] - A[4,1] - A[4,4]
# Constraints
R = []
for i in 1:4
    append!(R, [one(M) - sum(A[i,:])])
    append!(R, [one(M) - sum(A[:,i])])
end
append!(R, [one(M) - sum(A[1:2,1:2])])
append!(R, [one(M) - sum(A[1:2,3:4])])
append!(R, [one(M) - sum(A[3:4,1:2])])
append!(R, [one(M) - sum(A[3:4,3:4])])

# Semidefinite relaxation
println("Building model...")
op_eq=R
val,model,_,_ = npa(p, 2, min=false, op_eq=R)
println("Termination status ", termination_status(model))
println("Optimal value is   ", objective_value(model))

println("Building model...")
model = pcpop(p, 2, equalities=R)
set_optimizer(model, Mosek.Optimizer)
println("Optimizing...")
optimize!(model)
println("Termination status ", termination_status(model))
println("Optimal value is   ", objective_value(model))