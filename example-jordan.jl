include("traceGrobner.jl")

# CHSH setting
@pcmonoid M a b c d
@comms [a, b] [c, d]
Unipotent.([a,b,c,d])
build(M)

f = a*c + a*d + b*c - b*d

# Optimimzation
k = 1
diagonalize=true
Γ, C, A, b  = npa_dual(f, k, rm=true)
model, P, blkD = jordan_reduce(C, A, b, verbose=true, complex=true, diagonalize=diagonalize)
println("Termination status: ", termination_status(model))
println("Objective value: ", objective_value(model))

# # Example statepop 7.2.1
# @pcmonoid M a[2,0] b[2,0]
# Unipotent.(a)
# Unipotent.(b)
# @comms a b
# build(M)

# TM = make_trace_monoid(M, 6, tracial=false)
# p  = (state(a[1]*b[2], TM) + state(a[2]*b[1], TM))^2 
# p += (state(a[1]*b[1], TM) - state(a[2]*b[2], TM))^2
# basis = trace_monomials(TM, 0:3, tracial=false)
# model = tpop(p, TM, basis, tracial=false)

# # Jordan reduction
# C, A, b = write_canonical(model)
# model_red, P, blkD = jordan_reduce(C, A, b; verbose=true)

# # Compare optimal solutions
# set_optimizer(model, Mosek.Optimizer)
# set_silent!(model)
# optimize!(model)
# println(termination_status(model))
# println(objective_value(model))

# set_optimizer(model_red, Mosek.Optimizer)
# set_silent!(model_red)
# optimize!(model_red)
# println(termination_status(model_red))
# println(objective_value(model_red))