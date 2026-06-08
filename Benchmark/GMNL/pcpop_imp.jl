include("../../traceGrobner.jl")
using Combinatorics

function gmnl(n,k; primal=true, canonical=true)

    @pcmonoid M A[2*n,0]
    # Set variables to projectors
    Projector.(M.vertices)
    # Set commutation relations
    A = reshape(A, n,2)
    for (i, j) in combinations(1:n, 2)
        @comms A[i, :] A[j, :]
    end
    # Build the monoid
    build(M)
    start_setup = time()
    # Objective function
    obj = prod(A[:,1])
    # Linear constraints on the moments
    tr_eq = [[A[i,2]*A[(i+1),1],0] for i in 1:(n-1)]
    push!(tr_eq,[A[n,2]*A[1,1],0])
    push!(tr_eq,[prod([1-A[i,2] for i in 1:n]),0])
    # Optimization of the semidefinite relaxation
    model,_ = pcpop!(obj, k; tr_eq=tr_eq, min=false,optimize=false,primal=primal,canonical=canonical)
    stop_setup = time()
    set_optimizer(model, Mosek.Optimizer)
    start_solve = time()
    optimize!(model)
    stop_solve = time()
    elapsed_setup = stop_setup - start_setup
    elapsed_solve = stop_solve - start_solve
    return elapsed_setup, elapsed_solve
end

function avg_time(total_runs,n,k;primal=true,canonical=true)
    # Hot start
    gmnl(n,k; primal=primal, canonical=canonical)
    gmnl(n,k; primal=primal, canonical=canonical)

    total_setup_time = 0.0
    total_solve_time = 0.0
    for i in 1:total_runs
        setup_time, solve_time = gmnl(n,k; primal=primal, canonical=canonical)
        total_setup_time += setup_time
        total_solve_time += solve_time
    end
    avg_setup_time = total_setup_time / total_runs
    avg_solve_time = total_solve_time / total_runs
    return avg_setup_time, avg_solve_time
end