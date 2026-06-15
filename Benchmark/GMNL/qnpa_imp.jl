using QuantumNPA,MosekTools,JuMP

function gmnl(n,k;primal=false,optimize=true)
    start_setup = time()
    A = [projector([i],1,1:2) for i in 1:n]
    A = reduce(vcat, [v' for v in A])
    # Objective function
    obj = prod(A[:,1])
    # Linear constraints on the moments
    tr_eq = Vector{Any}([A[i,2]*A[(i+1)] for i in 1:(n-1)])
    push!(tr_eq,A[n,2]*A[1,1])
    push!(tr_eq,prod([Id-A[i,2] for i in 1:n]))
    # Level of semidefinite relaxation
    if !primal
    model = npa2jump_d(obj, k,eq=tr_eq; sense=:maximize,solver=Mosek.Optimizer)
    else
        model = npa2jump(obj, k,eq=tr_eq; sense=:maximize,solver=Mosek.Optimizer)
    end
    stop_setup = time()
    if optimize
        start_solve = time()
        optimize!(model)
        stop_solve = time()
        elapsed_setup = stop_setup - start_setup
        elapsed_solve = stop_solve - start_solve
        return elapsed_setup, elapsed_solve
    else
        elapsed_setup = stop_setup - start_setup
        return elapsed_setup,0.0
    end
end

function avg_time(total_runs,n,k;optimize=true,primal=false)
    total_setup_time = 0.0
    total_solve_time = 0.0
    gmnl(1,1; optimize=false, primal=primal)
    gmnl(1,1; optimize=false, primal=primal)

    for i in 1:total_runs
        setup_time, solve_time = gmnl(n,k; optimize=optimize, primal=primal)
        total_setup_time += setup_time
        total_solve_time += solve_time
    end
    avg_setup_time = total_setup_time / total_runs
    avg_solve_time = total_solve_time / total_runs
    return avg_setup_time, avg_solve_time
end