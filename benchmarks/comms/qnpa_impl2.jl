using QuantumNPA,MosekTools,JuMP,Graphs

function party_vec_gen(n)
    party_D = Dict([])
    G=SimpleGraph(n)
    [add_edge!(G, i, mod(i,n)+1) for i in 1:n]
    mc = maximal_cliques(complement(G))
    for i in 1:n
        party_D[i] = [j for j in 1:length(mc) if i in mc[j]]
    end
    return party_D
end

function polygon_bell(n,d;primal=false,optimize=true)
    setup_start = time()
    party_D = party_vec_gen(n)
    V=[]
    party_n = Dict([v=>0 for (k,v) in party_D])
    for i in 1:n
        st = party_n[party_D[i]]+1
        en = party_n[party_D[i]]+2
        push!(V, projector(party_D[i], 1, st:en))
        party_n[party_D[i]] += 2
    end
    # Objective function (projector form: dichotomic = 1 - 2*projector)
    obj = sum([let (a, b) = (V[i], V[i+1]); (Id-2*a[1])*(Id-2*b[1]) + (Id-2*a[1])*(Id-2*b[2]) + (Id-2*a[2])*(Id-2*b[1]) - (Id-2*a[2])*(Id-2*b[2]) end for i in 1:(n-1)])
    a=V[n];b=V[1];obj += (Id-2*a[1])*(Id-2*b[1]) + (Id-2*a[1])*(Id-2*b[2]) + (Id-2*a[2])*(Id-2*b[1]) - (Id-2*a[2])*(Id-2*b[2])
    # Optimization of the semidefinite relaxation
    if !primal
        model = npa2jump_d(obj, d; sense=:maximize,solver=Mosek.Optimizer)
    else
        model = npa2jump(obj, d; sense=:maximize,solver=Mosek.Optimizer)
    end
    stop_setup = time()
    if !optimize
        return nothing, stop_setup - setup_start, 0.0
    end
    start_solve = time()
    set_optimizer(model, Mosek.Optimizer)
    unset_silent(model)
    optimize!(model)
    stop_solve = time()
    if termination_status(model) != MOI.OPTIMAL
        @warn "Optimization terminated with status: $(termination_status(model))"
    end
    elapsed_setup = stop_setup - setup_start
    elapsed_solve = stop_solve - start_solve
    return objective_value(model), elapsed_setup, elapsed_solve
end

function avg_time(total_runs,n,k;primal=false,optimize=true)
    # hot run
    polygon_bell(n,1; primal=primal)
    polygon_bell(n,1; primal=primal)

    # Actual timer starts here
    total_setup_time = 0.0
    total_solve_time = 0.0
    for i in 1:total_runs
        _,setup_time, solve_time = polygon_bell(n,k; primal=primal,optimize=optimize)
        total_setup_time += setup_time
        total_solve_time += solve_time
    end
    avg_setup_time = total_setup_time / total_runs
    avg_solve_time = total_solve_time / total_runs
    return avg_setup_time, avg_solve_time
end