using QuantumNPA,MosekTools,JuMP

function chsh(k;primal=false)
    start_setup = time()
    @dichotomic A1 A2 B1 B2;
    # Objective function
    S = A1*(B1 + B2) + A2*(B1 - B2)
    if !primal
        model = npa2jump_d(S, k; sense=:maximize,solver=Mosek.Optimizer)
    else
        model = npa2jump(S, k; sense=:maximize,solver=Mosek.Optimizer)
    end
    stop_setup = time()
    start_solve = time()
    optimize!(model)
    stop_solve = time()
    elapsed_setup = stop_setup - start_setup
    elapsed_solve = stop_solve - start_solve
    return elapsed_setup, elapsed_solve
end

function avg_time(total_runs,k;primal=false)
    # hot run
    chsh(k; primal=primal)
    chsh(k; primal=primal)
    
    # Actual timer starts here
    total_setup_time = 0.0
    total_solve_time = 0.0
    for i in 1:total_runs
        setup_time, solve_time = chsh(k; primal=primal)
        total_setup_time += setup_time
        total_solve_time += solve_time
    end
    avg_setup_time = total_setup_time / total_runs
    avg_solve_time = total_solve_time / total_runs
    return avg_setup_time, avg_solve_time
end