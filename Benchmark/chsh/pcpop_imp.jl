include("../../traceGrobner.jl")

function chsh(k;primal=false,canonical=true)
    @pcmonoid M a[2,0] b[2,0]
    Unipotent.(M.vertices)
    @comms a b
    build(M)
    start_setup = time()

    # Objective function
    obj = a[1]*b[1] + a[1]*b[2] + a[2]*b[1] - a[2]*b[2]
    # Optimization of the semidefinite relaxation
    model,_ = pcpop!(obj, k; min=false,optimize=false,primal=primal,canonical=canonical)
    stop_setup = time()
    start_solve = time()
    set_optimizer(model, Mosek.Optimizer)
    optimize!(model)
    stop_solve = time()
    elapsed_setup = stop_setup - start_setup
    elapsed_solve = stop_solve - start_solve
    return elapsed_setup, elapsed_solve
end

function avg_time(total_runs,k;primal=true,canonical=true)
    # hot run
    chsh(k; primal=primal, canonical=canonical)
    chsh(k; primal=primal, canonical=canonical)

    # Actual timer starts here
    total_setup_time = 0.0
    total_solve_time = 0.0
    for i in 1:total_runs
        setup_time, solve_time = chsh(k; primal=primal, canonical=canonical)
        total_setup_time += setup_time
        total_solve_time += solve_time
    end
    avg_setup_time = total_setup_time / total_runs
    avg_solve_time = total_solve_time / total_runs
    return avg_setup_time, avg_solve_time
end