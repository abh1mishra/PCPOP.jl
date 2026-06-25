include("../../traceGrobner.jl")

function polygon_bell(n,d;primal=false,canonical=true)
    setup_start = time()
    @pcmonoid M V[2*n,0]
    Projector.(M.vertices)
    V = reshape(V, n, 2)
    for i in 1:n-1
        @comms V[i,:] V[i+1,:]
    end
    @comms V[n,:] V[1,:]
    build(M)
    # Objective function (projector form: dichotomic = 1 - 2*projector)
    obj = sum([let (a, b) = (V[i,:], V[i+1,:]); (1-2*a[1])*(1-2*b[1]) + (1-2*a[1])*(1-2*b[2]) + (1-2*a[2])*(1-2*b[1]) - (1-2*a[2])*(1-2*b[2]) end for i in 1:(n-1)])
    obj += (1-2*V[n,1])*(1-2*V[1,1]) + (1-2*V[n,1])*(1-2*V[1,2]) + (1-2*V[n,2])*(1-2*V[1,1]) - (1-2*V[n,2])*(1-2*V[1,2])
    # Optimization of the semidefinite relaxation
    model,_ = pcpop!(obj, d; min=false,optimize=false,primal=primal,canonical=canonical)
    stop_setup = time()
    start_solve = time()
    set_optimizer(model, Mosek.Optimizer)
    # set_silent(model)
    optimize!(model)
    stop_solve = time()

    if termination_status(model) != MOI.OPTIMAL
        error("Optimization failed with status: $(termination_status(model))")
    end
    elapsed_setup = stop_setup - setup_start
    elapsed_solve = stop_solve - start_solve
    return objective_value(model), elapsed_setup, elapsed_solve
end

function avg_time(total_runs,n,k;primal=false,canonical=true)
    # hot run
    polygon_bell(n,1; primal=primal, canonical=canonical)
    polygon_bell(n,1; primal=primal, canonical=canonical)

    # Actual timer starts here
    total_setup_time = 0.0
    total_solve_time = 0.0
    for i in 1:total_runs
        _,setup_time, solve_time = polygon_bell(n,k; primal=primal, canonical=canonical)
        total_setup_time += setup_time
        total_solve_time += solve_time
    end
    avg_setup_time = total_setup_time / total_runs
    avg_solve_time = total_solve_time / total_runs
    return avg_setup_time, avg_solve_time
end