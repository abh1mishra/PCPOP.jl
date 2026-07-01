function exact_solve(model::JuMP.Model, filename::AbstractString; T=Float64)
    # Write JuMP model to SDPA-sparse .dat-s
    file = filename * ".dat-s"
    JuMP.write_to_file(model, file)
    # Read SDPA-sparse into a CLRS model
    read_model = ClusteredLowRankSolver.sdpa_sparse_to_problem(file; T=T)
    # Solve with ClusteredLowRankSolver
    model_clrs = ClusteredLowRankSolver.solvesdp(read_model)

    # round the solution
    FF, g = find_field(model_clrs[2], model_clrs[3])
    status, problem, esol = exact_solution(read_model, model_clrs[2], model_clrs[3]; FF, g)

    return problem, esol
end

function npa_exact(obj, level;
    min=true,
    op_eq = 0, 
    op_ge = 0,
    tr_eq = 0,
    tr_ge = 0,
    lvl_lm=0,
    list_vars=[],
    cyclic=false,
    normalize=true,
    optimizer=Mosek.Optimizer,
    model_flags=[],
    rm=false)

    model = npa(obj, level; min=min;
                            op_eq = op_eq, 
                            op_ge = op_ge,
                            tr_eq = tr_eq,
                            tr_ge = tr_ge,
                            lvl_lm=lvl_lm,
                            list_vars=list_vars,
                            cyclic=cyclic,
                            normalize=normalize,
                            optimizer=optimizer,
                            model_flags=model_flags,
                            rm=rm)

    problem, esol = exact_solve(model, "saved_model") 
    
end