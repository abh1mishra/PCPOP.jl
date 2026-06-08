# Taken from QuantumNPA.jl with some modifications

function LinearAlgebra.dot(A::SparseMatrixCSC,
                           B::Symmetric{<:JuMP._MA.AbstractMutable,
                                        Matrix{<:JuMP._MA.AbstractMutable}})
    acc = zero(eltype(B))

    for j in 1:size(A, 2)
        for k in nzrange(A, j)
            add_to_expression!(acc, nonzeros(A)[k], B[rowvals(A)[k], j])
        end
    end

    return acc
end

function LinearAlgebra.dot(A::Symmetric{<:JuMP._MA.AbstractMutable,
                                        Matrix{<:JuMP._MA.AbstractMutable}},
                           B::SparseMatrixCSC)
    return dot(B, A)
end

function sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end

    return matrix
end


function npa_moment(operators::Vector;cPoly=1)
    N = length(operators)
    iops = collect(enumerate(operators))
    moment_D = Dict{AbstractMonomial,Matrix{Float64}}([])
    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(real_rep(conj(x)*cPoly*y))
            for (m,c) in p
                if haskey(moment_D, m)
                    sym_add!(moment_D[m], i, j, c)
                else
                    moment_D[m]= sym_add!(spzeros(N, N), i, j, c)
                end
            end
        end
    end

    return moment_D
end


function cyclic_npa_moment(operators::Vector;cPoly=1)
    N = length(operators)
    iops = collect(enumerate(operators))
    moment_D = Dict{AbstractMonomial,Matrix{Float64}}([])
    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(conj(x)*cPoly*y)
            for (m,c) in p
                if m isa Variable
                    println("error ",x, "  ", y, "  ", cPoly)
                end
                m_,m__ = cyclic_reduce(m), cyclic_reduce(m') 
                if haskey(moment_D, m_)
                    sym_add!(moment_D[m_], i, j, c)
                elseif haskey(moment_D, m__)
                    sym_add!(moment_D[m__], i, j, c)
                else
                    moment_D[m_] = sym_add!(spzeros(N, N), i, j, c)
                end
            end
        end
    end

    return moment_D
end

function npa_dual(obj, level;
    min=false,
    op_eq = [], 
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    lvl_lm=-1,
    list_vars=[],
    tracial=false,
    normalize=true,
    optimizer=Mosek.Optimizer,
    model_flags=[],
    rm=false)

    ops,ops_principal = basis_gen(obj,level,op_eq,op_ge,tr_eq,tr_ge,list_vars,lvl_lm)
    println("Number of operators in the principal moment matrix: ", length(ops_principal))
    model=Model(optimizer)
    for (flag,val) in model_flags
        set_optimizer_attribute(model,flag,val)
    end

    Id = one(first(ops_principal))

    pair_mat_ops = [(Id,ops_principal)]

    pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_ge])

    pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_eq])
    pair_mat_ops=vcat(pair_mat_ops,[(-g,ops) for g in op_eq])


    if !isempty(tr_ge)
        for (g, val) in tr_ge
            tr_ge_poly = g-val
            pair_mat_ops=vcat(pair_mat_ops,[(tr_ge_poly,[Id])])
        end
    end

    if !isempty(tr_eq)
        for (g, val) in tr_eq
            tr_eq_poly = g-val
            pair_mat_ops=vcat(pair_mat_ops,[(tr_eq_poly,[Id])])
            pair_mat_ops=vcat(pair_mat_ops,[(-tr_eq_poly,[Id])])
        end
    end

    LMI = [tracial ? cyclic_npa_moment(ops_i;cPoly=g) : npa_moment(ops_i;cPoly=g) for (g, ops_i) in pair_mat_ops]
    tsize = [length(i) for (j,i) in pair_mat_ops]

    Zs=[@variable(model, [1:l, 1:l], PSD) for l in tsize]


    min ? s=-1 : s=1
    tracial ? obj_poly = Polynomial(obj) : obj_poly = real_rep(Polynomial(obj))
    tracial ? Id = cyclic_reduce(Id) : Id = Id
    if !is_number(obj)

        obj_p = Main.coefficient(obj_poly,Id) + sum(s*LinearAlgebra.dot(Zs[i],LMI[i][Id]) for i in 1:length(LMI))        
        min ? @objective(model, Max, obj_p) : @objective(model, Min, obj_p)
    end
    println("Adding constraints")
    mons = union!(Set(), [keys(lmi) for lmi in LMI]...)
    for m in mons
        if m != Id
            constr = sum((s*LinearAlgebra.dot(Zs[i], LMI[i][m]) for i in 1:length(LMI) if haskey(LMI[i], m)), init=0.0) + Main.coefficient(obj_poly, m)
            @constraint(model, constr == 0)
        end
    end
    println("Constraints added")
    if rm
        return model, LMI, Zs
    end
    # Solve the model
    optimize!(model)
    # Check the optimization status
    if termination_status(model) != MOI.OPTIMAL
        @warn "The optimization problem is $(termination_status(model))."
    end
    # Extract the optimal value of the objective function
    optimal_value = objective_value(model)

    # Return the optimal value and the dictionary of optimal variables
    # return optimal_value, optimal_vars, model

    if is_number(obj)
        @warn "The objective function is a constant, it is a feasibility check"
        return termination_status(model),model,ops_principal,unique_mons,PM
    end

    return optimal_value, model, LMI,Zs

end
