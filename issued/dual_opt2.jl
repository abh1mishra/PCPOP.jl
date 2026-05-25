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

function upscale_poly_mat(p::Polynomial, tsize::Int, offset::Int)
    new_coeffs = typeof(p.coeffs[1])[]
    
    for c in p.coeffs
        n, m = size(c)
        # Create an inflated sparse matrix
        new_mat = spzeros(eltype(c), tsize, tsize)
        new_mat[offset+1:offset+n, offset+1:offset+m] = c
        push!(new_coeffs, new_mat)
    end

    # Return a new polynomial with identical monomials but inflated coefficient matrices
    return Polynomial(p.monomials, new_coeffs, p.monoid)
end


function npa_moment(operators::Vector;cPoly=1)
    N = length(operators)
    iops = collect(enumerate(operators))
    monoid = first(operators).monoid
    moment_mat = Polynomial{Matrix{Float64}}(monoid)
    unique_mons=[]
    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(real_rep(conj(x)*cPoly*y))

            for (m,c) in p
                m_i=findfirst(x->(x==m),unique_mons)
                if m_i !== nothing
                    sym_add!(Main.coefficient(moment_mat,m), i, j, c)
                else
                    push!(moment_mat.monomials, m)
                    push!(moment_mat.coeffs, sym_add!(spzeros(N, N), i, j, c))
                    push!(unique_mons,m)
                end
            end
        end
    end

    return moment_mat
end


function npa_dual(obj, level;
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

    # Delete redundant polynomials(polynomials which are numbers, so k*Id) from the op_ge and op_eq and warn for incompatible polynomials in op_ge and op_eq
    op_ge!=0 && (op_ge=sanity_check_op_ge(op_ge))
    op_eq!=0 && (op_eq=sanity_check_op_eq(op_eq))
    if all(is_number.(vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])))
        @warn "All the input polynomials are constants. The optimization will be trivial."
        isempty(list_vars) && throw(ArgumentError("The list of variables is empty. Please provide a non-empty list of variables."))
    end
    # at this pont, op_ge and op_eq constaints non-trivial polynomials or zero.
    if lvl_lm==0
        ops, ops_principal = get_monomials(obj,level; op_eq = op_eq, op_ge = op_ge, tr_eq = tr_eq, tr_ge = tr_ge,list_vars=list_vars)
    else
        if !isempty(list_vars)
            ops_principal= mons_at_level(list_vars, level)
            ops= mons_at_level(list_vars, lvl_lm)
        else
            pols=vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])
            vars=unique_array(union([variables(g) for g in pols]...))
            ops_principal = mons_at_level(vars, level)
            ops = mons_at_level(vars, lvl_lm)
        end
    end
    println("Number of operators in the principal moment matrix: ", length(ops_principal))
    model=Model(optimizer)
    for (flag,val) in model_flags
        set_optimizer_attribute(model,flag,val)
    end

    Id = one(first(ops_principal))
    pair_mat_ops = [(Id,ops_principal)]
    if op_ge != 0
        pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_ge])
    end

    if op_eq != 0
        pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_eq])
        pair_mat_ops=vcat(pair_mat_ops,[(-g,ops) for g in op_eq])
    end


    # if normalize
    #     push!(tr_eq, (Id, 1.0))
    # end

    if tr_ge != 0
        for (g, val) in tr_ge
            tr_ge_poly = g-val
            pair_mat_ops=vcat(pair_mat_ops,[(tr_eq_poly,[Id])])
        end
    end

    if tr_eq != 0
        for (g, val) in tr_eq
            tr_eq_poly = g-val
            pair_mat_ops=vcat(pair_mat_ops,[(tr_eq_poly,[Id])])
            pair_mat_ops=vcat(pair_mat_ops,[(-tr_eq_poly,[Id])])
        end
    end
    println("LMI begin")
    LMI = [npa_moment(ops_i;cPoly=g) for (g, ops_i) in pair_mat_ops]
    tsize = sum(length(i) for (j,i) in pair_mat_ops)
    upscaled_LMI = []
    offset = 0
    for i in 1:length(LMI)
        push!(upscaled_LMI,upscale_poly_mat(LMI[i], tsize,offset))
        offset += length(pair_mat_ops[i][2])
    end
    println("F begin")

    F = sum(upscaled_LMI)
    Z=@variable(model, [1:tsize, 1:tsize], PSD)
    println("F ends")

    obj_poly=real_rep(Polynomial(obj))
    if !is_number(obj)
        obj_p=0
        if cyclic
            for (m,c) in Polynomial(obj)
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i==nothing
                    throw(ArgumentError("level not enough"))
                end
                obj_p+=c*unique_vars[m_i]
            end
        else
            obj_p=Main.coefficient(obj_poly,Id) + LinearAlgebra.dot(Z,Main.coefficient(F,Id))
        end
        
        min ? @objective(model, Max, obj_p) : @objective(model, Min, obj_p)
    end
    println("Adding constraints")
    mons=monomials(F)
    for (i,m) in enumerate(mons)
        println("Adding constraint for monomial ", i/length(mons))
        if m!=Id
            constr = LinearAlgebra.dot(Z,Main.coefficient(F,m)) + Main.coefficient(obj_poly,m)
            @constraint(model, constr == 0)
        end
    end
    println("Constraints added")
    if rm
        return model, F, Z
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

    return optimal_value, model, F,Z

end
