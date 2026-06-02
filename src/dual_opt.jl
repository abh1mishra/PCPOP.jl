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

function cyclic_npa_moments_block_dual!(list_monomials::Vector{M},A,B,tsize;cPoly=1,unique_mons=[],unique_pos=[],offset=0,extra_zeros=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials
    
            # Initialize the row of the A matrix for constraints
            Ai=spzeros(tsize * tsize)
            
            if cPoly == 1
                m = list_monomials[i]'* list_monomials[j]
                m1,m2 = cyclic_reduce(m),cyclic_reduce(m')
                if m1==0 || m2==0
                    Ai[(i-1)*tsize + j] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    Ai[(i-1)*tsize + j] = 1.0
                    upi,upj = unique_pos[m_i]
                    Ai[(upi-1)*tsize + upj] = -1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                else
                    push!(unique_mons, m1)
                    push!(unique_pos, (i,j))
                end
            else
                monomial_product = list_monomials[i]'* cPoly * list_monomials[j]
                if monomial_product == 0
                    Ai[(offset+i-1)*tsize+offset+j] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                Ai[(offset+i-1)*tsize + offset+j] = 1.0
                for (m,c) in monomial_product
                    m1,m2 = cyclic_reduce(m),cyclic_reduce(m')
                    mi = findfirst(x->(x==m1 || x==m2),unique_mons)
                    upi,upj = unique_pos[mi]
                    Ai[(upi-1)*tsize + upj] = -c
                end
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                Ai=spzeros(tsize * tsize)
                Ai[(offset+i-1)*tsize+ j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end

        for i in 1:num_monomials
            for j in offset+num_monomials+1:tsize
                Ai=spzeros(tsize * tsize)
                Ai[(offset+i-1)*tsize + j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if cPoly == 1
        return unique_mons,unique_pos
    end

end

function npa_moments_block_dual!(list_monomials::Vector{M},A,B,tsize;cPoly=1,unique_mons=[],unique_pos=[],offset=0,extra_zeros=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials

            # Initialize the row of the A matrix for constraints
            Ai=spzeros(tsize*tsize)

            if cPoly == 1
                m = real_rep(list_monomials[i]'* list_monomials[j])
                if m == 0
                    Ai[(i-1)*tsize + j] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                m_i=findfirst(x->(x==m),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    Ai[(i-1)*tsize + j] = 1.0
                    upi,upj = unique_pos[m_i]
                    Ai[(upi-1)*tsize + upj] = -1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                else
                    push!(unique_mons, m)
                    push!(unique_pos, (i,j))
                end
            else
                monomial_product = real_rep(Polynomial(list_monomials[i]'* cPoly * list_monomials[j]))
                # assume the monomials PM exists and m and m' are in pm.
                if monomial_product == 0
                    Ai[(offset+i-1)*tsize + (offset+j)] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                Ai[(offset+i-1)*tsize + (offset+j)] = 1.0
                for (m,c) in monomial_product
                    m_i=findfirst(x->(x==m),unique_mons)
                    upi,upj = unique_pos[m_i]
                    Ai[(upi-1)*tsize + upj] = -c
                end
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                Ai=spzeros(tsize*tsize)
                Ai[(offset+i-1)*tsize + j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end

        for i in 1:num_monomials
            for j in offset+num_monomials+1:tsize
                Ai=spzeros(tsize*tsize)
                Ai[(offset+i-1)*tsize + j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if cPoly == 1
        return unique_mons,unique_pos
    end
end
function npa_canonical(obj, level;
    min=true,
    op_eq = [], 
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    lvl_lm=-1,
    list_vars=[],
    cyclic=false,
    normalize=true,
    optimizer=Mosek.Optimizer,
    model_flags=[],
    rm=false,
    extra_zeros=false,
    )

    # Delete redundant polynomials(polynomials which are numbers, so k*Id) from the op_ge and op_eq and warn for incompatible polynomials in op_ge and op_eq
    !isempty(op_ge) && (op_ge=sanity_check_op_ge(op_ge))
    !isempty(op_eq) && (op_eq=sanity_check_op_eq(op_eq))
    if all(is_number.(vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])))
        @warn "All the input polynomials are constants. The optimization will be trivial."
        isempty(list_vars) && throw(ArgumentError("The list of variables is empty. Please provide a non-empty list of variables."))
    end
    # at this pont, op_ge and op_eq constaints non-trivial polynomials or zero.
    if lvl_lm==-1
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
    A=Vector{SparseMatrixCSC{Float64, Int64}}([])
    B=Float64[]
    println("op_ge",length(op_ge))
    tsize=sum([[length(ops_principal)];[length(ops) for i in 1:length(op_ge)];[2*length(ops) for i in 1:length(op_eq)];[1 for i in length(tr_ge)]])

    println("Size of the PSD variable: ", tsize, "x", tsize)
    X = @variable(model, [1:tsize, 1:tsize], PSD)
    unique_mons,unique_pos = cyclic ? cyclic_npa_moments_block_dual!(ops_principal,A,B,tsize; extra_zeros=extra_zeros) : npa_moments_block_dual!(ops_principal,A,B,tsize; extra_zeros=extra_zeros)
    offset = length(ops_principal)
    println("Done building PM")


    for i in 1:length(op_ge)
        if cyclic
            cyclic_npa_moments_block_dual!(ops,A,B,tsize; cPoly=op_ge[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
        else
            npa_moments_block_dual!(ops,A,B,tsize; cPoly=op_ge[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
        end
        offset += length(ops)
    end

    println("Done building LMI")

    for i in 1:length(op_eq)
        if cyclic
            cyclic_npa_moments_block_dual!(ops,A,B,tsize; cPoly=op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            cyclic_npa_moments_block_dual!(ops,A,B,tsize; cPoly=-op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops) 

        else
            npa_moments_block_dual!(ops,A,B,tsize; cPoly=op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            npa_moments_block_dual!(ops,A,B,tsize; cPoly=-op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops)
        end
    end

    # Add the constraints for the principal moment matrix

    for i in 1:length(tr_eq)
        Ai=spzeros(tsize * tsize)
        if cyclic
            for (m,c) in Polynomial(tr_eq[i][1])
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                upi,upj = unique_pos[m_i]
                Ai[(upi-1)*tsize + upj] = c
            end

        else
            tr_eq_poly=real_rep(Polynomial(Polynomial(tr_eq[i][1])))
            for (m,c) in tr_eq_poly
                mi = findfirst(x->x==m,unique_mons)
                upi,upj = unique_pos[mi]
                Ai[(upi-1)*tsize + upj] = c
            end
        end
        push!(A, Ai)
        push!(B, tr_eq[i][2])
    end

    for i in 1:length(tr_ge)

        Ai=spzeros(tsize * tsize)
        Ai[offset*tsize + offset+1] = -1.0

        if cyclic
            for (m,c) in Polynomial(tr_ge_p)
                m1,m2 = (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                upi,upj = unique_pos[m_i]
                Ai[(upi-1)*tsize + upj] = c
            end
        else
            tr_ge_poly=real_rep(Polynomial(tr_ge_p))
            for (m,c) in tr_ge_poly
                m_i = findfirst(x->x==m,unique_mons)
                upi,upj = unique_pos[m_i]
                Ai[(upi-1)*tsize + upj] = c
            end
        end
        offset += 1
        push!(A, Ai)
        push!(B, tr_ge[i][2])
    end

    println("Done building trace constraints")

    if normalize
        id_elem=one(first(ops_principal))
        Ai=spzeros(tsize * tsize)
        if cyclic
            id_elem=cyclic_reduce(id_elem)
            id_i=findfirst(x->x==id_elem,unique_mons)
        else
            id_i = findfirst(x->x==id_elem,unique_mons)
        end
        upi,upj= unique_pos[id_i]   
        Ai[(upi-1)*tsize + upj] = 1.0
        push!(A, Ai)
        push!(B, 1.0)
    end
    println("Done building normalization constraint")
    C = spzeros(tsize * tsize)
    if !is_number(obj)
        if cyclic
            for (m,c) in Polynomial(obj)
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i==nothing
                    throw(ArgumentError("level not enough"))
                end
                upi,upj = unique_pos[m_i]
                C[(upi-1)*tsize + upj] = c
            end
        else
            obj_poly=real_rep(Polynomial(Polynomial(obj)))
            for (m,c) in obj_poly
                mi = findfirst(x->x==m,unique_mons)
                if mi==nothing
                    throw(ArgumentError("level not enough"))
                end
                upi,upj = unique_pos[mi]
                C[(upi-1)*tsize + upj] = c
            end
        end
        
    end
    println("Done building objective")
    if rm
        A = vcat([r' for r in A]...)
        return Dict(zip(unique_mons, unique_pos)),C,A,B
    end

    # Impose constraints AX=B
    lenA=length(A)
    Xvec=vec(X)
    for i in 1:lenA
        println(100*i/lenA, "%")
        Ai = A[i]
        bi = B[i]
        cons = dot(Ai, Xvec)
        @constraint(model, cons == bi)
    end

    # Impose the objective
    if !is_number(obj)
        obj_cons = dot(C, Xvec)
        min ? @objective(model, Min, obj_cons) : @objective(model, Max, obj_cons)
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
        return termination_status(model),model,ops_principal,Dict(zip(unique_mons, unique_pos))
    end

    return optimal_value, model,ops_principal,Dict(zip(unique_mons, unique_pos))

end
