using JuMP,Mosek,MosekTools

function mons_at_levelint(list_vars::Vector{Variable},level::Int)
    if isempty(list_vars)
        @error "The list of variables is empty. Please provide a non-empty list of variables."
    end
    if level==0
        return [one(prod(list_vars))]
    end
    Id=one(list_vars[1])
    if level==0
        return [Id]
    end
    mons=AbstractMonomial[Id]
    union!(mons,list_vars)
    if level==1
        return mons
    end
    temp_=copy(list_vars)
    for _ in 2:level
        temp_=union(monomials.(filter(res->(res!=Id && !(typeof(res)<:Number) ),kron(temp_,list_vars,1)))...)
        union!(mons,temp_)
    end
    return mons
end
function mons_at_level(list_vars::Vector{Variable},level::String)
    if isempty(list_vars)
        @error "The list of variables is empty. Please provide a non-empty list of variables."
    end
    if level=="0"
        return [one(prod(list_vars))]
    end
    Id=one(list_vars[1])
    list_vars=unique(list_vars)
    lvl_array=[]
    lvl_str=split(level,"+")
    total_types=Set([])
    type_var_dict=Dict{String,Vector{Variable}}()
    for i in lvl_str
        # if tryparse(Int, i) !== nothing
        #     push!(lvl_array,parse(Int,i))
        # else
        #     i_types=split(i,"*")
        #     push!(total_types,i_types...)
        #     push!(lvl_array,i_types)
        # end
        i_types=split(i,"*")
        li_arr=[]
        for j in i_types
            if tryparse(Int, j) === nothing
                push!(total_types,j)
                push!(li_arr, j)
            else
                push!(li_arr,parse(Int,j))
            end
        end
        push!(lvl_array, li_arr)
    end

    for i in total_types
        na_nu=string.(split(i,"["))
        if(typeof(list_vars[1].parent_monoid[]))<:GraphProductMonoid
            if length(na_nu)==1
                type_var_dict[i]=filter(x->extract_string_before_number(x.name)==i,list_vars)
            else
                rng_arr=parse_range(na_nu[2][1:end-1])
                type_var_dict[i]=filter(x->extract_string_before_number(x.name)==na_nu[1] && extract_index(x.name) in rng_arr,list_vars)
            end
        else
            if length(na_nu)==1
                type_var_dict[i]=filter(x->x.parent_monoid[].name==i,list_vars)
            else
                rng_arr=parse_range(na_nu[2][1:end-1])
                type_var_dict[i]=filter(x->extract_string_before_number(x.name)==na_nu[1] && extract_index(x.name) in rng_arr,list_vars)
            end
        end
    end
    # mons always has Id
    mons=AbstractMonomial[Id]

    for l in lvl_array
        mons_l=[]
        for l_ in l
            if l_ isa Int
                push!(mons_l,mons_at_levelint(list_vars,l_))
            else
                push!(mons_l,type_var_dict[l_])
            end
        end
        if mons_l==[]
            continue
        end
        if length(mons_l)==1
            union!(mons,monomials.(mons_l[1])...)
            continue
        end
        union!(mons,monomials.(kron(mons_l...))...)
    end
    return monomials(sum(mons))
end


function mons_at_level(list_vars::Vector{Variable},level::Int)
    return mons_at_level(list_vars, string(level))
end

function mons_at_level(p::Polynomial,level::String)
    if isempty(variables(p))
        return [one(p.monoid)]
    end
    return mons_at_level(variables(p),level)
end

function mons_at_level(p::Polynomial,level::Int)
    if isempty(variables(p))
        return [one(p.monoid)]
    end
    return mons_at_level(variables(p),level)
end

function mons_at_level(M::AbstractMonoid, k::Int)
    return mons_at_level(variables(M), k)
end

function get_monomials(obj, level; 
    op_eq = [], 
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    list_vars=[])
    # for a given level of the localising moment matrices, it returns the operators at that level and the ones of the principal moment matix (which will be larger)
    if !isempty(list_vars)
        ops= mons_at_level(list_vars, level)
    else
        pols=vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])
        vars=unique_array(union([variables(g) for g in pols]...))
        ops = mons_at_level(vars, level)
    end
    if isempty(op_ge) && isempty(op_eq)
        return ops, ops
    end
    deg = Int(ceil(maximum([[degree(g) for g in op_ge];[degree(g) for g in op_eq]])/2))
    extra_vars= unique_array(union([[variables(g) for g in op_ge];[variables(g) for g in op_eq]]...))
    ops_add = mons_at_level(extra_vars, deg)
    ops_principal = unique_array([ops_add[o]*ops[p] 
            for o in 1:length(ops_add) for p in 1:length(ops)])
    return ops, ops_principal
end


function cyclic_npa_moments_block!(list_monomials::Vector{M},X,tsize,model;cPoly=1,unique_mons=[],unique_pos=[],offset=0,extra_zeros=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials
                
            if cPoly == 1
                m = list_monomials[i]'* list_monomials[j]
                m1,m2 = cyclic_reduce(m),cyclic_reduce(m')
                if m1==0 || m2==0
                    @constraint(model, X[i,j] == 0)
                    continue
                end
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    upi,upj = unique_pos[m_i]
                    @constraint(model, X[i,j] - X[upi,upj] == 0)

                else
                    push!(unique_mons, m1)
                    push!(unique_pos, (i,j))
                end
            else
                monomial_product = list_monomials[i]'* cPoly * list_monomials[j]
                if monomial_product == 0
                    @constraint(model, X[offset+i,offset+j] == 0)
                    continue
                end
                PMpoly = 0
                for (m,c) in monomial_product
                    m1,m2 = cyclic_reduce(m),cyclic_reduce(m')
                    mi = findfirst(x->(x==m1 || x==m2),unique_mons)
                    if mi === nothing
                        throw(ArgumentError("level not enough"))
                    end
                    upi,upj = unique_pos[mi]
                    PMpoly += c*X[upi,upj]
                end
                @constraint(model, X[offset+i,offset+j] - PMpoly == 0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                @constraint(model, X[offset+i , j] == 0.0)
            end
        end

        for i in 1:num_monomials
            for j in offset+num_monomials+1:tsize
                @constraint(model, X[offset+i , j] == 0.0)
            end
        end
    end

    if cPoly == 1
        return unique_mons,unique_pos
    end

end

function npa_moments_block!(list_monomials::Vector{M},X,tsize,model;cPoly=1,unique_mons=[],unique_pos=[],offset=0,extra_zeros=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials

            if cPoly == 1
                m = real_rep(list_monomials[i]'* list_monomials[j])
                if m == 0
                    @constraint(model, X[i,j] == 0)
                    continue
                end
                m_i=findfirst(x->(x==m),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    upi,upj = unique_pos[m_i]
                    @constraint(model, X[i,j] - X[upi,upj] == 0)
                else
                    push!(unique_mons, m)
                    push!(unique_pos, (i,j))
                end
            else
                monomial_product = real_rep(Polynomial(list_monomials[i]'* cPoly * list_monomials[j]))
                # assume the monomials PM exists and m and m' are in pm.
                if monomial_product == 0
                    @constraint(model, X[offset+i,offset+j] == 0)
                    continue
                end
                PMpoly = 0
                for (m,c) in monomial_product
                    m_i=findfirst(x->(x==m),unique_mons)
                    upi,upj = unique_pos[m_i]
                    PMpoly += c*X[upi,upj]
                end
                @constraint(model, X[offset+i,offset+j] - PMpoly == 0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                Ai=spzeros(tsize*tsize)
                @constraint(model,X[offset+i , j] == 0.0)
            end
        end

        for i in 1:num_monomials
            for j in offset+num_monomials+1:tsize
                @constraint(model,X[offset+i , j] == 0.0)
            end
        end
    end

    if cPoly == 1
        return unique_mons,unique_pos
    end
end
function npa(obj, level;
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
    extra_zeros=false
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
    tsize=sum([[length(ops_principal)];[length(ops) for i in 1:length(op_ge)];[2*length(ops) for i in 1:length(op_eq)]])

    println("Size of the PSD variable: ", tsize, "x", tsize)
    X = @variable(model, [1:tsize, 1:tsize], PSD)
    unique_mons,unique_pos = cyclic ? cyclic_npa_moments_block!(ops_principal,X,tsize,model; extra_zeros=extra_zeros) : npa_moments_block!(ops_principal,X,tsize,model; extra_zeros=extra_zeros)
    offset = length(ops_principal)
    println("Done building PM")

    for i in 1:length(op_ge)
        if cyclic
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=op_ge[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
        else
            npa_moments_block!(ops,X,tsize,model; cPoly=op_ge[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
        end
        offset += length(ops)
    end

    println("Done building LMI")

    for i in 1:length(op_eq)
        if cyclic
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=-op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops) 

        else
            npa_moments_block!(ops,X,tsize,model; cPoly=op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            npa_moments_block!(ops,X,tsize,model; cPoly=-op_eq[i], unique_mons=unique_mons, unique_pos=unique_pos, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops)
        end
    end

    # Add the constraints for the principal moment matrix

    for i in 1:length(tr_eq)
        tr_eq_p = 0
        if cyclic
            for (m,c) in Polynomial(tr_eq[i][1])
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                upi,upj = unique_pos[m_i]
                tr_eq_p += c*X[upi,upj]
            end

        else
            tr_eq_poly=real_rep(Polynomial(Polynomial(tr_eq[i][1])))
            for (m,c) in tr_eq_poly
                mi = findfirst(x->x==m,unique_mons)
                upi,upj = unique_pos[mi]
                tr_eq_p += c*X[upi,upj]
            end
        end
        @constraint(model, tr_eq_p - tr_eq[i][2] == 0)
    end

    for i in 1:length(tr_ge)
        tr_ge_p=0
        if cyclic
            for (m,c) in Polynomial(tr_ge[i][1])
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                tr_ge_p+=c*X[unique_pos[m_i]...]
            end
        else
            tr_ge_poly=real_rep(Polynomial(Polynomial(tr_ge[i][1])))
            for (m,c) in tr_ge_poly
                m_i = findfirst(x->x==m,unique_mons)
                tr_ge_p+=c*X[unique_pos[m_i]...]
            end
        end
        @constraint(model, tr_ge_p >= tr_ge[i][2])
    end

    println("Done building trace constraints")

    if normalize
        id_elem=one(first(ops_principal))
        if cyclic
            id_elem=cyclic_reduce(id_elem)
            id_i=findfirst(x->x==id_elem,unique_mons)
        else
            id_i = findfirst(x->x==id_elem,unique_mons)
        end
        upi,upj= unique_pos[id_i]   
        @constraint(model, X[upi,upj] == 1)
    end
    println("Done building normalization constraint")
    obj_p = 0
    if !is_number(obj)
        if cyclic
            for (m,c) in Polynomial(obj)
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i==nothing
                    throw(ArgumentError("level not enough"))
                end
                upi,upj = unique_pos[m_i]
                obj_p += c*X[upi,upj]
            end
        else
            obj_poly=real_rep(Polynomial(Polynomial(obj)))
            for (m,c) in obj_poly
                mi = findfirst(x->x==m,unique_mons)
                if mi==nothing
                    throw(ArgumentError("level not enough"))
                end
                upi,upj = unique_pos[mi]
                obj_p += c*X[upi,upj]
            end
        end
        
    end
    println("Done building objective")
    if rm
        return Dict(zip(unique_mons, unique_pos)),model,ops_principal
    end


    # Impose the objective
    if !is_number(obj)
        min ? @objective(model, Min, obj_p) : @objective(model, Max, obj_p)
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
