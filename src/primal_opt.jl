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
    op_eq = 0, 
    op_ge = 0,
    tr_eq = 0,
    tr_ge = 0,
    list_vars=[])
    # for a given level of the localising moment matrices, it returns the operators at that level and the ones of the principal moment matix (which will be larger)
    if !isempty(list_vars)
        ops= mons_at_level(list_vars, level)
    else
        pols=vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])
        vars=unique_array(union([variables(g) for g in pols]...))
        ops = mons_at_level(vars, level)
    end
    if op_ge==0
        return ops, ops
    end
    deg = Int(ceil(maximum([degree(g) for g in op_ge])/2))
    extra_vars= unique_array(union([variables(g) for g in op_ge]...))
    ops_add = mons_at_level(extra_vars, deg)
    ops_principal = unique_array([ops_add[o]*ops[p] 
            for o in 1:length(ops_add) for p in 1:length(ops)])
    return ops, ops_principal
end


function cyclic_npa_moments_block(list_monomials::Vector{M},model;cPoly=1,unique_mons=[],unique_vars=[],eq=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables
    moments_matrix = Matrix{JuMP.AffExpr}(undef, num_monomials, num_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in 1:num_monomials
            
            # Compute the product of the monomials
            monomial_product = list_monomials[i]'* cPoly * list_monomials[j]
            # println("1",typeof(monomial_product))
            # Initialize the JuMP variable for the matrix entry
            moments_matrix[i, j] = 0.0
            # Check if the product already exists in the dictionary
            for (m,c) in monomial_product
                # println("2",typeof(m))
                m1,m2 =  (cyclic_reduce(m),cyclic_reduce(m'))
                if m1==0 || m2==0
                    continue
                end
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    moments_matrix[i, j] += c*unique_vars[m_i]
                else
                    # Create a new JuMP variable
                    new_var = @variable(model)
                    # Store the new variable in the dictionary
                    push!(unique_mons, m1)
                    push!(unique_vars, new_var)
                    # Use the new variable in the matrix
                    moments_matrix[i, j] += c*new_var
                end
            end
        end
    end

    if !eq
        @constraint(model, moments_matrix >=0,PSDCone())
    else
        @constraint(model, moments_matrix .== 0)
    end

    return moments_matrix, unique_mons, unique_vars
end

function npa_moments_block(list_monomials::Vector{M},model;cPoly=1,unique_mons=[],unique_vars=[],eq=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables
    moments_matrix = Matrix{JuMP.AffExpr}(undef, num_monomials, num_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in 1:num_monomials
            # Compute the product of the monomials
            monomial_product = real_rep(Polynomial(list_monomials[i]'* cPoly * list_monomials[j]))
            # Initialize the JuMP variable for the matrix entry
            moments_matrix[i, j] = 0.0
            for (m,c) in monomial_product
                # Check if the product already exists in the dictionary
                m_i=findfirst(x->(x==m),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    moments_matrix[i, j] += c*unique_vars[m_i]
                else
                    # Create a new JuMP variable
                    new_var = @variable(model)
                    # Store the new variable in the dictionary
                    push!(unique_mons, m)
                    push!(unique_vars, new_var)
                    # Use the new variable in the matrix
                    moments_matrix[i, j] += c*new_var
                end
            end
        end
    end
    if !eq
        @constraint(model, moments_matrix in PSDCone())
    else
        @constraint(model, moments_matrix .== 0)
    end
    return moments_matrix, unique_mons, unique_vars
end
function npa(obj, level;
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
    principal_moments_matrix, unique_mons, unique_vars = cyclic ? cyclic_npa_moments_block(ops_principal,model) : npa_moments_block(ops_principal,model)
    # Add the constraints for the principal moment matrix

    if tr_eq!=0
        for i in 1:length(tr_eq)
            tr_eq_p=0
            if cyclic
                for (m,c) in Polynomial(tr_eq[i][1])
                    m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                    m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                    tr_eq_p+=c*unique_vars[m_i]
                end
            else
                tr_eq_poly=real_rep(Polynomial(Polynomial(tr_eq[i][1])))
                for (m,c) in tr_eq_poly
                    m_i=findfirst(x->x==m,unique_mons)
                    tr_eq_p+=c*unique_vars[m_i]
                end
            end
            @constraint(model, tr_eq_p == tr_eq[i][2])
        end
    end
    if tr_ge!=0
        for i in 1:length(tr_ge)
            tr_ge_p=0
            if cyclic
                for (m,c) in Polynomial(tr_ge[i][1])
                    m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                    m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                    tr_ge_p+=c*unique_vars[m_i]
                end
            else
                tr_ge_poly=real_rep(Polynomial(Polynomial(tr_ge[i][1])))
                for (m,c) in tr_ge_poly
                    m_i=findfirst(x->x==m,unique_mons)
                    tr_ge_p+=c*unique_vars[m_i]
                end
            end
            @constraint(model, tr_ge_p >= tr_ge[i][2])
        end
    end
    if op_ge!=0
        for i in 1:length(op_ge)
            if cyclic
                cyclic_npa_moments_block(ops,model; cPoly=op_ge[i],unique_mons=unique_mons, unique_vars=unique_vars) 
            else
                npa_moments_block(ops,model; cPoly=op_ge[i], unique_mons=unique_mons, unique_vars=unique_vars) 
            end
        end
    end
    if op_eq!=0
        for i in 1:length(op_eq)
            if cyclic
                cyclic_npa_moments_block(ops,model; cPoly=op_eq[i],unique_mons=unique_mons, unique_vars=unique_vars,eq=true) 

            else
                npa_moments_block(ops,model; cPoly=op_eq[i], unique_mons=unique_mons, unique_vars=unique_vars,eq=true) 

            end
        end
    end
    if normalize
        id_elem=one(first(ops_principal))
        if cyclic
            id_elem=cyclic_reduce(id_elem)
            @constraint(model,unique_vars[findfirst(check->check==id_elem,unique_mons)]==1.0)
        else
            @constraint(model,unique_vars[findfirst(check->check==id_elem,unique_mons)]==1.0)
        end
    end
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
            obj_poly=real_rep(Polynomial(Polynomial(obj)))
            for (m,c) in obj_poly
                m_i=findfirst(x->x==m,unique_mons)
                if m_i==nothing
                    throw(ArgumentError("level not enough"))
                end
                obj_p+=c*unique_vars[m_i]
            end
        end
        
        min ? @objective(model, Min, obj_p) : @objective(model, Max, obj_p)
    end
    if rm
        return model,Dict(zip(unique_mons,unique_vars)),principal_moments_matrix
    end
    # Solve the model
    optimize!(model)
    # Check the optimization status
    if termination_status(model) != MOI.OPTIMAL
        @warn "The optimization problem is $(termination_status(model))."
    end
    # Extract the optimal value of the objective function
    optimal_value = objective_value(model)
    # Extract the optimal values of the variables and store it as a dictionary
    optimal_vars = Dict{AbstractMonomial, Float64}()
    for i in 1:length(unique_mons)
        optimal_vars[unique_mons[i]] = value(unique_vars[i])
    end
    # Return the optimal value and the dictionary of optimal variables
    # return optimal_value, optimal_vars, model

    if is_number(obj)
        @warn "The objective function is a constant, it is a feasibility check"
        return termination_status(model),model,ops_principal,Dict(zip(unique_mons,unique_vars)),principal_moments_matrix
    end

    return optimal_value, model,ops_principal,Dict(zip(unique_mons,unique_vars)),principal_moments_matrix

end
