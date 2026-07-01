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


function basis_gen(obj,level,op_eq,op_ge,tr_eq,tr_ge,list_vars,lvl_lm)
        # Delete redundant polynomials(polynomials which are numbers, so k*Id) from the op_ge and op_eq and warn for incompatible polynomials in op_ge and op_eq
    !isempty(op_ge) && (op_ge=sanity_check_op_ge(op_ge))
    !isempty(op_eq) && (op_eq=sanity_check_op_eq(op_eq))
    if all(is_number.(vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])))
        @warn "All the input polynomials are constants. It is a feasibility problem."
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
    return ops, ops_principal
end