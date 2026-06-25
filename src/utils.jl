# num_to_var(n::Int64,M::Base.RefValue{Monoid})=M[].variables[n]
# num_to_var(x::Array{Int64},M::Base.RefValue{Monoid})=num_to_var.(x,M)
# var2num(x::Variable)=x.value
# var2num(x::Array{Variable})=var2num.(x)
reverse_dict(dict)=Dict(value => key for (key, value) in dict)
Base.sort(dict::Dict)=sort(collect(dict), by = x->x[1])
IndexMap_rev= reverse_dict(IndexMap)

"""
    check_exponents_consistency(monomial::Monomial{VariableType}) where VariableType

    Checks the consistency of exponents in a monomial.

    # Arguments
    - `monomial`: A monomial in a monoid.

    # Returns
    - `false` if the exponents of any variable in the monomial are not consistent.
    - The exponents of the monomial if all exponents are consistent.

    # Notes
    The function first extracts the variables from the monoid of the monomial. Then, for each variable, it calculates the exponents of the variable in each of its cliques. If the exponents are not all the same, it returns `false`. If it has checked all variables and all their exponents are consistent, it returns the exponents of the monomial.
"""
function check_exponents_consistency(w::GraphProductWord{T}) where {T<:AbstractMonomial}
    
    vertices=w.monoid.vertices
    exponents = T<:Variable ? Dict{Variable,Vector{Variable}}([]) : Dict{AbstractMonoid,Vector{AbstractMonomial}}([])
    for vertex in vertices
        exponent=[clique_projector(w.clique_words[i],[vertex]) for i in vertex.clique_indices]
        if length(unique(exponent))!=1
            return false
        end
        exponents[vertex]=exponent[1]
    end
    return exponents
end

function get_clique_indices(monoid::m,cliques::Vector{<:Vector}) where {m<:Union{AbstractMonoid,Variable}}
    return [clique_index for (clique_index,clique) in enumerate(cliques) if monoid in clique]
end

function merge(x::Vector{Vector{w}},y::Vector{Vector{w}}) where w<: Any
    n = length(x)
    result = Vector{w}[]
    @inbounds for i in 1:n
        push!(result, vcat(x[i], y[i]))
    end
    return result
end

function merge(x::Dict{w,Int},y::Dict{w,Int}) where w<:AbstractMonomial

    result=Dict{w,Int}(x)
    for (key,value) in y
        result[key]=get(result,key,0)+value
    end
    return result
end

# Iterative merge for Dict{w,Int}
function merge(x::Dict{w,Int}, y::Dict{w,Int}, z::Dict{w,Int}...) where w<:AbstractMonomial
    result = merge(x, y)
    for arg in z
        result = merge(result, arg)
    end
    return result
end

merge(x::Dict{w,Int}) where w<:AbstractMonomial = x

"""
    subsequence(A::Vector{V}, B::Vector{V}) where V

    Check if vector `A` is a subsequence of vector `B`. 

    # Arguments
    - `A::Vector{V}`: The vector to check if it is a subsequence.
    - `B::Vector{V}`: The vector in which to check for the subsequence.

    # Returns
    - If `A` is a subsequence of `B`, returns a tuple of vectors representing the left and right factors of `B` around `A`. 
    - If `A` is not a subsequence of `B`, returns `false`.

    # Example
    ```julia
    julia> subsequence([2, 3], [1, 2, 3, 4])
    ([(1,), (4,)])

    # Notes
    This function uses a sliding window approach to check for subsequences.
"""
Base.copy(w::NCWord)=NCWord(w.monoid,copy(w.word),w.monomial)
function subsequence(A::Vector,B::Vector)

    indices_to_check_in_B=1:length(B)-length(A)+1
    # for subsequence check, A==B[index:index+length(A)-1] for some index
    left_and_right_factors=[(B[1:index-1],B[index+length(A):end]) for index in indices_to_check_in_B if B[index:index+length(A)-1] == A]
    return length(left_and_right_factors)>0 ? left_and_right_factors : false
end

function check_base(x,y)

    x==0||y==0 && return nothing
    x.monoid!=y.monoid && throw("Different monoids")
end

function Base.copy(v::Vector{Vector{V}}) where V<:AbstractMonomial

    z=Vector{Vector{V}}(undef,length(v))
    for (i,j) in enumerate(v)
        z[i]=copy(j)
    end
    return z
end

Base.:!(x::Vector)=length(x)==0 ? true : false

vars_to_str(x::Vector)=join([i.name for i in x],",")


"""
    index_elements(clique_words::Vector{Vector{V}}, position::Symbol, clique_indices::Vector{Int}) where V

    Indexes elements in the clique words based on their position.

    # Arguments
    - `clique_words`: A vector of vectors, it is the clique_words of soem monomial.
    - `position`: A symbol representing the position to index, can only be :first or :last.
    - `clique_indices`: A vector of indices of cliques to be checked.

    # Returns
    - A vector of variables at the `position` in the `clique_words` of the monomial
"""
function index_elements(clique_words::Vector{Vector{m}}, position::Symbol, clique_indices::Vector{Int}) where m<:AbstractMonomial
    return Vector{m}([eval(:($position($clique_words[$clique_index]))) for clique_index in clique_indices if !isempty(clique_words[clique_index])])
end

"""
    is_edge(variable::V, clique_words::Vector{Vector{V}}, position::Symbol) where V

    Checks if a variable forms an edge in the clique words.

    # Arguments
    - `variable`: A variable to check.
    - `clique_words`: A vector of vectors, each representing a clique word.
    - `position`: A symbol representing the position to check, can be :first or :last.

    # Returns
    - `true` if the variable forms an edge, `false` otherwise.
"""
function is_edge(var::m, clique_words::Vector{Vector{m}}, position::Symbol) where m <: AbstractMonomial
    clique_indices=isa(var,Variable) ? var.clique_indices : var.monoid.clique_indices

    for clique_index in clique_indices
        if isempty(clique_words[clique_index]) || (eval(:($position($clique_words[$clique_index]))) != var)
            return false
        end
    end
    return true
end

"""
get_edge_variables(clique_words::Vector{Vector{V}}, position::Symbol, clique_indices::Vector{Int}) where V

Gets the edges in the clique words.

# Arguments
- `clique_words`: A vector of vectors, each representing a clique word.
- `position`: A symbol representing the position to check, can be :first or :last.
- `clique_indices`: A vector of indices of cliques.

# Returns
- A set of edges.
"""
function get_edge_variables(clique_words::Vector{Vector{m}}, position::Symbol, clique_indices::Vector{Int}) where m <: AbstractMonomial

    potential_edge_variables=Set(index_elements(clique_words, position, clique_indices))
    return filter(variable -> is_edge(variable, clique_words, position), potential_edge_variables)
end

function get_edge_variables(clique_words::Vector{Vector{m}}, position::Symbol) where m<: AbstractMonomial

    return get_edge_variables(clique_words, position, collect(1:length(clique_words)))
end
"""
    clique_projector(elements::Vector{V}, clique::Vector{V}) where V

Filters elements based on their presence in a given clique.

# Arguments
- `elements`: A vector of elements.
- `clique`: A vector representing a clique.

# Returns
- A vector containing only the elements of `elements` that are also in `clique`.

# Notes
The function uses the `filter` function to keep only the elements in `elements` that are also present in `clique`.
"""
clique_projector(word::Vector{W},monoids::Vector{M}) where {W<: AbstractMonomial,M<: AbstractMonoid} =filter(i->i.monoid in monoids, word)
clique_projector(word::Vector{Variable},monoids::Vector{Variable}) =filter(i->i in monoids, word)

function get_root_monomials(m1::AbstractMonomial,m2::AbstractMonomial)
    if isa(m1,Variable)
        m1=m1.monomial[]
    end
    if isa(m2,Variable)
        m2=m2.monomial[]
    end

    n1,n2=Base.RefValue(m1),Base.RefValue(m2)

    while true
        while true
            if n1[].monoid==n2[].monoid
                return n1[],n2[]
            end
            !isdefined(n2[].monoid.parent_monoid,:x) && break
            # isdefined(n2[].monomial,:x) ? (n2=n2[].monomial) : (n2[].monomial[]=words_to_monomial(n2[].monoid.parent_monoid[],[n2[]]);n2=n2[].monomial)
            n2=Base.RefValue(monomial(n2[]))
        end
        n2=Base.RefValue(m2)
        !isdefined(n1[].monoid.parent_monoid,:x) && break
        # isdefined(n1[].monomial,:x) ? (n1=n1[].monomial) : (n1[].monomial[]=words_to_monomial(n1[].monoid.parent_monoid[],[n1[]]);n1=n1[].monomial)
        n1=Base.RefValue(monomial(n1[]))
    end
    throw("No common root_monoid")
end
        
function raise_monomial(m::AbstractMonomial,M::AbstractMonoid)
    if isa(m,Variable)
        m=monomial(m)
        m.monoid==M && return m
    end
    while m.monoid!=M
        m=monomial(m)
        m.monoid==M && return m
    end
end


function get_root_polynomials(p::Polynomial,q::Polynomial)
    if p.monoid==q.monoid
        return (p,q)
    end
    f_m_p,f_m_q=(one(p.monoid),one(q.monoid))
    f_m_p_,f_m_q_=get_root_monomials(f_m_p,f_m_q)

    if f_m_p.monoid==f_m_p_.monoid && f_m_q.monoid!=f_m_q_.monoid
        Mon=p.monoid
        q_=Polynomial{q.coeff_type}(Mon)
        for i in 1:length(q.monomials)
            push!(q_.monomials,raise_monomial(q.monomials[i],Mon))
            push!(q_.coeffs,q.coeffs[i])
        end
        return (p,q_)
    elseif f_m_p.monoid!=f_m_p_.monoid && f_m_q.monoid==f_m_q_.monoid
        Mon=q.monoid
        p_=Polynomial{p.coeff_type}(Mon)
        for i in 1:length(p.monomials)
            push!(p_.monomials,raise_monomial(p.monomials[i],Mon))
            push!(p_.coeffs,p.coeffs[i])
        end
        return (p_,q)
    else
        Mon=f_m_p_.monoid
        p_=Polynomial(Mon)
        q_=Polynomial(Mon)
        for i in 1:length(p.monomials)
            push!(p_.monomials,raise_monomial(p.monomials[i],Mon))
            push!(p_.coeffs,p.coeffs[i])
        end
        for i in 1:length(q.monomials)
            push!(q_.monomials,raise_monomial(q.monomials[i],Mon))
            push!(q_.coeffs,q.coeffs[i])
        end
        return (p_,q_)
    end
end

get_root_polynomials(p::Polynomial,m::M) where M<:AbstractMonomial=get_root_polynomials(p,Polynomial(m))
get_root_polynomials(m::M,p::Polynomial) where M<:AbstractMonomial=get_root_polynomials(Polynomial(m),p)

function get_root_monomials(m::NCWord,n::AbstractMonomial)
    if n isa Variable && n.parent_monoid == m.monoid
        return (m,monomial(n))
    end
    check_ncword_is_poly(m) ? get_root_polynomials(ncword_to_poly(m),n) : Base.invoke(get_root_monomials,Tuple{AbstractMonomial,AbstractMonomial},m,n)
end

get_root_monomials(n::N,m::NCWord) where N<:AbstractMonomial = reverse(get_root_monomials(m,n))

function get_root_monomials(m::NCWord,n::NCWord)
    if m.monoid==n.monoid
        return (m,n)
    end
    check_ncword_is_poly(m) || check_ncword_is_poly(n) ? get_root_polynomials(ncword_to_poly(m),ncword_to_poly(n)) : Base.invoke(get_root_monomials,Tuple{AbstractMonomial,AbstractMonomial},m,n)
end

function is_identity(x)
    return x == Base.one(x)
end


function ncword_to_poly(m::NCWord)
    if m==0
        return Polynomial{Float64}(m.monoid)
    end
    p=Polynomial{Float64}(m.monoid)
    for i in 1:length(m.word)
        push!(p.monomials,NCWord(m.monoid,AbstractAlgebra.monomial(m.word,i)))
        push!(p.coeffs,m.word.coeffs[i])
    end
    return p
end

function check_ncword_is_poly(m::NCWord)
    if m==0
        return true
    end
    if (length(m.word.coeffs)==1) && (m.word.coeffs[1]==1)
        return false
    end
    return true
 end

function transform_ncword(w::NCWord)
    res=check_ncword_is_poly(w) ? ncword_to_poly(w) :  w
end

function sort_polynomial(p::Polynomial)
    # Zip the monomials and coefficients
    zipped = zip(p.monomials, p.coeffs)
    
    # Sort the zipped pairs based on the monomials
    sorted_zipped = sort(zipped, by = x -> x[1])
    
    # Unzip the sorted result into two vectors
    sorted_monomials, sorted_coeffs = unzip(sorted_zipped)
    
    # Create a new polynomial with the sorted monomials and coefficients
    q = Polynomial(p.monoid)
    return Polynomial(sorted_monomials, sorted_coeffs, p.monoid)
end


function comms_system(G::GraphProductMonoid,M::Vector{NCMonoid})
    for (i,j) in collect(combinations(M,2))
        if isempty(intersect(i.name,j.name))
            push!(G.commutations,(i,j))
        end
    end
end


function extract_string_before_number(s::String)
    for i in 1:length(s)

        if haskey(IndexMap_rev, s[i]) || isdigit(s[i])
            return s[1:i-1]
        end
    end
    return s  # Return the whole string if no number is found
end

function extract_index(s::String)
    for i in 1:length(s)

        if haskey(IndexMap_rev, s[i]) || isdigit(s[i])
            ind_str=s[i:end]
            res=tryparse(Int,ind_str)
            if res !== nothing
                return res
            else
                res= map(c -> IndexMap_rev[c],ind_str)
                return parse(Int,res)
            end
        end
    end
    throw("No index found in string")
end

function parse_range(s::String)
    res= split(s,":")
    ind_start = parse(Int,res[1])
    ind_end = parse(Int,res[2])
    return collect(ind_start : ind_end)
end
# function prmt_terms_coeff(t1::Term,t2::Term,S)
#     a,b=S(t1.coefficient)*t1.monomial,S(t2.coefficient)*t2.monomial
#     return a,b
# end
# Inbuilt unique() does not work with new types
function unique_array(G::AbstractArray{T}) where T
    result = T[]
    for x in G
        if all(y -> !(x == y), result)
            push!(result, x)
        end
    end
    return result
end

function polynomial_to_aaword(p::Polynomial)
    if p==0
        return zero(p.monoid.base_ring)
    end
    return sum([p.coeffs[i]*p.monomials[i].word for i in 1:length(p.monomials)])
end

conj_min(x::Number) = real(x)

conj_min(m::AbstractMonomial) = min(m, conj(m))

function conj_minmax(O::AbstractMonomial)
    Oc = conj(O)

    if O == Oc
        return (O, 0, O)
    elseif O < Oc
        return (O, 1, Oc)
    else
        return (Oc, -1, O)
    end
end

function sanity_check_op_ge(op_ge::Vector)
    for (i,p) in enumerate(op_ge)
        if is_number(p) && p isa Polynomial
            deleteat!(op_ge,i)

            if p isa Polynomial && coefficient(p,one(p.monoid)) < 0
                throw(ArgumentError("Negative constant term in op_ge is not allowed"))
            end
        end
    end
        if isempty(op_ge)
        return 0
    else
        return op_ge
    end
end

function sanity_check_op_eq(op_eq::Vector)
    for (i,p) in enumerate(op_eq)
        if is_number(p)
            deleteat!(op_eq,i)
            if p isa Polynomial && coefficient(p,one(p.monoid)) != 0
                throw(ArgumentError("non-zero constant in op_eq is not allowed"))
            end
        end
    end
    if isempty(op_eq)
        return 0
    else
        return op_eq
    end
end

function close_graph!(g,D,m,E)
    n = nv(g)
    h = SimpleDiGraph(n * 2)

    for e in edges(g)
        add_edge!(h, src(e), dst(e))
    end

    for e in edges(g)
        add_edge!(h, src(e) + n, dst(e) + n)
    end

    pairs = [(w[end],w[1]) for w in m if length(w) > 1]
    pairs = [(D[(p[1],E[p[1]])] , D[(p[2],1)]+n) for p in pairs]

    for pair in pairs
        add_edge!(h, pair...)
    end

    for pair in pairs
        rem_edge!(h, pair...)
        if !has_path(h, pair...)
            add_edge!(g, pair[1], pair[2]-n)
            add_edge!(h, pair...)
        end
    end
end


function trace_mons_reduce(mons::Set,tr_pols::Vector)
    tr_polsc = Vector{Polynomial}([])
    for i in 1:length(tr_pols)
        tr_polsi = Polynomial(tr_pols[i])
        if iszero(tr_polsi)
            throw(ArgumentError("$i -th polynomial is a constant polynomial, which is not allowed. It will be ignored."))
        end
        tr_polsci=0
        for (m,c) in tr_polsi
            m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
            if m1 in mons
                tr_polsci+=c*m1
            elseif m2 in mons
                tr_polsci+=c*m2
            else
                throw(ArgumentError("Monomial $m in tr_pols constraint is not present in any of the LMIs."))
            end
        end
        push!(tr_polsc,tr_polsci)
    end
    return tr_polsc
end

function trace_mons_reduce(mons::Set,tr_pols::Vector{<:AbstractVector})
    just_pols = [i[1] for i in tr_pols]
    return [[pol,tr_pols[i][2]] for (i,pol) in enumerate(trace_mons_reduce(mons,just_pols))]
end

function trace_mons_reduce(mons::Set,tr_pols::Vector{<:Tuple})
    just_pols = [i[1] for i in tr_pols]
    return [(pol,tr_pols[i][2]) for (i,pol) in enumerate(trace_mons_reduce(mons,just_pols))]
end