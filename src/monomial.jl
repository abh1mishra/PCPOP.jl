"""
    GraphProductWord <: AbstractMonomial

    Structure `GraphProductWord`

    #Fields:
    - `monoid` : Parent monoid of the monomial.
    - `clique_words` : list of projections to cliques as vector of variables.
    - `edge_l` : set of letters that can be pushed to the left.
    - `edge_r` : set of letters that can be pushed to the right.

    Edges are updated during multiplication.

    #Example: 
    Take the partially commuting monoid M = < a, b, c : ab = ba >.
    The maximal cliques (non-commuting sets) are {a,c} and {b,c}.
    The monomial m = a*b*c is stored as:
    - m.monoid = M
    - m.clique_words = [[a,c], [b,c]]
    - m.edge_l = {a, b}
    - m.edge_r = {c}
"""
struct GraphProductWord{T<:AbstractMonomial} <: AbstractMonomial
    monoid::GraphProductMonoid
    monomial::Base.RefValue{AbstractMonomial}
    clique_words::Vector{Vector{T}}
    edge_l::Set{T}
    edge_r::Set{T}
end


# Iterations
Base.length(m::GraphProductWord) = length(m.clique_words)
Base.iterate(m::GraphProductWord) = iterate(m.clique_words)
Base.iterate(m::GraphProductWord, state) = iterate(m.clique_words, state) 
Base.getindex(m::GraphProductWord, i::Int64) = m.clique_words[i]

Base.copy(m::GraphProductWord)=GraphProductWord(m.monoid,m.monomial,copy(m.clique_words),copy(edge_l),copy(edge_r))
# Arithmetic
Base.:/(m::GraphProductWord,n::GraphProductWord)=divides(n,m)
Base.:*(m::GraphProductWord,n::GraphProductWord)=multiply(m,n)

# General Arithmetic

Base.:/(m1::AbstractMonomial,m2::AbstractMonomial)=general_divide(m1,m2)
# Base.:/(v::Variable,m::AbstractMonomial)=general_divide(v,m)
# Base.:/(m::AbstractMonomial,v::Variable)=general_divide(m,v)

Base.:*(m1::AbstractMonomial,m2::AbstractMonomial)=general_mult(m1,m2)

Base.:(==)(m1::AbstractMonomial,m2::AbstractMonomial)=general_equals(m1,m2)

# Display
Base.show(io::IO, mime::MIME"text/plain", m::GraphProductWord) = show(io,m)
Base.show(io::IO, mime::MIME"text/print", m::GraphProductWord)=show(io,m)
Base.show(io::IO,m::GraphProductWord)=show_monomial(io,m)

# Hashing and comparison
Base.hash(m::GraphProductWord, h::UInt)=hash(m.monoid, hash(m.clique_words, hash(0x7d6979235cb005d0, h))) 
Base.:(==)(m::GraphProductWord,n::GraphProductWord)=m.monoid == n.monoid ? (m.clique_words==n.clique_words) : general_equals(m,n)
Base.:(==)(m::GraphProductWord,n::Number)=n==1 && one(m)==m
Base.:(==)(n::Number,m::GraphProductWord)=n==1 && one(m)==m
Base.isless(m::GraphProductWord,n::GraphProductWord)=less_or_not(m,n)
Base.:<(m::GraphProductWord,n::GraphProductWord)=less_or_not(m,n)



one(w::W) where W <:AbstractMonomial = one(w.monoid)
is_identity(x::X) where X<:AbstractMonomial=one(x)==x


function monomial(w::W) where W<:AbstractMonomial
    if isdefined(w.monomial,:x)
        return one(w.monomial[])*w.monomial[]
    end
    monoid=nothing
    if isa(w,Variable) && isa(w.parent_monoid[],GraphProductMonoid)
        try
            monoid=w.parent_monoid[]
        catch e
            throw(ArgumentError("Parent monoid not defined"))
        end
    elseif isa(w.monoid.parent_monoid[],GraphProductMonoid)
        try
            monoid=w.monoid.parent_monoid[]
        catch e
            throw(ArgumentError("Parent monoid not defined"))
        end
    else
        throw(ArgumentError("Cannot convert to monomial, parent monoid is not a GraphProductMonoid"))
    end
    w.monomial[]=words_to_monomial(monoid,[w])
    return w.monomial[]*one(w.monomial[])
end
"""
    variables(monomial::GraphProductWord)

    Extracts and sorts the variables from the clique words of a monomial.

    # Arguments
    - `monomial`: A GraphProductWord object.

    # Returns
    - A sorted list of unique variables present in the clique words of the monomial.
"""
variables(m::GraphProductWord)=Vector{Variable}(variables(m.clique_words))

monomials(m::AbstractMonomial)=[m]

function variables(m_words::Vector{Vector{W}}) where W<:AbstractMonomial
    return  union(vcat(variables.(union(m_words...))...))
end
variables(v::Variable)=[v]
# exponents(m::GraphProductWord)=[i[2] for i in Exponents(m)]

"""
    Exponents(m::GraphProductWord{V}) where V

    computes the exponent of variables in monomial `m`.

    #Arguments
    - `m` : monomial in monoid `M`.

    #Returns
    - dictionary of variables and their exponents.

    #Example
    Take the monoid M = < a, b, c : ab = ba >. The Exponents of m = a*b*a*b^2*c^3 is {a:3, b:3, c:3}
"""
function exponents(m::GraphProductWord{Variable})
    vertices=m.monoid.vertices
    expo=[]
    for vertex in vertices
        push!(expo,length(clique_projector(m.clique_words[vertex.clique_indices[1]],[vertex])))
    end
    return Dict{Variable,Int}(zip(vertices,expo))
    # return Dict(sort(Base.reduce(merge,[exponents(word) for word in expo],init=Dict{Variable,Int}())))
end

function exponents(m::GraphProductWord)
    vertices=m.monoid.vertices
    if is_identity(m)
        return merge([exponents(one(mo)) for mo in vertices]...)
    end
    expo=[]
    for vertex in vertices
        push!(expo,exponents.(clique_projector(m.clique_words[vertex.clique_indices[1]],[vertex]))...)
    end
    return merge(expo...)
    # return Dict(sort(Base.reduce(merge,[exponents(word) for word in expo],init=Dict{Variable,Int}())))
end
degree(x::Number) = !iszero(x) ? 0 : -Inf

"""
    degree(m::GraphProductWord{V}) where {V}

    Computes the degree (total number of variables) of monomial `m`.
"""

# degree(m::GraphProductWord)=sum(values(exponents(m)))
degree(m::GraphProductWord)=sum(values(exponents(m));init=0)



function show_monomial(io::IO, m::GraphProductWord)

    for j in m
        print(io,"($(join([k for k in j],",")))")
    end
end

function general_mult(m1,m2)
    x1,x2=get_root_monomials(m1,m2)
    return x1*x2
end

function general_divide(m1,m2)
    x1,x2=get_root_monomials(m1,m2)
    return x1/x2
end



function general_equals(m1,m2)
    m1==1 && m2==1 && return true
    x1,x2=get_root_monomials(m1,m2)
    return x1==x2
end





function multiply(m::GraphProductWord{W},n::GraphProductWord{W}) where W<: AbstractMonomial
    if m.monoid!=n.monoid
        return general_mult(m,n)
    end

    parent_monoid=m.monoid

    m_words = copy(m.clique_words)
    n_words = copy(n.clique_words)
    m_l=copy(m.edge_l)
    m_r=copy(m.edge_r)
    n_l=copy(n.edge_l)
    n_r=copy(n.edge_r)

    prods=product(m_r,n_l)

    middle_list=Vector{Union{W,Polynomial}}()

    for (i,j) in prods
        if i.monoid==j.monoid
            M=i.monoid
            clique_indices=M.clique_indices
            res=i*j
            push!(middle_list,res)
            for k in clique_indices
                pop!(m_words[k])
                popfirst!(n_words[k])
            end
        end
    end

    if(isempty(middle_list))
        res_clique_words=merge(m_words,n_words)
        m_l=get_edge_variables(res_clique_words,:first)
        n_r=get_edge_variables(res_clique_words,:last)
        return GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),
        res_clique_words,m_l,n_r)
    else
        m_l=get_edge_variables(m_words,:first)
        n_r=get_edge_variables(n_words,:last)
        m_r=get_edge_variables(m_words,:last)
        n_l=get_edge_variables(n_words,:first)
        m_=GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),m_words,m_l,m_r)
        n_=GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),n_words,n_l,n_r)
        middle_prod=prod(middle_list)
        return m_*middle_prod*n_
    end
end



function multiply(m::GraphProductWord{Variable},n::GraphProductWord{Variable})
    if m.monoid!=n.monoid
        return general_mult(m,n)
    end

    parent_monoid=m.monoid

    m_words = copy(m.clique_words)
    n_words = copy(n.clique_words)
    m_l=copy(m.edge_l)
    m_r=copy(m.edge_r)
    n_l=copy(n.edge_l)
    n_r=copy(n.edge_r)

    prods=product(m_r,n_l)

    middle_list=Vector{}()

    for (i,j) in prods
        if j in i.ortho_conj
            return Polynomial(parent_monoid)
        end
        if(i==j && i.mult_type[]==:Projector)
            clique_indices=i.clique_indices
            for index in clique_indices
                pop!(m_words[index])
                popfirst!(n_words[index])
            end
            push!(middle_list,i)
        end
        # unitary and unipotent
        if (i==j && i.mult_type[]==:Unipotent) || (i'==j && i.mult_type[]==:Unitary)
            clique_indices=i.clique_indices
            for index in clique_indices
                pop!(m_words[index])
                popfirst!(n_words[index])
            end
            push!(middle_list,one(i))
        end
    end

    if(isempty(middle_list))
        res_clique_words=merge(m_words,n_words)
        m_l=get_edge_variables(res_clique_words,:first)
        n_r=get_edge_variables(res_clique_words,:last)
        return GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),
        res_clique_words,m_l,n_r)
    else
        m_l=get_edge_variables(m_words,:first)
        n_r=get_edge_variables(n_words,:last)
        m_r=get_edge_variables(m_words,:last)
        n_l=get_edge_variables(n_words,:first)
        m_=GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),m_words,m_l,m_r)
        n_=GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),n_words,n_l,n_r)
        middle_prod=prod(middle_list)
        return m_*middle_prod*n_
    end
end


"""
    Base.conj(m::GraphProductWord{V}) where V

    Takes adjoint of monomial `m` reversing the words with involution:
    (a1 ... ak)* = ak* ... a1*

"""
function conjugate(m::GraphProductWord)

    is_identity(m) && return m
    edge_type=typeof(m.edge_l)
    edge_l = edge_type(Set(conj.(m.edge_r)))
    edge_r = edge_type(Set(conj.(m.edge_l)))
    cw_=typeof(m.clique_words)([])
    for words in m.clique_words
        if isempty(words)
            push!(cw_,words)
        else
            push!(cw_,conj.(reverse(words)))
        end
    end

    return GraphProductWord(m.monoid,Base.RefValue{AbstractMonomial}(),cw_,edge_l,edge_r)

end
Base.conj(m::GraphProductWord{W}) where W<:AbstractMonomial=conjugate(m)
Base.adjoint(m::GraphProductWord)=Base.conj(m)


function Base.conj(m::GraphProductWord{Variable})
    is_identity(m) && return m
    m.monoid.conj_type[] && return conjugate(m)
    w=monomial_to_word(m)
    return prod(conj.(reverse(w)))
end



"""
    less_or_not(m::GraphProductWord, n::GraphProductWord)

Determines if a monomial is less than or not equal to another monomial.

# Arguments
- `m`: A GraphProductWord object.
- `n`: Another GraphProductWord object.

# Returns
- A boolean value. Returns true if `m` is less than `n`, false otherwise.

# Notes
The function first checks if `m` and `n` belong to the same monoid. If not, it throws an error. It then checks if `m` is equal to `n`. If they are equal, it returns false. It then checks the degrees of `m` and `n`. If the degree of `m` is less than the degree of `n`, it returns true. If the degree of `m` is greater than the degree of `n`, it returns false. It then iterates over the cliques in the monoid. If the clique word of `m` is equal to the clique word of `n`, it continues to the next iteration. Otherwise, it returns true if the length of the clique word of `m` is less than the length of the clique word of `n`, or if they have the same length and the clique word of `m` is less than the clique word of `n`.
"""
function less_or_not(m::GraphProductWord,n::GraphProductWord)
    (m.monoid!=n.monoid) && m.monoid.parent_monoid[]!=n.monoid.parent_monoid[] && throw(ArgumentError("NOT YET IMPLEMENTED"))
    if(m.monoid!=n.monoid)
        vertices=m.monoid.parent_monoid[].vertices
        return findfirst(item -> item == m.monoid, vertices) < findfirst(item -> item == n.monoid, vertices)
    end
    (m==n) && return false
    degree(m)<degree(n) && return true
    degree(m)>degree(n) && return false

    for i in 1:length(m.monoid.cliques)
        m.clique_words[i]==n.clique_words[i] && continue
        tot_length_m=sum(degree.(m.clique_words[i]),init=0)
        tot_length_n=sum(degree.(n.clique_words[i]),init=0)
        return ((tot_length_m<tot_length_n) || ((tot_length_m==tot_length_n) && (vcat(m.clique_words[i])<vcat(n.clique_words[i])))) 
    end
end


function compare(m, n)

    if m < n
        return -1
    elseif m > n
        return 1
    else
        @assert m == n
        return 0
    end
end


"""
    divide(m::GraphProductWord{V}, n::GraphProductWord{V}; all=false) where V

    Divide monomial `m` by monomial `n`.

    # Arguments
    - `m::GraphProductWord{V}`: The divisor monomial.
    - `n::GraphProductWord{V}`: The dividend monomial.
    - `all::Bool`: If `true`, returns all possible pairs of factors. If `false`, returns the first pair of factors found.

    # Returns
    - If `m` is divisible by `n`, returns a tuple of monomials representing the left and right factors of `m` around `n`. 
      That is, (l,r) such that l*m*r = n.
    - If `m` is not divisible by `n`, returns an empty array.

    # Notes
    - This function checks if `n` is divisible by `m` by checking if each clique word of `m` is a subsequence of the corresponding clique word of `n`.
    - If `m` is divisible by `n`, it generates all potential factors of `n` by `m` and checks if each potential factor is reconstructible.
    - If a potential factor is reconstructible, it updates the edge variables of the factor and adds it to the result.
    - The time complexity is O(n^2), where n is the number of clique words in `n`.
    - The space complexity is O(n), due to the creation of `potential_factors` and `left_and_right_potential_monomial_factors`.
"""
function divide(m::GraphProductWord,n::GraphProductWord;all=false)

    m==n && return (one(m),one(n))

    # clique_potential_factors ith element contains all potential tuples of left_and_right_word_factors for ith clique of n
    clique_potential_factors=[]
    for (index,clique_word) in enumerate(n.clique_words)
        ith_left_and_right_word_factors=subsequence(m.clique_words[index],clique_word)
        ith_left_and_right_word_factors!=false ? push!(clique_potential_factors,ith_left_and_right_word_factors) : return []
    end

    #= potential_factors is a vector of tuples of tuples, 
    where ith element/tuple represents a potential factor containing tuples where
    each tuple j is a left_and_right_factor for jth clique of n
    =#
    potential_factors=collect(Base.Iterators.product(clique_potential_factors...))

    dummy_edge_set=(eltype(m.monoid.vertices)<:AbstractMonoid) ? Set{AbstractMonomial}() : Set{Variable}()
    left_and_right_potential_monomial_factors=[]
    for i in potential_factors
        left_monomial=GraphProductWord(n.monoid,Base.RefValue{AbstractMonomial}(),[j[1] for j in i ],copy(dummy_edge_set),copy(dummy_edge_set))
        right_monomial=GraphProductWord(n.monoid,Base.RefValue{AbstractMonomial}(),[j[2] for j in i ],copy(dummy_edge_set),copy(dummy_edge_set))
        left_and_right_potential_monomial_factors=push!(left_and_right_potential_monomial_factors,(left_monomial,right_monomial))
    end
    # left_and_right_potential_monomial_factors=[(GraphProductWord(n.monoid,n.monomial,[j[1] for j in i ],dummy_edge_set,dummy_edge_set),GraphProductWord(n.monoid,n.monomial,[j[2] for j in i ],dummy_edge_set,dummy_edge_set)) for i in potential_factors]
    result=[]
    for (left_potential_monomial_factor,right_potential_monomial_factor) in left_and_right_potential_monomial_factors
        if is_reconstructible(left_potential_monomial_factor) & is_reconstructible(right_potential_monomial_factor)
            left_factor=left_potential_monomial_factor
            right_factor=right_potential_monomial_factor
            clique_indices=collect(1:length(left_factor.clique_words))
            union!(left_factor.edge_l,get_edge_variables(left_factor.clique_words,:first,clique_indices))
            union!(right_factor.edge_r,get_edge_variables(right_factor.clique_words,:last,clique_indices))
            union!(left_factor.edge_r,get_edge_variables(left_factor.clique_words,:last,clique_indices))
            union!(right_factor.edge_l,get_edge_variables(right_factor.clique_words,:first,clique_indices))
            all ? push!(result,(left_factor,right_factor)) : return (left_factor,right_factor)
        end
    end
    return result
end
divides(a::GraphProductWord,b::GraphProductWord) =!isempty(divide(a,b))
Base.zero(m::GraphProductWord)=0
