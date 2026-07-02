
struct NCMonoid <: AbstractMonoid
    # Mandatory fields
    name::String
    parent_monoid::Base.RefValue{AbstractMonoid}
    is_built::Base.RefValue{Bool}

    # wrapper for AbstractAlgebra and technical fields
    base_ring::AbstractAlgebra.Generic.FreeAssAlgebra{Rational{BigInt}}
    vertices::Vector{Variable}
    abs_alg_vars::Vector{
        AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}},
    }
    var_dict::Dict{
        Variable,
        AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}},
    }
    rev_var_dict::Dict{
        AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}},
        Variable,
    }
    relations::Array{
        AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}},
        1,
    }
    g_basis::Array{AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}}, 1}
    # To manipulate current monoid within parent monoid
    clique_indices::Vector{Int}
    commutes_with::Vector{AbstractMonoid}
end
function NCMonoid(
    name::String,
    variables::Vector{V};
    parent_monoid::Base.RefValue{AbstractMonoid} = Base.RefValue{AbstractMonoid}(),
) where {V <: Variable}
    base_ring, vars=free_associative_algebra(QQ, [i.name for i in variables])
    var_dict=Dict(zip(variables, vars))
    rev_var_dict=Dict(zip(vars, variables))
    ncmonoid=NCMonoid(
        name,
        parent_monoid,
        Base.RefValue(false),
        base_ring,
        variables,
        vars,
        var_dict,
        rev_var_dict,
        Array{AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}}, 1}([]),
        Array{AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}}, 1}([]),
        Vector{Int}(),
        Vector{AbstractMonoid}(),
    )
    for (j, i) in enumerate(variables)
        i.parent_monoid[]=ncmonoid
        i.monomial[]=NCWord(ncmonoid, vars[j])
    end
    return ncmonoid
end

struct NCWord <: AbstractMonomial
    monoid::NCMonoid
    word::AbstractAlgebra.Generic.FreeAssAlgElem{Rational{BigInt}}
    monomial::Base.RefValue{AbstractMonomial}
    function NCWord(
        monoid::NCMonoid,
        word::AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}},
    )
        new(monoid, word, Base.RefValue{AbstractMonomial}())
    end
    function NCWord(
        monoid::NCMonoid,
        word::AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{Rational{BigInt}},
        m::Base.RefValue{AbstractMonomial},
    )
        new(monoid, word, m)
    end
end

function add_mult!(a, b, c)
    rel=a*b-c
    rel_=conj(rel)
    add_relations!([rel, rel_])
end

function add_relations!(v::Vector)
    for i in v
        monoid=i.monoid
        if i isa Polynomial
            union!(monoid.relations, [polynomial_to_aaword(i)])
        else
            union!(monoid.relations, [i.word])
        end
    end
end

function exponents(m::NCWord)
    vertices=m.monoid.vertices
    if is_identity(m)
        return Dict{Variable, Int}(zip(vertices, zeros(Int, length(vertices))))
    end
    e_w=exponent_word(m.word, 1)
    return Dict{Variable, Int}(
        zip(vertices, [count(==(i), e_w) for i in 1:length(vertices)]),
    )
end

variables(w::NCWord) = Vector{Variable}(
    union(
        [
            map(
                k->w.monoid.rev_var_dict[w.monoid.abs_alg_vars[k]],
                exponent_word(w.word, i),
            ) for i in 1:length(w.word)
        ]...,
    ),
)
degree(w::NCWord) = total_degree(w.word)

function Base.:(==)(M::NCWord, N::NCWord)
    (M.monoid!=N.monoid) && return (length(M)==length(N)==1) &&
           (
               AbstractAlgebra.monomial(M.word, 1)==1 &&
               AbstractAlgebra.monomial(N.word, 1)==1
           ) &&
           (M.word.coeffs[1]==N.word.coeffs[1])
    return (M.word==N.word)
end
Base.:(==)(M::NCWord, N::Number) = M.word==N
Base.:(==)(M::Number, N::NCWord) = M==N.word

Base.:-(v::AbstractMonomial) = -1*v

Base.:+(w1::NCWord, w2::NCWord) = Polynomial(w1)+Polynomial(w2)
Base.:+(x::Number, y::Variable) = x+monomial(y)
Base.:+(x::Variable, y::Number) = y+x
Base.:+(x::Number, m::NCWord) = x*one(m)+m
Base.:+(m::NCWord, x::Number) = x+m
Base.:+(w::NCWord, v::Variable) = w+monomial(v)
Base.:+(v::Variable, w::NCWord) = w+v
Base.:+(v1::Variable, v2::Variable) = monomial(v1)+monomial(v2)

Base.:-(w1::NCWord, w2::NCWord) = w1+(-w2)
Base.:-(x::Number, y::Variable) = x+(-y)
Base.:-(x::Variable, y::Number) = x+(-y)
Base.:-(x::Number, m::NCWord) = x+(-m)
Base.:-(m::NCWord, x::Number) = m+(-x)
Base.:-(w1::NCWord, w2::Variable) = w1+(-w2)
Base.:-(w1::Variable, w2::NCWord) = w1+(-w2)
Base.:-(v1::Variable, v2::Variable) = v1+(-v2)

function Base.:*(w1::NCWord, w2::NCWord)
    return w1.monoid==w2.monoid ? multiply_ncword(w1, w2) : general_mult(w1, w2)
end

# Base.:*(v::Variable,w::NCWord) = v.monomial[]*w
# Base.:*(w::NCWord,v::Variable) = w*v.monomial[]
# Base.:*(w::NCWord,m::AbstractMonomial)=check_ncword_is_poly(w) ? ncword_to_poly(w)*m : general_mult(w,m)
# Base.:*(m::AbstractMonomial,w::NCWord)=check_ncword_is_poly(w) ? ncword_to_poly(w)*m : general_mult(m,w)
Base.:*(x::Number, y::Variable) = x*monomial(y)
Base.:*(x::Variable, y::Number) = y*x
Base.:*(x::Number, m::NCWord) = NCWord(m.monoid, x*m.word)
Base.:*(m::NCWord, x::Number) = x*m

Base.:/(w::NCWord, n::Number) = NCWord(w.monoid, w.word/n)
Base.:/(n::Number, w::NCWord) = divide_ncword(n*one(w.monoid), w)
Base.:/(w1::NCWord, w2::NCWord) =
    (w1.monoid==w2.monoid) ? divide_ncword(w1, w2) : general_divide(w1, w2)
# function Base.:/(w1::NCWord,w2::NCWord) 
#     if w1.monoid != w2.monoid
#         println(w1.monoid == w2.monoid)
#         println("******************************************")
#         println(typeof(w1))
#         println("******************************************")
#         println(typeof(w2))
#     end
#     (w1.monoid==w2.monoid) ? divide_ncword(w1,w2) : general_divide(w1,w2)
# end
divide(w1::NCWord, w2::NCWord) = w1/w2
divide(w::NCWord, n::Number) = w/n
divide(n::Number, w::NCWord) = n/w

Base.:/(w::NCWord, v::Variable) =
    (w.monoid==v.parent_monoid[]) ? w/monomial(v) :
    throw(ArgumentError("Cannot divide words from different monoids"))
Base.:/(v::Variable, w::NCWord) =
    (w.monoid==v.parent_monoid[]) ? monomial(v)/w :
    throw(ArgumentError("Cannot divide words from different monoids"))

Base.hash(m::NCWord, h::UInt) = hash(m.monoid, hash(m.word, hash(0x7d6979235cb005d0, h)))
Base.one(nc::NCMonoid) = NCWord(nc, AbstractAlgebra.one(nc.base_ring))
Base.show(io::IO, mime::MIME"text/plain", m::NCWord) = show(io, m)
Base.show(io::IO, mime::MIME"text/print", m::NCWord) = show(io, m)
Base.show(io::IO, m::NCWord) = print(io, m.word)

function Base.show(io::IO, mime::MIME"text/plain", monoid::NCMonoid)
    print(
        io,
        "The Non-Commutative MONOID IS $(monoid.is_built[] ? "BUILT" : "NOT YET BUILT")\n 
*** Variables *** \n# of variables = $(length(monoid.vertices))\nvariables = $(vars_to_str(monoid.vertices))\n
*** Relations are ***",
    )
    for rel in monoid.relations
        print(io, "\n$rel")
    end
end

function Base.isless(w1::NCWord, w2::NCWord)
    (w1.monoid!=w2.monoid) &&
        w1.monoid.parent_monoid[]!=w2.monoid.parent_monoid[] &&
        throw(ArgumentError("NOT YET IMPLEMENTED"))
    if (w1.monoid==w2.monoid)
        w1.word<w2.word
    else
        vertices=w1.monoid.parent_monoid[].vertices
        return findfirst(item -> item == w1.monoid, vertices) <
               findfirst(item -> item == w2.monoid, vertices)
    end
end
Base.:<(w1::NCWord, w2::NCWord) = isless(w1, w2)

Base.length(m::NCWord) = length(m.word)
Base.iterate(m::NCWord) = iterate(ncword_to_poly(m))
Base.iterate(m::NCWord, state) = iterate(ncword_to_poly(m), state)
Base.getindex(m::NCWord, i::Int64) = collect(monomials(ncword_to_poly(m)))[i]

function Base.conj(w::NCWord)
    w_=w.word
    res_word=zero(w_)
    vars=w.monoid.vertices
    var_dict=w.monoid.var_dict
    for (index, coeff) in zip(1:length(w_), collect(coefficients(w_)))
        rev_expo_word=reverse(exponent_word(w_, index))
        res_word+=coeff'*prod(
            [var_dict[conj(vars[value])] for value in rev_expo_word],
            init = 1,
        )
    end
    return NCWord(w.monoid, res_word)
end
Base.conj(v::Vector{W}) where {W <: AbstractMonomial} = conj.(v)
Base.adjoint(w::NCWord) = conj(w)

function multiply_ncword(w1::NCWord, w2::NCWord)
    monoid=w1.monoid
    replacement_rules=monoid.g_basis
    res=(monoid.is_built[] && !isempty(replacement_rules)) ?
        normal_form(w1.word*w2.word, replacement_rules) : w1.word*w2.word
    return transform_ncword(NCWord(monoid, res))
end
function divide_ncword(w1::NCWord, w2::NCWord)
    ok, l, r = Generic.word_divides_leftmost(w1.word.exps[1], w2.word.exps[1])
    if !ok
        return false, (nothing, nothing)
    else
        isempty(l) ? (l=one(w1.monoid)) :
        l=one(w1.monoid)*prod([w1.monoid.vertices[i] for i in l])
        isempty(r) ? (r=one(w1.monoid)) :
        r=one(w1.monoid)*prod([w1.monoid.vertices[i] for i in r])
        return true, (l, r)
    end
end
function add_ncword(w1::NCWord, w2::NCWord)
    res=w1.word+w2.word
    return NCWord(w1.monoid, res)
end

function Polynomial(m::NCWord, n)
    return ncword_to_poly(m)*n
end
function Polynomial(m::NCWord)
    return ncword_to_poly(m)
end

Base.:+(p::Polynomial, n::NCWord) = p+ncword_to_poly(n)
Base.:+(n::NCWord, p::Polynomial) = p+n
Base.:-(p::Polynomial, n::NCWord) = p-ncword_to_poly(n)
Base.:-(n::NCWord, p::Polynomial) = ncword_to_poly(n)-p
Base.:*(p::Polynomial, n::NCWord) = p*ncword_to_poly(n)
Base.:*(n::NCWord, p::Polynomial) = ncword_to_poly(n)*p

function Base.isless(w1::NCWord, w2::GraphProductWord)
    w1==w2 && return false
    w1.monoid.parent_monoid[]!=w2.monoid.parent_monoid[] &&
        throw(ArgumentError("NOT YET IMPLEMENTED"))
    vertices==w1.monoid.parent_monoid[].vertices
    return findfirst(item -> item == w1.monoid, vertices) <
           findfirst(item -> item == w2.monoid, vertices)
end
Base.isless(w1::GraphProductWord, w2::NCWord) = !isless(w2, w1)
Base.:<(w1::NCWord, w2::GraphProductWord) = isless(w1, w2)
Base.:<(w1::GraphProductWord, w2::NCWord) = !isless(w2, w1)

function build(m::NCMonoid)
    push!(
        m.g_basis,
        AbstractAlgebra.groebner_basis(
            Array(m.relations),
            2*maximum(total_degree.(m.relations), init = 0),
        )...,
    )
    m.is_built[]=true
    nothing
end
function build(m::NCMonoid, degree::Int)
    push!(m.g_basis, AbstractAlgebra.groebner_basis(m.relations, degree)...)
    m.is_built[]=true
    nothing
end
