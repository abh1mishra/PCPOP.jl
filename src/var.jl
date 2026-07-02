const IndexMap = Dict{Char, Char}(
    '-' => '₋',
    '0' => '₀',
    '1' => '₁',
    '2' => '₂',
    '3' => '₃',
    '4' => '₄',
    '5' => '₅',
    '6' => '₆',
    '7' => '₇',
    '8' => '₈',
    '9' => '₉',
)
"""
    struct Variable <: AbstractMonomial

    Represents a variable in a monoid.

    # Fields
    - `name`: A string representing the name of the variable.
    - `value`: An integer representing the value of the variable.
    - `monoid`: A reference to the Monoid object that the variable belongs to.
    - `monomial`: A reference to the AbstractWord object that the variable belongs to.
    - `conj`: A reference to the conjugate of the variable. Can also be `nothing`.
    - `clique_indices`: A vector of integers representing the indices of the cliques that the variable belongs to.
    - `commutes_with`: A vector of Variable objects representing the variables that this variable commutes with.

    # Constructor
    The constructor takes a `name` and an optional `value` (default is 0). It initializes `monoid`, `monomial`, and `conj` as empty references, and `clique_indices` and `commutes_with` as empty vectors.
"""
struct Variable <: AbstractMonomial
    name::Union{String, Char}
    parent_monoid::Base.RefValue{AbstractMonoid}
    monomial::Base.RefValue{AbstractMonomial}
    conj::Base.RefValue{Union{Variable, Nothing}}
    mult_type::Base.RefValue{Symbol}
    ortho_conj::Vector{Variable}
    # To manipulate current monoid within parent monoid
    clique_indices::Vector{Int}
    commutes_with::Vector{Variable}
    function Variable(name)
        new(
            string(name),
            Base.RefValue{AbstractMonoid}(),
            Base.RefValue{AbstractMonomial}(),
            Base.RefValue{Union{Variable, Nothing}}(nothing),
            Base.RefValue{Symbol}(:Free),
            Vector{Variable}(),
            Vector{Int}(),
            Vector{Variable}(),
        )
    end
end

name(v::Variable) = v.name
# function name_base_indices(v::Variable)
#     splits = split(v.name, r"[\[,\]]\s*", keepempty = false)
#     if length(splits) == 1
#         return v.name, Int[]
#     else
#         return splits[1], parse.(Int, splits[2:end])
#     end
# end

Base.show(io::IO, v::Variable) = print(io, v.name)

Base.:(==)(v::Variable, b::Number) = false
Base.:(==)(a::Variable, b::Variable) = (a.name==b.name)

# Base.isless(a::Variable,b::Variable)=isless(a.name,b.name)
Base.hash(v::Variable, h::UInt) = hash(v.name, h)

word(v::Variable) = v.word[]

function Base.isless(a::Variable, b::Variable)
    (a.parent_monoid[]!=b.parent_monoid[]) &&
        a.parent_monoid[].parent_monoid[]!=b.parent_monoid[].parent_monoid[] &&
        throw(ArgumentError("NOT YET IMPLEMENTED"))
    if (a.parent_monoid[]!=b.parent_monoid[])
        vertices=a.parent_monoid[].parent_monoid[].vertices
        return findfirst(item -> item == a.parent_monoid[], vertices) <
               findfirst(item -> item == b.parent_monoid[], vertices)
    end
    return a.name<b.name
end
Base.:<(a::Variable, b::Variable) = isless(a, b)
# Base.:(/)(v::Variable,m::AbstractMonomial)=divide(monomial(v),m)
# Base.:(/)(m::AbstractMonomial,v::Variable)=divide(m,monomial(v))
# divides(u,v)=divides(monomial(u),monomial(v))
# divide(u,v;all=false)=divide(monomial(u),monomial(v), all=all)

Base.conj(v::Variable) = isnothing(v.conj[]) ? v : v.conj[]
Base.adjoint(v::Variable) = conj(v)
is_herm(v::Variable) = !isnothing(v.conj[]) ? false : true

Base.one(v::Variable) = one(v.parent_monoid[])
is_identity(v::V) where {V <: Variable} = false

exponents(v::Variable) = Tuple([v==i ? 1 : 0 for i in v.parent_monoid[].vertices])
degree(v::Variable) = 1

function looper(args)
    vars=[]
    exprs=[]
    for var in args
        if isa(var, Symbol)
            if (contains(string(var), "_"))
                var_conj=Symbol(string(var)*"_")
                push!(vars, var_conj)
                push!(vars, var)
                push!(exprs, :($(esc(var_conj)) = $Variable($"$var_conj")))
                push!(exprs, :($(esc(var)) = $Variable($"$var")))
                push!(exprs, :($(esc(var)).conj[] = $(esc(var_conj))))
                push!(exprs, :($(esc(var_conj)).conj[] = $(esc(var))))
            else
                push!(vars, var)
                push!(exprs, :($(esc(var)) = $Variable($"$var")))
            end
        else
            k=var.args[1]
            k_=Symbol(string(k)*"_")
            k__=Symbol(string(k)*"__")
            herm=var.args[2]
            non_herm=var.args[3]

            push!(
                exprs,
                :(
                    ($(esc(k_)), $(esc(k__))) =
                        $varArray($(string(k)), 1:($(esc(non_herm))), true)
                ),
            )
            push!(exprs, :($(esc(k)) = $varArray($(string(k)), 1:($(esc(herm))))))
            push!(vars, k)
            push!(vars, k_)
            push!(vars, k__)
        end
    end
    return vars, exprs
end

function varArray(s::String, indices, herm::Bool = false)
    if (herm)
        z=Variable[]
        z_=Variable[]
        for i in indices
            subs=map(k->IndexMap[k], string(i))
            push!(z, Variable("$s"*subs*"_"))
            push!(z_, Variable("$s"*subs*"__"))
        end
        [(z[i].conj[] = z_[i]; z_[i].conj[] = z[i]) for i in 1:length(z)]
        return z, z_
    else
        z=Variable[]
        for i in indices
            subs=map(k->IndexMap[k], string(i))
            push!(z, Variable("$s"*subs))
        end
        return z
    end
end

macro ncmonoid(M, args...)
    vars, exprs=looper(args)
    push!(
        exprs,
        :($(esc(M))=$NCMonoid($(string(M)), vcat($(Expr(:tuple, esc.(vars)...))...))),
    )
    return :($(foldr((x, y) -> :($x; $y), exprs, init = :())););
end

macro pcmonoid(M, args...)
    vars, exprs=looper(args)
    push!(
        exprs,
        :(
            $(esc(M))=$GraphProductMonoid(
                $(string(M)),
                vcat($(Expr(:tuple, esc.(vars)...))...),
            )
        ),
    )
    return :($(foldr((x, y) -> :($x; $y), exprs, init = :())););
end

macro subsystem(M, args...)
    monoids=[]
    exprs=[]
    systems=[i.args[1] for i in args]
    herm=[(i.args[2], i.args[3]) for i in args]
    for index in 1:length(args)
        vars=[]
        k=Symbol(lowercase(string(systems[index])))
        k_=Symbol(string(k)*"_")
        k__=Symbol(string(k)*"__")
        push!(
            exprs,
            :(
                ($(esc(k_)), $(esc(k__))) =
                    $varArray($(string(k)), $(1:herm[index][2]), true)
            ),
        )
        push!(exprs, :($(esc(k)) = $varArray($(string(k)), $(1:herm[index][1]))))
        push!(vars, k)
        push!(vars, k_)
        push!(vars, k__)
        push!(
            exprs,
            :(
                $(esc(systems[index]))=$NCMonoid(
                    $(string(systems[index])),
                    vcat($(Expr(:tuple, esc.(vars)...))...),
                )
            ),
        )
        push!(monoids, systems[index])
    end
    push!(
        exprs,
        :(
            $(esc(M))=$GraphProductMonoid(
                $(string(M)),
                vcat($(Expr(:tuple, esc.(monoids)...))...),
            )
        ),
    )
    push!(exprs, :($comms_system($(esc(M)), vcat($(Expr(:tuple, esc.(monoids)...))...))))
    return :($(foldr((x, y) -> :($x; $y), exprs, init = :())););
end

macro pcmonoid_simple(M)
    exprs=[]
    push!(
        exprs,
        :(
            ($isdefined($(esc(M)), 1)) &&
            $typeof($(esc(M))) <: GraphProductMonoid{Variable} ||
            throw("Provide a valid Monoid")
        ),
    )
    push!(
        exprs,
        :(
            ($isempty($(esc(M)).vertices) || $isempty($(esc(M)).commutations)) &&
            throw("Provide commutations")
        ),
    )
    push!(exprs, :($(esc(M))=$PCMonoid_generator($(esc(M)))))
    return :($(foldr((x, y) -> :($x; $y), exprs, init = :())););
end
function Projector(var::Variable)
    isa(var.parent_monoid[], GraphProductMonoid) ?
    (var.mult_type[] = :Projector; var'.mult_type[] = :Projector) :
    (add_relations!([var*var-var]); add_relations!([var'*var'-var']))
    nothing
end
function Unipotent(var::Variable)
    isa(var.parent_monoid[], GraphProductMonoid) ?
    (var.mult_type[] = :Unipotent; var'.mult_type[] = :Unipotent) :
    (add_relations!([var*var-1]), add_relations!([var'*var'-1]))

    nothing
end
function Unitary(var::Variable)
    is_herm(var) && throw("Provide non-hermitian variable")
    isa(var.parent_monoid[], GraphProductMonoid) ?
    (var.mult_type[] = :Unitary; (var').mult_type[] = :Unitary) :
    add_relations!([var*var'-1])

    nothing
end

variables(x) = Variable[]
