include("pcpop.jl")

struct TraceMonoid <: AbstractMonoid
    state_monoid::AbstractMonoid
    base_monoid::AbstractMonoid
    vertices_free::Vector
    vertices_states::Vector
    dict_free::Dict
    dict_states::Dict
    tracial::Bool
end

function TraceMonoid(; state_monoid, base_monoid, vertices_free, vertices_states, dict_free, dict_states, tracial)
    return TraceMonoid(state_monoid, base_monoid, vertices_free, vertices_states, dict_free, dict_states, tracial)
end

function Base.show(io::IO, ::MIME"text/plain", tm::TraceMonoid)
    println(io, "TraceMonoid with:")
    println(io, "  Base monoid: ", typeof(tm.base_monoid))
    println(io, "  State monoid: ", typeof(tm.state_monoid))
    println(io, "  Free generators: ", tm.vertices_free)
    println(io, "  State generators: ", tm.vertices_states)
    println(io, "  Tracial: ", tm.tracial)
end

# Display cyclic words
function Base.show(io::IO, w::CyclicWord)
    print(io, "[$(GraphProductWord(w.ref_word))]")
end

function Base.show(io::IO, ::MIME"text/plain", w::CyclicWord)
    print(io, "[$(GraphProductWord(w.ref_word))]")
end

# Coefficient CyclicWord
function coefficient(p::Polynomial, b::CyclicWord)
    i = findall(x -> x == b, cyclic_reduce.(p.monomials))
    if isempty(i)
        return 0
    else
        return sum(p.coeffs[i])
    end
end

# Trace Words
struct TraceWord{T<:AbstractMonomial} <: AbstractMonomial
    monoid::TraceMonoid
    free::T
    states::T
end

# Display Trace Words
function Base.show(io::IO, w::TraceWord)
    println(io, w.free, w.states)
end

Base.:*(m1::TraceWord,m2::TraceWord) = trace_mult(m1, m2)
Base.:(==)(m1::TraceWord,m2::TraceWord) = trace_equals(m1,m2)

function trace_mult(m1::TraceWord,m2::TraceWord)
    if m1.monoid !== m2.monoid
        error("These words do not belong to the same monoid!")
    else
        return TraceWord(
            m1.monoid,
            m1.free*m2.free,
            m1.states*m2.states
        )
    end
end

function trace_equals(m1::TraceWord,m2::TraceWord)
    m1.monoid == m2.monoid && m1.free == m2.free && m1.states == m2.states
end



function list_projectors(M)
    projectors = []
    for var in M.vertices
        if var.mult_type[] == :Projector
            push!(projectors, var)
        end
    end
    return projectors
end

function list_unitaries(M)
    unitaries = []
    for var in M.vertices
        if var.mult_type[] == :Unitary
            push!(unitaries, var)
        end
    end
    return unitaries
end

function list_unipotents(M)
    unipotents = []
    for var in M.vertices
        if var.mult_type[] == :Unipotent
            push!(unipotents, var)
        end
    end
    return unipotents
end

"""
    function make_trace_monoid(M::AbstractMonoid, k::Int; tracial=false)

    Build monoid of traces in ``M´´ up to degree ``k´´.

    #Input:
    - ``M::AbstractMonoid´´ : Base monoid
    - ``k::Integer´´ : Degree of traces
    - ``tracial::Boolean´´ : Tracial state

    #Output:
    - ``TM::AbstractMonoid´´ : Monoid of traces in M up to degree k.

    ``TM´´ has the variables in ``M´´ plus one commutative variable ``tr(w)´´
    for each monomial ``w´´ in ``M´´ with degree up to ``k´´
"""
function make_trace_monoid(M::AbstractMonoid, k::Int; statesymbol="ρ", monomialsymbol="μ", tracial=false)
    traces = mons_at_level(M, k)
    if tracial
        traces = unique([cyclic_reduce(t) for t in traces])
    end
    num_m = length(M.vertices)
    num_t = length(traces)
    @eval @pcmonoid TM $(Symbol.(monomialsymbol, M.vertices)...) $(Symbol.(statesymbol, traces)...)
    μ = TM.vertices[1:num_m]
    ρ = TM.vertices[num_m+1:end]
    dict_monomials = Dict{AbstractMonomial, Variable}(M.vertices[i] => μ[i] for i in 1:num_m)
    dict_traces = Dict{AbstractMonomial, Variable}(traces[i] => ρ[i] for i in 1:num_t)
    for (m1, m2) in M.commutations
        @comms [dict_monomials[m1], dict_monomials[m2]]
    end
    @comms [μ, ρ]
    @comms [ρ, ρ]

    Projector.([dict_monomials[m] for m in list_projectors(M)])
    Unipotent.([dict_monomials[m] for m in list_unipotents(M)])
    Unitary.([dict_monomials[m] for m in list_unitaries(M)])

    # add conjugation relations 
    # ρ[w]' = ρ[w']
    # μ[a]' = μ[a']
    for a in M.vertices
        dict_monomials[a].conj[] = dict_monomials[a']
    end

    for t in traces
        dict_traces[t].conj[] = dict_traces[t']
    end

    build(TM)

    return TraceMonoid(
        state_monoid = TM,
        base_monoid = M,
        vertices_free = μ,
        vertices_states = ρ,
        dict_free = dict_monomials,
        dict_states = dict_traces,
        tracial = tracial
    )
end

"""

    Return trace monomials in ``M´´ up to degree ``k´´

    trace_monomials(M::AbstractMonoid, k::Int; tracial=false)

    #Input:
    - ``M::AbstractMonoid´´ : Base Monoid
    - ``k::Int´´ : Degree trace monomials

    #Output
    - Trace monomials in ``M´´ up to degree ``k´´.

    Example:

    M = < a, b>

    Trace monomials degree k = 1:

    a, b, ρa, ρb

    Trace monomials degree k = 2:

    aa, ab, ba, bb, ρaa, ρab, ρba, ρbb, 
    aρa, aρb, bρa, bρb, ρaρa, ρaρb, ρbρb

"""
function trace_monomials(TM::TraceMonoid, k::Int; tracial=false)
    return mons_at_level(TM.state_monoid, k)
end

function trace_monomials(TM::TraceMonoid, k::Int, t::Int; tracial=false)
    return vcat(mons_at_level(TM.vertices_free, k), mons_at_level(TM.vertices_states, t))
end

degree(m::CyclicWord) = degree(m.ref_word)
monomial_to_word(m::CyclicWord) = monomial_to_word(m.ref_word)

using Combinatorics: partitions
function pure_trace_monomials(TM::TraceMonoid, k::Int; tracial=false)
    all_monomials = []
    base_monomials = mons_at_level(TM.base_monoid.vertices, k)
    if tracial
        base_monomials = unique(cyclic_reduce.(base_monomials))
    end
    dict_monomials = Dict(n => [m for m in base_monomials if (degree(m) == n)] for n in 0:k)
    # ρ(w1) ... ρ(wl) where |w1| + ... + |wl| = k
    for n in 1:k
        for part in partitions(k, n)
            for w_I in Iterators.product([dict_monomials[i] for i in part]...)
                append!(all_monomials, [prod([state(w, TM) for w in w_I])])
            end
        end
    end
    all_monomials = unique([clean_one(m, TM) for m in all_monomials])
    return all_monomials
end

function trace_monomials(TM::TraceMonoid, k::UnitRange{Int64}; tracial=false, pure=false)
    union([trace_monomials(TM, n, tracial=tracial, pure=pure) for n in k]...)
end

function trace_monomials(TM::TraceMonoid, k::Int64; tracial=false, pure=false)
    if pure
        return pure_trace_monomials(TM, k, tracial=tracial)   
    else
        base_monomials = mons_at_level(TM.base_monoid.vertices, k)
        if tracial
            base_monomials = cyclic_reduce.(base_monomials)
        end
        dict_monomials = Dict(n => [m for m in base_monomials if (degree(m) == n)] for n in 0:k)
        dict_states = Dict(n => pure_trace_monomials(TM, n, tracial=tracial) for n in 0:k)
    end
    
    all_monomials = []
    for n in 0:k
        append!(all_monomials, [state_embedding(w0, TM)*w1 for w0 in dict_monomials[n] for w1 in dict_states[k-n]])
    end

    all_monomials = unique([clean_one(m, TM) for m in all_monomials])
    return all_monomials
end

# Transform TraceWord to AbstractWord
function state_to_word(word::TraceWord)
    return word.free*word.states
end

function free_part(word::AbstractMonomial, TM::TraceMonoid)
    word_free = [w for w in monomial_to_word(word) if w in TM.vertices_free]
    word_free = words_to_monomial(TM.state_monoid, word_free)
end

function state_part(word::AbstractMonomial, TM::TraceMonoid)
    word_states = [w for w in monomial_to_word(word) if w in TM.vertices_states]
    word_states = words_to_monomial(TM.state_monoid, word_states)
end

# Transform AbstractWord to TraceWord
function word_to_state(word::AbstractMonomial, TM::TraceMonoid)
    return TraceWord(
        TM,
        free_part(word, TM),
        state_part(word, TM)
        )
end

# Embed word to StateWord
function state_embedding(word::AbstractMonomial, TM::TraceMonoid)
    if word.monoid == TM.state_monoid
        return word
    elseif word.monoid !== TM.base_monoid
        error("That word does not belong to the base monoid!")
    end

    dict_free = TM.dict_free
    base_word =  words_to_monomial(TM.state_monoid, [dict_free[a] for a in monomial_to_word(word)])
end

function state_embedding(p::Polynomial, TM::TraceMonoid)
    sum([c*state_embedding(m, TM) for (c, m) in zip(p.coeffs, p.monomials)])
end

"""
    state_projection(word::TraceWord)
    state_projection(word::AbstractMonomial, TM::TraceMonoid)

    w0 ρw1 ... ρwn --> ρw0 ρw1 ... ρwn

    State projection of a Trace Word in a Trace Monoid
"""
function state_projection(word::AbstractMonomial, TM::TraceMonoid)
    dict_free = Dict(value => key for (key, value) in TM.dict_free)
    base_word =  words_to_monomial(TM.base_monoid, [dict_free[a] for a in monomial_to_word(free_part(word, TM))])
    if base_word == one(TM.base_monoid)
        return state_part(word, TM)
    elseif TM.tracial
        return TM.dict_states[cyclic_reduce(base_word)]*state_part(word, TM)
    else
        return TM.dict_states[base_word]*state_part(word, TM)
    end
end

function state_projection(word::TraceWord)
    TM = word.monoid
    dict_free = Dict(value => key for (key, value) in TM.dict_free)
    base_word =  words_to_monomial(TM.base_monoid, [dict_free[a] for a in monomial_to_word(word.free)])

    if base_word == one(TM.base_monoid)
        return word
    elseif TM.tracial
        return TraceWord(
        TM,
        one(TM.state_monoid),
        TM.dict_states[cyclic_reduce(base_word)]*word.states
        )
    else
        return TraceWord(
        TM,
        one(TM.state_monoid),
        TM.dict_states[base_word]*word.states
        )
    end
end

function state_projection(p::Polynomial, TM::TraceMonoid)
    return sum([c*state_projection(m, TM) for (c, m) in zip(p.coeffs, p.monomials)])
end

function state_projection(p::Polynomial)
    return sum([c*state_projection(m) for (c, m) in zip(p.coeffs, p.monomials)])
end

"""
    Get the formal state expectation value of a monomial.

    state(m::AbstractMonomial, TM::TraceMonoid)
    state(m::Polynomial, TM::TraceMonoid)

    # INPUT
    - m : monomial in monoid M
    - TM : Trace monoid with base monoid M

    # OUTPUT
    - ρ(μm) : state monomial in TM, where μm projects m into TM
"""
function state(m::AbstractMonomial, TM::TraceMonoid)
    state_projection(state_embedding(m, TM), TM)
end

function state(m::Variable, TM::TraceMonoid)
    state_projection(state_embedding(monomial(m), TM), TM)
end

function state(p::Polynomial, TM::TraceMonoid)
    return sum([c*state(m, TM) for (c, m) in zip(p.coeffs, p.monomials)])
end

"""
    Trace polynomial optimization problem in trace monoid `TΜ`.
    
    tpop(poly::Polynomial, TM::TraceMonoid, k::Int; equalities = [], truncate = "degree", tracial=false)

    #Arguments:
    - `poly` : trace polynomial `p` in trace monoid `TΜ` with degree d.
    - `TM` : trace monoid
    - `k` : level of relaxation
            when not specified set to smallest size d÷2.
    - `equalities` : list with polynomial constraints.
    
    # Output:
    - JuMP model (unsolved) with SDP relaxation POP.
    
    Moment relaxation Γ(u,v) = [u*v] cyclic class
    [w] = [σ(w)] for any cycle σ ∈ Sn 
    
    Max sum_i p(i) Γ(i)
    s.t. Γ(1) = 1
         Γ in PSD(k)
    
    Dual of sum of square decomposition:
    
    Min t
    s.t. [t - p] in [TSOS(k)]    
"""
function tpop(poly::Polynomial, TM::TraceMonoid, basis_psd; equalities = [], truncate = "degree", tracial=false)
    if isempty(equalities)  
        matrix_psd = [clean_one(state_projection(x'*y, TM), TM) for x in basis_psd, y in basis_psd]
    else
        max_degree = maximum([degree(g) for g in equalities])
        if truncate == "degree"  
            truncate = max_degree
        elseif truncate < max_degree
            throw(ArgumentError("Truncation degree $(truncate) expected at least constraints degree $(max_degree)"))
        end
        # Extend equalities to state polynomials
        append!(equalities, [state_projection(r, TM) for r in equalities])
        grobner_truncated = macaulay_grobner(equalities, truncate)
        matrix_psd = [state_projection(reduce_grobner(Polynomial(x'*y), grobner_truncated), TM) for x in basis_psd, y in basis_psd]
        matrix_psd = reduce_duplicates(matrix_psd)
    end   
    basis_constraints = StarAlgebras.Basis{UInt16}(
              unique(matrix_psd),
            )
    Γ = [basis_constraints[m] for m in matrix_psd]

    sos_model = JuMP.Model()
    JuMP.@variable sos_model t
    JuMP.@objective sos_model Min t
    n = length(basis_psd)
    P = JuMP.@variable sos_model P[1:n, 1:n] Symmetric
    JuMP.@constraint sos_model P in PSDCone()

    objective = state_projection(t*one(TM.state_monoid)-poly, TM)
    for (idx, b) in enumerate(basis_constraints)
        if typeof(b) <: AbstractMonomial
            c = coefficient(objective, b)
            JuMP.@constraint sos_model LinearAlgebra.dot(P, Γ .== idx) == c
        elseif typeof(b) <: Polynomial
            for (b_coeff, b_monomial) in zip(b.coeffs, b.monomials)
                c = b_coeff*coefficient(objective, b_monomial)
                JuMP.@constraint sos_model LinearAlgebra.dot(P, Γ .== idx) == c
            end
        end
    end

    return sos_model
end

# TODO: ρ() == 1
function clean_one(m::AbstractMonomial, TM::TraceMonoid)
    #id = state_projection
    #m_word = monomial_to_word(m)
    #filter!(x -> x !== TM.vertices_states[1], m_word)
    #return words_to_monomial(TM.state_monoid, m_word)
    return m
end

