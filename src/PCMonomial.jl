"""
    PCMonomial <: AbstractMonomial

    Primary data structure for partially-commutative monomials replacing GraphProductMonoid{Variable}.
    
    This is the optimized implementation that stores variables as UInt32 Variable references,
    providing significant performance improvements for arithmetic operations.

    #Fields:
    - `monoid` : Partially-commutative monoid of the monomial.
    - `monomial` : Reference to parent monomial in nested structure.
    - `clique_words` : normal form of a PC Monomial, Vector of vectors of UInt32.
    - `edge_l` : Set of UInt32 indices for left edge variables.
    - `edge_r` : Set of UInt32 indices for right edge variables.

    #Notes:
    - For nested GraphProductMonoid structures (e.g., GraphProductMonoid{NCMonoid}),
      use GraphProductWord{T} instead.
    - PCMonomial is the default monomial type produced by one(::GraphProductMonoid{Variable}).
"""
struct PCMonomial <: AbstractMonomial
    monoid::GraphProductMonoid{Variable}
    monomial::Base.RefValue{AbstractMonomial}
    clique_words::Vector{Vector{UInt32}}
    edge_l::Set{UInt32}
    edge_r::Set{UInt32}
end

# ============================================================================
# Index Mapping Functions
# ============================================================================

"""
    var_to_index(var) -> UInt32

Convert a Variable to its UInt32 index in monoid.vertices.
Uses cached dictionary for O(1) lookup if available (after build()).
Falls back to linear search during construction or if dictionary not available.
"""
@inline function var_to_index(var::Variable)::UInt32
    # Fast path: use cached dictionary (O(1) lookup)
    monoid = var.parent_monoid[]
    dict = monoid.var_index_dict[]
    if dict !== nothing
        return dict[var]
    end
    # Fallback: linear search (used during build() before dictionary is created)
    @inbounds for i in 1:length(monoid.vertices)
        monoid.vertices[i] === var && return UInt32(i)
    end
    # Final fallback with equality check
    return UInt32(findfirst(==(var), monoid.vertices))
end

"""
    index_to_var(monoid, idx) -> Variable

Convert a UInt32 index back to a Variable.
"""
@inline function index_to_var(monoid::GraphProductMonoid{Variable}, idx::UInt32)::Variable
    @inbounds return monoid.vertices[idx]
end

"""
    build_var_index_dict(monoid) -> Dict{Variable, UInt32}

Build a dictionary mapping Variables to their indices for O(1) lookup.
Should be called once during build() and stored in the monoid.
Note: This is now called automatically in build() and cached in monoid.var_index_dict[].
"""
function build_var_index_dict(monoid::GraphProductMonoid{Variable})
    return Dict{Variable, UInt32}(v => UInt32(i) for (i, v) in enumerate(monoid.vertices))
end

"""
    build_conj_indices(monoid) -> Vector{UInt32}

Build conjugate index mapping. conj_indices[i] = index of conj(vertices[i]).
If variable is hermitian, conj_indices[i] = i.
Note: This is now called automatically in build() and cached in monoid.conj_indices[].
"""
function build_conj_indices(monoid::GraphProductMonoid{Variable})
    # Use cached version if available
    cached = monoid.conj_indices[]
    if cached !== nothing
        return cached
    end

    n = length(monoid.vertices)
    conj_indices = Vector{UInt32}(undef, n)

    # Use var_index_dict for O(1) lookup if available
    dict = monoid.var_index_dict[]

    @inbounds for i in 1:n
        v = monoid.vertices[i]
        cv = conj(v)
        if cv === v || cv == v
            conj_indices[i] = UInt32(i)
        elseif dict !== nothing
            # O(1) lookup using dictionary
            conj_indices[i] = dict[cv]
        else
            # Fallback: linear search
            for j in 1:n
                if monoid.vertices[j] === cv || monoid.vertices[j] == cv
                    conj_indices[i] = UInt32(j)
                    break
                end
            end
        end
    end
    return conj_indices
end

# ============================================================================
# Conversion Functions: GraphProductWord{Variable} ↔ PCMonomial
# ============================================================================

"""
    PCMonomial(m::GraphProductWord{Variable}) -> PCMonomial

Convert a GraphProductWord{Variable} to PCMonomial.
"""
function PCMonomial(m::GraphProductWord{Variable})
    monoid = m.monoid
    n_cliques = length(m.clique_words)

    # Convert clique_words
    clique_words = Vector{Vector{UInt32}}(undef, n_cliques)
    @inbounds for i in 1:n_cliques
        cw = m.clique_words[i]
        ci = var_to_index.(cw)
        clique_words[i] = ci
    end

    # Convert edge sets
    edge_l = Set{UInt32}(var_to_index(v) for v in m.edge_l)
    edge_r = Set{UInt32}(var_to_index(v) for v in m.edge_r)

    return PCMonomial(
        monoid,
        Base.RefValue{AbstractMonomial}(),
        clique_words,
        edge_l,
        edge_r,
    )
end

"""
    GraphProductWord(m::PCMonomial) -> GraphProductWord{Variable}

Convert a PCMonomial back to GraphProductWord{Variable}.
"""
function GraphProductWord(m::PCMonomial)
    monoid = m.monoid
    n_cliques = length(m.clique_words)

    # Convert clique_words back to clique_words
    clique_words = Vector{Vector{Variable}}(undef, n_cliques)
    @inbounds for i in 1:n_cliques
        ci = m.clique_words[i]
        cw = Vector{Variable}(undef, length(ci))
        for j in eachindex(ci)
            cw[j] = index_to_var(monoid, ci[j])
        end
        clique_words[i] = cw
    end

    # Convert edge sets
    edge_l = Set{Variable}(index_to_var(monoid, idx) for idx in m.edge_l)
    edge_r = Set{Variable}(index_to_var(monoid, idx) for idx in m.edge_r)

    return GraphProductWord(
        monoid,
        Base.RefValue{AbstractMonomial}(),
        clique_words,
        edge_l,
        edge_r,
    )
end

# ============================================================================
# Basic Operations
# ============================================================================

Base.length(m::PCMonomial) = length(m.clique_words)

# Iteration yields Variables for user convenience
function Base.iterate(m::PCMonomial)
    isempty(m.clique_words) && return nothing
    return (get_clique_word(m, 1), 1)
end

function Base.iterate(m::PCMonomial, state::Int)
    state >= length(m.clique_words) && return nothing
    return (get_clique_word(m, state + 1), state + 1)
end

"""
    get_clique_word(m::PCMonomial, i::Int) -> Vector{Variable}

Get the i-th clique word as a Vector{Variable}.
"""
function get_clique_word(m::PCMonomial, i::Int)
    @inbounds ci = m.clique_words[i]
    return [index_to_var(m.monoid, idx) for idx in ci]
end

Base.getindex(m::PCMonomial, i::Int) = get_clique_word(m, i)

function Base.copy(m::PCMonomial)
    PCMonomial(
        m.monoid,
        Base.RefValue{AbstractMonomial}(),
        [copy(ci) for ci in m.clique_words],
        copy(m.edge_l),
        copy(m.edge_r),
    )
end

# ============================================================================
# Identity
# ============================================================================

"""
    one(M::GraphProductMonoid{Variable}) -> PCMonomial

Returns the identity element of the monoid as a PCMonomial.
"""
function Base.one(M::GraphProductMonoid{Variable})
    n_cliques = length(M.cliques)
    return PCMonomial(
        M,
        Base.RefValue{AbstractMonomial}(),
        [UInt32[] for _ in 1:n_cliques],
        Set{UInt32}(),
        Set{UInt32}(),
    )
end

Base.one(m::PCMonomial) = one(m.monoid)

function is_identity(m::PCMonomial)
    @inbounds for ci in m.clique_words
        !isempty(ci) && return false
    end
    return true
end

# ============================================================================
# Equality
# ============================================================================

function Base.:(==)(m::PCMonomial, n::PCMonomial)
    m.monoid != n.monoid && return false
    return m.clique_words == n.clique_words
end

Base.:(==)(m::PCMonomial, n::Number) = n == 1 && is_identity(m)
Base.:(==)(n::Number, m::PCMonomial) = m == n

# ============================================================================
# Hashing and Comparison
# ============================================================================

Base.hash(m::PCMonomial, h::UInt) =
    hash(m.monoid, hash(m.clique_words, hash(0x8a7b9c3d4e5f6012, h)))

function less_or_not(m::PCMonomial, n::PCMonomial)
    m.monoid != n.monoid &&
        throw(ArgumentError("Cannot compare monomials from different monoids"))
    m == n && return false

    deg_m = degree(m)
    deg_n = degree(n)
    deg_m < deg_n && return true
    deg_m > deg_n && return false

    @inbounds for i in 1:length(m.clique_words)
        m.clique_words[i] == n.clique_words[i] && continue
        len_m = length(m.clique_words[i])
        len_n = length(n.clique_words[i])
        len_m < len_n && return true
        len_m > len_n && return false
        return m.clique_words[i] < n.clique_words[i]
    end
    return false
end

Base.isless(m::PCMonomial, n::PCMonomial) = less_or_not(m, n)
Base.:<(m::PCMonomial, n::PCMonomial) = less_or_not(m, n)

# ============================================================================
# Display - Shows Variables not indices
# ============================================================================

function Base.show(io::IO, m::PCMonomial)
    if all(isempty.(m.clique_words))
        print(io, "Id")
        return
    end
    show_level = m.monoid.show_level[]
    if show_level == 1
        for i in 1:length(m.clique_words)
            vars = get_clique_word(m, i)
            print(io, "($(join(vars, ",")))")
        end
    elseif show_level == 2
        for i in 1:length(m.clique_words)
            vars = get_clique_word(m, i)
            if !isempty(vars)
                print(io, "($(join(vars, ",")))")
            end
        end
    else
        mw = monomial_to_word(m)
        print(io, join(mw))
    end
end

Base.show(io::IO, mime::MIME"text/plain", m::PCMonomial) = show(io, m)
Base.show(io::IO, mime::MIME"text/print", m::PCMonomial) = show(io, m)

# ============================================================================
# Fast Merge (Core Performance Improvement)
# ============================================================================

"""
    merge_clique_words(x, y) -> Vector{Vector{UInt32}}

Fast merge of clique words using UInt32 arrays.
"""
function merge_clique_words(x::Vector{Vector{UInt32}}, y::Vector{Vector{UInt32}})
    n = length(x)
    result = Vector{Vector{UInt32}}(undef, n)
    @inbounds for i in 1:n
        xi, yi = x[i], y[i]
        len_x, len_y = length(xi), length(yi)
        ri = Vector{UInt32}(undef, len_x + len_y)
        copyto!(ri, 1, xi, 1, len_x)
        copyto!(ri, len_x + 1, yi, 1, len_y)
        result[i] = ri
    end
    return result
end

"""
    merge_clique_words!(x, y) -> Vector{Vector{UInt32}}

In-place merge - appends y to x. Use when x is already a copy.
"""
function merge_clique_words!(x::Vector{Vector{UInt32}}, y::Vector{Vector{UInt32}})
    @inbounds for i in eachindex(x)
        append!(x[i], y[i])
    end
    return x
end

# ============================================================================
# Edge Variable Computation
# ============================================================================

"""
    get_edge_indices(clique_words, position, clique_idx_list, monoid) -> Set{UInt32}

Get edge variable indices from clique_words.
"""
function get_edge_indices(
    clique_words::Vector{Vector{UInt32}},
    position::Symbol,
    clique_idx_list::Vector{Int},
    monoid::GraphProductMonoid{Variable},
)
    result = Set{UInt32}()

    @inbounds for ci in clique_idx_list
        isempty(clique_words[ci]) && continue
        idx = position == :first ? clique_words[ci][1] : clique_words[ci][end]

        # Check if this is truly an edge (appears at position in all its cliques)
        var = monoid.vertices[idx]
        is_edge = true
        for var_ci in var.clique_indices
            if isempty(clique_words[var_ci])
                is_edge = false
                break
            end
            check_idx =
                position == :first ? clique_words[var_ci][1] : clique_words[var_ci][end]
            if check_idx != idx
                is_edge = false
                break
            end
        end
        is_edge && push!(result, idx)
    end
    return result
end

function get_edge_indices(
    clique_words::Vector{Vector{UInt32}},
    position::Symbol,
    monoid::GraphProductMonoid{Variable},
)
    return get_edge_indices(clique_words, position, collect(1:length(clique_words)), monoid)
end

# ============================================================================
# Multiplication
# ============================================================================

"""
    simple_multiply(m::PCMonomial, n::PCMonomial) -> PCMonomial

Simple multiplication - just merges clique indices. Used when monoid has no relations.
"""
function simple_multiply(m::PCMonomial, n::PCMonomial)
    m.monoid != n.monoid &&
        throw(ArgumentError("Cannot multiply monomials from different monoids"))
    monoid = m.monoid
    res_clique_words = merge_clique_words(m.clique_words, n.clique_words)
    # edge_l = get_edge_indices(res_clique_words, :first, monoid)
    # edge_r = get_edge_indices(res_clique_words, :last, monoid)
    # return PCMonomial(m.monoid, Base.RefValue{AbstractMonomial}(), 
    #                   res_clique_words, edge_l, edge_r)
    return PCMonomial(
        m.monoid,
        Base.RefValue{AbstractMonomial}(),
        res_clique_words,
        Set{UInt32}(),
        Set{UInt32}(),
    )
end

"""
    multiply(m::PCMonomial, n::PCMonomial) -> PCMonomial

Full multiplication with Projector/Unipotent/Unitary handling.
"""
function multiply(m::PCMonomial, n::PCMonomial)
    # m.monoid != n.monoid && throw(ArgumentError("Cannot multiply monomials from different monoids"))
    if m.monoid!=n.monoid
        return general_mult(m, n)
    end
    monoid = m.monoid
    if !monoid.has_relations[]
        return simple_multiply(m, n)
    end

    m_indices = [copy(ci) for ci in m.clique_words]
    n_indices = [copy(ci) for ci in n.clique_words]
    m_l = copy(m.edge_l)
    m_r = copy(m.edge_r)
    n_l = copy(n.edge_l)
    n_r = copy(n.edge_r)

    # Check for edge interactions
    middle_list=Vector{}()
    zero_result = false
    for i_idx in m_r
        for j_idx in n_l
            i_var = monoid.vertices[i_idx]
            j_var = monoid.vertices[j_idx]

            # Check orthogonality
            if j_var in i_var.ortho_conj
                zero_result = true
                break
            end
            # Projector: P * P = P
            if i_idx == j_idx && i_var.mult_type[] == :Projector
                var_clique_indices = i_var.clique_indices
                for ci in var_clique_indices
                    pop!(m_indices[ci])
                    popfirst!(n_indices[ci])
                end
                push!(middle_list, i_var)
            end

            # Unipotent: U * U = 1
            if i_idx == j_idx && i_var.mult_type[] == :Unipotent
                var_clique_indices = i_var.clique_indices
                for ci in var_clique_indices
                    pop!(m_indices[ci])
                    popfirst!(n_indices[ci])
                end
                push!(middle_list, one(i_var))
            end

            # Unitary: U * U† = 1
            if i_var.mult_type[] == :Unitary
                conj_i = conj(i_var)
                if j_var == conj_i
                    var_clique_indicesi = i_var.clique_indices
                    var_clique_indicesj = j_var.clique_indices
                    for ci in var_clique_indicesi
                        pop!(m_indices[ci])
                    end
                    for ci in var_clique_indicesj
                        popfirst!(n_indices[ci])
                    end
                    push!(middle_list, one(i_var))
                end
            end
        end
        zero_result && break
    end

    # Return zero polynomial if orthogonal
    zero_result && return Polynomial(monoid)

    if isempty(middle_list)
        merge_clique_words!(m_indices, n_indices)
        new_edge_l = get_edge_indices(m_indices, :first, monoid)
        new_edge_r = get_edge_indices(m_indices, :last, monoid)
        return PCMonomial(
            monoid,
            Base.RefValue{AbstractMonomial}(),
            m_indices,
            new_edge_l,
            new_edge_r,
        )
    else
        # Recompute edges after modifications
        new_m_edge_l = get_edge_indices(m_indices, :first, monoid)
        new_m_edge_r = get_edge_indices(m_indices, :last, monoid)
        new_n_edge_l = get_edge_indices(n_indices, :first, monoid)
        new_n_edge_r = get_edge_indices(n_indices, :last, monoid)

        m_ = PCMonomial(
            monoid,
            Base.RefValue{AbstractMonomial}(),
            m_indices,
            new_m_edge_l,
            new_m_edge_r,
        )
        n_ = PCMonomial(
            monoid,
            Base.RefValue{AbstractMonomial}(),
            n_indices,
            new_n_edge_l,
            new_n_edge_r,
        )

        middle_prod=prod(middle_list)

        return m_ * middle_prod * n_
    end
end

Base.:*(m::PCMonomial, n::PCMonomial) = multiply(m, n)

"""
    idx_to_monomial(monoid, idx) -> PCMonomial

Create a monomial from a single variable index.
"""
function idx_to_monomial(monoid::GraphProductMonoid{Variable}, idx::UInt32)
    var = monoid.vertices[idx]
    return words_to_monomial(monoid, [var])
end

# ============================================================================
# Conjugation
# ============================================================================

"""
    conjugate(m::PCMonomial, conj_indices::Vector{UInt32}) -> PCMonomial

Conjugation using precomputed conjugate indices.
"""
function conjugate(m::PCMonomial, conj_indices::Vector{UInt32})
    is_identity(m) && return m

    n_cliques = length(m.clique_words)
    new_clique_words = Vector{Vector{UInt32}}(undef, n_cliques)

    @inbounds for i in 1:n_cliques
        ci = m.clique_words[i]
        if isempty(ci)
            new_clique_words[i] = UInt32[]
        else
            # Reverse and conjugate each index
            new_ci = Vector{UInt32}(undef, length(ci))
            for j in eachindex(ci)
                new_ci[j] = conj_indices[ci[end - j + 1]]
            end
            new_clique_words[i] = new_ci
        end
    end

    # Swap and conjugate edges
    new_edge_l = Set{UInt32}(conj_indices[idx] for idx in m.edge_r)
    new_edge_r = Set{UInt32}(conj_indices[idx] for idx in m.edge_l)

    return PCMonomial(
        m.monoid,
        Base.RefValue{AbstractMonomial}(),
        new_clique_words,
        new_edge_l,
        new_edge_r,
    )
end

function conjugate(m::PCMonomial)
    is_identity(m) && return m

    # Use cached conj_indices from monoid (built during build())
    # Falls back to building on demand if not cached
    conj_indices = build_conj_indices(m.monoid)
    return conjugate(m, conj_indices)
end

function Base.conj(m::PCMonomial)
    is_identity(m) && return m
    m.monoid.conj_type[] && return conjugate(m)
    w=monomial_to_word(m)
    return prod(conj.(reverse(w)))
end
Base.adjoint(m::PCMonomial) = conj(m)

# ============================================================================
# Division
# ============================================================================

"""
    divide(m::PCMonomial, n::PCMonomial; all=false)

Division: find (l, r) such that l * m * r = n.
Delegates to GraphProductWord division for correctness.
"""
function divide(m::PCMonomial, n::PCMonomial; all = false)
    m == n && return true, (one(m.monoid), one(n.monoid))

    # Convert to GraphProductWord, use existing division, convert back
    m_gpw = GraphProductWord(m)
    n_gpw = GraphProductWord(n)

    ok, result_gpw = divide(m_gpw, n_gpw; all = all)
    if !ok
        return false, (nothing, nothing)
    end
    # Convert results back to PCMonomial
    if all
        return true, [(PCMonomial(l), PCMonomial(r)) for (l, r) in result_gpw]
    else
        (l, r) = result_gpw
        return true, (PCMonomial(l), PCMonomial(r))
    end
end

Base.:/(m::PCMonomial, n::PCMonomial) = divides(m, n)
divides(a::PCMonomial, b::PCMonomial) = (divide(a, b))[1]

# ============================================================================
# Degree and Exponents
# ============================================================================

function degree(m::PCMonomial)
    return sum(values(exponents(m)); init = 0)
end

function exponents(m::PCMonomial)
    monoid = m.monoid
    result = Dict{Variable, Int}()

    for (var_idx, var) in enumerate(monoid.vertices)
        if !isempty(var.clique_indices)
            # Count occurrences in first clique this var belongs to
            ci = var.clique_indices[1]
            cnt = Base.count(==(UInt32(var_idx)), m.clique_words[ci])
            result[var] = cnt
        end
    end

    return result
end

function variables(m::PCMonomial)
    return [var for (var, cnt) in exponents(m) if cnt > 0]
end

monomials(m::PCMonomial) = [m]

# ============================================================================
# monomial() function - find parent monomial in nested structure
# ============================================================================

function monomial(m::PCMonomial)
    if isdefined(m.monomial, :x)
        return m.monomial[]
    end

    monoid = m.monoid
    if isa(monoid.parent_monoid[], GraphProductMonoid)
        parent_monoid = monoid.parent_monoid[]
        m.monomial[] = words_to_monomial(parent_monoid, [m])
        return m.monomial[]
    else
        throw(
            ArgumentError(
                "Cannot convert to monomial, parent monoid is not a GraphProductMonoid",
            ),
        )
    end
end

# ============================================================================
# Zero
# ============================================================================

Base.zero(m::PCMonomial) = 0

# ============================================================================
# Interoperability with GraphProductWord{Variable}
# ============================================================================

# Allow multiplication between PCMonomial and GraphProductWord{Variable}
Base.:*(m::PCMonomial, n::GraphProductWord{Variable}) = m * PCMonomial(n)
Base.:*(m::GraphProductWord{Variable}, n::PCMonomial) = PCMonomial(m) * n

# Equality between PCMonomial and GraphProductWord{Variable}
Base.:(==)(m::PCMonomial, n::GraphProductWord{Variable}) = m == PCMonomial(n)
Base.:(==)(m::GraphProductWord{Variable}, n::PCMonomial) = PCMonomial(m) == n

# Comparison
Base.isless(m::PCMonomial, n::GraphProductWord{Variable}) = m < PCMonomial(n)
Base.isless(m::GraphProductWord{Variable}, n::PCMonomial) = PCMonomial(m) < n

# ============================================================================
# Interoperability with Variable
# ============================================================================

Base.:*(v::Variable, m::PCMonomial) = monomial(v) * m
Base.:*(m::PCMonomial, v::Variable) = m * monomial(v)
