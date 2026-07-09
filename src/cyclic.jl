struct CyclicWord <: AbstractMonomial
    monoid::AbstractMonoid
    exponents::Dict
    ref_word::PCMonomial
    reduced_word::PCMonomial
end

function cyclic_reduce(m::PCMonomial)
    monoid = m.monoid

    m_words = [copy(ci) for ci in m.clique_words]
    m_l = copy(m.edge_l)  # Set{UInt32}
    m_r = copy(m.edge_r)  # Set{UInt32}

    conj_indices = build_conj_indices(monoid)

    for (i_idx, j_idx) in product(m_r, m_l)
        i_var = monoid.vertices[i_idx]
        j_var = monoid.vertices[j_idx]

        if j_var in i_var.ortho_conj
            return 0
        end

        # projector
        if i_idx == j_idx &&
           i_var.mult_type[] == :Projector &&
           length(m_words[i_var.clique_indices[1]]) != 1
            for index in i_var.clique_indices
                pop!(m_words[index])
            end
            delete!(m_r, i_idx)
        end

        # unitary and unipotent
        conj_i_idx = conj_indices[i_idx]
        if (i_idx == j_idx && i_var.mult_type[] == :Unipotent) ||
           (conj_i_idx == j_idx && i_var.mult_type[] == :Unitary)
            if length(m_words[i_var.clique_indices[1]]) == 1
                continue
            end
            for index in i_var.clique_indices
                pop!(m_words[index])
                popfirst!(m_words[index])
            end
            delete!(m_r, i_idx)
            delete!(m_l, j_idx)
            union!(m_l, get_edge_indices(m_words, :first, i_var.clique_indices, monoid))
            union!(m_r, get_edge_indices(m_words, :last, i_var.clique_indices, monoid))
            return cyclic_reduce(
                PCMonomial(monoid, Base.RefValue{AbstractMonomial}(), m_words, m_l, m_r),
            )
        end
    end

    reduced_word = PCMonomial(monoid, Base.RefValue{AbstractMonomial}(), m_words, m_l, m_r)
    exp = exponents(reduced_word)

    return CyclicWord(monoid, exp, m, reduced_word)
end

function cyclic_reduce(poly::Polynomial{C_T, PCMonomial}) where {C_T}
    if iszero(poly)
        return Polynomial(Vector{CyclicWord}([]), Vector{C_T}([]), poly.monoid)
    end
    r=Polynomial(cyclic_reduce(first(poly.monomials)), first(poly.coeffs))
    for i in 2:length(poly.monomials)
        c_mon=cyclic_reduce(poly.monomials[i])
        index=findfirst(x->x==c_mon, r.monomials)
        if index===nothing
            push!(r.monomials, c_mon)
            push!(r.coeffs, poly.coeffs[i])
        else
            r.coeffs[index] += poly.coeffs[i]
        end
    end
    return r
end

function Base.:(==)(c1::CyclicWord, c2::CyclicWord)
    (c1.exponents!=c2.exponents || c1.monoid!=c2.monoid) && return false
    nvars = length(c1.monoid.vertices)
    m1 = c1.reduced_word
    m2 = c2.reduced_word
    for i in 2:nvars
        if m1^i/m2
            return true
        end
    end
    return false
end

Base.:(==)(c1::CyclicWord, c2::GraphProductWord{Variable}) = c1==cyclic_reduce(c2)
Base.:(==)(c1::GraphProductWord{Variable}, c2::CyclicWord) = cyclic_reduce(c1)==c2
Base.:(==)(c1::CyclicWord, c2::PCMonomial) = c1==cyclic_reduce(c2)
Base.:(==)(c1::PCMonomial, c2::CyclicWord) = cyclic_reduce(c1)==c2

conjugate(m::CyclicWord) = cyclic_reduce(m.ref_word')
Base.conj(m::CyclicWord) = conjugate(m)
Base.adjoint(m::CyclicWord) = Base.conj(m)

function Base.hash(t::CyclicWord, h::UInt)
    return hash(t.exponents, hash(t.monoid, hash(0x23269960ff982ff6, h)))
end

Base.zero(p::Polynomial{C_T, CyclicWord}) where {C_T} =
    Polynomial(CyclicWord[], C_T[], p.monoid)

function add_poly(
    p::Polynomial{C_T_1, CyclicWord},
    q::Polynomial{C_T_2, CyclicWord},
) where {C_T_1, C_T_2}
    T = MA.promote_operation(+, C_T_1, C_T_2)
    r = Polynomial(CyclicWord[], T[], p.monoid)  # Create a new polynomial with the same monoid
    for (m, c) in p
        push!(r.monomials, m)
        push!(r.coeffs, c)
    end
    for (m, c) in q
        index = findfirst(x -> x == m, r.monomials)
        if index === nothing
            push!(r.monomials, m)
            push!(r.coeffs, c)
        else
            r.coeffs[index] += c
            if r.coeffs[index] == 0
                deleteat!(r.monomials, index)
                deleteat!(r.coeffs, index)
            end
        end
    end
    return r
end

function coefficient(p::Polynomial{C_T, CyclicWord}, m::CyclicWord) where {C_T}
    index = findfirst(x -> x == m, p.monomials)
    return index === nothing ? zero(C_T) : p.coeffs[index]
end

Base.:+(p::Polynomial{C1, CyclicWord}, q::Polynomial{C2, CyclicWord}) where {C1, C2} =
    add_poly(p, q)
Base.:+(
    p::Polynomial{C1, CyclicWord},
    q::Polynomial{C2, M},
) where {C1, C2, M <: AbstractMonomial} = add_poly(p, cyclic_reduce(q))
Base.:+(
    q::Polynomial{C1, M},
    p::Polynomial{C2, CyclicWord},
) where {C1, C2, M <: AbstractMonomial} = add_poly(p, cyclic_reduce(q))

Base.:-(p::Polynomial{C1, CyclicWord}, q::Polynomial{C2, CyclicWord}) where {C1, C2} =
    add_poly(p, -q)
Base.:-(
    p::Polynomial{C1, CyclicWord},
    q::Polynomial{C2, M},
) where {C1, C2, M <: AbstractMonomial} = add_poly(p, -cyclic_reduce(q))
Base.:-(
    q::Polynomial{C1, M},
    p::Polynomial{C2, CyclicWord},
) where {C1, C2, M <: AbstractMonomial} = add_poly(cyclic_reduce(q), -p)

# NCWord must not mix with CyclicWord
Base.:+(p::Polynomial{CT, CyclicWord}, n::NCWord) where {CT} =
    error("cannot mix NCWord with CyclicWord")
Base.:+(n::NCWord, p::Polynomial{CT, CyclicWord}) where {CT} =
    error("cannot mix NCWord with CyclicWord")
Base.:-(p::Polynomial{CT, CyclicWord}, n::NCWord) where {CT} =
    error("cannot mix NCWord with CyclicWord")
Base.:-(n::NCWord, p::Polynomial{CT, CyclicWord}) where {CT} =
    error("cannot mix NCWord with CyclicWord")

Base.:+(p::Polynomial{CT, CyclicWord}, x::AbstractMonomial) where {CT} =
    x isa CyclicWord ? add_poly(p, Polynomial(x)) :
    add_poly(p, Polynomial(cyclic_reduce(x)))
Base.:+(x::AbstractMonomial, p::Polynomial{CT, CyclicWord}) where {CT} =
    x isa CyclicWord ? add_poly(Polynomial(x), p) :
    add_poly(Polynomial(cyclic_reduce(x)), p)
Base.:-(p::Polynomial{CT, CyclicWord}, x::AbstractMonomial) where {CT} =
    x isa CyclicWord ? add_poly(p, -Polynomial(x)) :
    add_poly(p, -Polynomial(cyclic_reduce(x)))
Base.:-(x::AbstractMonomial, p::Polynomial{CT, CyclicWord}) where {CT} =
    x isa CyclicWord ? add_poly(Polynomial(x), -p) :
    add_poly(Polynomial(cyclic_reduce(x)), -p)
Base.:+(p::Polynomial{CT, CyclicWord}, x::Number) where {CT} =
    add_poly(p, Polynomial([cyclic_reduce(one(p.monoid))], [x], p.monoid))
Base.:+(x::Number, p::Polynomial{CT, CyclicWord}) where {CT} =
    add_poly(Polynomial([cyclic_reduce(one(p.monoid))], [x], p.monoid), p)
Base.:+(x::Number, m::CyclicWord) = x+Polynomial(m)
Base.:+(m::CyclicWord, x::Number) = x+Polynomial(m)
Base.:-(m::CyclicWord, x::Number) = (-x)+Polynomial(m)
Base.:-(x::Number, m::CyclicWord) = x+(-Polynomial(m))

# Assuming m is coming from cyclic_reduce of some monomial
Base.iszero(m::CyclicWord) = false
