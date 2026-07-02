using RowEchelon: rref

"""
    macaulay_matrix(F, T)

    Compute Macaulay matrix of F with respect to T.

    # Input:
    - `F` : list of polynomials (fi)_i.
    - `T` : list of monomials (tj)_j.

    # Output:
    - `A` : Matrix with A_ij coefficient of term tj in fi.

    macaulay_matrix(F)

    Compute Macaulay matrix of F with respect to suppF.
"""
function macaulay_matrix(F, T)
    A = Matrix(vcat([[coefficient(f, t) for t in T]' for f in F]...))
end

function macaulay_matrix(F)
    T = reverse(sort(union([f.monomials for f in F]...)))
    A = macaulay_matrix(F, T)
    return A, T
end

"""
    echelon_reduce(F)
    
    Remove dependencies leading terms from reduced echelon form of Macaulay matrix.
    
    # Input:
    - `F` : list of polynomials.
    
    # Output:
    - `row_polynomials` : Row polynomials of the reduced echelon form of Mac(F, suppF).
"""
function echelon_reduce(F)
    A, T = macaulay_matrix(F)
    A_rref = rref(A)
    row_polynomials = unique([Ai'*T for Ai in eachrow(A_rref)])
    row_polynomials = filter(x -> !(x == 0), row_polynomials)
end

"""
    shifts(G, d)
    
    Construct all shifts of G up to degree d
   
    F = {l*g*r : (l,g,r) in M x G x M and deg(l*g*r) <= d}.
    
    # Input:
    - `G` : list of polynomials.
    
    # Output:
    - `F` : list shifts of G up to degree d.
    
    Monomials at PCPOP level `k` are introduced in `shifts(G, 2*k)`.
"""
function shifts(G, d)
    F = []
    for g in G
        deg_g = degree(g.monomials[1])
        for deg_l in 0:(d - deg_g)
            #   monomials_l = monomials(g.monoid, deg_l)
            monomials_l = mons_at_level(g, deg_l)
            for deg_r in 0:(d - deg_g - deg_l)
                #   monomials_r = monomials(g.monoid, deg_r)
                monomials_r = mons_at_level(g, deg_r)
                F = append!(F, [l*g*r for l in monomials_l for r in monomials_r])
            end
        end
    end
    return F
end

"""
    self_reduce(G)
    
    Reduce a Gröbner basis removing redundant polynomials.
    
    # Input:
    - `G` : Gröbner basis
    
    # Output:
    - Reduced Gröbner basis of G
"""

"""
    reduce(f, g)
    
    Reduce polynomial `f` with respect to `g`.

    # Input:
    -`f` : polynomial to be reduced
           f = sum_i a_i f_i
    -`g` : polynomial reduction
           g = sum_j b_j g_j
    Alternatively, `G` list of polynomials.

    If `f_i = l*g_1*r` replace `g_1` by `(g - g_1)`

    # Output:
    - Polynomial `f` after those replacements

    Output depends on they choice of factorization f_i = l*g_1*r.
    Grobner bases produce confluent replacements.
"""
function reduce_grobner(f::Polynomial, g::Polynomial)
    leading_monomial = g.monomials[1]
    leading_coefficient = g.coeffs[1]
    monomial_pointer=1
    while monomial_pointer <= length(f.monomials)
        m = f.monomials[monomial_pointer]
        ok, (l, r) = divide(m, leading_monomial)
        if ok
            f = f - (f.coeffs[monomial_pointer]/leading_coefficient)*l*g*r
        else
            monomial_pointer += 1
        end
    end
    # for (m,c) in f
    #     ok,factor = divide(m,leading_monomial)
    #     if ok
    #         f = f - (c/leading_coefficient)*factor[1]*g*factor[2]
    #     end
    # end
    return f
end

function reduce_grobner(f::Polynomial, G)
    for g in G
        f = reduce_grobner(f, g)
    end
    return f
end

function reduce_grobner(F::Vector, G)
    F_reduced = copy(F)
    return filter(x -> !(x==0), unique([reduce_grobner(f, G) for f in F_reduced]))
end

function self_reduce(G)
    G_reduced = []
    G_remains = copy(G)
    while !isempty(G_remains)
        g = pop!(G_remains)
        G_remains = reduce_grobner(G_remains, g)
        G_reduced = filter(x -> !(x==0), reduce_grobner(G_reduced, g))
        G_reduced = append!(G_reduced, [g])
    end
    return reverse(filter(x -> !(x==0), G_reduced))
end

"""
    macaulay_grobner(G, d)
    
    Truncated Gröbner basis of G up to degree d.
    
    # Input:
    - `G` : list of polynomials.
    - `d` : truncation degree.
    
    # Output:
    - `H` : reduced Gröbner basis of G truncated at degree d.
    
    Monomials at PCPOP level `k` are reduced in `macaulay_grobner(G, 2*k)`.
"""
function macaulay_grobner(G, d)
    F = shifts(G, d)
    H = echelon_reduce(F)
    H = clear(H)
    H = self_reduce(H)
end

# Debuggin
function clear(h::Polynomial)
    ind = findall(x->!iszero(x), h.coeffs)
    if !isempty(ind)
        return sum(h.coeffs[ind] .* h.monomials[ind])
    else
        return 0
    end
end

function clear(H)
    filter(x->!iszero(x), clear.(H))
end
