using LinearAlgebra
using SparseArrays
using MosekTools: Mosek
using JuMP

include("symmetries.jl")
include("grobner.jl")

# Compatibility issue dot() in julia1.11
dot_vars(C, iv::SparseVector{<:Number}) =
    sum(iv[i] * C[i] for i in iv.nzind)



# Remove unidentified duplicates
function reduce_duplicates(matrix::Matrix{T}) where T
    seen = T[]
    reps = T[]

    # Find all unique elements and build map from duplicates to rep
    for x in matrix
        found = findfirst(y -> x == y, seen)
        if isnothing(found)
            push!(seen, x)
            push!(reps, x)
        else
            push!(reps, seen[found])
        end
    end

    # Replace duplicates in the matrix with the first representative
    reduced = similar(matrix)
    i = 1
    for row in axes(matrix, 1), col in axes(matrix, 2)
        reduced[row, col] = reps[i]
        i += 1
    end

    return reduced
end


function variables(M::AbstractMonoid)
    return union(variables.(M.vertices)...)
end

function tracial_reduce(word::AbstractMonomial; tracial=false)
    if tracial
        return cyclic_reduce(word)
    else
        return word
    end
end

function tracial_reduce(p::Polynomial; tracial=false)
    sum([c*tracial_reduce(m, tracial=tracial) for (c, m) in zip(p.coeffs, p.monomials)])
end


"""
    Polynomial optimization problem in partially commutative monoid `Μ`.
    
    #Arguments:
    - `poly` : polynomial `p` in monoid `Μ` with degree d.
    - `k` : level of relaxation
            when not specified set to smallest size d÷2.
    - `equalities` : list with polynomial constraints.
    
    Additionally, for symmetry reduction:
    - `G` : Group of symmetries on monoid `Μ`.
    - `action` : of the group `G` on monomials in `Μ`.
    - `diagonalize` : block diagonal SDP with symmetries.
    
    Alternatively, wedderburn decomposition:
    - `wedderburn` : SymbolicWedderburn.WedderburnDecomposition.
    - `basis_psd`  : spanning subspace for SDP relaxation. 
    
    # Output:
    - JuMP model (unsolved) with SDP relaxation POP.
    
    Moment relaxation Γ(u,v) = u*v:
    
    Max sum_i p(i) Γ(i)
    s.t. Γ(1) = 1
         Γ in PSD(k)
    
    Dual of sum of square decomposition:
    
    Min t
    s.t. t - p in SOS(k)    
"""
function pcpop(poly::Polynomial, k::Int; equalities = [], inequalities=[], truncate = "degree", tracial=false)
    M = poly.monoid
    basis_psd = mons_at_level(M, k)
    cores_psd = union([Polynomial(one(M))], Polynomial.(inequalities))
    cores = union(monomials.(cores_psd)...)
    if isempty(equalities)  
        matrix_psd = Dict(s=>[tracial_reduce(x'*s*y, tracial=tracial) for x in basis_psd, y in basis_psd] for s in cores)
    else
        max_degree = maximum([degree(g) for g in equalities])
        if truncate == "degree"  
            truncate = max_degree
        elseif truncate < max_degree
            throw(ArgumentError("Truncation degree $(truncate) expected at least constraints degree $(max_degree)"))
        end
        grobner_truncated = macaulay_grobner(equalities, truncate)
        matrix_psd = Dict(s=>tracial_reduce.(reduce_duplicates([reduce_grobner(Polynomial(x'*s*y), grobner_truncated) for x in basis_psd, y in basis_psd], tracial=tracial)) for s in cores)
    end   
    basis_constraints = StarAlgebras.Basis{UInt16}(
              unique(union([collect(values(matrix_psd[s])) for s in cores]...)),
            )
    Γ = Dict(s=> [(basis_constraints[m]) for m in matrix_psd[s]] for s in cores)

    sos_model = JuMP.Model()
    JuMP.@variable sos_model t
    JuMP.@objective sos_model Min t
    n = length(basis_psd)
    P = Dict(s=> JuMP.@variable sos_model [1:n, 1:n] in PSDCone() for (i,s) in enumerate(cores_psd))

    #  base_name=Symbol("P$i")
    objective = t*one(M)-poly
    for (idx, b) in enumerate(basis_constraints)
        if typeof(b) <: AbstractMonomial
            c = coefficient(objective, b)
            JuMP.@constraint sos_model sum(cs*LinearAlgebra.dot(P[s], Γ[ms] .== idx) for s in cores_psd for (ms, cs) in s) == c
        elseif typeof(b) <: Polynomial
            for (b_coeff, b_monomial) in zip(b.coeffs, b.monomials)
                c = b_coeff*coefficient(objective, b_monomial)
                JuMP.@constraint sos_model sum(cs*LinearAlgebra.dot(P[s], Γ[ms] .== idx) for s in cores_psd for (ms, cs) in s) == c
            end
        end
    end

    return sos_model
end

function pcpop(poly::Polynomial; equalities=[], inequalities=[], truncate="degree", tracial=false)
    pcpop(poly, Int(ceil(degree(poly)/2)), equalities=equalities, inequalities=inequalities, truncate=truncate, tracial=tracial)
end

function pcpop(poly::Polynomial, k::Int, G::GroupsCore.Group, action::SymbolicWedderburn.Action; diagonalize=false)
    M = poly.monoid
    basis_psd = mons_at_level(M,k)
    basis_constraints = StarAlgebras.Basis{UInt16}(
              unique_array([x'*y for x in basis_psd, y in basis_psd]),
    )
    if diagonalize
        wedderburn = SymbolicWedderburn.WedderburnDecomposition(
            Float64,
            G,
            action,
            basis_constraints,
            basis_psd,
            semisimple=true,
        )
        
        return pcpop(poly, wedderburn, basis_psd)
    else
        
    tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
    invariant_vs = invariant_vectors(tbl, action, basis_constraints)
    
    Γ = [basis_constraints[x'*y] for x in basis_psd, y in basis_psd]

    sos_model = JuMP.Model()
    JuMP.@variable sos_model t
    JuMP.@objective sos_model Min t
    n = length(basis_psd)
    P = JuMP.@variable sos_model P[1:n, 1:n] Symmetric
    JuMP.@constraint sos_model P in PSDCone()

    # preallocating
    Γ_orb = similar(Γ, Float64)

    C = coefficient(t*one(M)-poly, basis_constraints)

    for iv in invariant_vs
        c = dot_vars(C, iv)
        # average Γs into Γ_orb with weights given by iv
        Γ_orb = invariant_constraint!(Γ_orb, Γ, iv)
        JuMP.@constraint sos_model dot(Γ_orb, P) == c
    end

    return sos_model
    end
end

## TODO : symmetries + equalities
##function pcpop(poly::Polynomial, k::Int, G::GroupsCore.Group, action::SymbolicWedderburn.Action; equalities=[], truncate="degree", diagonalize=false)
##    M = poly.monoid
##    basis_psd = mons_at_level(M, 1:k)
##    if isempty(equalities)  
##        matrix_psd = [x'*y for x in basis_psd, y in basis_psd]
##    else
##        max_degree = maximum([degree(g) for g in equalities])
##        if truncate == "degree"  
##            truncate = max_degree
##        elseif truncate < max_degree
##            throw(ArgumentError("Truncation degree $(truncate) expected at least constraints degree $(max_degree)"))
##        end
##        grobner_truncated = macaulay_grobner(equalities, truncate)
##        matrix_psd = [reduce_grobner(Polynomial(x'*y), grobner_truncated) for x in basis_psd, y in basis_psd]
##        matrix_psd = reduce_duplicates(matrix_psd)
##    end   
##    basis_constraints = StarAlgebras.Basis{UInt16}(
##              unique(matrix_psd),
##            )
##    if diagonalize
##        wedderburn = SymbolicWedderburn.WedderburnDecomposition(
##            Float64,
##            G,
##            action,
##            basis_constraints,
##            basis_psd,
##            semisimple=true,
##        )
##        
##        return pcpop(poly, wedderburn, basis_psd)
##    else
##        
##    tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
##    # ERROR when basis_constraints has polynomials!
##    # separate monomials and then linear combination?
##    invariant_vs = invariant_vectors(tbl, action, basis_constraints)
##    
##    Γ = [basis_constraints[m] for m in matrix_psd]
##
##    sos_model = JuMP.Model()
##    JuMP.@variable sos_model t
##    JuMP.@objective sos_model Min t
##    n = length(basis_psd)
##    P = JuMP.@variable sos_model P[1:n, 1:n] Symmetric
##    JuMP.@constraint sos_model P in PSDCone()
##
##    # preallocating
##    Γ_orb = similar(Γ, Float64)
##
##    C = coefficient(t*one(M)-poly, basis_constraints)
##    objective = t*one(M)-poly
##    for (idx, b) in enumerate(basis_constraints)
##        if typeof(b) <: AbstractMonomial
##            c = coefficient(objective, b)
##            JuMP.@constraint sos_model LinearAlgebra.dot(P, Γ .== idx) == c
##        elseif typeof(b) <: Polynomial
##            for (b_coeff, b_monomial) in zip(b.coeffs, b.monomials)
##                c = b_coeff*coefficient(objective, b_monomial)
##                JuMP.@constraint sos_model LinearAlgebra.dot(P, Γ .== idx) == c
##            end
##        end
##    end
##    for iv in invariant_vs
##        c = dot_vars(C, iv)
##        # average Γs into Γ_orb with weights given by iv
##        Γ_orb = invariant_constraint!(Γ_orb, Γ, iv)
##        JuMP.@constraint sos_model dot(Γ_orb, P) == c
##    end
##
##    return sos_model
##    end
##end

function invariant_constraint!(
    Γ_orb::AbstractMatrix{<:AbstractFloat},
    Γ::Matrix{<:Integer},
    invariant_vec::SparseVector,
)
    Γ_orb .= zero(eltype(Γ_orb))
    for i in eachindex(Γ)
        if Γ[i] ∈ SparseArrays.nonzeroinds(invariant_vec)
            Γ_orb[i] += invariant_vec[Γ[i]]
        end
    end
    return Γ_orb
end

function pcpop(
    poly::Polynomial,
    wedderburn::SymbolicWedderburn.WedderburnDecomposition,
    basis_psd;
)
    model = JuMP.Model()
    M = poly.monoid
    
    Γ = let basis_constraints = SymbolicWedderburn.basis(wedderburn)
        [basis_constraints[x'*y] for x in basis_psd, y in basis_psd]
    end

    JuMP.@variable model t
    JuMP.@objective model Min t
    psds = map(SymbolicWedderburn.direct_summands(wedderburn)) do ds
        dim = size(ds, 1)
        P = JuMP.@variable model [1:dim, 1:dim] Symmetric
        JuMP.@constraint model P in PSDCone()
        return P
    end

    # preallocating
    # Γπs = zeros.(eltype(wedderburn), size.(psds))
    Γ_orb = similar(Γ, eltype(wedderburn))

    C = coefficient(
        t*one(M)-poly,
        SymbolicWedderburn.basis(wedderburn),
    )
    for iv in invariant_vectors(wedderburn)
        c = dot_vars(C, iv)
        Γ_orb = invariant_constraint!(Γ_orb, Γ, iv)
        # Γπs = SymbolicWedderburn.diagonalize!(Γπs, Γ_orb, wedderburn)
        Γπs = SymbolicWedderburn.diagonalize(Γ_orb, wedderburn)

        JuMP.@constraint model sum(
            dot(Γπ, Pπ) for (Γπ, Pπ) in zip(Γπs, psds) if !iszero(Γπ)
        ) == c
    end
    return model
end

"""
    Trace polynomial optimization problem in partially commutative monoid `Μ`.
    
    #Arguments:
    - `poly` : polynomial `p` in monoid `Μ` with degree d.
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
function tpop(poly::Polynomial, k::Int; equalities = [], truncate = "degree", tracial=false)
    M = poly.monoid
    TM = make_trace_monoid(M, 2*k, tracial=tracial)
    basis_psd = trace_monomials(TM, k, tracial=tracial)
    if isempty(equalities)  
        matrix_psd = [state_projection(x'*y, TM) for x in basis_psd, y in basis_psd]
    else
        max_degree = maximum([degree(g) for g in equalities])
        if truncate == "degree"  
            truncate = max_degree
        elseif truncate < max_degree
            throw(ArgumentError("Truncation degree $(truncate) expected at least constraints degree $(max_degree)"))
        end
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

    objective = state_embedding(t*one(M)-poly, TM)
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

function tpop(poly::Polynomial, k::Int, t::Int; equalities = [], truncate = "degree", tracial=false)
    M = poly.monoid
    TM = make_trace_monoid(M, 2*k, tracial=tracial)
    basis_psd = trace_monomials(TM, k, t, tracial=tracial)
    if isempty(equalities)  
        matrix_psd = [state_projection(x'*y, TM) for x in basis_psd, y in basis_psd]
    else
        max_degree = maximum([degree(g) for g in equalities])
        if truncate == "degree"  
            truncate = max_degree
        elseif truncate < max_degree
            throw(ArgumentError("Truncation degree $(truncate) expected at least constraints degree $(max_degree)"))
        end
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

    objective = state_embedding(t*one(M)-poly, TM)
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