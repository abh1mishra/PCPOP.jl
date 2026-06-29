
"""
    Polynomial optimization problem in partially commutative monoid `Μ`.
    
    #Arguments:
    - `poly` : polynomial `p` in monoid `Μ` with degree d.
    - `k` : level of relaxation
            when not specified set to smallest size d÷2.
        Alternatively, accept basis for PSD relaxation

    Additional keyword arguments
    - `equalities` : list with polynomial constraints.
    - `equalities` : list of operator equality constraints.
    - `inequalities` : list of operator inequality constraints.
    - `moments` : list of moment equality constraints.
    - `normalize` : Boolean value for normalization L(1) = 1.
    - `tracial` : Boolean value for tracial condition L(ab) = L(ba).
    - `localize` : Boolean value to implement equalities
        true: pairs of inequalities r = 0  adds  r ≥ 0  &&  -r ≥ 0
        false: substitution rules r1 - r2 = 0  adds  r1 -> r2
    - `truncate` : Int degree truncation Gröbner bases for substitution rules.
    
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
function pcpop(p::Polynomial, k::Int; equalities=[], inequalities=[], moments=[], normalize=true, tracial=false, localize=false, truncate="degree")
    basis = mons_at_level(p.monoid, k)
    return pcpop(p, basis, equalities=equalities, 
                            inequalities=inequalities, 
                            moments=moments, 
                            normalize=normalize, 
                            tracial=tracial, 
                            localize=localize,
                            truncate=truncate)
end

function pcpop(p::Polynomial, basis_psd; equalities=[], inequalities=[], moments=[], normalize=true, tracial=false, localize=false, truncate="degree")
    if localize
        inequalities = union(inequalities, equalities, (-1).*equalities)
        equalities = []
    end
    println("Number of operators in the principal moment matrix: ", length(basis_psd))
    cores_psd = union([Polynomial(one(p.monoid))], Polynomial.(inequalities))
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
        matrix_psd = Dict(s=>tracial_reduce.(reduce_duplicates([reduce_grobner(Polynomial(x'*s*y), grobner_truncated) for x in basis_psd, y in basis_psd]), tracial=tracial) for s in cores)
    end   
    basis_constraints = StarAlgebras.Basis{UInt64}(
              unique(union([union(matrix_psd[s]) for s in cores]...)),
            )
    Γ = Dict(s=> [(basis_constraints[m]) for m in matrix_psd[s]] for s in cores)

    sos_model = JuMP.Model()

    if normalize
        moments = union([(one(p.monoid), 1)], moments)
    end
    nm = length(moments)
    JuMP.@variable sos_model t[1:nm]
    JuMP.@objective sos_model Min sum([moments[i][2]*t[i] for i in 1:nm])
    n = length(basis_psd)
    P = Dict(s=> JuMP.@variable sos_model [1:n, 1:n] in PSDCone() for (i,s) in enumerate(cores_psd))


    objective = sum([moments[i][1]*t[i] for i in 1:nm])-p
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