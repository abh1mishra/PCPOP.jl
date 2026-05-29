
############################
##  OPTIMIZATION METHODS  ##
############################

function pcpop(p::Polynomial, k::Int;
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    normalize = true,
    tracial = false,
    solver = Mosek.Optimizer,
    optimize = false,
    reduce = false,
    block_diag = false,
    primal = true)

    basis = mons_at_level(p.monoid, k)
    return pcpop(p, basis, op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           solver = solver,
                           optimize = optimize,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
end

"""
    Polynomial optimization problem in partially commutative monoid `Μ`.
    
    pcpop(p, k)
    pcpop(p, basis)

    #Arguments:
    - `p` : polynomial `p` in monoid `Μ`.
    - `k` : level of the hierarchy for the SDP relaxation.
            (alternatively `basis` of relaxation subspace)

    Additional keyword arguments
    - `op_eq` : list of polynomial equality constraints.
    - `op_ge` : list of polynomial inequality constraints.
    - `tr_eq` : list of moment equality constraints.
    - `tr_ge` : list of moment inequality constraints.
    - `normalize` : Boolean value for normalization L(1) = 1.
    - `tracial` : Boolean value for tracial condition L(ab) = L(ba).
    - `solver` : SDP solver for optimization.
    - `optimize` : Boolean value for optimizing the model.
    - `reduce` : Boolean value for Jordan reduction.
    - `block_diag` : Boolean value for block diagonalization of the SDP.
    - `primal` : Boolean value for primal (moment) or dual (sum of squares) implementation.

    # Output:
    - `model` : JuMP model SDP relaxation of PCPOP.
    
    Primal moment relaxation Γ(u,v) = u*v:
    
    Max ∑ p(w) Γ(w)
    s.t. Γ(1) = 1
         Γ in PSD(k)
    
    Dual of sum of square decomposition:
    
    Min t
    s.t. t - p in SOS(k)    
"""
function pcpop(p::Polynomial, basis;
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    normalize = true,
    tracial = false,
    solver = Mosek.Optimizer,
    optimize = false,
    reduce = false,
    block_diag = false,
    primal = true)

    if reduce
        Γ, C, A, b = npa_canonical(p, basis, 
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
        model = jordan_reduce(C, A, b, complex=true, diagonalize=block_diag)
    elseif primal
        model = npa(p, basis, 
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
    else
        model = sos(p, basis, 
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
    end

    set_solver(model, solver)

    if optimize
        optimize!(model)
    end

    return model
end

"""
    Dual sum of squares implementation.
"""
function sos(p::Polynomial, basis; 
            op_eq = [], 
            op_ge = [], 
            tr_eq = [],
            tr_ge = [],
            normalize = true, 
            tracial = false)

    op_ge = union(op_ge, op_eq, -op_eq)
    cores_psd = union([Polynomial(one(p.monoid))], Polynomial.(op_ge))
    cores = union(monomials.(cores_psd)...)
        
    matrix_psd = Dict(s => [tracial_reduce(x'*s*y, tracial=tracial) for x in basis_psd, y in basis_psd] for s in cores)
        
    basis_constraints = StarAlgebras.Basis{UInt16}(
              unique(union([union(matrix_psd[s]) for s in cores]...)),
            )
    Γ = Dict(s=> [(basis_constraints[m]) for m in matrix_psd[s]] for s in cores)

    model = JuMP.Model()

    if normalize
        tr_eq = union([(one(p.monoid), 1)], tr_eq)
    end
    n = length(tr_eq)
    m = length(tr_ge)
    JuMP.@variable model t[1:n+m]
    JuMP.@objective model Min sum([tr_eq[i][2]*t[i] for i in 1:n]) +
                              sum([tr_ge[i][2]*t[i] for i in n+1:m])  - p
    @constraint model t[n+1:n+m] .>= 0
    N = length(basis_psd)

    P = Dict(s=> JuMP.@variable model [1:N, 1:N] in PSDCone() for (i,s) in enumerate(cores_psd))

    for (idx, b) in enumerate(basis_constraints)
        if typeof(b) <: AbstractMonomial
            c = coefficient(objective, b)
            JuMP.@constraint sos_model sum(cs*LinearAlgebra.dot(P[s], Γ[ms] .== idx) for s in cores_psd for (ms, cs) in s) == c
        elseif typeof(b) <: Polynomial
            for (b_coeff, b_monomial) in zip(b.coeffs, b.monomials)
                c = b_coeff*coefficient(objective, b_monomial)
                JuMP.@constraint model sum(cs*LinearAlgebra.dot(P[s], Γ[ms] .== idx) for s in cores_psd for (ms, cs) in s) == c
            end
        end
    end

    return model
end

"""
    Primal moment implementation.
"""
function moments(p::Polynomial, basis; 
            op_eq = [], 
            op_ge = [], 
            tr_eq = [],
            tr_ge = [],
            normalize = true, 
            tracial = false,
            solver = Mosek.Optimizer,
            optimize = false)

    cores_zero = Polynomial.(op_eq)
    cores_psd = union([Polynomial(one(p.monoid))], Polynomial.(op_ge))
    if isempty(cores_zero)
        cores = union(monomials.(cores_psd)...)
    else
        cores = union(union(monomials.(cores_psd)...), union(monomials.(cores_zero)...)) 
    end
    matrix_psd = Dict(s => [tracial_reduce(x'*s*y, tracial=tracial) for x in basis_psd, y in basis_psd] for s in cores)

    basis_constraints = StarAlgebras.Basis{UInt16}(
        unique(union([union(matrix_psd[s]) for s in cores]...)),
        )

    Γ = Dict(s=> [(basis_constraints[m]) for m in matrix_psd[s]] for s in cores)

    if normalize
        tr_eq = union([(Polynomial(one(p.monoid)), 1)], tr_eq)
    end

    model = JuMP.Model()
    JuMP.@variable model y[1:length(basis_constraints)]

    # Objective function
    p = state_projection(p, TM)
    if tracial
        p = cyclic_reduce(p)
        moments = [(cyclic_reduce(Polynomial(m[1])), m[2]) for m in moments]
    end
    objective = sum(c*y[basis_constraints[m]] for (m,c) in p)
    JuMP.@objective model Max objective

    # Moment constraints
    for (s,b) in tr_eq
        JuMP.@constraint model sum(c*y[basis_constraints[m]] for (m,c) in s) == b
    end

    for (s,b) in tr_ge
        JuMP.@constraint model sum(c*y[basis_constraints[m]] for (m,c) in s) >= b
    end

    # Constraints
    P = Dict(s=>[y[i] for i in Γ[s]] for s in cores)
    # Positivity Γ(y) ≥ 0   
    for s in cores_psd
        JuMP.@constraint model sum(c*P[m] for (m, c) in s) in PSDCone()
    end
    # Zero Γ(y) = 0
    for s in cores_zero
        JuMP.@constraint model sum(c*P[m] for (m, c) in s) .== 0
    end

    return model   
end