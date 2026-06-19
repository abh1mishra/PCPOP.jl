
############################
##  OPTIMIZATION METHODS  ##
############################

function pcpop!(p, k;
    min=false,
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    lvl_lm=-1,
    list_vars=[],
    normalize = true,
    tracial = false,
    solver = Mosek.Optimizer,
    optimize = true,
    reduce = false,
    block_diag = false,
    model_flags=[],
    primal = false,
    canonical=true,
    extra_zeros=false,
    silent = true)

    basis,basis_principal = basis_gen(p,k,op_eq,op_ge,tr_eq,tr_ge,list_vars,lvl_lm)
    return pcpop!(p, basis,basis_principal;
                           min = min,
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           solver = solver,
                           optimize = optimize,
                           reduce = reduce,
                           block_diag = block_diag,
                           model_flags = [],
                           primal = primal,
                           canonical = canonical,
                           extra_zeros = extra_zeros,
                           silent=silent)
end

"""
    Polynomial optimization problem in partially commutative monoid `Μ`.
    
    pcpop!(p, k)
    pcpop!(p, basis)

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
function pcpop!(p, basis, basis_principal;
    min = false,
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    normalize = true,
    tracial = false,
    solver = Mosek.Optimizer,
    optimize = true,
    reduce = false,
    block_diag = false,
    model_flags=[],
    primal = false,
    canonical=true,
    extra_zeros=false,
    silent = true)

    if is_number(p)
        @warn "The objective function is a constant, it is a feasibility check"
    end
    if reduce
        Γ, C, A, b = npa_canonical(p, basis, basis_principal;
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           extra_zeros = extra_zeros)
        model,_,_ = jordan_reduce(C, A, b, complex=true, diagonalize=block_diag)
    elseif primal && canonical
        model,Γ,X = npa(p, basis, basis_principal;
            min = min,
            op_eq = op_eq,
            op_ge = op_ge,
            tr_eq = tr_eq,
            tr_ge = tr_ge,
            normalize = normalize,
            tracial = tracial,
            extra_zeros = extra_zeros)
    elseif primal && !canonical
         model,Γ,PM = npa_nc(p, basis, basis_principal;
            min = min,
            op_eq = op_eq,
            op_ge = op_ge,
            tr_eq = tr_eq,
            tr_ge = tr_ge,
            normalize = normalize,
            tracial = tracial)
    else
        model,Γ,Zs = npa_dual(p, basis, basis_principal;
            min = min,
            op_eq = op_eq,
            op_ge = op_ge,
            tr_eq = tr_eq,
            tr_ge = tr_ge,
            normalize = normalize,
            tracial = tracial)
        # model = sos(p, basis;
        #                     min = min,
        #                    op_eq = op_eq,
        #                    op_ge = op_ge,
        #                    tr_eq = tr_eq,
        #                    tr_ge = tr_ge,
        #                    normalize = normalize,
        #                    tracial = tracial,
        #                    block_diag = block_diag)
    end

    if optimize
        set_optimizer(model, solver)
        silent ? set_silent(model) : Nothing
        for (flag,val) in model_flags
            set_optimizer_attribute(model,flag,val)
        end
        optimize!(model)
        return objective_value(model), model,basis, basis_principal
    end
    return model, basis, basis_principal
end

"""
    Dual sum of squares implementation.
"""
function sos(p::Polynomial, basis; 
            min=false,
            op_eq = [], 
            op_ge = [], 
            tr_eq = [],
            tr_ge = [],
            normalize = true, 
            tracial = false,
            block_diag = false)

    op_ge = union(op_ge, op_eq, -op_eq)
    cores_psd = union([Polynomial(one(p.monoid))], Polynomial.(op_ge))
    cores = union(monomials.(cores_psd)...)
        
    matrix_psd = Dict(s => [tracial_reduce(x'*s*y, tracial=tracial) for x in basis, y in basis] for s in cores)
        
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
    if m > 0
        JuMP.@objective model Min sum([tr_eq[i][2]*t[i] for i in 1:n]) + sum([tr_ge[i][2]*t[i] for i in n+1:m])
        @constraint model t[n+1:n+m] .>= 0
    else
        JuMP.@objective model Min sum([tr_eq[i][2]*t[i] for i in 1:n])
    end
    N = length(basis_psd)

    P = Dict(s=> JuMP.@variable model [1:N, 1:N] in PSDCone() for (i,s) in enumerate(cores_psd))
    # Objective function sum of squares decomposition
    objective = sum([tr_eq[i][1]*t[i] for i in 1:n]) + (-1)*p
    if m > 0
        objective += sum([tr_ge[i][1]*t[i] for i in n+1:m])
    end
    for (idx, b) in enumerate(basis_constraints)
        if typeof(b) <: AbstractMonomial
            c = coefficient(objective, b)
            JuMP.@constraint model sum(cs*LinearAlgebra.dot(P[s], Γ[ms] .== idx) for s in cores_psd for (ms, cs) in s) == c
        elseif typeof(b) <: Polynomial
            for (b_coeff, b_monomial) in zip(b.coeffs, b.monomials)
                c = b_coeff*coefficient(objective, b_monomial)
                JuMP.@constraint model sum(cs*LinearAlgebra.dot(P[s], Γ[ms] .== idx) for s in cores_psd for (ms, cs) in s) == c
            end
        end
    end

    return model
end

function sos_block_diag(p::Polynomial, basis; 
    op_eq = [], 
    op_ge = [], 
    tr_eq = [],
    tr_ge = [],
    normalize = true, 
    tracial = false,
    verbose = false,
    epsilon = 1e-8,
    complex = false)

op_ge = union(op_ge, op_eq, -op_eq)
cores_psd = union([Polynomial(one(p.monoid))], Polynomial.(op_ge))
cores = union(monomials.(cores_psd)...)

matrix_psd = Dict(s => [tracial_reduce(x'*s*y, tracial=tracial) for x in basis, y in basis] for s in cores)

basis_constraints = StarAlgebras.Basis{UInt16}(
      unique(union([union(matrix_psd[s]) for s in cores]...)),
    )
Γ = Dict(s=> [(basis_constraints[m]) for m in matrix_psd[s]] for s in cores)
basis_psd = Dict(s => StarAlgebras.Basis{UInt16}(
    union(sum(c*matrix_psd[m] for (m,c) in s))) for s in cores_psd)
Γ_psd = Dict(s => [basis_psd[s][m] for m in basis_psd[s]] for s in cores_psd)

model = JuMP.Model()

if normalize
tr_eq = union([(one(p.monoid), 1)], tr_eq)
end

n = length(tr_eq)
m = length(tr_ge)
JuMP.@variable model t[1:n+m]
if m > 0
    JuMP.@objective model Min sum([tr_eq[i][2]*t[i] for i in 1:n]) + sum([tr_ge[i][2]*t[i] for i in n+1:m])
    @constraint model t[n+1:n+m] .>= 0
else
    JuMP.@objective model Min sum([tr_eq[i][2]*t[i] for i in 1:n])
end
N = length(basis)

# Objective function sum of squares decomposition
objective = sum([tr_eq[i][1]*t[i] for i in 1:n]) + (-1)*p
if m > 0
    objective += sum([tr_ge[i][1]*t[i] for i in n+1:m])
end

blocks = Dict()
parts = Dict()

for s in cores_psd

    Ps = SDPSymmetryReduction.Partition(Γs)
    blkD =
        SDPSymmetryReduction.blockDiagonalize(
            Ps,
            verbose,
            epsilon=epsilon,
            complex=complex
            )

    if blkD === nothing && !complex
        blkD =
            SDPSymmetryReduction.blockDiagonalize(
                Ps,
                verbose,
                epsilon=epsilon,
                complex=true
            )
    end

    @assert blkD !== nothing

    blocks[s] = blkD
    parts[s] = Dict(i => (Ps .== i)  for i in 1:Ps.nparts)
end

P = Dict(s => [JuMP.@variable(model, [1:i, 1:i] in PSDCone()) for i in blocks[s].blkSizes] for s in cores_psd)

    for (k,b) in enumerate(basis_constraints)

        c = coefficient(objective,b)

        lhs = AffExpr()

        for s in cores_psd
            for m in basis_psd[s]
                cm = coefficient(m, b)
                i = basis_psd[s][m]
                lhs += cm*sum(blocks[s].blks[i])
            end
        end

        @constraint(model, lhs == c)
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
    matrix_psd = Dict(s => [tracial_reduce(x'*s*y, tracial=tracial) for x in basis, y in basis] for s in cores)

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

function moments_block_diag(
    p::Polynomial, basis; 
    op_eq = [], 
    op_ge = [], 
    tr_eq = [],
    tr_ge = [],
    normalize = true, 
    tracial = false,
    verbose = false,
    epsilon = 1e-8,
    complex = false)

    cores_zero = Polynomial.(op_eq)
    cores_psd = union([Polynomial(one(p.monoid))], Polynomial.(op_ge))
    if isempty(cores_zero)
        cores = union(monomials.(cores_psd)...)
    else
        cores = union(union(monomials.(cores_psd)...), union(monomials.(cores_zero)...)) 
    end

    matrix_psd = Dict(s => [tracial_reduce(x'*s*y, tracial=tracial) for x in basis, y in basis] for s in cores)

    basis_constraints = StarAlgebras.Basis{UInt16}(
      unique(union([union(matrix_psd[s]) for s in cores]...)),
    )
    Γ = Dict(s=> [(basis_constraints[m]) for m in matrix_psd[s]] for s in cores)

    if normalize
        tr_eq = union([(one(p.monoid), 1)], tr_eq)
    end

    model = JuMP.Model()
    JuMP.@variable model y[1:length(basis_constraints)]

    # Objective function
    JuMP.@objective model Max sum([c*y[basis_constraints[m]] for (m,c) in p])

    # Moment constraints
    for (s,b) in tr_eq
        JuMP.@constraint model sum(c*y[basis_constraints[m]] for (m,c) in s) == b
    end
    for (s,b) in tr_ge
        JuMP.@constraint model sum(c*y[basis_constraints[m]] for (m,c) in s) >= b
    end

    # PSD constraints
    for s in cores_psd
        basis_psd = StarAlgebras.Basis{UInt16}(
            union(sum(c*matrix_psd[m] for (m,c) in s)))
        Γ_psd = [basis_psd[s][m] for m in basis_psd[s]]
        P_psd = SDPSymmetryReduction.Partition(Γ_psd)
        blkD = SDPSymmetryReduction.blockDiagonalize(P_psd, verbose, epsilon=epsilon, complex=complex)

        psdBlocks = sum(blkD.blks[i] .* sum(ci*y[basis_constraints[mi]]) for (mi, ci) in basis_psd[i] for i = 1:P.nparts)
        for blk in psdBlocks
            if size(blk, 1) > 1
                blk = realify(blk;complex=complex)
                @constraint(model, blk in PSDCone())
            else
                blk = realify(blk;complex=complex)
                @constraint(model, blk .>= 0)
            end
        end
    end

    # Zero Γ(y) = 0
    for s in cores_zero
        JuMP.@constraint model sum(c*y[x'*m*y] for (m, c) in s for x in basis for y in basis) .== 0
    end

    return model
end