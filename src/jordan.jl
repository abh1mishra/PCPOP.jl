using SDPSymmetryReduction
using SparseArrays: sparse, spzeros, SparseMatrixCSC
using StatsBase: sample

"""

    jordan_reduce(C, A, b; verbose=false, complex=false)

    Numerical reduction of the semidefinite program in canonical form

    sup ⟨C, x⟩
    s.t. Ax = b
         Γ ≥ 0

    Where `x´ is the vector with all variables in the PSD variable `Γ´.

    We call the implementation in SDPSymmetryReducction [Brosch & de Klerk 2020]
    or the Jordan algebra reduction method [Permenter & Parrilo 2016]

"""
function jordan_reduce(
    C,
    A,
    b;
    verbose = false,
    complex = false,
    epsilon = Base.rtoldefault(Float64),
    diagonalize = true,
    solver = default_solver(),
    optimize = false,
)
    if diagonalize
        return jordan_reduce_diagonal(
            C,
            A,
            b;
            verbose = verbose,
            complex = complex,
            epsilon = epsilon,
        )
    else
        P = admissible_subspace(
            SDPSymmetryReduction.Partition{UInt64},
            C,
            A,
            b;
            verbose = verbose,
        )
        PMat = hcat([sparse(vec(P.matrix .== i)) for i in 1:P.nparts]...)
        newA = A * PMat
        newB = b
        newC = C' * PMat
        model = Model()
        x = @variable(model, x[1:P.nparts])
        @objective model Max newC*x
        @constraint model newA*x .== newB
        @constraint model sum((P.matrix .== i)*x[i] for i in 1:P.nparts) in PSDCone()
        # Optimize
        if optimize
            set_solver(model, solver)
            optimize!(model)
        end
        return model, P, nothing
    end
end

function jordan_reduce_diagonal(
    C,
    A,
    b;
    verbose = false,
    complex = false,
    epsilon = Base.rtoldefault(Float64),
    solver = default_solver(),
    optimize = false,
)
    # Optimal invariant subspace
    P, blkD =
        block_diagonal(C, A, b, verbose = verbose, complex = complex, epsilon = epsilon)
    PMat = hcat([sparse(vec(P.matrix .== i)) for i in 1:P.nparts]...)
    newA = A * PMat
    newB = b
    newC = C' * PMat
    # Reduced model
    model = Model()
    x = @variable(model, x[1:P.nparts])
    # Linear constraints
    for i in 1:size(newA, 1)
        constraint_i = AffExpr(0)
        for j in 1:P.nparts
            add_to_expression!(constraint_i, newA[i, j], x[j])
        end
        @constraint(model, constraint_i == newB[i])
    end
    # Objective
    obj = AffExpr(0)
    for i in 1:P.nparts
        add_to_expression!(obj, newC[i], x[i])
    end
    @objective(model, Max, obj)
    # PSD constraints
    psdBlocks = sum(blkD.blks[i] .* x[i] for i in 1:P.nparts)
    for blk in psdBlocks
        if size(blk, 1) > 1
            blk = realify(blk; complex = complex)
            @constraint(model, blk in PSDCone())
        else
            blk = realify(blk; complex = complex)
            @constraint(model, blk .>= 0)
        end
    end
    # Optimize
    if optimize
        set_optimizer(model, solver)
        optimize!(model)
    end

    return model, P, blkD
end

# Block diagonalization through optimal admissible subspace
function block_diagonal(
    C,
    A,
    b;
    verbose = false,
    complex = false,
    epsilon = Base.rtoldefault(Float64),
)
    P = admissible_subspace(
        SDPSymmetryReduction.Partition{UInt32},
        C,
        A,
        b;
        verbose = verbose,
    )
    blkD = blockDiagonalize(P, verbose, complex = complex, epsilon = epsilon)
    return P, blkD
end

function realify(M::AbstractMatrix; complex = false)
    if !complex
        return M
    else
        realM = real(M)
        imagM = imag(M)
        return [realM -imagM; imagM realM]
    end
end

function cyclic_npa_moments_block_dual!(
    list_monomials::Vector{M},
    A,
    B,
    tsize;
    cPoly = 1,
    unique_mons = [],
    unique_pos = [],
    offset = 0,
    extra_zeros = false,
) where {M <: AbstractMonomial}

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials

            # Initialize the row of the A matrix for constraints
            Ai=spzeros(tsize * tsize)

            if cPoly == 1
                m = list_monomials[i]' * list_monomials[j]
                m1, m2 = cyclic_reduce(m), cyclic_reduce(m')
                if m1==0 || m2==0
                    Ai[(i - 1) * tsize + j] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                m_i=findfirst(x->(x==m1 || x==m2), unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    Ai[(i - 1) * tsize + j] = 1.0
                    upi, upj = unique_pos[m_i]
                    Ai[(upi - 1) * tsize + upj] = -1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                else
                    push!(unique_mons, m1)
                    push!(unique_pos, (i, j))
                end
            else
                monomial_product = list_monomials[i]' * cPoly * list_monomials[j]
                if monomial_product == 0
                    Ai[(offset + i - 1) * tsize + offset + j] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                Ai[(offset + i - 1) * tsize + offset + j] = 1.0
                for (m, c) in monomial_product
                    m1, m2 = cyclic_reduce(m), cyclic_reduce(m')
                    mi = findfirst(x->(x==m1 || x==m2), unique_mons)
                    upi, upj = unique_pos[mi]
                    Ai[(upi - 1) * tsize + upj] = -c
                end
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                Ai=spzeros(tsize * tsize)
                Ai[(offset + i - 1) * tsize + j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end

        for i in 1:num_monomials
            for j in (offset + num_monomials + 1):tsize
                Ai=spzeros(tsize * tsize)
                Ai[(offset + i - 1) * tsize + j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if cPoly == 1
        return unique_mons, unique_pos
    end
end

function npa_moments_block_dual!(
    list_monomials::Vector{M},
    A,
    B,
    tsize;
    cPoly = 1,
    unique_mons = [],
    unique_pos = [],
    offset = 0,
    extra_zeros = false,
) where {M <: AbstractMonomial}

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials

            # Initialize the row of the A matrix for constraints
            Ai=spzeros(tsize*tsize)

            if cPoly == 1
                m = real_rep(list_monomials[i]' * list_monomials[j])
                if m == 0
                    Ai[(i - 1) * tsize + j] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                m_i=findfirst(x->(x==m), unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    Ai[(i - 1) * tsize + j] = 1.0
                    upi, upj = unique_pos[m_i]
                    Ai[(upi - 1) * tsize + upj] = -1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                else
                    push!(unique_mons, m)
                    push!(unique_pos, (i, j))
                end
            else
                monomial_product =
                    real_rep(Polynomial(list_monomials[i]' * cPoly * list_monomials[j]))
                # assume the monomials PM exists and m and m' are in pm.
                if monomial_product == 0
                    Ai[(offset + i - 1) * tsize + (offset + j)] = 1.0
                    push!(A, Ai)
                    push!(B, 0.0)
                    continue
                end
                Ai[(offset + i - 1) * tsize + (offset + j)] = 1.0
                for (m, c) in monomial_product
                    m_i=findfirst(x->(x==m), unique_mons)
                    upi, upj = unique_pos[m_i]
                    Ai[(upi - 1) * tsize + upj] = -c
                end
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                Ai=spzeros(tsize*tsize)
                Ai[(offset + i - 1) * tsize + j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end

        for i in 1:num_monomials
            for j in (offset + num_monomials + 1):tsize
                Ai=spzeros(tsize*tsize)
                Ai[(offset + i - 1) * tsize + j] = 1.0
                push!(A, Ai)
                push!(B, 0.0)
            end
        end
    end

    if cPoly == 1
        return unique_mons, unique_pos
    end
end

"""
npa_canonical exports model,D,C,A,X,B 

where npa_canonical accepts the same arguments as npa, use rm=true flag to get the desired return values, otherwise it will continue to optimize.

model is the JuMP model, D is a dictionary of the unique monomials in the moment matrix and their positions,
C is the matrix for the objective function(so call @objective @objective(model,Max,dot(C,X)) )
X is the single largest PSD variable of form Blockdiagonal(Prinicpal_moment_matrix, localising_moment_matrix1, localising_moment_matrix2,...),
A is vector of the sparse matrices for the linear constraints, so every element in A is the size of X
B is the vector of the right hand side of the linear constraints, so length of B is the same as length of A.

In order to get vectorized form, call vec(A),Vec(X)...

Also, by default, the off diagonal terms are not set to zero. This shall include a lot of additional constraints, so use extra_zeros=true flag to set them to zero.
extra_zeros=true will include additional rows in A. I didn't see any significant change in SDP solution by including these zeros but the setup time increased dramatically. 
"""

function npa_canonical(
    obj,
    ops,
    ops_principal;
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    tracial = false,
    normalize = true,
    extra_zeros = false,
)
    println("Number of operators in the principal moment matrix: ", length(ops_principal))
    model=Model()
    A=Vector{SparseMatrixCSC{Float64, Int64}}([])
    B=Float64[]
    println("op_ge", length(op_ge))
    tsize=sum(
        [
            [length(ops_principal)];
            [length(ops) for i in 1:length(op_ge)];
            [2*length(ops) for i in 1:length(op_eq)];
            [1 for i in length(tr_ge)]
        ],
    )

    println("Size of the PSD variable: ", tsize, "x", tsize)
    X = @variable(model, [1:tsize, 1:tsize], PSD)
    unique_mons, unique_pos =
        tracial ?
        cyclic_npa_moments_block_dual!(
            ops_principal,
            A,
            B,
            tsize;
            extra_zeros = extra_zeros,
        ) : npa_moments_block_dual!(ops_principal, A, B, tsize; extra_zeros = extra_zeros)
    offset = length(ops_principal)
    println("Done building PM")

    for i in 1:length(op_ge)
        if tracial
            cyclic_npa_moments_block_dual!(
                ops,
                A,
                B,
                tsize;
                cPoly = op_ge[i],
                unique_mons = unique_mons,
                unique_pos = unique_pos,
                offset = offset,
                extra_zeros = extra_zeros,
            )
        else
            npa_moments_block_dual!(
                ops,
                A,
                B,
                tsize;
                cPoly = op_ge[i],
                unique_mons = unique_mons,
                unique_pos = unique_pos,
                offset = offset,
                extra_zeros = extra_zeros,
            )
        end
        offset += length(ops)
    end

    println("Done building LMI")

    for i in 1:length(op_eq)
        if tracial
            cyclic_npa_moments_block_dual!(
                ops,
                A,
                B,
                tsize;
                cPoly = op_eq[i],
                unique_mons = unique_mons,
                unique_pos = unique_pos,
                offset = offset,
                extra_zeros = extra_zeros,
            )
            offset += length(ops)
            cyclic_npa_moments_block_dual!(
                ops,
                A,
                B,
                tsize;
                cPoly = -op_eq[i],
                unique_mons = unique_mons,
                unique_pos = unique_pos,
                offset = offset,
                extra_zeros = extra_zeros,
            )
            offset += length(ops)

        else
            npa_moments_block_dual!(
                ops,
                A,
                B,
                tsize;
                cPoly = op_eq[i],
                unique_mons = unique_mons,
                unique_pos = unique_pos,
                offset = offset,
                extra_zeros = extra_zeros,
            )
            offset += length(ops)
            npa_moments_block_dual!(
                ops,
                A,
                B,
                tsize;
                cPoly = -op_eq[i],
                unique_mons = unique_mons,
                unique_pos = unique_pos,
                offset = offset,
                extra_zeros = extra_zeros,
            )
            offset += length(ops)
        end
    end

    # Add the constraints for the principal moment matrix

    for i in 1:length(tr_eq)
        Ai=spzeros(tsize * tsize)
        if tracial
            for (m, c) in Polynomial(tr_eq[i][1])
                m1, m2 = (cyclic_reduce(m), cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2), unique_mons)
                upi, upj = unique_pos[m_i]
                Ai[(upi - 1) * tsize + upj] = c
            end

        else
            tr_eq_poly=real_rep(Polynomial(Polynomial(tr_eq[i][1])))
            for (m, c) in tr_eq_poly
                mi = findfirst(x->x==m, unique_mons)
                upi, upj = unique_pos[mi]
                Ai[(upi - 1) * tsize + upj] = c
            end
        end
        push!(A, Ai)
        push!(B, tr_eq[i][2])
    end

    for i in 1:length(tr_ge)
        Ai=spzeros(tsize * tsize)
        Ai[offset * tsize + offset + 1] = -1.0

        if tracial
            for (m, c) in Polynomial(tr_ge[i][1])
                m1, m2 = (cyclic_reduce(m), cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2), unique_mons)
                upi, upj = unique_pos[m_i]
                Ai[(upi - 1) * tsize + upj] = c
            end
        else
            tr_ge_poly=real_rep(Polynomial(tr_ge[i][1]))
            for (m, c) in tr_ge_poly
                m_i = findfirst(x->x==m, unique_mons)
                upi, upj = unique_pos[m_i]
                Ai[(upi - 1) * tsize + upj] = c
            end
        end
        offset += 1
        push!(A, Ai)
        push!(B, tr_ge[i][2])
    end

    println("Done building trace constraints")

    if normalize
        id_elem=one(first(ops_principal))
        Ai=spzeros(tsize * tsize)
        if tracial
            id_elem=cyclic_reduce(id_elem)
            id_i=findfirst(x->x==id_elem, unique_mons)
        else
            id_i = findfirst(x->x==id_elem, unique_mons)
        end
        upi, upj = unique_pos[id_i]
        Ai[(upi - 1) * tsize + upj] = 1.0
        push!(A, Ai)
        push!(B, 1.0)
    end
    println("Done building normalization constraint")
    C = spzeros(tsize * tsize)
    if !is_number(obj)
        if tracial
            for (m, c) in Polynomial(obj)
                m1, m2 = (cyclic_reduce(m), cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2), unique_mons)
                if m_i==nothing
                    throw(ArgumentError("level not enough"))
                end
                upi, upj = unique_pos[m_i]
                C[(upi - 1) * tsize + upj] = c
            end
        else
            obj_poly=real_rep(Polynomial(Polynomial(obj)))
            for (m, c) in obj_poly
                mi = findfirst(x->x==m, unique_mons)
                if mi==nothing
                    throw(ArgumentError("level not enough"))
                end
                upi, upj = unique_pos[mi]
                C[(upi - 1) * tsize + upj] = c
            end
        end
    end
    println("Done building objective")
    A = vcat([r' for r in A]...)
    return Dict(zip(unique_mons, unique_pos)), C, A, B
end

"""
    
    write_canonical(model::Model)

    Write a JuMP.Model in canonical form (C, A, b)

    sup ⟨C, X⟩
    s.t. AX = b
         X ≥ 0

    (!) For now only one PSD variable and EqualTo constraints.
    (!) Missing the constraints identifying equivalent moments.

    See MOI documentation p.333 for more information about dualization
    http://jump.dev/MathOptInterface.jl/stable/MathOptInterface.pdf

"""

function write_canonical(model::Model)
    # Get the model data
    vars = all_variables(model)
    obj = objective_function(model)
    C = [JuMP.coefficient(obj, v) for v in vars]

    A = []
    b = Float64[]
    for cref in all_constraints(model, include_variable_in_set_constraints = false)
        # get function and set for constraint cref
        func = constraint_object(cref).func
        set = constraint_object(cref).set
        if func isa JuMP.GenericAffExpr && set isa MOI.EqualTo
            row = [JuMP.coefficient(func, v) for v in vars]
            push!(A, row)
            push!(b, set.value)
        end
    end

    return C, hcat(A...), b
end

function vec_sparse(M::SparseMatrixCSC)
    I, J, V = findnz(M)

    # Convert 2D indices (I, J) to 1D linear indices
    linear_indices = I .+ (J .- 1) .* size(M, 1)

    # Build the true SparseVector
    v = sparsevec(linear_indices, V, length(M))
    return v
end
