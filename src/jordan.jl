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
function jordan_reduce(C, A, b; verbose=false, complex=false)
    # Optimal invariant subspace
    P, blkD = block_diagonal(C, A, b, verbose=verbose, complex=complex)
    PMat = hcat([sparse(vec(P.P .== i)) for i = 1:P.n]...)
    newA = A * PMat
    newB = b
    newC = C' * PMat
    # Reduced model
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => !verbose))
    x = @variable(model, x[1:P.n])
    # Linear constraints
    for i in 1:size(newA, 1)
        constraint_i = AffExpr(0)
            for j in 1:P.n
                add_to_expression!(constraint_i, newA[i, j], x[j])
            end
        @constraint(model, constraint_i == newB[i])
    end
    # Objective
    obj = AffExpr(0)
    for i in 1:P.n
        add_to_expression!(obj, newC[i], x[i])
    end
    @objective(model, Max, obj)
    # PSD constraints
    psdBlocks = sum(blkD.blks[i] .* x[i] for i = 1:P.n)
    for blk in psdBlocks
        if size(blk, 1) > 1
            blk = realify(blk;complex=complex)
            @constraint(model, blk in PSDCone())
        else
            blk = realify(blk;complex=complex)
            @constraint(model, blk .>= 0)
        end
    end
    # Optimize
    optimize!(model)  
    if verbose
        @show termination_status(model)
        @show objective_value(model)
    end    
    return model, P, blkD
end

# Block diagonalization through optimal admissible subspace
function block_diagonal(C, A, b; verbose=false, complex=false)
    P = admPartSubspace(C, A, b, verbose)
    blkD = blockDiagonalize(P, complex=complex)
    return P, blkD
end

function realify(M::AbstractMatrix; complex=false)
    if !complex
        return M
    else
        realM = real(M)
        imagM = imag(M)
        return [realM -imagM; imagM realM]
    end
end

"""
    
    write_canonical(model::Model)

    Write a JuMP.Model in canonical form (C, A, b)

    sup ⟨C, X⟩
    s.t. AX = b
         X,X_i ≥ 0
X=@variable(model, X[1:2,1:2], PSD)
 --- IGNORE ---
    (!) For now only one PSD variable and EqualTo constraints.
    (!) Missing the constraints identifying equivalent moments.

    See MOI documentation p.333 for more information about dualization
    http://jump.dev/MathOptInterface.jl/stable/MathOptInterface.pdf

"""

"""
npa_dual exports model,D,C,A,X,B 

where npa_dual accepts the same arguments as npa, use rm=true flag to get te desired return values, otherwise it will continue to optmize.

model is the JuMP model, D is a dictionary of the unique monomials in the moment matrix and their positions,
C is the matrix for the objective function(so call @objective @objective(model,Max,dot(C,X)) )
X is the single largest PSD variable of form Blockdiagonal(Prinicpal_moment_matrix, localising_moment_matrix1, localising_moment_matrix2,...),
A is vector of the sparse matrices for the linear constraints, so every element in A is the size of X
B is the vector of the right hand side of the linear constraints, so length of B is the same as length of A.

In order to get vectorized form, call vec(A),Vec(X)...

Also, by default, the off diagonal terms are not set to zero. This shall include a lot of additional constraints, so use extra_zeros=true flag to set them to zero.
extra_zeros=true will include additional rows in A. I didn't see any significant change in SDP solution by including these zeros but the setup time increased dramatically. 
"""


function write_canonical(model::Model)
    # Get the model data
    vars = all_variables(model)                     
    obj = objective_function(model)                
    C = [JuMP.coefficient(obj, v) for v in vars]

    A = []
    b = Float64[]
    for cref in all_constraints(model, include_variable_in_set_constraints=false)
        # get function and set for constraint cref
        func = constraint_object(cref).func
        set  = constraint_object(cref).set
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

# # Example statepop 7.2.1
# @pcmonoid M a[2,0] b[2,0]
# Unipotent.(a)
# Unipotent.(b)
# @comms a b
# build(M)

# TM = make_trace_monoid(M, 6, tracial=false)
# p  = (state(a[1]*b[2], TM) + state(a[2]*b[1], TM))^2 
# p += (state(a[1]*b[1], TM) - state(a[2]*b[2], TM))^2
# basis = trace_monomials(TM, 0:3, tracial=false)
# model = tpop(p, TM, basis, tracial=false)

# # Jordan reduction
# C, A, b = write_canonical(model)
# model_red, P, blkD = jordan_reduce(C, A, b; verbose=true)

# # Compare optimal solutions
# set_optimizer(model, Mosek.Optimizer)
# set_silent!(model)
# optimize!(model)
# println(termination_status(model))
# println(objective_value(model))

# set_optimizer(model_red, Mosek.Optimizer)
# set_silent!(model_red)
# optimize!(model_red)
# println(termination_status(model_red))
# println(objective_value(model_red))