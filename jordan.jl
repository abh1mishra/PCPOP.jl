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
            @constraint(model, blk in PSDCone())
        else
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
    blkD = blockDiagonalize(P, complex)
    return P, blkD
end


"""
    
    write_canonical(model::Model)

    Write a JuMP.Model in canonical form (C, A, b)

    sup ⟨C, x⟩
    s.t. Ax = b
         Γ ≥ 0

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

# Example statepop 7.2.1
@pcmonoid M a[2,0] b[2,0]
Unipotent.(a)
Unipotent.(b)
@comms a b
build(M)

TM = make_trace_monoid(M, 6, tracial=false)
p  = (state(a[1]*b[2], TM) + state(a[2]*b[1], TM))^2 
p += (state(a[1]*b[1], TM) - state(a[2]*b[2], TM))^2
basis = trace_monomials(TM, 0:3, tracial=false)
model = tpop(p, TM, basis, tracial=false)

# Jordan reduction
C, A, b = write_canonical(model)
model_red, P, blkD = jordan_reduce(C, A, b; verbose=true)

# Compare optimal solutions
set_optimizer(model, Mosek.Optimizer)
set_silent!(model)
optimize!(model)
println(termination_status(model))
println(objective_value(model))

set_optimizer(model_red, Mosek.Optimizer)
set_silent!(model_red)
optimize!(model_red)
println(termination_status(model_red))
println(objective_value(model_red))