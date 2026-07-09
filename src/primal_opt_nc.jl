function cyclic_npa_moments_block_nc(
    list_monomials::Vector{M},
    model;
    cPoly = 1,
    mons_D = Dict{AbstractMonomial,AbstractJuMPScalar}([]),
    eq = false,
    progress = false,
) where {M <: AbstractMonomial}

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables
    moments_matrix = Matrix{JuMP.AffExpr}(undef, num_monomials, num_monomials)

    if progress
        prog_obj = Progress(
            num_monomials^2,
            desc = "Creation of moment matrix for polynomial $(cPoly)";
            showspeed = true,
            output = stdout,
            dt = 0.0,
        )
    end

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in 1:num_monomials
            if progress
                next!(prog_obj)
            end

            # Compute the product of the monomials
            monomial_product = list_monomials[i]' * cPoly * list_monomials[j]

            moments_matrix[i, j] = 0.0
            # Check if the product already exists in the dictionary
            for (m, c) in monomial_product
                m1, m2 = (cyclic_reduce(m), cyclic_reduce(m'))
                if m1==0 || m2==0
                    continue
                end
                if haskey(mons_D, m1) || haskey(mons_D, m2)
                    # Use the existing JuMP variable
                    var = haskey(mons_D, m1) ? mons_D[m1] : mons_D[m2]
                    # Use the existing JuMP variable
                    moments_matrix[i, j] += c*var
                else
                    # Create a new JuMP variable
                    new_var = @variable(model)
                    # Store the new variable in the dictionary
                    mons_D[m1] = new_var
                    # Use the new variable in the matrix
                    moments_matrix[i, j] += c*new_var
                end
            end
        end
    end

    if !eq
        @constraint(model, moments_matrix >= 0, PSDCone())
    else
        @constraint(model, moments_matrix .== 0)
    end

    return moments_matrix, mons_D
end

function npa_moments_block_nc(
    list_monomials::Vector{M},
    model;
    cPoly = 1,
    mons_D = Dict{AbstractMonomial,AbstractJuMPScalar}([]),
    eq = false,
    progress = false,
) where {M <: AbstractMonomial}

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables
    moments_matrix = Matrix{JuMP.AffExpr}(undef, num_monomials, num_monomials)

    if progress
        prog_obj = Progress(
            num_monomials*(num_monomials+1)÷2,
            desc = "Creation of moment matrix for polynomial $(cPoly)";
            showspeed = true,
            output = stdout,
            dt = 0.0,
        )
    end

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials
            if progress
                next!(prog_obj)
            end
            # Compute the product of the monomials
            monomial_product =
                real_rep(Polynomial(list_monomials[i]' * cPoly * list_monomials[j]))
            # Initialize the JuMP variable for the matrix entry
            moments_matrix[i, j] = 0.0
            for (m, c) in monomial_product
                # Check if the product already exists in the dictionary

                if haskey(mons_D, m)
                    # Use the existing JuMP variable
                    moments_matrix[i, j] += c*mons_D[m]
                else
                    # Create a new JuMP variable
                    new_var = @variable(model)
                    # Store the new variable in the dictionary
                    mons_D[m] = new_var
                    # Use the new variable in the matrix
                    moments_matrix[i, j] += c*new_var
                end
            end
            if i != j
                moments_matrix[j, i] = moments_matrix[i, j]
            end
        end
    end
    if !eq
        @constraint(model, moments_matrix in PSDCone())
    else
        @constraint(model, moments_matrix .== 0)
    end
    return moments_matrix, mons_D
end

function npa_nc(
    obj,
    ops,
    ops_principal;
    min = true,
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    tracial = false,
    normalize = true,
    progress = false,
)
    model=Model()

    if eltype(ops_principal) <: AbstractMonomial
        principal_moments_matrix, mons_D =
            tracial ?
            cyclic_npa_moments_block_nc(ops_principal, model; progress = progress) :
            npa_moments_block_nc(ops_principal, model; progress = progress)
    else
        PMS = []
        mons_D = Dict{AbstractMonomial,AbstractJuMPScalar}([])
        for (i, ops_principal_i) in enumerate(ops_principal)
            principal_moments_matrix_i, mons_D =
                tracial ?
                cyclic_npa_moments_block_nc(ops_principal_i, model; progress = progress,mons_D=mons_D) :
                npa_moments_block_nc(ops_principal_i, model; progress = progress,mons_D=mons_D)
            push!(PMS, principal_moments_matrix_i)
        end
    end
    # Add the constraints for the principal moment matrix
    if !isempty(tr_eq)
        for i in 1:length(tr_eq)
            tr_eq_p=0
            if tracial
                for (m, c) in Polynomial(tr_eq[i][1])
                    m1, m2 = (cyclic_reduce(m), cyclic_reduce(m'))
                    var = haskey(mons_D, m1) ? mons_D[m1] : mons_D[m2]
                    tr_eq_p += c*var
                end
            else
                tr_eq_poly=real_rep(Polynomial(Polynomial(tr_eq[i][1])))
                for (m, c) in tr_eq_poly
                    var = mons_D[m]
                    tr_eq_p+=c*var
                end
            end
            @constraint(model, tr_eq_p == tr_eq[i][2])
        end
    end
    if !isempty(tr_ge)
        for i in 1:length(tr_ge)
            tr_ge_p=0
            if tracial
                for (m, c) in Polynomial(tr_ge[i][1])
                    m1, m2 = (cyclic_reduce(m), cyclic_reduce(m'))
                    var = haskey(mons_D, m1) ? mons_D[m1] : mons_D[m2]
                    tr_ge_p += c*var
                end
            else
                tr_ge_poly=real_rep(Polynomial(Polynomial(tr_ge[i][1])))
                for (m, c) in tr_ge_poly
                    var = mons_D[m]
                    tr_ge_p+=c*var
                end
            end
            @constraint(model, tr_ge_p >= tr_ge[i][2])
        end
    end
    if !isempty(op_ge)
        for i in 1:length(op_ge)
            if tracial
                cyclic_npa_moments_block_nc(
                    ops,
                    model;
                    cPoly = op_ge[i],
                    mons_D = mons_D,
                    progress = progress,
                )
            else
                npa_moments_block_nc(
                    ops,
                    model;
                    cPoly = op_ge[i],
                    mons_D = mons_D,
                    progress = progress,
                )
            end
        end
    end
    if !isempty(op_eq)
        for i in 1:length(op_eq)
            if tracial
                cyclic_npa_moments_block_nc(
                    ops,
                    model;
                    cPoly = op_eq[i],
                    mons_D = mons_D,
                    eq = true,
                    progress = progress,
                )

            else
                npa_moments_block_nc(
                    ops,
                    model;
                    cPoly = op_eq[i],
                    mons_D = mons_D,
                    eq = true,
                    progress = progress,
                )
            end
        end
    end
    if normalize
        eltype(ops_principal)<:AbstractMonomial ? id_elem=one(first(ops_principal)) : id_elem=one(prod([first(i) for i in ops_principal]))
        if tracial
            id_elem=cyclic_reduce(id_elem)
            @constraint(model, mons_D[id_elem]==1.0)
        else
            @constraint(model, mons_D[id_elem]==1.0)
        end
    end
    if !is_number(obj)
        obj_p=0
        if tracial
            for (m, c) in Polynomial(obj)
                m1, m2 = (cyclic_reduce(m), cyclic_reduce(m'))
                var = haskey(mons_D, m1) ? mons_D[m1] : mons_D[m2]
                obj_p += c*var
            end
        else
            obj_poly=real_rep(Polynomial(Polynomial(obj)))
            for (m, c) in obj_poly
                var = mons_D[m]
                obj_p+=c*var
            end
        end

        min ? @objective(model, Min, obj_p) : @objective(model, Max, obj_p)
    end

    if eltype(ops_principal) <: AbstractMonomial
        return model, mons_D, principal_moments_matrix
    else
        return model, mons_D, PMS
    end
end
