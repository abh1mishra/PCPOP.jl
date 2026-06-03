function cyclic_npa_moments_block_nc(list_monomials::Vector{M},model;cPoly=1,unique_mons=[],unique_vars=[],eq=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables
    moments_matrix = Matrix{JuMP.AffExpr}(undef, num_monomials, num_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in 1:num_monomials
            
            # Compute the product of the monomials
            monomial_product = list_monomials[i]'* cPoly * list_monomials[j]
            # println("1",typeof(monomial_product))
            # Initialize the JuMP variable for the matrix entry
            moments_matrix[i, j] = 0.0
            # Check if the product already exists in the dictionary
            for (m,c) in monomial_product
                # println("2",typeof(m))
                m1,m2 =  (cyclic_reduce(m),cyclic_reduce(m'))
                if m1==0 || m2==0
                    continue
                end
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    moments_matrix[i, j] += c*unique_vars[m_i]
                else
                    # Create a new JuMP variable
                    new_var = @variable(model)
                    # Store the new variable in the dictionary
                    push!(unique_mons, m1)
                    push!(unique_vars, new_var)
                    # Use the new variable in the matrix
                    moments_matrix[i, j] += c*new_var
                end
            end
        end
    end

    if !eq
        @constraint(model, moments_matrix >=0,PSDCone())
    else
        @constraint(model, moments_matrix .== 0)
    end

    return moments_matrix, unique_mons, unique_vars
end

function npa_moments_block_nc(list_monomials::Vector{M},model;cPoly=1,unique_mons=[],unique_vars=[],eq=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables
    moments_matrix = Matrix{JuMP.AffExpr}(undef, num_monomials, num_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials
            # Compute the product of the monomials
            monomial_product = real_rep(Polynomial(list_monomials[i]'* cPoly * list_monomials[j]))
            # Initialize the JuMP variable for the matrix entry
            moments_matrix[i, j] = 0.0
            for (m,c) in monomial_product
                # Check if the product already exists in the dictionary
                m_i=findfirst(x->(x==m),unique_mons)
                if m_i !== nothing
                    # Use the existing JuMP variable
                    moments_matrix[i, j] += c*unique_vars[m_i]
                else
                    # Create a new JuMP variable
                    new_var = @variable(model)
                    # Store the new variable in the dictionary
                    push!(unique_mons, m)
                    push!(unique_vars, new_var)
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
    return moments_matrix, unique_mons, unique_vars
end

function npa_nc(obj, ops,ops_principal;
    min=true,
    op_eq = [], 
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    tracial=false,
    normalize=true)

    println("Number of operators in the principal moment matrix and LMI: ", length(ops_principal)," ",length(ops))
    model=Model()

    principal_moments_matrix, unique_mons, unique_vars = tracial ? cyclic_npa_moments_block_nc(ops_principal,model) : npa_moments_block_nc(ops_principal,model)
    # Add the constraints for the principal moment matrix

    if !isempty(tr_eq)
        for i in 1:length(tr_eq)
            tr_eq_p=0
            if tracial
                for (m,c) in Polynomial(tr_eq[i][1])
                    m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                    m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                    tr_eq_p+=c*unique_vars[m_i]
                end
            else
                tr_eq_poly=real_rep(Polynomial(Polynomial(tr_eq[i][1])))
                for (m,c) in tr_eq_poly
                    m_i=findfirst(x->x==m,unique_mons)
                    tr_eq_p+=c*unique_vars[m_i]
                end
            end
            @constraint(model, tr_eq_p == tr_eq[i][2])
        end
    end
    if !isempty(tr_ge)
        for i in 1:length(tr_ge)
            tr_ge_p=0
            if tracial
                for (m,c) in Polynomial(tr_ge[i][1])
                    m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                    m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                    tr_ge_p+=c*unique_vars[m_i]
                end
            else
                tr_ge_poly=real_rep(Polynomial(Polynomial(tr_ge[i][1])))
                for (m,c) in tr_ge_poly
                    m_i=findfirst(x->x==m,unique_mons)
                    tr_ge_p+=c*unique_vars[m_i]
                end
            end
            @constraint(model, tr_ge_p >= tr_ge[i][2])
        end
    end
    if !isempty(op_ge)
        for i in 1:length(op_ge)
            if tracial
                cyclic_npa_moments_block_nc(ops,model; cPoly=op_ge[i],unique_mons=unique_mons, unique_vars=unique_vars) 
            else
                npa_moments_block_nc(ops,model; cPoly=op_ge[i], unique_mons=unique_mons, unique_vars=unique_vars) 
            end
        end
    end
    if !isempty(op_eq)
        for i in 1:length(op_eq)
            if tracial
                cyclic_npa_moments_block_nc(ops,model; cPoly=op_eq[i],unique_mons=unique_mons, unique_vars=unique_vars,eq=true) 

            else
                npa_moments_block_nc(ops,model; cPoly=op_eq[i], unique_mons=unique_mons, unique_vars=unique_vars,eq=true) 

            end
        end
    end
    if normalize
        id_elem=one(first(ops_principal))
        if tracial
            id_elem=cyclic_reduce(id_elem)
            @constraint(model,unique_vars[findfirst(check->check==id_elem,unique_mons)]==1.0)
        else
            @constraint(model,unique_vars[findfirst(check->check==id_elem,unique_mons)]==1.0)
        end
    end
    if !is_number(obj)
        obj_p=0
        if tracial
            for (m,c) in Polynomial(obj)
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                m_i=findfirst(x->(x==m1 || x==m2),unique_mons)
                if m_i==nothing
                    throw(ArgumentError("level not enough"))
                end
                obj_p+=c*unique_vars[m_i]
            end
        else
            obj_poly=real_rep(Polynomial(Polynomial(obj)))
            for (m,c) in obj_poly
                m_i=findfirst(x->x==m,unique_mons)
                if m_i==nothing
                    throw(ArgumentError("level not enough"))
                end
                obj_p+=c*unique_vars[m_i]
            end
        end
        
        min ? @objective(model, Min, obj_p) : @objective(model, Max, obj_p)
    end

    return model,Dict(zip(unique_mons,unique_vars)),principal_moments_matrix

end