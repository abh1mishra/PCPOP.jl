function sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end

    return matrix
end

function npa_block(list_monomials;cPoly=1)
    num_monomials = length(list_monomials)
    moments_D = Dict{AbstractMonomial,SparseMatrixCSC{Float64, Int64}}([])
    N = num_monomials
    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials
            # Compute the product of the monomials
            monomial_product = real_rep(Polynomial(list_monomials[i]'* cPoly * list_monomials[j]))
            # Initialize the JuMP variable for the matrix entry
            for (m,c) in monomial_product
                # Check if the product already exists in the dictionary
                if haskey(moments_D,m)
                    sym_add!(moments_D[m], i, j, c)
                else
                    moments_D[m] = sym_add!(spzeros(N, N), i, j, c)
                end
            end
        end
    end
    return moments_D
end