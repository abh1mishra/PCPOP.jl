function cyclic_npa_moments_block!(list_monomials::Vector{M},X,tsize,model;cPoly=1,mons_pos_D=Dict([]),offset=0,extra_zeros=false) where M<:AbstractMonomial

    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Initialize the matrix of JuMP variables

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials
            println("i: ", i, " j: ", j)
            if cPoly == 1
                m = list_monomials[i]'* list_monomials[j]
                m1,m2 = cyclic_reduce(m),cyclic_reduce(m')
                if m1==0 || m2==0
                    @constraint(model, X[i,j] == 0)
                    continue
                end
                if  haskey(mons_pos_D,m1) || haskey(mons_pos_D,m2)
                    # Use the existing JuMP variable
                    upi,upj = haskey(mons_pos_D,m1) ? mons_pos_D[m1] : mons_pos_D[m2]
                    @constraint(model, X[i,j] - X[upi,upj] == 0)

                else
                    mons_pos_D[m1]= (i,j)
                end
            else
                monomial_product = list_monomials[i]'* cPoly * list_monomials[j]
                if monomial_product == 0
                    @constraint(model, X[offset+i,offset+j] == 0)
                    continue
                end
                PMpoly = 0
                for (m,c) in monomial_product
                    m1,m2 = cyclic_reduce(m),cyclic_reduce(m')
                    if !haskey(mons_pos_D,m1) && !haskey(mons_pos_D,m2)
                        throw(ArgumentError("level not enough"))
                    end
                    upi,upj = haskey(mons_pos_D,m1) ? mons_pos_D[m1] : mons_pos_D[m2]
                    PMpoly += c*X[upi,upj]
                end
                @constraint(model, X[offset+i,offset+j] - PMpoly == 0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                @constraint(model, X[offset+i , j] == 0.0)
            end
        end

        for i in 1:num_monomials
            for j in offset+num_monomials+1:tsize
                @constraint(model, X[offset+i , j] == 0.0)
            end
        end
    end

    if cPoly == 1
        return mons_pos_D
    end

end

function npa_moments_block!(list_monomials::Vector{M},X,tsize,model;cPoly=1,mons_pos_D=Dict([]),offset=0,extra_zeros=false) where M<:AbstractMonomial
    # Get the number of monomials
    num_monomials = length(list_monomials)

    # Iterate over the list of monomials to fill the matrix
    for i in 1:num_monomials
        for j in i:num_monomials

            if cPoly == 1
                m = real_rep(list_monomials[i]'* list_monomials[j])
                if m == 0
                    @constraint(model, X[i,j] == 0)
                    continue
                end
                if haskey(mons_pos_D,m)
                    # Use the existing JuMP variable
                    upi,upj = mons_pos_D[m]
                    @constraint(model, X[i,j] - X[upi,upj] == 0)
                else
                    mons_pos_D[m] = (i,j)
                end
            else
                monomial_product = real_rep(Polynomial(list_monomials[i]'* cPoly * list_monomials[j]))
                # assume the monomials PM exists and m and m' are in pm.
                if monomial_product == 0
                    @constraint(model, X[offset+i,offset+j] == 0)
                    continue
                end
                PMpoly = 0
                for (m,c) in monomial_product
                    upi,upj = mons_pos_D[m]
                    PMpoly += c*X[upi,upj]
                end
                @constraint(model, X[offset+i,offset+j] - PMpoly == 0)
            end
        end
    end

    if extra_zeros
        for i in 1:num_monomials
            for j in 1:offset
                @constraint(model,X[offset+i , j] == 0.0)
            end
        end

        for i in 1:num_monomials
            for j in offset+num_monomials+1:tsize
                @constraint(model,X[offset+i , j] == 0.0)
            end
        end
    end

    if cPoly == 1
        return mons_pos_D
    end
end
function npa(obj, ops, ops_principal;
    min=true,
    op_eq = [], 
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    tracial=false,
    normalize=true,
    extra_zeros=false
    )

    println("Number of operators in the principal moment matrix: ", length(ops_principal))

    model=Model()
    tsize=sum([[length(ops_principal)];[length(ops) for i in 1:length(op_ge)];[2*length(ops) for i in 1:length(op_eq)]])

    X = @variable(model, [1:tsize, 1:tsize], PSD)
    mons_pos_D = tracial ? cyclic_npa_moments_block!(ops_principal,X,tsize,model; extra_zeros=extra_zeros) : npa_moments_block!(ops_principal,X,tsize,model; extra_zeros=extra_zeros)
    offset = length(ops_principal)
    for i in 1:length(op_ge)
        if tracial
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=op_ge[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
        else
            npa_moments_block!(ops,X,tsize,model; cPoly=op_ge[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
        end
        offset += length(ops)
    end


    for i in 1:length(op_eq)
        if tracial
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            cyclic_npa_moments_block!(ops,X,tsize,model; cPoly=-op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops) 

        else
            npa_moments_block!(ops,X,tsize,model; cPoly=op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros) 
            offset += length(ops)
            npa_moments_block!(ops,X,tsize,model; cPoly=-op_eq[i], mons_pos_D=mons_pos_D, offset=offset, extra_zeros=extra_zeros)
            offset += length(ops)
        end
    end

    # Add the constraints for the principal moment matrix

    for i in 1:length(tr_eq)
        tr_eq_p = 0
        if tracial
            for (m,c) in Polynomial(tr_eq[i][1])
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                upi,upj = haskey(mons_pos_D,m1) ? mons_pos_D[m1] : mons_pos_D[m2]
                tr_eq_p += c*X[upi,upj]
            end

        else
            tr_eq_poly=real_rep(Polynomial(tr_eq[i][1]))
            for (m,c) in tr_eq_poly
                upi,upj = mons_pos_D[m]
                tr_eq_p += c*X[upi,upj]
            end
        end
        try
            @constraint(model, tr_eq_p - tr_eq[i][2] == 0)
        catch e
            println("error",tr_eq_p, " ", tr_eq[i][2])
        end
    end

    for i in 1:length(tr_ge)
        tr_ge_p=0
        if tracial
            for (m,c) in Polynomial(tr_ge[i][1])
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                upi,upj = haskey(mons_pos_D,m1) ? mons_pos_D[m1] : mons_pos_D[m2]
                tr_ge_p+=c*X[upi,upj]
            end
        else
            tr_ge_poly=real_rep(Polynomial(Polynomial(tr_ge[i][1])))
            for (m,c) in tr_ge_poly
                upi,upj = mons_pos_D[m]
                tr_ge_p+=c*X[upi,upj]
            end
        end
        @constraint(model, tr_ge_p >= tr_ge[i][2])
    end


    if normalize
        id_elem=one(first(ops_principal))
        if tracial
            id_elem=cyclic_reduce(id_elem)
        end
        upi,upj = mons_pos_D[id_elem]   
        @constraint(model, X[upi,upj] == 1)
    end
    obj_p = 0
    if !is_number(obj)
        if tracial
            for (m,c) in Polynomial(obj)
                m1,m2= (cyclic_reduce(m),cyclic_reduce(m'))
                if !haskey(mons_pos_D,m1) && !haskey(mons_pos_D,m2)
                    throw(ArgumentError("level not enough"))
                end
                upi,upj = haskey(mons_pos_D,m1) ? mons_pos_D[m1] : mons_pos_D[m2]
                obj_p += c*X[upi,upj]
            end
        else
            obj_poly=real_rep(Polynomial(obj))
            for (m,c) in obj_poly
                if !haskey(mons_pos_D,m)
                    throw(ArgumentError("level not enough"))
                end
                upi,upj = mons_pos_D[m]
                obj_p += c*X[upi,upj]
            end
        end
        min ? @objective(model, Min, obj_p) : @objective(model, Max, obj_p)
    end

    mons_pos_D = Dict([k => X[v...] for (k,v) in mons_pos_D])
    return model,mons_pos_D,X    
end
