# Taken from QuantumNPA.jl with some modifications

function LinearAlgebra.dot(A::SparseMatrixCSC,
                           B::Symmetric{<:JuMP._MA.AbstractMutable,
                                        Matrix{<:JuMP._MA.AbstractMutable}})
    acc = zero(eltype(B))

    for j in 1:size(A, 2)
        for k in nzrange(A, j)
            add_to_expression!(acc, nonzeros(A)[k], B[rowvals(A)[k], j])
        end
    end

    return acc
end

function LinearAlgebra.dot(A::Symmetric{<:JuMP._MA.AbstractMutable,
                                        Matrix{<:JuMP._MA.AbstractMutable}},
                           B::SparseMatrixCSC)
    return dot(B, A)
end

function sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end

    return matrix
end


function npa_moment(operators::Vector,Zsi::Symmetric{VariableRef, Matrix{VariableRef}};cPoly=1)
    N = length(operators)
    iops = collect(enumerate(operators))
    moment_D = Dict{AbstractMonomial,AbstractJuMPScalar}([])
    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(real_rep(conj(x)*cPoly*y))
            for (m,c) in p
                term = c*Zsi[i, j]
                if i != j
                    term += c*Zsi[j, i]
                end
                
                if haskey(moment_D, m)
                    moment_D[m] += term
                else
                    moment_D[m] = term
                end
            end
        end
    end

    return moment_D
end


function cyclic_npa_moment(operators::Vector,Zsi::Symmetric{VariableRef, Matrix{VariableRef}};cPoly=1)
    N = length(operators)
    iops = collect(enumerate(operators))
    moment_D = Dict{AbstractMonomial,AbstractJuMPScalar}([])
    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(conj(x)*cPoly*y)
            for (m,c) in p
                if m isa Variable
                    println("error ",x, "  ", y, "  ", cPoly)
                end
                m_,m__ = cyclic_reduce(m), cyclic_reduce(m')
                
                term = c*Zsi[i, j]
                if i != j
                    term += c*Zsi[j, i]
                end
                
                if haskey(moment_D, m_)
                    moment_D[m_] += term
                elseif haskey(moment_D, m__)
                    moment_D[m__] += term
                else
                    moment_D[m_] = term
                end
            end
        end
    end

    return moment_D
end

function obj_gen(obj_poly, LMI, Id, s)
    obj_p = Main.coefficient(obj_poly,Id) + sum((s*LMI[i][Id] for i in 1:length(LMI) if haskey(LMI[i], Id)), init=0.0)        
    return obj_p
end

function cons_gen!(LMI, obj_poly, Id, s,model)
    mons = union!(Set(), [keys(lmi) for lmi in LMI]...)
    for m in mons
        if m != Id
            constr = sum((s*LMI[i][m] for i in 1:length(LMI) if haskey(LMI[i], m)), init=0.0) + Main.coefficient(obj_poly, m)
            @constraint(model, constr == 0)
        end
    end
end

function npa_dual(obj, ops,ops_principal;
    min=false,
    op_eq = [], 
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    tracial=false,
    normalize=true)
    println("npa_dual in action")
    println("Number of operators in the principal moment matrix: ", length(ops_principal))
    model=Model()

    Id = one(first(ops_principal))

    pair_mat_ops = [(Id,ops_principal)]

    pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_ge])

    pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_eq])
    pair_mat_ops=vcat(pair_mat_ops,[(-g,ops) for g in op_eq])


    tsize = [length(i) for (j,i) in pair_mat_ops]

    # Process the tracial constraints 
    if !isempty(tr_ge)
        tr_ge = [[real_rep(Polynomial(tr_ge[1])),tr_ge[2]]]
    end

    if !isempty(tr_eq)
        tr_eq = [[real_rep(Polynomial(tr_eq[1])),tr_eq[2]]]
    end

    Zs = [@variable(model, [1:l, 1:l], PSD) for l in tsize]
    Zmeq = [@variable(model) for i in tr_eq]
    Zmge = [@variable(model, lower_bound=0) for i in tr_ge]

    LMI = [tracial ? cyclic_npa_moment(ops_i,Zs[i];cPoly=g) : npa_moment(ops_i,Zs[i];cPoly=g) for (i,(g, ops_i)) in enumerate(pair_mat_ops)]

    min ? s=-1 : s=1
    tracial ? Id = cyclic_reduce(Id) : Id = Id

    obj_poly = obj*Id
    obj_poly += sum(Zmeq[i]*tr_eq[i][2] for i in 1:length(tr_eq); init=0)
    obj_poly += sum(Zmge[i]*tr_ge[i][2] for i in 1:length(tr_ge); init=0)

    S = -obj
    S += sum(tr_eq[i][1] for i in 1:length(tr_eq); init=0)
    S += sum(tr_ge[i][1] for i in 1:length(tr_ge); init=0)

    tracial ? obj_poly = Polynomial(obj_poly) : obj_poly = real_rep(Polynomial(obj_poly))

    obj_p = obj_gen(obj_poly, LMI, Id, s,[tr_eq;tr_ge])
    min ? @objective(model, Max, obj_p) : @objective(model, Min, obj_p)

    cons_gen!(LMI, obj_poly, Id, s,model,[tr_eq;tr_ge])

    return model, LMI,Zs

end
