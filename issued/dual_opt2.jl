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


    if !isempty(tr_ge)
        for (g, val) in tr_ge
            tr_ge_poly = g-val
            pair_mat_ops=vcat(pair_mat_ops,[(tr_ge_poly,[Id])])
        end
    end

    if !isempty(tr_eq)
        for (g, val) in tr_eq
            tr_eq_poly = g-val
            pair_mat_ops=vcat(pair_mat_ops,[(tr_eq_poly,[Id])])
            pair_mat_ops=vcat(pair_mat_ops,[(-tr_eq_poly,[Id])])
        end
    end

    tsize = [length(i) for (j,i) in pair_mat_ops]

    Zs=[@variable(model, [1:l, 1:l], PSD) for l in tsize]

    LMI = [tracial ? cyclic_npa_moment(ops_i,Zs[i];cPoly=g) : npa_moment(ops_i,Zs[i];cPoly=g) for (i,(g, ops_i)) in enumerate(pair_mat_ops)]

    min ? s=-1 : s=1
    tracial ? Id = cyclic_reduce(Id) : Id = Id

    obj_poly = obj*Id
    tracial ? obj_poly = Polynomial(obj_poly) : obj_poly = real_rep(Polynomial(obj_poly))

    obj_p = obj_gen(obj_poly, LMI, Id, s)
    min ? @objective(model, Max, obj_p) : @objective(model, Min, obj_p)

    cons_gen!(LMI, obj_poly, Id, s,model)

    return model, LMI,Zs

end
