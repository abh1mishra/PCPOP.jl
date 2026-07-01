function model_new_obj(model,S,V,mons,LMI,s)
    delete.(model,V)
    V=cons_gen!(LMI, S, s,model,mons)
    return model,V
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


function cons_gen!(LMI, obj_poly, s,model,mons)
    V = Vector([])
    for m in mons
        constr = sum((s*LMI[i][m] for i in 1:length(LMI) if haskey(LMI[i], m)), init=0.0) - coefficient(obj_poly, m)
        push!(V, @constraint(model, constr == 0))
    end
    return V
end

function npa_dual(obj, ops,ops_principal;
    min=false,
    op_eq = [], 
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    tracial=false,
    normalize=true,
    change_objective=false)

    model=Model()

    Id = one(first(ops_principal))

    pair_mat_ops = [(Id,ops_principal)]

    pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_ge])

    pair_mat_ops=vcat(pair_mat_ops,[(g,ops) for g in op_eq])
    pair_mat_ops=vcat(pair_mat_ops,[(-g,ops) for g in op_eq])


    tsize = [length(i) for (j,i) in pair_mat_ops]

    if normalize
        push!(tr_eq,[Id,1.0])
    end
    # Process the tracial constraints 


    Zs = [@variable(model, [1:l, 1:l], PSD) for l in tsize]
    Zmeq = [@variable(model) for i in tr_eq]
    Zmge = [@variable(model, lower_bound=0) for i in tr_ge]

    LMI = [tracial ? cyclic_npa_moment(ops_i,Zs[i];cPoly=g) : npa_moment(ops_i,Zs[i];cPoly=g) for (i,(g, ops_i)) in enumerate(pair_mat_ops)]

    min ? s=-1 : s=1
    mons = union!(Set(), [keys(lmi) for lmi in LMI]...)
    if !isempty(tr_ge)
        if tracial
            tr_ge = trace_mons_reduce(mons, tr_ge)
        else
            tr_ge = [[real_rep(Polynomial(tr_ge[i][1])),tr_ge[i][2]] for i in 1:length(tr_ge)]
        end
    end

    if !isempty(tr_eq)
        if tracial            
            tr_eq = trace_mons_reduce(mons, tr_eq)
        else
            tr_eq = [[real_rep(Polynomial(tr_eq[i][1])),tr_eq[i][2]] for i in 1:length(tr_eq)]
        end
    end
    obj_poly = sum(Zmeq[i]*tr_eq[i][2] for i in 1:length(tr_eq); init=0)
    obj_poly += sum(-s*Zmge[i]*tr_ge[i][2] for i in 1:length(tr_ge); init=0)

    if !(obj isa Number) && !iszero(Polynomial(obj))
        if tracial
            obj = trace_mons_reduce(mons,[obj])[1]
        else
            obj = real_rep(Polynomial(obj))
        end
    end
    
    S = -obj
    S += sum(Zmeq[i]*tr_eq[i][1] for i in 1:length(tr_eq); init=0)
    S += sum(-s*Zmge[i]*tr_ge[i][1] for i in 1:length(tr_ge); init=0)

    min ? @objective(model, Max, obj_poly) : @objective(model, Min, obj_poly)

    V=cons_gen!(LMI, S, s,model,mons)

    if change_objective
        return model,S,V,mons,LMI
    end

    return model, LMI[1],Zs

end
