

############################
##  OPTIMIZATION METHODS  ##
############################

function pcpop(p::Polynomial, k::Int;
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    normalize = true,
    tracial = false,
    solver = Mosek.Optimizer,
    optimize = false,
    reduce = false,
    block_diag = false,
    primal = true)

    basis = mons_at_level(p.monoid, k)
    return pcpop(p, basis, op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           solver = solver,
                           optimize = optimize,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
end

function pcpop(p::Polynomial, basis;
    op_eq = [],
    op_ge = [],
    tr_eq = [],
    tr_ge = [],
    normalize = true,
    tracial = false,
    solver = Mosek.Optimizer,
    optimize = false,
    reduce = false,
    block_diag = false,
    primal = true)

    if reduce
        Γ, C, A, b = npa_canonical(p, basis, 
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           solver = solver,
                           optimize = optimize,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
        model = jordan_reduce(C, A, b, complex=true, diagonalize=block_diag)
    elseif primal
        model = npa(p, basis, 
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           solver = solver,
                           optimize = optimize,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
    else
        model = sos(p, basis, 
                           op_eq = op_eq,
                           op_ge = op_ge,
                           tr_eq = tr_eq,
                           tr_ge = tr_ge,
                           normalize = normalize,
                           tracial = tracial,
                           solver = solver,
                           optimize = optimize,
                           reduce = reduce,
                           block_diag = block_diag,
                           primal = primal)
    end

    if optimize
        optimize!(model)
    end

    return model
end