@testset "Jordan reduction" begin
    diagonalize = true # whether to diagonalize the moment matrix or not, this is a parameter of the jordan reduction method and it will not change the values of the plot (they use diagonalization in their code)
    complex = true
    G=0.8  # value of the guessing probability
    @pcmonoid M a b c d e f
    ρ=[b, c, d]
    σ=a
    PA=[e, f]
    Projector.([e, f])
    build(M)
    PA=[e f; 1-e 1-f]
    ge = [σ-1/3*ρ[i] for i in 1:3]
    ge=vcat(ge, [i-i*i for i in ρ])
    ge_constr = [[-σ, -G]]
    eq_constr = [[ρ[x], 1] for x in 1:3]
    wit =
        -(2*PA[1, 1]-1)*ρ[1]-(2*PA[1, 2]-1)*ρ[1]-(2*PA[1, 1]-1)*ρ[2]+(2*PA[1, 2]-1)*ρ[2]+(
            2*PA[1, 1]-1
        )*ρ[3]
    level = 2 # this is the level of the localising matrices and it will be enough to retrieve the values of the plot (they use level 3 of the principal moment matrix)
    D, C, A, b=npa_canonical(
        wit,
        level;
        op_ge = ge,
        tr_eq = eq_constr,
        tr_ge = ge_constr,
        cyclic = true,
        min = false,
        normalize = false,
        rm = true,
    )
    model, P, blkD =
        jordan_reduce(C, A, b, verbose = true, complex = complex, diagonalize = diagonalize)
end

@testset "routed bell" begin
    @pcmonoid M UA[2, 0] UBS[2, 0] UBL[2, 0]
    Unipotent.(M.vertices)
    @comms UA UBS
    @comms UA UBL

    # SRQ constraint, joint measuribility of BL
    @comms UBL

    build(M)

    # C_s = 2\sqrt(2) constraint

    Cs = UA[1]*UBS[1] + UA[1]*UBS[2] + UA[2]*UBS[1] - UA[2]*UBS[2]
    tr_eq = [[Cs, 2*sqrt(2)]]

    θ = 1/2
    level=2

    obj = tan(θ)*UA[1]*UBL[1] + UA[1]*UBL[2] + UA[2]*UBL[1] - tan(θ)*UA[2]*UBL[2]
    Γ, Cvec, Amat, Bvec=npa_canonical(obj, level; tr_eq = tr_eq, min = false, rm = true)
    diagonalize=true
    complex=true
    model_red, P, blkD = jordan_reduce(
        Cvec,
        Amat,
        Bvec,
        verbose = true,
        complex = complex,
        diagonalize = diagonalize,
    )
end

@testset "routed bell2" begin
    @pcmonoid M A[4, 0] BS[4, 0] BL[4, 0]
    Projector.(M.vertices)
    @comms A BS
    @comms A BL
    @ortho A[1] A[2]
    @ortho A[3] A[4]
    @ortho BS[1] BS[2]
    @ortho BS[3] BS[4]
    @ortho BL[1] BL[2]
    @ortho BL[3] BL[4]
    # SRQ constraint, joint measuribility of BL
    @comms BL

    build(M)
    PBS = [BS[1] BS[3]; BS[2] BS[4]; 1-BS[1]-BS[2] 1-BS[3]-BS[4]]
    PBL = [BL[1] BL[3]; BL[2] BL[4]; 1-BL[1]-BL[2] 1-BL[3]-BL[4]]
    PA = [A[1] A[3]; A[2] A[4]; 1-A[1]-A[2] 1-A[3]-A[4]]

    level = "2+A*A*A+BS*BS*BS+BL*BL*BL+A*A*BL+A*A*BS+BS*BS*BL"
    # level=2

    Γ, Cmat, Amat, Xmat, Bvec =
        npa_dual(0, level; min = false, rm = true, list_vars = M.vertices)
    ηl = @variable(model)
    ηs = 0.98
    ηa=ηs
    Γ = Dict([k=>Xmat[v...] for (k, v) in Γ])
    [
        @constraint(model, Γ[PA[a, x] * PBS[b, y]] == ηa*ηs*chsh_correlations(a, b, x, y))
        for a in 1:2 for x in 1:2 for b in 1:2 for y in 1:2
    ]
    id = one(A[1])
    [
        @constraint(model, Γ[PA[a, x] * id] == ηa*chsh_correlations(a, -1, x, 1)) for
        a in 1:2 for x in 1:2
    ]

    [
        @constraint(model, Γ[PBS[b, y] * id] == ηs*chsh_correlations(-1, b, 1, y)) for
        b in 1:2 for y in 1:2
    ]

    [
        @constraint(model, Γ[PA[a, x] * PBL[b, y]] == ηa*ηl*chsh_correlations(a, b, x, y))
        for a in 1:2 for x in 1:2 for b in 1:2 for y in 1:2
    ]

    [
        @constraint(model, Γ[PBL[b, y] * id] == ηl*chsh_correlations(-1, b, 1, y)) for
        b in 1:2 for y in 1:2
    ]

    @objective(model, Max, ηl)
    optimize!(model)
    println("Optimal value is ", value(ηl))
end

@testset "routed bell 3" begin
    @subsystem M A[2, 0] B[4, 0]
    add_relations!([i^3-i for i in a])
    add_relations!([i^3-i for i in b])
    bs = b[1:2]
    bl = b[3:4]

    # SRQ constraint, joint measuribility of BL
    add_relations!([bl[1]*bl[2]-bl[2]*bl[1]])

    build(M)

    level = "3+A*A*A*A+b[1:2]*b[1:2]*b[1:2]*b[1:2]+b[3:4]*b[3:4]*b[3:4]*b[3:4]+A*A*b[1:2]*b[1:2]+A*A*b[3:4]*b[3:4]+b[1:2]*b[1:2]*b[3:4]*b[3:4]"
    # level = 3
    # model, Γ, pm = npa(0,level;min=false,rm=true,list_vars = [a;b],optimizer=SDPA.Optimizer,model_flags=[("Mode",SDPA.PARAMETER_STABLE_BUT_SLOW)])
    Γ, model, Xmat, pm = npa(0, level; min = false, rm = true, list_vars = [a; b])

    ηl = @variable(model)
    ηs = 0.98
    ηa=ηs
    Γ = Dict([k=>Xmat[v...] for (k, v) in Γ])
    # AxBsy
    [
        @constraint(
            model,
            Γ[a[x] * bs[y]] == ηa*ηs*chsh_correlations(1, 1, x, y; obs = true)
        ) for x in 1:2 for y in 1:2
    ]#4

    #Ax^2Bsy
    [
        @constraint(
            model,
            Γ[a[x] ^ 2 * bs[y]] == ηa*ηs*chsh_correlations(1, 1, -1, y; obs = true)
        ) for x in 1:2 for y in 1:2
    ]#4
    #AxBsy^2
    [
        @constraint(
            model,
            Γ[a[x] * bs[y] ^ 2] == ηa*ηs*chsh_correlations(1, 1, x, -1; obs = true)
        ) for x in 1:2 for y in 1:2
    ]#4
    #Ax^2Bsy^2
    [@constraint(model, Γ[a[x] ^ 2 * bs[y] ^ 2] == ηa*ηs) for x in 1:2 for y in 1:2]#4

    id = monomial(one(monomial(a[1])))

    #Ax
    [
        @constraint(model, Γ[a[x] * id] == ηa*chsh_correlations(1, -1, x, -1; obs = true))
        for x in 1:2
    ]#2
    #Ax^2
    [@constraint(model, Γ[a[x] ^ 2 * id] == ηa) for x in 1:2]#2

    #Bsy
    [
        @constraint(
            model,
            Γ[bs[y] * id] == ηs*chsh_correlations(-1, -1, -1, y; obs = true)
        ) for y in 1:2
    ]#2
    #Bsy^2
    [@constraint(model, Γ[bs[y] ^ 2 * id] == ηs) for y in 1:2]#2

    #AxBly
    [
        @constraint(
            model,
            Γ[a[x] * bl[y]] == ηa*ηl*chsh_correlations(1, 1, x, y; obs = true)
        ) for x in 1:2 for y in 1:2
    ]#4

    #Ax^2Bly
    [
        @constraint(
            model,
            Γ[a[x] ^ 2 * bl[y]] == ηa*ηl*chsh_correlations(1, 1, -1, y; obs = true)
        ) for x in 1:2 for y in 1:2
    ]#4
    #AxBly^2
    [
        @constraint(
            model,
            Γ[a[x] * bl[y] ^ 2] == ηa*ηl*chsh_correlations(1, 1, x, -1; obs = true)
        ) for x in 1:2 for y in 1:2
    ]#4
    #Ax^2Bly^2
    [@constraint(model, Γ[a[x] ^ 2 * bl[y] ^ 2] == ηa*ηl) for x in 1:2 for y in 1:2]#4

    #Bly
    [
        @constraint(
            model,
            Γ[bl[y] * id] == ηl*chsh_correlations(-1, -1, -1, y; obs = true)
        ) for y in 1:2
    ]#2
    #Bly^2
    [@constraint(model, Γ[bl[y] ^ 2 * id] == ηl) for y in 1:2]#2

    @objective(model, Max, ηl)
    optimize!(model)
    println("Optimal value is ", value(ηl))
end
