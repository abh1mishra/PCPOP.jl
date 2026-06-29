@testset "Polynomial Optimization  " begin

    @pcmonoid M a b c d

    @comms [a, b] [c, d]
    Projector.([a,b,c,d])

    build(M)

    f = a*c + a*d + b*c - b*d
    
    val,sos_model,_ = pcpop(f, 1)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)

    val,sos_model,_ = pcpop(f, 2)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)
    @test abs(2-objective_value(sos_model)) <= 1e-6
    
end

@testset "Symmetry Reduction       " begin

    @pcmonoid M a b c d

    @comms [a, b] [c, d]
    Projector.([a,b,c,d])

    build(M)

    f = a*c + a*d + b*c - b*d
    
    action = OnLetters()
    p1 = PG.perm"(1,3)(2,4)" #"
    p2 = PG.perm"(1,4)(2,3)" #"
    G = PG.PermGroup(p1, p2)

    sa_model = pcpop(f, 2, G, action)
    set_optimizer(sa_model, Mosek.Optimizer)
    set_silent(sa_model)
    optimize!(sa_model);
    @test termination_status(sa_model) == MOI.TerminationStatusCode(1)
    @test abs(2-objective_value(sa_model)) <= 1e-6
    dual_vars = all_constraints(sa_model, AffExpr, MOI.EqualTo{Float64})
    @test length(dual_vars) == 13
    
end

@testset "Block Diagonalization    " begin
   
    @pcmonoid M a b c d

    @comms [a, b] [c, d]
    Projector.([a,b,c,d])

    build(M)

    f = a*c + a*d + b*c - b*d
    
    action = OnLetters()
    p1 = PG.perm"(1,3)(2,4)" #"
    p2 = PG.perm"(1,4)(2,3)" #"
    G = PG.PermGroup(p1, p2)
    T = Float64
    
    monomial_basis = mons_at_level(f, 2)
    basis_constraints = StarAlgebras.Basis{UInt16}(
              unique([x*y for x in monomial_basis for y in monomial_basis]),
    )
    
    wedderburn = SymbolicWedderburn.WedderburnDecomposition(
            T,
            G,
            action,
            basis_constraints,
            monomial_basis,
            semisimple=true,
        )
    
    U  = wedderburn.Uπs
    SU =[size(u,1) for u in U]
    @test length(U) == 4
    @test SU == [5,3,3,2]
    
    blk_model = pcpop(f, 2, G, action, diagonalize=true)
    set_optimizer(blk_model, Mosek.Optimizer)
    set_silent(blk_model)
    optimize!(blk_model);
    @test termination_status(blk_model) == MOI.TerminationStatusCode(1)
    @test abs(2-objective_value(blk_model)) <= 1e-6
    dual_vars = all_constraints(blk_model, AffExpr, MOI.EqualTo{Float64})
    @test length(dual_vars) == 13
    sdp_vars = all_constraints(blk_model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)
    @test 1 + sum([s*(s+1)/2 for s in SU]) == length(all_variables(blk_model))
    
end