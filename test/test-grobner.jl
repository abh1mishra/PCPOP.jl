@testset "Macaulay Reduction       " begin
    @pcmonoid M a b c
    build(M)

    G = [a*b - b*a, b*c - c*b]

    @test reduce_grobner(a*b*c - a*c*b, b*c - c*b) == 0
    @test reduce_grobner(a*b*c - a*c*b, G) == 0
    @test reduce_grobner([a*b*c - a*c*b, b*c - c*b], b*c - c*b) == []
    @test reduce_grobner(G, G) == []
    @test G == self_reduce(G)

    H = [a*b - b*a, b*c - c*b, a*b*c - a*c*b]

    @test reduce_grobner(H, G) == []
    @test 1.0*G == self_reduce(1.0*H)
end

@testset "Macaulay Optimization    " begin
    @pcmonoid M a b c d
    Projector.([a, b, c, d])
    build(M)

    G = [a*c - c*a, a*d - d*a, b*c - c*b, b*d - d*b]

    f = a*c + a*d + b*c - b*d

    model = pcpop(f)
    set_optimizer(model, default_solver())
    set_silent(model)
    optimize!(model)

    @test termination_status(model) == MOI.TerminationStatusCode(2)

    model_localized = pcpop(f, equalities = G)
    set_optimizer(model_localized, default_solver())
    set_silent(model_localized)
    optimize!(model_localized)

    @test termination_status(model_localized) == MOI.TerminationStatusCode(1)

    # Compare with partially commuting setting
    @pcmonoid M a b c d

    @comms [a, b] [c, d]
    Projector.([a, b, c, d])

    build(M)

    f = a*c + a*d + b*c - b*d

    model_pc = pcpop(f, 1)
    set_optimizer(model_pc, default_solver())
    set_silent(model_pc)
    optimize!(model_pc)
    @test termination_status(model_pc) == MOI.TerminationStatusCode(1)
    @test abs(objective_value(model_pc) - objective_value(model_localized)) <= 1e-6

    # Second level ncpop + Macaulay
    @pcmonoid M a b c d
    Projector.([a, b, c, d])
    build(M)
    G = [a*c - c*a, a*d - d*a, b*c - c*b, b*d - d*b]
    f = a*c + a*d + b*c - b*d
    global model_localized = pcpop(f, 2, equalities = G)
    set_optimizer(model_localized, default_solver())
    set_silent(model_localized)
    optimize!(model_localized)
    @test termination_status(model_localized) == MOI.TerminationStatusCode(1)

    # Second level pcpop
    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Projector.([a, b, c, d])
    build(M)
    f = a*c + a*d + b*c - b*d
    global model_pc = pcpop(f, 2)
    set_optimizer(model_pc, default_solver())
    set_silent(model_pc)
    optimize!(model_pc)
    @test termination_status(model_pc) == MOI.TerminationStatusCode(1)
    @test objective_value(model_pc) <= 2+1e-6
    @test abs(objective_value(model_pc) - objective_value(model_localized)) <= 1e-6
end
