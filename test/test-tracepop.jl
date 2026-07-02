@testset "Trace Monoids            " begin
    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Projector.([a, b, c, d])
    build(M)

    TM = make_trace_monoid(M, 2)

    μ = TM.vertices_free
    ρ = TM.vertices_states
    w = μ[1]*ρ[2]*ρ[7]
    ρw = state_projection(w, TM)
    rev_dict_μ=Dict([j=>i for (i, j) in TM.dict_free])

    @test free_part(w, TM) == μ[1]
    @test state_part(w, TM) == ρ[2]*ρ[7]
    @test w == free_part(w, TM)*state_part(w, TM) == state_part(w, TM)*free_part(w, TM)
    @test ρw == TM.dict_states[monomial(rev_dict_μ[μ[1]])]*ρ[2]*ρ[7]
    @test word_to_state(ρw, TM) == state_projection(word_to_state(w, TM))

    f = μ[1] - TM.dict_states[monomial(rev_dict_μ[μ[1]])]
    g = μ[2]*μ[1] - TM.dict_states[rev_dict_μ[μ[2]] * rev_dict_μ[μ[1]]]

    @test state_projection(f, TM) == 0
    @test state_projection(g, TM) !== 0

    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Projector.([a, b, c, d])
    build(M)

    TM = make_trace_monoid(M, 2, tracial = true)

    μ = TM.vertices_free
    ρ = TM.vertices_states
    w = μ[1]*ρ[2]*ρ[7]
    ρw = state_projection(w, TM)
    rev_dict_μ=Dict([j=>i for (i, j) in TM.dict_free])
    m2=TM.dict_states[cyclic_reduce(monomial(rev_dict_μ[μ[1]]))]
    m6=TM.dict_states[cyclic_reduce(rev_dict_μ[μ[2]]*rev_dict_μ[μ[1]])]

    @test free_part(w, TM) == μ[1]
    @test state_part(w, TM) == ρ[2]*ρ[7]
    @test w == free_part(w, TM)*state_part(w, TM) == state_part(w, TM)*free_part(w, TM)
    @test ρw == m2*ρ[2]*ρ[7]
    @test word_to_state(ρw, TM) == state_projection(word_to_state(w, TM))

    @test state_projection(μ[1]*μ[2], TM) ==
          state_projection(μ[2]*μ[1], TM) ==
          TM.dict_states[cyclic_reduce(rev_dict_μ[μ[2]]*rev_dict_μ[μ[1]])]
    res=m2*m6
    @test state_projection(μ[1]*μ[2]*m2, TM) == state_projection(μ[2]*m2*μ[1], TM) == res

    f = μ[1] - m2
    g = μ[2]*μ[1] - m6

    @test state_projection(f, TM) == 0
    @test state_projection(g, TM) == 0
end

@testset "Trace Optimization       " begin
    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Projector.([a, b, c, d])
    build(M)

    f = a*c + a*d + b*c - b*d

    val, sos_model, _ = pcpop(f, 1; tracial = true)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)

    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Projector.([a, b, c, d])
    build(M)

    f = a*c + a*d + b*c - b*d

    val, sos_model, _ = pcpop(f, 2; tracial = true)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)
    @test abs(objective_value(sos_model) - 2) <= 1e-6
end

@testset "Statepop Optimization    " begin
    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Projector.([a, b, c, d])
    build(M)

    f = a*c + a*d + b*c - b*d

    sos_model = tpop(f, 1, 0)
    set_optimizer(sos_model, default_solver())
    set_silent(sos_model)
    optimize!(sos_model)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)

    f = a*c + a*d + b*c - b*d

    sos_model = tpop(f, 2, 1)
    set_optimizer(sos_model, default_solver())
    set_silent(sos_model)
    optimize!(sos_model)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)
end

@testset "Tracepop Optimization    " begin
    @pcmonoid M a b c d
    @comms [a, b] [c, d]
    Unipotent.([a, b, c, d])
    build(M)

    f = a*c + a*d + b*c - b*d

    sos_model = tpop(f, 1, 0, tracial = true)
    set_optimizer(sos_model, default_solver())
    set_silent(sos_model)
    optimize!(sos_model)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)

    f = a*c + a*d + b*c - b*d

    sos_model = tpop(f, 1, tracial = true)
    set_optimizer(sos_model, default_solver())
    set_silent(sos_model)
    optimize!(sos_model)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)

    f = a*c + a*d + b*c - b*d

    sos_model = tpop(f, 2, 1, tracial = true)
    set_optimizer(sos_model, default_solver())
    set_silent(sos_model)
    optimize!(sos_model)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)
end

@testset "Examples in Klep et al.  " begin
    # Example 6.1 f == 0
    @pcmonoid M a b c
    Projector(a)
    Projector(b)
    Projector(c)
    build(M)

    TM = make_trace_monoid(M, 6, tracial = true)

    f = - state(a*b*c, TM) - state(a*b, TM)*state(c, TM)

    basis_psd = trace_monomials(TM, 0:3, tracial = true)

    sos_model = tpop(f, TM, basis_psd, tracial = true)
    set_optimizer(sos_model, default_solver())
    set_silent(sos_model)
    optimize!(sos_model)
    @test termination_status(sos_model) == MOI.TerminationStatusCode(1)
    @test abs(-objective_value(sos_model) - (-1/32)) <= 1e-6
end

@testset "Statepop Klep et al.     " begin

    #Example 7.2.1 quadratic Bell inequality
    @pcmonoid M a[2, 0] b[2, 0]
    Unipotent.(a)
    Unipotent.(b)
    @comms a b
    build(M)

    TM = make_trace_monoid(M, 6, tracial = false)

    p = (state(a[1]*b[2], TM) + state(a[2]*b[1], TM))^2
    p += (state(a[1]*b[1], TM) - state(a[2]*b[2], TM))^2

    basis = trace_monomials(TM, 0:3, tracial = false)

    model = tpop(p, TM, basis, tracial = false)
    set_optimizer(model, default_solver())
    set_silent(model)
    optimize!(model)
    @test termination_status(model) == MOI.TerminationStatusCode(1)
    @test abs(objective_value(model) - 4) <= 1e-6
end
