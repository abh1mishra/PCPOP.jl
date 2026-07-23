@testset "Standard NPA             " begin
    @pcmonoid M a[2, 0] b[2, 0]
    @comms a b
    Unipotent.([a; b])
    build(M)
    chsh=a[1]*(b[1]+b[2])+a[2]*(b[1]-b[2])
    level=2
    optimal_value, _, _=pcpop(chsh, level; min = false)
    @test abs(optimal_value-2*sqrt(2)) < 1e-6
end

@testset "Standard NPA Cyclic      " begin
    @pcmonoid M a b c d e
    @comms [a, b] [c, d]
    Unipotent.([a, b, c, d])
    Projector(e)
    build(M)
    chsh = a*(c+d) + b*(c-d)
    chsh=e*chsh
    level=2
    tr_eq=[(e, 1.0)]
    optimal_value, _, _=pcpop(
        chsh,
        level;
        min = false,
        tr_eq = tr_eq,
        tracial = true,
        normalize = false,
        primal = false,
    )
    @test abs(optimal_value-2*sqrt(2)) < 1e-6
end
