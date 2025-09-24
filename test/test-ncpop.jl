@testset "ncpop" begin
    @ncmonoid M x1 x2
    Projector(x1)
    build(M)

    obj=x1*x2+x2*x1
    op_ge=[x2-x2^2+0.5]
    level=2
    @test abs(npa(obj,level;op_ge=op_ge)[1]+0.7499999767923403) < 1e-6
end