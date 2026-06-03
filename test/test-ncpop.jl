@testset "ncpop" begin
    # To define the non-commutative variables X1,X2 and Non-commutative monoid M
    @ncmonoid M X1 X2
    #  To define the equality constraints
    Projector(X1)
    build(M)

    # To define the objective function, inequality constraints and level of relaxation
    obj=X1*X2+X2*X1
    op_ge=[X2-X2^2+0.5]
    level=2

    # To solve the non-commutative polynomial optimization problem
    ov,model,_=pcpop!(obj,level;min=true,op_ge=op_ge)
    @test abs(ov+0.7499999767923403) < 1e-6
end