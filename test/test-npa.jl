@testset "Standard NPA" begin
    @subsystem M A[2,0] B[2,0]
    [add_mult!(i,i,1) for i in [a;b]]
    build(M)
    chsh=a[1]*(b[1]+b[2])+a[2]*(b[1]-b[2])
    level=2
    optimal_value,_,_=npa(chsh,level;min=false)
    @test abs(optimal_value-2*sqrt(2)) < 1e-6
end
@testset "Standard NPA with Cyclic Reduction" begin
    @pcmonoid M a b c d e
    @comms [a,b] [c,d]
    Unipotent.([a,b,c,d])
    Projector(e)
    build(M)
    chsh = a*(c+d)+ b*(c-d)
    chsh=e*chsh
    level=2
    tr_eq=[(e,1.0)]
    optimal_value,_,_=npa(chsh,level;min=false,tr_eq=tr_eq,cyclic=true,normalize=false)
    @test abs(optimal_value-2*sqrt(2)) < 1e-6
end

@testset "Standard NPA with Cyclic Reduction and Unipotent" begin
    G=0.8  # value of the guessing probability
    @pcmonoid M a b c d e f
    ρ=[b,c,d]
    σ=a
    PA=[e,f]
    Projector.([e,f])
    build(M)
    PA=[e f;1-e 1-f]
    ge = [σ-1/3*ρ[i] for i in 1:3]
    ge=vcat(ge,[i-i*i for i in ρ])
    ge_constr = [[-σ,-G]]
    eq_constr = [[ρ[x],1] for x in 1:3]
    wit = -(2*PA[1,1]-1)*ρ[1]-(2*PA[1,2]-1)*ρ[1]-(2*PA[1,1]-1)*ρ[2]+(2*PA[1,2]-1)*ρ[2]+(2*PA[1,1]-1)*ρ[3]
    level = 2 # this is the level of the localising matrices and it will be enough to retrieve the values of the plot (they use level 3 of the principal moment matrix)
    ov,_,_=npa(wit,level;op_ge=ge,tr_eq=eq_constr,tr_ge=ge_constr,cyclic=true,min=false,normalize=false)
    @test abs(ov-4.26552739737409) < 1e-4
end