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
    ov,_=pcpop(wit,level;op_ge=ge,tr_eq=eq_constr,tr_ge=ge_constr,tracial=true,min=false,normalize=false,block_diag=diagonalize,reduce=true)
    println("Optimal value is ", ov)
end

@testset "routed bell" begin
    diagonalize = true # whether to diagonalize the moment matrix or not, this is a parameter of the jordan reduction method and it will not change the values of the plot (they use diagonalization in their code)
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
    ov, model, _=pcpop(obj, level; tr_eq = tr_eq, min = false,reduce=true,block_diag=diagonalize) 
end

