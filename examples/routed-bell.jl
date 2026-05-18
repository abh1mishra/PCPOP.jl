using LinearAlgebra
include("../traceGrobner.jl")
include("/home/abhishek-mishra/Dropbox/repo/QCMaterial/SemiDIQKD/cgb_PnM/bb84_corr.jl")

function bb84_correlations(a, b, x, y)
    # Shared state |\phi^+> = (|00> + |11>) / sqrt(2)
    state = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    rho = state * state'
    Id = I(2)
    # Projectors for Z, X
    # x,y = 1 -> Z basis
    P_Z1 = [1.0 0.0; 0.0 0.0]  # Outcome 1 (+1)
    P_Z2 = [0.0 0.0; 0.0 1.0]  # Outcome 2 (-1)
    
    # x,y = 2 -> X basis
    P_X1 = [1.0 1.0; 1.0 1.0] / 2.0  # Outcome 1 (+1)
    P_X2 = [1.0 -1.0; -1.0 1.0] / 2.0 # Outcome 2 (-1)
    
    # Map index x,y -> Z or X; a,b -> 1 or 2
    proj = Dict(
        (1, 1) => P_Z1, (1, 2) => P_Z2,
        (2, 1) => P_X1, (2, 2) => P_X2
    )
    
    if a == -1
        P_B = proj[(y, b)]
        P_B = kron(Id, P_B)
        return (real(tr(rho * P_B)))
    end
    if b == -1
        P_A = proj[(x, a)]
        P_A = kron(P_A, Id)
        return (real(tr(rho * P_A)))
    end
    # Alice and Bob's projectors
    P_A = proj[(x, a)]
    P_B = proj[(y, b)]
    
    # Joint projector P_A ⊗ P_B
    P_AB = kron(P_A, P_B)

    # Probability = tr(rho * P_AB)
    return real(tr(rho * P_AB))
end

function chsh_correlations(a, b, x, y)
    # Shared state |\phi^+> = (|00> + |11>) / sqrt(2)
    state = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    rho = state * state'
    Id = I(2)

    # Alice's observables (Z, X)
    P_Z1 = [1.0 0.0; 0.0 0.0]  # x=1, a=1 (+1)
    P_Z2 = [0.0 0.0; 0.0 1.0]  # x=1, a=2 (-1)
    P_X1 = [1.0 1.0; 1.0 1.0] / 2.0  # x=2, a=1 (+1)
    P_X2 = [1.0 -1.0; -1.0 1.0] / 2.0 # x=2, a=2 (-1)

    proj_A = Dict(
        (1, 1) => P_Z1, (1, 2) => P_Z2,
        (2, 1) => P_X1, (2, 2) => P_X2
    )

    # Bob's observables: B1 = (Z+X)/sqrt(2), B2 = (Z-X)/sqrt(2)
    Z = [1.0 0.0; 0.0 -1.0]
    X = [0.0 1.0; 1.0 0.0]
    B1 = (Z + X) / sqrt(2)
    B2 = (Z - X) / sqrt(2)

    PB_1_plus = (Id + B1) / 2.0   # y=1, b=1 (+1)
    PB_1_minus = (Id - B1) / 2.0  # y=1, b=2 (-1)
    PB_2_plus = (Id + B2) / 2.0   # y=2, b=1 (+1)
    PB_2_minus = (Id - B2) / 2.0  # y=2, b=2 (-1)

    proj_B = Dict(
        (1, 1) => PB_1_plus, (1, 2) => PB_1_minus,
        (2, 1) => PB_2_plus, (2, 2) => PB_2_minus
    )

    if a == -1
        P_B = proj_B[(y, b)]
        P_B = kron(Id, P_B)
        return real(tr(rho * P_B))
    end
    if b == -1
        P_A = proj_A[(x, a)]
        P_A = kron(P_A, Id)
        return real(tr(rho * P_A))
    end

    # Joint projector P_A ⊗ P_B
    P_A = proj_A[(x, a)]
    P_B = proj_B[(y, b)]
    P_AB = kron(P_A, P_B)

    return real(tr(rho * P_AB))
end

####### Example 1 fig 5 in https://arxiv.org/abs/2310.07484
@pcmonoid M A[2,0] BS[2,0] BL[2,0]
Projector.(M.vertices)
@comms A BS
@comms A BL

# SRQ constraint, joint measuribility of BL
@comms BL

build(M)

UA = [a-(1-a) for a in A]
UBS = [b-(1-b) for b in BS]
UBL = [b-(1-b) for b in BL]

# C_s = 2\sqrt(2) constraint

Cs = UA[1]*UBS[1] + UA[1]*UBS[2] + UA[2]*UBS[1] - UA[2]*UBS[2]
tr_eq = [[Cs,2*sqrt(2)]]

θ = 0.0
level=2

obj = tan(θ)*UA[1]*UBL[1] + UA[1]*UBL[2] + UA[2]*UBL[1] - tan(θ)*UA[2]*UBL[2]

ov,_,_=npa(obj,level;tr_eq=tr_eq,min=false)
println("Optimal value is ", ov)

##### Example 2 fig 10 in https://arxiv.org/abs/2310.07484, UB on \eta_l in CHSH with inefficient detectors
@pcmonoid M A[4,0] BS[4,0] BL[4,0]
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

model, Γ, pm = npa(0,level;min=false,rm=true,list_vars = M.vertices)
ηl = @variable(model)
ηs = 0.98
ηa=ηs
[@constraint(model,Γ[PA[a,x]*PBS[b,y]] == ηa*ηs*chsh_correlations(a,b,x,y)) for a in 1:2 for x in 1:2 for b in 1:2 for y in 1:2]
id = one(A[1])
@constraint(model,Γ[PA[1,1]*id] == ηa*chsh_correlations(1,-1,1,1))
@constraint(model,Γ[PA[1,2]*id] == ηa*chsh_correlations(1,-1,2,1))
@constraint(model,Γ[PA[2,1]*id] == ηa*chsh_correlations(2,-1,1,1))
@constraint(model,Γ[PA[2,2]*id] == ηa*chsh_correlations(2,-1,2,1))

@constraint(model,Γ[PBS[1,1]*id] == ηs*chsh_correlations(-1,1,1,1))
@constraint(model,Γ[PBS[2,1]*id] == ηs*chsh_correlations(-1,2,1,1))
@constraint(model,Γ[PBS[1,2]*id] == ηs*chsh_correlations(-1,1,1,2))
@constraint(model,Γ[PBS[2,2]*id] == ηs*chsh_correlations(-1,2,1,2))


[@constraint(model,Γ[PA[a,x]*PBL[b,y]] == ηa*ηl*chsh_correlations(a,b,x,y)) for a in 1:2 for x in 1:2 for b in 1:2 for y in 1:2]

@constraint(model,Γ[PBL[1,1]*id] == ηl*chsh_correlations(-1,1,1,1))
@constraint(model,Γ[PBL[2,1]*id] == ηl*chsh_correlations(-1,2,1,1))
@constraint(model,Γ[PBL[1,2]*id] == ηl*chsh_correlations(-1,1,1,2))
@constraint(model,Γ[PBL[2,2]*id] == ηl*chsh_correlations(-1,2,1,2))

@objective(model, Max, ηl)
optimize!(model)
println("Optimal value is ", value(ηl))