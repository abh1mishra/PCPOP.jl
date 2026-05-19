using LinearAlgebra
include("../traceGrobner.jl")
include("/home/abhishek-mishra/Dropbox/repo/QCMaterial/SemiDIQKD/cgb_PnM/bb84_corr.jl")

function bb84_correlations(a, b, x, y; obs=false)
    # Shared state |\phi^+> = (|00> + |11>) / sqrt(2)
    state = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    rho = state * state'
    Id = I(2)

    if obs
        Z = [1.0 0.0; 0.0 -1.0]
        X = [0.0 1.0; 1.0 0.0]
        obs_dict = Dict(1 => Z, 2 => X)
        
        if x != -1 && y != -1
            O_A = obs_dict[x]
            O_B = obs_dict[y]
            return real(tr(rho * kron(O_A, O_B)))
        elseif x != -1
            O_A = obs_dict[x]
            return real(tr(rho * kron(O_A, Id)))
        elseif y != -1
            O_B = obs_dict[y]
            return real(tr(rho * kron(Id, O_B)))
        else
            return 1.0 # tr(rho * Id \otimes Id)
        end
    end

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

function chsh_correlations(a, b, x, y; obs=false)
    # Shared state |\phi^+> = (|00> + |11>) / sqrt(2)
    state = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    rho = state * state'
    Id = I(2)
    Z = [1.0 0.0; 0.0 -1.0]
    X = [0.0 1.0; 1.0 0.0]

    if obs
        # Alice's observables (Z, X)
        obs_A = Dict(1 => Z, 2 => X)
        
        # Bob's observables: B1 = (Z+X)/sqrt(2), B2 = (Z-X)/sqrt(2)
        B1 = (Z + X) / sqrt(2)
        B2 = (Z - X) / sqrt(2)
        obs_B = Dict(1 => B1, 2 => B2)

        if x != -1 && y != -1
            O_A = obs_A[x]
            O_B = obs_B[y]
            return real(tr(rho * kron(O_A, O_B)))
        elseif x != -1
            O_A = obs_A[x]
            return real(tr(rho * kron(O_A, Id)))
        elseif y != -1
            O_B = obs_B[y]
            return real(tr(rho * kron(Id, O_B)))
        else
            return 1.0
        end
    end

    # Alice's observables (Z, X)
    P_Z1 = [1.0 0.0; 0.0 0.0]  # x=1, a=1 (+1)
    P_Z2 = [0.0 0.0; 0.0 1.0]  # x=1, a=2 (-1)
    P_X1 = [1.0 1.0; 1.0 1.0] / 2.0  # x=2, a=1 (+1)
    P_X2 = [1.0 -1.0; -1.0 1.0] / 2.0 # x=2, a=1 (+1)

    proj_A = Dict(
        (1, 1) => P_Z1, (1, 2) => P_Z2,
        (2, 1) => P_X1, (2, 2) => P_X2
    )

    # Bob's observables: B1 = (Z+X)/sqrt(2), B2 = (Z-X)/sqrt(2)
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
@pcmonoid M UA[2,0] UBS[2,0] UBL[2,0]
Unipotent.(M.vertices)
@comms UA UBS
@comms UA UBL

# SRQ constraint, joint measuribility of BL
@comms UBL

build(M)

# C_s = 2\sqrt(2) constraint

Cs = UA[1]*UBS[1] + UA[1]*UBS[2] + UA[2]*UBS[1] - UA[2]*UBS[2]
tr_eq = [[Cs,2*sqrt(2)]]

θ = pi/4
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
[@constraint(model,Γ[PA[a,x]*id] == ηa*chsh_correlations(a,-1,x,1)) for a in 1:2 for x in 1:2]

[@constraint(model,Γ[PBS[b,y]*id] == ηs*chsh_correlations(-1,b,1,y)) for b in 1:2 for y in 1:2]

[@constraint(model,Γ[PA[a,x]*PBL[b,y]] == ηa*ηl*chsh_correlations(a,b,x,y)) for a in 1:2 for x in 1:2 for b in 1:2 for y in 1:2]

[@constraint(model,Γ[PBL[b,y]*id] == ηl*chsh_correlations(-1,b,1,y)) for b in 1:2 for y in 1:2]

@objective(model, Max, ηl)
optimize!(model)
println("Optimal value is ", value(ηl))


##### Example 3 using trichotomic observables fig 10 in https://arxiv.org/abs/2310.07484, UB on \eta_l in CHSH with inefficient detectors
@subsystem M A[2,0] B[4,0]
add_relations!([i^3-i for i in a])
add_relations!([i^3-i for i in b])
bs = b[1:2]
bl = b[3:4]

# SRQ constraint, joint measuribility of BL
add_relations!([bl[1]*bl[2]-bl[2]*bl[1]])

build(M)

level = "3+a*a*a*a+b[1:2]*b[1:2]*b[1:2]*b[1:2]+b[3:4]*b[3:4]*b[3:4]*b[3:4]+a*a*b[1:2]*b[1:2]+a*a*b[3:4]*b[3:4]+b[1:2]*b[1:2]*b[3:4]*b[3:4]"
model, Γ, pm = npa(0,level;min=false,rm=true,list_vars = [a;b])
ηl = @variable(model)
ηs = 0.98
ηa=ηs
# AxBsy
[@constraint(model,Γ[a[x]*bs[y]] == ηa*ηs*chsh_correlations(1,1,x,y;obs=true)) for x in 1:2 for y in 1:2]#4

#Ax^2Bsy
[@constraint(model,Γ[a[x]^2*bs[y]] == ηa*ηs*chsh_correlations(1,1,-1,y;obs=true)) for x in 1:2 for y in 1:2]#4
#AxBsy^2
[@constraint(model,Γ[a[x]*bs[y]^2] == ηa*ηs*chsh_correlations(1,1,x,-1;obs=true)) for x in 1:2 for y in 1:2]#4
#Ax^2Bsy^2
[@constraint(model,Γ[a[x]^2*bs[y]^2] == ηa*ηs) for x in 1:2 for y in 1:2]#4

id = monomial(one(monomial(a[1])))

#Ax
[@constraint(model,Γ[a[x]*id] == ηa*chsh_correlations(1,-1,x,-1;obs=true)) for x in 1:2]#2
#Ax^2
[@constraint(model,Γ[a[x]^2*id] == ηa) for x in 1:2]#2

#Bsy
[@constraint(model,Γ[bs[y]*id] == ηs*chsh_correlations(-1,-1,-1,y;obs=true)) for y in 1:2]#2
#Bsy^2
[@constraint(model,Γ[bs[y]^2*id] == ηs) for y in 1:2]#2

#AxBly
[@constraint(model,Γ[a[x]*bl[y]] == ηa*ηl*chsh_correlations(1,1,x,y;obs=true)) for x in 1:2 for y in 1:2]#4

#Ax^2Bly
[@constraint(model,Γ[a[x]^2*bl[y]] == ηa*ηl*chsh_correlations(1,1,-1,y;obs=true)) for x in 1:2 for y in 1:2]#4
#AxBly^2
[@constraint(model,Γ[a[x]*bl[y]^2] == ηa*ηl*chsh_correlations(1,1,x,-1;obs=true)) for x in 1:2 for y in 1:2]#4
#Ax^2Bly^2
[@constraint(model,Γ[a[x]^2*bl[y]^2] == ηa*ηl) for x in 1:2 for y in 1:2]#4

#Bly
[@constraint(model,Γ[bl[y]*id] == ηl*chsh_correlations(-1,-1,-1,y;obs=true)) for y in 1:2]#2
#Bly^2
[@constraint(model,Γ[bl[y]^2*id] == ηl) for y in 1:2]#2

@objective(model, Max, ηl)
optimize!(model)
println("Optimal value is ", value(ηl))