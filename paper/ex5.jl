using cgb

println("==== Case n=2 ====")
function run_n2()
    @pcmonoid M A[2,0] B[2,0]
    Projector.([A;B])
    @comms A B
    build(M)
    obj = A[1]*B[1]
    
    # n=2 constraints analog to Eq 2 (standard Hardy)
    tr_eq = [
        [A[2] * B[1], 0],
        [A[1] * B[2], 0],
        [(1-A[2]) * (1-B[2]), 0]
    ]
    level = 2
    ov, model, var_dict, principal_matrix = npa(obj, level; tr_eq=tr_eq, min=false)
    println("Objective value (n=2): ", ov)
end
run_n2()

println("\n==== Case n=3 ====")
function run_n3(level)
    # Case n=3: Three parties A, B, C
    # Each party has two measurement settings: U (index 1) and D (index 2)
    # Both are projectors
    
    @pcmonoid M A[2,0] B[2,0] C[2,0]
    Projector.([A;B;C])

    # All parties are spatially separated, so all observables commute
    @comms A B
    @comms B C
    @comms A C
    build(M)

    # Objective: P(+1,+1,+1 | U1, U2, U3) = p > 0
    # In code: A[1]=U_A, B[1]=U_B, C[1]=U_C
    obj = A[1]*B[1]*C[1]

    # Constraints from Eq (2) of arXiv:2311.07266:
    # For i=1,2,...,n: P(+1,+1 | D_i, U_{i+1}) = 0 with n+1 ≡ 1
    # For n=3:
    #   i=1: P(+1,+1 | D_1, U_2) = 0  →  A[2] * B[1] = 0
    #   i=2: P(+1,+1 | D_2, U_3) = 0  →  B[2] * C[1] = 0
    #   i=3: P(+1,+1 | D_3, U_1) = 0  →  C[2] * A[1] = 0  (cyclic: 3+1 ≡ 1)
    # And: P(-1,-1,-1 | D_1, D_2, D_3) = 0  →  (1-A[2])*(1-B[2])*(1-C[2]) = 0
    tr_eq = [
        [A[2] * B[1], 0],                    # P(+1,+1 | D_1, U_2) = 0
        [B[2] * C[1], 0],                    # P(+1,+1 | D_2, U_3) = 0
        [C[2] * A[1], 0],                    # P(+1,+1 | D_3, U_1) = 0
        [(1-A[2]) * (1-B[2]) * (1-C[2]), 0]  # P(-1,-1,-1 | D_1, D_2, D_3) = 0
    ]

    ov, model, var_dict, principal_matrix = npa(obj, level; tr_eq=tr_eq, min=false)
    println("Objective value (n=3): ", ov)
end

println("\n==== Case n=4 ====")
function run_n4(level)
    # Case n=4: Four parties A, B, C, D
    # Each party has two measurement settings: U (index 1) and D (index 2)
    
    @pcmonoid M A[2,0] B[2,0] C[2,0] D[2,0]
    Projector.([A;B;C;D])

    # All parties are spatially separated, so all observables commute
    @comms A B
    @comms B C
    @comms C D
    @comms D A
    @comms A C
    @comms B D

    build(M)

    # Objective: P(+1,+1,+1,+1 | U1, U2, U3, U4) = p > 0
    obj = A[1]*B[1]*C[1]*D[1]

    # Constraints from Eq (2) of arXiv:2311.07266:
    # For i=1,2,3,4: P(+1,+1 | D_i, U_{i+1}) = 0 with 4+1 ≡ 1
    #   i=1: P(+1,+1 | D_1, U_2) = 0  →  A[2] * B[1] = 0
    #   i=2: P(+1,+1 | D_2, U_3) = 0  →  B[2] * C[1] = 0
    #   i=3: P(+1,+1 | D_3, U_4) = 0  →  C[2] * D[1] = 0
    #   i=4: P(+1,+1 | D_4, U_1) = 0  →  D[2] * A[1] = 0  (cyclic: 4+1 ≡ 1)
    # And: P(-1,-1,-1,-1 | D_1, D_2, D_3, D_4) = 0
    tr_eq = [
        [A[2] * B[1], 0],                              # P(+1,+1 | D_1, U_2) = 0
        [B[2] * C[1], 0],                              # P(+1,+1 | D_2, U_3) = 0
        [C[2] * D[1], 0],                              # P(+1,+1 | D_3, U_4) = 0
        [D[2] * A[1], 0],                              # P(+1,+1 | D_4, U_1) = 0
        [(1-A[2]) * (1-B[2]) * (1-C[2]) * (1-D[2]), 0] # P(-1,-1,-1,-1 | D's) = 0
    ]

    ov, model, var_dict, principal_matrix = npa(obj, level; tr_eq=tr_eq, min=false)
    println("Objective value (n=4): ", ov)
end
