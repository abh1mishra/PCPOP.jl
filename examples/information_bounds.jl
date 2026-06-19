include("traceGrobner.jl")
# PROJECTORS

# Parameters
G = 0.8

# Build the monoid
@pcmonoid M ρ[3,0] σ A[2,0]
Projector.(A)
build(M)
PA=[A[1] A[2];1-A[1] 1-A[2]]

# Conditions on the operators
op_ge = vcat([σ-(1/3)*r for r in ρ], [r-r^2 for r in ρ]) 

# Conditions on the moments
tr_ge = [[-σ,-G]]
tr_eq = [[ρ[x],1] for x in 1:3]

# Optimization of the semidefinite relaxation
obj = -(2*PA[1,1]-1)*ρ[1] - (2*PA[1,2]-1)*ρ[1] - (2*PA[1,1]-1)*ρ[2]
obj+=  (2*PA[1,2]-1)*ρ[2] + (2*PA[1,1]-1)*ρ[3]
val, model, _ = pcpop!(obj, 2; min=false,
                 op_ge=op_ge,
                 tr_eq=tr_eq,
                 tr_ge=tr_ge,
                 tracial=true,
                 primal=false,
                 normalize=false)
println("Optimal value is ", val)

# UNITARIES

G = 0.8
# Build the monoid
@pcmonoid M ρ[3,0] σ a[2,0]
Unipotent.(a)
build(M)
# Conditions on the operators
op_ge = vcat([σ-(1 /3)*r for r in ρ], [r-r^2 for r in ρ]) 
# Conditions on the moments
tr_ge = [[-σ,-G]]
tr_eq = [[ρ[x],1] for x in 1:3]
# Optimization of the semidefinite relaxation
obj = -a[1]*ρ[1] - a[2]*ρ[1] - a[1]*ρ[2] +  a[2]*ρ[2] + a[1]*ρ[3]
val, model, _ = pcpop!(obj, 2; min=false,
                 op_ge=op_ge,
                 tr_eq=tr_eq,
                 tr_ge=tr_ge,
                 tracial=true,
                 normalize=false)
println("Optimal value is ", val)