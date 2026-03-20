include("./traceGrobner.jl")
include("./src/pcpop.jl")

function index_moments(basis)
    basis_constraints = StarAlgebras.Basis{UInt16}(
              unique([x*y for x in basis for y in basis]),
    )
    M = [basis_constraints[x*y] for x in basis, y in basis]
    return M, basis_constraints
end

###########
# Setting #
###########

cbA=cliqueBase()

@hermitian cbA AL[1:2] AS[1:2] BL[1:2] BS[1:2]

@comms [AL; AS] [BL; BS]
@comms AL
@projector cbA.v

build(cbA)

f = +(cbA.v...)*(1/2) + 1

k = 2
basisA = monomials([AS; AL; BS; BL], k)

momentsA, basis_momentsA = index_moments(basisA)

model = JuMP.Model()
JuMP.@variable model η
JuMP.@objective model Max η

n = length(basisA)
P = JuMP.@variable model P[1:n, 1:n] Symmetric
JuMP.@constraint model P in PSDCone()

objective = η*f
for (idx, b) in enumerate(basis_momentsA)
    c = MP.coefficient(objective, b)
    JuMP.@constraint model LinearAlgebra.dot(P, momentsA .== idx) == c
end

# commutations, moments
set_optimizer(model, Mosek.Optimizer)
optimize!(model)
@show objective_value(model)