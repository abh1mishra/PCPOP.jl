include("traceGrobner.jl")
include("src/pcpop.jl")

###########
# Setting #
###########

cb=cliqueBase()

@hermitian cb a b c d

@comms [a, b] [c, d]
@projector [a,b,c,d]

build(cb)

f = a*c + a*d + b*c - b*d
    
#######################
# SYMBOLIC WEDDERBURN #
#######################

# Moment relaxation
k = 2
monomial_basis = monomials(cb, k)
basis_constraints = StarAlgebras.Basis{UInt16}(
              unique([x*y for x in monomial_basis for y in monomial_basis]),
    )

# Symmetries
action = OnLetters()

g = PG.perm"(1,2)" #"
ag = SymbolicWedderburn.action(action, g, a)

p1 = PG.perm"(1,3)(2,4)" #"
p2 = PG.perm"(1,4)(2,3)" #"
G = PG.PermGroup(p1, p2)
T = Float64

tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
ehom = SymbolicWedderburn.CachedExtensionHomomorphism(
        G,
        action,
        monomial_basis
    )

ψ = SymbolicWedderburn.action_character(ehom, tbl)

inv = invariant_vectors(tbl, action, basis_constraints)

sa_basis = SymbolicWedderburn.symmetry_adapted_basis(
            T,
            G,
            action,
            monomial_basis,
            semisimple = false,
        )

ss_basis = SymbolicWedderburn.symmetry_adapted_basis(
            T,
            G,
            action,
            monomial_basis,
            semisimple = true,
        )

wedderburn = SymbolicWedderburn.WedderburnDecomposition(
            T,
            G,
            action,
            basis_constraints,
            monomial_basis,
            semisimple=true,
        )

##########
# PC-POP #
##########

# SOS relaxation
k = 2

# Unsymmetrized
sos_model = pcpop(f, k)
set_optimizer(sos_model, Mosek.Optimizer)
set_silent(sos_model)
optimize!(sos_model)
@show termination_status(sos_model)
@show objective_value(sos_model)

# Symmetrized
sa_model = pcpop(f, k, G, action)
set_optimizer(sa_model, Mosek.Optimizer)
set_silent(sa_model)
optimize!(sa_model)
@show termination_status(sa_model)
@show objective_value(sa_model)

# Block diagonal
blk_model = pcpop(f, k, G, action, diagonalize=true)
set_optimizer(blk_model, Mosek.Optimizer)
set_silent(blk_model)
optimize!(blk_model)
@show termination_status(blk_model)
@show objective_value(blk_model)