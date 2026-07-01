module PCPOP

export @pcmonoid, @ncmonoid, @subsystem, @comms,@ortho, @pcmonoid_simple, GraphProductMonoid, npa, Projector, Unipotent, Unitary,
add_relations!, add_mult!,build,monomial,Polynomial,
get_monomials,mons_at_level, npa_moments_block, cyclic_npa_moments_block,
cyclic_reduce,real_rep,variables,pcpop!, pcpop, tpop, obj_gen, cons_gen!,model_new_obj, npa_dual, basis_gen,make_trace_monoid,state,trace_monomials,
default_solver, mosek_available

using Graphs,Combinatorics,Base.Iterators,
AbstractAlgebra ,LinearAlgebra, JuMP,SparseArrays, Mosek, MosekTools, Clarabel, StatsBase, ClusteredLowRankSolver
import MutableArithmetics as MA
include("GM.jl")
include("var.jl")
include("monoid.jl")
include("monomial.jl")
include("PCMonomial.jl")  # PCMonomial - fast implementation for GraphProductMonoid{Variable}
include("polynomial.jl")
include("NCMonoid.jl")
include("utils.jl")
include("representations.jl")
include("cyclic.jl")
include("npa_level.jl")
include("primal_opt.jl")
include("primal_opt_nc.jl")
include("dual_opt.jl")
include("jordan.jl")
include("symmetries.jl")
include("grobner.jl")
include("pcpop_sos.jl")
include("tpop.jl")
include("optimization.jl")

const package_types=Dict([NCMonoid=>NCWord,
GraphProductMonoid=>GraphProductWord,
Variable=>PCMonomial,
AbstractMonoid=>AbstractMonomial,
GraphProductMonoid{Variable}=>PCMonomial,
TraceMonoid=>TraceWord])


end
