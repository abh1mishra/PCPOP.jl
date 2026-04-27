using Graphs,Combinatorics,Base.Iterators,AbstractAlgebra , JuMP, Mosek,MosekTools, StatsBase,NautyGraphs
import MutableArithmetics as MA
include("src/GM.jl")
include("src/var.jl")
include("src/monoid.jl")
include("src/monomial.jl")
include("src/PCMonomial.jl")  # PCMonomial - fast implementation for GraphProductMonoid{Variable}
include("src/polynomial.jl")
include("src/NCMonoid.jl")
include("src/utils.jl")
include("src/representations.jl")
include("src/cyclic.jl")
include("src/primal_opt.jl")
include("src/pcpop.jl")
include("src/grobner.jl")
include("src/tpop.jl")

const package_types=Dict([NCMonoid=>NCWord,
GraphProductMonoid=>GraphProductWord,
Variable=>PCMonomial,
AbstractMonoid=>AbstractMonomial,
GraphProductMonoid{Variable}=>PCMonomial,
TraceMonoid=>TraceWord])