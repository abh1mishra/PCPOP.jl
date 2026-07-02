using PCPOP
using Test
using JuMP
using Mosek, MosekTools
using SymbolicWedderburn
using SymbolicWedderburn: StarAlgebras

# The test suite also exercises internal (non-exported) functions:
using PCPOP:
    OnLetters,
    PG,
    Polynomial,
    add_comms!,
    coefficient,
    divide,
    free_part,
    make_trace_monoid,
    monomial,
    monomial_to_foata,
    state,
    state_part,
    state_projection,
    tpop,
    trace_monomials,
    word_to_state

include("test-ncpop.jl")
include("test-cyclic.jl")
include("test-monomials.jl")
include("test-optimization.jl")
include("test-npa.jl")
include("test-tracepop.jl")
