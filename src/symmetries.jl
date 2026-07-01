using SymbolicWedderburn
using AbstractPermutations
import PermutationGroups   # only used via the `PG` alias below; `import` avoids a `degree` name clash
using GroupsCore

using SymbolicWedderburn.StarAlgebras

const PG = PermutationGroups
const AP = AbstractPermutations

"""
    Structure OnLetters <: SymbolicWedderburn.ByPermutations.
    Permutes letters in `word` with permutation `perm`.

    #Example: (12)(abca) = bacb.
"""
struct OnLetters <: SymbolicWedderburn.ByPermutations end
function SymbolicWedderburn.action(
    ::OnLetters,
    perm::PG.Perms.AbstractPermutation,
    # perm::AbstractPermutations.AbstractPermutation,
    word::Variable,
)

    vertices = monomial(word).monoid.vertices
    return vertices[findfirst(isequal(word), vertices)^perm]
end

function SymbolicWedderburn.action(
    ::OnLetters,
    perm::PG.Perms.AbstractPermutation,
    # perm::AbstractPermutations.AbstractPermutation,
    word::Vector{Variable},
)
    return [SymbolicWedderburn.action(OnLetters(), perm, letter) for letter in word]
end

function SymbolicWedderburn.action(
    ::OnLetters,
    perm::PG.Perms.AbstractPermutation,
    # perm::AbstractPermutations.AbstractPermutation,
    word::AbstractMonomial,
)
    if word == one(word)
        return word
    else
        ncword = SymbolicWedderburn.action(OnLetters(), perm, monomial_to_word(word))
        return words_to_monomial(word.monoid, ncword)
    end
end

function SymbolicWedderburn.action(
    ::OnLetters,
    perm::PG.Perms.AbstractPermutation,
    # perm::AbstractPermutations.AbstractPermutation,
    word::Polynomial,
)
    return sum([coef*SymbolicWedderburn.action(OnLetters(), perm, mono) for (coef, mono) in zip(word.coeffs, word.monomials)])
end