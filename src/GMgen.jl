# ============================================================================
# PCMonoid Generator - Creates nested GraphProductMonoid{NCMonoid} from GraphProductMonoid{Variable}
# ============================================================================

"""
    PCMonoid_generator(GM::GraphProductMonoid)

Generate a nested GraphProductMonoid{NCMonoid} from a GraphProductMonoid{Variable}.
Groups variables by their clique membership patterns and creates NCMonoid subsystems.
"""
function PCMonoid_generator(GM::GraphProductMonoid)
    vars = GM.vertices
    comms = GM.commutations
    name = GM.name
    graph = SimpleGraph(length(vars))

    for (i, j) in comms
        add_edge!(graph, Edge(findfirst(isequal(i), vars), findfirst(isequal(j), vars)))
    end
    graph = complement(graph)
    cliques = maximal_cliques(graph)
    cliques = [[vars[index] for index in clique] for clique in maximal_cliques(graph)]
    types = Dict{Vector{Vector{Variable}}, Vector{Variable}}()
    for var in vars
        k = filter(i -> var in i, cliques)
        if haskey(types, k)
            push!(types[k], var)
        else
            types[k] = [var]
        end
    end
    monoids = Vector{NCMonoid}([])
    for (index, (key, value)) in enumerate(types)
        ncmon = NCMonoid("$name$index", value)
        push!(monoids, ncmon)
    end
    
    monoid_comms = Vector{Tuple{NCMonoid, NCMonoid}}()
    for (i, j) in combinations(monoids, 2)
        edges = map(((l, r),) -> (findfirst(isequal(l), vars), findfirst(isequal(r), vars)), product(i.vertices, j.vertices))
        any(k -> has_edge(graph, k), edges) || push!(monoid_comms, (i, j))
    end
    res = GraphProductMonoid{NCMonoid}(name, monoids; commutations=monoid_comms)
    build(res)
    return res
end