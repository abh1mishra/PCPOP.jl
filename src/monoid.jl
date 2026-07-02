"""
    GraphProductMonoid

    #Notes
    A mutable struct representing a partially commutative monoid `M`.
    The `dependency graph` describes the commutation relations among variables, where connected variables don't commute and disconnected variables do commute. For the cliques, we use the maximnal cliques for now.
    
    #Fields:
    - `variables` : variables in `M`.
    - `num_variables` : number of variables.
    - `commutations` : pairs of commuting variables.
    - `cliques` : maximal cliques (totally connected sets) of dependency graph.
    - `mult_rules` : dictionary (a, b) => c where a*b = c.
    - `projectors` : dictionary (a, a) => a.
    - `orthogonal_pairs` : set of unordered pairs {a,b} such a*b = 0.
    - `is_built` : boolean indicating if the monoid is built.

    The monoid `M` should be built by calling `build(M)` before using it for computations.

    Example usage:
    ```julia
    GraphProductMonoid("M", [A, B, C, BC])
    @comms A B C
    ```
    Alternate usage:
    ```julia
    @pcmonoid M A B C BC
    @comms A B C
    ```
"""
struct GraphProductMonoid{T <: Union{AbstractMonoid, Variable}} <: AbstractMonoid
    # mandatory fields
    name::String
    parent_monoid::Base.RefValue{AbstractMonoid}
    is_built::Base.RefValue{Bool}
    conj_type::Base.RefValue{Bool}

    # technical fields
    vertices::Vector{T}
    commutations::Vector{Tuple{T, T}}
    cliques::Vector{Vector{T}}

    # To manipulate current monoid within parent monoid
    commutes_with::Vector{AbstractMonoid}
    clique_indices::Vector{Int}

    # Indicates if the monoid has implicit algebraic relations (e.g., a^2 = 1, a*b = 0, etc.)
    has_relations::Base.RefValue{Bool}

    # Cached index mappings for GraphProductMonoid{Variable} (PCMonomial optimization)
    # For non-Variable monoids, these remain empty/nothing
    var_index_dict::Base.RefValue{Union{Nothing, Dict{Variable, UInt32}}}
    conj_indices::Base.RefValue{Union{Nothing, Vector{UInt32}}}

    #  show fields, show_level can be 1,2,3 to control the amount of information printed when show is called on the monoid
    show_level::Base.RefValue{Int}

    function GraphProductMonoid{T}(
        name::String,
        vertices::Vector{T};
        parent_monoid::Base.RefValue{AbstractMonoid} = Base.RefValue{AbstractMonoid}(),
        commutations::Vector{Tuple{T, T}} = Vector{Tuple{T, T}}([]),
        empty_cliques = Vector{Vector{T}}([]),
        commutes_with = Vector{T}([]),
        clique_indices = Vector{Int}([]),
    ) where {T <: Union{AbstractMonoid, Variable}}
        res=new{T}(
            name,
            parent_monoid,
            Base.RefValue(false),
            Base.RefValue(false),
            vertices,
            commutations,
            empty_cliques,
            commutes_with,
            clique_indices,
            Base.RefValue(false),
            Base.RefValue{Union{Nothing, Dict{Variable, UInt32}}}(nothing),
            Base.RefValue{Union{Nothing, Vector{UInt32}}}(nothing),
            Base.RefValue{Int}(1),
        )
        for i in vertices
            i.parent_monoid[]=res
        end
        return res
    end
end
GraphProductMonoid(
    name::String,
    vertices::Vector{T},
) where {T <: Union{AbstractMonoid, Variable}} = GraphProductMonoid{T}(name, vertices)

function Base.show(io::IO, mime::MIME"text/plain", monoid::GraphProductMonoid)
    print(
        io,
        "THE MONOID IS $(monoid.is_built[] ? "BUILT" : "NOT YET BUILT")\n 
*** Variables *** \n# of variables = $(length(monoid.vertices))\nvariables = $(vars_to_str(monoid.vertices))\n
*** Commutation structure *** \ncommutations = ",
    )
    display(mime, monoid.commutations)
    print(io, "\nmaximal cliques = $(vars_to_str.(monoid.cliques))")
end

# A GraphProductMonoid is identidfied by its name and type of representation i.e graph or clique normal form
Base.:(==)(M::GraphProductMonoid, N::GraphProductMonoid) =
    (M.name==N.name) && (M.vertices==N.vertices) && (M.commutations==N.commutations)

"""
    add_elems!(M::GraphProductMonoid, m_arr::Vector{m}) where {m <: Union{AbstractMonoid, Variable}}

    # Notes
    Adds elements to the monoid `M` and updates their parent monoid reference.

    # Arguments
    - `M`: The GraphProductMonoid to which elements are added.
    - `m_arr`: A vector of elements (either AbstractMonoid or Variable) to be added.

    # Returns
    - None. The function modifies the `M` in-place.
"""
function add_elems!(
    M::GraphProductMonoid,
    m_arr::Vector{m},
) where {m <: Union{AbstractMonoid, Variable}}
    push!(M.vertices, m_arr...)
    for i in m_arr
        i.parent_monoid[]=M
    end
end

"""
    add_comms!(comms)

    # Notes
    Adds commutations to a given set.

    # Arguments
    - `comms`: A set of commutations. It can be a vector of Variables or a vector of vectors.

    # Returns
    - None. The function modifies the `comms` set in-place.
"""
function add_comms!(comms)
    (length(comms)==1 && isa(comms[1], Vector)) ? comms=comms[1] : nothing
    length(comms)==1 &&
        isa(comms[1], Union{Variable, AbstractMonoid}) &&
        throw("Provide commutations")

    z=[(isa(i, Union{Variable, AbstractMonoid})) ? [i] : i for i in comms]
    z=collect(combinations(z, 2))
    z=[vec(collect(product(i...))) for i in z]
    z=union(z...)
    add_comms!(z)
end

"""
    macro comms(monoids...)

    #Notes
    A macro to add commutations between monoids.

    # Arguments
    - `monoids...`: A variable number of monoid arguments.

    # Returns
    - None. The macro modifies the commutations in-place.

    Example:
    ```julia
    @pcmonoid M1 a b c
    @pcmonoid M2 x y z
    @pcmonoid M3 p q r
    @comms M1 M2 M3
    ```
    This will add commutations between M1, M2, and M3.
"""
macro comms(monoids...)
    return :($add_comms!($(Expr(:tuple, esc.(monoids)...))))
end

"""
    add_comms!(comms)

    # Notes
    Adds commutations to a given set.

    # Arguments
    - `comms`: A set of commutations. It can be a vector of Variables or a vector of vectors.

    # Returns
    - None. The function modifies the `comms` set in-place.
"""
function add_comms!(
    comms::Vector{Tuple{m, n}},
) where {m <: Union{Variable, AbstractMonoid}, n <: Union{Variable, AbstractMonoid}}
    M=isempty(comms) ? throw("Provide commutations") : first(first(comms)).parent_monoid[]
    M.is_built[] && throw("Operators already built")
    if m==Variable
        extra_comms=Vector{Tuple{m, m}}([])
        for (i, j) in comms
            (!is_herm(i)) && push!(extra_comms, (i', j))
            (!is_herm(j)) && push!(extra_comms, (i, j'))
            (!is_herm(i) && !is_herm(j)) && push!(extra_comms, (i', j'))
        end
        comms=unique(union(comms, extra_comms))
    end
    push!(M.commutations, comms...)
end

function add_orthos!(orthos)
    (length(orthos)==1 && isa(orthos[1], Vector)) ? orthos=orthos[1] : nothing
    length(orthos)==1 && isa(orthos[1], Variable) && throw("Provide orthogonal pairs")

    z=[(isa(i, Variable)) ? [i] : i for i in orthos]
    z=collect(combinations(z, 2))
    z=[vec(collect(product(i...))) for i in z]
    z=union(z...)
    add_orthos!(z)
end

"""
    macro ortho(vars...)

    #Notes
    A macro to add orthogonal pairs of variables. Works the same way as `@comms` macro, but for orthogonal pairs.

    # Arguments
    - `vars...`: A variable number of variable arguments.

    # Returns
    - None. The macro modifies the orthogonal pairs in-place.

    Example:
    ```julia
    @pcmonoid M a b c
    @ortho a b
    ```
    This will add the orthogonal pair (a, b) to the monoid M, so a*b = 0. 

    ```julia
    @pcmonoid M A[2,0] B[2,0] C[2,0]
    @ortho A[2,0] B[2,0]
    ```
    This will add the orthogonal pairs for each a in A and b in B to the monoid M.
"""
macro ortho(vars...)
    return :($add_orthos!($(Expr(:tuple, esc.(vars)...))))
end

function add_orthos!(orthos::Vector{Tuple{m, m}}) where {m <: Variable}
    M=isempty(orthos) ? throw("Provide orthogonal pairs") :
      first(first(orthos)).parent_monoid[]
    M.is_built[] && throw("Operators already built")

    for (i, j) in orthos
        if i.parent_monoid[] != j.parent_monoid[]
            throw("Variables belong to different monoids")
        end
        if typeof(i.parent_monoid[]) == NCMonoid
            add_relations!([i*j, j'*i'])
        else
            push!(i.ortho_conj, j)
            push!(j'.ortho_conj, i')
        end
    end
end

"""
    add_comms!(monoids::Vector{T}, comms::Matrix{Int64}) where T <: Union{Variable, AbstractMonoid}

    # Notes
    Adds commutations between monoids or variables based on a matrix of commutation relations.

    # Arguments
    - `monoids`: A vector of Variable objects.
    - `comms`: A matrix of integers representing commutation relations. If `comms[i, j]` is 1, a commutation is added between `monoids[i]` and `monoids[j]`.

    # Returns
    - None. The function modifies the `monoids` vector in-place.

"""
function add_comms!(
    elems::Vector{T},
    comms::Matrix{Int64},
) where {T <: Union{Variable, AbstractMonoid}}
    # Check for dimension mismatch between the number of elements and the size of the commutation matrix
    ((n=length(elems))!=size(comms)[1]) ? throw("Dimension mismatch") : nothing

    # Creates a vector of tuples representing the pairs of elements that commute based on the commutation matrix
    z=Vector{Tuple{T, T}}()

    # Iterate through the upper triangle of the commutation matrix to find pairs of elements that commute
    for i in 1:n
        for j in i:n
            if comms[i, j]==1
                push!(z, (elems[i], elems[j]))
            end
        end
    end

    # Add the commutations to the monoid using the add_comms! function
    add_comms!(z)
end

"""
    build(M::GraphProductMonoid{V})

    # Notes
    Builds the monoid `M` and sets `is_built` as true. If "M" is a GraphProduct Monoid of other monoids, it builds the child monoids as well.

    # Arguments
    - `M`: The GraphProductMonoid to be built.
    # Returns
    - None. The function modifies the `M` in-place. 

    Example usage:
    ```julia
    @pcmonoid M a b c
    @comms a b c
    build(M)
    ```
    WARNING: It is not possible to change the fields of a built monoid.
"""
function build(parent_monoid::GraphProductMonoid)

    # checking irregularities in input phase
    parent_monoid.is_built[] && throw("Operators already built")
    num_monoids=length(unique!(parent_monoid.vertices))
    num_monoids==0 && throw("No variables")

    # Ignores duplicate commutations
    unique!(parent_monoid.commutations)

    # building the independency graph
    graph=SimpleGraph(num_monoids)
    for (i, j) in parent_monoid.commutations
        add_edge!(
            graph,
            Edge(
                findfirst(isequal(i), parent_monoid.vertices),
                findfirst(isequal(j), parent_monoid.vertices),
            ),
        )
    end

    # taking the complement of the graph to get the dependency graph
    graph=complement(graph)

    # finding the maximal cliques of the dependency graph
    cliques=sort([
        [parent_monoid.vertices[j] for j in clique] for clique in maximal_cliques(graph)
    ])

    #  updating the variables with the clique_indices they belong to and the variables they commute with
    for child_monoid in parent_monoid.vertices
        push!(child_monoid.clique_indices, get_clique_indices(child_monoid, cliques)...)
        push!(
            child_monoid.commutes_with,
            setdiff(
                parent_monoid.vertices,
                union(cliques[child_monoid.clique_indices]...),
            )...,
        )
    end

    push!(parent_monoid.cliques, cliques...)

    # special handling for GraphProductMonoid{Variable} (PCMonomial optimization)
    if isa(parent_monoid.vertices[1], Variable)
        # Build the var_index_dict for O(1) variable lookup (PCMonomial optimization)
        parent_monoid.var_index_dict[] = Dict{Variable, UInt32}(
            v => UInt32(i) for (i, v) in enumerate(parent_monoid.vertices)
        )

        # Set the monomial field of each child variable to its corresponding PCMonomial representation, helps in multiplication of variables in GraphProductMonoid{Variable} 
        for child_monoid in parent_monoid.vertices
            child_monoid.monomial[]=words_to_monomial(parent_monoid, [child_monoid])
            if child_monoid.mult_type[] != :Free || !isempty(child_monoid.ortho_conj)
                parent_monoid.has_relations[]=true
            end
        end

        # Check if the the cliques in the dependency graph are closed under conjugation, if so, then taking conjugate just conjugates the cliques words : tuple sof non-commutative words.
        all(
            i -> (
                Set(i.commutes_with)==Set(conj(i).commutes_with) &&
                !(i in conj(i).commutes_with)
            ),
            parent_monoid.vertices,
        ) && (parent_monoid.conj_type[]=true)

        # Build the conj_indices mapping for O(1) conjugate lookup (PCMonomial optimization)
        n = length(parent_monoid.vertices)
        conj_indices = Vector{UInt32}(undef, n)
        var_index_dict = parent_monoid.var_index_dict[]
        @inbounds for i in 1:n
            v = parent_monoid.vertices[i]
            cv = conj(v)
            if cv === v || cv == v
                conj_indices[i] = UInt32(i)
            else
                # Use dictionary for O(1) lookup
                conj_indices[i] = var_index_dict[cv]
            end
        end
        parent_monoid.conj_indices[] = conj_indices
    end

    # if the parent monoid is a GraphProductMonoid of other monoids, build the child monoids as well
    if isa(parent_monoid.vertices[1], AbstractMonoid)
        build.(parent_monoid.vertices)
        parent_monoid.conj_type[]=true
    end

    #  Monoid now finally built
    parent_monoid.is_built[]=true
end

"""
    reset(M::Monoid{V})

    Resets the monoid `M` clearing all its fields.
"""
function reset(monoid::GraphProductMonoid)
    monoid.vertices=Vector{AbstractMonoid}()
    monoid.commutations=Vector{Tuple{AbstractMonoid, AbstractMonoid}}()
    monoid.cliques=Vector{Vector{AbstractMonoid}}()
    monoid.is_built[]=false
end

# one() for GraphProductMonoid{Variable} returns PCMonomial (fast path)
# Defined in PCMonomial.jl
"""
    # Notes
    one() for nested GraphProductMonoid (e.g., GraphProductMonoid{NCMonoid}) returns GraphProductWord

    # Arguments
    - `M`: The GraphProductMonoid for which the identity element is requested.

    # Returns
    - The identity element of the monoid, which is a GraphProductWord.

    # Example usage:
    ```juliaulia
    @pcmonoid M1 a b c
    @pcmonoid M2 x y z
    @pcmonoid M3 p q r
    GPM = GraphProductMonoid("GPM", [M1, M2, M3])
    identity_element = one(GPM)
    ```
    This will return the identity element of the GraphProductMonoid `GPM`.
"""
function Base.one(M::GraphProductMonoid{T}) where {T <: AbstractMonoid}
    t = package_types[T]
    return GraphProductWord(
        M,
        Base.RefValue{AbstractMonomial}(),
        Vector{Vector{t}}([[] for i in 1:length(M.cliques)]),
        Set{t}(),
        Set{t}(),
    )
end
