#= This file implements important functions used in represneting monomials or changing their representations.
1. word_to_monomial
2. monomial_to_graph
3. monomial_to_word
4. is_reconstructible
5. graph_to_foata
6. monomial_to_foata
=#


"""
    word_to_monomial(monoid::Monoid, word::Vector{V}) where V

    Converts a word to a monomial.

    # Arguments:
    - `monoid`: The monoid to which the word belongs.
    - `word`: The word to be converted.

    # Returns:
    - `monomial``: The converted pcword as monomial.
"""
function words_to_monomial(monoid::GraphProductMonoid, words::Vector{AW};edge_l=Set{AW}(),edge_r=Set{AW}()) where AW<:AbstractMonomial
    
    monomial = one(monoid)
    words=filter(i->!is_identity(i),words)
    isempty(words) && return monomial
    # for single variable
    if length(words) == 1
        word= words[1]
        clique_indices=isa(word,Variable) ? word.clique_indices : word.monoid.clique_indices
        for clique_index in clique_indices
            push!(monomial.clique_words[clique_index], word)
        end
        push!(monomial.edge_l, word)
        push!(monomial.edge_r, word)
        return monomial
    end
    for (clique_index, clique) in enumerate(monoid.cliques)
        push!(monomial.clique_words[clique_index], clique_projector(words, clique)...)
    end

    isempty(edge_l) ? union!(monomial.edge_l,get_edge_variables(monomial.clique_words,:first,collect(1:length(monoid.cliques))) ) : union!(monomial.edge_l,edge_l)
    isempty(edge_r) ? union!(monomial.edge_r,get_edge_variables(monomial.clique_words,:last,collect(1:length(monoid.cliques))) ) : union!(monomial.edge_r,edge_r)
    return monomial
end


"""
    monomial_to_graph(monomial::Monomial)

    Converts a monomial into a graph representation.

    # Arguments
    - `monomial`: A monomial in a monoid.

    # Returns
    - `graph`: A SimpleDiGraph where each node represents a variable in the monomial and each edge represents a pair of variables that do not commute.
    - `node_dict`: A dictionary mapping each variable and its occurrence in the monomial to a unique integer.

    # Notes
    The function first checks if a variable occurs differently in cliques. If it does, the function returns false. Otherwise, it generates labels for each variable and its occurrence in the monomial, labels the words in the monomial based on the cliques and the generated labels, finds pairs of labels that do not commute, and finally creates a graph from these pairs.
"""
function monomial_to_graph(monomial::GraphProductWord)

    cliques = monomial.monoid.cliques
    exponents = check_exponents_consistency(monomial)

    # check if a variable occurs differently in cliques
    (exponents == false) && return false

    node_dict = generate_labels(exponents)
    clique_word_node = label_clique_words(monomial.clique_words, cliques, node_dict)

    # find the pairs of labels which don't commute
    non_commuting_pairs = [(node[j], node[j+1]) for node in clique_word_node for j in 1:length(node)-1]

    graph = SimpleDiGraph(length(node_dict))
    E=Edge.(non_commuting_pairs)
    for e in E
          add_edge!(graph, e)
    end
    for e in E
        rem_edge!(graph,e)
        if !has_path(graph,e.src,e.dst)
            add_edge!(graph,e)
         end
    end
    return graph, node_dict,exponents
end


"""
    generate_labels(exponents::Vector{Pair{V, Int}}) where V

    Generates labels for each variable and its occurrence in a monomial.

    # Arguments
    - `exponents`: A vector of (variable=>it's degree/exponent)  degree is the number of times they occur in the monomial.

    # Returns
    - `node_dict`: A dictionary mapping each variable and its occurrence in the monomial to a unique integer.

    # Notes
    The function iterates over each variable and its degree in the `exponents` dictionary and for each occurrence of the variable, it assigns a unique integer.
"""
function generate_labels(exponents::Dict)

    # all the variables and their occurrences in the monomial are given incremental numbers to be used in graph
    #  i.e if check_exponents_consistency is {a=>2,b=>3,c=>1} then node_dict={(a,1)=>1,(a,2)=>2,(b,1)=>3,(b,2)=>4,(b,3)=>5,(c,1)=>6}
    incremental_counter = 0
    node_dict = eltype(keys(exponents))<:Variable ? Dict{Tuple{Variable,Int64},Int64}() : Dict{Tuple{AbstractMonomial,Int64},Int64}()

    for (monoid, words) in exponents
        for count in 1:length(words)
            incremental_counter = node_dict[(words[count], count)] = incremental_counter + 1
        end
    end
    return node_dict
end


"""
    label_clique_words(clique_words::Vector{Vector{V}}, cliques::Vector{Vector{V}}, node_dict::Dict{Tuple{V,Int64},Int64}) where V

    Labels the words in a monomial based on the cliques and a dictionary of nodes.

    # Arguments
    - `clique_words`: A vector of vectors representing the words in the monomial.
    - `cliques`: A vector of vectors representing the cliques in the monomial.
    - `node_dict`: A dictionary mapping each variable and its occurrence in the monomial to a unique integer.

    # Returns
    - `result`: A vector of vectors where each inner vector represents a word in the monomial and each element in the inner vector is the label of a variable in the word.

    # Notes
    The function iterates over each clique and word in the monomial, finds the indices where each variable in the clique occurs in the word, and updates the corresponding elements in the word with the node label for the variable.
"""
function label_clique_words(clique_words::Vector, cliques::Vector, node_dict::Dict)
    
    result = Vector{Vector{Int64}}()
    for (clique, word) in zip(cliques, clique_words)
        labelled_word = Array{Int64}(undef, length(word))
        for monoid in clique
            #  for each variable in clique, find the indices where it occurs in word and then use variable and the indices to update labelled_word with the node label for the corresponding variable 
            variable_indices = eltype(eltype(cliques))<: Variable ? findall(x -> x==monoid, word) : findall(x -> x.monoid==monoid, word)
            num_variable_occurrences = length(variable_indices)
            if(num_variable_occurrences != 0)
                [labelled_word[index] = node_dict[(word[index], occurrence)] for (index, occurrence) in zip(variable_indices, 1:num_variable_occurrences)]
            end
        end
        push!(result, labelled_word)
    end
    return result
end


"""
    is_reconstructible(monomial::Monomial{V}) where V

    Checks if a monomial is reconstructible.

    # Arguments
    - `monomial`: A monomial in a monoid.

    # Returns
    - `Boolean`: Returns true if the monomial is reconstructible, false otherwise.

    # Notes
    The function first converts the monomial into a graph. If the conversion fails, it returns false. Otherwise, it checks if the graph is cyclic. If the graph is cyclic, it returns false, indicating that the monomial is not reconstructible. If the graph is not cyclic, it returns true, indicating that the monomial is reconstructible.
"""
function is_reconstructible(monomial::GraphProductWord)

    graph_conversion_result = monomial_to_graph(monomial)
    if graph_conversion_result==false
        return false
    else
        graph, node_dict,exponents = graph_conversion_result
        return !is_cyclic(graph)
    end
end


"""
    clique_words_to_labelled_word(node_dict::Dict{Tuple{V,Int64},Int64}, monomial::Monomial{V}) where V

    Converts the clique words of a monomial to a labelled word based on a dictionary of nodes.

    # Arguments
    - `node_dict`: A dictionary mapping each variable and its occurrence in the monomial to a unique integer.
    - `monomial`: A monomial in a monoid.

    # Returns
    - `labelled_word`: A vector where each element is the label of a variable in the word.

    # Notes
    The function first labels the clique words in the monomial based on the cliques and the node dictionary. Then, it iterates over each labelled clique word and adds the labels that are not already in the word to the word.
"""
function clique_words_to_labelled_word(node_dict::Dict,exponents, monomial::GraphProductWord)
    
    labelled_word = Vector{Int}()
    labelled_clique_words = label_clique_words(monomial.clique_words, monomial.monoid.cliques, node_dict)
    for labelled_clique_word in labelled_clique_words
        union!(labelled_word, labelled_clique_word)
    end
    return labelled_word
end


"""
    monomial_to_word(monomial::Monomial{V}) where V

    Converts a monomial to a word.

    # Arguments
    - `monomial`: A monomial in a monoid.

    # Returns
    - `word`: A vector where each element is a variable in the word.

    # Notes
    The function first converts the monomial into a graph and a dictionary of nodes, and then converts the clique words in the monomial to a labelled word based on the dictionary of nodes. Then, it corrects the labelled word to the correct word by commuting variables in the word which can commute or are in the wrong order. Finally, it converts the labels in the corrected word back to variables using the reversed dictionary of nodes.
"""
function monomial_to_word(monomial::GraphProductWord)
    
    graph, node_dict,exponents = monomial_to_graph(monomial)
    reversed_node_dict = reverse_dict(node_dict)
    labelled_word = clique_words_to_labelled_word(node_dict,exponents, monomial)

    # Correct the labelled word to the correct word by commuting variables in the word which can commute or are in the wrong order
    for i in 2:length(labelled_word)
        for j in i-1:-1:1
            if !has_path(graph, labelled_word[j], labelled_word[j+1])
                temp = labelled_word[j]
                labelled_word[j] = labelled_word[j+1]
                labelled_word[j+1] = temp
            end
        end
        if labelled_word[i] == 0
            labelled_word[i] = labelled_word[i-1]
        end
    end
    word = [reversed_node_dict[i][1] for i in labelled_word]
    return word
end


"""
    graph_to_foata(graph::SimpleDiGraph, node_dict::Dict{Tuple{V,Int64},Int64}) where V

    Converts a graph to a Foata normal form.

    # Arguments
    - `graph`: A directed graph.
    - `node_dict`: A dictionary of nodes in the graph.

    # Returns
    - `foata_form`: A vector of monomials representing the Foata normal form of the graph.

    # Notes
    The function first reverses the dictionary of nodes and the graph. Then, it finds the source nodes in the graph and stores them in a dictionary. It then iterates over the dictionary of nodes and for each node that is not a source node, it calculates its shortest distance from the source nodes and stores it in the dictionary. Finally, it converts the dictionary to a vector of monomials representing the Foata normal form of the graph.
"""
function graph_to_foata(graph::SimpleDiGraph, node_dict::Dict)
    
    reversed_node_dict = reverse_dict(node_dict)
    graph = reverse(graph)
    distance_dict = Dict{Int64,Vector{eltype(keys(node_dict))}}()
    source_nodes = findall(i -> i == 0, indegree(graph))
    distance_dict[0] = getindex.(Ref(reversed_node_dict), source_nodes)
    for (node, index) in node_dict
        if !(index in source_nodes)
            distance = maximum([length(a_star(graph, source_node, index)) for source_node in source_nodes])
            (haskey(distance_dict, distance)) ? push!(distance_dict[distance], node) : distance_dict[distance] = [node]
        end
    end
    levels = length(distance_dict)
    # foata_form = Vector{GraphProductWord}(undef, levels)
    foata_form = Vector{Vector{Base.tuple_type_head(eltype(keys(node_dict)))}}(undef, levels)
    
    # for i in 0:levels-1
    #     foata_form[levels-i] = length(distance_dict[i]) == 1 ? monomial(distance_dict[i][1][1]) : prod([node[1] for node in distance_dict[i]])
    # end

    for i in 0:levels-1
        foata_form[levels-i] = length(distance_dict[i]) == 1 ? [distance_dict[i][1][1]] : sort([node[1] for node in distance_dict[i]])
    end
    return foata_form
end


monomial_to_foata(m::GraphProductWord)=graph_to_foata(monomial_to_graph(m)[1:2]...)
