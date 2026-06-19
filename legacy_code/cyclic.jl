# Legacy version for GraphProductWord{Variable} - delegates to PCMonomial version
function cyclic_reduce(m::GraphProductWord{V}) where V<: Variable
    return cyclic_reduce(PCMonomial(m))
end

#= Legacy cyclic_reduce for GraphProductWord{Variable} - kept for reference
function cyclic_reduce_legacy(m::GraphProductWord{V}) where V<: Variable
    parent_monoid=m.monoid

    m_words = copy(m.clique_words)
    m_l=copy(m.edge_l)
    m_r=copy(m.edge_r)
    prod=product(m_r,m_l)
    for (i,j) in prod
        if j in i.ortho_conj
            return 0
        end
        # projector
        if i==j && i.mult_type[]==:Projector && length(m_words[i.clique_indices[1]])!=1
            for index in i.clique_indices
                pop!(m_words[index])
            end
            delete!(m_r,i)
        end
        # unitary and unipotent
        if (i==j && i.mult_type[]==:Unipotent) || (i'==j && i.mult_type[]==:Unitary)
            if length(m_words[i.clique_indices[1]])==1
                continue
            end
            for index in i.clique_indices
                pop!(m_words[index])
                popfirst!(m_words[index])
            end
            delete!(m_r,i)
            delete!(m_l,j)
            union!(m_l,get_edge_variables(m_words,:first,i.clique_indices))
            union!(m_r,get_edge_variables(m_words,:last,i.clique_indices))
            return cyclic_reduce(GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),m_words,m_l,m_r))
        end
    end
    reduced_word=GraphProductWord(parent_monoid,Base.RefValue{AbstractMonomial}(),m_words,m_l,m_r)
    graph, node_dict,exponents=monomial_to_graph(reduced_word)
    for word in m_words
        if length(word)>1
            first_letter_label=node_dict[word[1],1]
            last_letter_label= node_dict[word[end],length(exponents[word[end]])]
            add_edge!(graph,last_letter_label,first_letter_label)
        end
    end
    node_dict=Dict([(value,key) for (key,value) in node_dict])
    exponents=Dict([key=>length(value) for (key,value) in exponents])
    return CyclicWord(parent_monoid,graph,exponents,node_dict,Base.RefValue(reduced_word))
end
=#
