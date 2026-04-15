struct CyclicWord <: AbstractMonomial
    monoid::AbstractMonoid
    graph::NautyDiGraph{UInt64}
    exponents::Dict
    node_dict::Dict
    ref_word::PCMonomial
end

function cyclic_reduce(m::PCMonomial)
    monoid = m.monoid
    
    m_words = [copy(ci) for ci in m.clique_words]
    m_l = copy(m.edge_l)  # Set{UInt32}
    m_r = copy(m.edge_r)  # Set{UInt32}
    
    conj_indices = build_conj_indices(monoid)
    
    for (i_idx, j_idx) in product(m_r, m_l)
        i_var = monoid.vertices[i_idx]
        j_var = monoid.vertices[j_idx]
        
        if j_var in i_var.ortho_conj
            return 0
        end
        
        # projector
        if i_idx == j_idx && i_var.mult_type[] == :Projector && length(m_words[i_var.clique_indices[1]]) != 1
            for index in i_var.clique_indices
                pop!(m_words[index])
            end
            delete!(m_r, i_idx)
        end
        
        # unitary and unipotent
        conj_i_idx = conj_indices[i_idx]
        if (i_idx == j_idx && i_var.mult_type[] == :Unipotent) || (conj_i_idx == j_idx && i_var.mult_type[] == :Unitary)
            if length(m_words[i_var.clique_indices[1]]) == 1
                continue
            end
            for index in i_var.clique_indices
                pop!(m_words[index])
                popfirst!(m_words[index])
            end
            delete!(m_r, i_idx)
            delete!(m_l, j_idx)
            union!(m_l, get_edge_indices(m_words, :first, i_var.clique_indices, monoid))
            union!(m_r, get_edge_indices(m_words, :last, i_var.clique_indices, monoid))
            return cyclic_reduce(PCMonomial(monoid, Base.RefValue{AbstractMonomial}(), m_words, m_l, m_r))
        end
    end
    
    reduced_word = PCMonomial(monoid, Base.RefValue{AbstractMonomial}(), m_words, m_l, m_r)
    
    # Build graph (returns UInt32-based node_dict and exponents)
    graph, node_dict_idx, exponents_idx = monomial_to_graph(reduced_word)
    
    # Add cyclic edges
    for word in m_words
        if length(word) > 1
            first_label = node_dict_idx[(word[1], 1)]
            last_label = node_dict_idx[(word[end], exponents_idx[word[end]])]
            add_edge!(graph, last_label, first_label)
        end
    end
    
    # Convert to Variable-based format for CyclicWord
    # node_dict = Dict(label => (monoid.vertices[var_idx], occ) for ((var_idx, occ), label) in node_dict_idx)
    node_dict = Dict(label => var_idx for ((var_idx, occ), label) in node_dict_idx)

    label_n = Vector{UInt32}(undef, length(node_dict))
    for (k, v) in node_dict
        label_n[k] = v
    end

    exponents = Dict(monoid.vertices[var_idx] => count for (var_idx, count) in exponents_idx)


    graph_nauty = DenseNautyGraph{true}(collect(edges(graph)); vertex_labels=label_n)
    return CyclicWord(monoid, graph_nauty, exponents, node_dict, reduced_word)
end

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

function cyclic_reduce(poly::Polynomial)
    r=Polynomial(cyclic_reduce(first(poly.monomials)),first(poly.coeffs))
    for i in 2:length(poly.monomials)
        c_mon=cyclic_reduce(poly.monomials[i])
        index=findfirst(x->x==c_mon, r.monomials)
        if index===nothing
            push!(r.monomials, c_mon)
            push!(r.coeffs, poly.coeffs[i])
        else
            r.coeffs[index] += poly.coeffs[i]
        end
    end
    return r
end

function Base.:(==)(c1::CyclicWord,c2::CyclicWord)
    (c1.exponents!=c2.exponents || c1.monoid!=c2.monoid) && return false
    #edges_1=Set([(c1.node_dict[e.src],c1.node_dict[e.dst]) for e in edges(c1.graph)])
    #edges_2=Set([(c2.node_dict[e.src],c2.node_dict[e.dst]) for e in edges(c2.graph)])
    # edges_1=[(c1.node_dict[e.src][1],c1.node_dict[e.dst][1]) for e in edges(c1.graph)]
    # edges_2=[(c2.node_dict[e.src][1],c2.node_dict[e.dst][1]) for e in edges(c2.graph)]
    # (length(edges_1)!=length(edges_2) || countmap(edges_1)!=countmap(edges_2)) && return false
    return NautyGraphs.is_isomorphic(c1.graph,c2.graph)
end

Base.:(==)(c1::CyclicWord,c2::GraphProductWord{Variable}) = c1==cyclic_reduce(c2)
Base.:(==)(c1::GraphProductWord{Variable},c2::CyclicWord) = cyclic_reduce(c1)==c2
Base.:(==)(c1::CyclicWord,c2::PCMonomial) = c1==cyclic_reduce(c2)
Base.:(==)(c1::PCMonomial,c2::CyclicWord) = cyclic_reduce(c1)==c2

function Base.hash(t::CyclicWord, h::UInt)
    edges_t=Set([(t.node_dict[e.src],t.node_dict[e.dst]) for e in edges(t.graph)])
    return hash(t.exponents,hash(edges_t,hash(t.monoid,hash(0x23269960ff982ff6, h))))
end