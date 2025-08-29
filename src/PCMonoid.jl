function PCMonoid_generator(GM::GraphProductMonoid)
    vars=GM.vertices
    comms=GM.commutations
    name=GM.name
    graph=SimpleGraph(length(vars))

    for (i, j) in comms
        add_edge!(graph,Edge(findfirst(isequal(i),vars),findfirst(isequal(j),vars)))
    end
    graph=complement(graph)
    cliques=maximal_cliques(graph)
    cliques=[[vars[index] for index in clique] for clique in maximal_cliques(graph)]
    types=Dict{Vector{Vector{Variable}},Vector{Variable}}()
    for var in vars
        k=filter(i->var in i,cliques)
        if haskey(types,k)
            push!(types[k],var)
        else
            types[k]=[var]
        end
    end
    monoids=Vector{NCMonoid}([])
    # var_dict=Dict{}()
    # subsystems=Dict{NCMonoid,Vector{Int}}()
    for (index,(key,value)) in enumerate(types)
        ncmon=NCMonoid("$name$index",value)
        # for (i,j) in enumerate(value)
        #     var_dict[j]=ncmon.vertices[i]
        #     if haskey(relations,j)
        #         if relations[j]==:Unipotent
        #             ncmon.relations[(ncmon.vertices[i],ncmon.vertices[i])]=Tuple([])
        #         elseif relations[j]==:Projector
        #             ncmon.relations[(ncmon.vertices[i],ncmon.vertices[i])]=Tuple([ncmon.vertices[i]])
        #         end
        #     end
        # end
        # subsystems[ncmon]=value
        push!(monoids,ncmon)
    end
    
    monoid_comms=Vector{Tuple{NCMonoid,NCMonoid}}()
    for (i,j) in combinations(monoids,2)
       edges=map(((l,r), )-> ( findfirst(isequal(l),vars) , findfirst(isequal(r),vars) ),product(i.vertices,j.vertices))
       any(k->has_edge(graph,k),edges) || push!(monoid_comms,(i,j))
    end
    res=GraphProductMonoid{NCMonoid}(name,monoids;commutations=monoid_comms)
    build(res)
    return res
end
# name_creator(x::Vector{Vector{Int}})=join(join.(x,"_"),"?")
