abstract type AbstractFactorableSubgraph end
abstract type DominatorSubgraph <: AbstractFactorableSubgraph end
abstract type PostDominatorSubgraph <: AbstractFactorableSubgraph end

struct Edge{T<:Integer}
    edge_value::SymbolicT
    top_vertex::T
    bott_vertex::T
    reachable_vars::BitVector
    reachable_roots::BitVector
end

edge_value(edge::Edge) = edge.edge_value
top_vertex(edge::Edge) = edge.top_vertex
bott_vertex(edge::Edge) = edge.bott_vertex
reachable_vars(edge::Edge) = edge.reachable_vars
reachable_roots(edge::Edge) = edge.reachable_roots
vertices(edge::Edge) = (top_vertex(edge), bott_vertex(edge))
times_used(edge::Edge) = sum(reachable_roots(edge)) * sum(reachable_vars(edge))

struct DerivativeGraph{T<:Integer}
    symbols::Vector{SymbolicT} # postorder -> symbolic expr
    definitions::IdDict{SymbolicT, T} # symbolic expr -> postorder
    roots::Vector{SymbolicT}
    vars::Vector{SymbolicT}
    varset::Set{SymbolicT} # for fast checking if expr is var
    var_idx_to_postorder::Vector{T}
    postorder_to_var_idx::IdDict{T,T}
    root_idx_to_postorder::Vector{T}
    postorder_to_root_idx::IdDict{T,T}
    parent_edges::Dict{T, Vector{Edge{T}}}
    child_edges::Dict{T, Vector{Edge{T}}}
    doms::Vector{Union{Nothing, T}}
    pdoms::Vector{Union{Nothing, T}}
    dom_masks::Vector{BitVector}
    pdom_masks::Vector{BitVector}
end

function DerivativeGraph(roots::AbstractVector{SymbolicT}, vars::AbstractVector{SymbolicT}, idx_type::Type=Int32)
    dg = DerivativeGraph{idx_type}(
        SymbolicT[],
        IdDict{SymbolicT, idx_type}(),
        roots,
        vars,
        Set(vars),
        Vector{idx_type}(undef, length(vars)),
        IdDict{idx_type, idx_type}(),
        Vector{idx_type}(undef, length(roots)),
        IdDict{idx_type, idx_type}(),
        Dict{idx_type, Vector{Edge{idx_type}}}(),
        Dict{idx_type, Vector{Edge{idx_type}}}(),
        Union{Nothing, idx_type}[],
        Union{Nothing, idx_type}[],
        BitVector[],
        BitVector[]
    )

    populate_dergraph!(dg)

    sizehint!(dg.doms, length(dg.symbols))
    sizehint!(dg.pdoms, length(dg.symbols))
    sizehint!(dg.dom_masks, length(dg.symbols))
    sizehint!(dg.pdom_masks, length(dg.symbols))

    append!(dg.doms, _get_dominators(dg))
    append!(dg.pdoms, _get_postdominators(dg))
    append!(dg.dom_masks, calculate_dominance_mask(dg.doms))
    append!(dg.pdom_masks, calculate_dominance_mask(dg.pdoms))

    return dg
end

Base.eachindex(dg::DerivativeGraph{T}) where {T} = T(1):T(length(dg.symbols)) # iterator over postorder indices with type T
Base.length(dg::DerivativeGraph) = length(dg.symbols)

root_postorders(dg::DerivativeGraph) = dg.root_idx_to_postorder

# these return references to the internal data structures, not for external use
parent_edges(dg::DerivativeGraph{T}, node::T) where {T} = dg.parent_edges[node]
child_edges(dg::DerivativeGraph{T}, node::T) where {T} = dg.child_edges[node]
parent_nodes(dg::DerivativeGraph{T}, node::T) where {T} = top_vertex.(parent_edges(dg, node))
child_nodes(dg::DerivativeGraph{T}, node::T) where {T} = bott_vertex.(child_edges(dg, node))

function hasedge(dg::DerivativeGraph{T}, edge::Edge{T}) where {T}
    is_child_edge = edge in child_edges(dg, top_vertex(edge))
    is_parent_edge = edge in parent_edges(dg, bott_vertex(edge))

    @assert is_child_edge == is_parent_edge "edge is only in one of the child/parent edge lists"

    return is_child_edge && is_parent_edge
end

function rem_edge!(dg::DerivativeGraph{T}, edge::Edge{T}) where {T}
    @assert hasedge(dg, edge) "edge is not in the graph"
    
    top_vert, bott_vert = vertices(edge)
    top_child_edges = child_edges(dg, top_vert)
    bott_parent_edges = parent_edges(dg, bott_vert)
    deleteat!(top_child_edges, findfirst(isequal(edge), top_child_edges))
    deleteat!(bott_parent_edges, findfirst(isequal(edge), bott_parent_edges))

    return nothing
end

function add_edge!(dg::DerivativeGraph{T}, top_vertex::T, bott_vertex::T, edge_value::SymbolicT) where {T}
    var_reachability = reachable_vars(dg, bott_vertex)
    root_reachability = reachable_roots(dg, bott_vertex)
    new_edge = Edge{T}(edge_value, top_vertex, bott_vertex, var_reachability, root_reachability)

    add_edge!(dg, new_edge)

    return nothing
end

function add_edge!(dg::DerivativeGraph{T}, edge::Edge{T}) where {T}
    top_vertex, bott_vertex = vertices(edge)
    push!(dg.child_edges[top_vertex], edge)
    push!(dg.parent_edges[bott_vertex], edge)

    return nothing
end

function reachable_roots(dg::DerivativeGraph{T}, node::T) where {T}
    edges = dg.parent_edges[node]
    roots_mask = falses(length(dg.roots))
    for edge in edges
        roots_mask .|= reachable_roots(edge)
    end

    if haskey(dg.postorder_to_root_idx, node)
        roots_mask[dg.postorder_to_root_idx[node]] = 1
    end

    return roots_mask
end

function reachable_vars(dg::DerivativeGraph{T}, node::T) where {T}
    edges = dg.child_edges[node]
    vars_mask = falses(length(dg.vars))
    for edge in edges
        vars_mask .|= reachable_vars(edge)
    end

    if haskey(dg.postorder_to_var_idx, node)
        vars_mask[dg.postorder_to_var_idx[node]] = 1
    end

    return vars_mask
end

# handles multiplication of >2 arguments
function nary_derivative_idx(expr::SymbolicT, arg_idx::Integer)
    @match expr begin
        BSImpl.AddMul(; coeff, dict, variant) => begin
            if variant == SymbolicUtils.AddMulVariant.ADD
                return COMMON_ONE
            else
                args = copy(parent(arguments(expr)))
                args[arg_idx] = COMMON_ONE
                return prod(args)
            end
        end
        _ => return derivative_idx(expr, arg_idx)
    end
end

function populate_dergraph!(dg::DerivativeGraph)
    for (root_idx, root) in enumerate(dg.roots)
        local post_idx
        if root in dg.varset
            post_idx = populate_dergraph_var!(dg, root, root_idx)
        elseif iscall(root)
            post_idx = populate_dergraph!(dg, root, root_idx)
        end

        dg.root_idx_to_postorder[root_idx] = post_idx
        dg.postorder_to_root_idx[post_idx] = root_idx
    end
end

function populate_dergraph!(dg::DerivativeGraph{T}, expr::SymbolicT, root_idx::Integer) where {T}
    haskey(dg.definitions, expr) && return populate_root_reachabilities!(dg, dg.definitions[expr], root_idx)

    # assume iscall(expr)
    args = parent(arguments(expr))
    arg_idx_to_post_idx = IdDict{Int, T}()
    sizehint!(arg_idx_to_post_idx, length(args)) # faster?
    for (arg_idx, arg) in enumerate(args)
        if arg in dg.varset
            arg_idx_to_post_idx[arg_idx] = populate_dergraph_var!(dg, arg, root_idx)
        elseif iscall(arg)
            arg_idx_to_post_idx[arg_idx] = populate_dergraph!(dg, arg, root_idx)
        end
    end

    push!(dg.symbols, expr)
    post_idx::T = length(dg.symbols)
    dg.definitions[expr] = post_idx
    dg.child_edges[post_idx] = Edge{T}[]
    dg.parent_edges[post_idx] = Edge{T}[]

    for (arg_idx, arg_post_idx) in arg_idx_to_post_idx
        partial_der = nary_derivative_idx(expr, arg_idx)
        arg_reachable_vars = reachable_vars(dg, arg_post_idx)
        arg_reachable_roots = reachable_roots(dg, arg_post_idx)
        arg_reachable_roots[root_idx] = 1
        new_edge = Edge{T}(partial_der, post_idx, arg_post_idx, arg_reachable_vars, arg_reachable_roots)
        push!(dg.child_edges[post_idx], new_edge)
        push!(dg.parent_edges[arg_post_idx], new_edge)
    end

    return post_idx
end


function populate_root_reachabilities!(dg::DerivativeGraph{T}, node::T, root_idx::Integer) where {T}
    for child_edge in dg.child_edges[node]
        child_edge.reachable_roots[root_idx] = 1
        populate_root_reachabilities!(dg, child_edge.bott_vertex, root_idx)
    end

    return node
end

function populate_dergraph_var!(dg::DerivativeGraph{T}, var::SymbolicT, root_idx::Integer) where {T}
    haskey(dg.definitions, var) && return populate_root_reachabilities!(dg, dg.definitions[var], root_idx)

    push!(dg.symbols, var)
    post_idx::T = length(dg.symbols) # postorder number
    dg.definitions[var] = post_idx
    var_idx = findfirst(isequal(var), dg.vars)
    dg.var_idx_to_postorder[var_idx] = post_idx
    dg.postorder_to_var_idx[post_idx] = var_idx
    dg.child_edges[post_idx] = Edge{T}[]
    dg.parent_edges[post_idx] = Edge{T}[]

    return post_idx
end

function generate_graph(dg::DerivativeGraph{T}) where {T}
    g = SymbolicUtils.OrderedDiGraph{T}(length(dg.symbols))
    for node in 1:length(dg.symbols)
        for edge in dg.child_edges[node]
            SymbolicUtils.Graphs.add_edge!(g, edge.top_vertex, edge.bott_vertex)
        end
    end

    return g
end

function _get_dominators(dg::DerivativeGraph{T}) where {T}
    doms = Vector{Union{Nothing, T}}(undef, length(dg))
    root_idxs = root_postorders(dg)
    doms[root_idxs] = root_idxs

    function get_common_parent(a::T, b::T)::Union{Nothing, T}
        # move a and b up the graph through their immediate dominators until they meet
        while a != b
            (isnothing(a) || isnothing(b)) && return nothing
            !(a < b && a != doms[a]) && !(b < a && b != doms[b]) && return nothing
            while !isnothing(a) && a < b && a != doms[a]
                a = doms[a]
            end
            while !isnothing(b) && b < a && b != doms[b]
                b = doms[b]
            end
        end
        return a
    end

    changed = true # keeps track of when changes stop happening
    while changed
        changed = false
        for node in reverse(eachindex(dg))
            node in root_idxs && continue
            
            parents = parent_nodes(dg, node)
            
            new_idom = first(parents)

            for parent in parents
                parent == first(parents) && continue
                if isassigned(doms, parent)
                    new_idom = get_common_parent(parent, new_idom)
                    isnothing(new_idom) && break
                end
            end

            if doms[node] != new_idom
                doms[node] = new_idom
                changed = true
            end
        end
    end


    return doms
end

function _get_postdominators(dg::DerivativeGraph{T}) where {T}
    pdoms = Vector{Union{Nothing, T}}(undef, length(dg))
    pdoms[dg.var_idx_to_postorder] = dg.var_idx_to_postorder

    function get_common_child(a::T, b::T)::Union{Nothing, T}
        # move a and b up the graph through their immediate dominators until they meet
        (isnothing(pdoms[a]) || isnothing(pdoms[b])) && return nothing
        while a != b
            !(a > b && a != pdoms[a]) && !(b > a && b != pdoms[b]) && return nothing
            while a > b && a != pdoms[a]
                a = pdoms[a]
            end
            while b > a && b != pdoms[b]
                b = pdoms[b]
            end
        end
        return a
    end

    changed = true # keeps track of when changes stop happening
    while changed
        changed = false
        for node in eachindex(dg)
            node in dg.var_idx_to_postorder && continue

            children = child_nodes(dg, node)
            if isempty(children)
                pdoms[node] = nothing
                continue
            end
            new_pidom = first(children)

            for child_idx in 2:length(children)
                if isassigned(pdoms, children[child_idx])
                    new_pidom = get_common_child(children[child_idx], new_pidom)
                    isnothing(new_pidom) && break
                end
            end

            if pdoms[node] != new_pidom
                pdoms[node] = new_pidom
                changed = true
            end
        end
    end


    return pdoms
end

mutable struct FactorableSubgraph{T<:Integer, S<:AbstractFactorableSubgraph}
    subgraph_value::SymbolicT
    top_vertex::T
    bott_vertex::T
    reachable_vars::BitVector
    reachable_roots::BitVector
    edges::Set{Edge{T}}

    function FactorableSubgraph{T, DominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector) where {T<:Integer}
        new{T, DominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, Set{Edge{T}}())
    end

    function FactorableSubgraph{T, PostDominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector) where {T<:Integer}
        new{T, PostDominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, Set{Edge{T}}())
    end
end

# constructs a factorable subgraph and calculates its edges and subgraph value
function FactorableSubgraph{T, S}(dg::DerivativeGraph{T}, top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector) where {T<:Integer, S<:AbstractFactorableSubgraph}
    sub = FactorableSubgraph{T, S}(top_vertex, bott_vertex, reachable_vars, reachable_roots)
    populate_subgraph_edges!(dg, sub)
    return sub
end

# also works for postdominators
function calculate_dominance_mask(dominators::Vector{T}) where {T}
    dom_mask = Vector{BitVector}(undef, length(dominators))
    for i in eachindex(dominators)
        dom_mask[i] = falses(length(dominators))
    end
    for (node, dom) in enumerate(dominators)
        dom_mask[dom][node] = 1
    end

    return dom_mask
end

function select_forward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, DominatorSubgraph}) where {T}
    return dg.dom_masks[forward_vertex(sub)]
end

function select_forward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T}
    return dg.pdom_masks[forward_vertex(sub)]
end

function select_backward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, DominatorSubgraph}) where {T}
    return dg.pdom_masks[backward_vertex(sub)]
end

function select_backward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T}
    return dg.dom_masks[backward_vertex(sub)]
end

# forward is in the direction of dominated to dominating or postdominated to postdominating
forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = parent_edges(dg, edge.top_vertex)
forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = child_edges(dg, edge.bott_vertex)
forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, DominatorSubgraph}, node::T) where {T} = parent_edges(dg, node)
forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, PostDominatorSubgraph}, node::T) where {T} = child_edges(dg, node)

forward_vertex(::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = edge.top_vertex
forward_vertex(::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = edge.bott_vertex
backward_vertex(::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = edge.bott_vertex
backward_vertex(::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = edge.top_vertex
forward_vertex(sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = sub.top_vertex
forward_vertex(sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = sub.bott_vertex
backward_vertex(sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = sub.bott_vertex
backward_vertex(sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = sub.top_vertex

function subgraph_edges(sub::FactorableSubgraph)
    return sub.edges
end

function populate_subgraph_edges!(dg::DerivativeGraph{T}, sub::FactorableSubgraph) where {T}
    union!(sub.edges, _subgraph_edges(dg, sub, backward_vertex(sub), select_forward_dominance_mask(dg, sub)))
end

function _subgraph_edges(dg::DerivativeGraph{T}, sub::FactorableSubgraph, node::T, dom_mask::BitVector) where {T}
    visited = Set{Edge{T}}()
    sub_edges = Set{Edge{T}}()
    for edge in forward_edges(dg, sub, node)
        sub.subgraph_value += _subgraph_edges!(sub_edges, dg, sub, edge, visited, dom_mask)
    end

    return sub_edges
end

function _subgraph_edges!(edges::Set{Edge{T}}, dg::DerivativeGraph{T}, sub::FactorableSubgraph{T}, edge::Edge{T}, visited::Set{Edge{T}}, dom_mask::BitVector) where {T}
    edge in visited && return COMMON_ZERO
    
    push!(visited, edge)

    if forward_vertex(sub, edge) == forward_vertex(sub)
        push!(edges, edge)
        return edge.edge_value
    end

    !dom_mask[forward_vertex(sub, edge)] && return COMMON_ZERO

    push!(edges, edge)
    path_value = COMMON_ZERO
    for next_edge in forward_edges(dg, sub, edge)
        path_value += _subgraph_edges!(edges, dg, sub, next_edge, visited, dom_mask)
    end

    path_value *= edge.edge_value

    return path_value
end

function get_factorable_subgraphs(dg::DerivativeGraph{T}) where {T}
    dom_factorable_subgraphs = Set{FactorableSubgraph{T, DominatorSubgraph}}()
    for (dominated, dominating) in enumerate(dg.doms)
        if !isnothing(dominating) && length(child_edges(dg, dominating)) > 1 && length(parent_edges(dg, T(dominated))) > 1
            reachable_vars_mask = reachable_vars(dg, T(dominated))
            reachable_roots_mask = reachable_roots(dg, dominating)
            push!(dom_factorable_subgraphs, FactorableSubgraph{T, DominatorSubgraph}(dg, dominating, T(dominated), reachable_vars_mask, reachable_roots_mask))
        end
    end
    pdom_factorable_subgraphs = Set{FactorableSubgraph{T, PostDominatorSubgraph}}()
    for (postdominated, postdominating) in enumerate(dg.pdoms)
        if !isnothing(postdominating) && length(parent_edges(dg, postdominating)) > 1 && length(child_edges(dg, T(postdominated))) > 1
            reachable_vars_mask = reachable_vars(dg, postdominating)
            reachable_roots_mask = reachable_roots(dg, T(postdominated))
            push!(pdom_factorable_subgraphs, FactorableSubgraph{T, PostDominatorSubgraph}(dg, T(postdominated), postdominating, reachable_vars_mask, reachable_roots_mask))
        end
    end
    
    return union(dom_factorable_subgraphs, pdom_factorable_subgraphs)
end

function factor_subgraph!(dg::DerivativeGraph{T}, sub::FactorableSubgraph) where {T}
    sub_edges = subgraph_edges(sub)
    dom_mask = select_backward_dominance_mask(dg, sub)
    
    for edge in sub_edges
        if dom_mask[backward_vertex(sub, edge)] || backward_vertex(sub, edge) == backward_vertex(sub)
            rem_edge!(dg, edge)
        end
    end

    add_edge!(dg, sub.top_vertex, sub.bott_vertex, sub.subgraph_value)
end