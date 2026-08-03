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

# use loose definition of edge equality to allow for edge values changing
Base.:(==)(a::Edge, b::Edge) = a.top_vertex == b.top_vertex && a.bott_vertex == b.bott_vertex


"""
    $TYPEDEF

Representation of a derivative as a DAG. Nodes represent the expression tree, and edges represent partial derivatives.

Used for the D* AD algorithm.

Postorder numbers refer to the indices of nodes in the graph, where the lowest numbers are variables and the highest are roots. 
"""
struct DerivativeGraph{T<:Integer}
    symbols::Vector{SymbolicT} # postorder number -> symbolic expression
    definitions::IdDict{SymbolicT, T} # symbolic expression -> postorder number
    roots::Vector{SymbolicT} # root index -> root symbolic expression
    vars::Vector{SymbolicT} # variable index -> variable symbolic expression
    varset::Set{SymbolicT} # for fast checking if an expression is a variable
    var_idx_to_postorder::Vector{T}
    postorder_to_var_idx::IdDict{T,T}
    root_idx_to_postorder::Vector{T}
    postorder_to_root_idx::IdDict{T,T}
    parent_edges::Dict{T, Vector{Edge{T}}} # node -> parent edges
    child_edges::Dict{T, Vector{Edge{T}}} # node -> child edges
    doms::Vector{Union{Nothing, T}} # immediate dominators of nodes
    pdoms::Vector{Union{Nothing, T}} # immediate postdominators of nodes
    dom_masks::Vector{BitVector} # BitVectors representing which nodes each node dominates (e.g. dom_masks[dominator][dominated])
    pdom_masks::Vector{BitVector} # BitVectors representing which nodes each node postdominates (e.g. pdom_masks[postdominator][postdominated])
end

"""
    DerivativeGraph(roots::AbstractVector{SymbolicT}, vars::AbstractVector{SymbolicT}, idx_type::Type=Int32) -> DerivativeGraph

Creates and populates a `DerivativeGraph`.

# Arguments

- `roots::AbstractVector{SymbolicT}`: Expressions to take the derivative of
- `vars::AbstractVector{SymbolicT}`: Variables to take the derivative with respect to
- `idx_type::Type`: Integer type used to store indices in the graph
    (**Default**: `Int32`)
"""
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

    append!(dg.doms, get_dominators(dg))
    append!(dg.pdoms, get_postdominators(dg))
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

# handles terms with >2 arguments (e.g. multiplication of 3+ things)
function nary_derivative_idx(expr::SymbolicT, arg_idx::Integer)
    @match expr begin
        BSImpl.AddMul(; coeff, dict, variant) => begin
            if variant == SymbolicUtils.AddMulVariant.ADD
                return COMMON_ONE
            else
                args = copy(parent(arguments(expr)))
                args[arg_idx] = COMMON_ONE
                return SymbolicUtils.mul_worker(VartypeT, args)
            end
        end
        _ => return derivative_idx(expr, arg_idx)
    end
end

# called in `DerivativeGraph` constructor to recursively iterate through the graph to fill out edges + reachabilities
function populate_dergraph!(dg::DerivativeGraph)
    for (root_idx, root) in enumerate(dg.roots)
        local post_idx
        if root in dg.varset
            post_idx = populate_dergraph_var!(dg, root, root_idx)
        else
            post_idx = populate_dergraph!(dg, root, root_idx)
        end

        dg.root_idx_to_postorder[root_idx] = post_idx
        dg.postorder_to_root_idx[post_idx] = root_idx
    end
end

function populate_dergraph!(dg::DerivativeGraph{T}, expr::SymbolicT, root_idx::Integer) where {T}
    haskey(dg.definitions, expr) && return populate_root_reachabilities!(dg, dg.definitions[expr], root_idx)

    !iscall(expr) && return

    args = parent(arguments(expr))
    arg_idx_to_post_idx = Vector{T}(undef, length(args))
    for (arg_idx, arg) in enumerate(args)
        if arg in dg.varset
            arg_idx_to_post_idx[arg_idx] = populate_dergraph_var!(dg, arg, root_idx)
        elseif iscall(arg)
            arg_idx_to_post_idx[arg_idx] = populate_dergraph!(dg, arg, root_idx)
        else
            arg_idx_to_post_idx[arg_idx] = T(-1)
        end
    end

    push!(dg.symbols, expr)
    post_idx::T = length(dg.symbols)
    dg.definitions[expr] = post_idx
    dg.child_edges[post_idx] = Edge{T}[]
    dg.parent_edges[post_idx] = Edge{T}[]

    # add new edges
    for (arg_idx, arg_post_idx) in enumerate(arg_idx_to_post_idx)
        arg_post_idx == T(-1) && continue
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
        if !child_edge.reachable_roots[root_idx]
            child_edge.reachable_roots[root_idx] = 1
            populate_root_reachabilities!(dg, child_edge.bott_vertex, root_idx)
        end
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

# for visualization purposes - converts a `DerivativeGraph` into an `OrderedDiGraph`
function generate_graph(dg::DerivativeGraph{T}) where {T}
    g = SymbolicUtils.OrderedDiGraph{T}(length(dg.symbols))
    for node in 1:length(dg.symbols)
        for edge in dg.child_edges[node]
            SymbolicUtils.Graphs.add_edge!(g, edge.top_vertex, edge.bott_vertex)
        end
    end

    return g
end

# Follows the algorithm described in this paper: https://www.cs.tufts.edu/comp/150FP/archive/keith-cooper/dom14.pdf
function get_dominators(dg::DerivativeGraph{T}) where {T}
    doms = Vector{Union{Nothing, T}}(undef, length(dg))
    root_idxs = root_postorders(dg)
    doms[root_idxs] = root_idxs

    # moves two nodes up the graph until they meet
    function get_common_parent(a::T, b::T)::Union{Nothing, T}
        # move a and b up the graph through their immediate dominators until they meet
        (isnothing(doms[a]) || isnothing(doms[b])) && return nothing
        while a != b
            !(a < b && a != doms[a]) && !(b < a && b != doms[b]) && return nothing
            while a < b && a != doms[a]
                a = doms[a]
            end
            while b < a && b != doms[b]
                b = doms[b]
            end
        end
        return a
    end

    changed = true # keeps track of when changes stop happening
    while changed
        changed = false
        for node in reverse(eachindex(dg))
            parents = parent_nodes(dg, node)
            
            isempty(parents) && continue
            
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

function get_postdominators(dg::DerivativeGraph{T}) where {T}
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

# represents a subgraph of a DerivativeGraph, defined by either a dominator or postdominator node. can be factored into a single edge
mutable struct FactorableSubgraph{T<:Integer, S<:AbstractFactorableSubgraph}
    subgraph_value::SymbolicT
    top_vertex::T
    bott_vertex::T
    reachable_vars::BitVector
    reachable_roots::BitVector
    edges::Set{Edge{T}}
    subgraph_count::Int

    function FactorableSubgraph{T, DominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector) where {T<:Integer}
        new{T, DominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, Set{Edge{T}}(), sum(reachable_vars)*sum(reachable_roots))
    end

    function FactorableSubgraph{T, PostDominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector) where {T<:Integer}
        new{T, PostDominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, Set{Edge{T}}(), sum(reachable_vars)*sum(reachable_roots))
    end
end

# Loops through the immediate dominators to fill out the dominance masks. Also works for postdominators.
function calculate_dominance_mask(dominators::Vector{T}) where {T}
    dom_mask = Vector{BitVector}(undef, length(dominators))
    for i in eachindex(dominators)
        dom_mask[i] = falses(length(dominators))
    end
    for (node, dom) in enumerate(dominators)
        isnothing(dom) && continue
        # walk up through the dominators
        while true
            dom_mask[dom][node] = 1
            dom == dominators[dom] && break
            dom = dominators[dom]
        end
    end

    return dom_mask
end

# the following functions allow functions to treat Dominator and PostDominator subgraphs the same by using forward and backward instead of up and down the graph
# forward is in the direction of dominated to dominating or postdominated to postdominating. i.e. from the factor base to the factor node
select_forward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = dg.dom_masks[forward_vertex(sub)]
select_forward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = dg.pdom_masks[forward_vertex(sub)]
select_backward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = dg.pdom_masks[backward_vertex(sub)]
select_backward_dominance_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = dg.dom_masks[backward_vertex(sub)]

forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = parent_edges(dg, edge.top_vertex)
forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = child_edges(dg, edge.bott_vertex)
forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, DominatorSubgraph}, node::T) where {T} = parent_edges(dg, node)
forward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, PostDominatorSubgraph}, node::T) where {T} = child_edges(dg, node)

backward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = child_edges(dg, edge.bott_vertex)
backward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = parent_edges(dg, edge.top_vertex)
backward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, DominatorSubgraph}, node::T) where {T} = child_edges(dg, node)
backward_edges(dg::DerivativeGraph{T}, ::FactorableSubgraph{T, PostDominatorSubgraph}, node::T) where {T} = parent_edges(dg, node)

forward_vertex(::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = edge.top_vertex
forward_vertex(::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = edge.bott_vertex
backward_vertex(::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = edge.bott_vertex
backward_vertex(::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = edge.top_vertex
forward_vertex(sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = sub.top_vertex
forward_vertex(sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = sub.bott_vertex
backward_vertex(sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = sub.bott_vertex
backward_vertex(sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = sub.top_vertex

subgraph_edges(sub::FactorableSubgraph) = sub.edges

# NOT called in the constructor of `FactorableSubgraph`. Instead, delayed until right before factoring to account for changes to the `DerivativeGraph` from previous factoring
# Also calculates the subgraph's value
# Recurses through the subgraph to find all its edges
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

    # reached the end of the subgraph
    if forward_vertex(sub, edge) == forward_vertex(sub)
        push!(edges, edge)
        return edge.edge_value
    end

    !dom_mask[forward_vertex(sub, edge)] && return COMMON_ZERO # edge not in subgraph

    push!(edges, edge)
    path_value = COMMON_ZERO
    for next_edge in forward_edges(dg, sub, edge)
        path_value += _subgraph_edges!(edges, dg, sub, next_edge, visited, dom_mask)
    end

    path_value *= edge.edge_value

    return path_value
end

# the number of times a subgraph is used by all possible partial derivatives
subgraph_count(sub::FactorableSubgraph) = sub.subgraph_count

"""
    get_factorable_subgraphs(dg::DerivativeGraph{T}) where {T} -> BinaryHeap

Generates a heap of factorable subgraphs ordered with `FactorOrder` (smallest + most used are factored first).
"""
function get_factorable_subgraphs(dg::DerivativeGraph{T}) where {T}
    subs = DataStructures.BinaryHeap{Union{FactorableSubgraph{T, DominatorSubgraph}, FactorableSubgraph{T, PostDominatorSubgraph}}, FactorOrder}()
    for (dominated, dominating) in enumerate(dg.doms)
        if !isnothing(dominating) && length(child_edges(dg, dominating)) > 1 && length(parent_edges(dg, T(dominated))) > 1
            reachable_vars_mask = reachable_vars(dg, T(dominated))
            reachable_roots_mask = reachable_roots(dg, dominating)
            push!(subs, FactorableSubgraph{T, DominatorSubgraph}(dominating, T(dominated), reachable_vars_mask, reachable_roots_mask))
        end
    end

    for (postdominated, postdominating) in enumerate(dg.pdoms)
        if !isnothing(postdominating) && length(parent_edges(dg, postdominating)) > 1 && length(child_edges(dg, T(postdominated))) > 1
            reachable_vars_mask = reachable_vars(dg, postdominating)
            reachable_roots_mask = reachable_roots(dg, T(postdominated))
            push!(subs, FactorableSubgraph{T, PostDominatorSubgraph}(T(postdominated), postdominating, reachable_vars_mask, reachable_roots_mask))
        end
    end
    
    return subs
end

# factors a subgraph from dg, replacing it with a single edge (keeping original edges when necessary)
function factor_subgraph!(dg::DerivativeGraph{T}, sub::FactorableSubgraph) where {T}
    # TODO: add in case for branching subgraph (factoring creates new subgraph inside of another)
    # TODO: add complete checks if subgraph is still valid
    # check that the factor and factor base nodes are still a factor and factor base
    (length(backward_edges(dg, sub, forward_vertex(sub))) < 2 || length(forward_edges(dg, sub, backward_vertex(sub))) < 2) && return false

    populate_subgraph_edges!(dg, sub)
    sub_edges = subgraph_edges(sub)
    dom_mask = select_backward_dominance_mask(dg, sub)

    for edge in sub_edges
        # if edge is completely contained in the subgraph, remove it
        if dom_mask[backward_vertex(sub, edge)] || backward_vertex(sub, edge) == backward_vertex(sub)
            rem_edge!(dg, edge)
        end
    end

    # add final subgraph edge
    add_edge!(dg, Edge{T}(sub.subgraph_value, sub.top_vertex, sub.bott_vertex, sub.reachable_vars, sub.reachable_roots))

    return true
end

# ordering subgraphs should be factored in
struct FactorOrder <: Base.Order.Ordering
end

Base.lt(::FactorOrder, a, b) = factor_order(a, b)
Base.isless(::FactorOrder, a, b) = factor_order(a, b)

function factor_order(a::FactorableSubgraph, b::FactorableSubgraph)
    a_diff = abs(a.top_vertex - a.bott_vertex)
    b_diff = abs(b.top_vertex - b.bott_vertex)

    # factor the smaller subgraph first (guarantees that if a ⊂ b then a is factored first)
    # then, factor the more used subgraph first

    return a_diff < b_diff || (a_diff == b_diff && subgraph_count(a) > subgraph_count(b))
end

# Factor all subgraphs in the `DerivativeGraph`. This is the key step in the D* algorithm.
function factor_subgraphs!(dg::DerivativeGraph)
    subs = get_factorable_subgraphs(dg)

    while !isempty(subs)
        # factor the first subgraph according to `FactorOrder``
        sub = pop!(subs)
        factor_subgraph!(dg, sub)

        # update dominance after factoring
        # TODO: more efficiently propogate changes to other subgraphs?
        dg.doms .= get_dominators(dg)
        dg.pdoms .= get_postdominators(dg)
        dg.dom_masks .= calculate_dominance_mask(dg.doms)
        dg.pdom_masks .= calculate_dominance_mask(dg.pdoms)
    end
end

# evaluate a FULLY FACTORED DerivativeGraph
# currently only accounts for R1->R1 functions by evaluating a single path of edges
function evaluate_paths(dg::DerivativeGraph)
    if length(dg.roots) == 1 && length(dg.vars) == 1
        return evaluate_path(dg, only(child_edges(dg, dg.root_idx_to_postorder[1])), dg.var_idx_to_postorder[1])
    end
end

# evaluate the product of single linear path of edges (this represents applying the chain rule)
function evaluate_path(dg::DerivativeGraph{T}, edge::Edge{T}, goal::T) where T
    edge.bott_vertex == goal && return edge.edge_value

    result = COMMON_ZERO
    for child_edge in child_edges(dg, edge.bott_vertex)
        result += evaluate_path(dg, child_edge, goal)
    end

    return result * edge.edge_value
end

"""
    dstar_derivative(roots::AbstractVector{SymbolicT}, vars::AbstractVector{SymbolicT}) -> SymbolicT

Computes the derivative of `roots` w.r.t. `vars` using the D* automatic differentiation algorithm using the `DerivativeGraph` data structure.

(see [this paper](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/main-65.pdf) for more details on the algorithm)

# Arguments
- `roots::AbstractVector{SymbolicT}`: Vector of expressions to differentate
- `vars::AbstractVector{SymbolicT}`: Vector of variables to differentate w.r.t.
"""
function dstar_derivative(roots::AbstractVector{SymbolicT}, vars::AbstractVector{SymbolicT})
    dg = DerivativeGraph(roots, vars)
    factor_subgraphs!(dg)

    return evaluate_paths(dg)
end

dstar_derivative(roots::AbstractVector, vars::AbstractVector) = dstar_derivative(unwrap.(roots), unwrap.(vars))
dstar_derivative(root, var) = dstar_derivative([root], [var])
