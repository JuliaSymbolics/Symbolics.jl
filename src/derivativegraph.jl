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
Base.hash(e::Edge, h::UInt) = hash((e.top_vertex, e.bott_vertex), h)

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
    var_idx_to_postorder::Dict{Int,T}
    postorder_to_var_idx::Dict{T,T}
    root_idx_to_postorder::Dict{Int,T}
    postorder_to_root_idx::Dict{T,T}
    parent_edges::Dict{T, Vector{Edge{T}}} # node -> parent edges
    child_edges::Dict{T, Vector{Edge{T}}} # node -> child edges
    dirty_roots::BitVector # roots touched by factoring a subgraph that need to have doms/pdoms recomputed
    dirty_vars::BitVector # variables touched by factoring a subgraph that need to have doms/pdoms recomputed
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
        IdDict{Int, idx_type}(),
        IdDict{idx_type, idx_type}(),
        IdDict{Int, idx_type}(),
        IdDict{idx_type, idx_type}(),
        Dict{idx_type, Vector{Edge{idx_type}}}(),
        Dict{idx_type, Vector{Edge{idx_type}}}(),
        trues(length(roots)),
        trues(length(vars))
    )

    populate_dergraph!(dg)

    sizehint!(dg.var_idx_to_postorder, length(vars))
    sizehint!(dg.root_idx_to_postorder, length(roots))

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

function propogate_var_reachability(dg::DerivativeGraph{T}, node::T) where {T}
    _parent_edges = parent_edges(dg, node)
    isempty(_parent_edges) && return nothing

    new_reachability = reachable_vars(dg, node) # computes reachable_vars from union of child edge reachabilities

    # shrink parent reachabilities
    for edge in _parent_edges
        old_reachability = copy(edge.reachable_vars)
        edge.reachable_vars .&= new_reachability
        if edge.reachable_vars != old_reachability
            dg.dirty_vars .|= old_reachability .& .~edge.reachable_vars
            propogate_var_reachability(dg, top_vertex(edge))
        end
    end
end
propogate_var_reachability(dg::DerivativeGraph{T}, edge::Edge{T}) where {T} = propogate_var_reachability(dg, edge.top_vertex)

function propogate_root_reachability(dg::DerivativeGraph{T}, node::T) where {T}
    _child_edges = child_edges(dg, node)
    isempty(_child_edges) && return nothing

    new_reachability = reachable_roots(dg, node) # computes reachable_roots from union of parent edge reachabilities

    # shrink child reachabilities
    for edge in _child_edges
        old_reachability = copy(edge.reachable_roots)
        edge.reachable_roots .&= new_reachability
        if edge.reachable_roots != old_reachability
            dg.dirty_roots .|= old_reachability .& .~edge.reachable_roots
            propogate_root_reachability(dg, bott_vertex(edge))
        end
    end
end

propogate_root_reachability(dg::DerivativeGraph{T}, edge::Edge{T}) where {T} = propogate_root_reachability(dg, edge.bott_vertex)

function rem_edge!(dg::DerivativeGraph{T}, edge::Edge{T}) where {T}
    @assert hasedge(dg, edge) "edge is not in the graph"

    # removing this edge can only affect doms/pdoms for its roots+vars, so mark them as dirty for recomputation
    dg.dirty_roots .|= edge.reachable_roots
    dg.dirty_vars .|= edge.reachable_vars

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

    # if edge already exists, merge in roots and vars
    existing_idx = findfirst(isequal(edge), dg.child_edges[top_vertex])
    if !isnothing(existing_idx)
        existing = dg.child_edges[top_vertex][existing_idx]
        existing.reachable_roots .|= edge.reachable_roots
        existing.reachable_vars .|= edge.reachable_vars

        dg.dirty_roots .|= edge.reachable_roots
        dg.dirty_vars .|= edge.reachable_vars

        return nothing
    end

    push!(dg.child_edges[top_vertex], edge)
    push!(dg.parent_edges[bott_vertex], edge)

    # adding an edge can only affect doms/pdoms for its roots+vars, so mark them as dirty for recomputation
    dg.dirty_roots .|= edge.reachable_roots
    dg.dirty_vars .|= edge.reachable_vars

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

is_root_reachable(dg::DerivativeGraph{T}, node::T, root::Integer) where {T} = haskey(dg.root_idx_to_postorder, root) && (dg.root_idx_to_postorder[root] == node || any(e -> reachable_roots(e)[root], parent_edges(dg, node)))

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

is_var_reachable(dg::DerivativeGraph{T}, node::T, var::Integer) where {T} = dg.var_idx_to_postorder[var] == node || any(e -> reachable_vars(e)[var], child_edges(dg, node))

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
        _ => begin
            der = derivative_idx(expr, arg_idx)
            @assert !isnothing(der) """
            Unable to compute derivative of $expr w.r.t. argument $arg_idx.
            If this is a user-registered function, make sure to register its derivatives with `@register_derivative`.
            """
            return der
        end
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

        isnothing(post_idx) && continue

        dg.root_idx_to_postorder[root_idx] = post_idx
        dg.postorder_to_root_idx[post_idx] = root_idx
    end
end

function populate_dergraph!(dg::DerivativeGraph{T}, expr::SymbolicT, root_idx::Integer) where {T}
    haskey(dg.definitions, expr) && return populate_root_reachabilities!(dg, dg.definitions[expr], root_idx)

    !iscall(expr) && return nothing

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
        # figure out reachable vars from child edges
        arg_reachable_vars = reachable_vars(dg, arg_post_idx)
        # new edge, so only reachable by the given root
        arg_reachable_roots = falses(length(dg.roots))
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
    var_idxs = findall(isequal(var), dg.vars)
    for var_idx in var_idxs
        dg.var_idx_to_postorder[var_idx] = post_idx
    end
    dg.postorder_to_var_idx[post_idx] = first(var_idxs)
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
function get_dominators(dg::DerivativeGraph{T}, root::Integer) where {T}
    doms = Vector{Union{Nothing, T}}(undef, length(dg))
    root_idxs = collect(values(root_postorders(dg)))
    doms[root_idxs] = root_idxs

    # moves two nodes up the graph until they meet
    function get_common_parent(a::T, b::T)::Union{Nothing, T}
        # move a and b up the graph through their immediate dominators until they meet
        (isnothing(doms[a]) || isnothing(doms[b])) && return nothing
        while a != b
            !(a < b && a != doms[a]) && !(b < a && b != doms[b]) && return nothing
            while !isnothing(a) && a < b && a != doms[a]
                a = doms[a]
            end
            isnothing(a) && return nothing
            while !isnothing(b) && b < a && b != doms[b]
                b = doms[b]
            end
            isnothing(b) && return nothing
        end
        return a
    end

    changed = true # keeps track of when changes stop happening
    while changed
        changed = false
        for node in reverse(eachindex(dg))
            # skip over nodes not reachable from root
            if !is_root_reachable(dg, node, root)
                doms[node] = nothing
                continue
            end

            # filter parents by root
            parents = [top_vertex(e) for e in parent_edges(dg, node) if reachable_roots(e)[root]]
            
            if isempty(parents)
                doms[node] = node
                continue
            end
            
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

function get_postdominators(dg::DerivativeGraph{T}, var::Integer) where {T}
    pdoms = Vector{Union{Nothing, T}}(undef, length(dg))
    var_idxs = collect(values(dg.var_idx_to_postorder))
    pdoms[var_idxs] = var_idxs

    function get_common_child(a::T, b::T)::Union{Nothing, T}
        # move a and b up the graph through their immediate dominators until they meet
        (isnothing(pdoms[a]) || isnothing(pdoms[b])) && return nothing
        while a != b
            !(a > b && a != pdoms[a]) && !(b > a && b != pdoms[b]) && return nothing
            while !isnothing(a) && a > b && a != pdoms[a]
                a = pdoms[a]
            end
            isnothing(a) && return nothing
            while !isnothing(b) && b > a && b != pdoms[b]
                b = pdoms[b]
            end
            isnothing(b) && return nothing
        end
        return a
    end

    changed = true # keeps track of when changes stop happening
    while changed
        changed = false
        for node in eachindex(dg)
            if !is_var_reachable(dg, node, var)
                pdoms[node] = nothing
                continue
            end
            
            children = [bott_vertex(e) for e in child_edges(dg, node) if reachable_vars(e)[var]]

            if isempty(children)
                pdoms[node] = node
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
    dominance_mask::BitVector
    edges::Set{Edge{T}}
    dg::DerivativeGraph{T}
    times_used::Int

    function FactorableSubgraph{T, DominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector, dominance_mask::BitVector, dg::DerivativeGraph{T}) where {T<:Integer}
        new{T, DominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, dominance_mask, Set{Edge{T}}(), dg, sum(dominance_mask)*sum(reachable_vars))
    end

    function FactorableSubgraph{T, PostDominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector, dominance_mask::BitVector, dg::DerivativeGraph{T}) where {T<:Integer}
        new{T, PostDominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, dominance_mask, Set{Edge{T}}(), dg, sum(dominance_mask)*sum(reachable_roots))
    end
end

Base.show(io::IO, sub::FactorableSubgraph{T,S}) where {T,S} = print(io, "$S($(sub.top_vertex), $(sub.bott_vertex))")

Base.:(==)(::FactorableSubgraph{T, DominatorSubgraph}, ::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = false
Base.:(==)(::FactorableSubgraph{T, PostDominatorSubgraph}, ::FactorableSubgraph{T, DominatorSubgraph}) where {T} = false

Base.:(==)(a::FactorableSubgraph{T, S}, b::FactorableSubgraph{T, S}) where {T, S} = a.top_vertex == b.top_vertex && a.bott_vertex == b.bott_vertex
Base.hash(e::FactorableSubgraph{T, DominatorSubgraph}, h::UInt) where {T} = hash((e.top_vertex, e.bott_vertex, 0), h)
Base.hash(e::FactorableSubgraph{T, PostDominatorSubgraph}, h::UInt) where {T} = hash((e.top_vertex, e.bott_vertex, 1), h)

# Loops through the immediate dominators to fill out the dominance masks. Also works for postdominators.
function calculate_dominance_mask(dominators::Vector{T}) where {T}
    dom_mask = Vector{BitVector}(undef, length(dominators))
    for i in eachindex(dominators)
        dom_mask[i] = falses(length(dominators))
    end
    for (node, dom) in enumerate(dominators)
        # walk up through the dominators
        while !isnothing(dom)
            dom_mask[dom][node] = 1
            dom == dominators[dom] && break
            dom = dominators[dom]
        end
    end

    return dom_mask
end

# the following functions allow functions to treat Dominator and PostDominator subgraphs the same by using forward and backward instead of up and down the graph
# forward is in the direction of dominated to dominating or postdominated to postdominating. i.e. from the factor base to the factor node
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

# analogous to forward; dominance is the root masks for dominator subgraphs and variable masks for postdominator subgraphs (validated in get_factorable_subgraphs)
# analogous to backward; nondominance is the var reachability for dominator subgraphs and root reachability for postdominator subgraphs
# not the same thing as dominator/postdominator but related
dominance_mask(::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = reachable_roots(edge)
dominance_mask(::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = reachable_vars(edge)

nondominance_mask(::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = reachable_vars(edge)
nondominance_mask(::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = reachable_roots(edge)
nondominance_mask(sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = sub.reachable_vars
nondominance_mask(sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = sub.reachable_roots

# copies correct masks for constructing a new subgraph edge
sub_edge_reachable_roots(sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = copy(sub.dominance_mask)
sub_edge_reachable_roots(sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = copy(nondominance_mask(sub))
sub_edge_reachable_vars(sub::FactorableSubgraph{T, DominatorSubgraph}) where {T} = copy(nondominance_mask(sub))
sub_edge_reachable_vars(sub::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = copy(sub.dominance_mask)

subgraph_edges(sub::FactorableSubgraph) = sub.edges

# NOT called in the constructor of `FactorableSubgraph`. Instead, delayed until right before factoring to account for changes to the `DerivativeGraph` from previous factoring
# Also calculates the subgraph's value
# Recurses through the subgraph to find all its edges
function populate_subgraph_edges!(dg::DerivativeGraph{T}, sub::FactorableSubgraph) where {T}
    union!(sub.edges, _subgraph_edges(dg, sub, backward_vertex(sub)))
end

# finds edges with overlapping nondominance reachabilities for combining paths
# this prevents combining the same relationship twice, such as where one branch is from a previous factored subgraph and the others are part of a subgraph of a different type (dom/pdom) being factored between the same nodes
function find_edge_group(sub::FactorableSubgraph, edges::Vector{Edge{T}}) where {T}
    isempty(edges) && return BitVector()
    # start with the first edge in the group
    mask = copy(nondominance_mask(sub, edges[1]))
    group::BitVector = zeros(length(edges))
    group[1] = 1

    # loop through edges until the group is unchanged to prevent a scenario where the nondominance looks like {1}, {2,3}, {1,3}, where one loop wouldn't include edge 2
    changed = true
    while changed
        changed = false
        for (i,edge) in pairs(edges)
            group[i] && continue # edge not already in group
            nondom_mask = nondominance_mask(sub, edge)
            if any(nondom_mask .& mask) # any overlap
                changed = true
                mask .|= nondom_mask
                group[i] = 1
            end
        end
    end

    return group
end

function _subgraph_edges(dg::DerivativeGraph{T}, sub::FactorableSubgraph, node::T) where {T}
    cache = Dict{Edge{T}, SymbolicT}()
    sub_edges = Set{Edge{T}}()

    # discover all edges, but only sum those part of a group of reachabilities (explained in find_edge_group)
    edges = forward_edges(dg, sub, node)
    edge_group_mask = find_edge_group(sub, edges)
    for (in_group,edge) in zip(edge_group_mask,edges)
        sub.subgraph_value += _subgraph_edges!(sub_edges, dg, sub, edge, cache)*in_group
    end

    return sub_edges
end


function _subgraph_edges!(edges::Set{Edge{T}}, dg::DerivativeGraph{T}, sub::FactorableSubgraph, edge::Edge{T}, cache::Dict{Edge{T}, SymbolicT}) where {T}
    haskey(cache, edge) && return cache[edge]

    edge_dom_mask = dominance_mask(sub, edge)
    if !any(edge_dom_mask .& sub.dominance_mask) # edge not in subgraph (no overlap)
        cache[edge] = COMMON_ZERO
        return COMMON_ZERO
    end

    # reached the end of the subgraph
    if forward_vertex(sub, edge) == forward_vertex(sub)
        push!(edges, edge)
        cache[edge] = edge.edge_value
        return edge.edge_value
    end

    # discover all edges, but only sum those part of a group of reachabilities (explained in find_edge_group)
    push!(edges, edge)
    path_value = COMMON_ZERO
    next_edges = forward_edges(dg, sub, edge)
    edge_group_mask = find_edge_group(sub, next_edges)
    for (in_group,next_edge) in zip(edge_group_mask,next_edges)
        path_value += _subgraph_edges!(edges, dg, sub, next_edge, cache)*in_group
    end

    path_value *= edge.edge_value

    cache[edge] = path_value
    return path_value
end

# the number of times a subgraph is used by all possible partial derivatives
subgraph_count(sub::FactorableSubgraph) = sub.times_used

"""
    get_factorable_subgraphs(dg::DerivativeGraph{T};
        dom_cache=Dict{Int, Vector{Union{Nothing,T}}}(),
        pdom_cache=Dict{Int, Vector{Union{Nothing,T}}}()) where {T} -> BinaryHeap

Generates a heap of factorable subgraphs ordered with `FactorOrder` (smallest + most used are factored first).

`dom_cache`/`pdom_cache` map root/var index -> the `get_dominators`/`get_postdominators` result for that
root/var. When supplied (and reused across repeated calls, e.g. from `factor_subgraphs!`), a root/var's
dominators are only recomputed if `dg.dirty_roots`/`dg.dirty_vars` marks it as changed since it was last
cached
"""
function get_factorable_subgraphs(dg::DerivativeGraph{T};
        dom_cache::Dict{Int, Vector{Union{Nothing,T}}}=Dict{Int, Vector{Union{Nothing,T}}}(),
        pdom_cache::Dict{Int, Vector{Union{Nothing,T}}}=Dict{Int, Vector{Union{Nothing,T}}}()) where {T}
    subs = DataStructures.BinaryHeap{Union{FactorableSubgraph{T, DominatorSubgraph}, FactorableSubgraph{T, PostDominatorSubgraph}}, FactorOrder}()
    dom_pairs = Dict{Tuple{T,T}, BitVector}() # maps (dominated, dominating) pairs to the bitmask of all roots that reach the pair
    for root in keys(dg.root_idx_to_postorder)
        # only recompute doms if root is dirty
        (dg.dirty_roots[root] || !haskey(dom_cache, root)) && (dom_cache[root] = get_dominators(dg, root))

        doms = dom_cache[root]
        for (dominated, dominating) in pairs(doms)
            # check dominated node is a factor base (2+ parents)
            isnothing(dominating) && continue
            count(e -> reachable_roots(e)[root], parent_edges(dg, T(dominated))) > 1 || continue
            
            dom_pair = (dominated, dominating)
            if !haskey(dom_pairs, dom_pair)
                dom_pairs[dom_pair] = zeros(length(dg.roots))
            end
            dom_pairs[dom_pair][root] = 1
        end
    end

    # create subgraphs for each pair
    for ((dominated, dominating), root_mask) in dom_pairs
        reachable_vars_mask = reachable_vars(dg, T(dominated))
        reachable_roots_mask = reachable_roots(dg, dominating)
        push!(subs, FactorableSubgraph{T, DominatorSubgraph}(dominating, T(dominated), reachable_vars_mask, reachable_roots_mask, root_mask, dg))
    end

    # repeat the same process, but for postdominators and variables
    pdom_pairs = Dict{Tuple{T,T}, BitVector}()
    for var in keys(dg.var_idx_to_postorder)
        # only recompute pdoms if var is dirty
        (dg.dirty_vars[var] || !haskey(pdom_cache, var)) && (pdom_cache[var] = get_postdominators(dg, var))

        pdoms = pdom_cache[var]
        for (postdominated, postdominating) in pairs(pdoms)
            # check postdominated node is a factor base (2+ children)
            isnothing(postdominating) && continue
            count(e -> reachable_vars(e)[var], child_edges(dg, T(postdominated))) > 1 || continue

            pdom_pair = (postdominated, postdominating)
            if !haskey(pdom_pairs, pdom_pair)
                pdom_pairs[pdom_pair] = zeros(length(dg.vars))
            end
            pdom_pairs[pdom_pair][var] = 1
        end
    end

    for ((postdominated, postdominating), var_mask) in pdom_pairs
        reachable_vars_mask = reachable_vars(dg, postdominating)
        reachable_roots_mask = reachable_roots(dg, T(postdominated))
        push!(subs, FactorableSubgraph{T, PostDominatorSubgraph}(T(postdominated), postdominating, reachable_vars_mask, reachable_roots_mask, var_mask, dg))
    end

    # dom_cache/pdom_cache are now up to date with all roots/vars that were dirty coming in
    fill!(dg.dirty_roots, false)
    fill!(dg.dirty_vars, false)

    return subs
end

# checks for structural true dominance, masked by variable (used for pdom subgraphs)
# TODO: optimize
function is_dominator(dg::DerivativeGraph{T}, dominator::T, dominated::T, var_mask::BitVector, cache::Dict{Tuple{T,T}, Bool}=Dict{Tuple{T,T},Bool}()) where {T}
    dominator == dominated && return true
    cache_key = (dominator, dominated)
    haskey(cache, cache_key) && return cache[cache_key]
    cache[cache_key] = false

    next_edges = filter(e -> any(reachable_vars(e) .& var_mask), parent_edges(dg, dominated))

    isempty(next_edges) && return false
    for parent_edge in next_edges
        is_dominator(dg, dominator, top_vertex(parent_edge), var_mask, cache) || return false
    end

    cache[cache_key] = true
    return true
end

function is_postdominator(dg::DerivativeGraph{T}, postdominator::T, postdominated::T, root_mask::BitVector, cache::Dict{Tuple{T,T}, Bool}=Dict{Tuple{T,T},Bool}()) where {T}
    postdominator == postdominated && return true
    cache_key = (postdominator, postdominated)
    haskey(cache, cache_key) && return cache[cache_key]
    cache[cache_key] = false

    next_edges = filter(e -> any(reachable_roots(e) .& root_mask), child_edges(dg, postdominated))

    isempty(next_edges) && return false
    for child_edge in next_edges
        is_postdominator(dg, postdominator, bott_vertex(child_edge), root_mask, cache) || return false
    end

    cache[cache_key] = true
    return true
end

check_dominance(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T} = is_postdominator(dg, sub.bott_vertex, edge.bott_vertex, sub.dominance_mask)
check_dominance(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T} = is_dominator(dg, sub.top_vertex, edge.top_vertex, sub.dominance_mask)

# factors a subgraph from dg, replacing it with a single edge (keeping original edges when necessary)
function factor_subgraph!(dg::DerivativeGraph{T}, sub::FactorableSubgraph) where {T}
    # check that the factor and factor base nodes are still a factor and factor base
    (length(backward_edges(dg, sub, forward_vertex(sub))) < 2 || length(forward_edges(dg, sub, backward_vertex(sub))) < 2) && return false

    populate_subgraph_edges!(dg, sub)
    sub_edges = subgraph_edges(sub)

    # TODO: cut this down
    extras = Dict{Edge{T}, Tuple{BitVector, BitVector}}()
    for edge in sub_edges
        dom_extra = dominance_mask(sub, edge) .& .~sub.dominance_mask
        nondom_extra = nondominance_mask(sub, edge) .& .~nondominance_mask(sub)

        if !any(dom_extra) && !any(nondom_extra)
            check_dominance(dg, sub, edge) || return false
        end

        extras[edge] = (dom_extra, nondom_extra)
    end

    for edge in sub_edges
        # for comparison to determine dirty roots+vars
        old_roots = copy(edge.reachable_roots)
        old_vars = copy(edge.reachable_vars)

        dom_extra, nondom_extra = extras[edge]
        if !any(dom_extra) && !any(nondom_extra)
            # edge is completely contained in subgraph
            rem_edge!(dg, edge)
        elseif any(dom_extra) && !any(nondom_extra)
            # connected to another root/var, so just narrow scope to that root/var
            dominance_mask(sub, edge) .= dom_extra
        elseif any(nondom_extra) && !any(dom_extra)
            # connected to another root/var, so just narrow scope to that root/var
            nondominance_mask(sub, edge) .= nondom_extra
        end

        dg.dirty_roots .|= old_roots .!= edge.reachable_roots
        dg.dirty_vars .|= old_vars .!= edge.reachable_vars
    end

    # add new subgraph edge, using the dominance_mask for roots/vars as applicable
    sub_edge = Edge{T}(sub.subgraph_value, sub.top_vertex, sub.bott_vertex, sub_edge_reachable_vars(sub), sub_edge_reachable_roots(sub))
    add_edge!(dg, sub_edge)

    # propagate changes to reachability only after fully factoring subgraph
    propogate_root_reachability(dg, sub_edge)
    propogate_var_reachability(dg, sub_edge)

    for edge in sub_edges
        propogate_root_reachability(dg, edge)
        propogate_var_reachability(dg, edge)
    end
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
function factor_subgraphs!(dg::DerivativeGraph{T}) where {T}
    dom_cache = Dict{Int, Vector{Union{Nothing,T}}}()
    pdom_cache = Dict{Int, Vector{Union{Nothing,T}}}()
    subs = get_factorable_subgraphs(dg; dom_cache, pdom_cache)

    # subgraphs that factor_subgraph! refused to factor (ambiguous classification) so we
    # don't re-propose and re-refuse the same (top_vertex, bott_vertex, kind) forever
    rejected = Set{FactorableSubgraph}()

    while !isempty(subs)
        # factor the first subgraph according to `FactorOrder``
        sub = pop!(subs)
        sub in rejected && continue

        factor_subgraph!(dg, sub) || push!(rejected, sub)
        subs = get_factorable_subgraphs(dg; dom_cache, pdom_cache)
    end
end

# evaluate the derivative of root w.r.t. var using a fully factored DerivativeGraph
function evaluate_path(dg::DerivativeGraph, root::Integer, var::Integer, cache::Vector{Dict{Edge,SymbolicT}})
    haskey(dg.root_idx_to_postorder, root) || return COMMON_ZERO
    haskey(dg.var_idx_to_postorder, var) || return COMMON_ZERO
    root_postorder = dg.root_idx_to_postorder[root]
    var_postorder = dg.var_idx_to_postorder[var]

    root_postorder == var_postorder && return COMMON_ONE

    next_edges = filter(e -> reachable_vars(e)[var], dg.child_edges[root_postorder])
    isempty(next_edges) && return COMMON_ZERO # path from root to var does not exist
    @assert length(next_edges) == 1 "Error in graph factoring. There is >1 path from root to var."

    return evaluate_path(dg, first(next_edges), var, cache)
end

function evaluate_path(dg::DerivativeGraph, edge::Edge, var::Integer, cache::Vector{Dict{Edge,SymbolicT}})
    edge.bott_vertex == dg.var_idx_to_postorder[var] && return edge.edge_value # reached var
    haskey(cache[var], edge) && return cache[var][edge]

    next_edges = filter(e -> reachable_vars(e)[var], dg.child_edges[edge.bott_vertex])
    isempty(next_edges) && return COMMON_ZERO # path from root to var does not exist
    @assert length(next_edges) == 1 "Error in graph factoring. There is >1 path from root to var."

    result = evaluate_path(dg, first(next_edges), var, cache) * edge.edge_value
    cache[var][edge] = result

    return result
end

"""
$(SIGNATURES)

Computes the Jacobian of `roots` w.r.t. `vars` using the D* automatic differentiation algorithm using the [`DerivativeGraph`](@ref) data structure.

(see [this paper](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/main-65.pdf) for more details on the algorithm)

Mostly the same usage as [`jacobian`](@ref). More limited in input expressions (doesn't support `ifelse` or nested differentials), but asymptotically faster for large Rn->Rm expressions.

# Arguments

- `roots::AbstractVector`: Vector of expressions to differentate or array-type symbolic expression (e.g. function registered with `@register_array_symbolic`)
- `vars::AbstractVector`: Vector of variables to differentate w.r.t. or single array-type variable (e.g. `@variables x[1:4]`)
"""
function dstar_jacobian(roots::AbstractVector, vars::AbstractVector{SymbolicT})
    roots isa Arr && (roots = scalarize(unwrap(roots)))
    roots isa AbstractVector{Num} && (roots = unwrap.(roots))
    dg = DerivativeGraph(roots, vars)
    factor_subgraphs!(dg)

    result = Matrix{SymbolicT}(undef, length(roots), length(vars))
    cache = [Dict{Edge,SymbolicT}() for _ in eachindex(vars)]

    for root in eachindex(roots)
        for var in eachindex(vars)
            result[root, var] = evaluate_path(dg, root, var, cache)
        end
    end

    return result
end

function dstar_jacobian(roots, vars)
    # input validation copied from `jacobian`
    roots = vec(scalarize(roots))
    if roots isa Vector{Num}
        roots = unwrap.(roots)::Vector{SymbolicT}
    elseif roots isa Vector{SymbolicT}
    else
        roots = roots::Vector{eltype(roots)}
    end
    # Suboptimal, but prevents wrong results on Arr for now. Arr resulting from a symbolic function will fail on this due to unknown size.
    vars = vec(scalarize(vars))
    if vars isa Vector{Num}
        vars = unwrap.(vars)::Vector{SymbolicT}
    elseif vars isa Vector{SymbolicT}
    else
        error("This should not happen! `vars` must be convertible to Vector{SymbolicT}. \nReceived vars = $vars")
    end
    return dstar_jacobian(roots, vars)
end

"""
$(SIGNATURES)

Computes the derivative of `root` w.r.t. `var` using the D* differentiation algorithm.

Mostly the same usage as [`derivative`](@ref), but more limited in input expressions (doesn't support `ifelse` or nested differentials).

Wrapper for R1->R1 case of `dstar_jacobian`. See [`dstar_jacobian`](@ref) for more information.

# Arguments
- `root`: Expression to differentate
- `var`: Variable to differentate w.r.t.
"""
dstar_derivative(root::Union{Num,SymbolicT}, var::Union{Num,SymbolicT}) = only(dstar_jacobian([root], [var]))
