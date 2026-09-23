abstract type AbstractFactorableSubgraph end
abstract type DominatorSubgraph <: AbstractFactorableSubgraph end
abstract type PostDominatorSubgraph <: AbstractFactorableSubgraph end

# `mutable` so edges have object identity: two edges between the same vertices
# represent different path sets (e.g. a factored subgraph edge alongside a
# preexisting edge) and must be distinguishable, like FastDifferentiation's `PathEdge`
mutable struct Edge{T<:Integer}
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

# scratch buffers reused across graph traversals (subgraph reachability walks,
# factor-base bypass checks, shared-path detection, subgraph edge collection) to
# avoid allocating per call. Buffers are grown lazily to `length(dg.symbols)`.
struct DGScratch{T<:Integer}
    seen::Vector{BitVector} # two reachability seen-sets (forward and backward walks are live simultaneously)
    stack::Vector{T} # DFS stack
    pairs::Set{Tuple{Int,Int}} # has_shared_path pair set
    counts::Dict{T,T} # populate_subgraph_edges! in-degree counts
    vterms::Dict{T,Vector{SymbolicT}} # populate_subgraph_edges! per-node sum terms
end

DGScratch{T}() where {T} = DGScratch{T}(BitVector[BitVector(), BitVector(), BitVector()], T[],
    Set{Tuple{Int,Int}}(), Dict{T,T}(), Dict{T,Vector{SymbolicT}}())

# grow `buf` to at least `n` bits and reset it to all false
function _seen_buf!(buf::BitVector, n::Integer)
    length(buf) < n && resize!(buf, n)
    fill!(buf, false)
    return buf
end

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
    scratch::DGScratch{T} # reusable traversal buffers
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
        trues(length(vars)),
        DGScratch{idx_type}()
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

function hasedge(dg::DerivativeGraph{T}, edge::Edge{T}) where {T}
    is_child_edge = edge in child_edges(dg, top_vertex(edge))
    is_parent_edge = edge in parent_edges(dg, bott_vertex(edge))

    @assert is_child_edge == is_parent_edge "edge is only in one of the child/parent edge lists"

    return is_child_edge && is_parent_edge
end

function propagate_var_reachability(dg::DerivativeGraph{T}, node::T) where {T}
    _parent_edges = parent_edges(dg, node)
    isempty(_parent_edges) && return nothing

    new_reachability = reachable_vars(dg, node) # computes reachable_vars from union of child edge reachabilities

    # shrink parent reachabilities
    for edge in _parent_edges
        old_reachability = copy(edge.reachable_vars)
        edge.reachable_vars .&= new_reachability
        if edge.reachable_vars != old_reachability
            dg.dirty_vars .|= old_reachability .& .~edge.reachable_vars
            propagate_var_reachability(dg, top_vertex(edge))
        end
    end
end
propagate_var_reachability(dg::DerivativeGraph{T}, edge::Edge{T}) where {T} = propagate_var_reachability(dg, edge.top_vertex)

function propagate_root_reachability(dg::DerivativeGraph{T}, node::T) where {T}
    _child_edges = child_edges(dg, node)
    isempty(_child_edges) && return nothing

    new_reachability = reachable_roots(dg, node) # computes reachable_roots from union of parent edge reachabilities

    # shrink child reachabilities
    for edge in _child_edges
        old_reachability = copy(edge.reachable_roots)
        edge.reachable_roots .&= new_reachability
        if edge.reachable_roots != old_reachability
            dg.dirty_roots .|= old_reachability .& .~edge.reachable_roots
            propagate_root_reachability(dg, bott_vertex(edge))
        end
    end
end

propagate_root_reachability(dg::DerivativeGraph{T}, edge::Edge{T}) where {T} = propagate_root_reachability(dg, edge.bott_vertex)

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

    # parallel edges between the same vertices are allowed: they are distinct
    # objects covering disjoint sets of paths, so never merge them
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
                # `arguments` materializes `dict` entries in order, preceded by
                # `coeff` when it isn't 1; drop this argument's dict entry and
                # rebuild via the fast `Mul` ctor rather than re-canonicalizing
                # all factors through `mul_worker`
                args = parent(arguments(expr))
                key_pos = arg_idx - (length(args) - length(dict))
                newdict = copy(dict)
                for (i, k) in enumerate(keys(dict))
                    i == key_pos && (delete!(newdict, k); break)
                end
                return SymbolicUtils.Mul{VartypeT}(coeff, newdict; type = symtype(expr), shape = shape(expr))
            end
        end
        _ => begin
            der = derivative_idx(expr, arg_idx)
            isnothing(der) && throw(DerivativeNotDefinedError(expr, arg_idx))
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
    # `ifelse` conditions are treated as piecewise-constant, matching
    # `expand_derivatives` (`D(ifelse(c,a,b)) == ifelse(c,D(a),D(b))`). The
    # condition is excluded from the graph entirely so non-differentiable
    # subterms (comparisons) are never traversed.
    op = operation(expr)
    cond_idx = op === ifelse || op === ifelse_eager || op === ifelse_branching ? 1 : 0
    for arg_idx in reverse(eachindex(args))
        arg = args[arg_idx]
        if arg_idx == cond_idx
            arg_idx_to_post_idx[arg_idx] = T(-1)
        elseif arg in dg.varset
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

    # add new edges; partial derivatives of identical arguments are summed into a
    # single edge (the total derivative w.r.t. that argument)
    partial_ders = Dict{T, SymbolicT}()
    reachable_masks = Dict{T, BitVector}()
    for (arg_idx, arg_post_idx) in enumerate(arg_idx_to_post_idx)
        arg_post_idx == T(-1) && continue
        arg_reachable_vars = get!(() -> reachable_vars(dg, arg_post_idx), reachable_masks, arg_post_idx)
        if any(arg_reachable_vars)
            existing = get(partial_ders, arg_post_idx, nothing)
            partial_ders[arg_post_idx] = isnothing(existing) ? nary_derivative_idx(expr, arg_idx) : existing + nary_derivative_idx(expr, arg_idx)
        else
            # the edge can never reach a var, so its value can never pass a
            # mask check and is never read; don't bother computing the partial
            partial_ders[arg_post_idx] = COMMON_ZERO
        end
    end
    for (arg_post_idx, partial_der) in partial_ders
        # figure out reachable vars from child edges
        arg_reachable_vars = reachable_masks[arg_post_idx]
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

# Follows the algorithm described in this paper: https://www.cs.tufts.edu/comp/150FP/archive/keith-cooper/dom14.pdf
function get_dominators(dg::DerivativeGraph{T}, root::Integer) where {T}
    doms = Vector{Union{Nothing, T}}(undef, length(dg))
    fill!(doms, nothing)
    for ri in values(root_postorders(dg))
        doms[ri] = ri
    end
    # roots whose expressions aren't calls have no postorder entry
    haskey(dg.root_idx_to_postorder, root) || return doms

    # nodes reachable from `root`, computed once up front rather than via
    # is_root_reachable inside the fixpoint loop
    reach = _seen_buf!(dg.scratch.seen[3], length(dg.symbols))
    stack = empty!(dg.scratch.stack)
    root_node = dg.root_idx_to_postorder[root]
    reach[root_node] = true
    push!(stack, root_node)
    while !isempty(stack)
        node = pop!(stack)
        for e in child_edges(dg, node)
            reachable_roots(e)[root] || continue
            b = bott_vertex(e)
            reach[b] && continue
            reach[b] = true
            push!(stack, b)
        end
    end

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
            if !reach[node]
                doms[node] = nothing
                continue
            end

            # intersect over parents reachable from root, without allocating
            new_idom::Union{Nothing, T} = nothing
            first_parent = true
            for e in parent_edges(dg, node)
                reachable_roots(e)[root] || continue
                parent = top_vertex(e)
                if first_parent
                    new_idom = parent
                    first_parent = false
                elseif isassigned(doms, parent)
                    new_idom = get_common_parent(parent, new_idom)
                    isnothing(new_idom) && break
                end
            end

            if first_parent
                doms[node] = node
                continue
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
    for vi in values(dg.var_idx_to_postorder)
        pdoms[vi] = vi
    end

    # nodes that can reach `var`, computed once up front rather than via
    # is_var_reachable inside the fixpoint loop
    reach = _seen_buf!(dg.scratch.seen[3], length(dg.symbols))
    stack = empty!(dg.scratch.stack)
    var_node = dg.var_idx_to_postorder[var]
    reach[var_node] = true
    push!(stack, var_node)
    while !isempty(stack)
        node = pop!(stack)
        for e in parent_edges(dg, node)
            reachable_vars(e)[var] || continue
            t = top_vertex(e)
            reach[t] && continue
            reach[t] = true
            push!(stack, t)
        end
    end

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
            if !reach[node]
                pdoms[node] = nothing
                continue
            end

            # intersect over children that can reach var, without allocating
            new_pidom::Union{Nothing, T} = nothing
            first_child = true
            for e in child_edges(dg, node)
                reachable_vars(e)[var] || continue
                child = bott_vertex(e)
                if first_child
                    new_pidom = child
                    first_child = false
                elseif isassigned(pdoms, child)
                    new_pidom = get_common_child(child, new_pidom)
                    isnothing(new_pidom) && break
                end
            end

            if first_child
                pdoms[node] = node
                continue
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
    edges::OrderedCollections.OrderedSet{Edge{T}}
    dg::DerivativeGraph{T}
    times_used::Int

    function FactorableSubgraph{T, DominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector, dominance_mask::BitVector, dg::DerivativeGraph{T}) where {T<:Integer}
        new{T, DominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, dominance_mask, OrderedCollections.OrderedSet{Edge{T}}(), dg, sum(dominance_mask)*sum(reachable_vars))
    end

    function FactorableSubgraph{T, PostDominatorSubgraph}(top_vertex::T, bott_vertex::T, reachable_vars::BitVector, reachable_roots::BitVector, dominance_mask::BitVector, dg::DerivativeGraph{T}) where {T<:Integer}
        new{T, PostDominatorSubgraph}(COMMON_ZERO, top_vertex, bott_vertex, reachable_vars, reachable_roots, dominance_mask, OrderedCollections.OrderedSet{Edge{T}}(), dg, sum(dominance_mask)*sum(reachable_roots))
    end
end

Base.show(io::IO, sub::FactorableSubgraph{T,S}) where {T,S} = print(io, "$S($(sub.top_vertex), $(sub.bott_vertex))")

Base.:(==)(::FactorableSubgraph{T, DominatorSubgraph}, ::FactorableSubgraph{T, PostDominatorSubgraph}) where {T} = false
Base.:(==)(::FactorableSubgraph{T, PostDominatorSubgraph}, ::FactorableSubgraph{T, DominatorSubgraph}) where {T} = false

Base.:(==)(a::FactorableSubgraph{T, S}, b::FactorableSubgraph{T, S}) where {T, S} = a.top_vertex == b.top_vertex && a.bott_vertex == b.bott_vertex
Base.hash(e::FactorableSubgraph{T, DominatorSubgraph}, h::UInt) where {T} = hash((e.top_vertex, e.bott_vertex, 0), h)
Base.hash(e::FactorableSubgraph{T, PostDominatorSubgraph}, h::UInt) where {T} = hash((e.top_vertex, e.bott_vertex, 1), h)

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

# When factoring a subgraph, an edge whose reachability extends outside the
# subgraph in the dominance direction must keep serving those outside paths.
# `outside_edge` creates a parallel edge carrying `dom_extra` (the dominance
# reachability outside the subgraph) together with the edge's full
# nondominance reachability; the original edge is then narrowed to the
# in-subgraph dominance part and the outside nondominance part. This is the
# analogue of FastDifferentiation's `add_non_dom_edges!`/`reset_edge_masks!`.
outside_edge(sub::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}, dom_extra::BitVector) where {T} =
    Edge{T}(edge.edge_value, edge.top_vertex, edge.bott_vertex, copy(edge.reachable_vars), dom_extra)
outside_edge(sub::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}, dom_extra::BitVector) where {T} =
    Edge{T}(edge.edge_value, edge.top_vertex, edge.bott_vertex, dom_extra, copy(edge.reachable_roots))

subgraph_edges(sub::FactorableSubgraph) = sub.edges

# all nodes reachable from `start` moving through `sub` along in-subgraph edges
# only; `forward` chooses the factor-base-to-factor-node direction
function _subgraph_reachable(dg::DerivativeGraph{T}, sub::FactorableSubgraph, start::T, forward::Bool) where {T}
    seen = _seen_buf!(dg.scratch.seen[forward ? 1 : 2], length(dg.symbols))
    seen[start] = true
    stack = empty!(dg.scratch.stack)
    push!(stack, start)
    while !isempty(stack)
        node = pop!(stack)
        for e in (forward ? forward_edges(dg, sub, node) : backward_edges(dg, sub, node))
            test_edge(sub, e) || continue
            next = forward ? forward_vertex(sub, e) : backward_vertex(sub, e)
            seen[next] && continue
            seen[next] = true
            push!(stack, next)
        end
    end
    return seen
end

# accumulate `sum` into `node`'s path-product terms and propagate their total
# forward once every in-subgraph backward edge has contributed (cf.
# FastDifferentiation's `_evaluate_branching_subgraph`). Terms are collected per
# node and summed once with `add_worker` rather than building an intermediate
# `Add` per contributing edge.
function _vertex_sum!(dg::DerivativeGraph{T}, sub::FactorableSubgraph, sum::SymbolicT, node::T, counts::Dict{T,T}, vertex_terms::Dict{T,Vector{SymbolicT}}) where {T}
    push!(get!(()->SymbolicT[], vertex_terms, node), sum)
    (counts[node] -= 1) == 0 || return
    total = SymbolicUtils.add_worker(VartypeT, vertex_terms[node])
    for e in forward_edges(dg, sub, node)
        e in sub.edges || continue
        _vertex_sum!(dg, sub, total * e.edge_value, forward_vertex(sub, e), counts, vertex_terms)
    end
end

# NOT called in the constructor of `FactorableSubgraph`. Instead, delayed until right before factoring to account for changes to the `DerivativeGraph` from previous factoring
# collects all edges on in-subgraph paths between the factor base and factor
# node, and computes the subgraph's value as the sum of all such path products
# via a vertex-sum DP (correct whether or not in-subgraph paths rejoin)
function populate_subgraph_edges!(dg::DerivativeGraph{T}, sub::FactorableSubgraph) where {T}
    empty!(sub.edges)
    fwd_ok = _subgraph_reachable(dg, sub, backward_vertex(sub), true)
    bwd_ok = _subgraph_reachable(dg, sub, forward_vertex(sub), false)
    for node in T.(findall(fwd_ok))
        for e in forward_edges(dg, sub, node)
            test_edge(sub, e) && bwd_ok[forward_vertex(sub, e)] && push!(sub.edges, e)
        end
    end

    counts = empty!(dg.scratch.counts)
    for e in sub.edges
        node = forward_vertex(sub, e)
        counts[node] = get(counts, node, zero(T)) + one(T)
    end
    counts[backward_vertex(sub)] = one(T)
    vertex_terms = empty!(dg.scratch.vterms)
    _vertex_sum!(dg, sub, COMMON_ONE, backward_vertex(sub), counts, vertex_terms)
    sub.subgraph_value = haskey(vertex_terms, forward_vertex(sub)) ?
        SymbolicUtils.add_worker(VartypeT, vertex_terms[forward_vertex(sub)]) : COMMON_ZERO
    return nothing
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
                dom_pairs[dom_pair] = falses(length(dg.roots))
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
    pdom_pairs = empty!(dom_pairs)
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
                pdom_pairs[pdom_pair] = falses(length(dg.vars))
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

# nondominance-direction reachability of `edge` along paths that bypass the factor
# base: for a dominator subgraph, the vars reachable from the edge's lower
# endpoint without passing through the dominated node (the paper's `b pdom e.1`
# test); for a postdominator subgraph, the roots reachable from the edge's upper
# endpoint without passing through the postdominated node (`b dom e.2`)
function bypass_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, DominatorSubgraph}, edge::Edge{T}) where {T}
    reach = falses(length(dg.vars))
    seen = _seen_buf!(dg.scratch.seen[1], length(dg.symbols))
    stack = empty!(dg.scratch.stack)
    push!(stack, edge.bott_vertex)
    while !isempty(stack)
        node = pop!(stack)
        (node == sub.bott_vertex || seen[node]) && continue
        seen[node] = true
        var_idx = get(dg.postorder_to_var_idx, node, nothing)
        var_idx !== nothing && (reach[var_idx] = true)
        for child_edge in child_edges(dg, node)
            push!(stack, child_edge.bott_vertex)
        end
    end
    return reach
end

function bypass_mask(dg::DerivativeGraph{T}, sub::FactorableSubgraph{T, PostDominatorSubgraph}, edge::Edge{T}) where {T}
    reach = falses(length(dg.roots))
    seen = _seen_buf!(dg.scratch.seen[1], length(dg.symbols))
    stack = empty!(dg.scratch.stack)
    push!(stack, edge.top_vertex)
    while !isempty(stack)
        node = pop!(stack)
        (node == sub.top_vertex || seen[node]) && continue
        seen[node] = true
        root_idx = get(dg.postorder_to_root_idx, node, nothing)
        root_idx !== nothing && (reach[root_idx] = true)
        for parent_edge in parent_edges(dg, node)
            push!(stack, parent_edge.top_vertex)
        end
    end
    return reach
end

# `a ⊆ b` for equal-length BitVectors, without allocating (broadcast `.<=`
# materializes a temporary BitVector and `test_edge` is called in inner loops)
@inline function _mask_subset(a::BitVector, b::BitVector)
    ac, bc = a.chunks, b.chunks
    @inbounds for i in eachindex(ac)
        iszero(ac[i] & ~bc[i]) || return false
    end
    return true
end

# an edge is on a valid path within `sub` iff its reachability covers all of the
# subgraph's dominance and nondominance masks. Every in-subgraph edge must serve
# the full (dominance x nondominance) pair set, or the single factored edge could
# not represent the subgraph's value for all of the pairs it claims.
test_edge(sub::FactorableSubgraph, edge::Edge) =
    _mask_subset(sub.dominance_mask, dominance_mask(sub, edge)) &&
    _mask_subset(nondominance_mask(sub), nondominance_mask(sub, edge))

# the unique valid next edge on an in-subgraph path from `edge` toward the factor
# node, or `nothing` if the path ends or branches
function next_valid_edge(dg::DerivativeGraph, sub::FactorableSubgraph, edge::Edge)
    next = nothing
    for e in forward_edges(dg, sub, edge)
        test_edge(sub, e) || continue
        next === nothing || return nothing
        next = e
    end
    return next
end

# whether `start_edge` lies on an unbroken, unbranched in-subgraph path from the
# factor base to the factor node
function isa_connected_path(dg::DerivativeGraph, sub::FactorableSubgraph, start_edge::Edge)
    test_edge(sub, start_edge) || return false
    edge = start_edge
    while forward_vertex(sub, edge) != forward_vertex(sub)
        edge = next_valid_edge(dg, sub, edge)
        edge === nothing && return false
    end
    return true
end

# whether any two of `edges` share a single (dominance, nondominance) pair.
# Parallel edges with disjoint coverage are one path split across edge objects,
# not distinct paths, and must not be treated as a factorable branch.
function has_shared_path(sub::FactorableSubgraph, edges::Vector{Edge{T}}) where {T}
    seen = empty!(sub.dg.scratch.pairs)
    sub_nondom = nondominance_mask(sub)
    for e in edges
        for d in findall(dominance_mask(sub, e))
            sub.dominance_mask[d] || continue
            for n in findall(nondominance_mask(sub, e))
                sub_nondom[n] || continue
                (d, n) in seen && return true
                push!(seen, (d, n))
            end
        end
    end
    return false
end

# whether `sub` is still a factorable subgraph: prior factoring may have deleted
# or narrowed its edges (cf. FastDifferentiation's `subgraph_exists`)
function subgraph_exists(dg::DerivativeGraph, sub::FactorableSubgraph)
    fwd = forward_edges(dg, sub, backward_vertex(sub))
    bwd = backward_edges(dg, sub, forward_vertex(sub))
    (has_shared_path(sub, fwd) && count(e -> test_edge(sub, e), bwd) >= 2) || return false
    return count(e -> isa_connected_path(dg, sub, e), fwd) >= 2
end

# factors a subgraph from dg, replacing it with a single edge (keeping original edges when necessary)
function factor_subgraph!(dg::DerivativeGraph{T}, sub::FactorableSubgraph) where {T}
    # check that the factor and factor base nodes are still a factor and factor base
    subgraph_exists(dg, sub) || return false

    populate_subgraph_edges!(dg, sub)
    sub_edges = subgraph_edges(sub)

    for edge in sub_edges
        # for comparison to determine dirty roots+vars
        old_roots = copy(edge.reachable_roots)
        old_vars = copy(edge.reachable_vars)

        # dominance-direction reachability outside the subgraph is split off onto a
        # parallel edge (cf. FastDifferentiation's `add_non_dom_edges!`)
        dom_extra = dominance_mask(sub, edge) .& .~sub.dominance_mask
        any(dom_extra) && add_edge!(dg, outside_edge(sub, edge, dom_extra))

        # the original edge keeps its in-subgraph dominance reachability and only
        # the nondominance reachability of paths that bypass the factor base or
        # leave the subgraph entirely (cf. FastDifferentiation's `reset_edge_masks!`)
        dominance_mask(sub, edge) .&= sub.dominance_mask
        nondominance_mask(sub, edge) .&= bypass_mask(dg, sub, edge) .| .~nondominance_mask(sub)

        # once either mask is empty the edge is completely absorbed by the subgraph edge
        (!any(edge.reachable_roots) || !any(edge.reachable_vars)) && rem_edge!(dg, edge)

        dg.dirty_roots .|= old_roots .!= edge.reachable_roots
        dg.dirty_vars .|= old_vars .!= edge.reachable_vars
    end

    # add new subgraph edge, using the dominance_mask for roots/vars as applicable
    sub_edge = Edge{T}(sub.subgraph_value, sub.top_vertex, sub.bott_vertex, sub_edge_reachable_vars(sub), sub_edge_reachable_roots(sub))
    add_edge!(dg, sub_edge)

    # propagate changes to reachability only after fully factoring subgraph
    propagate_root_reachability(dg, sub_edge)
    propagate_var_reachability(dg, sub_edge)

    for edge in sub_edges
        propagate_root_reachability(dg, edge)
        propagate_var_reachability(dg, edge)
    end
    return true
end

# ordering subgraphs should be factored in
struct FactorOrder <: Base.Order.Ordering
end

Base.Order.lt(::FactorOrder, a, b) = factor_order(a, b)
Base.isless(::FactorOrder, a, b) = factor_order(a, b)

function factor_order(a::FactorableSubgraph, b::FactorableSubgraph)
    a_diff = abs(a.top_vertex - a.bott_vertex)
    b_diff = abs(b.top_vertex - b.bott_vertex)

    # factor the smaller subgraph first (guarantees that if a ⊂ b then a is factored first)
    # then, factor the more used subgraph first

    return a_diff < b_diff || (a_diff == b_diff && subgraph_count(a) > subgraph_count(b))
end

# Factor all subgraphs in the `DerivativeGraph`. This is the key step in the D* algorithm.
# Each heap of candidate subgraphs is drained completely (cf. FastDifferentiation's
# single-pass `factor!`): `subgraph_exists` validates proposals lazily against the
# current edges, so stale proposals are skipped. A proposal's masks remain valid while
# stale because factoring preserves the graph's (root, var) pair set. The heap is only
# recomputed after a full drain that made progress, to pick up candidates created or
# resurrected by earlier factorings.
function factor_subgraphs!(dg::DerivativeGraph{T}) where {T}
    dom_cache = Dict{Int, Vector{Union{Nothing,T}}}()
    pdom_cache = Dict{Int, Vector{Union{Nothing,T}}}()
    subs = get_factorable_subgraphs(dg; dom_cache, pdom_cache)
    rejected = Set{Tuple{T,T,DataType,BitVector}}()
    factored = Set{Tuple{T,T,DataType,BitVector}}()

    while !isempty(subs)
        made_progress = false
        while !isempty(subs)
            # factor the first subgraph according to `FactorOrder`
            sub = pop!(subs)

            # an identical proposal was already factored or rejected: any remaining
            # parallel edges cover disjoint (root, var) pairs, so re-factoring
            # cannot make progress
            signature = (sub.top_vertex, sub.bott_vertex, typeof(sub), copy(sub.dominance_mask))
            (signature in factored || signature in rejected) && continue
            if factor_subgraph!(dg, sub)
                push!(factored, signature)
                made_progress = true
            else
                push!(rejected, signature)
            end
        end
        made_progress || break
        subs = get_factorable_subgraphs(dg; dom_cache, pdom_cache)
    end
end

# evaluate the derivative of root w.r.t. var using a fully factored DerivativeGraph
function evaluate_path(dg::DerivativeGraph{T}, root::Integer, var::Integer, cache::Vector{Dict{Edge{T},SymbolicT}}) where {T}
    haskey(dg.root_idx_to_postorder, root) || return COMMON_ZERO
    haskey(dg.var_idx_to_postorder, var) || return COMMON_ZERO
    root_postorder = dg.root_idx_to_postorder[root]
    var_postorder = dg.var_idx_to_postorder[var]

    root_postorder == var_postorder && return COMMON_ONE

    # sum the products of all paths from root to var (the factored graph is a
    # sum-of-products representation; parallel edges are distinct summands)
    terms = SymbolicT[]
    for e in dg.child_edges[root_postorder]
        (reachable_roots(e)[root] && reachable_vars(e)[var]) || continue
        push!(terms, evaluate_path(dg, e, root, var, cache))
    end
    return isempty(terms) ? COMMON_ZERO : SymbolicUtils.add_worker(VartypeT, terms)
end

function evaluate_path(dg::DerivativeGraph{T}, edge::Edge{T}, root::Integer, var::Integer, cache::Vector{Dict{Edge{T},SymbolicT}}) where {T}
    edge.bott_vertex == dg.var_idx_to_postorder[var] && return edge.edge_value # reached var
    haskey(cache[var], edge) && return cache[var][edge]

    terms = SymbolicT[]
    for e in dg.child_edges[edge.bott_vertex]
        (reachable_roots(e)[root] && reachable_vars(e)[var]) || continue
        push!(terms, evaluate_path(dg, e, root, var, cache))
    end
    result = (isempty(terms) ? COMMON_ZERO : SymbolicUtils.add_worker(VartypeT, terms)) * edge.edge_value
    cache[var][edge] = result

    return result
end

# if `ex` is a branch-guard partial `f(c, 1, 0)`/`f(c, 0, 1)` for an
# `ifelse`-family op `f`, returns `(f, c, is_then_branch)`; else `nothing`
function _ifelse_guard(ex)
    iscall(ex) || return nothing
    f = operation(ex)
    (f === ifelse || f === ifelse_eager || f === ifelse_branching) || return nothing
    args = SymbolicUtils.arguments(ex)
    _isone(args[2]) && _iszero(args[3]) && return (f, args[1], true)
    _iszero(args[2]) && _isone(args[3]) && return (f, args[1], false)
    return nothing
end

# if `ex`'s numerator product contains a branch-guard factor `f(c,1,0)` or
# `f(c,0,1)`, strips it and returns `(f, c, is_then_branch, rest)`; else
# `nothing`. Guards only ever appear as product factors — including `Div`
# numerators, since they multiply but never divide.
function _strip_ifelse_guard(ex)
    g = _ifelse_guard(ex)
    g === nothing || return (g..., COMMON_ONE)
    @match ex begin
        BSImpl.AddMul(; coeff, dict, variant) && if variant == SymbolicUtils.AddMulVariant.MUL end => begin
            for k in keys(dict)
                g = _ifelse_guard(k)
                g === nothing && continue
                rest = copy(dict)
                delete!(rest, k)
                return (g..., SymbolicUtils.Mul{VartypeT}(coeff, rest; type = symtype(ex), shape = shape(ex)))
            end
            return nothing
        end
        BSImpl.Div(; num, den, simplified) => begin
            r = _strip_ifelse_guard(num)
            r === nothing && return nothing
            f, c, is_then, stripped = r
            return (f, c, is_then, SymbolicUtils.Div{VartypeT}(stripped, den, simplified; type = symtype(ex), shape = shape(ex)))
        end
        _ => nothing
    end
end

# Branch-guard edge values multiply the whole downstream path product, so raw
# output contains `f(c,1,0)*X`-style factors. Folding them into
# `f(c,X,0)`/`f(c,0,X)` keeps the condition a genuine select — untaken branches
# are never evaluated and a dead-zone `0*Inf`/`0*NaN` cannot poison the result —
# matching `expand_derivatives`' `ifelse(c, Da, Db)` form.
function _fold_ifelse_guards(ex)
    r = _strip_ifelse_guard(ex)
    if r !== nothing
        f, c, is_then, rest = r
        folded = _fold_ifelse_guards(rest)
        return is_then ? f(c, folded, COMMON_ZERO) : f(c, COMMON_ZERO, folded)
    end
    @match ex begin
        BSImpl.AddMul(; coeff, dict, variant) => begin
            newdict = empty(dict)
            changed = false
            for (k, v) in dict
                k2 = _fold_ifelse_guards(k)
                changed |= k2 !== k
                newdict[k2] = v
            end
            changed || return ex
            if variant == SymbolicUtils.AddMulVariant.ADD
                # f(c,A,0)+f(c,0,B) selects exactly one of A,B — merge into
                # f(c,A,B) to match the `ifelse(c, Da, Db)` form
                for (k, v) in collect(newdict)
                    haskey(newdict, k) || continue
                    iscall(k) || continue
                    f = operation(k)
                    (f === ifelse || f === ifelse_eager || f === ifelse_branching) || continue
                    args = SymbolicUtils.arguments(k)
                    _iszero(args[2]) || continue
                    for (k2, v2) in collect(newdict)
                        (k2 === k || !haskey(newdict, k2) || !iscall(k2)) && continue
                        operation(k2) === f || continue
                        args2 = SymbolicUtils.arguments(k2)
                        (isequal(args2[1], args[1]) && _iszero(args2[3])) || continue
                        # v*f(c,0,B)+v2*f(c,A,0) = f(c,v2*A,v*B)
                        merged = f(args[1], v2 * args2[2], v * args[3])
                        delete!(newdict, k)
                        delete!(newdict, k2)
                        newdict[merged] = get(newdict, merged, 0) + 1
                        break
                    end
                end
                isempty(newdict) && return coeff
                return SymbolicUtils.Add{VartypeT}(coeff, newdict; type = symtype(ex), shape = shape(ex))
            end
            return SymbolicUtils.Mul{VartypeT}(coeff, newdict; type = symtype(ex), shape = shape(ex))
        end
        BSImpl.Div(; num, den, simplified) => begin
            n2 = _fold_ifelse_guards(num)
            d2 = _fold_ifelse_guards(den)
            (n2 === num && d2 === den) && return ex
            return SymbolicUtils.Div{VartypeT}(n2, d2, simplified; type = symtype(ex), shape = shape(ex))
        end
        _ => ex
    end
end

"""
$(SIGNATURES)

Computes the Jacobian of `roots` w.r.t. `vars` using the D* automatic differentiation algorithm using the [`DerivativeGraph`](@ref) data structure.

(see [this paper](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/main-65.pdf) for more details on the algorithm)

Mostly the same usage as [`jacobian`](@ref). More limited in input expressions (doesn't support nested differentials), but asymptotically faster for large Rn->Rm expressions.

# Arguments

- `roots::AbstractVector`: Vector of expressions to differentate or array-type symbolic expression (e.g. function registered with `@register_array_symbolic`)
- `vars::AbstractVector`: Vector of variables to differentate w.r.t. or single array-type variable (e.g. `@variables x[1:4]`)
"""
function dstar_jacobian(roots::AbstractVector, vars::AbstractVector{SymbolicT})
    roots isa Arr && (roots = scalarize(unwrap(roots)))
    roots isa AbstractVector{Num} && (roots = unwrap.(roots))

    # deduplicate roots; a duplicate root would populate the graph with an
    # identical subgraph and can cause redundant factoring proposals
    unique_roots = similar(roots, 0)
    root_map = Vector{Int}(undef, length(roots))
    seen = Dict{SymbolicT, Int}()
    for (root_idx, root) in enumerate(roots)
        unique_idx = get(seen, root, 0)
        if unique_idx == 0
            push!(unique_roots, root)
            unique_idx = length(unique_roots)
            seen[root] = unique_idx
        end
        root_map[root_idx] = unique_idx
    end

    # deduplicate vars the same way; a duplicate var maps several var indices to
    # the same node, and only one gets its reachability bit propagated
    unique_vars = similar(vars, 0)
    var_map = Vector{Int}(undef, length(vars))
    empty!(seen)
    for (var_idx, var) in enumerate(vars)
        unique_idx = get(seen, var, 0)
        if unique_idx == 0
            push!(unique_vars, var)
            unique_idx = length(unique_vars)
            seen[var] = unique_idx
        end
        var_map[var_idx] = unique_idx
    end

    dg = DerivativeGraph(unique_roots, unique_vars)
    factor_subgraphs!(dg)

    result = Matrix{SymbolicT}(undef, length(unique_roots), length(unique_vars))
    idx_type = keytype(dg.child_edges)
    cache = [Dict{Edge{idx_type},SymbolicT}() for _ in eachindex(unique_vars)]

    has_ifelse = any(dg.symbols) do s
        iscall(s) || return false
        o = operation(s)
        o === ifelse || o === ifelse_eager || o === ifelse_branching
    end

    for root in eachindex(unique_roots)
        for var in eachindex(unique_vars)
            r = evaluate_path(dg, root, var, cache)
            result[root, var] = has_ifelse ? _fold_ifelse_guards(r) : r
        end
    end

    return result[root_map, var_map]
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
    _res = dstar_jacobian(roots, vars)
    res = similar(_res, Num)
    map!(Num, res, _res)
    return res
end

"""
$(SIGNATURES)

Computes the derivative of `root` w.r.t. `var` using the D* differentiation algorithm.

Mostly the same usage as [`derivative`](@ref), but more limited in input expressions (doesn't support nested differentials).

Wrapper for R1->R1 case of `dstar_jacobian`. See [`dstar_jacobian`](@ref) for more information.

# Arguments
- `root`: Expression to differentate
- `var`: Variable to differentate w.r.t.
"""
dstar_derivative(root::Union{Num,SymbolicT}, var::Union{Num,SymbolicT}) = Num(only(dstar_jacobian(unwrap.([root]), unwrap.([var]))))
