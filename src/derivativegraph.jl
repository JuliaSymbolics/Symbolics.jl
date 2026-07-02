struct Edge{T<:Integer}
    value
    top_vertex::T
    bott_vertex::T
    reachable_vars::BitVector
    reachable_roots::BitVector
end

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
        Dict{idx_type, Vector{Edge{idx_type}}}()
    )

    populate_dergraph!(dg::DerivativeGraph{idx_type})

    return dg
end

Base.eachindex(dg::DerivativeGraph{T}) where {T} = T(1):T(length(dg.symbols)) # iterator over postorder indices with type T


function reachable_roots(dg::DerivativeGraph{T}, node::T) where {T}
    edges = dg.parent_edges[node]
    roots_mask = falses(length(dg.roots))
    for edge in edges
        roots_mask .|= edge.reachable_roots
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
        vars_mask .|= edge.reachable_vars
    end

    if haskey(dg.postorder_to_var_idx, node)
        vars_mask[dg.postorder_to_var_idx[node]] = 1
    end

    return vars_mask
end

function populate_dergraph!(dg::DerivativeGraph)
    for (root_idx, root) in enumerate(dg.roots)
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
    haskey(dg.definitions, expr) && populate_root_reachabilities(dg, dg.definitions[expr], root_idx)

    # assume iscall(expr)
    args = parent(arguments(expr))
    arg_idx_to_post_idx = IdDict{T, T}()
    sizehint!(arg_idx_to_post_idx, length(args)) # faster?
    for (arg_idx, arg) in enumerate(args)
        if arg in dg.varset
            arg_idx_to_post_idx[T(arg_idx)] = populate_dergraph_var!(dg, arg, root_idx)
        elseif iscall(arg)
            arg_idx_to_post_idx[T(arg_idx)] = populate_dergraph!(dg, arg, root_idx)
        end
    end

    push!(dg.symbols, expr)
    post_idx::T = length(dg.symbols)
    dg.definitions[expr] = post_idx
    dg.child_edges[post_idx] = Edge{T}[]
    dg.parent_edges[post_idx] = Edge{T}[]

    for (arg_idx, arg_post_idx) in arg_idx_to_post_idx
        partial_der = derivative_idx(expr, arg_idx)
        arg_reachable_vars = reachable_vars(dg, arg_post_idx)
        arg_reachable_roots = reachable_roots(dg, arg_post_idx)
        arg_reachable_roots[root_idx] = 1
        new_edge = Edge{T}(partial_der, post_idx, arg_post_idx, arg_reachable_vars, arg_reachable_roots)
        push!(dg.child_edges[post_idx], new_edge)
        push!(dg.parent_edges[arg_post_idx], new_edge)
    end

    return post_idx
end

function populate_root_reachabilities!(dg::DerivativeGraph, node::Integer, root_idx::Integer)
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