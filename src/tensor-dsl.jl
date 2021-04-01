### @arrayop
#

struct ArrayOp
    op         # High-level operation
    output_idx # output indices
    expr       # Used in pattern matching
    f          # When called with Symbolic arguments, gives the RHS of a tensor expression as a Term
               # Useful to infer eltype
end

macro arrayop(name, output_idx, expr, options...)
    @assert output_idx.head == :tuple
    oidxs = filter(x->x isa Symbol, output_idx)
    idxs = union(oidxs, find_indices(expr))
    fbody = position_args(deepcopy(expr))
    oftype(x,T) = :($x::$T)
    quote
        let
            @syms $(map(x->oftype(x, Int), idxs)...)

            $ArrayOp($(esc(name)),
                     $output_idx,
                     $expr,
                     (__args__...,) -> $fbody)
        end
    end
end

# Find all symbolic indices in expr
function find_indices(expr, idxs=[])
    !(expr isa Expr) && return idxs
    if expr.head == :ref
        return append!(idxs, expr.args[2:end])
    elseif expr.head == :call && expr.args[1] == :getindex || expr.args[1] == getindex
        return append!(idxs, filter(x->x isa Symbol, expr.args[3:end]))
    else
        foreach(x->find_indices(x, idxs), expr.args)
        return idxs
    end
end

struct AxisOf
    A
    dim
end

Base.get(a::AxisOf) = axes(a.A, a.dim)

function idx_to_axes(expr, dict=Dict{Sym, Vector}())
    if istree(expr)
        if operation(expr) === (getindex)
            args = arguments(expr)
            for (axis, sym) in enumerate(@views args[2:end])
                !(sym isa Sym) && continue
                axesvec = Base.get!(() -> [], dict, sym)
                push!(axesvec, AxisOf(car(args), axis))
            end
        else
            foreach(ex->idx_to_axes(ex, dict), arguments(expr))
        end
    end
    dict
end


function shape_propagate(output_idx, expr)
    matches = idx_to_axes(expr)
    for (sym, ms) in matches
        @assert !isempty(ms) "dimension of $sym is unknown"
        to_check = filter(m->!isnothing(shape(m.A)), ms)
        # Only check known dimensions. It may be "known symbolically"
        isempty(to_check) && continue
        reference = axes(first(to_check).A, first(to_check).dim)
        for i in 2:length(ms)
            m = ms[i]
            s=shape(m.A)
            if s !== Unknown()
                if !isequal(axes(m.A, m.dim), reference)
                    "expected axes($(m.A), $(m.dim)) = $(reference)" |> DimensionMismatch |> throw
                end
            end
        end
    end

    map(output_idx) do i
        mi = matches[i]
        @assert !isempty(mi)
        get(first(mi))
    end
end

# replace _1 with __args__[1]
function position_args(expr, arrs=[])
    !(expr isa Expr) && return expr
    if expr.head == :ref
        name = expr.args[1]
        if startswith(string(name), "_")
            i = parse(Int, replace(string(name), "_" => ""))
            return Expr(:ref,
                        :(__args__[$i]),
                        expr.args[2:end]...)
        end
    elseif expr.head == :call &&
        expr.args[1] == :getindex ||
        expr.args[1] == getindex

        name = expr.args[2]
        if startswith(string(name), "_")
            i = parse(Int, replace(string(name), "_" => ""))
            return Expr(:call,
                        expr.args[1],
                        :(__args__[$i]),
                        expr.args[3:end]...)
        end
    end

    return Expr(expr.head, map(position_args, expr.args)...)
end
