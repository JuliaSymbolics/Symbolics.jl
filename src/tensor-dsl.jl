### @arrayop
#

struct ArrayOp
    op         # High-level operation
    output_idx # output indices
    expr       # Used in pattern matching
               # Useful to infer eltype
    arity
    reduce
end

positional(i) = Term{Array}(positional, [i])

macro arrayop(name, output_idx, expr, reduce=+)
    @assert output_idx.head == :tuple
    oidxs = filter(x->x isa Symbol, output_idx.args)
    iidxs, arity = find_indices(expr)
    idxs = union(oidxs, iidxs)
    fbody = call2term(deepcopy(expr))
    oftype(x,T) = :($x::$T)
    positionals = [:($(Symbol("_$i")) = $(Sym{Array}(Symbol("_$i"))))
                   for i in 1:arity]
    quote
        let
            @syms $(map(x->oftype(x, Int), idxs)...)
            $(positionals...)

            $ArrayOp($(esc(name)),
                     $output_idx,
                     $fbody,
                     $arity,
                     $reduce)
        end
    end
end

get_pos(x::Symbol) = startswith(string(x), "_") ? parse(Int, replace(string(x), "_"=>"")) : 0
# Find all symbolic indices in expr
function find_indices(expr, idxs=[], arity=Ref{Int}(0))
    !(expr isa Expr) && return idxs, arity[]
    if expr.head == :ref
        arity[] = max(arity[], get_pos(expr.args[1]))
        return append!(idxs, filter(x->x isa Symbol, expr.args[2:end])), arity[]
    elseif expr.head == :call && expr.args[1] == :getindex || expr.args[1] == getindex
        arity[] = max(arity[], get_pos(expr.args[2]))
        return append!(idxs, filter(x->x isa Symbol, expr.args[3:end])), max(arity, get_pos(expr.args[1])), arity[]
    else
        foreach(x->find_indices(x, idxs, arity), expr.args)
        return idxs, arity[]
    end
end

# replace _1 with __args__[1]
function call2term(expr, arrs=[])
    !(expr isa Expr) && return expr
    if expr.head == :call
        return Expr(:call, term, map(call2term, expr.args)...)
    end

    return Expr(expr.head, map(call2term, expr.args)...)
end


### Shape propagate

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

makeposvars(n) = [Sym{Array}(Symbol("_$i")) for i in 1:n]

function instantiate_rhs(aop::ArrayOp, args...)
    substitute(aop.expr, Dict(makeposvars(aop.arity) .=> args))
end

function propagate_shape(aop::ArrayOp, args...)
    output_idx = aop.output_idx
    expr = instantiate_rhs(aop, args...)

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

# TODO: have fallback
function propagate_eltype(aop::ArrayOp, args...)
    eltype(symtype(instantiate_rhs(aop, args...)))
end

function propagate_ndims(aop::ArrayOp, args...)
    length(aop.output_idx)
end
