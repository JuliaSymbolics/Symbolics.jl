const _ = Sym{AbstractArray}(:_)

struct TensorOp
    output_idx::Tuple
    expr
end

struct AxisOf
    A
    dim
end

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

function shape_propagate(t::TensorOp)
    matches = idx_to_axes(t.expr)
    for (sym, ms) in matches
        @assert !isempty(ms) "dimension of $sym is unknown"
        to_check = filter(m->!isnothing(shape(m.A)), ms)
        # Only check known dimensions. It may be "known symbolically"
        isempty(to_check) && continue
        reference = axes(first(to_check).A, first(to_check).dim)
        for i in 2:length(ms)
            m = ms[i]
            @maybe s=shape(m.A) begin
                @assert isequal(axes(m.A, m.dim), reference)
            end
        end
    end

    map(t.output_idx) do i
        mi = matches[i]
        @assert !isempty(mi)
        first(mi)
    end
end
