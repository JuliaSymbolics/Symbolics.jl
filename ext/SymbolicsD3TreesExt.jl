module SymbolicsD3TreesExt

using Symbolics
using SymbolicUtils
using D3Trees
using AbstractTrees

function sym_nodelabel(x)
    if x isa SymbolicUtils.BasicSymbolic
        if SymbolicUtils.iscall(x)
            string(SymbolicUtils.operation(x))
        else
            string(x)
        end
    else
        string(x)
    end
end

function D3Trees.D3Tree(t::Symbolics.Num; kwargs...)
    symbolic = Symbolics.unwrap(t)
    nodes = collect(AbstractTrees.PreOrderDFS(symbolic))
    merged_kwargs = merge(
        (
            text = sym_nodelabel.(nodes),
            tooltip = string.(nodes),
            init_expand = 3
        ),
        kwargs
    )
    D3Tree(symbolic; merged_kwargs...)
end

end # module SymbolicsD3TreesExt