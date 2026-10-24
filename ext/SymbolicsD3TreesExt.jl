module SymbolicsD3TreesExt

using Symbolics
using SymbolicUtils
using D3Trees
using D3Trees.AbstractTrees

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

    merged_kwargs = merge(
        (
            detect_repeat = false,
            text = sym_nodelabel, tooltip = string,
            lazy_expand_after_depth = 8, lazy_subtree_depth = 4, init_expand = 8
        ),
        kwargs
    )

    D3Tree(symbolic; merged_kwargs...)
end

end # module SymbolicsD3TreesExt