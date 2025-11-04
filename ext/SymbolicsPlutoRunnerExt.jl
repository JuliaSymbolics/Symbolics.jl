module SymbolicsPlutoRunnerExt

import Symbolics

const WARNED_LATEXIFY = Ref(false)

function Symbolics.warn_load_latexify(::Nothing)
    Base.get_extension(Symbolics, :SymbolicsLatexifyExt) === nothing || return
    WARNED_LATEXIFY[] && return
    @warn """
    Attempting to print a symbolic expression in a Pluto notebook. Please run \
    `import Latexify` to enable pretty-printing of symbolic expressions.
    """
    WARNED_LATEXIFY[] = true
    return nothing
end

end
