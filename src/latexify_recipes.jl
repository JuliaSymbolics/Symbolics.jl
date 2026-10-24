"""
    _toexpr_metadata(O, ::Type{Ctx}, val; latexwrapper) -> Union{Nothing, Any}

Hook for customizing Latexify output based on metadata. Return `nothing` to fall back
to the default behavior.

Defined in `Symbolics` (not in the Latexify extension) so downstream packages can add
methods via `import Symbolics` without reaching into the extension with
`Base.get_extension`. The Latexify extension calls this hook while building LaTeX
expressions and provides the `_toexpr_plain` and `default_latex_wrapper` helpers for
composing output inside a method.
"""
_toexpr_metadata(::Any, ::DataType, @nospecialize(val); kwargs...) = nothing

"""
    _toexpr_op(op, args; latexwrapper) -> Union{Nothing, Any}

Hook for customizing Latexify output based on the operation of a call. Return `nothing`
to fall back to the default behavior.

Defined in `Symbolics` (not in the Latexify extension) so downstream packages can add
methods via `import Symbolics` without reaching into the extension with
`Base.get_extension`.
"""
_toexpr_op(@nospecialize(op), args; kwargs...) = nothing
