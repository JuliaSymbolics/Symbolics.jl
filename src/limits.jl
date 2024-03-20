"""
    limit(expr, var, h[, side::Symbol])

Compute the limit of `expr` as `var` approaches `h`.

`side` indicates the direction from which `var` approaches `h`. It may be one of `:left`,
`:right`, or `:both`. If `side` is `:both` and the two sides do not align, an error is
thrown. Side defaults to `:both` for finite `h`, `:left` for `h = Inf`, and `:right` for
`h = -Inf`.

`expr` must be compoesed of `log`, `exp`, constants, and the rational opperators `+`, `-`,
`*`, and `/`. This limitation may eventually be relaxed.

!!! warning
    Because symbolic limit computation is undecidable, this function necessarily employs
    heuristics and may occasionally return wrong answers. Nevertheless, please report wrong
    answers as issues as we aim to have heuristics that produce correct answers in all
    practical cases.
"""
limit(expr, var, h, side...) = SymbolicLimits.limit(expr, var, h, side...)[1]
