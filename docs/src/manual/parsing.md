# Parsing Julia Expressions to Symbolic Expressions

Julia expressions such as `:(y - x)` are fundamentally different from symbolic
expressions, as they do not have an algebra defined on them. Thus, it can be
very helpful when building domain-specific languages (DSLs) and parsing files
to convert from Julia expressions to Symbolics.jl expressions for further
manipulation. Towards this end is the `parse_expr_to_symbolic` which performs
the parsing.

!!! warn
    Take the limitations mentioned in the `parse_expr_to_symbolic` docstrings
    seriously! Because Julia expressions contain no symbolic metadata, there
    is limited information and thus the parsing requires heuristics to work. 

```@docs
parse_expr_to_symbolic
@parse_expr_to_symbolic
```