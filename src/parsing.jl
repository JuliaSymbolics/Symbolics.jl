"""
```julia
parse_expr_to_symbolic(ex, mod::Module)
```

Applies the `parse_expr_to_symbolic` function in the current module, i.e.
`parse_expr_to_symbolic(ex, mod)` where `mod` is the module of the function
caller.

## Arguments

* `ex`: the expression to parse
* `mod`: the module to apply the parsing in. See the limitations section for details.

## Example

```julia
ex = :(y(t) ~ x(t))
parse_expr_to_symbolic(ex,Main) # gives the symbolic expression `y(t) ~ x(t)` in empty Main

# Now do a whole system

ex = [:(y ~ x)
      :(y ~ -2x + 3 / z)
      :(z ~ 2)]
eqs = parse_expr_to_symbolic.(ex, (Main,))

@variables x y z
ex = [y ~ x
      y ~ -2x + 3 / z
      z ~ 2]
all(isequal.(eqs,ex)) # true
```
## Limitations

### Symbolic-ness Tied to Environment Definitions

The parsing to a symbolic expression has to be able to recognize the difference between
functions, numbers, and globals defined within one's Julia environment and those that
are to be made symbolic. The way this functionality handles this problem is that it
does not define anything as symbolic that is already defined in the chosen `mod` module.
For example, `f(x,y)` will have `f` as non-symbolic if the function `f` (named `f`)
is defined in `mod`, i.e. if `isdefined(mod,:f)` is true. When the symbol is defined, it
will be replaced by its value. Notably, this means that the parsing behavior changes
depending on the environment that it is applied.

For example:

```julia
parse_expr_to_symbolic(:(x - y),@__MODULE__) # x - y
x = 2.0
parse_expr_to_symbolic(:(x - y),@__MODULE__) # 2.0 - y
```

This is required to detect that standard functions like `-` are functions instead of
symbolic symbols. For safety, one should create anonymous modules or other sub-environments
to ensure no stray variables are defined.

### Metadata is Blank

Because all the variables defined by the expressions are not defined with the standard
`@variables`, there is no metadata that is or can be associated with any of the generated
variables. Instead, they all have blank metadata, but are defined in the `Real` domain.
Thus, the variables which come out of this parsing may not evaluate as equal to a symbolic
variable defined elsewhere.
"""
function parse_expr_to_symbolic end

parse_expr_to_symbolic(x::Number, mod::Module) = x
function parse_expr_to_symbolic(x::Symbol, mod::Module)
  if isdefined(mod, x)
    getfield(mod, x)
  else
    (@variables $x)[1]
  end
end
function parse_expr_to_symbolic(ex, mod::Module)
    if ex.head == :call
        if isdefined(mod, ex.args[1])
            return getfield(mod,ex.args[1])(parse_expr_to_symbolic.(ex.args[2:end],(mod,))...)
        else
            x = parse_expr_to_symbolic(ex.args[1], mod)
            ys = parse_expr_to_symbolic.(ex.args[2:end],(mod,))
            return Term{Real}(x,[ys...])
        end
    elseif ex.head == :ref
        return ex
    end
end

"""
```julia
@parse_expr_to_symbolic ex
```

Applies the `parse_expr_to_symbolic` function in the current module, i.e.
`parse_expr_to_symbolic(ex, mod)` where `mod` is the module of the function
caller.

## Arguments

* `ex`: the expression to parse

## Example

```julia
ex = :(y(t) ~ x(t))
@parse_expr_to_symbolic ex # gives the symbolic expression `y(t) ~ x(t)`
```

## Limitations

The same limitations apply as for the function `parse_expr_to_symbolic`.
See its docstring for more details.
"""
macro parse_expr_to_symbolic(ex)
    :(parse_expr_to_symbolic($ex, @__MODULE__))
end
