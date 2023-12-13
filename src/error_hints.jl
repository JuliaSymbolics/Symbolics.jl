const FLOAT_NUM_CONVERSION_ERROR = """

An implicit conversion of symbolic Num into a Float64 was encountered. Common causes for this error include:

1. A substitution resulted in numerical values, but is not automatically converted to numerical values. For example:

```julia
@variables x y
result = [1.0]
A = [x x^2
     y+x 2x + y]
sub = substitute(A, Dict([x=>2.0, y=>2]))
result[1] = sub[1]
```

The reason is because the result of `sub` is still symbolic:

```julia
julia> sub = substitute(A, Dict([x=>2.0, y=>2]))
2×2 Matrix{Num}:
 2.0  4.0
 4.0  6.0
```

Notice that the return is a `Matrix{Num}`, not a `Matrix{Float64}`. To fix this, ensure that the type is converted, i.e.
`sub = Symbolics.unwrap.(substitute(A, Dict([x=>2.0, y=>2])))`, or `result[1] = Symbolics.unwrap(sub[1])`.

2. Captured symbolic values inside of registered functions. An example of this version of the error is the following:

```julia
@variables x y
ff(x) = x + y
@register_symbolic ff(x)
result = [0.0]
result[1] = eval(build_function(ff(x),[x,y]))(1.0) # Error
```

Notice that the is the fact that the generated function call `eval(build_function(ff(x),x))(1.0)` returns a symbolic
value (`1.0 + y`), not a numerical value. This is because the value `y` is a global symbol enclosed into `ff`. To fix
this, ensure that any symbolic value that is used in registered functions is passed in, i.e. :

```julia
ff(x) = x + y
@register_symbolic ff(x)
```
"""

Base.Experimental.register_error_hint(MethodError) do io, e, args, kwargs
  if e isa MethodError && Num in args && e.f <: Number
      println(io, FLOAT_NUM_CONVERSION_ERROR)
  end
end