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

const SYMBOLIC_BOOLEAN_CONTROL_FLOW = """
    Symbolic expression found in a control flow expression that requires a boolean. This can occur
    for example in a function like:

    ```julia
    f(x) = if x > 1
        x^2
    else
        2x
    end
    ```

    Notice that if `x` is symbolic, then `x > 1` is a symbolic expression and `if x > 1`
    does not know which branch to take without knowing the value of `x`, and thus we
    cannot know the mathematical expression. Thus instead of using `if/else`, use the
    function `ifelse(x>1, x^2, 2x)`.

    Sometimes there are functions that are fundamentally incompatible with symbolic tracing.
    For example:

    ```julia 
    function factorial(x)
        out = x
        while x > 1
            x -= 1
            out *= x
        end
        out
    end
    ```

    In this function, if `x=5` then the expression is ` x*(x-1)*(x-2)*(x-3)*(x-4)`,
    while if `x=3` then ` x*(x-1)*(x-2)`. Thus the symbolic expression is ill-defined,
    and Symbolics is not compatible with this function. If you have this case, you
    may want to register the symbolic function as a primitive, i.e.

    ```julia
    @register_symbolic f(x)
    ```

    then makes `f(x)` work by not expanding this function.

    For more information, see https://docs.sciml.ai/Symbolics/stable/manual/faq/ and
    https://docs.sciml.ai/Symbolics/stable/manual/functions/ .
"""

const SYMBOLIC_ARRAY_BOUNDS_ERROR = """
    A symbolic variable was used in array bounds or range construction. This commonly occurs when:
    
    1. Trying to create arrays with symbolic dimensions:
    
    ```julia
    @variables n
    @variables X[1:n]  # Error: n is symbolic, array size must be concrete
    ```
    
    To fix this, use concrete (numerical) values for array bounds:
    
    ```julia
    @variables X[1:10]  # Works: concrete size
    ```
    
    If you need variable-sized arrays, consider using a vector of variables or working 
    with symbolic indexing at runtime rather than at variable declaration time.
    
    For more information about array variables, see https://docs.sciml.ai/Symbolics/stable/manual/variables/
"""

Base.Experimental.register_error_hint(TypeError) do io, e, args, kwargs
    if e.expected === Bool && typeof(e.got) <: Union{Symbolics.Arr, Num, Symbolics.BasicSymbolic}
        println(io, "\n")
        println(io, SYMBOLIC_BOOLEAN_CONTROL_FLOW)
        # Also show array bounds hint as this is a common case
        println(io, "\n") 
        println(io, SYMBOLIC_ARRAY_BOUNDS_ERROR)
    end
end
