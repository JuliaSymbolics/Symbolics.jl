# Frequently Asked Questions

## Limits of Symbolic Computation

### Transforming my function to a symbolic equation has failed. What do I do?

If you see the error:

```
ERROR: TypeError: non-boolean (Num) used in boolean context
```

this is likely coming from an algorithm which cannot be traced into a purely
symbolic algorithm. Many numerical solvers for instance have this property. It
shows up when you're doing something like if `x < tol`. If x is a number, then
this is true or false. If x is symbol, then it's `x < tol`, so Julia just cannot
know how many iterations to do and throws an error.

This shows up in adaptive algorithms, for example:

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

The number of iterations this algorithm runs for is dependent on the value of
`x`, and so there is no static representation of the algorithm. If `x` is 5,
then it's `out = x*(x-1)*(x-2)*(x-3)*(x-4)`, while if `x` is 3 then it's
`out = x*(x-1)*(x-2)`. Thus it should be no surprise that:

```julia
@variables x
factorial(x)
```

fails. It's not that there is anything wrong with this code, but it's not going
to work because fundamentally this is not a symbolically-representable algorithm.

The space of algorithms which can be turned into symbolic algorithms is what we
call quasi-static, that is, there is a way to represent the algorithm as static.
Loops are allowed, but the amount of loop iterations should not require that you
know the value of the symbol `x`. If the algorithm is quasi-static, then Symbolics.jl
tracing will produce the static form of the code, unrolling the operations, and
generating a flat representation of the algorithm.

#### What can be done?

If you need to represent this function `f` symbolically, then you'll need to make
sure it's not traced and instead is directly represented in the underlying
computational graph. Just like how `sqrt(x)` symbolically does not try to
represent the underlying algorithm, this needs to be done to your `f`. This is
done by doing `@register_symbolic f(x)`. If you need to define things like derivatives to
`f`, then [the function registration documentation](@ref function_registration).
