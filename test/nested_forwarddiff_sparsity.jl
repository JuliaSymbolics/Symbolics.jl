using ForwardDiff, SparseArrays, Symbolics, PreallocationTools, Test
# Test Nesting https://discourse.julialang.org/t/preallocationtools-jl-with-nested-forwarddiff-and-sparsity-pattern-detection-errors/107897

function foo(x, cache)
    d = get_tmp(cache, x)

    d[:] = x

    0.5 * x'*x
end

function residual(r, x, cache)
    function foo_wrap(x)
        foo(x, cache)
    end

    r[:] = ForwardDiff.gradient(foo_wrap, x)
end

cache = DiffCache(zeros(2))
pattern = Symbolics.jacobian_sparsity((r, x) -> residual(r, x, cache), zeros(2), zeros(2))
@test pattern == sparse([1 0
                         0 1])
