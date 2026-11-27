using Symbolics
using Symbolics: InverseFourierTransform, unwrap
using Symbolics.SymbolicUtils: iscall, operation, arguments, isconst, unwrap_const
using Test

@variables t w a d
@syms f(t)

# `h` stands for `f` in the numerical checks: smooth, decaying, neither even nor odd.
h(u) = exp(-u^2) * (1 + u / 2)

# Simpson's rule on [-L, L]; the integrands below decay like exp(-u^2/9).
function quad(g; L = 60.0, n = 24_000)
    du = 2L / n
    return du / 3 * sum((i == 0 || i == n ? 1 : iseven(i) ? 2 : 4) * g(-L + i * du) for i in 0:n)
end

# `(N, β)` with the transform N ∫ g(u) exp(iβku) du evaluated at k.
function kernel(::Type{T}, convention) where {T}
    N, β = Dict(:angular => (1.0, -1.0), :ordinary => (1.0, -2π), :unitary => (1 / sqrt(2π), -1.0))[convention]
    T === FourierTransform && return N, β
    return (convention === :angular ? 1 / 2π : N), -β
end
transform_num(T, convention, g, k) = ((N, β) = kernel(T, convention); N * quad(u -> g(u) * cis(β * k * u)))

# Value of `ex` with the symbols in `vals` and `f` replaced by `h`, computing every transform
# term by quadrature.
function evalnum(ex, vals)
    ex = unwrap(ex)
    isconst(ex) && return unwrap_const(ex)
    iscall(ex) || return vals[ex]
    op, args = operation(ex), arguments(ex)
    if (op === real || op === imag) && Symbolics.is_fourier_term(args[1])
        F = operation(args[1])
        @assert isequal(only(arguments(args[1])), unwrap(f(F.var)))
        return op(transform_num(typeof(F), F.convention, h, evalnum(F.point, vals)))
    end
    return op(map(x -> evalnum(x, vals), args)...)
end
evalnum(z::Complex{Num}, vals) = evalnum(real(z), vals) + im * evalnum(imag(z), vals)

@testset "rules agree with quadrature of the definition" begin
    # forward in `t`, evaluated at `w`; inverse in `w`, evaluated at `t`
    for (T, x, k) in ((FourierTransform, t, w), (InverseFourierTransform, w, t)),
            convention in (:angular, :ordinary, :unitary)
        transform = T === FourierTransform ? fourier : inverse_fourier
        cases = [
            (f(x), u -> h(u)),
            (3f(x) - 2a * f(x) / 5, u -> (3 - 2 * 0.7 / 5) * h(u)),
            (f(2x + 1), u -> h(2u + 1)),
            (f(-x), u -> h(-u)),
            (f(x - d), u -> h(u - 0.3)),
            (f(-x / 3 + d) + f(0.5x) / a, u -> h(-u / 3 + 0.3) + h(0.5u) / 0.7),
            (Complex{Num}(Num(f(x)), Num(f(x + 1))), u -> h(u) + im * h(u + 1)),
        ]
        for (ex, g) in cases
            result = transform(ex, x => k; convention)
            @test result isa Complex{Num}
            for kv in (-1.3, 0.0, 0.4, 2.1)
                vals = Dict(unwrap(a) => 0.7, unwrap(d) => 0.3, unwrap(k) => kv)
                @test evalnum(result, vals) ≈ transform_num(T, convention, g, kv) atol = 1.0e-10
            end
        end
    end
end

@testset "constants" begin
    F = FourierTransform(t => w)
    @test isequal(fourier(0, t => w), zero(Complex{Num}))
    @test isequal(fourier(Num(0), t => w), zero(Complex{Num}))
    @test isequal(fourier(a, t => w), a * F(Num(1)))
    @test isequal(fourier([f(t), 0], t => w), [F(f(t)), zero(Complex{Num})])
end

@testset "everything else stays an exact unevaluated term" begin
    F = FourierTransform(t => w)
    @syms s z
    I = Integral(s in Symbolics.DomainSets.ClosedInterval(0, t^2))(s)
    for ex in (
            f(t)^2, f(t) * f(2t), exp(-t^2), t * f(t), f(t^2), f(a * t), f(t + sign(t)),
            f(t) * I, f(t) / (1 + t^2), f(t + z), t, f(Inf * t), f(1.0e-320 * t),
        )
        @test isequal(fourier(ex, t => w), F(ex))
    end
    @test isequal(fourier(f(t)^2 + f(t - 1), t => w), F(f(t)^2) + fourier(f(t - 1), t => w))
    # the integrand is expanded first, and an operator term in it is not taken as constant
    inner = real(FourierTransform(d => t)(f(2d)))
    @test Symbolics.hasnode(unwrap(real(fourier(inner * f(t), t => w)))) do x
        Symbolics.is_fourier_term(x) && isequal(operation(x).point, unwrap(t * (1 // 2)))
    end
end

@testset "numeric overflow" begin
    # exact: integer arithmetic does not wrap around
    m = typemin(Int)
    @test isequal(fourier(f(m * t), t => w), abs(1 // big(m)) * FourierTransform(t => w * (1 // big(m)))(f(t)))
    @test isequal(fourier(f(2^62 * (2t + a)), t => w), fourier(f(big(2)^63 * t + big(2)^62 * a), t => w))
    # floats are computed without underflow or overflow, and a result a Float64 cannot hold
    # leaves the term unevaluated
    small = fourier(f(1.0e-170t + 1.0e-170), t => 1.0e-170)
    @test Symbolics.hasnode(z -> iscall(z) && operation(z) === sin, unwrap(imag(small)))
    for (ex, k) in ((f(t + 1.0e300), 1.0e300), (f(t + 1.0e-200), 1.0e-200), (f(1.0e200t + 1.0e-200a), w))
        @test isequal(fourier(ex, t => k), FourierTransform(t => k)(ex))
    end
    @test Symbolics._hasbigfloat(real(fourier(f(t + big(1) / 3), t => w)))
end

@testset "operator construction" begin
    @test isequal(FourierTransform(t => w), FourierTransform(t => w))
    @test hash(FourierTransform(t => w)) == hash(FourierTransform(t => w))
    @test FourierTransform(t => w) != InverseFourierTransform(t => w)
    @test !isequal(FourierTransform(t => w), FourierTransform(t => w; convention = :unitary))
    @test isequal(FourierTransform(t => w)(f(t)), FourierTransform(t => w)(f(t)))
    @test hash(FourierTransform(t => w)(f(t))) == hash(FourierTransform(t => w)(f(t)))
    @test repr(FourierTransform(t => w)) == "FourierTransform(t => w)"
    @test repr(InverseFourierTransform(w => t; convention = :ordinary)) == "InverseFourierTransform(w => t; convention = :ordinary)"
    @syms z
    for p in (t => t, t => 2t + w, t => im, t => z, 2t => w, t => real(FourierTransform(d => w)(f(d))))
        @test_throws ArgumentError FourierTransform(p)
    end
    @test_throws ArgumentError FourierTransform(t => w; convention = (1, -1))
    @test_throws ArgumentError FourierTransform(t => w)(1)
end

@testset "safety of unevaluated terms" begin
    T = FourierTransform(t => w)(f(t))
    Ta = FourierTransform(t => w)(f(a * t))
    # derivatives stay unevaluated, never 0 or a chain-rule product
    for (D, ex) in ((Differential(w), real(T)), (Differential(w), imag(T)), (Differential(a), real(Ta)), (Differential(t), real(T)))
        @test isequal(expand_derivatives(D(ex)), D(ex))
    end
    bare = unwrap(T)
    @test isequal(expand_derivatives(Differential(w)(bare)), Differential(w)(bare))

    # the sparsity pattern includes the evaluation point and the parameters
    @test Symbolics.jacobian_sparsity([real(T)], [t, w])[1, 2]
    @test all(Symbolics.jacobian_sparsity([imag(Ta)], [a, w]))

    # substitution replaces a whole term, and does not enter one
    @test isequal(substitute(real(T) + w, Dict(w => 1)), real(T) + 1)
    @test isequal(substitute_in_deriv(real(T), Dict(w => 1, t => 2)), real(T))
    @test isequal(substitute(real(T) + 1, Dict(real(T) => 2)), Num(3))
    G = Integral(t in Symbolics.DomainSets.ClosedInterval(0, a))(real(T))
    @test Symbolics.hasnode(isequal(unwrap(real(T))), unwrap(expand_derivatives(Differential(a)(G))))
    # the Leibniz rule would substitute the integration variable into a transform term
    @syms s
    for G in (
            Integral(w in Symbolics.DomainSets.ClosedInterval(0, a))(real(T)),
            Integral(s in Symbolics.DomainSets.ClosedInterval(0, a))(real(FourierTransform(t => w)(s * f(t)))),
        )
        @test isequal(expand_derivatives(Differential(a)(G)), Differential(a)(G))
    end
    # transform terms are opaque to the Hessian sparsity analysis
    @variables x y
    H = Symbolics.hessian_sparsity(real(FourierTransform(t => w)(a * f(t))) * x + y^2, [w, x, y, a, d])
    # rows and columns of `w` and `a` (in the transform term) are marked against every variable
    # that occurs; `d` does not occur, and `x`, `y` pair only as the rest of the expression says
    @test Matrix(H) == Bool[1 1 1 1 0; 1 0 0 1 0; 1 0 1 1 0; 1 1 1 1 0; 0 0 0 0 0]
    @test !Symbolics.isaffine(real(T), [w])

    # generated code cannot evaluate a transform term
    @test_throws UndefVarError build_function(real(T), w; expression = Val(false))(1.0)
end
