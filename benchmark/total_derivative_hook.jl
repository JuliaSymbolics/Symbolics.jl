using Symbolics, BenchmarkTools, Printf

@variables x t u(t) p[1:400] z[1:30]
Dx = Differential(x)
Dt = Differential(t)
simple = Dx(sin(x)^2 + x^3 + exp(x))
time_expr = Dt(sin(u)^2 + u^3 + exp(u))
params(n) = Dx(x * sum(sin(p[i]) + exp(p[i]^2) for i in 1:n))
expr40, expr400 = params(40), params(400)
zs = collect(z)
eqs = [sin(z[i]) + z[i + 1]^2 + exp(z[i + 2]) + p[i] * z[i] for i in 1:28]

cases = [
    "simple Dx" => @benchmarkable(expand_derivatives($simple)),
    "Dt with u(t)" => @benchmarkable(expand_derivatives($time_expr)),
    "28x30 jacobian" => @benchmarkable(Symbolics.jacobian($eqs, $zs)),
    "Dx, 40 parameters" => @benchmarkable(expand_derivatives($expr40)),
    "Dx, 400 parameters" => @benchmarkable(expand_derivatives($expr400))
]
println("Julia ", VERSION, "; source ", pathof(Symbolics))
println("Warm derivative caches; median of up to 10000 samples, evals=10 (3 seconds/case max).")
for (label, bench) in cases
    run(bench; samples = 10, evals = 1) # compile and populate all caches
    trial = run(bench; samples = 10000, evals = 10, seconds = 3)
    estimate = median(trial)
    @printf("%s: %.3f us; %d bytes; %d allocations\n",
        label, estimate.time / 1000, estimate.memory, estimate.allocs)
end
