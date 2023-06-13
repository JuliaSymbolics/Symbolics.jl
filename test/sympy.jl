if get(ENV, "CI", nothing) !== nothing
    # copied and modified from https://github.com/tkf/IPython.jl/blob/master/test/install_dependencies.jl


    # Adding Pkg in test/REQUIRE would be an error in 0.6.  Using
    # Project.toml still has some gotchas.  So:
    Pkg = Base.require(Base.PkgId(Base.UUID(0x44cfe95a1eb252eab672e2afdf69b78f), "Pkg"))

    # Let PyCall.jl use Python interpreter from Conda.jl
    # See: https://github.com/JuliaPy/PyCall.jl
    ENV["PYTHON"] = ""
    Pkg.build("PyCall")
end

using Test
using SymPy
using Symbolics

@variables t p(t) x y(t)
expr = x * p + (x^2 - 1 + y) * (p + 2t)
sexpr = symbolics_to_sympy(expr)
sp = symbolics_to_sympy(p)

@test SymPy.simplify(symbolics_to_sympy(Symbolics.solve_for(expr, p))) == SymPy.solve(sexpr, sp)[1]

@variables k
sk = symbols("k", positive = True)
@test symbolics_to_sympy(k) !== sk
@test symbolics_to_sympy(k, Dict()) != sk
@test symbolics_to_sympy(k, Dict("k"=>sk)) == sk
@test symbolics_to_sympy(k, Dict("x"=>sk)) != sk
@test symbolics_to_sympy(k, [sk]) == sk
@test symbolics_to_sympy(k, (sk,)) == sk
