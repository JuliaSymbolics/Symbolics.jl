using Documenter, Symbolics, SymbolicUtils

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

# Make sure that plots don't throw a bunch of warnings / errors!
ENV["GKSwstype"] = "100"
ENV["JULIA_DEBUG"] = "Documenter"
using Plots

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(
    sitename="Symbolics.jl",
    authors="Chris Rackauckas",
    modules=[Symbolics,SymbolicUtils],
    clean=true, doctest=false, linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :cross_references, :linkcheck],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
                             mathengine = mathengine,
                             canonical="https://docs.sciml.ai/Symbolics/stable/"),
    pages=[
        "Home" => "index.md",
        "getting_started.md",
        "Tutorials" => Any[
            "tutorials/auto_parallel.md",
            "tutorials/converting_to_C.md",
            "tutorials/perturbation.md"
        ],
        "Manual" => Any[
            "Basics" => [
                "manual/variables.md",
                "manual/expression_manipulation.md",
                "manual/faq.md",
                "manual/external.md",
                "manual/arrays.md",
                "manual/build_function.md",
                "manual/functions.md",
                "manual/parsing.md",
                "manual/io.md",
                "manual/types.md",
            ],

            "Algebra" => [
                "manual/solver.md",
                "manual/groebner.md",
            ],

            "Calculus" => [
                "manual/derivatives.md",
                "manual/integration.md",
                "manual/ode.md",
                "manual/taylor.md",
                "manual/sparsity_detection.md",
                "manual/limits.md",
            ],  
        ],
        "Comparison Against SymPy" => "comparison.md",
    ]
)

deploydocs(
   repo = "github.com/JuliaSymbolics/Symbolics.jl.git";
   push_preview = true
)
