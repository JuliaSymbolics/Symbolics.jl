using Documenter, Symbolics, SymbolicUtils

# Make sure that plots don't throw a bunch of warnings / errors!
ENV["GKSwstype"] = "100"
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
    clean=true,doctest=false,
    strict=[
        :doctest,
        :linkcheck,
        :parse_error,
        :example_block,
        # Other available options are
        # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
    ],
    format = Documenter.HTML(#analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"]),
                             mathengine = mathengine,
                             #canonical="https://mtk.sciml.ai/stable/"),
    pages=[
        "Home" => "index.md",
        "tutorials/symbolic_functions.md",
        "Tutorials" => Any[
            "tutorials/auto_parallel.md",
            "tutorials/converting_to_C.md"
        ],
        "Examples" => Any[
            "examples/perturbation.md"
        ],
        "Manual" => Any[
            "manual/variables.md",
            "manual/expression_manipulation.md",
            "manual/derivatives.md",
            "manual/groebner.md",
            "manual/arrays.md",
            "manual/build_function.md",
            "manual/functions.md",
            "manual/io.md",
            "manual/sparsity_detection.md",
            "manual/types.md",
            "manual/faq.md"
        ],
        "Comparison Against SymPy" => "comparison.md",
    ]
)

deploydocs(
   repo = "github.com/JuliaSymbolics/Symbolics.jl.git";
   push_preview = true
)
