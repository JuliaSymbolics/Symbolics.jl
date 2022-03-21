using Documenter, Symbolics, SymbolicUtils

makedocs(
    sitename="Symbolics.jl",
    authors="Chris Rackauckas",
    modules=[Symbolics,SymbolicUtils],
    clean=true,doctest=false,
    format = Documenter.HTML(#analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"]),
                             #canonical="https://mtk.sciml.ai/stable/"),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorials/symbolic_functions.md",
            "tutorials/auto_parallel.md",
            "tutorials/converting_to_C.md"
        ],
        "Manual" => Any[
            "manual/variables.md",
            "manual/expression_manipulation.md",
            "manual/derivatives.md",
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
