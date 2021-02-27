using Documenter, Symbolics

makedocs(
    sitename="Symbolics.jl",
    authors="Chris Rackauckas",
    modules=[Symbolics],
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
            "variables.md",
            "expression_manipulation.md",
            "derivatives.md",
            "build_function.md",
            "functions.md",
            "io.md",
            "sparsity_detection.md"
        ],
        "Comparison Against SymPy" => "comparison.md",
    ]
)

deploydocs(
   repo = "github.com/JuliaSymbolics/Symbolics.jl.git";
   push_preview = true
)
