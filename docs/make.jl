using Documenter, RandomMatrixDistributions

makedocs(
    sitename = "RandomMatrixDistributions.jl",
    authors = "Damian Pavlyshyn",
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Implementation" => "implementation.md",
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/damian-t-p/RandomMatrixDistributions.jl.git",
)
