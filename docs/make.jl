using Documenter, RandomMatrixDistributions

makedocs(
    sitename = "RandomMatrixDistributions.jl",
    authors = "Damian Pavlyshyn",
    pages = [
        "Home" => "index.md",
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/damian-t-p/RandomMatrixDistributions.jl.git",
)
