using Documenter, BGEN

makedocs(
    format = Documenter.HTML(),
    sitename = "BGEN.jl",
    authors = "Seyoon Ko",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/BGEN.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    devbranch = "main"
)
