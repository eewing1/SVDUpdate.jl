using SVDUpdate
using Documenter

DocMeta.setdocmeta!(SVDUpdate, :DocTestSetup, :(using SVDUpdate); recursive=true)

makedocs(;
    modules=[SVDUpdate],
    authors="Evan Ewing <eewing@mit.edu> and contributors",
    repo="https://github.com/eewing1/SVDUpdate.jl/blob/{commit}{path}#{line}",
    sitename="SVDUpdate.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://eewing1.github.io/SVDUpdate.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/eewing1/SVDUpdate.jl",
    devbranch="main",
)
