using Symlens
using Documenter

DocMeta.setdocmeta!(Symlens, :DocTestSetup, :(using Symlens); recursive=true)

makedocs(;
    modules=[Symlens],
    authors="Yilun Guan",
    repo="https://github.com/guanyilun/laughing-succotash/blob/{commit}{path}#{line}",
    sitename="Symlens.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://guanyilun.github.io/Symlens.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/guanyilun/laughing-succotash",
    devbranch="master",
)
