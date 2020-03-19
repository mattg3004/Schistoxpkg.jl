using Documenter, Schistoxpkg

makedocs(;
    modules=[Schistoxpkg],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/mattg3004/Schistoxpkg.jl/blob/{commit}{path}#L{line}",
    sitename="Schistoxpkg.jl",
    authors="Matthew",
    assets=String[],
)

deploydocs(;
    repo="github.com/mattg3004/Schistoxpkg.jl",
)
