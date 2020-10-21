using Documenter, Schistoxpkg

makedocs(
    modules=[Schistoxpkg],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Functions" => Any[
           "man/functions.md",
         ],
         "Parameters" => Any[
           "man/parameters.md",
         ]
    ],
    sitename="Schistoxpkg.jl",
    authors="Matthew",
)

deploydocs(
    repo="github.com/mattg3004/Schistoxpkg.jl.git",
    target = "build"
)
