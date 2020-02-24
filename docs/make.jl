using Documenter, CurvilinearCalculus

makedocs(;
    modules=[CurvilinearCalculus],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/fgerick/CurvilinearCalculus.jl/blob/{commit}{path}#L{line}",
    sitename="CurvilinearCalculus.jl",
    authors="Felix Gerick <felixgerick@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/fgerick/CurvilinearCalculus.jl",
)
