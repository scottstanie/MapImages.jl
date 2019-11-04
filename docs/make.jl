using Documenter, MapImages

makedocs(;
    modules=[MapImages],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/scottstanie/MapImages.jl/blob/{commit}{path}#L{line}",
    sitename="MapImages.jl",
    authors="scott <scott.stanie@gmail.com>",
    assets=String[],
)

deploydocs(
    repo = "github.com/scottstanie/MapImages.jl.git",
)
