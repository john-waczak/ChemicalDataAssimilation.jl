using ChemicalDataAssimilation
using Documenter

DocMeta.setdocmeta!(ChemicalDataAssimilation, :DocTestSetup, :(using ChemicalDataAssimilation); recursive=true)

makedocs(;
    modules=[ChemicalDataAssimilation],
    authors="John Waczak <john.louis.waczak@gmail.com>",
    repo="https://github.com/john-waczak/ChemicalDataAssimilation.jl/blob/{commit}{path}#{line}",
    sitename="ChemicalDataAssimilation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://john-waczak.github.io/ChemicalDataAssimilation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/john-waczak/ChemicalDataAssimilation.jl",
    devbranch="main",
)
