using Pkg
using Documenter
using SequenceVariation
using Revise

# see https://github.com/tlienart/LiveServer.jl/issues/140#issuecomment-1271591251
Revise.revise()

makedocs(;
    checkdocs=:exports,
    linkcheck=true,
    sitename="SequenceVariation.jl",
    format=Documenter.HTML(),
    modules=[SequenceVariation],
    pages=[
        "Home" => "index.md",
        "Working with haplotypes" => "haplotypes.md",
        "Working with variations" => "variations.md",
        "Comparing variations" => "compare.md",
        "API Reference" => "api.md",
    ],
    authors=replace(
        join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => ""
    ) * ", The BioJulia Organisation, and other contributors.",
)

deploydocs(;
    repo="github.com/BioJulia/SequenceVariation.jl.git",
    devbranch="master",
    push_preview=true,
)
