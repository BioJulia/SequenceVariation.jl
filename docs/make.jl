using Pkg
using Documenter
using SequenceVariation

makedocs(;
    checkdocs = :exports,
    linkcheck = true,
    sitename = "SequenceVariation.jl",
    format = Documenter.HTML(),
    modules = [SequenceVariation],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/SequenceVariation.jl.git",
    devbranch = "master",
    push_preview = true,
)
