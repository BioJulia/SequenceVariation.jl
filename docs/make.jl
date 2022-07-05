using Documenter
using SequenceVariation

makedocs(
    sitename = "SequenceVariation",
    format = Documenter.HTML(),
    modules = [SequenceVariation]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
