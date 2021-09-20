using Documenter
using MOSLab

makedocs(
    sitename = "MOSLab",
    format = Documenter.HTML(),
    modules = [MOSLab]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/Rapos0/MOSLab.jl.git"
)
