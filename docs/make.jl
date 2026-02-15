using Documenter
using sfeee

makedocs(
    sitename = "sfeee",
    format = Documenter.HTML(),
    modules = [sfeee]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
