using Documenter

push!(LOAD_PATH, "../src/")
include("../src/NLS_module.jl")
using .NLS

makedocs(
    sitename = "Documentation NLS",  # Name of your documentation site
    modules = [NLS],                 # Modules to document
    format = Documenter.HTML(),            # Output format (HTML by default)
    pages = Any[
        "Home" => "index.md"              # Main page of the docs
    ]
)
