# When building docs, set DOCS so package code can detect we're in a docs
# build and avoid running top-level side-effects (packages should guard
# heavy side-effects with `if get(ENV, "DOCS", "") != "true"`).
ENV["DOCS"] = "true"

using Documenter
using MixSwitchStochEpi

makedocs(
    sitename = "MixSwitchStochEpi.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [MixSwitchStochEpi],
    pages = [
        "Home" => "index.md",
                "API Reference" => [
            "All API" => "api_all.md",
            "Activity Scoring" => "api/activity_scoring.md",
            "Contact Models" => "api/contact_models.md",
            "Epidemic Dynamics" => "api/epidemic_dynamics.md",
            "Utilities" => "api/utilities.md",
        ],
    ]
)