using Documenter
using PCPOP

DocMeta.setdocmeta!(PCPOP, :DocTestSetup, :(using PCPOP); recursive=true)

makedocs(
    sitename = "PCPOP.jl",
    modules  = [PCPOP],
    authors  = "Abhishek Mishra, Moisés Bermejo Morán",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://abh1mishra.github.io/PCPOP.jl",
        edit_link  = "main",
    ),
    repo = Documenter.Remotes.GitHub("abh1mishra", "PCPOP.jl"),
    pages = [
        "Home" => "index.md",
        "API"  => "api.md",
    ],
    # Keep the build from failing on undocumented exports / broken links while
    # docstring coverage is still incomplete. Tighten (remove) once complete.
    warnonly = true,
)

deploydocs(
    repo = "github.com/abh1mishra/PCPOP.jl.git",
    devbranch = "main",
    push_preview = true,
)
