using Documenter
using FRANZ

makedocs(
    sitename = "FRANZ",
    modules = [FRANZ],
    format = Documenter.HTML(
        edit_link = "main",
    ),
)

deploydocs(
    repo = "github.com/leonardromano/FRANZ.jl",
    devbranch = "main",
    push_preview = true,
)
