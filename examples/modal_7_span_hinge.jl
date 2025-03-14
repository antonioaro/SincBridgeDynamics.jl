using SincBridgeDynamics
using CairoMakie, LaTeXStrings

# set theme for plots
set_theme!(merge(sincdyn_theme, theme_latexfonts()))

# bridge data
L = [16, 6.2, 28.4, 6.2, 16, 6.2, 28.4, 6.2, 16, 6.2, 28.4, 6.2, 16]
E = 57.4e9
I = 0.245
μ = 5376.7

# support conditions
bc = ["pinned", "hinge", "roller", "roller", "hinge", "hinge", "roller", "roller", "hinge", "hinge",
    "roller", "roller", "hinge", "pinned"]

# bridge model
bridge = init_bridge_model(E, I, L, μ, bc)

# lower and upper frequencies
FREQB = 1.0
FREQE = 3.0

# modal analysis
modes = modal(bridge, FREQB, FREQE)

# plot
iλ = 1

fig = Figure(size=(600, 200))

axs = Axis(fig[1, 1],
    xlabel=latexstring("x \\, \\mathrm{[m]}"),
    ylabel=latexstring("\\phi "),
    yticks=WilkinsonTicks(3)
)

fn = round(modes.ω[iλ] / 2 / π, digits=2)
axs.title = latexstring("f_{$(iλ)}=$(fn)\\,\\mathrm{Hz}")

x = [0; cumsum(L)]
vlines!(axs, x[[1, 3, 4, 7, 8, 11, 12, 14]], color=:lightgray)
vlines!(axs, x[[2, 5, 6, 9, 10, 13]], color=:lightgray, linestyle=:dash)
lines!(axs, modes.x, modes.Φ[iλ])

fig_tit = "./figs/bridge-5-mode-$(iλ).pdf"
save(fig_tit, fig)

fig


