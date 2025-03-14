using SincBridgeDynamics
using CairoMakie, LaTeXStrings

# set theme for plots
set_theme!(merge(sincdyn_theme, theme_latexfonts()))

# bridge data
L = [43.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0,
    44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 43.0]
E = 30.0e9
I = 18.875
μ = 45860.0

# support conditions
bc = repeat(["roller"], 29)
bc[1] = "pinned"
bc[end] = "pinned"

# bridge model
bridge = init_bridge_model(E, I, L, μ, bc)

# lower and upper frequencies
FREQB = 1.0
FREQE = 3.0

# modal analysis
modes = modal(bridge, FREQB, FREQE)

# plot
iλ = 2

fig = Figure(size=(600, 200))

axs = Axis(fig[1, 1],
    xlabel=latexstring("x \\, \\mathrm{[m]}"),
    ylabel=latexstring("\\phi "),
    yticks=WilkinsonTicks(3)
)

fn = round(modes.ω[iλ] / 2 / π, digits=2)
axs.title = latexstring("f_{$(iλ)}=$(fn)\\,\\mathrm{Hz}")

x = [0; cumsum(L)]
vlines!(axs, x, color=:lightgray)
lines!(axs, modes.x, modes.Φ[iλ])

fig_tit = "./figs/bridge-6-mode-$(iλ).pdf"
save(fig_tit, fig)

fig


