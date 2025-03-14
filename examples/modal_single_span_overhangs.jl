using SincBridgeDynamics
using CairoMakie, LaTeXStrings

# set theme for plots
set_theme!(merge(sincdyn_theme, theme_latexfonts()))

#----------------------
# single-span
#----------------------

# bridge data
E=3.09e10
I=8.59
L=[40.0]
μ=34610.0

# support conditions
bc = ["pinned", "pinned"]

# bridge model
bridge = init_bridge_model(E, I, L, μ, bc)

# lower and upper frequencies
FREQB = 2.0
FREQE = 30.0

# modal analysis
modes = modal(bridge, FREQB, FREQE; maxit=150)

#----------------------
# simply supported with overhangs
#----------------------

# bridge data
L=[0.8, 40.0, 0.8]

# support conditions
bc = ["free", "roller", "roller", "free"]

# bridge model
bridge_ = init_bridge_model(E, I, L, μ, bc)

# modal analysis
modes_ = modal(bridge_, FREQB, FREQE)

# plot
iλ = 3

fig = Figure(size=(600, 200))

axs = Axis(fig[1, 1],
    xlabel=latexstring("x \\, \\mathrm{[m]}"),
    ylabel=latexstring("\\phi "),
    yticks=WilkinsonTicks(3)
)

fn = round(modes.ω[iλ] / 2 / π, digits=2)
axs.title = latexstring("f_{$(iλ)}=$(fn)\\,\\mathrm{Hz}")

# plot
x = [0; cumsum(L)]
vlines!(axs, x[[2,3]], color=:lightgray)

lines!(axs, modes_.x, modes_.Φ[iλ], label="Overhangs")
lines!(axs, modes.x .+ 0.8, modes.Φ[iλ], linestyle=:dash, label="Simply")

ylims!(axs, -0.002, 0.002)
axislegend(; position=:lb)

fig_tit = "./figs/bridge-2-mode-$(iλ).pdf"
save(fig_tit, fig)

fig