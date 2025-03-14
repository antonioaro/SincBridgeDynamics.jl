using SincBridgeDynamics
using CairoMakie, LaTeXStrings

# set theme for plots
set_theme!(merge(sincdyn_theme, theme_latexfonts()))

#----------------------
# 3-span fixed
#----------------------

# bridge data
L = [17.4, 20.3, 19.3]
E = 34e9
I = 1.513
μ = 18680

# support conditions
bc = ["fixed", "roller", "roller", "fixed"]

# bridge model
bridge = init_bridge_model(E, I, L, μ, bc)

# lower and upper frequencies
FREQB = 1.0
FREQE = 15.0

# modal analysis
modes = modal(bridge, FREQB, FREQE)

#----------------------
# 3-span oversails
#----------------------

# bridge data
L = [2.0, 17.4, 20.3, 19.3, 2.0]

# support conditions
bc = ["fixed", "roller", "roller", "roller", "roller", "fixed"]

# bridge model
bridge_ = init_bridge_model(E, I, L, μ, bc)

# modal analysis
modes_ = modal(bridge_, FREQB, FREQE)

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

# plot
x = [0; cumsum(L)]
vlines!(axs, x[[2,3,4,5]], color=:lightgray)
vlines!(axs, x[[1, end]], color=:lightgray, linestyle=:dash)

lines!(axs, modes_.x, -modes_.Φ[iλ], label="Oversails")
lines!(axs, modes.x .+ 2, -modes.Φ[iλ], linestyle=:dash, label="Fixed")

ylims!(axs, -0.0025, 0.0025)
axislegend(; position=:lb)

fig_tit = "./figs/bridge-4-mode-$(iλ).pdf"
save(fig_tit, fig)

fig