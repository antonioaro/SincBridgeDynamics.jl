using SincBridgeDynamics
using CairoMakie, LaTeXStrings

# set theme for plots
set_theme!(merge(sincdyn_theme, theme_latexfonts()))

# bridge data
L = [16, 6.2, 28.4, 6.2, 16, 6.2, 28.4, 6.2, 16, 6.2, 28.4, 6.2, 16]
E = 57.4e9
I = 0.245
μ = 5376.7
ζ =0.01

# support conditions
bc = ["pinned", "hinge", "roller", "roller", "hinge", "hinge", "roller", "roller", "hinge", "hinge",
    "roller", "roller", "hinge", "pinned"]

# bridge model
bridge = init_bridge_model(E, I, L, μ, bc, ζ)

# minimum modes
MINMODE = 3

# lower and upper frequencies
FREQB = 0.1
FREQE = 3.0

# modal analysis
modes = modal(bridge, FREQB, FREQE)

# train models
train = init_train_model(TRAIN_MODELS[1])
Vt = 40

# solution to normal coordinate
iλ = 1
fn = modes.ω[iλ] / 2 / π

# time step and vector
Δt, NSTEPS = init_time_analysis(FREQB, FREQE, bridge.Lb, Vt, train.leng, MINMODE)

# preallocate
Y = zeros(NSTEPS)
∂²Y = zeros(NSTEPS)

Yt = zeros(NSTEPS)
∂²Yt = zeros(NSTEPS)

# algorithm
alg = HighSampling(; Δt=Δt)

# modal solution
_solver!(alg, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

# energy percentile
freq, E = energy_percentile(Yt, Δt)
_, ∂²E = energy_percentile(∂²Yt, Δt)

iE99 = findfirst(∂²E .> 0.99)

#------------------
# plot
#------------------

fig = Figure()

axs = Axis(fig[1, 1],
    xlabel=latexstring("f/f_$(iλ) "),
    ylabel=latexstring(" E_q ")
)

# natural frequency
vlines!(axs, 1, color=:lightgray)
text!(axs, 1 - 0.01, 105, text=latexstring("f_$(iλ)"), color=:black, rotation=π / 2)

# load frequency
fV = modes.λ[iλ, 1] * Vt / bridge.Le[1] / 2 / pi
vlines!(axs, fV / fn, color=:lightgray)
text!(axs, fV / fn - 0.01, 105, text=latexstring("f_V"), color=:black, rotation=π / 2)

# fE99
vlines!(axs, freq[iE99] ./ fn, color=:red)
text!(axs, freq[iE99] ./ fn + 0.1, 105, text=latexstring("f_{E_{99}}"), color=:red, rotation=π / 2)

# energy percentile
lines!(axs, freq ./ fn, E .* 100, label=latexstring("Y_{$(iλ)}(t)"))
lines!(axs, freq ./ fn, ∂²E .* 100, label=latexstring("\\ddot{Y}_{$(iλ)}(t)"))

#
limits!(axs, 0, 2, 0.0, 120)
axislegend(axs; position=:rb)

fig_tit = "./figs/bridge5_eq_$(iλ).pdf"
save(fig_tit, fig)

fig