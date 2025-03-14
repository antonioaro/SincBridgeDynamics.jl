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
ζ = 0.01

# support conditions
bc = repeat(["roller"], 29)
bc[1] = "pinned"
bc[end] = "pinned"

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
train = init_train_model(TRAIN_MODELS[2]; distributed=true)
Vt = 80

# solution to normal coordinate
iλ = 2
fn = modes.ω[iλ] / 2 / π

# time step and vector
Δt, NSTEPS = init_time_analysis(FREQB, FREQE, bridge.Lb, Vt, train.leng, MINMODE)

# preallocate
Y = zeros(NSTEPS)
∂²Y = zeros(NSTEPS)

Yt = zeros(NSTEPS)
∂²Yt = zeros(NSTEPS)

#------------------
# plot
#------------------

fig = Figure()

axs = Axis(fig[1, 1],
    xlabel=latexstring("f/f_$(iλ) "),
    ylabel=latexstring("|\\ddot{Y}_$(iλ)(ιω)| \\, \\mathrm{[m/s^2/Hz]}"),
    yscale=log10
)

# natural frequency
vlines!(axs, 1, color=:lightgray)
text!(axs, 1 - 0.01, 1, text=latexstring("f_$(iλ)"), color=:black, rotation=π / 2)

# load frequency
fV = modes.λ[1, iλ] * Vt / bridge.Le[1] / 2 / pi
vlines!(axs, fV / fn, color=:lightgray)
text!(axs, fV / fn - 0.01, 1, text=latexstring("f_V"), color=:black, rotation=π / 2)

# bogie frequency
I = argmax(train.DTS)
fb = Vt / train.ν[I]
vlines!(axs, fb / fn, color=:lightgray)
text!(axs, fb / fn - 0.01, 1, text=latexstring("f_b"), color=:black, rotation=π / 2)

# band limit
Ω = _band_limit(bridge, modes, train, Vt, iλ)
vlines!(axs, Ω / fn, color=:lightgray)
text!(axs, Ω / fn - 0.01, 1, text=latexstring("Ω"), color=:black, rotation=π / 2)

#------------------
# high-sampling
#------------------

# algorithm
alg = HighSampling(; Δt=Δt)

# modal solution
_solver!(alg, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

# fourier transform
freq, h∂²Yt = rffts(∂²Yt, Δt)
l2 = lines!(axs, freq ./ fn, abs.(h∂²Yt), color=Cycled(2), label="High sampling")

# ------------------
# time discrete analysis
# ------------------

# algorithm
alg = TimeDiscrete(; Δt=Δt)

# modal solution
_solver!(alg, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

# fourier transform
freq, h∂²Yt = rffts(∂²Yt, Δt)
l1 = lines!(axs, freq ./ fn, abs.(h∂²Yt), color=Cycled(1), label="Resampled")

# limits
limits!(axs, 0, 5, 1e-7, 1e1)
axislegend(axs, [l2, l1], ["Resampled", "High sampling"], position=:lb)

# save
fig_tit = "./figs/bridge6_ha_$(iλ).pdf"
save(fig_tit, fig)

fig