using SincBridgeDynamics
using CairoMakie, LaTeXStrings

# set theme for plots
set_theme!(merge(sincdyn_theme, theme_latexfonts()))

# bridge data
E = 3.09e10
I = 8.59
L = [0.8, 40.0, 0.8]
μ = 34610.0
ζ = 0.01

# support conditions
bc = ["free", "roller", "roller", "free"]

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
Vt = 80

# solution to normal coordinate
iλ = 1
fn = modes.ω[iλ] / 2 / π

# time step and vector
Δt, NSTEPS = init_time_analysis(FREQB, FREQE, bridge.Lb, Vt, train.leng, MINMODE)
t = collect(range(0.0, step=Δt, stop=Δt * (NSTEPS - 1)))

# preallocate
Y = zeros(NSTEPS)
∂²Y = zeros(NSTEPS)

Yt = zeros(NSTEPS)
∂²Yt = zeros(NSTEPS)

#------------------
# plot
#------------------

# time history of displacement -------------
fig1 = Figure()
ax1 = Axis(fig1[1, 1],
    xlabel=latexstring("t \\, \\mathrm{[s]}"),
    ylabel=latexstring("Y_{$(iλ)}(t) \\, \\mathrm{[m]}"))

# time history of acceleration -------------
fig2 = Figure()
ax2 = Axis(fig2[1, 1],
    xlabel=latexstring("t \\, \\mathrm{[s]}"),
    ylabel=latexstring("\\ddot{Y}_{$(iλ)}(t) \\, \\mathrm{[m/s^2]}"))


# frequency content of displacement -------------
fig3 = Figure()

ax3 = Axis(fig3[1, 1],
    xlabel=latexstring("f/f_$(iλ) "),
    ylabel=latexstring("|Y_$(iλ)(ιω)| \\, \\mathrm{[m/Hz]}"),
    yscale=log10
)

# natural frequency
vlines!(ax3, 1, color=:lightgray)
text!(ax3, 1 - 0.01, 2e-8, text=latexstring("f_$(iλ)"), color=:black, rotation=π / 2)

# load frequency
fV = modes.λ[iλ, 1] * Vt / bridge.Le[1] / 2 / pi
vlines!(ax3, fV / fn, color=:lightgray)
text!(ax3, fV / fn - 0.01, 2e-8, text=latexstring("f_V"), color=:black, rotation=π / 2)

# band limit
Ω = _band_limit(bridge, modes, train, Vt, iλ)
vlines!(ax3, Ω / fn, color=:lightgray)
text!(ax3, Ω / fn - 0.01, 2e-8, text=latexstring("Ω"), color=:black, rotation=π / 2)

# frequency content of acceleration -------------
fig4 = Figure()

ax4 = Axis(fig4[1, 1],
    xlabel=latexstring("f/f_$(iλ) "),
    ylabel=latexstring("|\\ddot{Y}_$(iλ)(ιω)| \\, \\mathrm{[m/s^2/Hz]}"),
    yscale=log10
)

# natural frequency
vlines!(ax4, 1, color=:lightgray)
text!(ax4, 1 - 0.01, 1e-5, text=latexstring("f_$(iλ)"), color=:black, rotation=π / 2)

# load frequency
fV = modes.λ[iλ, 1] * Vt / bridge.Le[1] / 2 / pi
vlines!(ax4, fV / fn, color=:lightgray)
text!(ax4, fV / fn - 0.01, 1e-5, text=latexstring("f_V"), color=:black, rotation=π / 2)

# band limit
Ω = _band_limit(bridge, modes, train, Vt, iλ)
vlines!(ax4, Ω / fn, color=:lightgray)
text!(ax4, Ω / fn - 0.01, 1e-5, text=latexstring("Ω"), color=:black, rotation=π / 2)

# ------------------
# low sampling
# ------------------

# algorithm
alg = HighSampling(; Δt=1 / (2.56 * Ω))

# modal solution
_solver!(alg, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

# plot
t_ = collect(range(0.0, step=alg.Δt, stop=alg.Δt * (NSTEPS - 1)))

scatter!(ax1, t_, -Yt, label="Low sampling")
scatter!(ax2, t_, ∂²Yt, label="Low sampling")

# ------------------
# time discrete analysis
# ------------------

# algorithm
alg = TimeDiscrete(; Δt=Δt)

# modal solution
_solver!(alg, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

# fourier transform
freq, hYtd = rffts(Yt, Δt)
freq, h∂²Ytd = rffts(∂²Yt, Δt)

# plot
lines!(ax1, t, -Yt, label="Resampled")
lines!(ax2, t, ∂²Yt, label="Resampled")

#------------------
# high-sampling
#------------------

# algorithm
alg = HighSampling(; Δt=Δt)

# modal solution
_solver!(alg, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

# fourier transform
freq, hYt = rffts(Yt, Δt)
freq, h∂²Yt = rffts(∂²Yt, Δt)

# plot
lines!(ax1, t, -Yt, label="High sampling", linestyle=:dash)
lines!(ax2, t, ∂²Yt, label="High sampling", linestyle=:dash)

l2 = lines!(ax3, freq ./ fn, abs.(hYt), color=Cycled(2), label="High sampling")
l2 = lines!(ax4, freq ./ fn, abs.(h∂²Yt), color=Cycled(2), label="High sampling")

l1 = lines!(ax3, freq ./ fn, abs.(hYtd), color=Cycled(1), label="Resampled")
l1 = lines!(ax4, freq ./ fn, abs.(h∂²Ytd), color=Cycled(1), label="Resampled")

# limits
limits!(ax1, 0, 2, -1.75e-8, 0.5e-8)
limits!(ax2, 0, 2, -2.5e-6, 1.5e-6)
limits!(ax3, 0, 5, 1e-15, 1e-7)
limits!(ax4, 0, 5, 1e-11, 1e-4)

# legend
axislegend(ax1, position=:rb)
axislegend(ax2, position=:rb)
axislegend(ax3, [l1, l2], ["Resampled", "High sampling"], position=:rt)
axislegend(ax4, [l1, l2], ["Resampled", "High sampling"], position=:rt)

# save
fig_tit = "./figs/bridge2_u_$(iλ).pdf"
save(fig_tit, fig1)

fig_tit = "./figs/bridge2_a_$(iλ).pdf"
save(fig_tit, fig2)

fig_tit = "./figs/bridge2_hu_$(iλ).pdf"
save(fig_tit, fig3)

fig_tit = "./figs/bridge2_ha_$(iλ).pdf"
save(fig_tit, fig4)

fig4