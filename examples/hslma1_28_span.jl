using SincBridgeDynamics
using FFTW
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

# post_point
xB = [1082.50]

# minimum modes
MINMODE = 3

# lower and upper frequencies
FREQB = 0.1
FREQE = 30.0

# modal analysis
modes = modal(bridge, FREQB, FREQE; maxit=200)

# train models
train = init_train_model(TRAIN_MODELS[2]; distributed=true)
Vt = 80

# time step and vector
Δt, NSTEPS = init_time_analysis(FREQB, FREQE, bridge.Lb, Vt, train.leng, MINMODE)
t = collect(range(0.0, step=Δt, stop=Δt * (NSTEPS - 1)))

# figures
fig1 = Figure()
ax1 = Axis(fig1[1, 1],
    xlabel=latexstring("t \\, \\mathrm{[s]}"),
    ylabel=latexstring("u(t) \\, \\mathrm{[m]}"))

fig2 = Figure()
ax2 = Axis(fig2[1, 1],
    xlabel=latexstring("t \\, \\mathrm{[s]}"),
    ylabel=latexstring("\\ddot{u}(t) \\, \\mathrm{[m/s^2]}"))

fig3 = Figure()
ax3 = Axis(fig3[1, 1],
    xlabel=latexstring("f \\, \\mathrm{[Hz]}"),
    ylabel=latexstring("|u(ιω)| \\, \\mathrm{[m/Hz]}"),
    yscale=log10)

fig4 = Figure()
ax4 = Axis(fig4[1, 1],
    xlabel=latexstring("f \\, \\mathrm{[Hz]}"),
    ylabel=latexstring("|\\ddot{u}(ιω)| \\, \\mathrm{[m/s^2/Hz]}"),
    yscale=log10)

# ------------------
# time discrete analysis
# ------------------

# algorithm
alg = TimeDiscrete(; Δt=Δt)

# normal coordinate
@time U, A = msup(alg, bridge, modes, train, MINMODE, FREQE, Vt, NSTEPS, xB)

# plot
lines!(ax1, t, -U[:], label="Resampled")
lines!(ax2, t, -A[:], label="Resampled")

# fourier transform
freq, hU = rffts(U[:], Δt)
_, hA = rffts(A[:], Δt)

lines!(ax3, freq, abs.(hU), label="Resampled")
lines!(ax4, freq, abs.(hA), label="Resampled")

# ------------------
# high-sampling
# ------------------

# algorithm
alg = HighSampling(; Δt=Δt)

# normal coordinate
@time U, A = msup(alg, bridge, modes, train, MINMODE, FREQE, Vt, NSTEPS, xB)

# plot
lines!(ax1, t, -U[:], label="High sampling", linestyle=:dash)
lines!(ax2, t, -A[:], label="High sampling", linestyle=:dash)

# fourier transform
freq, hU = rffts(U[:], Δt)
_, hA = rffts(A[:], Δt)

lines!(ax3, freq, abs.(hU), label="High sampling", linestyle=:dash)
lines!(ax4, freq, abs.(hA), label="High sampling", linestyle=:dash)

# lims
limits!(ax1, 5, 28, -0.002, 0.001)
limits!(ax2, 5, 28, -0.1, 0.1)
limits!(ax3, 0, 30, 1e-10, 1e-2)
limits!(ax4, 0, 30, 1e-7, 1e-0)

# legend
axislegend(ax1, position=:lb,)
axislegend(ax2, position=:lb,)
axislegend(ax3, position=:rt,)
axislegend(ax4, position=:rt,)

# save 
fig_tit = "./figs/bridge6_disp.pdf"
save(fig_tit, fig1)

fig_tit = "./figs/bridge6_acc.pdf"
save(fig_tit, fig2)

fig_tit = "./figs/bridge6_hu.pdf"
save(fig_tit, fig3)

fig_tit = "./figs/bridge6_ha.pdf"
save(fig_tit, fig4)

fig4