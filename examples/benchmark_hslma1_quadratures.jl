using SincBridgeDynamics
using StructuralDynamicsODESolvers
using BenchmarkTools
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
FREQB = 1.0
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

# ------------------
# trapezoidal
# ------------------

# algorithm
alg = Numerical(Δt, Trapezoidal)

# solver
@time U, A = msup(alg, bridge, modes, train, MINMODE, FREQE, Vt, NSTEPS, xB)

# plot
lines!(ax1, t, -U[:], label="Newmark")
lines!(ax2, t, -A[:], label="Newmark")

# ------------------
# bathe
# ------------------

# algorithm
alg = Numerical(Δt, Bathe)

# solver
@time U, A = msup(alg, bridge, modes, train, MINMODE, FREQE, Vt, NSTEPS, xB)

# plot
lines!(ax1, t, -U[:], label="Bathe")
lines!(ax2, t, -A[:], label="Bathe")

# lims
limits!(ax1, 18, 20, -0.002, 0.0005)
limits!(ax2, 18, 20, -0.1, 0.1)

# legend
axislegend(ax1, position=:rb,)
axislegend(ax2, position=:rb,)

# save 
fig_tit = "./figs/benchmark_u.pdf"
save(fig_tit, fig1)

fig_tit = "./figs/benchmark_d2u.pdf"
save(fig_tit, fig2)

fig2