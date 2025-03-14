using RailBridgeDynamics
using CairoMakie, LaTeXStrings

# set theme for plots
set_theme!(merge(sincdyn_theme, theme_latexfonts()))

# train models
train = init_train_model(TRAIN_MODELS[2]; distributed=true)

# ------------------
# figure
# ------------------

fig1 = Figure()

ax1 = Axis(fig1[1, 1],
    xlabel=latexstring("ν \\, \\mathrm{[m]}"),
    ylabel=latexstring("G(ν) \\, \\mathrm{[kN]}")
)

lines!(ax1, train.ν, train.DTS .* 1e-3)
limits!(3, 30, 0, nothing)

fig_tit = "./figs/dts_a1.pdf"
save(fig_tit, fig1)


fig2 = Figure()

ax2 = Axis(fig2[1, 1],
    xlabel=latexstring("x \\, \\mathrm{[m]}"),
    ylabel=latexstring("P \\, \\mathrm{[kN]}")
)

stem!(ax2, train.axle, train.load .* 1e-3)
limits!(nothing, nothing, 0, nothing)

fig_tit = "./figs/stem_a1.pdf"
save(fig_tit, fig2)

fig2