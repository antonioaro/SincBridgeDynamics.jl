sincdyn_theme = Theme(
    size=(600, 300),
    Axis=(
        xgridvisible=false,
        ygridvisible=false,
        xlabelsize=22,
        ylabelsize=22,
        xtickalign=1,
        xticksize=10,
        ytickalign=1,
        yticksize=10,
        xticklabelsize=16,
        yticklabelsize=16,
        # yticks=WilkinsonTicks(5; k_min=5)
    ),
    Legend=(framecolor=(:white, 0.0),
        nbanks=2,
        backgroundcolor=(:white, 0.85),
        labelsize=16
    ),
    fontsize=18
)