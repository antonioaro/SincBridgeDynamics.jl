module SincBridgeDynamics

using LinearAlgebra
using NonlinearEigenproblems
using SparseArrays
using SLEEF
using StaticArrays
using FFTW
using StructuralDynamicsODESolvers
using CairoMakie, LaTeXStrings

# add https://github.com/ONSAS/StructuralDynamicsODESolvers.jl

include("interface.jl")
export init_bridge_model

include("modal.jl")
include("constraints.jl")
export modal

include("damping.jl")
export ec1damping

include("plots.jl")
export sincdyn_theme

include("train_load.jl")
export init_train_model

include("TRAIN_MODELS.jl")
export TRAIN_MODELS

include("utils.jl")
export init_time_analysis, energy_percentile, rffts

include("time_discrete.jl")
export HighSampling, TimeDiscrete, _band_limit

include("solvers.jl")
export _solver!, msup

include("quadratures.jl")
export Numerical

end
