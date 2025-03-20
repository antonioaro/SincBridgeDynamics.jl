# SincBridgeDynamics.jl

![Julia](https://img.shields.io/badge/Julia-1.9+-blue)  
📡 Simulation of bridge dynamics using sinc-based methods.

## 🚀 Installation

You can install this package in Julia using:

```julia
using Pkg
Pkg.add(url="https://github.com/ONSAS/StructuralDynamicsODESolvers.jl")
Pkg.add(url="https://github.com/antonioaro/SincBridgeDynamics.jl")
```

## 📖 Usage

Import the package in Julia:

```julia
using SincBridgeDynamics
```

Example usage:

```julia
# bridge data
E=3.09e10
I=8.59
L=[40.0]
μ=34610.0
ζ = 0.01

# support conditions
bc = ["pinned", "pinned"]

# post_point
xB = [20.0]

# bridge model
bridge = init_bridge_model(E, I, L, μ, bc, ζ)

# minimum modes
MINMODE = 3

# lower and upper frequencies
FREQB = 2.0
FREQE = 30.0

# modal analysis
modes = modal(bridge, FREQB, FREQE; maxit=150)

# train models
train = init_train_model(TRAIN_MODELS[2]; distributed=true)
Vt = 80

# time step and vector
Δt, NSTEPS = init_time_analysis(FREQB, FREQE, bridge.Lb, Vt, train.leng, MINMODE)
t = collect(range(0.0, step=Δt, stop=Δt * (NSTEPS - 1)))

# algorithm
alg = TimeDiscrete(; Δt=Δt)

# solver
U, A = msup(alg, bridge, modes, train, MINMODE, FREQE, Vt, NSTEPS, xB)
```

## 🛠 Development

If you want to contribute, clone the repository and use Julia's development mode:

```julia
Pkg.dev("SincBridgeDynamics")
```

To run tests (not available):

```julia
using Pkg
Pkg.test("SincBridgeDynamics")
```

## 📜 License

This project is licensed under the MIT License.

---
🤖 **Built with Julia ❤️**

