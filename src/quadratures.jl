struct Numerical{N,T} <: AbstractSolver
    Δt::N
    integration_scheme::T
end

Numerical(Δt::N, integration_scheme::T) where {N,T} = Numerical{N,T}(Δt, integration_scheme)

# time integration scheme
function _solver(alg::Numerical{N,T}, NSTEPS::Int, Vt, Y, ∂²Y, modes, bridge, iλ; σ=0.0) where {N,T}

    Δt = alg.Δt

    # modal force
    x = (range(0.0, step=Δt, stop=Δt * (NSTEPS - 1)) .- σ) .* Vt
    @views fΨ = _modeigen(modes.λ[iλ, :], modes.v[:, iλ], bridge.Le, x)
    fΨ .*= 1 / modes.Mn[iλ]

    # numerical integration
    m = 1.0
    c = 2.0 * m * bridge.ζ * modes.ω[iλ]
    k = modes.ω[iλ]^2

    M = m * ones(1, 1)
    C = c * ones(1, 1)
    K = k * ones(1, 1)
    R = [[fΨ[i]] for i in eachindex(fΨ)]

    X = nothing # state constraints are ignored
    B = ones(1, 1)
    sys = SecondOrderConstrainedLinearControlContinuousSystem(M, C, K, B, X, R)

    U₀ = zeros(1)
    V₀ = zeros(1)

    ivp_free = InitialValueProblem(sys, (U₀, V₀))

    alg_ = alg.integration_scheme(Δt)
    sol = solve(ivp_free, alg_; NSTEPS=NSTEPS - 1)

    return displacements(sol, 1), accelerations(sol, 1)
end