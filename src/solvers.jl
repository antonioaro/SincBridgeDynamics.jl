"""
    msup(alg, bridge, modes, train, MINMODE, FREQE, Vt, NSTEPS, xB) -> (U, A)

Performs mode superposition analysis for a moving load problem.

# Parameters
- `alg`: Solver algorithm for the modal analysis.
- `bridge`: Bridge model.
- `modes`: Modal properties of the bridge.
- `train`: Train model.
- `MINMODE::Int`: Minimum mode index to be considered in the analysis.
- `FREQE::Float64`: Upper frequency limit for mode superposition.
- `Vt::Float64`: Train velocity.
- `NSTEPS::Int`: Number of time steps.
- `xB::Vector{Float64}`: Positions where responses are computed.

# Returns
- `U::Matrix{Float64}`: Displacement response at `xB` locations over time.
- `A::Matrix{Float64}`: Acceleration response at `xB` locations over time.
"""
function msup(alg, bridge, modes, train, MINMODE, FREQE, Vt, NSTEPS, xB)

    U = zeros(NSTEPS, length(xB))
    A = zeros(NSTEPS, length(xB))

    Y = zeros(NSTEPS)
    ∂²Y = zeros(NSTEPS)

    Yt = zeros(NSTEPS)
    ∂²Yt = zeros(NSTEPS)

    # loop over mode shapes
    for iλ in eachindex(modes.ω)
        if modes.ω[iλ] .<= (2.0 * π * FREQE) || iλ <= MINMODE

            # clean buffer
            Yt .= 0.0
            ∂²Yt .= 0.0

            Y .= 0.0
            ∂²Y .= 0.0

            # modal solution            
            _solver!(alg, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

            # mode superposition
            _modesupeporsition!(U, A, Yt, ∂²Yt, bridge, modes, xB, iλ)
        end
    end

    return U, A
end

# normal coordinate
function _solver!(alg::AbstractSolver, bridge, modes, train, Yt, ∂²Yt, Y, ∂²Y, Vt, NSTEPS, iλ)

    # distributed load
    x_sl = train.distributed ? SVector(0.0, 0.6, 0.6) : SVector(0.0)
    p_sl = train.distributed ? SVector(0.25, 0.5, 0.25) : SVector(1.0)

    # time step for discrete analysis
    if alg isa Analytic
        if alg.discrete
            alg.band_limit = _band_limit(bridge, modes, train, Vt, iλ)
            _init_discrete_analysis!(alg, NSTEPS)

            # FFTW plan
            rfft_plan, irfft_plan, x, X = _fftw_plans(Yt, alg.NSTEPS, NSTEPS, alg.new_size)
        end
    end

    # loop over axle load
    Yt .= 0.0
    ∂²Yt .= 0.0

    @fastmath @inbounds @simd for i_axle in 1:train.naxl
        σ = train.axle[i_axle] / Vt
        p_axle = train.load[i_axle]

        for i_sl in eachindex(x_sl)
            σ += x_sl[i_sl] / Vt

            # solution
            Y, ∂²Y = _solver(alg, NSTEPS, Vt, Y, ∂²Y, modes, bridge, iλ; σ)

            # sum
            @views @. Yt += Y * p_axle * p_sl[i_sl]
            @views @. ∂²Yt += ∂²Y * p_axle * p_sl[i_sl]
        end
    end

    if alg isa Analytic
        if alg.discrete
            _fft_resample!(Yt, alg.NSTEPS, NSTEPS, alg.new_size, x, X, rfft_plan, irfft_plan)
            _fft_resample!(∂²Yt, alg.NSTEPS, NSTEPS, alg.new_size, x, X, rfft_plan, irfft_plan)
        end
    end

    nothing
end

function _modesupeporsition!(U, A, Yt, ∂²Yt, bridge, modes, post_point, iλ)

    # mode shape
    @views Ψ = _modeigen(modes.λ[iλ, :], modes.v[:, iλ], bridge.Le, post_point)

    # mode superposition
    @views mul!(U, Yt, Ψ', 1.0, 1.0)
    @views mul!(A, ∂²Yt, Ψ', 1.0, 1.0)

    nothing
end

function _fft_resample!(Y::AbstractVector{T}, n::Int, NSTEPS::Int, new_size, x, X, rfft_plan, irfft_plan) where {T}

    # do nothing
    if n == NSTEPS
        return x
    end

    # Clear real and complex buffer
    x .= 0.0
    X .= 0.0

    # rfft: Apply real-to-complex FFT to the real-valued input x
    nfft = n ÷ 2 + 1
    @views mul!(X[1:nfft], rfft_plan, Y[1:n])

    # irfft: Apply complex-to-real inverse FFT to the complex buffer X
    @views mul!(x, irfft_plan, X)
    @views @. Y = x[1:NSTEPS] * new_size / n

    nothing
end

function _fftw_plans(arr::AbstractVector, N, NSTEPS, new_size)

    # Real and Complex buffer for FFT results
    x = Vector{Float64}(undef, new_size)
    X = Vector{ComplexF64}(undef, new_size ÷ 2 + 1)

    # FFTW plans
    rfft_plan = plan_rfft(@view arr[1:N]; flags=FFTW.ESTIMATE)
    irfft_plan = plan_irfft(X, new_size; flags=FFTW.ESTIMATE)

    return rfft_plan, irfft_plan, x, X
end