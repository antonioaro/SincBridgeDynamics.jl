abstract type AbstractSolver end

mutable struct Analytic{N} <: AbstractSolver
    Δt::N
    discrete::Bool
    NSTEPS::Int
    new_size::Int
    band_limit::N
    Δt_::N
end

TimeDiscrete(; Δt::N) where {N} = Analytic(Δt, true, 0, 0, 0.0, 0.0)
HighSampling(; Δt::N) where {N} = Analytic{N}(Δt, false, 0, 0, 0.0, 0.0)

function _band_limit(bridge, modes, train, Vt, iλ)

    # natural frequency
    fn = modes.ω[iλ] / 2 / π

    # single load period
    fl = modes.λ[iλ, 1] * Vt / bridge.Le[1] / 2 / pi

    # bogie passaging frequency
    I = argmax(train.DTS)
    fb = train.naxl > 1 ? Vt / train.ν[I] : 0.0

    return 2 * max(fn, fl, fb) * 1.1
end

function _init_discrete_analysis!(alg, NSTEPS)
    alg.Δt_ = 1 / (alg.band_limit * 2.56)
    alg.NSTEPS = nextpow(2, (NSTEPS * alg.Δt) / alg.Δt_)
    alg.new_size = Int(cld(alg.NSTEPS * alg.Δt_, alg.Δt))

    nothing
end

function _solver(alg::Analytic{N}, NSTEPS::Int, Vt, Y, ∂²Y, modes, bridge, iλ; σ=0.0) where {N}
    Δt = alg.discrete ? alg.Δt_ : alg.Δt
    NSTEPS = alg.discrete ? alg.NSTEPS : NSTEPS

    λ, v, ω, Mn = modes.λ, modes.v, modes.ω, modes.Mn
    ζ, L, Ne = bridge.ζ, bridge.Le, bridge.Ne

    ωD = ω[iλ] * √(1 - ζ^2)
    scale = 1.0 / (ωD * Mn[iλ])

    x₀ = 0.0
    ∑A⁻ = ∑B⁻ = ∑∂²A⁻ = ∑∂²B⁻ = 0.0
    Y .= 0.0
    ∂²Y .= 0.0

    @inbounds for i in eachindex(L)
        α₁, α₂, α₃, α₄ = @view v[(i-1)*4+1:(i-1)*4+4, iλ]

        t₀ = x₀ / Vt + σ
        t⁺ = t₀ + L[i] / Vt

        A₀, B₀, ∂²A₀, ∂²B₀ = _integral_terms(λ[iλ, i], ω[iλ], ζ, α₁, α₂, α₃, α₄, L[i], Vt, 0.0, t₀)

        n₀, n⁺ = Int(cld(t₀, Δt)), Int(cld(t⁺, Δt))

        @simd for j in n₀:n⁺
            t = j * Δt
            A⁺, B⁺, ∂²A⁺, ∂²B⁺, Pn = _integral_terms(λ[iλ, i], ω[iλ], ζ, α₁, α₂, α₃, α₄, L[i], Vt, t - t₀, t₀)

            sin_, cos_ = SLEEF.sincos_fast(ωD * t)
            exp_ = SLEEF.exp(-ζ * ω[iλ] * t)

            Y[j+1] = ((A⁺ - A₀ + ∑A⁻) * sin_ - (B⁺ - B₀ + ∑B⁻) * cos_) * exp_ * scale
            ∂²Y[j+1] = (((∂²A⁺ - ∂²A₀ + ∑∂²A⁻) * sin_ + (∂²B⁺ - ∂²B₀ + ∑∂²B⁻) * cos_) * exp_ + ωD * Pn) * scale
        end

        if i == Ne
            @simd for j in n⁺:NSTEPS-1
                t = j * Δt
                sin_, cos_ = SLEEF.sincos_fast(ωD * t)
                exp_ = SLEEF.exp(-ζ * ω[iλ] * t)

                Y[j+1] = ((A⁺ - A₀ + ∑A⁻) * sin_ - (B⁺ - B₀ + ∑B⁻) * cos_) * exp_ * scale
                ∂²Y[j+1] = (((∂²A⁺ - ∂²A₀ + ∑∂²A⁻) * sin_ + (∂²B⁺ - ∂²B₀ + ∑∂²B⁻) * cos_) * exp_ + ωD * Pn) * scale
            end
        end

        x₀ += L[i]

        A⁺, B⁺, ∂²A⁺, ∂²B⁺, Pn = _integral_terms(λ[iλ, i], ω[iλ], ζ, α₁, α₂, α₃, α₄, L[i], Vt, t⁺ - t₀, t₀)

        ∑A⁻ += A⁺ - A₀
        ∑B⁻ += B⁺ - B₀
        ∑∂²A⁻ += ∂²A⁺ - ∂²A₀
        ∑∂²B⁻ += ∂²B⁺ - ∂²B₀
    end

    return Y, ∂²Y
end

function _integral_terms(λ, ω, ζ, α₁, α₂, α₃, α₄, L, Vt, t, T)

    # integration limit
    t̃ = t < L / Vt ? t : L / Vt

    # resonance frequency
    ωD = ω * √(1 - ζ^2)

    # constants
    A₀ = λ * Vt / L
    C = ζ * ω

    # trigonometric and hiperbolic functions
    sin⁻, cos⁻ = SLEEF.sincos_fast(A₀ * t̃ - ωD * (t̃ + T))
    sin⁺, cos⁺ = SLEEF.sincos_fast(A₀ * t̃ + ωD * (t̃ + T))
    sin_, cos_ = SLEEF.sincos_fast(ωD * (t̃ + T))
    sina_, cosa_ = SLEEF.sincos_fast(A₀ * t̃)

    cosh_ = cosh(A₀ * t̃)
    sinh_ = sinh(A₀ * t̃)

    # integral terms
    A₁ = (C * sin⁺ - (A₀ + ωD) * cos⁺) / (C^2 + (A₀ + ωD)^2)
    A₂ = (C * sin⁻ - (A₀ - ωD) * cos⁻) / (C^2 + (A₀ - ωD)^2)
    A₃ = (C * cos⁺ + (A₀ + ωD) * sin⁺) / (C^2 + (A₀ + ωD)^2)
    A₄ = (C * cos⁻ + (A₀ - ωD) * sin⁻) / (C^2 + (A₀ - ωD)^2)
    A₅ = ωD * (A₀^2 + C^2 + ωD^2) * sin_ - C * (A₀^2 - C^2 - ωD^2) * cos_
    A₆ = A₀ * ((A₀^2 - C^2 + ωD^2) * cos_ - 2 * C * ωD * sin_)
    A₇ = ωD * (A₀^2 + C^2 + ωD^2) * cos_ + C * (A₀^2 - C^2 - ωD^2) * sin_
    A₈ = A₀ * ((A₀^2 - C^2 + ωD^2) * sin_ + 2 * C * ωD * cos_)

    D = (A₀^2 + C^2 + ωD^2)^2 - (2 * A₀ * C)^2
    E = SLEEF.exp(C * (t̃ + T))

    I₁ = (A₁ + A₂) / 2
    I₂ = (A₃ + A₄) / 2
    I₃ = (A₅ * sinh_ + A₆ * cosh_) / D
    I₄ = (A₅ * cosh_ + A₆ * sinh_) / D

    J₁ = (-A₃ + A₄) / 2
    J₂ = (A₁ - A₂) / 2
    J₃ = (-A₇ * sinh_ + A₈ * cosh_) / D
    J₄ = (-A₇ * cosh_ + A₈ * sinh_) / D

    An = E * (α₁ * I₁ + α₂ * I₂ + α₃ * I₃ + α₄ * I₄)
    Bn = E * (α₁ * J₁ + α₂ * J₂ + α₃ * J₃ + α₄ * J₄)

    ∂²An = (-ωD^2 + C^2) * An - 2.0 * ωD * C * Bn
    ∂²Bn = -2.0 * ωD * C * An + (ωD^2 - C^2) * Bn

    Pn = α₁ * sina_ + α₂ * cosa_ + α₃ * sinh_ + α₄ * cosh_

    return An, Bn, ∂²An, ∂²Bn, Pn
end
