"""
    init_time_analysis(FREQB, FREQE, Lb, Vt, Lt, Nλ) -> (Δt, NSTEPS)

Initializes the time discretization parameters.

# Parameters
- `FREQB::Float64`: Lower bound of the frequency range.
- `FREQE::Float64`: Upper bound of the frequency range.
- `Lb::Float64`: Total length of the bridge.
- `Vt::Float64`: Train velocity.
- `Lt::Float64`: Train length.
- `Nλ::Int`: Number of modes shapes.

# Returns
- `Δt::Float64`: Time step for numerical integration.
- `NSTEPS::Int`: TNumber of time steps.
"""
function init_time_analysis(FREQB, FREQE, Lb, Vt, Lt, Nλ)
    nst1 = 20.0
    nst2 = 1.0

    TB = 1 / FREQB
    TE = 1 / FREQE

    Δt = minimum([TE / nst1, Lb / (Vt * Nλ * nst2)])
    tSpan = (Lb + Lt) / Vt + 10.0 * TB
    NSTEPS = Int(cld(tSpan , Δt))

    return Δt, NSTEPS
end

function energy_percentile(x, Δt)
    f, X= rffts(x, Δt)
    X = cumsum(abs.(X) .^ 2)
    X ./= maximum(X)

    return f, X
end

function rffts(x, Δt)
    X = rfft(x) * Δt
    nfft = length(X) * 2 - 1
    f = rfftfreq(nfft, 1 / Δt)

    return f, X
end