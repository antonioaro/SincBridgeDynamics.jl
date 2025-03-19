struct Modal{N}
    λ::Matrix{N}
    v::Matrix{N}
    ω::Vector{N}
    x::Vector{N}
    Φ::Vector{Vector{N}}
    Mn::Vector{N}
end

function Modal(; λ::Matrix{N}, v::Matrix{N}, ω::Vector{N}, x::Vector{N}, Φ::Vector{Vector{N}}, Mn::Vector{N}) where {N}
    return Modal{N}(λ, ω, v, x, Φ, Mn)
end

"""
    modal(bridge, FREQB, FREQE; neigs=Inf, maxit=100) -> Modal

Performs modal analysis of a multi-span beam within a frequency range defined by `FREQB` (lower frequency) and `FREQE` (upper frequency).

# Parameters
- `bridge`: Object representing the bridge structure.
- `FREQB::Float64`: Lower bound of the frequency range.
- `FREQE::Float64`: Upper bound of the frequency range.
- `neigs::Int=Inf`: Maximum number of modes to compute. By default, all modes within the frequency range are computed.
- `maxit::Int=100`: Maximum number of iterations for solving the eigenvalue problem.

# Returns
A `Modal` object containing:
- `λ`: Eigenvalues.
- `v`: Eigenvectors.
- `ω`: Natural frequencies.
- `x`: Spatial coordinate.
- `Φ`: Mass-normalized mode shapes.
- `Mn`: Modal mass.

"""
function modal(bridge, FREQB, FREQE; neigs=Inf, maxit=100)

    Nx = 101#31

    x = Vector{Float64}(undef, Nx * bridge.Ne)
    x₀ = 0.0
    for i in 1:bridge.Ne
        x₁ = x₀ + bridge.Le[i]
        x[Nx*(i-1)+1:Nx*i] .= LinRange(x₀, x₁, Nx)
        x₀ = x₁
    end

    # nonlinear eigenvalue problem (λ,v)
    λ, v, ω = _modaleigs(bridge.E, bridge.I, bridge.Le, bridge.μ, bridge.bc, FREQB, FREQE, neigs, maxit)
    nλ = size(λ, 1)

    # modal mass and eigenvectors
    Mn = Vector{Float64}(undef, nλ)
    Φ = Vector{Vector{Float64}}(undef, nλ)

    @inbounds for iλ = 1:nλ
        # mode shape
        @views Φ[iλ] = _modeigen(λ[iλ, :], v[:, iλ], bridge.Le, x)
        # modal mass
        Mn[iλ] = _modalmass(λ[iλ, :], v[:, iλ], bridge.Le, bridge.μ)
        # mass normalisation
        Φ[iλ] .*= 1 / √(Mn[iλ])
    end

    return Modal(λ, v, ω, x, Φ, Mn)
end

function _modaleigs(E, I, L, μ, bc, FREQB, FREQE, neigs, maxit)

    # number of elements
    ne = length(L)

    # normalised parameter
    K̃ = @. E * I / μ / L^4
    K̃₀, k₀ = findmin(K̃)

    Λ = @. √(√(K̃₀ / K̃))
    Γ = @. E * I / (E * I[k₀])

    # frequency range
    ωB = 2.0 * π * FREQB
    λB = √(ωB / √(K̃₀))

    ωE = 2.0 * π * FREQE
    λE = √(ωE / √(K̃₀))

    # number of eigenvector
    nv = 4 * ne

    # vector of matrices, function an system matrix
    AA = Vector{Matrix{Float64}}()
    fii = Vector{Function}()
    T = spzeros(nv, nv)

    # loop over elements
    irow = 0

    @inbounds for ie = 1:ne
        bc_type1 = boundary_conditions[bc[ie]]
        bc_type2 = boundary_conditions[bc[ie+1]]

        jcol = 4 * (ie - 1)
        irow = boundary_condition(bc_type1, AA, fii, T, irow, jcol, Λ[ie], L[ie], Γ[ie], 0.0)
        irow = boundary_condition(bc_type2, AA, fii, T, irow, jcol, Λ[ie], L[ie], Γ[ie], 1.0)
    end

    # group unique fii
    λ_ = π / 6
    val = [func(λ_) for func in fii]

    # ia = unique(i -> round(val[i], sigdigits=3), eachindex(val))
    ia = unique(i -> val[i], eachindex(val))
    fi = fii[ia]

    ic = [findall(x -> x == v, val) for v in val[ia]]

    A = Vector{Matrix{Float64}}()
    for k in eachindex(ic)
        push!(A, sum(AA[ic[k]]))
    end

    # minum eigenvalue, center and scaling factor
    center = λB + (λE - λB) / 2
    scaling = center - λE

    # solve  (set neigs=inf to find all eigenvalue with maxit)
    nep = SPMF_NEP(A, fi)
    (λ, v) = iar(nep, σ=center, γ=scaling, v=ones(nv), tol=1e-9, neigs=neigs, maxit=maxit)

    # find real(λ) > 0 
    i = findall(real(λ) .> λB)
    λ = real(λ[i])
    v = real(v[:, i])

    # sort
    p = sortperm(λ)
    λ = λ[p]
    v = v[:, p]

    # element eigenvalue
    λₑ = λ * Λ'

    # natural frequency
    ω = λ .^ 2 * √(K̃₀)

    return λₑ, v, ω
end

function _modeigen(λ, v, L, x)
    Ne = length(L)
    Ψ = zeros(Float64, length(x))
    x₀ = 0.0

    @inbounds @simd for i in 1:Ne
        x₁ = x₀ + L[i]

        α1 = v[(i-1)*4+1]
        α2 = v[(i-1)*4+2]
        α3 = v[(i-1)*4+3]
        α4 = v[(i-1)*4+4]

        @inbounds for j in eachindex(x)
            if x[j] ≥ x₀ && x[j] ≤ x₁
                xe = λ[i] / L[i] * (x[j] - x₀)
                s, c = SLEEF.sincos_fast(xe)
                Ψ[j] = α1 * s + α2 * c + α3 * sinh(xe) + α4 * cosh(xe)
            end
        end
        x₀ = x₁
    end

    return Ψ
end


function _modalmass(λ, v, L, μ)
    Mn = 0.0
    for i in eachindex(L)
        l = (i - 1) * 4
        Me = _element_mass(L[i], λ[i], v[l+1:l+4], μ)
        Mn += Me
    end

    return Mn
end

function _element_mass(L, λ, v, μ)

    s, c = SLEEF.sincos_fast(λ)
    sh = sinh(λ)
    ch = cosh(λ)

    α1, α2, α3, α4 = v

    Mn = 2.0 * λ * (α1^2 + α2^2 - α3^2 + α4^2) +
         4.0 * (α1 * α2 * (1.0 - c^2) - α3 * α4 * (1.0 - ch^2)) +
         4.0 * ch * ((α2 * α3 - α1 * α4) * c + (α1 * α3 + α2 * α4) * s) +
         2.0 * (α2^2 - α1^2) * s * c +
         2.0 * sh * (2.0 * (-α1 * α3 + α2 * α4) * c + (α3^2 + α4^2) * ch + 2.0 * (α1 * α4 + α2 * α3) * s)

    Mn *= L * μ / (4.0 * λ)

    return Mn
end