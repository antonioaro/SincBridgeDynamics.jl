struct Bridge
    E::Float64
    I::Vector{Float64}
    Ne::Int
    Le::Vector{Float64}
    μ::Float64
    ζ::Float64
    Lb::Float64
    bc::Vector{String}
end

"""
    init_bridge_model(E, I, Le, μ, bc, ζ=0.0) -> Bridge

Initializes a bridge model with the given properties.

# Parameters
- `E::Float64`: Young's modulus.
- `I::Union{Float64, Vector{Float64}}`: Moment of inertia, either as a single value (applied to all elements) or as a vector specifying values per element.
- `Le::Union{Float64, Vector{Float64}}`: Element lengths, either as a single value (for uniform elements) or as a vector.
- `μ::Float64`: Linear mass density of the bridge.
- `bc::Vector{String}`: Boundary conditions: pinned, fixed, free, roller, joint, hinge.
- `ζ::Float64=0.0`: Structural damping ratio (default: `0.0`).

# Returns
A `Bridge` object containing the discretized representation of the bridge.
"""
function init_bridge_model(E, I, Le, μ, bc, ζ=0.0)

    Le = Le isa AbstractVector ? Le : [Le]

    # number of elements    
    Ne = length(Le)

    # bridge length
    Lb = sum(Le)

    # inertia
    I = I isa AbstractVector ? I : fill(I, Ne)

    return Bridge(E, I, Ne, Le, μ, ζ, Lb, bc)
end
