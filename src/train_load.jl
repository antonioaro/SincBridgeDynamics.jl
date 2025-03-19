struct TrainModel
    type::String
    axle::Vector{Float64}
    load::Vector{Float64}
    naxl::Int
    leng::Float64
    ν::Vector{Float64}
    DTS::Vector{Float64}    
    distributed::Bool
end

"""
    init_train_model(train; distributed=false) -> TrainModel

Initializes a train model for dynamic analysis.

# Parameters
- `train::Tuple`: A tuple containing:
- `type`: Train type identifier.
- `axle::Vector{Float64}`: Positions of the train's axles.
- `load::Vector{Float64}`: Load applied by each axle.
- `distributed::Bool=false`: If `true`, the train load is considered as a distributed load.

# Returns
A `TrainModel` object containing:
- `type`: Train type.
- `axle`: Axle positions.
- `load`: Axle loads.
- `naxl`: Number of axles.
- `leng`: Train length.
- `λ, G`: Train wavelength and Dyanamic Train Signature.
- `distributed`: Whether the train load is distributed.

"""
function init_train_model(train; distributed=false)
    type = train[1]
    axle = train[2]
    load = train[3]

    naxl = length(axle)
    leng = axle[end]
    ν, G = _DTS(axle, load)

    return TrainModel(type, axle, load, naxl, leng, ν, G, distributed)
end

function _DTS(train_axle, train_load)
    ν = range(3.0, step=0.01, stop=35)
    G = zeros(length(ν))
    ζ = 0.0

    naxl = length(train_axle)
    G_ = zeros(naxl)
    train_axle .-= train_axle[1]

    two_pi = 2π

    for (iν, νi) in pairs(ν)
        G_ .= 0.0

        for iaxl in 1:naxl
            A, B = 0.0, 0.0

            @inbounds for jaxl in 1:iaxl
                δ = train_axle[jaxl] / νi
                exp_term = exp(-two_pi * ζ * δ)
                cos_term, sin_term = cospi(2 * δ), sinpi(2 * δ)

                A += train_load[jaxl] * cos_term * exp_term
                B += train_load[jaxl] * sin_term * exp_term
            end

            G_[iaxl] = hypot(A, B)
        end

        G[iν] = maximum(G_)
    end

    return ν, G
end