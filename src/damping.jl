function ec1damping(material::Vector{String}, L::Vector{Float64}, deck_type::String)

    nM = length(material)
    leng = maximum(L)

    ζ = zeros(Float64, nM)

    for iM in 1:nM
        if material[iM] in ["Composite steel-concrete", "Steel"]
            if leng <= 20.0
                ζ[iM] = 0.5 + 0.125 * (20.0 - leng)
            else
                ζ[iM] = 0.5
            end
        elseif material[iM] in ["Prestressed concrete"]
            if leng <= 20.0
                ζ[iM] = 1.0 + 0.07 * (20.0 - leng)
            else
                ζ[iM] = 1.0
            end
        elseif material[iM] in ["Reinforced concrete"] ||
               deck_type in ["Filler beam", "Filler beam reinforced"]
            if leng <= 20.0
                ζ[iM] = 1.5 + 0.07 * (20.0 - leng)
            else
                ζ[iM] = 1.5
            end
        end
    end

    return minimum(ζ) / 100.0
end