using Statistics

function to_mixing_ratio(Xnd,  Mnd)
    return Xnd ./ Mnd
end


function to_mixing_ratio(Xnd::AbstractMatrix, Mnd::AbstractVector)
    @assert size(Xnd, 2) == size(Mnd,1)
    Xmr = copy(Xnd)

    for j ∈ axes(Xnd, 2)
        Xmr[:,j] .= Xmr[:,j] ./ Mnd[j]
    end

    return Xmr
end


function get_reasonable_mr_units(X)
    mean_val = mean(X)

    units = ["%","ppm", "ppb", "ppt", "ppq"]
    scale = [1e-2, 1e-6, 1e-9, 1e-12, 1e-15]
    multiplier = [1e2, 1e6, 1e9, 1e12, 1e15]

    # compute abs difference between mean and scales
    idx_min = argmin(abs.(mean_val .- scale))
    min_diff = minimum(abs.(mean_val .- scale))

    if min_diff < 0.1 && min_diff > 100
        return "mixing ratio", 1.0
    else
        return units[idx_min], multiplier[idx_min]
    end
end
