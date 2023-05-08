
# need to filter out any NaN's that appear


# assume that all measurements are of direct state variables, i.e. concentrations
function JObs(u::Vector{Float64}, idx_meas::Vector{Int})
    return I(length(idx_meas))
end
