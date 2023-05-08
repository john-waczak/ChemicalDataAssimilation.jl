function Rmat_nonan(idx_t::Int, idx_meas::BitVector, meas_ϵ::AbstractArray{Float64}; fudge_fac::Float64=1.0)
    return diagm((fudge_fac .* meas_ϵ[idx_meas,idx_t]) .^ 2)
end

function Rinv_nonan(idx_t::Int, idx_meas::BitVector, meas_ϵ::AbstractArray{Float64}; fudge_fac::Float64=1.0)
    return diagm(1 ./ ((fudge_fac .* meas_ϵ[idx_meas,idx_t]) .^ 2))
end


function get_idxs_not_nans(meas_ϵ_vec)
#    return [!isnan(meas_ϵ_vec[idx]) for idx ∈ 1:length(meas_ϵ_vec)]
    return .!(isnan.(meas_ϵ_vec))
end




# need to filter out any NaN's that appear
function Obs(u::Vector{Float64}, idx_meas::Vector{Int})
    return u[idx_meas]
end


# assume that all measurements are of direct state variables, i.e. concentrations
function JObs(u_h::Vector{Float64}, uₐ::Vector{Float64}, idx_meas::Vector{Int})
    res = zeros(length(u_h), length(uₐ))
    for i ∈ 1:length(idx_meas)
        res[i, idx_meas[i]] = 1
    end

    return res
end



function getTimeIndex(t, ts)
    return argmin(abs.(ts .- t))
end


# define helper function to update the analysis result
# function update_uₐ!(uₐ, ua_new, t_now, ts)
#     idx_t = getTimeIndex(t_now, ts)
#     uₐ[:,idx_t] .= ua_new
# end

function update_uₐ!(uₐ, ua_new, idx_t)
    uₐ[:,idx_t] .= ua_new
end




# tell Zygote to ignore our update function
Zygote.@nograd update_uₐ!

