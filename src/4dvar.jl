function Rmat(idx_t::Int, meas_ϵ::AbstractArray{Float64}; fudge_fac::Float64=1.0)
    return diagm((fudge_fac .* meas_ϵ[:,idx_t]) .^ 2)
end

function Rinv(idx_t::Int, meas_ϵ::AbstractArray{Float64}; fudge_fac::Float64=1.0)
    return diagm(1 ./ ((fudge_fac .* meas_ϵ[:,idx_t]) .^ 2))
end



# assume we already have constant defined with
# indices of the measurements
function ObsOpMeas(u, idx_meas)   # <--- observation operator for direct chemical measurements
    return u[idx_meas]
end



