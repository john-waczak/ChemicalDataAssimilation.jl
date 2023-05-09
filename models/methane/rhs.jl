function rhs!(du, u, p, t)
    # set everything to sero
    du .= 0.0

    # get the time index
    tval,idx_t = findmin(x -> abs.(x.- t), ts)

    # get the current ro2_ratio
    ro2_ratio = sum(u[idx_ro2])
    ro2_ratio = ro2_ratio/RO2ᵢ

    # set up product temporary value:
    prod_temp = 1.0

    # update derivatives
    @inbounds for i ∈ 1:length(derivatives)
        prod_temp = 1.0
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end

    prod_temp = 1.0

    @inbounds for i ∈ 1:length(derivatives_ro2)
        prod_temp = 1.0
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_ro2[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end
end

