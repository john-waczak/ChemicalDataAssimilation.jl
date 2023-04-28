struct DerivativeTerm
    idx_du::Int           # index of derivative term du[i] to update
    idxs_in::Vector{Int}  # indices of reactants that appear in rate
    idx_k::Int            # reaction rate coefficient index
    prefac::Int           # ± 1
    needs_ro2::Bool       # whether or not we need the ro2 ratio
end


function DerivativeTerms(rxn::CollisionReaction)
    dts = DerivativeTerm[]

    # loop over reactants
    for idx_du ∈ rxn.idxs_in
        push!(dts, DerivativeTerm(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            -1,
            rxn.needs_ro2
        ))
    end

    # loop over products
    for idx_du ∈ rxn.idxs_out
        push!(dts, DerivativeTerm(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            1,
            rxn.needs_ro2
        ))
    end

    return dts
end


function DerivativeTerms(rxn::Photodissociation)
    dts = DerivativeTerm[]

    # loop over reactants
    for idx_du ∈ rxn.idxs_in
        push!(dts, DerivativeTerm(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            -1,
            false
        ))
    end

    # loop over products
    for idx_du ∈ rxn.idxs_out
        push!(dts, DerivativeTerm(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            1,
            false
        ))
    end

    return dts

end



function get_time_index(t::Float64, Δt_step::Float64)
    return round(Int, t/Δt_step)
end


function get_step_time(t::Float64, Δt_step::Float64)
    return round(Int, t/Δt_step) * Δt_step
end


function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::DerivativeTerm,
                            ro2_ratio::Float64,
                            df_rrate_coeffs_mech::DataFrame,
                            Δt_step::Float64
                            )

    du[deriv_term.idx_du] += deriv_term.prefac * df_rrate_coeffs_mech[idx_t, "k_$(deriv_term.idx_k)"] * prod(u[deriv_term.idxs_in])


    # if we need to, multiply by ratio to use actual RO2 value
    if deriv_term.needs_ro2
        du[deriv_term.idx_du] *= ro2_ratio
    end
end


function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::DerivativeTerm,
                            ro2_ratio::Float64,
                            K_matrix::Matrix{Float64},
                            Δt_step::Float64
                            )

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[idx_t, deriv_term.idx_k] * prod(u[deriv_term.idxs_in])


    # if we need to, multiply by ratio to use actual RO2 value
    if deriv_term.needs_ro2
        du[deriv_term.idx_du] *= ro2_ratio
    end
end
