abstract type RxnJacobian end

struct JacobianTerm <: RxnJacobian
    # duᵢ/duⱼ
    i::Int
    j::Int
    idxs_in::Vector{Int}
    idx_k::Int
    prefac::Int
end

struct JacobianTermRO2 <: RxnJacobian
    i::Int
    j::Int
    idxs_in::Vector{Int}
    idx_k::Int
    prefac::Int
end



function JacobianTerms(drxn::DerivativeTerm)
    jac_terms = JacobianTerm[]

    i = drxn.idx_du
    # there will be a term for each i ∈
    for j ∈ drxn.idxs_in
        push!(jac_terms, JacobianTerm(
            i,
            j,
            [idx for idx ∈ drxn.idxs_in if idx != j],  # all the terms left after differentiation w.r.t uⱼ
            drxn.idx_k,
            drxn.prefac
        ))
    end

    return jac_terms
end

function JacobianTerms(drxn::DerivativeTermRO2)
    jac_terms = JacobianTermRO2[]

    i = drxn.idx_du
    # there will be a term for each i ∈
    for j ∈ drxn.idxs_in
        push!(jac_terms, JacobianTermRO2(
            i,
            j,
            [idx for idx ∈ drxn.idxs_in if idx != j],  # all the terms left after differentiation w.r.t uⱼ
            drxn.idx_k,
            drxn.prefac
        ))
    end

    return jac_terms
end



function update_jacobian!(idx_t::Int,
                          Jac,
                          u::Vector{Float64},
                          jac_term::JacobianTerm,
                          ro2_ratio::Float64,
                          K_matrix::Matrix{Float64},
                          Δt_step::Float64,
                          prod_temp::Float64
                          )
    for idx ∈ jac_term.idxs_in
        prod_temp *= u[idx]
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[idx_t, jac_term.idx_k] * prod_temp
end


function update_jacobian!(idx_t::Int,
                          Jac,
                          u::Vector{Float64},
                          jac_term::JacobianTermRO2,
                          ro2_ratio::Float64,
                          K_matrix::Matrix{Float64},
                          Δt_step::Float64,
                          prod_temp::Float64
                          )
    for idx ∈ jac_term.idxs_in
        prod_temp *= u[idx]
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[idx_t, jac_term.idx_k] * prod_temp * ro2_ratio
end


