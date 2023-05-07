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




jac_func="""
function jac!(Jac, u, p, t)
    # set everything to sero
    Jac .= 0.0

    # get the time index
    tval,idx_t = findmin(x -> abs.(x.- t), ts)

    # get the current ro2_ratio
    ro2_ratio = sum(u₀[idx_ro2])
    ro2_ratio = ro2_ratio/RO2ᵢ

    # set up product temporary value:
    prod_temp = 1.0

    @inbounds for i ∈ 1:length(jacobian_terms)
        update_jacobian!(
            idx_t,
            Jac,
            u,
            jacobian_terms[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end

    prod_temp = 1.0

    @inbounds for i ∈ 1:length(jacobian_terms_ro2)
        update_jacobian!(
            idx_t,
            Jac,
            u,
            jacobian_terms_ro2[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end
end
"""


function write_jac_func(;model_name::String="mcm")
    outpath = "./models/$(model_name)/jacobian.jl"

    # if it already exists, remove it so we can recreate it
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end

    open(outpath, "w") do f
        println(f, jac_func)
    end

end


function generate_jac_prototype(jacobian_terms::Vector{JacobianTerm}, jacobian_terms_ro2::Vector{JacobianTermRO2})
    idx_pairs = []
    for jac_term ∈ jacobian_terms
        push!(idx_pairs, (jac_term.i, jac_term.j))
    end
    for jac_term ∈ jacobian_terms_ro2
        push!(idx_pairs, (jac_term.i, jac_term.j))
    end
    idx_pairs = unique(idx_pairs)

    I = [idx_pair[1] for idx_pair ∈ idx_pairs]
    J = [idx_pair[2] for idx_pair ∈ idx_pairs]
    V = zeros(size(I))

    Jac = sparse(I,J,V)

    return Jac
end

