abstract type RxnDerivative end

struct DerivativeTerm <: RxnDerivative
    idx_du::Int           # index of derivative term du[i] to update
    idxs_in::Vector{Int}  # indices of reactants that appear in rate
    idx_k::Int            # reaction rate coefficient index
    prefac::Int           # ± 1
end

struct DerivativeTermRO2 <: RxnDerivative
    idx_du::Int           # index of derivative term du[i] to update
    idxs_in::Vector{Int}  # indices of reactants that appear in rate
    idx_k::Int            # reaction rate coefficient index
    prefac::Int           # ± 1
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
        ))
    end

    # loop over products
    for idx_du ∈ rxn.idxs_out
        push!(dts, DerivativeTerm(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            1,
        ))
    end

    return dts
end


function DerivativeTerms(rxn::CollisionReactionRO2)
    dts = DerivativeTermRO2[]

    # loop over reactants
    for idx_du ∈ rxn.idxs_in
        push!(dts, DerivativeTermRO2(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            -1,
        ))
    end

    # loop over products
    for idx_du ∈ rxn.idxs_out
        push!(dts, DerivativeTermRO2(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            1,
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
        ))
    end

    # loop over products
    for idx_du ∈ rxn.idxs_out
        push!(dts, DerivativeTerm(
            idx_du,
            rxn.idxs_in,
            rxn.idx_k,
            1,
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
                            K_matrix::Matrix{Float64},
                            Δt_step::Float64,
                            prod_temp::Float64
                            )

    # prod_temp = u[deriv_term.idxs_in[1]]
    # for i ∈ 2:size(deriv_term.idxs_in,1)
    for i ∈ axes(deriv_term.idxs_in,1)
        prod_temp *= u[deriv_term.idxs_in[i]]
    end


    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[idx_t, deriv_term.idx_k] * prod_temp

    # du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[idx_t, deriv_term.idx_k] * prod(u[deriv_term.idxs_in])
end

function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::DerivativeTermRO2,
                            ro2_ratio::Float64,
                            K_matrix::Matrix{Float64},
                            Δt_step::Float64,
                            prod_temp::Float64
                            )

    # prod_temp = u[deriv_term.idxs_in[1]]
    # for i ∈ 2:size(deriv_term.idxs_in,1)
    for i ∈ axes(deriv_term.idxs_in,1)
        prod_temp *= u[deriv_term.idxs_in[i]]
    end

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[idx_t, deriv_term.idx_k] * prod_temp * ro2_ratio

    # du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[idx_t, deriv_term.idx_k] * prod(u[deriv_term.idxs_in]) * ro2_ratio
end



rhs_func = """
function rhs!(du, u, p, t)
    # set everything to sero
    du .= 0.0

    # get the time index
    tval,idx_t = findmin(x -> abs.(x.- t), ts)

    # get the current ro2_ratio
    ro2_ratio = sum(u[idx_ro2])
    ro2_ratio = ro2_ratio/RO2ᵢ
    ro2_ratio = ro2_ratio > 0.0 ? ro2_ratio : 1.0  # make sure we don't have issues with ratio of 0

    # set up product temporary value:
    prod_temp = 1.0

    # update derivatives
    @inbounds for i ∈ 1:length(derivatives)
        prod_temp = 1.0  # <-- need to start fresh for each derivative
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
        prod_temp = 1.0  # <-- need to start fresh for each derivative
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
"""

function write_rhs_func(;model_name::String="mcm")
    outpath = "./models/$(model_name)/rhs.jl"

    # if it already exists, remove it so we can recreate it
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end

    open(outpath, "w") do f
        println(f, rhs_func)
    end

end

