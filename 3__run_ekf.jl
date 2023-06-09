ENV["GKSwstype"] = 100

println("Setting Up Julia Environment")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# investigate possible uses of StaticArrays.jl, specifically the SizedArray{} decorator
using ChemicalDataAssimilation
using DelimitedFiles, CSV, DataFrames  # file reading
using Plots, StatsPlots # visualization
using Statistics # for mean, std, and var functions
using DifferentialEquations
using Sundials  # For CVODE_BDF()
using Zygote, ForwardDiff
using SciMLSensitivity
using ProgressMeter
using BenchmarkTools
#using Optimization
#using OptimizationOptimJL  # for ADAM
using LinearAlgebra
using Measurements
using SparseArrays
using ArgParse



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--mechanism_path"
            help = "Path to mechanism `.fac` file specifying the chemical mechanism to be used."
            arg_type = String
            default = "mechanism-files/extracted/alkanes/methane.fac"
        "--model_name"
            help = "Name for the resulting model used in output paths"
            arg_type = String
            default = "methane"
        "--time_step"
            help = "The time step used during integration of mechanism (in minutes)."
            arg_type = Float64
            default = 15.0
        "--try_solve"
            help = "Whether or not to precompile solvers by calling once."
            action = :store_true
        "--fudge_fac"
            help = "A fudge factor for manipulating scale of measurement uncertainties"
            arg_type = Float64
            default = 0.5  # 1.0
        "--epsilon"
            help = "Estimated background uncertainty for diagonal of B matrix, i.e. uncertainty in initial condition"
            arg_type = Float64
            default = 1.0  # 1.0
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)

    @assert isfile(parsed_args[:mechanism_path]) "Supplied mechanism path does not exist"


    # make sure that the datapath and outpath exist
    if !ispath("models/$(parsed_args[:model_name])")
        println("$(parsed_args[:model_name]) directory does not exist in `./models`. Creating now...")
        mkpath(parsed_args[:datapath])
    end


    @assert ispath("data/no_ap/number_densities.csv") "Can not find  data/no_ap/number_densities.csv"
    @assert ispath("data/w_ap/number_densities.csv") "Can not find  data/w_ap/number_densities.csv"
    @assert ispath("data/no_ap/number_densities_ϵ.csv") "Can not find  data/no_ap/number_densities_ϵ.csv"
    @assert ispath("data/w_ap/number_densities_ϵ.csv") "Can not find  data/w_ap/number_densities_ϵ.csv"
    @assert ispath("models/$(parsed_args[:model_name])/4dvar/u0.csv")  "Can not find models/$(parsed_args[:model_name])/4dvar/u0.csv"



    return parsed_args
end



# 0. parse arguments and set up output directory
println("Parsing command line arguments...")
parsed_args = parse_commandline()
mechpath = parsed_args[:mechanism_path]
model_name = parsed_args[:model_name]


if !isdir("models/$model_name/EKF")
    mkpath("models/$model_name/EKF")
end




# --------------------------------------------------------------------------------------------------------------------------
# 1. read in dataframes
# -------------------------------------------------------------------------------------------------------------------------
println("Loading Data into DataFrames...")

fac_dict = read_fac_file(mechpath)
df_species = CSV.File("models/$model_name/species.csv") |> DataFrame;
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_params_ϵ = CSV.File("models/$model_name/state_parameters_ϵ.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;
df_number_densities_ϵ = CSV.File("models/$model_name/number_densities_ϵ.csv") |> DataFrame;

# --------------------------------------------------------------------------------------------------------------------------
# 2. Generate Initial Conditions
# --------------------------------------------------------------------------------------------------------------------------

df_u₀ = CSV.File("models/$model_name/4dvar/u0.csv") |> DataFrame
u₀ = df_u₀.u0
@assert typeof(u₀) == Vector{Float64}


# --------------------------------------------------------------------------------------------------------------------------
# 3. Get RO₂ indices and initial value
# --------------------------------------------------------------------------------------------------------------------------
include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices
ro2_sum = sum(u₀[idx_ro2])


# --------------------------------------------------------------------------------------------------------------------------
# 4. Get reaction rate coefficients
# --------------------------------------------------------------------------------------------------------------------------

df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame;

# make sure every entry is a Float64
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end



# --------------------------------------------------------------------------------------------------------------------------
# 5. Get measurement indices
# --------------------------------------------------------------------------------------------------------------------------

const idx_meas::Vector{Int} = Int[]

is_meas_in_mechanism = [spec ∈ df_species[!, "MCM Name"] for spec ∈ names(df_number_densities[:, Not([:t, :w_ap])])]

df_nd_to_use = df_number_densities[:, Not([:t, :w_ap])]
df_nd_to_use = df_nd_to_use[:, is_meas_in_mechanism]

df_nd_to_use_ϵ = df_number_densities_ϵ[:, Not([:t, :w_ap])]
df_nd_to_use_ϵ = df_nd_to_use_ϵ[:, is_meas_in_mechanism]

for spec_name ∈ names(df_nd_to_use)
    println(spec_name)
    push!(idx_meas, df_species[df_species[!, "MCM Name"] .== spec_name, :idx_species][1])
end


# --------------------------------------------------------------------------------------------------------------------------
# 6. Setup all other constants
# --------------------------------------------------------------------------------------------------------------------------

const Δt_step::Float64 = parsed_args[:time_step]  # time step in minutes
const rxns::Array{ChemicalDataAssimilation.Reaction} = ChemicalDataAssimilation.Reaction[]
const derivatives::Vector{ChemicalDataAssimilation.DerivativeTerm} = ChemicalDataAssimilation.DerivativeTerm[]
const derivatives_ro2::Vector{ChemicalDataAssimilation.DerivativeTermRO2} = ChemicalDataAssimilation.DerivativeTermRO2[]
const jacobian_terms::Vector{ChemicalDataAssimilation.JacobianTerm} = ChemicalDataAssimilation.JacobianTerm[]
const jacobian_terms_ro2::Vector{ChemicalDataAssimilation.JacobianTermRO2} = ChemicalDataAssimilation.JacobianTermRO2[]
const RO2ᵢ::Float64 =ro2_sum > 0.0 ? ro2_sum : 1.0  # make sure we have at least "1 particle" in RO2 sum to prevent division issues
const K_matrix::Matrix{Float64} = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts::Vector{Float64} = df_rrate_coeffs_mech.t
#const fudge_fac::Float64 = 0.5  # for measurement uncertainty
const fudge_fac::Float64 = parsed_args[:fudge_fac]  # for measurement uncertainty
const meas_ϵ::Matrix{Float64} = Matrix(df_nd_to_use_ϵ)'
const W::Matrix{Float64} = Matrix(df_nd_to_use)'
const P::Matrix{Float64} = zeros(nrow(df_species), nrow(df_species))
const P_diag::Matrix{Float64} = zeros(nrow(df_species), length(ts)) # i.e. P_diag[i] == P[i,i]
const Q::Matrix{Float64} = zeros(size(P))
const uₐ::Matrix{Float64} = zeros(length(u₀), length(ts))


const tmin::Float64 = minimum(ts)
const tmax::Float64 = maximum(ts)
const tol::Float64 = 1e-3
const tspan = (tmin, tmax)

const ϵ::Float64 = parsed_args[:epsilon]
const ϵ_min::Float64 = 1e-12


const try_solve::Bool = parsed_args[:try_solve]

# --------------------------------------------------------------------------------------------------------------------------
# 7. Generate reactions, derivatives, jacobians
# --------------------------------------------------------------------------------------------------------------------------

# create vector of reaction objects
species, reactions = parse_rxns(fac_dict["reaction_definitions"])

@showprogress for i ∈ 1:length(reactions)
    try
        push!(rxns, parse_rxn(reactions[i], i, df_species))
    catch e
        println(reactions[i])
        println("\n")
        println(e)
        break
    end
end

println("num reactions: ", length(rxns))
@assert length(rxns) == length(reactions)

# create derivative list
@showprogress for rxn ∈ rxns
    dts = DerivativeTerms(rxn)
    if eltype(dts) <: ChemicalDataAssimilation.DerivativeTerm
        for dt ∈ dts
            push!(derivatives, dt)
        end
    else
        for dt ∈ dts
            push!(derivatives_ro2, dt)
        end
    end
end

println("num derivative terms: ", length(derivatives) + length(derivatives_ro2))


# create jacobian list
@showprogress for drxn ∈ derivatives
    j_terms = JacobianTerms(drxn)
    for j_term ∈ j_terms
        push!(jacobian_terms, j_term)
    end
end

@showprogress for drxn ∈ derivatives_ro2
    j_terms = JacobianTerms(drxn)
    for j_term ∈ j_terms
        push!(jacobian_terms_ro2, j_term)
    end
end


println("num jacobian terms: ", size(jacobian_terms,1) + size(jacobian_terms_ro2,1))



# --------------------------------------------------------------------------------------------------------------------------
# 8. Get rhs and jacobian functions, generate jacobian prototype sparse matrix
# --------------------------------------------------------------------------------------------------------------------------
include("models/$model_name/rhs.jl")
include("models/$model_name/jacobian.jl")

# evaluate each once to precompile
du = copy(u₀)
rhs!(du, u₀, nothing, -180.0)

jac_prototype = generate_jac_prototype(jacobian_terms, jacobian_terms_ro2)

println("jacobian non-zero percentage: $(length(nonzeros(jac_prototype))/length(jac_prototype)*100)%")

jac!(jac_prototype, u₀, nothing, 1.0)



# --------------------------------------------------------------------------------------------------------------------------
# 9. Set up ODE Defaults
# --------------------------------------------------------------------------------------------------------------------------

fun = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)

if try_solve
  sol = solve(ode_prob, CVODE_BDF(); saveat=15.0, reltol=tol, abstol=tol);
end

# test that fetching R matrix works
is_meas_not_nan = get_idxs_not_nans(meas_ϵ[:,end])
@benchmark Rmat_nonan(1, is_meas_not_nan, meas_ϵ; fudge_fac=fudge_fac)

# test observation operator
u_h = Obs(u₀, idx_meas[is_meas_not_nan])



# --------------------------------------------------------------------------------------------------------------------------
# 10. Initialize Covariance Matrices
# --------------------------------------------------------------------------------------------------------------------------

for i ∈ 1:length(u₀)
    P[i,i] = (ϵ * u₀[i])^2 + ϵ_min^2
    P_diag[i,1] = P[i,i]
end

# initially, set Q to match P
Q .= P


# set the first value to the background estimate
uₐ[:,1] .= u₀


function model_forward(u_now, t_now)
    _prob = remake(ode_prob, u0=u_now, tspan=(t_now, t_now+Δt_step))
    solve(_prob, CVODE_BDF(), reltol=tol, abstol=tol, dense=false,save_everystep=false,save_start=false, sensealg=QuadratureAdjoint())[:,end]
end


u_test, DM_test = Zygote.withjacobian(model_forward, u₀, ts[1])

# --------------------------------------------------------------------------------------------------------------------------
# 11. Assimilation Loop
# --------------------------------------------------------------------------------------------------------------------------

const u_now::Vector{Float64} = zeros(size(uₐ[:,1]))
const DM::Matrix{Float64} = DM_test[1]


@showprogress for k ∈ 1:length(ts)-1  # because we always update the *next* value
    # k=2

    # collect current model estimate
    u_now .= uₐ[:,k]  # should preallocate this

    # --------------------------
    # Forecast Step
    # --------------------------

    # run model forward one Δt
    u_next, DM_tup = Zygote.withjacobian(model_forward, u_now, ts[k])  # <-- can't make this mutating

    DM .= DM_tup[1]  # result is a 1-element tuple, so we index it

    # collect observations
    local is_meas_not_nan = get_idxs_not_nans(W[:,k+1])
    local idx_meas_nonan = idx_meas[is_meas_not_nan]
    local u_h = Obs(u_next, idx_meas_nonan)


    # update loop for the mode covariance matrix, Q
    if k > 1
        Q .= 0.0
        for j ∈ axes(Q,1)
            # define error growth rate
            eg = abs(u_next[j]-u_now[j])/max(1.0, u_next[j])  # note: we should actually do the sens. analysis here

            # clip the growth rate to 1%-50% range
            eg = max(min(eg, 0.5), 0.01)


            # if we have a measurement, further manipulate the growth rate
            if j ∈ idx_meas_nonan
                u_h_j = u_h[idx_meas_nonan .== j][1]
                diff = abs(u_next[j] - u_h_j)
                ratio = diff/u_h_j

                if ratio > 1
                    eg = max(eg, 0.5)
                elseif ratio > 0.5
                    eg = max(eg, 0.2)
                elseif ratio > 0.2
                    eg = max(eg, 0.1)
                else
                    continue
                end
            end

            # finally, update Q
            Q[j,j] = eg * u_next[j]^2 + ϵ_min^2
        end
    end

    # update the background covariance matrix
    P .= DM*P*DM' + Q

    # --------------------------
    # Analysis Step
    # --------------------------

    DH = JObs(u_h, u_next, idx_meas_nonan)
    R = Rmat_nonan(k+1, is_meas_not_nan, meas_ϵ; fudge_fac=fudge_fac)
    denom = DH*P*DH' + R

    #Kalman = P*DH'/denom
    Kalman = P*DH'* inv(denom)  # <-- this seems slightly faster for now

    u_next .= u_next .+ Kalman*(W[is_meas_not_nan, k+1] - u_h)

    P .= (I(length(u_next)) - Kalman*DH)*P


    # do some clipping to make sure we stay reasonable
    for j ∈ axes(P,1)
        P[j,j] = min(max(0, u_next[j]^2), P[j,j])
        P[j,j] = max(0.05*u_next[j]^2, P[j,j])
    end


    # filter negative values to zero
    u_next[u_next .≤ 0.0] .= 0.0

    # update the analysis vector
    update_uₐ!(uₐ, u_next, k+1)
    P_diag[:,k+1] .= [P[i,i] for i ∈ axes(P,1)]
end


# --------------------------------------------------------------------------------------------------------------------------
# 12. Post Process
# --------------------------------------------------------------------------------------------------------------------------


# combine output w/ uncertainty from diagonal
uₐ_nd = uₐ .± sqrt.(P_diag)

# convert final output into mixing ratios
M = df_params.M .± (fudge_fac .* df_params_ϵ.M)

uₐ_mr = to_mixing_ratio(uₐ_nd, M)

# chop off values and uncertainties for easier plotting
ua_mr_vals = Measurements.value.(uₐ_mr)
ua_mr_ϵ = Measurements.uncertainty.(uₐ_mr)


# save to output file
df_ekf = DataFrame()
df_ekf_ϵ = DataFrame()

@showprogress for i ∈ axes(ua_mr_vals, 1)
    df_ekf[!, df_species[i, "MCM Name"]] = ua_mr_vals[i,:]
    df_ekf_ϵ[!, df_species[i, "MCM Name"]] = ua_mr_ϵ[i,:]
end

df_ekf[!, :times] = ts
df_ekf_ϵ[!, :times] = ts

CSV.write("models/$model_name/EKF/ekf_output.csv", df_ekf)
CSV.write("models/$model_name/EKF/ekf_ϵ_output.csv", df_ekf_ϵ)


# combine measurements with uncertainties
W_mr = W .± (fudge_fac .* meas_ϵ)
W_mr = to_mixing_ratio(W_mr, M)

W_mr_val = Measurements.value.(W_mr)
W_mr_ϵ = Measurements.uncertainty.(W_mr)

# save measurements to csv files for final output
size(W_mr_val)
idx_meas

df_w = DataFrame()
df_w_ϵ = DataFrame()
@showprogress for i ∈ axes(W_mr_val, 1)
    df_w[!, df_species[idx_meas[i], "MCM Name"]] = W_mr_val[i,:]
    df_w_ϵ[!, df_species[idx_meas[i], "MCM Name"]] = W_mr_ϵ[i,:]
end

CSV.write("models/$model_name/EKF/ekf_measurements.csv", df_w)
CSV.write("models/$model_name/EKF/ekf_measurements_ϵ.csv", df_w_ϵ)



# --------------------------------------------------------------------------------------------------------------------------
# 13. Plots
# --------------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------------
# 14. Compute Lifetime
# --------------------------------------------------------------------------------------------------------------------------

# we write the lifetime for species i as

# τᵢ = ∑ⱼτᵢⱼ where j are all reactions for which i is a reactant

# the form for each τᵢⱼ depends on the specific reaction type

# Photodissociation reaction:
# X + hν ⟶ products
# Ẋ = -jX  [molecules/cm³/s]
# τ = 1\j  [s]

# Collision reaction:
# X + ∑ⱼYⱼ ⟶ products
# Ẋ = -kX⋅ΠⱼYⱼ   [molecules/cm³/s]
# τ = 1\(kΠⱼYⱼ)

# Collision reaction w/ RO2 (i.e. all termolecular reactions) will look the same as above since M is included already inside of our computed k.

derivatives
derivatives_ro2

# now we need to combine

ua_nd =  Measurements.value.(uₐ_nd)
size(ua_nd)

τs = copy(ua_nd)  # preallocate matrix to hold values
ℓ_mat = zeros(size(ua_nd))  # loss rate
#ℓ = 1.0

@showprogress for d ∈ 1:length(derivatives)
    derivative = derivatives[d]
    if derivative.prefac < 0.0 # i.e. if it's negative so that we have a reactant not product
        for idx_t ∈ axes(K_matrix,1)
            ℓ = K_matrix[idx_t,derivative.idx_k]
            for i ∈ derivative.idxs_in
                ℓ  *= ua_nd[i, idx_t]
            end

            ℓ_mat[derivative.idx_du, idx_t] += ℓ
        end
    end
end

@showprogress for d ∈ 1:length(derivatives_ro2)
    derivative = derivatives_ro2[d]
    if derivative.prefac < 0.0 # i.e. if it's negative so that we have a reactant not product
        for idx_t ∈ axes(K_matrix,1)
            ℓ = K_matrix[idx_t,derivative.idx_k]
            for i ∈ derivative.idxs_in
                ℓ  *= ua_nd[i, idx_t]
            end

            ℓ_mat[derivative.idx_du, idx_t] += ℓ
        end
    end
end


for j ∈ axes(τs, 2), i ∈ axes(τs,1)
    if isinf(τs[i,j]/ℓ_mat[i,j] ) || isnan(τs[i,j]/ℓ_mat[i,j] )
        println("idx: ", (i,j), "\tu:\t", τs[i,j], "\tℓ:\t",ℓ_mat[i,j], "\tτ:\t", τs[i,j]/ℓ_mat[i,j] )
    end

    τs[i,j] = τs[i,j] / ℓ_mat[i,j]
end



# generate lifetime dataframes for output
df_τs = DataFrame()
@showprogress for i ∈ axes(τs, 1)
    df_τs[!, df_species[i, "MCM Name"]] = τs[i,:]
end

df_τs.t = df_params.t

CSV.write("models/$model_name/EKF/lifetimes.csv", df_τs)


df_τs_means = DataFrame()
@showprogress for i ∈ axes(τs, 1)
    df_τs_means[!, df_species[i, "MCM Name"]] = [mean(τs[i,:])]
end

df_τs_means

τs_means = Matrix(df_τs_means)
idx_sort = sortperm(τs_means, dims=2, rev=true)

spec_name = []
mean_lifetime = []

for idx ∈ idx_sort
    push!(spec_name, names(df_τs_means)[idx])
    push!(mean_lifetime,df_τs_means[1, idx] )
end

df_τs_means_sorted = DataFrame(:species => spec_name, :τ_seconds => mean_lifetime)

df_τs_means_sorted.τ_minutes = (df_τs_means_sorted.τ_seconds ./ 60)
df_τs_means_sorted.τ_hours = (df_τs_means_sorted.τ_minutes ./ 60)
df_τs_means_sorted.τ_days = (df_τs_means_sorted.τ_hours ./ 24)
df_τs_means_sorted.τ_weeks = (df_τs_means_sorted.τ_days ./7)
df_τs_means_sorted.τ_years = (df_τs_means_sorted.τ_days ./365)


CSV.write("models/$model_name/EKF/mean_lifetimes.csv", df_τs_means_sorted[7:end,:])
