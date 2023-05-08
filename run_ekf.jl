using ChemicalDataAssimilation

# file reading
using DelimitedFiles, CSV, DataFrames

# visualization
using Plots, StatsPlots

# for mean, std, and var functions
using Statistics

# for ODE problem and integration
using DifferentialEquations
using Sundials

# for autodiff of ODE problem solutions
using Zygote, ForwardDiff
using SciMLSensitivity

# convenience packages
using ProgressMeter
using BenchmarkTools

# for optimization of 4dvar cost function
using Optimization
using OptimizationOptimJL

# for generating specialized matrix types
using LinearAlgebra


# --------------------------------------------------------------------------------------------------------------------------
# Setup paths
# --------------------------------------------------------------------------------------------------------------------------

mechpath = "mechanism-files/extracted/alkanes/methane.fac"
model_name = "methane"

# mechpath = "mechanism-files/extracted/full/mcm_subset.fac"
# model_name = "mcm_full"

@assert ispath(mechpath) == true

if !ispath("models/$model_name")
    mkpath("models/$model_name")
end

fac_dict = read_fac_file(mechpath)

const Δt_step = 15.0  # time step in minutes

df_species = CSV.File("models/$model_name/species.csv") |> DataFrame;
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;
df_number_densities_ϵ = CSV.File("models/$model_name/number_densities_ϵ.csv") |> DataFrame;

# df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame
# df_rrate_coeffs = CSV.File("./models/$model_name/rrate_coeffs.csv") |> DataFrame;


# --------------------------------------------------------------------------------------------------------------------------
# Generate Initial Conditions
# --------------------------------------------------------------------------------------------------------------------------

df_u₀ = CSV.File("models/$model_name/4dvar/u0.csv") |> DataFrame
u₀ = df_u₀.u₀
@assert typeof(u₀) == Vector{Float64}


# --------------------------------------------------------------------------------------------------------------------------
# Get RO₂ indices and initial value
# --------------------------------------------------------------------------------------------------------------------------
include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices

ro2_sum = sum(u₀[idx_ro2])
const RO2ᵢ =ro2_sum > 0 ? ro2_sum : 1.0  # make sure we have at least "1 particle"



# --------------------------------------------------------------------------------------------------------------------------
# Get reaction rate coefficients
# --------------------------------------------------------------------------------------------------------------------------

df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame;
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end


# --------------------------------------------------------------------------------------------------------------------------
# Get measurement indices
# --------------------------------------------------------------------------------------------------------------------------

is_meas_in_mechanism = [spec ∈ df_species[!, "MCM Name"] for spec ∈ names(df_number_densities[:, Not([:t, :w_ap])])]

df_nd_to_use = df_number_densities[:, Not([:t, :w_ap])]
df_nd_to_use = df_nd_to_use[:, is_meas_in_mechanism]

df_nd_to_use_ϵ = df_number_densities[:, Not([:t, :w_ap])]
df_nd_to_use_ϵ = df_nd_to_use_ϵ[:, is_meas_in_mechanism]

const idx_meas = Int[]
for spec_name ∈ names(df_nd_to_use)
    println(spec_name)
    push!(idx_meas, df_species[df_species[!, "MCM Name"] .== spec_name, :idx_species][1])
end

idx_meas


# --------------------------------------------------------------------------------------------------------------------------
# Generate reactions, derivatives, jacobians
# --------------------------------------------------------------------------------------------------------------------------

# create vector of reaction objects
species, reactions = parse_rxns(fac_dict["reaction_definitions"])
const rxns = ChemicalDataAssimilation.Reaction[]
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
const derivatives = ChemicalDataAssimilation.DerivativeTerm[]
const derivatives_ro2 = ChemicalDataAssimilation.DerivativeTermRO2[]
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
const jacobian_terms = ChemicalDataAssimilation.JacobianTerm[]
const jacobian_terms_ro2 = ChemicalDataAssimilation.JacobianTermRO2[]
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
# Generate lookup table for reaction rate coefficients and time values
# --------------------------------------------------------------------------------------------------------------------------

const K_matrix = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts = df_rrate_coeffs_mech.t

@assert size(K_matrix, 1) == length(ts)
@assert size(K_matrix, 2) == length(rxns)


# --------------------------------------------------------------------------------------------------------------------------
# Get rhs and jacobian functions, generate jacobian prototype sparse matrix
# --------------------------------------------------------------------------------------------------------------------------
include("models/$model_name/rhs.jl")
include("models/$model_name/jacobian.jl")

# evaluate each once to precompile
du = copy(u₀)
rhs!(du, u₀, nothing, -180.0)

jac_prototype = generate_jac_prototype(jacobian_terms, jacobian_terms_ro2)
# jac_prototype=zeros(size(u₀,1),size(u₀,1))

jac!(jac_prototype, u₀, nothing, 1.0)




# --------------------------------------------------------------------------------------------------------------------------
#  Set up ODE Defaults
# --------------------------------------------------------------------------------------------------------------------------

tmin = minimum(ts)
tmax = maximum(ts)
#tmax = 0.0  # only go until t=0 for 4dVar
tspan = (tmin, tmax)
tol = 1e-3

#ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(rhs!, u₀, tspan)
#fun = ODEFunction(rhs!; jac=jac!, jac_prototype=Jac)
fun = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)


sol = solve(
    ode_prob,
    CVODE_BDF();
    saveat=15.0,
    reltol=tol,
    abstol=tol,
);


# define R
meas_ϵ = collect(df_number_densities_ϵ[1, Not([measurements_to_ignore..., :t, :w_ap])])
fudge_fac = 0.1
#fudge_fac = 1.0
const R = diagm( (fudge_fac .* meas_ϵ) .^ 2)
const Rinv = inv(R)



# define W matrix of observations


## verify that each time series includes NaN's and that observation operator ignores NaN values...

# define Jacobian of Observation Operator

# define Q

const Q = 0.0 * I(nrow(df_species))  # we're not including process noise for now

# define covariance matrix and set it to some initial values
const P = Symmetric(zeros(nrow(df_species), nrow(df_species)))

fudge_prefac = 0.1
for i ∈ 1:length(u₀)
    P[i,i] = (fudge_prefac * u₀[i])^2
end


# pre-allocate vector to contain diagonal entries (we will think about full matrix serialization later)
const P_diag = zeros(nrow(df_species), length(ts)) # i.e. P_diag[i] == P[i,i]
P_diag[:,1] .= [P[i,i] for i ∈ 1:length(u₀)]



# preallocate kalman gain matrix
const Kₖ = zeros(nrow(df_species), length(idx_meas))

# A×X=B     --> X = A\B
# A == X×B  --> X = A/B  (we'll use this one I think)


# preallocate reaction rate coefficients and times
const K_matrix = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts = df_rrate_coeffs_mech.t


# include rhs and jacobian
include("models/$model_name/rhs.jl")
include("models/$model_name/jacobian.jl")


# define ODE Problem
fun = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)

# integrator = init(
#     ode_prob,
#     CVODE_BDF();
#     # saveat=15.0,
#     dense=false,
# #    save_everystep=false,
# #    save_start=false,
# #    save_end=false,
#     reltol=tol,
#     abstol=tol
# )


# # step!(integrator,
# #       Δt_step,
# #       true
# #       )

# # integrator.sol

# initialize the analysis matrix
uₐ = zeros(length(u₀), length(ts))

# set the first value to the background estimate
uₐ[:,1] .= u₀


function getTimeIndex(t, Δt)
    return Int(round(t/Δt)) + 1  # since arrays are 1-indexed
end

# define helper function to update the analysis result
function update_uₐ!(uₐ, ua_new, t_now, Δt)
    idx_t = getTimeIndex(t_now, Δt)
    uₐ[:,idx_t] .= ua_new
end

# tell Zygote to ignore our update function
Zygote.@nograd update_uₐ!


function model_forward!(u_now, t_now)
    _prob = remake(ode_prob, u0=u_now, tspan=(t_now, t_now+Δt_step))
    solve(_prob, CVODE_BDF(), reltol=tol, abstol=tol, dense=false,save_everystep=false,save_start=false, sensealg=QuadratureAdjoint())[:,end]
end


@benchmark Zygote.withjacobian(model_forward!, u₀, ts[1])


@benchmark solve(ode_prob, CVODE_BDF(linear_solver=:Dense), saveat=15.0, reltol=tol, abstol=tol)

# @benchmark solve(ode_prob, ImplicitMidpoint(), dt=Δt_step, saveat=15.0, reltol=tol, abstol=tol, verbose=false)
# @benchmark solve(ode_prob, FBDF(), saveat=15.0, reltol=tol, abstol=tol, verbose=false)
# @benchmark solve(ode_prob, CVODE_BDF(linear_solver=:LapackDense), saveat=15.0, reltol=tol, abstol=tol)
# @benchmark solve(ode_prob, CVODE_BDF(linear_solver=:KLU), saveat=15.0, reltol=tol, abstol=tol)


idx_meas
