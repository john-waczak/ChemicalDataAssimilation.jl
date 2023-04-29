using ChemicalDataAssimilation
using Plots
using Statistics
using DelimitedFiles, CSV, DataFrames
using ProgressMeter
using BenchmarkTools

using Statistics
using DifferentialEquations
using Sundials

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

df_species = CSV.File("models/$model_name/species.csv") |> DataFrame
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame
df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame
df_rrate_coeffs = CSV.File("./models/$model_name/rrate_coeffs.csv") |> DataFrame

# grab ro2 indices
include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices


# generate initial conditions
init_path = "./mechanism-files/initial_concentrations/full.txt"
isfile(init_path)

init_dict = generate_init_dict(init_path, df_params.M[1])
u₀    = zeros(Float64, nrow(df_species))

names(df_species)

for (key, val) ∈ init_dict
    try
        println("$(key): $(val)")
        idx = df_species[df_species[!, "MCM Name"] .== key, :idx_species][1]

        println("\tidx: ", idx)
        u₀[idx] = val
    catch e
        println("\t$key not in mechanism")
    end
end

# need to update this to use measurements

# set up initial RO2 value using u₀
ro2_sum = sum(u₀[idx_ro2])
const RO2ᵢ =ro2_sum > 0 ? ro2_sum : 1.0  # make sure we have at least "1 particle"

# now we can load in the reaction rate coefficients
df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end


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


# convert reaction rates into matrix for faster lookup
const K_matrix = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts = df_rrate_coeffs_mech.t

# define values used for update
idx_t = 0
tval = -180.0
ro2_ratio = 1.0

function rhs!(du, u, p, t)
    # set everything to sero
    du .= 0.0

    # get the time index
    tval,idx_t = findmin(x -> abs.(x.- t), ts)

    # get the current ro2_ratio
    ro2_ratio = sum(u₀[idx_ro2])
    ro2_ratio = ro2_ratio/RO2ᵢ

    # update derivatives
    @inbounds for i ∈ 1:length(derivatives)
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives[i],
            ro2_ratio,
            K_matrix,
            Δt_step
        )
    end

    @inbounds for i ∈ 1:length(derivatives_ro2)
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_ro2[i],
            ro2_ratio,
            K_matrix,
            Δt_step
        )
    end
end


# this should preallocate the rhs function
du = copy(u₀)
@benchmark rhs!(du, u₀,nothing, 1.0)


tmin = minimum(ts)
tmax = maximum(ts)

tspan = (tmin, tmax)

ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(rhs!, u₀, tspan)


tol = 1e-6
@benchmark solve(ode_prob,
                 CVODE_BDF();
                 saveat=15.0,
                 reltol=tol,
                 abstol=tol
                 )


