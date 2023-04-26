using ChemicalDataAssimilation
using Plots
using Statistics
using DelimitedFiles, CSV, DataFrames
using ProgressMeter

mechpath = "mechanism-files/extracted/alkanes/methane.fac"
model_name = "methane"

@assert ispath(mechpath) == true

if !ispath("models/$model_name")
    mkpath("models/$model_name")
end


fac_dict = read_fac_file(mechpath)

fac_dict["generic_rate_coefficients"]
fac_dict["complex_rate_coefficients"]
fac_dict["peroxy_radicals"]
fac_dict["reaction_definitions"]

const Δt_step = 15.0  # time step in minutes

# 1. generate species indices, names, etc...
generate_species_df("models/names.csv", fac_dict; model_name=model_name)
df_species = CSV.File("models/$model_name/species.csv") |> DataFrame


# # 2. load in measurements
# lab_measurements = readdlm("models/measurement_names.csv", ',')
# idx_meas = Bool.(lab_measurements[2,:])
# mcm_name = lab_measurements[3, idx_meas]


# . generate lookup table for M, O2, N2, H2O
generate_densities("data/no_ap/number_densities.csv",
                   "data/w_ap/number_densities.csv";
                   model_name=model_name
                   )

df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame


# 2. Generate lookup table for Photolysis Rates
generate_photolysis_rates("data/no_ap/photolysis_rates.csv",
                          "data/w_ap/photolysis_rates.csv";
                          model_name=model_name,
                          Δt_step=Δt_step
                          )
df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame


# 3. Generate lookup table for generic/complex rate coefficients

# this will include *both* generic and complex
rate_list = generate_rrates_funcs(fac_dict; model_name=model_name)
include("./models/$(model_name)/rrates.jl")

generate_rrates(fac_dict,
                df_params,
                rate_list;
                model_name=model_name
                )

include("./models/$model_name/rrate_coeffs.jl")
df_rrate_coeffs = CSV.File("./models/$model_name/rrate_coeffs.csv") |> DataFrame



# generate list reaction, reaction rate coefficients
fac_dict["reaction_definitions"]



# now let's define reaction structs
# when we create the reactions, the reaction rate function should parsed to replace expression names with their dataframe
# counterparts, i.e.
#     KRO2HO2 ⟶ df_rrate_coeffs.KRO2HO2[idx_time]
#     T       ⟶ df_params.temperature
# where idx_time is fetched via the current time

# we need to check if any of the actual reaction rates used in the equations are repeated anywhere... ← it appears they don't!


# we should define a function, get_row_index(t::Float64) to return the index given an input time... NOTE: will this be differentiable?

# the reaction rate function in our structs should be like
# rrate_coeff(idx_time, df_params, df_rrate_coeffs, df_photolysis, etc...)


# we can then write a generic function to evaluate the reaction rate coefficient for the
# abstract reaction type, i.e.
# function reaction_rate_coeff(rxn::Reaction)




# 4. generate final lookup table for *all* reaction rate coefficients






