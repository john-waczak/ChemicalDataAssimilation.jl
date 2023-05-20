using ChemicalDataAssimilation
using Plots
using Statistics
using DelimitedFiles, CSV, DataFrames
using ProgressMeter

mechpath = "mechanism-files/extracted/alkanes/methane.fac"
model_name = "methane"

# mechpath = "mechanism-files/extracted/full/mcm_subset.fac"
# model_name = "mcm_full"


@assert ispath(mechpath) == true

if !ispath("models/$model_name")
    mkpath("models/$model_name")
end


fac_dict = read_fac_file(mechpath)

fac_dict["generic_rate_coefficients"]
fac_dict["complex_rate_coefficients"]
fac_dict["peroxy_radicals"]
fac_dict["reaction_definitions"]

const Δt_step::Float64 = 15.0  # time step in minutes

# 1. generate species indices, names, etc...
generate_species_df("models/names.csv", fac_dict; model_name=model_name)
df_species = CSV.File("models/$model_name/species.csv") |> DataFrame



# . generate lookup table for M, O2, N2, H2O
generate_densities("data/no_ap/number_densities.csv",
                   "data/w_ap/number_densities.csv";
                   model_name=model_name
                   )

df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame


# 2. Generate lookup table for Photolysis Rates
use_recalculated_photo_rates = true
if use_recalculated_photo_rates
    df_photolysis = CSV.File("data/photolysis_rates_corrected.csv") |> DataFrame
else
    generate_photolysis_rates("data/no_ap/photolysis_rates.csv",
                              "data/w_ap/photolysis_rates.csv";
                              model_name=model_name,
                              Δt_step=Δt_step
                              )
    df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame
end


# photolysis_fudge_fac = 1e6
# for colname ∈ names(df_photolysis)
#     df_photolysis[!,colname] .= df_photolysis[!,colname] .* photolysis_fudge_fac
# end

# df_photolysis


describe(df_photolysis)


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


# generate indices for ro2 sum
generate_ro2(fac_dict,
             model_name=model_name
             )
include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices



# generate sane initial conditions
init_path = "./mechanism-files/initial_concentrations/full.txt"
isfile(init_path)

init_dict = generate_init_dict(init_path, df_params.M[1])
u₀    = zeros(Float64, nrow(df_species))

names(df_species)

for (key, val) ∈ init_dict
    try
        println("$(key): $(val)")
        idx = df_species[df_species[!, "MCM Name"] .== key, :idx_species][1]
        #idx = findfirst(x -> x == key, species)
        println("idx: ", idx)
        u₀[idx] = val
    catch e
        println("\t$key not in mechanism")
    end
end

u₀



# combine parameters into one long tuple
df_species[idx_ro2,:]
ro2_sum = sum(u₀[idx_ro2])
const RO2ᵢ =ro2_sum > 0 ? ro2_sum : 1.0  # make sure we have at least "1 particle"


# now we should be able to generate dataframe with all reaction rate coefficients
generate_rrates_mechanism(fac_dict,
                          rate_list;
                          model_name=model_name
                          )

include("./models/$model_name/rrates_mechanism.jl")
df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end


# now let's create a file holding all of the reaction structs in a vector that we can include

fac_dict["reaction_definitions"][1]
species, reactions = parse_rxns(fac_dict["reaction_definitions"])
size(reactions)


rxns = generate_reaction_list(fac_dict, df_species)


# generate stoichiometry matrix for later visualization
spec_list, N = generate_stoich_mat(
    fac_dict,
    model_name=model_name
)
