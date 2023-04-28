using ChemicalDataAssimilation
using Plots
using Statistics
using DelimitedFiles, CSV, DataFrames
using ProgressMeter

# mechpath = "mechanism-files/extracted/alkanes/methane.fac"
# model_name = "methane"

mechpath = "mechanism-files/extracted/full/mcm_subset.fac"
model_name = "mcm_full"


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



# now let's create a file holding all of the reaction structs in a vector that we can include

fac_dict["reaction_definitions"][1]
species, reactions = parse_rxns(fac_dict["reaction_definitions"])
size(reactions)

for i ∈ 1:length(reactions)
    reactants, reactants_stoich, products, products_stoich, rrate_string = reactions[i]
    if nothing ∈ products
        println(i, "\t", products)
    end
end


rxns = ChemicalDataAssimilation.Reaction[]
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



# generate stoichiometry matrix for later visualization



# now let's define reaction structs
# when we create the reactions, the reaction rate function should parsed to replace expression names with their dataframe
# counterparts, i.e.
#     KRO2HO2 ⟶ df_rrate_coeffs.KRO2HO2[idx_time]
#     T       ⟶ df_params.temperature
# where idx_time is fetched via the current time

# we need to check if any of the actual reaction rates used in the equations are repeated anywhere... ← it appears they don't!


# we should define a function, get_row_index(t::Float64) to return the index given an input time... NOTE: will this be differentiable?

# the reaction rate function in our structs should be like
# rrate_coeff(idx_time, df_params, df_rrate_coeffs, df_photolysis, ro2_ratio, etc...)


# we can then write a generic function to evaluate the reaction rate coefficient for the
# abstract reaction type, i.e.
# function reaction_rate_coeff(rxn::Reaction)




# 4. generate final lookup table for *all* reaction rate coefficients



# generate file with constant vector of reactions
# i.e.
# const rxns = Reaction[...]

species, reactions = parse_rxns(fac_dict["reaction_definitions"])
unique_species = [spec for spec ∈ species if spec != nothing]

for i ∈ 1:length(reactions)
    reactants, reactants_stoich, products, products_stoich, rrate_string = reactions[i]
    if occursin("J<", rxn_rate_string) || occursin("J <", rxn_rate_string)
        #generate photodissociation
    else
        # generate collision
    end
end
