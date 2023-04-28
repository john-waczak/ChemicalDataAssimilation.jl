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


rxns = generate_reaction_list(fac_dict, df_species)

# generate stoichiometry matrix for later visualization
spec_list, N = generate_stoich_mat(
    fac_dict,
    model_name=model_name
)


# generate derivative list

rxns[1].idxs_in
rxns[1].idxs_out
rxns[1].idx_k
rxns[1].needs_ro2

derivatives = ChemicalDataAssimilation.DerivativeTerm[]
@showprogress for rxn ∈ rxns
    dts = DerivativeTerms(rxn)
    for dt ∈ dts
        push!(derivatives, dt)
    end
end

length(derivatives)


prod(u₀[rxns[1].idxs_in])

# generate ODE RHS function


# we should define a function, get_row_index(t::Float64) to return the index given an input time... NOTE: will this be differentiable?


du = copy(u₀)

using BenchmarkTools

@benchmark update_derivative!(1,
                   du,
                   u₀,
                   derivatives[1],
                   10.0,
                   df_rrate_coeffs_mech,
                   Δt_step,
                   )



@benchmark 1.0 * @view df_rrate_coeffs_mech[1, "k_1"]



K_matrix = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
typeof(K_matrix)
@benchmark 1.0 * K_matrix[1,1]


# so it's clearly worth it to use a matrix


ttest = 1.1

function rhs1(t)
    _, idx_t = findmin(x -> abs.(x.- ttest), df_rrate_coeffs_mech.t )
    ro2_ratio = sum(u₀[idx_ro2])/RO2ᵢ
    for derivative ∈ derivatives
        update_derivative!(
            idx_t,
            du,
            u₀,
            derivative,
            ro2_ratio,
            df_rrate_coeffs_mech,
            Δt_step
        )
    end
end


function rhs2(t)
    _, idx_t = findmin(x -> abs.(x.- ttest), df_rrate_coeffs_mech.t )
    ro2_ratio = sum(u₀[idx_ro2])/RO2ᵢ
    for derivative ∈ derivatives
        update_derivative!(
            idx_t,
            du,
            u₀,
            derivative,
            ro2_ratio,
            K_matrix,
            Δt_step
        )
    end
end

@benchmark rhs1(1.0)
@benchmark rhs2(1.0)




df_rrate_coeffs_mech


step_time = round(ttest/Δt_step) * Δt_step
df_rrate_coeffs_mech.t


findmin(x->abs(x-t), a)
