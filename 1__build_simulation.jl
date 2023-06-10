ENV["GKSwstype"] = 100

println("Setting Up Julia Environment")
using Pkg
Pkg.activate(".")
Pkg.instantiate()
println("\t...Finished")

using ChemicalDataAssimilation
using Plots
using Statistics
using DelimitedFiles, CSV, DataFrames
using ProgressMeter
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
        "--use_updated_photolysis"
            help = "Whether or not to use updated photolysis"
            action = :store_true
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)

    @assert isfile(parsed_args[:mechanism_path]) "Supplied mechanism path does not exist"


    # make sure that the datapath and outpath exist
    if !ispath("models/$(parsed_args[:model_name])")
        println("$(parsed_args[:model_name]) directory does not exist in `./models`. Creating now...")
        mkpath("models/$(parsed_args[:model_name])")
    end


    @assert ispath("data/no_ap/number_densities.csv") "Can not find  data/no_ap/number_densities.csv"
    @assert ispath("data/w_ap/number_densities.csv") "Can not find  data/w_ap/number_densities.csv"
    @assert ispath("data/no_ap/number_densities_ϵ.csv") "Can not find  data/no_ap/number_densities_ϵ.csv"
    @assert ispath("data/w_ap/number_densities_ϵ.csv") "Can not find  data/w_ap/number_densities_ϵ.csv"


    return parsed_args
end


# parse arguments and set up output directory
parsed_args = parse_commandline()
mechpath = parsed_args[:mechanism_path]
model_name = parsed_args[:model_name]

# set time step
Δt_step = parsed_args[:time_step]  # time step in minutes



# generate dictionary with .fac file information
fac_dict = read_fac_file(mechpath)

# 1. generate species indices, names, etc...
generate_species_df("models/names.csv", fac_dict; model_name=model_name)
df_species = CSV.File("models/$model_name/species.csv") |> DataFrame


# 2. generate lookup table for M, O2, N2, H2O
generate_densities("data/no_ap/number_densities.csv",
                    "data/w_ap/number_densities.csv";
                    model_name=model_name
                    )
generate_densities("data/no_ap/number_densities_ϵ.csv",
                    "data/w_ap/number_densities_ϵ.csv";
                    model_name=model_name
                    )

df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame



# 3. Generate lookup table for Photolysis Rates
use_recalculated_photo_rates = parsed_args[:use_updated_photolysis]
if use_recalculated_photo_rates
    include("1b__build_new_photolysis_rates.jl")  # has the code to generate re-fitted photo rates
    df_photolysis = CSV.File("data/photolysis_rates_corrected.csv") |> DataFrame
else
    generate_photolysis_rates("data/no_ap/photolysis_rates.csv",
                              "data/w_ap/photolysis_rates.csv";
                              model_name=model_name,
                              Δt_step=Δt_step
                              )
    df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame
end


# 4. Generate lookup table for generic/complex rate coefficients

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


# 5. generate indices for ro2 sum
generate_ro2(fac_dict,
              model_name=model_name
              )
include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices



# 6. generate sane initial conditions
init_path = "./mechanism-files/initial_concentrations/full.txt"
@assert isfile(init_path) "Cant find ./mechanism-files/initial_concentrations/full.txt"


init_dict = generate_init_dict(init_path, df_params.M[1])
u₀ = zeros(Float64, nrow(df_species))

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


# 7. Generate reaction rate coefficients lookup table

# combine parameters into one long tuple
ro2_sum = sum(u₀[idx_ro2])
RO2ᵢ =ro2_sum > 0 ? ro2_sum : 1.0  # make sure we have at least "1 particle" or ratio trick will break


# now we should be able to generate dataframe with all reaction rate coefficients
generate_rrates_mechanism(fac_dict,
                          rate_list;
                          model_name=model_name
                          )

# load file and generate CSV
include("./models/$model_name/rrates_mechanism.jl")

# read in resulting CSV and fix any non-float data
df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end

CSV.write("./models/$model_name/rrate_coeffs_mech.csv", df_rrate_coeffs_mech)


# 8. generate stoichiometry matrix for later visualization
spec_list, N = generate_stoich_mat(
    fac_dict,
    model_name=model_name
);


# 9. generate rhs func
write_rhs_func(model_name=model_name)
include("models/$model_name/rhs.jl")


# 10. generate jacobian func
write_jac_func(model_name=model_name)
include("models/$model_name/jacobian.jl")





