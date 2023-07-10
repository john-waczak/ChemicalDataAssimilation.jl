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
        "--data_basepath"
            help = "Path to data files to be used for testing"
            arg_type = String
            default = "data/intertek-emergency-testing"
        "--collection_id"
            help = "Name of collection to analyze"
            arg_type = String
            default = "empty"
        "--unc_ext"
            help = "Extension for uncertainty files."
            arg_type = String
            default = "_std"
        "--mechanism_path"
            help = "Path to mechanism `.fac` file specifying the chemical mechanism to be used."
            arg_type = String
            default = "mechanism-files/extracted/alkanes/methane.fac"
        "--model_name"
            help = "Name for the resulting model used in output paths"
            arg_type = String
            default = "empty_methane"
        "--time_step"
            help = "The time step used during integration of mechanism (in minutes)."
            arg_type = Float64
            default = 15.0
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)

    @assert ispath(parsed_args[:data_basepath])
    @assert ispath(joinpath(parsed_args[:data_basepath], "number_densities", parsed_args[:collection_id]))
    @assert isfile(parsed_args[:mechanism_path]) "Supplied mechanism path does not exist"


    # make sure that the outpath exists
    if !ispath("models/$(parsed_args[:model_name])")
        println("$(parsed_args[:model_name]) directory does not exist in `./models`. Creating now...")
        mkpath("models/$(parsed_args[:model_name])")
    end

    return parsed_args
end


# parse arguments and set up output directory
parsed_args = parse_commandline()

data_basepath = parsed_args[:data_basepath]
model_name = parsed_args[:model_name]
mechanism_path = parsed_args[:mechanism_path]
collection_id = parsed_args[:collection_id]
unc_ext = parsed_args[:unc_ext]
Δt_step = parsed_args[:time_step]  # time step in minutes





# generate dictionary with .fac file information
fac_dict = read_fac_file(mechanism_path)

# 1. generate species indices, names, etc...
println("Species: ")
println(size(generate_species(fac_dict)))


names_path = joinpath(data_basepath, "names.csv")
@assert ispath(names_path)
generate_species_df(names_path, fac_dict; model_name=model_name)
df_species = CSV.File("models/$model_name/species.csv") |> DataFrame


# 2. generate lookup table for M, O2, N2, H2O
generate_densities(
    joinpath(data_basepath, "number_densities", collection_id, "number_densities.csv"),
    joinpath(data_basepath, "number_densities", collection_id, "number_densities"*unc_ext*".csv"),
    model_name=model_name
)


df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame



# 3. Generate lookup table for Photolysis Rates
generate_photolysis_rates(
    joinpath(data_basepath, "rates", collection_id, "photolysis_rates.csv"),
    model_name=model_name,
    Δt_step=Δt_step
)

df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame


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
# println("idx ro2: ", idx_ro2)
# println("type: ", typeof(idx_ro2), " eltype: ", eltype(idx_ro2))


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





