ENV["GKSwstype"] = 100

println("Setting Up Julia Environment")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using ChemicalDataAssimilation
using DelimitedFiles, CSV, DataFrames
using Plots, StatsPlots
using Statistics
using ProgressMeter
using SparseArrays
using Measurements
using ArgParse


# # set random number seed
# rng = StableRNG(42)



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--model_name", "-n"
            help = "Name for the resulting model used in output paths"
            arg_type = String
            default = "methane"
        "--time_step", "-t"
            help = "The time step used during integration of mechanism (in minutes)."
            arg_type = Float64
            default = 15.0
        "--restart"
            help = "Whether or not to restart 4d_var from previous fitresult"
            action = :store_true
        "--fudge_fac", "-f"
            help = "A fudge factor for manipulating scale of measurement uncertainties"
            arg_type = Float64
            default = 0.5
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)
    @assert ispath("models/$(parsed_args[:model_name])/4dvar/u0.csv")  "Can not find models/$(parsed_args[:model_name])/4dvar/u0.csv"
    return parsed_args
end



# 0. parse arguments and set up output directory
println("Parsing command line arguments...")
parsed_args = parse_commandline()
model_name = parsed_args[:model_name]
Δt = parsed_args[:time_step]
fudge_fac = parsed_args[:fudge_fac]

# 1. load in dataframes
df_species = CSV.File("models/$model_name/species.csv") |> DataFrame;
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_params_ϵ = CSV.File("models/$model_name/state_parameters_ϵ.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;
df_number_densities_ϵ = CSV.File("models/$model_name/number_densities_ϵ.csv") |> DataFrame;
df_sol = CSV.File("models/$model_name/4dvar/u0_integrated.csv") |> DataFrame


# 2. generate measurement indices
idx_meas::Vector{Int} = Int[]
measurements_to_ignore = [:C2H6]  # skip any with nans
is_meas_in_mechanism = [spec ∈ df_species[!, "MCM Name"] for spec ∈ names(df_number_densities[:, Not([measurements_to_ignore..., :t, :w_ap])])]
idx_ts_to_use = df_number_densities.t .≤ 0.0
df_nd_to_use = df_number_densities[idx_ts_to_use, Not([measurements_to_ignore..., :t, :w_ap])]
df_nd_to_use_ϵ = df_number_densities_ϵ[idx_ts_to_use, Not([measurements_to_ignore..., :t, :w_ap])]

for spec_name ∈ names(df_nd_to_use)
    println(spec_name)
    push!(idx_meas, df_species[df_species[!, "MCM Name"] .== spec_name, :idx_species][1])
end

# 3. generate trajectory for measurements
u_sol = Matrix(df_sol)'
times = df_number_densities.t[df_number_densities.t .≤ 0]

spec_names = df_species[idx_meas, "MCM Name"]
u_meas = u_sol[idx_meas,:]


# 4. convert number density to mixing ratio
M = df_params[df_params.t .≤ 0.0, :M] .± (fudge_fac .* df_params_ϵ[df_params.t .≤ 0.0, :M])  # <-- total number density
u_mr = to_mixing_ratio(u_meas , Measurements.value.(M))
W =  Matrix(df_nd_to_use)' .±  (fudge_fac .* Matrix(df_nd_to_use_ϵ)')
W_mr = to_mixing_ratio(W, M)

# # convert to ppb
# u_ppb = u_mr .* 1e9
# W_ppb = W_mr .* 1e9

# u_mr[1,1]

# --------------------------------------------------------------------------------------------------------------------------
#  Generate plots
# --------------------------------------------------------------------------------------------------------------------------
println("Generating plots...")
@showprogress for i ∈ 1:length(idx_meas)
    plot_spec_name = df_species[idx_meas[i], "MCM Name"]

    units, unit_mult = get_reasonable_mr_units(u_mr[i,:])

    plot(times,
         u_mr[i,:] .* unit_mult,
         label="4dVar",
         title=plot_spec_name,
         lw=3
         )

    scatter!(times,
             Measurements.value.(W_mr[i,:]) .* unit_mult,
             yerror=Measurements.uncertainty.(W_mr[i,:]) .* unit_mult,
             label="Measurements",
             )
    xlabel!("time [minutes]")
    ylabel!("concentration [$(units)]")

    savefig("models/$model_name/4dvar/$(plot_spec_name).png")
    savefig("models/$model_name/4dvar/$(plot_spec_name).pdf")
end

println("...All done!")
