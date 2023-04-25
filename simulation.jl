using ChemicalDataAssimilation
using Plots
using Statistics
using DelimitedFiles, CSV, DataFrames


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


Δt_step = 15.0  # time step in minutes

# 1. generate lookup table for M, O2, N2, H2O

df_number_densities_no_ap = CSV.File("data/no_ap/number_densities.csv") |> DataFrame
df_number_densities_no_ap.w_ap = [false for _ ∈ 1:nrow(df_number_densities_no_ap)]

df_number_densities_w_ap = CSV.File("data/w_ap/number_densities.csv") |> DataFrame
df_number_densities_w_ap.w_ap = [false for _ ∈ 1:nrow(df_number_densities_w_ap)]

df_number_densities = vcat(df_number_densities_no_ap[1:end-1,:], df_number_densities_w_ap)

# 2. Generate lookup table for Photolysis Rates

# join into one dataframe & add column with flag for no ap vs w ap
df_photolysis_no_ap = CSV.File("data/no_ap/photolysis_rates.csv") |> DataFrame
df_photolysis_no_ap.w_ap = [false for _ ∈ 1:nrow(df_photolysis_no_ap)]

df_photolysis_w_ap = CSV.File("data/w_ap/photolysis_rates.csv") |> DataFrame
df_photolysis_w_ap.w_ap = [false for _ ∈ 1:nrow(df_photolysis_w_ap)]

df_photolysis = vcat(df_photolysis_no_ap[1:end-1,:], df_photolysis_w_ap)
df_photolysis.t = df_photolysis.Δt .* Δt_step
select!(df_photolysis, Not(:Δt))

# get rid of the other dataframes
df_photolysis_no_ap = nothing
df_photolysis_w_ap = nothing
GC.gc()

# 3. Generate lookup table for generic rate coefficients

# this will include *both* generic and complex
df_general_rate_coeffs = DataFrame()
# copy in the times
df_general_rate_coeffs.t = df_photolysis.t
df_general_rate_coeffs.w_ap = df_photolysis.w_ap

# generate functions for generic and complex rrate coeffs
rate_list = generate_rrates(fac_dict; model_name=model_name)
include("./models/$(model_name)/rrates.jl")

for rate ∈ rate_list
    rrate_string = "$(rate).(df_number_densities.temperature, df_number_densities.pressure, df_number_densities.M, df_number_densities.O2, df_number_densities.N2,df_number_densities.H2O)"

    df_general_rate_coeffs[!, rate] =  eval(Meta.parse(rrate_string))
end

df_general_rate_coeffs

plot(df_general_rate_coeffs.t, df_general_rate_coeffs[!, "KRO2NO"])


# 4. generate lookup table for complex rate coefficients

# 5. generate final lookup table for *all* reaction rate coefficients

# 6. generate lookup table for MCM name → MCM index → SmilesString/IUPAC name

df_names = CSV.File("models/names.csv") |> DataFrame

chemical_formula = []
for row ∈ eachrow(df_names)
    if !ismissing(row["InChl Identifier"])
        if occursin("/", row["InChl Identifier"])
            push!(chemical_formula, split(row["InChl Identifier"], "/")[2])
        else
            push!(chemical_formula, missing)
        end
    else
        push!(chemical_formula, missing)
    end
end

df_names.chemical_formula = chemical_formula

df_names[1, "InChl Identifier"]

# load in measurements
lab_measurements = readdlm("models/measurement_names.csv", ',')
idx_meas = Bool.(lab_measurements[2,:])


mcm_name = lab_measurements[3, idx_meas]


# 7. Generate lookup table for all measured
