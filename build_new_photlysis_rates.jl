using CSV, DataFrames
using Plots
using LinearAlgebra
using Statistics
using ProgressMeter

using DataInterpolations
using QuadGK
using LaTeXStrings

mechpath = "mechanism-files/extracted/alkanes/methane.fac"
model_name = "methane"

# mechpath = "mechanism-files/extracted/full/mcm_subset.fac"
# model_name = "mcm_full"

@assert ispath(mechpath) == true

if !ispath("models/$model_name")
    mkpath("models/$model_name")
end


# load in number densities to get correct number of rows
# for output dataframe
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;

N = nrow(df_number_densities)


spec_path = "./data/spectra"
@assert ispath(spec_path) == true

# generate list of files
spec_files = []
for (root, dirs, files) ∈ walkdir(spec_path)
    for file ∈ files
        if endswith(file, ".txt")
            push!(spec_files, joinpath(root, file))
        end
    end
end

println("\n")
println(spec_files)
println("\n")


function read_oceanoptics(path)
    f = readlines(path)
    idx_wavelengths = findfirst(x -> x == ">>>>>Begin Spectral Data<<<<<", f) + 1 # find anything between } and ;

    f[idx_wavelengths:end]
    λs = [parse(Float64, split(line, "\t")[1]) for line ∈ f[idx_wavelengths:end]]
    Is = [parse(Float64, split(line, "\t")[2]) for line ∈ f[idx_wavelengths:end]]

    Is[Is .< 0.0] .= 0.0
    Is[Is .== -0.0] .= 0.0

    return Is, λs
end

# loop through all files and compute the mean spectrum
Is_list = []
λs_list = []

@showprogress for spec_f ∈ spec_files
    I, λ = read_oceanoptics(spec_f)
    push!(Is_list, I)
    push!(λs_list, λ)
end


Is = mean(Is_list)
λs = mean(λs_list)

# plot a sample spectrum
plot(λs, Is,
     xlabel="wavelength [nm]",
     ylabel="Irradiance [μW⋅cm⁻²⋅nm]",
     lw=2,
     label="",
     )


# let's convert to Photons/s/cm^2/nm using E_photon = hf = hc/λ

h_plank = 6.62607015e-34 * 1e6 # μJ Hz-1
c = 299792458.0e9 # nm / s

E_photon = (h_plank*c) ./ λs

Is

Is_photon = Is ./ E_photon

p = plot(λs, Is_photon,
         xlabel="wavelength [nm]",
         ylabel="Irradiance [Photons⋅s⁻¹cm⁻²⋅nm]",
         lw=2,
         label="",
         )

savefig("mean_irradiance.png")


λmin = minimum(λs)
λmax = maximum(λs)


# load in σ and Φ data
df_σs = CSV.File(download("https://raw.githubusercontent.com/john-waczak/Photolysis.jl/main/mcm/photolysis_%CF%83.csv")) |> DataFrame
df_Φs = CSV.File(download("https://raw.githubusercontent.com/john-waczak/Photolysis.jl/main/mcm/photolysis_%CE%A6.csv")) |> DataFrame

@assert nrow(df_σs) == length(Is)
@assert nrow(df_Φs) == length(Is)


# J = ∫dλ I(λ,T) σ(λ,T) Φ(λ,T)

# -----------------------------------------------------------------
# go through each of the photolysis rates
# -----------------------------------------------------------------

df_photolysis = DataFrame()

@showprogress for i ∈ 1:56
    try
        df_photolysis[:, "J_$i"] = get_J(i)[1] .* ones(N)
    catch e
        continue
    end
end

df_photolysis

# write the output to files
CSV.write("./data/photolysis_rates_corrected.csv", df_photolysis)


