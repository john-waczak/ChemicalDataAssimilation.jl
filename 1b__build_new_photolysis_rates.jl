# This file is to be included within 1__build_simulation.jl

using CSV, DataFrames
using Plots
using LinearAlgebra
using Statistics
using ProgressMeter
using DataInterpolations
using QuadGK
using LaTeXStrings



N = nrow(df_number_densities)

spec_path = "./data/spectra"
@assert (ispath(spec_path) == true) "Can't find spectra in ./data/spectra"

if !ispath("models/photolysis")
    mkpath("models/photolysis")
end


# generate list of files
spec_files = []
for (root, dirs, files) ∈ walkdir(spec_path)
    for file ∈ files
        if endswith(file, ".txt")
            push!(spec_files, joinpath(root, file))
        end
    end
end


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

savefig("models/photolysis/sample_spec.png")
savefig("models/photolysis/sample_spec.pdf")

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

savefig("models/photolysis/mean_irradiance.png")
savefig("models/photolysis/mean_irradiance.pdf")


# compute wavelength boundaries
λmin = minimum(λs)
λmax = maximum(λs)


# load in σ and Φ data
df_σs = CSV.File(download("https://raw.githubusercontent.com/john-waczak/Photolysis.jl/main/mcm/photolysis_%CF%83.csv")) |> DataFrame
df_Φs = CSV.File(download("https://raw.githubusercontent.com/john-waczak/Photolysis.jl/main/mcm/photolysis_%CE%A6.csv")) |> DataFrame

# save cross section and quantum yield data for later viewing:
CSV.write("models/photolysis/cross_sections.csv", df_σs)
CSV.write("models/photolysis/quantum_yields.csv", df_Φs)

@assert nrow(df_σs) == length(Is)
@assert nrow(df_Φs) == length(Is)


# J = ∫dλ I(λ,T) σ(λ,T) Φ(λ,T)

# -----------------------------------------------------------------
# go through each of the photolysis rates
# -----------------------------------------------------------------


function get_J(idx_reaction)
    integrand = ConstantInterpolation(Is_photon .* df_σs[:, "σ_$idx_reaction"] .* df_Φs[:, "Φ_$idx_reaction"], λs)
    J, J_error = quadgk(integrand, λmin, λmax)
end


plot(
    λs,
    log10.(Is_photon),
    ylabel="log10 Irradiance or σ",
    lw=3,
    color=:royalblue,
    label="I(λ,T)",
    legend=:inside,
)

plot!(
    λs,
    log10.(df_σs[:, :σ_1]),
    lw=3,
    color=:orange,
    #ylabel="log10 σ",
    label="σ(λ,T)",
)


plot!(
    twinx(),
    λs,
    df_Φs[:,:Φ_1],
    lw=3,
    color=:gray,
    ylabel="Quantum Yield",
    label="Φ(λ,T)",
    legend=:right
)

title!("J₁ = $(get_J(1)[1])")
savefig("models/photolysis/J1_photolysis_summary.png")
savefig("models/photolysis/J1_photolysis_summary.pdf")



df_photolysis = DataFrame()

@showprogress for i ∈ 1:56
    try
        df_photolysis[:, "J_$i"] = get_J(i)[1] .* ones(N)
    catch e
        # not all numbers in this range have corresponding photolysis rates. Just ignore the missing ones
        continue
    end
end


# write the output to files
CSV.write("./data/photolysis_rates_corrected.csv", df_photolysis)


