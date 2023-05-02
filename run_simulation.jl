using ChemicalDataAssimilation
using Plots
using Statistics
using DelimitedFiles, CSV, DataFrames
using ProgressMeter
using BenchmarkTools

using Statistics
using DifferentialEquations
using Sundials

mechpath = "mechanism-files/extracted/alkanes/methane.fac"
model_name = "methane"

# mechpath = "mechanism-files/extracted/full/mcm_subset.fac"
# model_name = "mcm_full"

@assert ispath(mechpath) == true

if !ispath("models/$model_name")
    mkpath("models/$model_name")
end

fac_dict = read_fac_file(mechpath)

const Δt_step = 15.0  # time step in minutes

df_species = CSV.File("models/$model_name/species.csv") |> DataFrame
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame
df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame
df_rrate_coeffs = CSV.File("./models/$model_name/rrate_coeffs.csv") |> DataFrame



# grab ro2 indices
include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices


# generate initial conditions
init_path = "./mechanism-files/initial_concentrations/full.txt"
isfile(init_path)

init_dict = generate_init_dict(init_path, df_params.M[1])
u₀    = zeros(Float64, nrow(df_species))

names(df_species)

for (key, val) ∈ init_dict
    try
        println("$(key): $(val)")
        idx = df_species[df_species[!, "MCM Name"] .== key, :idx_species][1]

        println("\tidx: ", idx)
        u₀[idx] = val
    catch e
        println("\t$key not in mechanism")
    end
end

# need to update this to use measurements
df_nd_init = df_number_densities[1, Not([:C2H6, :t])]

for name ∈ names(df_nd_init)
    if name ∈ df_species[:, "MCM Name"]
        idx = df_species[df_species[:, "MCM Name"] .== name, :].idx_species[1]
        println("Old val: ", u₀[idx])
        u₀[idx] = df_nd_init[name]
        println("New val: ", u₀[idx])
    end
end



# set up initial RO2 value using u₀
ro2_sum = sum(u₀[idx_ro2])
const RO2ᵢ =ro2_sum > 0 ? ro2_sum : 1.0  # make sure we have at least "1 particle"

# now we can load in the reaction rate coefficients
df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end


# create vector of reaction objects
species, reactions = parse_rxns(fac_dict["reaction_definitions"])
const rxns = ChemicalDataAssimilation.Reaction[]
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


# create derivative list
const derivatives = ChemicalDataAssimilation.DerivativeTerm[]
const derivatives_ro2 = ChemicalDataAssimilation.DerivativeTermRO2[]
@showprogress for rxn ∈ rxns
    dts = DerivativeTerms(rxn)
    if eltype(dts) <: ChemicalDataAssimilation.DerivativeTerm
        for dt ∈ dts
            push!(derivatives, dt)
        end
    else
        for dt ∈ dts
            push!(derivatives_ro2, dt)
        end
    end
end


# create jacobian list
const jacobian_terms = ChemicalDataAssimilation.JacobianTerm[]
const jacobian_terms_ro2 = ChemicalDataAssimilation.JacobianTermRO2[]
@showprogress for drxn ∈ derivatives
    j_terms = JacobianTerms(drxn)
    for j_term ∈ j_terms
        push!(jacobian_terms, j_term)
    end
end

@showprogress for drxn ∈ derivatives_ro2
    j_terms = JacobianTerms(drxn)
    for j_term ∈ j_terms
        push!(jacobian_terms_ro2, j_term)
    end
end




# convert reaction rates into matrix for faster lookup
const K_matrix = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts = df_rrate_coeffs_mech.t

# define values used for update
# idx_t = 0
# tval = -180.0
# ro2_ratio = 1.0


function rhs!(du, u, p, t)
    # set everything to sero
    du .= 0.0

    # get the time index
    tval,idx_t = findmin(x -> abs.(x.- t), ts)

    # get the current ro2_ratio
    ro2_ratio = sum(u[idx_ro2])
    ro2_ratio = ro2_ratio/RO2ᵢ

    # set up product temporary value:
    prod_temp = 1.0

    # update derivatives
    @inbounds for i ∈ 1:length(derivatives)
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end

    prod_temp = 1.0

    @inbounds for i ∈ 1:length(derivatives_ro2)
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_ro2[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end
end


function jac!(Jac, u, p, t)
    # set everything to sero
    Jac .= 0.0

    # get the time index
    tval,idx_t = findmin(x -> abs.(x.- t), ts)

    # get the current ro2_ratio
    ro2_ratio = sum(u₀[idx_ro2])
    ro2_ratio = ro2_ratio/RO2ᵢ

    # set up product temporary value:
    prod_temp = 1.0

    @inbounds for i ∈ 1:length(jacobian_terms)
        update_jacobian!(
            idx_t,
            Jac,
            u,
            jacobian_terms[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end

    prod_temp = 1.0

    @inbounds for i ∈ 1:length(jacobian_terms_ro2)
        update_jacobian!(
            idx_t,
            Jac,
            u,
            jacobian_terms_ro2[i],
            ro2_ratio,
            K_matrix,
            Δt_step,
            prod_temp
        )
    end
end



# this should preallocate the rhs function
du = copy(u₀)
du
rhs!(du, u₀, nothing, -180.0)

du

@benchmark rhs!(du, u₀,nothing, 1.0)  # got it down to 3 allocations!

# set up jacobian prototype as a sparse matrix
idx_pairs = []
for jac_term ∈ jacobian_terms
    push!(idx_pairs, (jac_term.i, jac_term.j))
end
for jac_term ∈ jacobian_terms_ro2
    push!(idx_pairs, (jac_term.i, jac_term.j))
end
idx_pairs = unique(idx_pairs)
size(idx_pairs)
size(derivatives)

I = [idx_pair[1] for idx_pair ∈ idx_pairs]
J = [idx_pair[2] for idx_pair ∈ idx_pairs]
V = zeros(size(I))

using SparseArrays
Jac = sparse(I,J,V)

df_species

@benchmark jac!(Jac, u₀, nothing, 1.0)

# set up jacobian prototype as a non sparse matrix
Jac = zeros(size(u₀,1), size(u₀,1))

@benchmark jac!(Jac, u₀, nothing, 1.0)



tmin = minimum(ts)
tmax = maximum(ts)
#tmax = 0.0

tspan = (tmin, tmax)

# add a bump to each element so we start off nonzero
u₀ .+= 1e-10
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(rhs!, u₀, tspan)

fun = ODEFunction(rhs!; jac=jac!, jac_prototype=Jac)
ode_prob2 = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)


tol = 1e-3
@benchmark solve(ode_prob,
                 CVODE_BDF();
                 saveat=15.0,
                 reltol=tol,
                 abstol=tol,
                 )

@benchmark solve(ode_prob2,
                 CVODE_BDF();
                 saveat=15.0,
                 reltol=tol,
                 abstol=tol,
                 )




# sol = solve(
#     ode_prob,
#     CVODE_BDF();
#     saveat=15.0,
#     reltol=tol,
#     abstol=tol,
# );

sol = solve(
    ode_prob2,
    CVODE_BDF();
    saveat=15.0,
    reltol=tol,
    abstol=tol,
);


size(sol)

df_species
speciesiwant = df_species[1:9, "MCM Name"]
#speciesiwant = ["CO", "O3", "NO2", "APINENE", "H2", "CH4"]
# speciesiwant = ["H2O2"]
#speciesiwant = ["HO2"]


size(sol)


plotting_species_df = df_species[[name ∈ speciesiwant for name ∈ df_species[!, "MCM Name"]], :]

names(plotting_species_df)

p = plot()
i = 1
for row ∈ eachrow(plotting_species_df)
    try
        # plot!(p, sol.t ./ (60*24), sol[row.idx_species, :] .* (nₘ/m_init), label=row.formula, xlabel="t [days]", ylabel="concentration [molecules/cc]", legend=:outertopright)
        plot!(p, sol.t ./ (60), sol[row.idx_species, :] ./ df_params.M, label=row.formula, xlabel="t [days]", ylabel="concentration [mixing ratio]", legend=:outertopright)
    catch e
        println(e)
        println(row)
    end
end

display(p)


