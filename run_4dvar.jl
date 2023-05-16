using ChemicalDataAssimilation

# file reading
using DelimitedFiles, CSV, DataFrames

# visualization
using Plots, StatsPlots

# for mean, std, and var functions
using Statistics

# for ODE problem and integration
using DifferentialEquations
using Sundials

# for autodiff of ODE problem solutions
using Zygote, ForwardDiff
using SciMLSensitivity

# convenience packages
using ProgressMeter
using BenchmarkTools

# for optimization of 4dvar cost function
using Optimization
using OptimizationOptimJL
using OptimizationFlux

# for generating specialized matrix types
using LinearAlgebra


# --------------------------------------------------------------------------------------------------------------------------
# Setup paths
# -------------------------------------------------------------------------------------------------------------------------

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

df_species = CSV.File("models/$model_name/species.csv") |> DataFrame;
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_params_ϵ = CSV.File("models/$model_name/state_parameters_ϵ.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;
df_number_densities_ϵ = CSV.File("models/$model_name/number_densities_ϵ.csv") |> DataFrame;

# df_photolysis = CSV.File("models/$model_name/photolysis_rates.csv") |> DataFrame
# df_rrate_coeffs = CSV.File("./models/$model_name/rrate_coeffs.csv") |> DataFrame;


# --------------------------------------------------------------------------------------------------------------------------
# Generate Initial Conditions
# --------------------------------------------------------------------------------------------------------------------------

df_number_densities
measurements_to_ignore = [:C2H6]  # skip any with nans

init_path = "./mechanism-files/initial_concentrations/full.txt"
@assert isfile(init_path)

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

df_nd_init = df_number_densities[1, Not([measurements_to_ignore..., :t])]

for name ∈ names(df_nd_init)
    if name ∈ df_species[:, "MCM Name"]
        idx = df_species[df_species[:, "MCM Name"] .== name, :].idx_species[1]
        println("Old val: ", u₀[idx])
        u₀[idx] = df_nd_init[name]
        println("New val: ", u₀[idx])
    end
end


# --------------------------------------------------------------------------------------------------------------------------
# Get RO₂ indices and initial value
# --------------------------------------------------------------------------------------------------------------------------
include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices

ro2_sum = sum(u₀[idx_ro2])
const RO2ᵢ =ro2_sum > 0 ? ro2_sum : 1.0  # make sure we have at least "1 particle"


# --------------------------------------------------------------------------------------------------------------------------
# Get reaction rate coefficients
# --------------------------------------------------------------------------------------------------------------------------

df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame;
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end


# --------------------------------------------------------------------------------------------------------------------------
# Generate reactions, derivatives, jacobians
# --------------------------------------------------------------------------------------------------------------------------

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

println("num reactions: ", length(rxns))
@assert length(rxns) == length(reactions)

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

println("num derivative terms: ", length(derivatives) + length(derivatives_ro2))


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


println("num jacobian terms: ", size(jacobian_terms,1) + size(jacobian_terms_ro2,1))



# --------------------------------------------------------------------------------------------------------------------------
# Generate lookup table for reaction rate coefficients and time values
# --------------------------------------------------------------------------------------------------------------------------

const K_matrix = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts = df_rrate_coeffs_mech.t

@assert size(K_matrix, 1) == length(ts)
@assert size(K_matrix, 2) == length(rxns)


# --------------------------------------------------------------------------------------------------------------------------
# Get rhs and jacobian functions, generate jacobian prototype sparse matrix
# --------------------------------------------------------------------------------------------------------------------------
include("models/$model_name/rhs.jl")
include("models/$model_name/jacobian.jl")

# evaluate each once to precompile
du = copy(u₀)
rhs!(du, u₀, nothing, -180.0)

jac_prototype = generate_jac_prototype(jacobian_terms, jacobian_terms_ro2)

jac!(jac_prototype, u₀, nothing, 1.0)



# --------------------------------------------------------------------------------------------------------------------------
#  Set up ODE Defaults
# --------------------------------------------------------------------------------------------------------------------------

tmin = minimum(ts)
#tmax = maximum(ts)
tmax = 0.0  # only go until t=0 for 4dVar
tspan = (tmin, tmax)
tol = 1e-3

#ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(rhs!, u₀, tspan)
#fun = ODEFunction(rhs!; jac=jac!, jac_prototype=Jac)
fun = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)


sol = solve(
    ode_prob,
    CVODE_BDF();
    saveat=15.0,
    reltol=tol,
    abstol=tol,
);



# --------------------------------------------------------------------------------------------------------------------------
#  Set Up Loss Function for 4dVar
# --------------------------------------------------------------------------------------------------------------------------


measurements_to_ignore

idx_ts_to_use = ts .≤ 0.0
df_nd_to_use = df_number_densities[idx_ts_to_use, Not([measurements_to_ignore..., :t, :w_ap])]
df_species


# set up measurement covariance matrix
#meas_ϵ = collect(df_number_densities[1, Not([measurements_to_ignore..., :t, :w_ap])])
const meas_ϵ = Matrix(df_number_densities_ϵ[:, Not([measurements_to_ignore..., :t, :w_ap])])'

const fudge_fac = 1.0
#const fudge_fac = 0.5
#const fudge_fac = 0.25
@benchmark Rmat(1,meas_ϵ; fudge_fac=fudge_fac)
@benchmark Rinv(1, meas_ϵ; fudge_fac=fudge_fac)


# fudge_fac = 0.1
# fudge_fac = 1.0
# const R = diagm( (fudge_fac .* meas_ϵ) .^ 2)
# const Rinv = inv(R)

# σ_b = 0.1
# B = diagm( σ_b^2 .* u₀ .^ 2)


const idx_meas = Int[]
for spec_name ∈ names(df_nd_to_use)
    push!(idx_meas, df_species[df_species[!, "MCM Name"] .== spec_name, :idx_species][1])
end

size(sol)


W = collect(Matrix(df_nd_to_use)')  # make it match the shape of solution object which is (n_dims, n_times)
size(W)
@assert eltype(W) == Float64
@assert size(W,1) == length(idx_meas)

res = ObsOpMeas(sol, idx_meas)
W



minimum(u₀[u₀ .> 0])
# u0a = u₀ .+ 5e8  # what's reasonable here?
# u0a = u₀ .+ 1e7
# u0a = u₀ .+ 1e9
#u0a = u₀ .+ 1000
u0a = u₀ .+ 1.0
log_u0a = log.(u0a)


const u0b = copy(u0a) # i.e. "background guess"
# compute background covariance matrix
const B = diagm((fudge_fac .* (u₀ .+ 1)) .^2)
const Binv = inv(B)
println(Binv)


function loss(log_u0a)
    u0a = exp.(log_u0a)
    _prob = remake(ode_prob; u0=u0a)
    sol = solve(
        _prob,
        CVODE_BDF();
        saveat=Δt_step,
        reltol=tol,
        abstol=tol,
        sensealg=QuadratureAdjoint(),
        verbose=false
    )


    l = 0.0
    for j ∈ axes(W,2)
        l += 0.5 * ((W[:,j] .- ObsOpMeas(sol[:,j], idx_meas))' * Rinv(j, meas_ϵ; fudge_fac=fudge_fac) * (W[:,j] .- ObsOpMeas(sol[:,j], idx_meas)))[1]
    end

    # add in term for B
    # l += 0.5*(u0a-u0b)'*Binv*(u0a-u0b)

    return l
end


losses1 = []
losses2 = []
function callback(p, lossval)
    println("current loss: ", lossval)
    push!(losses1, lossval)
    false
end

function callback2(p, lossval)
    println("current loss: ", lossval)
    push!(losses2, lossval)
    false
end


# @benchmark loss(log_u0a)
loss(log_u0a)
Zygote.gradient(loss, log_u0a)


#ParameterHandling.value(u0a)
#ParameterHandling.flatten(u0a)

optf = OptimizationFunction((x,p)->loss(x), Optimization.AutoZygote())
#opt_prob = Optimization.OptimizationProblem(optf, u0a, lb=zeros(size(u₀)), ub=1e50*ones(size(u₀))) #, ub= [Inf for _ ∈ 1:length(u₀)])
opt_prob = Optimization.OptimizationProblem(optf, log_u0a)

opt_sol = solve(opt_prob,
                ADAM(0.1);
                maxiters=300,
                callback=callback)

# finish it off with the slower, but better BFGS
prob2 = remake(opt_prob, u0=opt_sol.u)
opt_sol = solve(prob2,
                Optim.BFGS(initial_stepnorm=0.01);
                callback=callback2,
                allow_f_increases=false,
                )

# see example here: https://docs.sciml.ai/DiffEqFlux/stable/examples/neural_ode/

u0a_final = exp.(opt_sol.u)


u0a_final

df_species
# save output to a file

if !isdir("models/$model_name/4dvar")
    mkpath("models/$model_name/4dvar")
end


iterations1 = 1:length(losses1)
iterations2 = (iterations1[end]+1):(iterations1[end]+length(losses2))
plot(iterations1, losses1, xlabel="iteration", ylabel="loss", title="4D-Var Training", lw=3, label="ADAM")

plot!(iterations2, losses2, lw=3, label="BFGS")
savefig("models/$model_name/4dvar/training_loss.png")
savefig("models/$model_name/4dvar/training_loss.pdf")



df_out = DataFrame(:u₀ => u0a_final)
CSV.write("models/$model_name/4dvar/u0.csv", df_out)


# generate solution with final values:
_prob = remake(ode_prob; u0=u0a_final)
sol = solve(
    _prob,
    CVODE_BDF();
    saveat=Δt_step,
    reltol=tol,
    abstol=tol,
    verbose=false
)


# try the whole thing again but starting with the final value from the assimilation
u0a_adjusted =  copy(Matrix(sol)[:,end])
for i ∈ 1:length(u₀)
    if u₀[i] != 0.0
        u0a_adjusted[i] = u₀[i]
    end
end

u0a_adjusted


opt_prob = Optimization.OptimizationProblem(optf, log.(u0a_adjusted))
opt_sol = solve(opt_prob,
                ADAM(0.1);
                maxiters=300,
                callback=callback)

# finish it off with the slower, but better BFGS
prob2 = remake(opt_prob, u0=opt_sol.u)
opt_sol = solve(prob2,
                Optim.BFGS(initial_stepnorm=0.01);
                callback=callback2,
                allow_f_increases=false,
                )

u0a_adjusted = exp.(opt_sol.u)
df_out = DataFrame(:u₀ => u0a_adjusted)
CSV.write("models/$model_name/4dvar/u0_adjusted.csv", df_out)


_prob = remake(ode_prob; u0=u0a_adjusted)
sol = solve(
    _prob,
    CVODE_BDF();
    saveat=Δt_step,
    reltol=tol,
    abstol=tol,
    verbose=false
)



using Measurements


u_sol = Matrix(sol)
times = tmin:Δt_step:tmax



spec_names = df_species[idx_meas, "MCM Name"]
u_meas = u_sol[idx_meas,:]
# collect total number density
M = df_params[df_params.t .≤ 0.0, :M]
M_ϵ = df_params_ϵ[df_params.t .≤ 0.0, :M]

# convert number density to mixing ratio
u_mr = copy(u_meas)
for j ∈ axes(u_mr,2)
    u_mr[:,j] .= u_mr[:,j] ./ M[j]
end

# convert to ppb
u_ppb = u_mr .* 1e9


# combine measurements w/ uncertainty
W_ϵ = similar(W)
for i ∈ 1:length(idx_meas)
    W_ϵ[i,:] .= [sqrt.(Rmat(t, meas_ϵ; fudge_fac=fudge_fac)[i,i]) for t ∈ 1:size(W,2)]
end

# convert to mixing ratio
W_mr = W .± W_ϵ
for j ∈ axes(W_mr, 2)
    W_mr[:,j] .= W_mr[:,j] ./ (M[j] ± M_ϵ[j])
end

# convert to ppb
W_ppb = W_mr .* 1e9



# # visualize the results:
# @showprogress for i ∈ 1:length(idx_meas)
#     times = tmin:Δt_step:tmax
#     plot_spec_name = df_species[idx_meas[i], "MCM Name"]
#     plot(times,
#          sol[i,:],
#          label="4dVar",
#          title=plot_spec_name,
#          lw=3
#         )

#     scatter!(times, W[i,:],
#             yerror=[sqrt(Rmat(t,meas_ϵ;fudge_fac=fudge_fac)[i,i]) for t ∈ 1:size(W,2)],
#             label="Measurements",
#             )
#     xlabel!("time [minutes]")
#     ylabel!("concentration [molecules/cm³]")

#     savefig("models/$model_name/4dvar/$(plot_spec_name).png")
#     savefig("models/$model_name/4dvar/$(plot_spec_name).pdf")
# end


# visualize the results:
@showprogress for i ∈ 1:length(idx_meas)
    plot_spec_name = df_species[idx_meas[i], "MCM Name"]
    plot(times,
         u_ppb[i,:],
         label="4dVar",
         title=plot_spec_name,
         lw=3
         )

    scatter!(times,
             Measurements.value.(W_ppb[i,:]),
             yerror=Measurements.uncertainty.(W_ppb[i,:]),
             label="Measurements",
             )
    xlabel!("time [minutes]")
    ylabel!("concentration [ppb]")

    savefig("models/$model_name/4dvar/$(plot_spec_name)_ppb.png")
    savefig("models/$model_name/4dvar/$(plot_spec_name)_ppb.pdf")
end


u0a_final .- u0a_adjusted


# idx_meas
# idx_not_meas = [idx for idx∈1:length(u₀) if !(idx∈idx_meas)]
# u0a_final[idx_not_meas] .- u₀[idx_not_meas]
# # u0a_final

# i=idx_meas[1]
# u0a_final[idx_meas[1]]
# println(mean([sqrt(Rmat(t,meas_ϵ;fudge_fac=fudge_fac)[i,i]) for t ∈ 1:size(W,2)]))





