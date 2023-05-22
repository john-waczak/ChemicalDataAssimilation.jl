using ChemicalDataAssimilation
using DelimitedFiles, CSV, DataFrames
using Plots, StatsPlots
using Statistics
using DifferentialEquations
using Sundials
using Zygote, ForwardDiff
using SciMLSensitivity
using ProgressMeter
using BenchmarkTools
using Optimization
using OptimizationOptimJL
using OptimizationFlux
using LinearAlgebra
using SparseArrays

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


df_species = CSV.File("models/$model_name/species.csv") |> DataFrame;
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_params_ϵ = CSV.File("models/$model_name/state_parameters_ϵ.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;
df_number_densities_ϵ = CSV.File("models/$model_name/number_densities_ϵ.csv") |> DataFrame;

# --------------------------------------------------------------------------------------------------------------------------
# Generate Initial Conditions
# --------------------------------------------------------------------------------------------------------------------------

measurements_to_ignore = [:C2H6]  # skip any with nans

init_path = "./mechanism-files/initial_concentrations/full.txt"
@assert isfile(init_path)

init_dict = generate_init_dict(init_path, df_params.M[1])
u₀    = zeros(Float64, nrow(df_species))

# first set values based on MCM defaults
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

# next update with initial values of our measured species
df_nd_init = df_number_densities[1, Not([measurements_to_ignore..., :t])]
for name ∈ names(df_nd_init)
    if name ∈ df_species[:, "MCM Name"]
        idx = df_species[df_species[:, "MCM Name"] .== name, :].idx_species[1]
        println("Old val: ", u₀[idx])
        u₀[idx] = df_nd_init[name]
        println("New val: ", u₀[idx])
    end
end

for i ∈ axes(u₀,1)
    println(df_species[i,"MCM Name"], ":\t\t\t", u₀[i])
end


# --------------------------------------------------------------------------------------------------------------------------
# Get RO₂ indices and initial value
# --------------------------------------------------------------------------------------------------------------------------

include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices
ro2_sum = sum(u₀[idx_ro2])



# --------------------------------------------------------------------------------------------------------------------------
# Get reaction rate coefficients
# --------------------------------------------------------------------------------------------------------------------------

df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame;

# make sure everything is a Float64
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end



# --------------------------------------------------------------------------------------------------------------------------
# Get measurement indices
# --------------------------------------------------------------------------------------------------------------------------

const idx_meas::Vector{Int} = Int[]

is_meas_in_mechanism = [spec ∈ df_species[!, "MCM Name"] for spec ∈ names(df_number_densities[:, Not([measurements_to_ignore..., :t, :w_ap])])]

idx_ts_to_use = df_params.t .≤ 0.0
df_nd_to_use = df_number_densities[idx_ts_to_use, Not([measurements_to_ignore..., :t, :w_ap])]
df_nd_to_use_ϵ = df_number_densities_ϵ[idx_ts_to_use, Not([measurements_to_ignore..., :t, :w_ap])]

for spec_name ∈ names(df_nd_to_use)
    println(spec_name)
    push!(idx_meas, df_species[df_species[!, "MCM Name"] .== spec_name, :idx_species][1])
end

idx_meas


# --------------------------------------------------------------------------------------------------------------------------
# Setup constants, preallocated vectors, preallocated arrays, etc...
# --------------------------------------------------------------------------------------------------------------------------
const Δt_step::Float64 = 15.0  # time step in minutes
const rxns::Array{ChemicalDataAssimilation.Reaction} = ChemicalDataAssimilation.Reaction[]
const derivatives::Vector{ChemicalDataAssimilation.DerivativeTerm} = ChemicalDataAssimilation.DerivativeTerm[]
const derivatives_ro2::Vector{ChemicalDataAssimilation.DerivativeTermRO2} = ChemicalDataAssimilation.DerivativeTermRO2[]
const jacobian_terms::Vector{ChemicalDataAssimilation.JacobianTerm} = ChemicalDataAssimilation.JacobianTerm[]
const jacobian_terms_ro2::Vector{ChemicalDataAssimilation.JacobianTermRO2} = ChemicalDataAssimilation.JacobianTermRO2[]
const RO2ᵢ::Float64 =ro2_sum > 0.0 ? ro2_sum : 1.0  # make sure we have at least "1 particle" in RO2 sum to prevent division issues
const K_matrix::Matrix{Float64} = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts::Vector{Float64} = df_rrate_coeffs_mech.t
const fudge_fac::Float64 = 0.5  # for measurement uncertainty
const meas_ϵ::Matrix{Float64} = Matrix(df_nd_to_use_ϵ)'
const W::Matrix{Float64} = Matrix(df_nd_to_use)'

const tmin::Float64 = minimum(ts)
const tmax::Float64 = 0.0 # maximum(ts)
const tol::Float64 = 1e-3
const tspan = (tmin, tmax)

const ϵ::Float64 = 1.0
const ϵ_min::Float64 = 1e-12


const try_solve::Bool = true
const use_background_cov::Bool = false


# --------------------------------------------------------------------------------------------------------------------------
# Generate reactions, derivatives, jacobians
# --------------------------------------------------------------------------------------------------------------------------

# create vector of reaction objects
species, reactions = parse_rxns(fac_dict["reaction_definitions"])

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
# Get rhs and jacobian functions, generate jacobian prototype sparse matrix
# --------------------------------------------------------------------------------------------------------------------------
include("models/$model_name/rhs.jl")
include("models/$model_name/jacobian.jl")

# evaluate each once to precompile
du = copy(u₀)
rhs!(du, u₀, nothing, -180.0)

jac_prototype = generate_jac_prototype(jacobian_terms, jacobian_terms_ro2)


println("jacobian non-zero percentage: $(length(nonzeros(jac_prototype))/length(jac_prototype)*100)%")

jac!(jac_prototype, u₀, nothing, 1.0)



# --------------------------------------------------------------------------------------------------------------------------
#  Set up ODE Defaults
# --------------------------------------------------------------------------------------------------------------------------

fun = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)

if try_solve
    sol = solve(ode_prob, CVODE_BDF(); saveat=15.0, reltol=tol, abstol=tol);
end



# --------------------------------------------------------------------------------------------------------------------------
#  Set Up Loss Function for 4dVar
# --------------------------------------------------------------------------------------------------------------------------

@benchmark Rmat(1,meas_ϵ; fudge_fac=fudge_fac)
@benchmark Rinv(1, meas_ϵ; fudge_fac=fudge_fac)


@assert eltype(W) == Float64
@assert size(W,1) == length(idx_meas)

# try out observation operator
if try_solve
    res = ObsOpMeas(sol[1], idx_meas)
end


# create initial analysis vector with small offset
u0a = u₀ .+ 1.0

# we optimize the logarithm instead to force
# concentrations to stay positive
log_u0a = log.(u0a)

# compute background covariance matrix
const u0b::Vector{Float64} = copy(u0a) # i.e. "background guess"
const B::Matrix{Float64} = diagm((ϵ .* (u₀)) .^2  .+ ϵ_min^2)
const Binv::Matrix{Float64} = inv(B)


function loss(log_u0a)
    # first convert back to number density
    u0a = exp.(log_u0a)

    # remake problem using current value
    _prob = remake(ode_prob; u0=u0a)

    # integrate the model forward
    sol = solve(
        _prob,
        CVODE_BDF();
        saveat=Δt_step,
        reltol=tol,
        abstol=tol,
        sensealg=QuadratureAdjoint(),
        verbose=false
    )

    # compute loss
    l = 0.0

    for j ∈ axes(W,2)
        l += 0.5 * ((W[:,j] .- ObsOpMeas(sol[j], idx_meas))' * Rinv(j, meas_ϵ; fudge_fac=fudge_fac) * (W[:,j] .- ObsOpMeas(sol[j], idx_meas)))[1]
    end

    # optionally, add additional loss term quantifying our belief in the inital condition vector
    if use_background_cov
        l += 0.5*(u0a-u0b)'*Binv*(u0a-u0b)
    end

    return l
end

loss(log_u0a)


if try_solve
    @benchmark loss(log_u0a)
end



# --------------------------------------------------------------------------------------------------------------------------
#  Solve the Optimization Problem
# --------------------------------------------------------------------------------------------------------------------------


losses1::Vector{Float64} = []
losses2::Vector{Float64} = []
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


if try_solve
    Zygote.gradient(loss, log_u0a)
end


# define the optimization function and declare we are using Zygote for AutoGrad
optf = OptimizationFunction((x,p)->loss(x), Optimization.AutoZygote())
opt_prob = Optimization.OptimizationProblem(optf, log_u0a)


# solve first with ADAM which is fast but can get stuck in local minimum
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

for i ∈ axes(u0a_final, 1)
    println(df_species[i, "MCM Name"], "\t\t", u0a_final[i])
end

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
df_out

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
u0a_adjusted =  copy(sol[end])
for i ∈ 1:length(u₀)
    if u₀[i] > 1.0
        println("\tcurrent:\t", u0a_adjusted[i])
        u0a_adjusted[i] = u₀[i]
        println("\tfinal:\t", u0a_adjusted[i])
    end
end


for i ∈ axes(u0a_adjusted, 1)
    println(i, "\t", df_species[i, "MCM Name"], "\t\t", u0a_adjusted[i])
end

# # do optimzation on these adjusted values...
# opt_prob = Optimization.OptimizationProblem(optf, log.(u0a_adjusted))
# opt_sol = solve(opt_prob,
#                 ADAM(0.1);
#                 maxiters=300,
#                 callback=callback)

# # finish it off with the slower, but better BFGS
# prob2 = remake(opt_prob, u0=opt_sol.u)
# opt_sol = solve(prob2,
#                 Optim.BFGS(initial_stepnorm=0.01);
#                 callback=callback2,
#                 allow_f_increases=false,
#                 )
#u0a_adjusted = exp.(opt_sol.u)

df_out = DataFrame(:u₀ => u0a_adjusted)
CSV.write("models/$model_name/4dvar/u0_adjusted.csv", df_out)


_prob = remake(ode_prob; u0=u0a_adjusted)
sol2 = solve(
    _prob,
    CVODE_BDF();
    saveat=Δt_step,
    reltol=tol,
    abstol=tol,
    verbose=false
)



# --------------------------------------------------------------------------------------------------------------------------
#  convert final solution to mixing ratios
# --------------------------------------------------------------------------------------------------------------------------

using Measurements

adjusted = false
if adjusted
    u_sol = Matrix(sol2)
else
    u_sol = Matrix(sol)
end

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



# --------------------------------------------------------------------------------------------------------------------------
#  Generate plots
# --------------------------------------------------------------------------------------------------------------------------

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

df_rrate_coeffs_mech[:, "k_41"]
