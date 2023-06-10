ENV["GKSwstype"] = 100

println("Setting Up Julia Environment")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

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
using ParameterHandling, EarlyStopping
using JSON
using StableRNGs
using ArgParse


# # set random number seed
# rng = StableRNG(42)



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
        "--restart"
            help = "Whether or not to restart 4d_var from previous fitresult"
            action = :store_true
        "--try_solve"
            help = "Whether or not to precompile solvers by calling once."
            action = :store_true
        "--use_background_cov"
            help = "Whether or not to use background covariance matrix in loss"
            action = :store_true
        "--fudge_fac"
            help = "A fudge factor for manipulating scale of measurement uncertainties"
            arg_type = Float64
            default = 0.5
        "--epsilon"
            help = "Estimated background uncertainty for diagonal of B matrix, i.e. uncertainty in initial condition"
            arg_type = Float64
            default = 0.5
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)

    @assert isfile(parsed_args[:mechanism_path]) "Supplied mechanism path does not exist"

    if !ispath("models/$(parsed_args[:model_name])")
        println("$(parsed_args[:model_name]) directory does not exist in `./models`. Creating now...")
        mkpath("models/$(parsed_args[:model_name])")
    end


    @assert ispath("data/no_ap/number_densities.csv") "Can not find  data/no_ap/number_densities.csv"
    @assert ispath("data/w_ap/number_densities.csv") "Can not find  data/w_ap/number_densities.csv"
    @assert ispath("data/no_ap/number_densities_ϵ.csv") "Can not find  data/no_ap/number_densities_ϵ.csv"
    @assert ispath("data/w_ap/number_densities_ϵ.csv") "Can not find  data/w_ap/number_densities_ϵ.csv"

    if parsed_args[:restart]
        @assert ispath("models/$(parsed_args[:model_name])/4dvar/u0.csv")  "Can not find models/$(parsed_args[:model_name])/4dvar/u0.csv"
    end



    return parsed_args
end



# 0. parse arguments and set up output directory
println("Parsing command line arguments...")
parsed_args = parse_commandline()
mechpath = parsed_args[:mechanism_path]
model_name = parsed_args[:model_name]
want_restart = parsed_args[:restart]

if !isdir("models/$model_name/4dvar")
    mkpath("models/$model_name/EKF")
end



# --------------------------------------------------------------------------------------------------------------------------
# 1. read in dataframes
# -------------------------------------------------------------------------------------------------------------------------
println("Loading Data into DataFrames...")

fac_dict = read_fac_file(mechpath)
df_species = CSV.File("models/$model_name/species.csv") |> DataFrame;
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_params_ϵ = CSV.File("models/$model_name/state_parameters_ϵ.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;
df_number_densities_ϵ = CSV.File("models/$model_name/number_densities_ϵ.csv") |> DataFrame;

# --------------------------------------------------------------------------------------------------------------------------
# 2. Generate Initial Conditions
# --------------------------------------------------------------------------------------------------------------------------
println("Loading initial concentrations...")
measurements_to_ignore = [:C2H6]  # skip any with nans
init_path = "./mechanism-files/initial_concentrations/full.txt"
@assert isfile(init_path)

init_dict = generate_init_dict(init_path, df_params.M[1])
u₀    = zeros(Float64, nrow(df_species))

if want_restart
    df_u0 = CSV.File("models/$model_name/4dvar/u0.csv") |> DataFrame
    u₀ = df_u0.u0
else
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
end

println("Initial concentrations")
println("---")
for i ∈ axes(u₀,1)
    println(df_species[i,"MCM Name"], ":\t\t\t", u₀[i])
end
println("---")



# --------------------------------------------------------------------------------------------------------------------------
# 3. Get RO₂ indices and initial value
# --------------------------------------------------------------------------------------------------------------------------
println("Computing initial RO2 value...")

include("./models/$model_name/ro2.jl")
idx_ro2  # this is the array w/ ro2 indices
ro2_sum = sum(u₀[idx_ro2])




# --------------------------------------------------------------------------------------------------------------------------
# 4. Get reaction rate coefficients
# --------------------------------------------------------------------------------------------------------------------------
println("Reading")

df_rrate_coeffs_mech = CSV.File("./models/$model_name/rrate_coeffs_mech.csv") |> DataFrame;

# make sure everything is a Float64
for i ∈ 3:ncol(df_rrate_coeffs_mech)
    if eltype(df_rrate_coeffs_mech[:,i]) != Float64
        println("Fixing ", names(df_rrate_coeffs_mech)[i])
        df_rrate_coeffs_mech[!,i] = parse.(Float64, df_rrate_coeffs_mech[:,i])
    end
end



# --------------------------------------------------------------------------------------------------------------------------
# 5. Get measurement indices
# --------------------------------------------------------------------------------------------------------------------------
println("Determining measurement indices...")
const idx_meas::Vector{Int} = Int[]

is_meas_in_mechanism = [spec ∈ df_species[!, "MCM Name"] for spec ∈ names(df_number_densities[:, Not([measurements_to_ignore..., :t, :w_ap])])]

idx_ts_to_use = df_params.t .≤ 0.0
df_nd_to_use = df_number_densities[idx_ts_to_use, Not([measurements_to_ignore..., :t, :w_ap])]
df_nd_to_use_ϵ = df_number_densities_ϵ[idx_ts_to_use, Not([measurements_to_ignore..., :t, :w_ap])]

for spec_name ∈ names(df_nd_to_use)
    println(spec_name)
    push!(idx_meas, df_species[df_species[!, "MCM Name"] .== spec_name, :idx_species][1])
end



# --------------------------------------------------------------------------------------------------------------------------
# 6. Setup constants, preallocated vectors, preallocated arrays, etc...
# --------------------------------------------------------------------------------------------------------------------------
println("Preallocating outputs...")
const Δt_step::Float64 = parsed_args[:time_step]  # time step in minutes
const rxns::Array{ChemicalDataAssimilation.Reaction} = ChemicalDataAssimilation.Reaction[]
const derivatives::Vector{ChemicalDataAssimilation.DerivativeTerm} = ChemicalDataAssimilation.DerivativeTerm[]
const derivatives_ro2::Vector{ChemicalDataAssimilation.DerivativeTermRO2} = ChemicalDataAssimilation.DerivativeTermRO2[]
const jacobian_terms::Vector{ChemicalDataAssimilation.JacobianTerm} = ChemicalDataAssimilation.JacobianTerm[]
const jacobian_terms_ro2::Vector{ChemicalDataAssimilation.JacobianTermRO2} = ChemicalDataAssimilation.JacobianTermRO2[]
const RO2ᵢ::Float64 =ro2_sum > 0.0 ? ro2_sum : 1.0  # make sure we have at least "1 particle" in RO2 sum to prevent division issues
const K_matrix::Matrix{Float64} = Matrix{Float64}(df_rrate_coeffs_mech[:, 3:end])
const ts::Vector{Float64} = df_rrate_coeffs_mech.t
const fudge_fac::Float64 = parsed_args[:fudge_fac]  # for measurement uncertainty
const meas_ϵ::Matrix{Float64} = Matrix(df_nd_to_use_ϵ)'
const W::Matrix{Float64} = Matrix(df_nd_to_use)'

const tmin::Float64 = minimum(ts)
const tmax::Float64 = 0.0 # maximum(ts)
const tol::Float64 = 1e-3
const tspan = (tmin, tmax)

const ϵ::Float64 = parsed_args[:epsilon]  # 1.0
const ϵ_min::Float64 = 1e-12

const try_solve::Bool = parsed_args[:try_solve]
const use_background_cov::Bool = parsed_args[:use_background_cov]



# --------------------------------------------------------------------------------------------------------------------------
# 7. Generate reactions, derivatives, jacobians
# --------------------------------------------------------------------------------------------------------------------------
println("Generating reaction lists, derivatives, jacobians...")

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
# 8. Get rhs and jacobian functions, generate jacobian prototype sparse matrix
# --------------------------------------------------------------------------------------------------------------------------
println("Fetching rhs and jac functions...")
include("models/$model_name/rhs.jl")
include("models/$model_name/jacobian.jl")

# evaluate each once to precompile
du = copy(u₀)
rhs!(du, u₀, nothing, -180.0)

jac_prototype = generate_jac_prototype(jacobian_terms, jacobian_terms_ro2)


println("jacobian non-zero percentage: $(length(nonzeros(jac_prototype))/length(jac_prototype)*100)%")

jac!(jac_prototype, u₀, nothing, 1.0)



# --------------------------------------------------------------------------------------------------------------------------
# 9. Set up ODE Defaults
# --------------------------------------------------------------------------------------------------------------------------
println("Setting up ODE defaults...")
fun = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)

if try_solve
    sol = solve(ode_prob, CVODE_BDF(); saveat=15.0, reltol=tol, abstol=tol);
end


# --------------------------------------------------------------------------------------------------------------------------
# 10. Set Up Loss Function for 4dVar
# --------------------------------------------------------------------------------------------------------------------------
println("Setting up loss function...")
@benchmark Rmat(1,meas_ϵ; fudge_fac=fudge_fac)
@benchmark Rinv(1, meas_ϵ; fudge_fac=fudge_fac)


@assert eltype(W) == Float64
@assert size(W,1) == length(idx_meas)


# try out observation operator
if try_solve
    res = ObsOpMeas(sol[1], idx_meas)
end


# create initial analysis vector with small offset
current_u0a = (;u0 = [positive(u+1.0) for u ∈ u₀ ])
u0a, unflatten = ParameterHandling.value_flatten(current_u0a)
# NOTE: `unflatten` both reconstructs the NamedTuple with our paramters and applies inverse transform to Positive


# compute background covariance matrix
const u0b::Vector{Float64} = copy(u0a) # i.e. "background guess"
const B::Matrix{Float64} = diagm((ϵ .* (u₀)) .^2  .+ ϵ_min^2)
const Binv::Matrix{Float64} = inv(B)


#function loss(log_u0a)
function loss(u0a)
    u0a_now = unflatten(u0a)

    # first convert back to number density
    # u0a = exp.(log_u0a)

    # remake problem using current value
    #_prob = remake(ode_prob; u0=u0a)

    _prob = remake(ode_prob; u0=u0a_now.u0)

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

if try_solve
    @benchmark loss(u0a)
end



# --------------------------------------------------------------------------------------------------------------------------
#  Solve the Optimization Problem
# --------------------------------------------------------------------------------------------------------------------------
println("Defining callback functions...")

# set up stoppers for each callback function
stopper1 = EarlyStopper(Disjunction(Warmup(Patience(n=5);n=3), Warmup(NumberSinceBest(n=10);n=10), TimeLimit(t=24.0)))

losses1::Vector{Float64} = []
losses2::Vector{Float64} = []


df_out = DataFrame(unflatten(u0a))
CSV.write("models/$model_name/4dvar/u0.csv", df_out)


function callback(u0a, lossval)
    println("current loss: ", lossval)
    df_out.u0 = unflatten(u0a).u0
    CSV.write("models/$model_name/4dvar/u0.csv", df_out)
    push!(losses1, lossval)
    done!(stopper1, lossval)  # return false unless we've matched a stopping criterion
end


if try_solve
    Zygote.gradient(loss, log_u0a)
end


# define the optimization function and declare we are using Zygote for AutoGrad
optf = OptimizationFunction((x,p)->loss(x), Optimization.AutoZygote())
opt_prob = Optimization.OptimizationProblem(optf, u0a)

println("First round of optimization:")
# solve first with ADAM which is fast but can get stuck in local minimum

method1 = ADAM(0.1)
method2 = BFGS(initial_stepnorm=0.01)

#method2=LBFGS()  # <-- took bad steps
if want_restart
    method1 = ADAM(0.005)
end


opt_sol = solve(opt_prob,
                method1;
                maxiters=300,
                callback=callback)
u0a = opt_sol.u
# finish it off with the slower, but better BFGS
prob2 = remake(opt_prob, u0=u0a)
stopper2 = EarlyStopper(Disjunction(Warmup(Patience(n=5);n=3), Warmup(NumberSinceBest(n=30);n=3), TimeLimit(t=24.0)))
function callback2(u0a, lossval)
    println("current loss: ", lossval)
    df_out.u0 = unflatten(u0a).u0
    CSV.write("models/$model_name/4dvar/u0.csv", df_out)
    push!(losses2, lossval)
    done!(stopper2, lossval)  # return false unless we've matched a stopping criterion
end

println("Second round of optimization:")
opt_sol = solve(prob2,
                method2,
                callback=callback2,
                allow_f_increases=false,
                )

# see example here: https://docs.sciml.ai/DiffEqFlux/stable/examples/neural_ode/

u0a_final = unflatten(opt_sol.u).u0

for i ∈ axes(u0a_final, 1)
    println(df_species[i, "MCM Name"], "\t\t", u0a_final[i])
end


# save output to a file

if !isdir("models/$model_name/4dvar")
    mkpath("models/$model_name/4dvar")
end

println("Plotting loss curve...")
iterations1 = 1:length(losses1)
iterations2 = (iterations1[end]):(iterations1[end]+length(losses2)-1)
plot(iterations1, losses1, xlabel="iteration", ylabel="loss", title="4D-Var Training", lw=3, label="ADAM")
plot!(iterations2, losses2, lw=3, label="BFGS")

savefig("models/$model_name/4dvar/training_loss.png")
savefig("models/$model_name/4dvar/training_loss.pdf")



df_out = DataFrame(:u0 => u0a_final)

CSV.write("models/$model_name/4dvar/u0.csv", df_out)



_prob = remake(ode_prob; u0=u0a_final)
sol = solve(
    _prob,
    CVODE_BDF();
    saveat=Δt_step,
    reltol=tol,
    abstol=tol,
    verbose=false
)


df_out = DataFrame()
@showprogress for i ∈ 1:nrow(df_species)
    df_out[!, df_species[i, "MCM Name"]] = sol[i,:]
end

df_out[1,:]
CSV.write("models/$model_name/4dvar/u0_integrated.csv", df_out)

