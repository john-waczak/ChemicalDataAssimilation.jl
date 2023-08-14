ENV["GKSwstype"] = 100

println("Setting Up Julia Environment")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# investigate possible uses of StaticArrays.jl, specifically the SizedArray{} decorator
using ChemicalDataAssimilation
using DelimitedFiles, CSV, DataFrames  # file reading
using Plots, StatsPlots # visualization
using Statistics # for mean, std, and var functions
using ProgressMeter
using BenchmarkTools
using LinearAlgebra
using Measurements
using SparseArrays
using ArgParse



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--model_name"
        help = "Name for the resulting model used in output paths"
        arg_type = String
        default = "empty_methane"
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)
    @assert ispath("models/$(parsed_args[:model_name])/EKF/ekf_output.csv")  "Can not find models/$(parsed_args[:model_name])/EKF/ekf_output.csv"
    return parsed_args
end




# 0. parse arguments and set up output directory
println("Parsing command line arguments...")
parsed_args = parse_commandline()
model_name = parsed_args[:model_name]


# --------------------------------------------------------------------------------------------------------------------------
# 1. read in dataframes
# -------------------------------------------------------------------------------------------------------------------------
println("Loading Data into DataFrames...")

df_species = CSV.File("models/$model_name/species.csv") |> DataFrame;
df_params = CSV.File("models/$model_name/state_parameters.csv") |> DataFrame
df_params_ϵ = CSV.File("models/$model_name/state_parameters_ϵ.csv") |> DataFrame
df_number_densities = CSV.File("models/$model_name/number_densities.csv") |> DataFrame;
df_number_densities_ϵ = CSV.File("models/$model_name/number_densities_ϵ.csv") |> DataFrame;


# --------------------------------------------------------------------------------------------------------------------------
# 2. get measurement indices
# --------------------------------------------------------------------------------------------------------------------------
idx_meas::Vector{Int} = Int[]

measurements_to_ignore = [:C2H6, :SO2]  # skip any with nans or negative values
is_meas_in_mechanism = [spec ∈ df_species[!, "MCM Name"] for spec ∈ names(df_number_densities[:, Not([measurements_to_ignore..., :t, :w_ap])])]

df_nd_to_use = df_number_densities[:, Not([measurements_to_ignore..., :t, :w_ap])]
df_nd_to_use = df_nd_to_use[:, is_meas_in_mechanism]

df_nd_to_use_ϵ = df_number_densities_ϵ[:, Not([measurements_to_ignore..., :t, :w_ap])]
df_nd_to_use_ϵ = df_nd_to_use_ϵ[:, is_meas_in_mechanism]

for spec_name ∈ names(df_nd_to_use)
    println(spec_name)
    push!(idx_meas, df_species[df_species[!, "MCM Name"] .== spec_name, :idx_species][1])
end

# --------------------------------------------------------------------------------------------------------------------------
# 3. read in EKF outputs
# --------------------------------------------------------------------------------------------------------------------------

df_ua = CSV.File("models/$model_name/EKF/ekf_output.csv") |> DataFrame
df_ua_ϵ = CSV.File("models/$model_name/EKF/ekf_ϵ_output.csv") |> DataFrame
df_w = CSV.File("models/$model_name/EKF/ekf_measurements.csv") |> DataFrame
df_w_ϵ = CSV.File("models/$model_name/EKF/ekf_measurements_ϵ.csv") |> DataFrame



ts = df_ua.times
ua_mr_vals = Matrix(df_ua[:, Not([:times])])'
ua_mr_ϵ = Matrix(df_ua_ϵ[:, Not([:times])])'
W_mr_val = Matrix(df_w)'
W_mr_ϵ = Matrix(df_w_ϵ)'




# --------------------------------------------------------------------------------------------------------------------------
# 4. Generate time-series plots & report
# --------------------------------------------------------------------------------------------------------------------------
df_species[3,:]
idx_0 = findfirst(x -> x == 0.0, ts)

# set output for quarto documents
report_basepath = joinpath("models/$model_name/EKF")
reportpath = joinpath(report_basepath, "summary.qmd")

# open report file
f = open(reportpath, "w")

model_name = "high_unprimed_methane"
model_name_pretty = join(uppercasefirst.(split(model_name, "_")[1:end-1]), " ")

# write header text
qmd_header = """
---
title: \"ActivePure Research Labs Chamber Assimilation\"
subtitle: \"Test: $(model_name_pretty)\"
date: today
format:
    pdf:
      keep-tex: true
reference-location: margin
citation-location: margin
---
"""
println(f, qmd_header)

# generate plots and summary info:
@showprogress for i ∈ 1:nrow(df_species)
    plot_spec_name = df_species[df_species.idx_species .== i, "MCM Name"][1]
    units, unit_mult = get_reasonable_mr_units(ua_mr_vals[i,:])

    mean_noap = round(mean(ua_mr_vals[i,2:idx_0] .* unit_mult), sigdigits=3)
    median_noap = round(median(ua_mr_vals[i,2:idx_0] .* unit_mult), sigdigits=3)
    std_noap = round(std(ua_mr_vals[i,2:idx_0] .* unit_mult), sigdigits=3)
    max_noap = round(maximum(ua_mr_vals[i,2:idx_0] .* unit_mult), sigdigits=3)
    min_noap = round(minimum(ua_mr_vals[i,2:idx_0] .* unit_mult), sigdigits=3)

    mean_wap = round(mean(ua_mr_vals[i,idx_0:end] .* unit_mult), sigdigits=3)
    median_wap = round(median(ua_mr_vals[i,idx_0:end] .* unit_mult), sigdigits=3)
    std_wap = round(std(ua_mr_vals[i,idx_0:end] .* unit_mult), sigdigits=3)
    max_wap = round(maximum(ua_mr_vals[i,idx_0:end] .* unit_mult), sigdigits=3)
    min_wap = round(minimum(ua_mr_vals[i,idx_0:end] .* unit_mult), sigdigits=3)

    p1 = plot(
        ts[2:idx_0],
        ua_mr_vals[i,2:idx_0] .* unit_mult,
        ribbon=ua_mr_ϵ[i,2:idx_0] .* unit_mult,
        xlabel="time [minutes]",
        ylabel="$plot_spec_name concentration [$(units)]",
        label="No ActivePure",
        lw=3,
        link=:y
    )
    plot!(
        ts[idx_0:end],
        ua_mr_vals[i,idx_0:end] .* unit_mult,
        ribbon=ua_mr_ϵ[i,idx_0:end] .* unit_mult,
        label="With ActivePure",
        lw=3,
    )

    if i ∈ idx_meas
        idx_to_use = findfirst(x->x==i, idx_meas)
        scatter!(ts[2:end],
                 W_mr_val[idx_to_use,2:end] .* unit_mult,
                 yerror=W_mr_ϵ[idx_to_use,2:end] .* unit_mult,
                 label="Measurements",
                 )
    end


    p2 = violin([plot_spec_name], ua_mr_vals[i,2:idx_0] .* unit_mult, side=:left, label="No ActivePure", ymirror=true)
    violin!([plot_spec_name], ua_mr_vals[i,idx_0:end] .* unit_mult , side=:right, label="With ActivePure")

    mag = 1.35
    plot(p1,p2, layout=grid(1,2,widths=[0.7,0.3]), size=(mag*600, mag*400), margins=5Plots.mm)

    savepath = "models/$model_name/EKF/$(plot_spec_name)"
    savepath_png = savepath * ".png"
    savepath_pdf = savepath * ".pdf"

    savefig(savepath_png)
    savefig(savepath_pdf)

    println(f, "\n")
    println(f, "# $(plot_spec_name)\n")
    println(f, "![$(plot_spec_name)]($(plot_spec_name).pdf)\n")
    # write out statistics for No AP vs w AP
    table_text = """
|   | mean | median | standard deviation | min | max | units |
|---|:----:|:------:|:------------------:|:---:|:---:|:-----:|
|Without ActivePure|$(mean_noap)|$(median_noap)|$(std_noap)|$(min_noap)|$(max_noap)|$(units)|
|With ActivePure|$(mean_wap)|$(median_wap)|$(std_wap)|$(min_wap)|$(max_wap)|$(units)|

"""
    println(f, table_text)
    println(f, "\n")

end



close(f)



# --------------------------------------------------------------------------------------------------------------------------
# 5. Generate lifetime plots
# --------------------------------------------------------------------------------------------------------------------------
df_τ = CSV.File("models/$model_name/EKF/lifetimes.csv") |> DataFrame
τs = Matrix(df_τ[!, Not([:t])])'

τs

idx_0 = findfirst(x -> x == 0.0, ts)

@showprogress for i ∈ 1:nrow(df_species)
    plot_spec_name = df_species[df_species.idx_species .== i, "MCM Name"][1]

    units, units_mult = get_reasonable_time_units(τs[i,:])

    try

        p1 = plot(
            ts[2:idx_0],
            τs[i,2:idx_0] .* units_mult,
            xlabel="time [minutes]",
            ylabel="$plot_spec_name lifetime [$(units)]",
            lw=3,
            label="No ActivePure",
        )


        plot!(
            ts[idx_0:end],
            τs[i,idx_0:end] .* units_mult,
            label="With ActivePure",
            lw=3,
        )

        p2 = violin([plot_spec_name], [τ * units_mult for τ ∈ τs[i,2:idx_0] if !isnan(τ) && !isinf(τ)] , side=:left, label="No ActivePure", ymirror=true)
        violin!([plot_spec_name], [τ * units_mult for τ ∈ τs[i,idx_0:end] if !isnan(τ) && !isinf(τ)] , side=:right, label="With ActivePure")

        mag = 1.35
        plot(p1,p2, layout=grid(1,2,widths=[0.7,0.3]), size=(mag*600, mag*400), margins=5Plots.mm)

        savefig("models/$model_name/EKF/$(plot_spec_name)_lifetime.png")
        savefig("models/$model_name/EKF/$(plot_spec_name)_lifetime.pdf")

    catch e
        println(plot_spec_name*" failed to plot!")
        println(e)
    end
end








# --------------------------------------------------------------------------------------------------------------------------
# Reaction Network Visualization
# --------------------------------------------------------------------------------------------------------------------------


# using SparseArrays
# using GraphRecipes

# N_raw = CSV.File("models/$model_name/N.csv") |> DataFrame
# N = sparse(N_raw.I, N_raw.J, N_raw.V)

# # now remove the old file from memory
# N_raw = nothing
# GC.gc()

# N[5,:]
# df_species[5,:]

# to construct graph, we loop over each reaction and add edges between products and reactants.

# # our graph is really a multigraph, i.e. a graph where we can have n-many edges between each pair of nodes. 

# graphplot([[1,1,2,2],[1,1,1],[1]], names="node_".*string.(1:3), nodeshape=:circle, self_edge_size=0.25)


# rxn = rxns[1]
# fieldnames(typeof(rxn))


# edges_mat = [[] for _ ∈ 1:nrow(df_species)]

# @showprogress for rxn ∈ rxns
#     for i ∈ rxn.idxs_in
#         for j ∈ rxn.idxs_out
#             push!(edges_mat[i], j)
#         end
#     end
# end

# graphplot(edges_mat,
#           names=df_species[!, "MCM Name"],
#           method=:chorddiagram,
# #          fontsize=10,
# #          nodeshape=:circle,
# #          nodesize=0.125,
#           size=(1000,1000)
#           )

# savefig("models/$model_name/network_vis.png")
# savefig("models/$model_name/network_vis.pdf")

# # ideas, instead have a marker for species and a separate marker for reactions. Then we can easily
# # visualize the inputs and outputs of each reaction.
# # we can then scale/color the reaction nodes by mean lifetime (or do a video of the lifetime)
# # but to do so we need to fix the positions of each node so that they don't change between frames
# # further we can color each of the species nodes by it's category i.e. source, scavenger, reactive species, reservoir....

# # further, we can use MetaGraphs.jl to include meta-information such as lifetimes, reaction rate coefficients, nice
# # latexified formula names, etc...



# # do a graph for each individual species including all those reactions which it is involved in.



