function generate_photolysis_rates(data_path::String; model_name::String="mcm", Δt_step::Float64=15.0)
    # if file already exists, delete it
    outpath = "./models/$(model_name)/photolysis_rates.csv"

    if isfile(outpath)
        rm(outpath)
    end

    # make sure output dir exists
    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end

    # load in data
    df_photolysis = CSV.File(data_path) |> DataFrame

    # add column to indicate active pure or not
    idx_wap = [dt .≥ 0.0 for dt ∈ df_photolysis.Δt]
    df_photolysis.w_ap = idx_wap

    # add in time coordinate
    df_photolysis.t = df_photolysis.Δt .* Δt_step
    select!(df_photolysis, Not(:Δt))

    CSV.write(outpath, df_photolysis)
end
