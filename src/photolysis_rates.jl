function generate_photolysis_rates(no_ap_path::String, w_ap_path::String; model_name::String="mcm", Δt_step::Float64=15.0)
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
    df_photolysis_no_ap = CSV.File(no_ap_path) |> DataFrame
    df_photolysis_w_ap = CSV.File(w_ap_path) |> DataFrame

    # add column to indicate active pure or not
    df_photolysis_no_ap.w_ap = [false for _ ∈ 1:nrow(df_photolysis_no_ap)]
    df_photolysis_w_ap.w_ap = [false for _ ∈ 1:nrow(df_photolysis_w_ap)]

    # join them together
    df_photolysis = vcat(df_photolysis_no_ap[1:end-1,:], df_photolysis_w_ap)

    # add in time coordinate
    df_photolysis.t = df_photolysis.Δt .* Δt_step
    select!(df_photolysis, Not(:Δt))

    CSV.write(outpath, df_photolysis)
end
