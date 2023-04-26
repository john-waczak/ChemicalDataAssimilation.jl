function generate_densities(no_ap_path::String, w_ap_path::String; model_name::String="mcm")
    # if file already exists, delete it
    outpath1 = "./models/$(model_name)/number_densities.csv"
    outpath2 = "./models/$(model_name)/state_parameters.csv"  # i.e. M, O2, N2, H2O + Temp, Pressure

    if isfile(outpath1) || isfile(outpath2)
        rm(outpath1)
        rm(outpath2)
    end

    # make sure output dir exists
    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end


    # load in data
    df_number_densities_no_ap = CSV.File(no_ap_path) |> DataFrame
    df_number_densities_w_ap = CSV.File(w_ap_path) |> DataFrame

    # add column to indicate active pure or not
    df_number_densities_no_ap.w_ap = [false for _ ∈ 1:nrow(df_number_densities_no_ap)]
    df_number_densities_w_ap.w_ap = [true for _ ∈ 1:nrow(df_number_densities_w_ap)]

    # join them together
    df_number_densities = vcat(df_number_densities_no_ap[1:end-1,:], df_number_densities_w_ap)


    # separate into state variables and measurements

    state_params = ["M", "O2", "N2", "H2O"]

    df_state = df_number_densities[:, vcat(state_params, "t", "temperature", "pressure", "w_ap")]
    df_nd = df_number_densities[:, Not(vcat(state_params, "temperature", "pressure"))]

    CSV.write(outpath1, df_nd)
    CSV.write(outpath2, df_state)
end
