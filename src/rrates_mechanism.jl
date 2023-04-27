function rxnrate_to_expression(rxnrate_string, idx_rate, rate_list)
    rxnrate_string = strip(rxnrate_string)

    rate_list_out = ["df_rrate_coeffs.$rate" for rate ∈ rate_list]
    rate_replace = Pair.(rate_list, rate_list_out)   # map variable to its initial value
    res = replace(rxnrate_string, rate_replace...)

    res = replace(res,
                  ";" => "",
                  "TEMP" => "df_params.temperature",
                  "EXP" => "exp",
                  "LOG10" => "log10",
                  "@" => "^",
                  "**" => "^",
                  "D-" => "e-",
                  "D+" => "e",
                  "D8" => "e8",
                  "D9" => "e9",
                  "D10" => "e10",
                  "D11" => "e11",
                  "D12" => "e12",
                  "D13" => "e13",
                  "D14" => "e14",
                  "D15" => "e15",
                  "D16" => "e16",
                  "D17" => "e17",
                  "D18" => "e18",
                  "D19" => "e19",
                  "J<" => "df_photolysis.J_",
                  "J <" => "df_photolysis.J_",  # weird space for some reason
                  ">" => "",
                  # replace number densities by their functions as well
                  "*M" => "*df_params.M",
                  "*O2" => "*df_params.O2",
                  "*N2" => "*df_params.N2",
                  "*H2O" => "*df_params.H2O",
                  "*RO2" => "*RO2ᵢ"
                  )

    res = "df_mech_rate_coeffs.k_$idx_rate = @. "*res*"*my_ones"

    return res
end


function generate_rrates_mechanism(fac_dict, rate_list; model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./models/$(model_name)/rrates_mechanism.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end

    species, reactions = parse_rxns(fac_dict["reaction_definitions"])
    unique_species = [spec for spec ∈ species if spec != nothing]

    open(outpath, "w") do f
        println(f, "my_ones = ones(nrow(df_params))\n\n")
        println(f, "df_mech_rate_coeffs = DataFrame()")
        println(f, "# copy in the times")
        println(f, "df_mech_rate_coeffs.t = df_params.t")
        println(f, "df_mech_rate_coeffs.w_ap = df_params.w_ap")

        for reaction ∈ reactions
        end

        for i ∈ 1:length(reactions)
            reactants, reactants_stoich, products, products_stoich, rrate_string = reactions[i]
            rrate_string = rxnrate_to_expression(rrate_string, i, rate_list)
            println(f, rrate_string)
        end

        println(f, "\n\n")
        println(f, "# save the file")
        println(f, "CSV.write(\"./models/$model_name/rrate_coeffs_mech.csv\", df_mech_rate_coeffs)")

    end

end


