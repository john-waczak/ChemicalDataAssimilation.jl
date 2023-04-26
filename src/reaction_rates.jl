function get_rrate_list(fac_dict::Dict)
    rrates = []
    grcs = fac_dict["generic_rate_coefficients"]
    crcs = fac_dict["complex_rate_coefficients"]

    for grc ∈ grcs
        lhs, rhs = split(grc, "=")
        push!(rrates, lhs)
    end

    for crc ∈ crcs
        lhs, rhs = split(crc, "=")
        push!(rrates, lhs)
    end

    return strip.(rrates)
end



#function rc_to_expression(rc_string, rate_list; params::NamedTuple = (T=291.483, P=1013.2))
function rc_to_expression(rc_string, rate_list)
    #param_names = join(String.(keys(params)), ",")
    #ps = "($(param_names))"
    ps = "(T, P, M, O2, N2, H2O)"


    lhs, rhs = split(rc_string, "=")
    lhs_out = strip(lhs)
    lhs = strip(lhs)*ps

    rhs = replace(rhs,
                   ";" => "",
                   "TEMP" => "T",
                   "EXP" => "exp",
                   "LOG10" => "log10",
                   "@" => "^",
                   "**" => "^",
                   "D-" => "e-",
                   "D+" => "e",
                   "D-" => "e-",
                   "D+" => "e",
                   "D7" => "e7",
                   "D17" => "e17",
                   # replace number densities by their functions as well
                   # "M" => "M$(ps)",
                   # "O2" => "O2$(ps)",
                   # "N2" => "N2$(ps)",
                   # "H2O" => "H2O$(ps)",
                  )

    if !isempty(rate_list)
        rate_list_out = [strip(rate)*ps for rate ∈ rate_list]
        rate_replace = Pair.(rate_list, rate_list_out)   # map variable to its initial value
        rhs = replace(rhs, rate_replace...)
    end


    # return valid julia code and original reaction rate name so we can replace it later
    return join([lhs, rhs], " = "), lhs_out
end


#function generate_rrates(fac_dict::Dict; model_name::String="mcm", params::NamedTuple=(T=291.483, P=1013.2))
function generate_rrates_funcs(fac_dict::Dict; model_name::String="mcm")
    outpath = "./models/$(model_name)/rrates.jl"
    outpath2 = "./models/$(model_name)/rate_names.txt"
    if isfile(outpath)
        rm(outpath)
    end

    if isfile(outpath2)
        rm(outpath2)
    end

    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end

    rate_names = get_rrate_list(fac_dict)
    idx_sort = sortperm([length(rate) for rate ∈ rate_names], rev=true)
    rate_names = rate_names[idx_sort]


    # write generic rate coefficients
    grcs = fac_dict["generic_rate_coefficients"]
    open(outpath, "a") do f
        println(f, "# -------------------------")
        println(f, "# generic rate coefficients")
        println(f, "# -------------------------")

        for grc ∈ grcs
            res, lhs = rc_to_expression(grc, rate_names)
            println(f, rstrip(res))
        end

        println(f, "\n\n")
    end


    # write complex rate coefficients
    crcs = fac_dict["complex_rate_coefficients"]
    open(outpath, "a") do f
        println(f, "# -------------------------")
        println(f, "# complex rate coefficients")
        println(f, "# -------------------------")

        for crc ∈ crcs
            #res, lhs = rc_to_expression(crc, rate_names; params)
            res, lhs = rc_to_expression(crc, rate_names)
            println(lhs)
            println(f, rstrip(res))
        end
    end

    # write the rate names to file
    open(outpath2, "a") do f
        for rate ∈ rate_names
            println(f, rate)
        end
    end

    return String.(rate_names)
end




function generate_rrates(fac_dict::Dict, df_params::DataFrame, rate_list; model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./models/$(model_name)/rrate_coeffs.jl"

    if isfile(outpath)
        rm(outpath)
    end

    # make sure output dir exists
    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end

    open(outpath, "w") do f

        println(f,"# this will include *both* generic and complex")
        println(f, "df_rate_coeffs = DataFrame()")
        println(f, "# copy in the times")
        println(f, "df_rate_coeffs.t = df_params.t")
        println(f, "df_rate_coeffs.w_ap = df_params.w_ap")

        for rate ∈ rate_list
            println(f, "df_rate_coeffs[!, \"$rate\"] = $(rate).(df_params.temperature, df_params.pressure, df_params.M, df_params.O2, df_params.N2, df_params.H2O)")
        end

        # save the file
        println(f, "\n\n")
        println(f, "CSV.write(\"./models/$model_name/rrate_coeffs.csv\", df_rate_coeffs)")
    end
end
