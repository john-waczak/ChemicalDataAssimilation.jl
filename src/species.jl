function generate_name_conversion_table(path_to_names_csv::String)
    df_names = CSV.File(path_to_names_csv) |> DataFrame

    # update with "nice" name for each. Eventually, we should update to auto-generate the LaTeX
    chemical_formula = []
    @showprogress for row ∈ eachrow(df_names)
        if !ismissing(row["InChl Identifier"])
            if occursin("/", row["InChl Identifier"])
                push!(chemical_formula, split(row["InChl Identifier"], "/")[2])
            else
                push!(chemical_formula, missing)
            end
        else
            push!(chemical_formula, missing)
        end
    end

    df_names.formula = chemical_formula

    return df_names
end


function generate_species(fac_dict::Dict)
    rxns = fac_dict["reaction_definitions"]
    species, reactions = parse_rxns(rxns)
    unique_species = [spec for spec ∈ species if spec != nothing]
    return String.(unique_species)  # want String, not SubString
end



function get_species_df(df_names::DataFrame, species::Vector{String})
    df_out = df_names[df_names[!, "MCM Name"] .∈ (species,), :]
    df_out.idx_species = ones(Int, nrow(df_out))

    for row ∈ eachrow(df_out)
        idx = findfirst(x -> x == row["MCM Name"], species)
        row.idx_species = idx
    end

    # sort by the idx
    sort!(df_out, :idx_species)

    # update formula names for broken species
    for row ∈ eachrow(df_out)
        if ismissing(row.formula)
            println("\t", row["MCM Name"], " is missing formula!")

            if row["MCM Name"] == "O"
                println("\tsetting O ⟶ O(3P)")
                row.formula = "O(3P)"
            elseif   row["MCM Name"] == "O1D"
                println("\tsetting O1D ⟶ O(1D)")
                row.formula = "O(1D)"
            elseif   row["MCM Name"] == "NA"
                println("\tsetting NA ⟶ HNO3(aq)")
                row.formula = "HNO3(aq)"
            elseif   row["MCM Name"] == "SA"
                println("\tsetting SA ⟶ H2SO4(aq)")
                row.formula = "H2SO4(aq)"
            elseif   row["MCM Name"] == "CL"
                println("\tsetting CL ⟶ Cl")
                row.formula = "Cl"
            else
                continue
            end
        end
    end

    # update common ones


    return df_out
end




# overload for building from file names
function generate_species_df(path_to_names_csv::String, fac_dict::Dict; model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./models/$(model_name)/species.csv"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./models/$(model_name)")
        mkdir("./models/$(model_name)")
    end

    df_names = generate_name_conversion_table(path_to_names_csv)
    species = generate_species(fac_dict)

    df_species = get_species_df(df_names, species)

    CSV.write(outpath, df_species)
end
