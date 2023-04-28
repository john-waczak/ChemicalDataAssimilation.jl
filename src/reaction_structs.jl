# Define Reaction Type Hierarchy
abstract type Reaction end

struct Photodissociation <: Reaction
    idxs_in::Vector{Int}
    idxs_out::Vector{Int}
    idx_k::Int
end

struct CollisionReaction <: Reaction
    idxs_in::Vector{Int}
    idxs_out::Vector{Int}
    idx_k::Int
    needs_ro2::Bool
end


function parse_rxn(reaction, idx_reaction, df_species)
    needs_ro2 = false
    is_photolysis = false

    reactants, reactants_stoich, products, products_stoich, rrate_string = reaction

    idx_reactants = get_spec_idx(String.(reactants), df_species)
    idx_products = products == [nothing] ? Int[] : get_spec_idx(String.(products), df_species)

    if occursin("*RO2", rrate_string)
        needs_ro2=true
    end

    if occursin("J<", rrate_string) || occursin("J <", rrate_string)
        is_photolysis = true
    end

    if is_photolysis
        return Photodissociation(
            idx_reactants,
            idx_products,
            idx_reaction
        )
    else
        return CollisionReaction(
            idx_reactants,
            idx_products,
            idx_reaction,
            needs_ro2
        )
    end
end



function generate_reaction_list(fac_dict::Dict, df_species::DataFrame)
    species, reactions = parse_rxns(fac_dict["reaction_definitions"])
    rxns = ChemicalDataAssimilation.Reaction[]
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
    return rxns
end



function reaction_rate_coeff(rxn::Reaction, time::Float64, ro2_ratio::Float64)
    return nothing
end


function rhs_update!(u, rxn::Reaction)
    return nothing
end
