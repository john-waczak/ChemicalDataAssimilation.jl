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


# abstract type Photodissociation <: Reaction end
# abstract type CollisionReaction <: Reaction end


# # Define Concrete Types
# struct Photodissociation11 <: Photodissociation
#     idx_r1::Int
#     idx_p1::Int

#     # reaction rate index
#     idx_k::Int
# end

# struct Photodissociation12 <: Photodissociation
#     idx_r1::Int
#     idx_p1::Int
#     idx_p2::Int

#     # reaction rate index
#     idx_k::Int
# end

# struct Photodissociation13 <: Photodissociation
#     idx_r1::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int

#     # reaction rate index
#     idx_k::Int
# end

# struct Photodissociation14 <: Photodissociation
#     idx_r1::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int
#     idx_p4::Int

#     # reaction rate index
#     idx_k::Int
# end

# struct Photodissociation15 <: Photodissociation
#     idx_r1::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int
#     idx_p4::Int
#     idx_p5::Int

#     # reaction rate index
#     idx_k::Int
# end





# struct Collision11 <: CollisionReaction
#     idx_r1::Int
#     idx_p1::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision12 <: CollisionReaction
#     idx_r1::Int
#     idx_p1::Int
#     idx_p2::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision13 <: CollisionReaction
#     idx_r1::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision14 <: CollisionReaction
#     idx_r1::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int
#     idx_p4::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end


# struct Collision20 <: CollisionReaction
#     idx_r1::Int
#     idx_r2::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision21 <: CollisionReaction
#     idx_r1::Int
#     idx_r2::Int
#     idx_p1::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision22 <: CollisionReaction
#     idx_r1::Int
#     idx_r2::Int
#     idx_p1::Int
#     idx_p2::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision23 <: CollisionReaction
#     idx_r1::Int
#     idx_r2::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision24 <: CollisionReaction
#     idx_r1::Int
#     idx_r2::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int
#     idx_p4::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end

# struct Collision25 <: CollisionReaction
#     idx_r1::Int
#     idx_r2::Int
#     idx_p1::Int
#     idx_p2::Int
#     idx_p3::Int
#     idx_p4::Int
#     idx_p5::Int

#     # reaction rate index
#     idx_k::Int

#     needs_ro2::Bool
# end




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




function reaction_rate_coeff(rxn::Reaction, time::Float64, ro2_ratio::Float64)
    return nothing
end


function rhs_update!(u, rxn::Reaction)
    return nothing
end
