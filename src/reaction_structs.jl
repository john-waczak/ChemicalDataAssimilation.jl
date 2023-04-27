# Define Reaction Type Hierarchy
abstract type Reaction end
abstract type Photodissociation <: Reaction end
abstract type CollisionReaction <: Reaction end


# Define Concrete Types
struct Photodissociation11
    idx_r1::Int
    idx_p1::Int

    # reaction rate index
    idx_k::Int
end

struct Photodissociation12
    idx_r1::Int
    idx_p1::Int
    idx_p2::Int

    # reaction rate index
    idx_k::Int
end

struct Collision11
    idx_r1::Int
    idx_p1::Int

    # reaction rate index
    idx_k::Int

    needs_ro2::Bool
end

struct Collision20
    idx_r1::Int
    idx_r2::Int

    # reaction rate index
    idx_k::Int

    needs_ro2::Bool
end

struct Collision21
    idx_r1::Int
    idx_r2::Int
    idx_p1::Int

    # reaction rate index
    idx_k::Int

    needs_ro2::Bool
end

struct Collision22
    idx_r1::Int
    idx_r2::Int
    idx_p1::Int
    idx_p2::Int

    # reaction rate index
    idx_k::Int

    needs_ro2::Bool
end



function reaction_rate_coeff(rxn::Reaction, time::Float64, ro2_ratio::Float64)
    return nothing
end




function rhs_update!(u, rxn::Reaction)
    return nothing
end
