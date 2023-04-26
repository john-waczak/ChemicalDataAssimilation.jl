abstract type Reaction end


# optionally, we could further reduce this to abstract types for these and
# then define concrete types based on the number of products

# abstract type Photodissociation <: Reaction end
# abstract type Bimolecular <: Reaction end
# abstract type Termolecular <: Reaction end

# struct Bimolecular22  a struct for a bimolecular reaction with 2 reactants and 2 products
# end





struct Photodissociation <: Reaction
    rate_coeff::Function
end


struct Bimolecular <: Reaction
    rate_coeff::Function
    needs_ro2::Bool
end


struct Termolecular <: Reaction
    rate_coeff::Function
    needs_ro2::Bool
end
