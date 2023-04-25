module ChemicalDataAssimilation

# Write your package code here.
include("parse_fac.jl")
include("reaction_rates.jl")

export read_fac_file, parse_rxns
export get_rrate_list, rc_to_expression, generate_rrates
end
