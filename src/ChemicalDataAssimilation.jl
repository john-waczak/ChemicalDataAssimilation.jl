module ChemicalDataAssimilation

using ProgressMeter
using CSV, DataFrames
using DelimitedFiles
include("parse_fac.jl")  # reading fac file
include("reaction_rates.jl")  # generate time series of generic/complex rate coeffs
include("species.jl")  # generate species list and name lookup
include("reference_measurements.jl")  #  generate time series of measurements/num densities
include("photolysis_rates.jl")  # generate time series of photolysis rates
include("reaction_structs.jl")  # define structs for storing reactions
include("initialize.jl")  # define method to sanely initialize concentrations
include("ro2.jl")  # define indices for ro2 sum
include("rrates_mechanism.jl")  # generate time series for actual reaction rate coeffs


#---
export read_fac_file, parse_rxns
#---
export get_rrate_list, rc_to_expression, generate_rrates_funcs, generate_rrates
#---
export generate_name_conversion_table, generate_species, generate_species_df
#---
export generate_densities
#---
export generate_photolysis_rates
#---
export Reaction, Photodissociation, CollisionReaction
export Photodissociation11, Photodissociation12
export Collision11, Collision20, Collision21, Collision22
#---
export rxnrate_to_expression
#---
export generate_init_dict
#---
export generate_ro2
#---
export generate_rrates_mechanism



end
