module ChemicalDataAssimilation

using ProgressMeter
using CSV, DataFrames
using DelimitedFiles
using SparseArrays
using LinearAlgebra

include("parse_fac.jl")  # reading fac file
include("reaction_rates.jl")  # generate time series of generic/complex rate coeffs
include("species.jl")  # generate species list and name lookup
include("reference_measurements.jl")  #  generate time series of measurements/num densities
include("photolysis_rates.jl")  # generate time series of photolysis rates
include("reaction_structs.jl")  # define structs for storing reactions
include("initialize.jl")  # define method to sanely initialize concentrations
include("ro2.jl")  # define indices for ro2 sum
include("rrates_mechanism.jl")  # generate time series for actual reaction rate coeffs
include("stoich_mats.jl")
include("derivative_structs.jl")
include("jacobian_struct.jl")
include("4dvar.jl")
include("EKF.jl")

#---
export read_fac_file, parse_rxns, get_spec_idx
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
export parse_rxn, generate_reaction_list
#---
export rxnrate_to_expression
#---
export generate_init_dict
#---
export generate_ro2
#---
export generate_rrates_mechanism
#---
export generate_stoich_mat, get_sparse_mat
#---
export DerivativeTerms, update_derivative!, write_rhs_func
#---
export JacobianTerms, update_jacobian!, write_jac_func, generate_jac_prototype
#---
export ObsOpMeas, Rmat, Rinv
#---
export JObs
end
