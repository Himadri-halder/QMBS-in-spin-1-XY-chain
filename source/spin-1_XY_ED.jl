using IterTools
using LinearAlgebra
using SparseArrays
using OrderedCollections
using Combinatorics
# using Arpack
using Random
using Polynomials
using LsqFit
using Base.Iterators:product
using Statistics
using JLD2,FileIO
using Printf
using Base.Threads

include("functions.jl")
include("lattice.jl")
include("basis.jl")
include("bimag_nematicneel_state.jl")
include("second_NN_states.jl")
include("spin-1_XY_hamiltonian.jl")
include("perturbative_corrections.jl")
include("spin-spin_correlation_matrix.jl")
include("spin_nematic_director_matrix.jl")