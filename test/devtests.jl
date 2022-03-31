using Revise
using Games
# using Test
# using HiGHS
using JuMP
using Gurobi
#using Ipopt
#using GLPK

include("tests.jl")

optimizer = optimizer_with_attributes(Gurobi.Optimizer)
# optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0)  #, "tol"=>1e-4)  #, "NonConvex"=>2)

example = Examples.three_users_mapping

player_set = example.player_set
utility = example.utility

# uc = utility_combs(player_set, utility)

# sh = shapley_value(player_set, utility)

ref_dist = Dict(
    zip(example.player_set, fill(0.0, length(example.player_set)))
)

mode = EnumMode(player_set, utility)

a = ref_in_core(mode, ref_dist, optimizer)