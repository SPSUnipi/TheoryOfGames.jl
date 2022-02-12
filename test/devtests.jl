using Games
# using Test
using HiGHS
using JuMP
using Gurobi
using Ipopt
using GLPK

include("tests.jl")

optimizer = optimizer_with_attributes(HiGHS.Optimizer)
# optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0)  #, "tol"=>1e-4)  #, "NonConvex"=>2)

example = Examples.three_users_mapping

player_set = example.player_set
utility = example.utility

# uc = utility_combs(player_set, utility)

# sh = shapley_value(player_set, utility)

a = nucleolus(player_set, utility, optimizer)