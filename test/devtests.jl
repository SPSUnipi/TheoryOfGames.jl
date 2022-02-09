using Games
# using Test
using Gurobi
using JuMP



utility(x) = 1 in x ? 1.0 : 0.0
player_set = [1, 2, 3]

utility_combs(player_set, utility)

shapley_value(player_set, utility)

a = least_core(player_set, utility, Gurobi.Optimizer)