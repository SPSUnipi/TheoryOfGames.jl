using Revise
using Games
# using Test
# using HiGHS
using JuMP
using Ipopt
#using GLPK

include("tests.jl")



function to_RobustMode(example)
    util_combs = utility_combs(example.player_set, example.utility)

    let util_combs=util_combs
        
        callback_benefit_by_coalition = (coal)->util_combs[Set(coal)]

        function callback_worst_coalition(profit_dist)
            min_surplus_combs = Dict(
                comb=>sum(Float64[profit_dist[c] for c in comb]) - util_combs[Set(comb)]
                for comb in keys(util_combs)
            )
            min_surplus, least_benefit_coal = findmin(min_surplus_combs)
            return least_benefit_coal, util_combs[Set(least_benefit_coal)], min_surplus
        end

        return Games.RobustMode(example.player_set, callback_benefit_by_coalition, callback_worst_coalition)
    end
end

optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0)
# optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0)  #, "tol"=>1e-4)  #, "NonConvex"=>2)

example = Examples.three_users_mapping

player_set = example.player_set
utility = example.utility

enum_mode = Games.EnumMode(example.player_set, example.utility)


val_enum = in_core(enum_mode, optimizer)

# uc = utility_combs(player_set, utility)

# sh = shapley_value(player_set, utility)

# ref_dist = Dict(
#     zip(example.player_set, fill(0.0, length(example.player_set)))
# )

# mode = EnumMode(player_set, utility)

# a = ref_in_core(mode, ref_dist, optimizer)

mode_example = to_RobustMode(example)

profit_distribution, min_surplus, history = least_core(mode_example, optimizer)


p_test = JuMP.Containers.DenseAxisArray(
    [8,8,8],
    [1,2,3]
)

in_core_check = verify_in_core(profit_distribution, enum_mode, optimizer)