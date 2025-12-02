using Revise
using TheoryOfGames
# using Test
# using HiGHS
using JuMP
using Ipopt
using YAML
using FileIO
using JLD2

include("tests.jl")


"Function to create EnumMode out of an example"
to_EnumMode(example) = EnumMode(example.player_set, example.utility)

"Function to create IterMode out of an example"
function to_IterMode(example)
    util_combs = utility_combs(example.player_set, example.utility)
    keys_no_grand_coalition = setdiff(keys(util_combs), [Set(), Set(example.player_set)])

    let  player_set=example.player_set, util_combs=util_combs, keys_no_grand_coalition=keys_no_grand_coalition
        
        callback_benefit_by_coalition = (coal)->util_combs[Set(coal)]

        function callback_worst_coalition(profit_dist; kwargs...)
            min_surplus_combs = Dict(
                comb=>sum(Float64[profit_dist[c] for c in comb]) - util_combs[Set(comb)]
                for comb in keys_no_grand_coalition
            )
            min_surplus, least_benefit_coal = findmin(min_surplus_combs)
            return [(
                least_profitable_coalition_status=JuMP.Containers.DenseAxisArray(
                    [pl in least_benefit_coal ? 1.0 : 0.0 for pl in player_set],
                    player_set,
                ),
                least_profitable_coalition=least_benefit_coal,
                coalition_benefit=util_combs[Set(least_benefit_coal)],
                min_surplus=min_surplus,
            )]
        end

        return TheoryOfGames.IterMode(example.player_set, callback_benefit_by_coalition, callback_worst_coalition)
    end
end

optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0)
# optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0)  #, "tol"=>1e-4)  #, "NonConvex"=>2)

example = Examples.three_users_mapping

player_set = example.player_set
utility = example.utility

enum_mode = TheoryOfGames.EnumMode(example.player_set, example.utility)

# save("test.jld2", enum_mode)

# loaded_enum = load("test.jld2", EnumMode())

# val_enum = in_core(enum_mode, optimizer)

# uc = utility_combs(player_set, utility)

# sh = shapley_value(player_set, utility)

# ref_dist = Dict(
#     zip(example.player_set, fill(0.0, length(example.player_set)))
# )

# mode = EnumMode(player_set, utility)

# a = ref_in_core(mode, ref_dist, optimizer)

mode_iter = to_IterMode(example)
# mode_enum = to_EnumMode(example)

# nucleolus(mode_enum, optimizer)

profit_distribution, min_surplus, history = least_core(mode_iter, optimizer, raw_outputs=true)

#specific_least_core_dist_enum, min_surplus_enum, model_dist_enum = var_least_core(mode_enum, optimizer; raw_outputs=true)
#profit_distribution, min_surplus, history, model_dist = in_core(mode_iter, optimizer; raw_outputs=true)
# verify_in_core(profit_distribution, mode_iter, optimizer)


# profit_distribution, min_surplus, history = least_core(mode_example, optimizer, raw_outputs=true)
#[9.249999967500331, 5.500000064999378, 9.24999996750029]

# p_test = JuMP.Containers.DenseAxisArray(
#     [8,8,8],
#     [1,2,3]
# )

# in_core_check = verify_in_core(profit_distribution, enum_mode, optimizer)