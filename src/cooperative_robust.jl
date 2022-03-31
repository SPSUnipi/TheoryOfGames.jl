"""
least_core(mode, player_set, utilities, optimizer; verbose)

Function to calculte the least core profit distribution for a game described by the utility function
and the grand coalition of player_set.
This function implements a robust-like approach by using a column generating approach
to iteratively add constraints to the master problem.
To do so, the callback function stored in the mode category is exploited

Inputs
------
mode : RobustMode
    Calculation mode that contains the reference to the callback function
    that given a profit distribution scheme, returns a tuple of
    the set of the coalition with the worst profit and its total benefit to be shared
    callback_worst_coalition accepts one argument (current profit sharing)
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
leastcore_dist : Dict
Dictionary of the fair distributions of the profits among the players
"""
function least_core(
        mode::RobustMode,
        optimizer;
        verbose=true
    )

    callback_worst_coalition = mode.callback_worst_coalition

    # total combinations
    comb_set = combinations(player_set)

    # get empty coalition
    empty_coal = empty_set(player_set)
    empty_val = utilities[empty_coal]

    # initialize JuMP model
    model_dist = Model(optimizer)

    @variable(model_dist, minimum(values(utilities)) <= profit_dist[u in player_set] <= maximum(values(utilities)))

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == utilities[Set(player_set)])

    # the gain of the worst group of the current iteration
    @variable(model_dist, minimum(values(utilities)) <= delta_worst <= maximum(values(utilities)))

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in comb_set; length(comb) > 0],
        sum([profit_dist[pl] for pl in comb]) >= utilities[Set(comb)] + delta_worst
    )

    # specify the objective maximize the benefit or remaining in the coalition
    @objective(model_dist, Max, delta_worst)

    # optimize
    optimize!(model_dist)

    lc_dist = Dict(zip(player_set, value.(profit_dist).data))

    return lc_dist
end