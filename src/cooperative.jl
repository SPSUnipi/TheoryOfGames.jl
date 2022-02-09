"""
    shapley_value(player_set, utility, mode; verbose)

Function to calculte the shapley value for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility : Function
    Function that describes the utility function for every coalition in player_set
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool
    When true, it shows a progress bar to describe the current execution status

Outputs
------
shapley_value : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function shapley_value(player_set::Vector, utility::Function, mode::EnumMode=EnumMode(); verbose=true)
    
    # get the combinations of utilities for every coalition
    utilities = utility_combs(player_set, utility, verbose=verbose)

    # initialize shapley dictionary
    type_output = eltype(values(utilities))
    sh_val = Dict{eltype(player_set), type_output}()

    n_pl = length(player_set)

    empty_coal = empty_set(player_set)

    for pl in player_set  # see https://en.wikipedia.org/wiki/Shapley_value
        sh_val[pl] = 1/n_pl*(
            sum(type_output[
                (utilities[union(Set(comb), pl)] - utilities[Set(comb)])/ binomial(n_pl-1, length(comb))
                for comb in combinations(setdiff(player_set, pl))
            ])  # general case except empty set
            + 
            (utilities[union(empty_coal, pl)] - utilities[empty_coal]) / binomial(n_pl-1, 0)  # case comb == empty_set
        )
    end

    return sh_val
end


"""
    least_core(player_set, utility, mode, optimizer; verbose)

Function to calculte the least core profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility : Function
    Function that describes the utility function for every coalition in player_set
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool
    When true, it shows a progress bar to describe the current execution status

Outputs
------
leastcore_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function least_core(player_set::Vector, utility::Function, optimizer, mode::EnumMode=EnumMode(); verbose=true)
    
    # get the combinations of utilities for every coalition
    utilities = utility_combs(player_set, utility, verbose=verbose)

    # total combinations
    comb_set = combinations(player_set)
    
    # get empty coalition
    empty_coal = empty_set(player_set)
    empty_val = utilities[empty_coal]
    type_val = empty_val

    # initialize JuMP model
    model_dist = Model(optimizer)

    @variable(model_dist, zero(empty_val) <= profit_dist[u in player_set] <= utilities[Set(player_set)])

    # the gain of the worst group of the current iteration
    @variable(model_dist, zero(empty_val) <= delta_worst <= utilities[Set(player_set)])

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in Set.(comb_set)],
        sum(type_val[profit_dist[pl] for pl in comb]) >= utilities[comb] + delta_worst
    )

    # specify the objective maximize the benefit or remaining in the coalition
    @objective(model_dist, Max, delta_worst)

    # optimize
    optimize!(model_dist)

    lc_dist = Dict(zip(player_set, value.(profit_dist).data))

    return lc_dist
end


"""
    nucleolus(player_set, utility, mode; verbose)

Function to calculte the nucleolus profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility : Function
    Function that describes the utility function for every coalition in player_set
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool
    When true, it shows a progress bar to describe the current execution status

Outputs
------
leastcore_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function nucleolus(player_set::Vector, utility::Function, optimizer, mode::EnumMode=EnumMode(); verbose=true)
    
    # get the combinations of utilities for every coalition
    utilities = utility_combs(player_set, utility, verbose=verbose)

    # total combinations
    comb_set = combinations(player_set)

    # get empty coalition
    empty_coal = empty_set(player_set)
    empty_val = utilities[empty_coalition]

    # initialize JuMP model
    # combinations with a fixed value of the variable
    fixed_combs = Dict{typeof(empty_coal), typeof(empty_val)}(empty_coal=>empty_val)

    # total combinations
    comb_set = combinations(user_agg_set)
    n_coalitions = number_coalitions(player_set)

    # Definition of JuMP model
    model_dist = Model(optimizer)

    @variable(model_dist, zero(empty_val) <= profit_dist[pl in player_set] <= utilities[Set(player_set)])

    # the gain of the worst group of each iteration
    @variable(model_dist, delta_worst >= zero(empty_val))

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in Set.(comb_set)],
        sum(type_val[profit_dist[pl] for pl in comb]) >= utilities[Set(comb)] + (
            (comb in keys(fixed_combs)) ? fixed_combs[comb] : delta_worst
        )
            
    )

    # specify the objective
    @objective(model_dist, Max, delta_worst)

    while length(fixed_combs) != n_coalitions

        ## Model definition

        

        # optimize the iteration
        optimize!(model_dist)
        #calculate the value of the current delta_worst
        delta_worst_val = value(delta_worst)

        excluded_combs = get_most_unhappy_coals(
            user_agg_set, value.(NPV_mod), calc_ref, shapley_agg_dict, delta_worst_val
        )

        for ex_comb in excluded_combs
            if !(ex_comb in keys(fixed_combs)) && (abs(dual(con_subset_profit[ex_comb])) > 1e-6)
                fixed_combs[ex_comb] = delta_worst_val
                # printfmtln("add {}", ex_comb)
            end
        end

        # printfmtln("delta: {:.2f}, n_combs: {:d}", delta_worst_val, length(excluded_combs))
    end

    return lc_dist
end
