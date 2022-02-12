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
function shapley_value(player_set, utility::Function, mode::EnumMode=EnumMode(); verbose=true)
    
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
function least_core(player_set, utility::Function, optimizer, mode::EnumMode=EnumMode(); verbose=true)
    
    # get the combinations of utilities for every coalition
    utilities = utility_combs(player_set, utility, verbose=verbose)

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
nucleolus_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function nucleolus(
        player_set,
        utility::Function,
        optimizer,
        mode::EnumMode=EnumMode();
        verbose=true,
        tol=1e-5
    )
    
    # get the combinations of utilities for every coalition
    utilities = utility_combs(player_set, utility, verbose=verbose)

    # total combinations
    comb_set = combinations(player_set)

    # initialize JuMP model

    # total combinations
    comb_set = combinations(player_set)
    n_coalitions = number_coalitions(player_set)

    # Definition of JuMP model
    model_dist = Model(optimizer)

    @variable(model_dist, minimum(values(utilities)) <= profit_dist[pl in player_set] <= maximum(values(utilities)))

    # the gain of the worst group of each iteration
    @variable(model_dist, minimum(values(utilities)) <= delta_worst <= maximum(values(utilities)))

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == utilities[Set(player_set)])

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in comb_set; length(comb) > 0],
        sum([profit_dist[pl] for pl in comb]) >= utilities[Set(comb)] + delta_worst
    )

    # specify the objective
    @objective(model_dist, Max, delta_worst)

    # if option verbose setup progresss bar
    pbar = verbose ? ProgressBar(1:n_coalitions) : nothing

    fixed_constraints = 1  # initialize constraint to the null value

    while fixed_constraints < n_coalitions

        # optimize the iteration
        optimize!(model_dist)

        # calculate the value of the current delta_worst
        delta_worst_val = value(delta_worst)

        # get dual variable of con_subset_profit
        dual_val = dual.(con_subset_profit)

        # when the variable is non-zero means that the constraint is binding
        # thus update the constraint to specify the corresponding value of delta_worst
        # as limit
        for comb in comb_set
            if dual_val[comb] > tol
                # get current rhs of the constraint
                pre_rhs = normalized_rhs(con_subset_profit[comb])

                # disable delta_worst for the specific constraint
                set_normalized_coefficient(con_subset_profit[comb], delta_worst, 0.0)

                # update rhs
                set_normalized_rhs(con_subset_profit[comb], pre_rhs + delta_worst_val)

                # update number of fixed constraints
                fixed_constraints += 1
            end
        end

        # if option verbose update progress bar
        if verbose
            pbar.current = fixed_constraints
            ProgressBars.display_progress(pbar)
        end
    end

    nuc_dist = Dict(zip(player_set, value.(profit_dist).data))

    return nuc_dist
end


"""
    specific_in_core(player_set, utility, dist_objective, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the core
and maximizes a specific distribution objective specified by dist_objective
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility : Function
    Function that describes the utility function for every coalition in player_set
dist_objective : Function
    Function to build the objective function of the profit distribution
    It shall a function with a single argument, which is the JuMP model and
    builds the objective functions using the JuMP @NLobjective or @objective commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function specific_in_core(
        player_set,
        utility::Function,
        dist_objective::Function,
        optimizer,
        mode::EnumMode=EnumMode(); 
        verbose=true,
        tol=1e-5
    )
    
    # get the combinations of utilities for every coalition
    utilities = utility_combs(player_set, utility, verbose=verbose)

    # total combinations
    comb_set = combinations(player_set)

    # initialize JuMP model
    model_dist = Model(optimizer)

    @variable(model_dist, minimum(values(utilities)) <= profit_dist[u in player_set] <= maximum(values(utilities)))

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == utilities[Set(player_set)])

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in comb_set; length(comb) > 0],
        sum([profit_dist[pl] for pl in comb]) >= utilities[Set(comb)] + 0.0 # the 0.0 is to remind 
                                            # that delta_worst shall be zero to belong to the core
    )

    # The objective function is calculated using the dist_objective function
    dist_objective(model_dist)

    # optimize
    optimize!(model_dist)

    lc_dist = Dict(zip(player_set, value.(profit_dist).data))

    return lc_dist
end


"""
    in_core(player_set, utility, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the core
for profit distribution for a game described by the utility function
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
in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function in_core(player_set, utility::Function, optimizer, mode::EnumMode=EnumMode(); verbose=true, tol=1e-5)
    
    # create feasibility problem
    in_core_objective(m) = @objective(m, Max, 0.0)

    return specific_in_core(
        player_set,
        utility,
        in_core_objective,  # the core is only a feasibility problem, thus specify constant obj. funciton
        optimizer,
        mode; 
        verbose=verbose,
        tol=tol
    )
end


"""
    var_core(player_set, utility, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the core
and minimizes the variance of the profit allocation among the plauers
for a game described by the utility function and the grand coalition of player_set.

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
var_core : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function var_core(player_set, utility::Function, optimizer, mode::EnumMode=EnumMode(); verbose=true, tol=1e-5)

    # create the objective function for the problem
    function var_objective(m)
        obj = sum(m[:profit_dist].^2)/length(m[:profit_dist]) - sum(m[:profit_dist]/length(m[:profit_dist]))^2
        println(typeof(obj))
        @objective(m,
        Min, 
        obj)
    end

    return specific_in_core(
        player_set,
        utility,
        var_objective,
        optimizer,
        mode; 
        verbose=verbose,
        tol=tol
    )
end
