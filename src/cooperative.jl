"""
    shapley_value(player_set, utilities, mode; verbose)

Function to calculte the shapley value for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary
    AbstractDictionary storing the utilities for every coalition.
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
shapley_value : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function shapley_value(
        player_set,
        utilities,
        mode::EnumMode=EnumMode();
        verbose=true,
    )

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
    shapley_value(player_set, utility_func, mode; verbose)

Function to calculte the shapley value for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility_func : Function
    Function to derive the utility by each coalition
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
shapley_value : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function shapley_value(
        player_set,
        utility_func::Function,
        mode::EnumMode=EnumMode();
        verbose=true,
    )

    utilities = utility_combs(player_set, utility_func, verbose=verbose)

    return shapley_value(player_set, utilities, mode, verbose=berbose)
end


"""
    least_core(player_set, utilities, optimizer, mode; verbose)

Function to calculte the least core profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary
    AbstractDictionary storing the utilities for every coalition.
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
leastcore_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function least_core(
        player_set,
        utilities,
        optimizer,
        mode::EnumMode=EnumMode();
        verbose=true
    )

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
    least_core(player_set, utility_func, optimizer, mode; verbose)

Function to calculte the least core profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility_func : Function
    Function to derive the utility by each coalition
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
leastcore_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function least_core(
        player_set,
        utility_func::Function,
        optimizer,
        mode::EnumMode=EnumMode();
        verbose=true
    )

    utilities = utility_combs(player_set, utility_func, verbose=verbose)

    return least_core(player_set, utilities, optimizer, mode, verbose=berbose)
end



"""
    nucleolus(player_set, utilities, mode; verbose, tol)

Function to calculte the nucleolus profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary
    AbstractDictionary storing the utilities for every coalition.
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
tol (optional, default 1e-5)
    Accepted tolerance for identifying the active constraints in the optimization
    procedure

Outputs
------
nucleolus_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function nucleolus(
        player_set,
        utilities,
        optimizer,
        mode::EnumMode=EnumMode();
        verbose=true,
        tol=1e-5
    )

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

    # optimize the iteration
    optimize!(model_dist)

    # last visited combination, to avoid repetitions
    last_visited_combs = [empty_set(player_set)]
    visited_nodes = []

    while (fixed_constraints < n_coalitions)

        # calculate the value of the current delta_worst
        delta_worst_val = value(delta_worst)

        # get dual variable of con_subset_profit
        dual_val = dual.(con_subset_profit)

        # set of current visited nodes to avoid infinite loops
        visited_combs = []

        # when the variable is non-zero means that the constraint is binding
        # thus update the constraint to specify the corresponding value of delta_worst
        # as limit
        for comb in comb_set
            if dual_val[comb] > tol

                # update visited combs
                push!(visited_combs, comb)

                if comb ∉ last_visited_combs
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
        end

        # optimize the iteration
        optimize!(model_dist)

        if Set(last_visited_combs) == Set(visited_nodes)
            # if previous and visited combs coincide, then
            # stop the iterations
            fixed_constraints = n_coalitions
            break
        else
            last_visited_combs = visited_nodes
        end

        # if option verbose update progress bar
        if verbose
            pbar.current = fixed_constraints
            ProgressBars.display_progress(pbar)
        end
    end

    # if option verbose update progress bar
    if verbose
        pbar.current = fixed_constraints
        ProgressBars.display_progress(pbar)
    end

    nuc_dist = Dict(zip(player_set, value.(profit_dist).data))

    return nuc_dist
end



"""
    nucleolus(player_set, utility_func, mode; verbose, tol)

Function to calculte the nucleolus profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility_func : Function
    Function to derive the utility by each coalition
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
tol (optional, default 1e-5)
    Accepted tolerance for identifying the active constraints in the optimization
    procedure

Outputs
------
nucleolus_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function nucleolus(
        player_set,
        utility_func::Function,
        optimizer,
        mode::EnumMode=EnumMode();
        verbose=true,
        tol=1e-5
    )

    utilities = utility_combs(player_set, utility_func, verbose=verbose)

    return nucleolus(player_set, utilities, optimizer, mode, verbose=verbose, tol=tol)
end


"""
    specific_in_core(player_set, utilities, dist_objective, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the core
and maximizes a specific distribution objective specified by dist_objective
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary
    AbstractDictionary storing the utilities for every coalition.
dist_objective : Function
    Function to build the objective function of the profit distribution
    It shall a function with two arguments, which are the JuMP model and the player set,
    and it shall build the objective functions using the JuMP @NLobjective or @objective
    commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function specific_in_core(
        player_set,
        utilities,
        dist_objective::Function,
        optimizer,
        mode::EnumMode=EnumMode();
        verbose=true,
    )

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
    dist_objective(model_dist, player_set)

    # optimize
    optimize!(model_dist)

    specific_in_core_dist = Dict(zip(player_set, value.(profit_dist).data))

    return specific_in_core_dist
end


"""
    specific_in_core(player_set, utilities, dist_objective, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the core
and maximizes a specific distribution objective specified by dist_objective
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility_func : Function
    Function to derive the utility by each coalition
dist_objective : Function
    Function to build the objective function of the profit distribution
    It shall a function with two arguments, which are the JuMP model and the player set,
    and it shall build the objective functions using the JuMP @NLobjective or @objective
    commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function specific_in_core(
        player_set,
        utility_func::Function,
        dist_objective::Function,
        optimizer,
        mode::EnumMode=EnumMode();
        verbose=true,
    )

    utilities = utility_combs(player_set, utility_func, verbose=verbose)

    return specific_in_core(player_set, utilities, dist_objective, optimizer, mode, verbose=verbose)
end


"""
    in_core(player_set, utility, optimizer, mode; verbose, utilities)

Function to calculte a stable profit distribution that belongs to the core
for profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary or Function
    AbstractDictionary storing the utilities for every coalition.
    or Function to derive the utility by each coalition
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function in_core(player_set, utilities, optimizer, mode::EnumMode=EnumMode(); kwargs...)
    
    # create feasibility problem
    in_core_objective(m, player_set) = @objective(m, Max, 0.0)

    return specific_in_core(
        player_set,
        utilities,
        in_core_objective,  # the core is only a feasibility problem, thus specify constant obj. funciton
        optimizer,
        mode; 
        kwargs...
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
utilities : AbstractDictionary or Function
    AbstractDictionary storing the utilities for every coalition.
    or Function to derive the utility by each coalition
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
var_core : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function var_core(player_set, utilities, optimizer, mode::EnumMode=EnumMode(); kwargs...)

    # create the objective function for the problem
    function var_objective(m, player_set)
        obj = sum(m[:profit_dist].^2)/length(m[:profit_dist]) - sum(m[:profit_dist]/length(m[:profit_dist]))^2
        @objective(m, Min, obj)
    end

    return specific_in_core(
        player_set,
        utilities,
        var_objective,
        optimizer,
        mode; 
        kwargs...
    )
end


"""
    ref_in_core(player_set, utilities, ref_dist, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the core
and minimizes the variance of the profit allocation among the players
with respect to a pre-defined reward distribution function
for a game described by the utility function and the grand coalition of player_set.


Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary or Function
    AbstractDictionary storing the utilities for every coalition.
    or Function to derive the utility by each coalition
ref_dist : AbstrctDict
    Reference distribution by player
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
norm : Any
    Normalization denominator for every player
    Default value nothing, hence for every player the Normalization
    factor is 1.0
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function ref_in_core(
        player_set,
        utilities,
        ref_dist,
        optimizer,
        mode::EnumMode=EnumMode();
        kwargs...
    )

    # create the objective function for the problem
    function var_objective(m, player_set)
        n_players = length(player_set)

        # initialize dictionary for the normalization denominator
        norm_den = isnothing(norm) ? Dict(zip(player_set, fill(1.0, n_players))) : norm

        obj = sum(
            ((m[:profit_dist][p] - ref_dist[p])/norm_den[p]).^2
            for p in player_set
        )
        @objective(m, Min, obj)
    end

    return specific_in_core(
        player_set,
        utilities,
        var_objective,
        optimizer,
        mode; 
        kwargs...
    )
end


"""
    specific_least_core(player_set, utilities, dist_objective, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes a specific objective for the profit allocation among the plauers,
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary
    AbstractDictionary storing the utilities for every coalition
dist_objective : Function
    Function to build the objective function of the profit distribution, after the least core
    is obtained
    It shall be a function with two arguments: the JuMP model of the problem and the player set,
    and the function shall build the desired objective functions using the 
    JuMP @NLobjective or @objective commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
utilities (optional, default nothing)

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function specific_least_core(
        player_set,
        utilities,
        dist_objective::Function,
        optimizer,
        mode::EnumMode=EnumMode(); 
        verbose=true,
    )

    # total combinations
    comb_set = combinations(player_set)

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

    # optimize and get the value of delta_worst
    optimize!(model_dist)

    # fix the value of delta_worst
    fix(delta_worst, value(delta_worst), force=true)

    # Update objective function according to the desired dist_objective function
    dist_objective(model_dist, player_set)

    # re-optimize to get the result according with dist_objective
    optimize!(model_dist)

    specific_least_core_dist = Dict(zip(player_set, value.(profit_dist).data))

    return specific_least_core_dist
end


"""
    specific_least_core(player_set, utility_func, dist_objective, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes a specific objective for the profit allocation among the plauers,
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
player_set : Vector
    Vector of the players of the coalition
utility_func : Function
    Function to derive the utility by each coalition
dist_objective : Function
    Function to build the objective function of the profit distribution, after the least core
    is obtained
    It shall be a function with two arguments: the JuMP model of the problem and the player set,
    and the function shall build the desired objective functions using the 
    JuMP @NLobjective or @objective commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
utilities (optional, default nothing)

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function specific_least_core(
        player_set,
        utility_func::Function,
        dist_objective::Function,
        optimizer,
        mode::EnumMode=EnumMode(); 
        verbose=true,
    )

    utilities = utility_combs(player_set, utility_func, verbose=verbose)

    return specific_least_core(player_set, utilities, dist_objective, optimizer, mode, verbose=verbose)
end


"""
    var_least_core(player_set, utilities, optimizer, mode; verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes the variance of the profit allocation among the players
for a game described by the utility function and the grand coalition of player_set.


Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary or Function
    AbstractDictionary storing the utilities for every coalition.
    or Function to derive the utility by each coalition
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function var_least_core(
        player_set,
        utilities,
        optimizer,
        mode::EnumMode=EnumMode(); 
        kwargs...
    )

    # create the objective function for the problem
    function var_objective(m, player_set)
        obj = sum(m[:profit_dist].^2)/length(m[:profit_dist]) - sum(m[:profit_dist]/length(m[:profit_dist]))^2
        @objective(m, Min, obj)
    end

    return specific_least_core(
        player_set,
        utilities,
        var_objective,
        optimizer,
        mode; 
        kwargs...
    )
end


"""
    ref_least_core(player_set, utilities, ref_dist, optimizer, mode; norm, verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes the variance of the profit allocation among the players
with respect to a pre-defined reward distribution function
for a game described by the utility function and the grand coalition of player_set.


Inputs
------
player_set : Vector
    Vector of the players of the coalition
    utilities : AbstractDictionary or Function
    AbstractDictionary storing the utilities for every coalition.
    or Function to derive the utility by each coalition.
ref_dist : AbstrctDict
    Reference distribution by player
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
norm (optional, default nothing)
    Normalization denominator for every player
    Default value nothing, hence for every player the Normalization
    factor is 1.0
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function ref_least_core(
        player_set,
        utilities,
        ref_dist,
        optimizer,
        mode::EnumMode=EnumMode();
        norm=nothing,
        kwargs...
    )

    # create the objective function for the problem
    function var_objective(m, player_set)
        n_players = length(player_set)

        # initialize dictionary for the normalization denominator
        norm_den = isnothing(norm) ? Dict(zip(player_set, fill(1.0, n_players))) : norm

        obj = sum(
            ((m[:profit_dist][p] - ref_dist[p])/norm_den[p]).^2
            for p in player_set
        )
        @objective(m, Min, obj)
    end

    return specific_least_core(
        player_set,
        utilities,
        var_objective,
        optimizer,
        mode;
        kwargs...
    )
end




"""
    ref_least_core(player_set, utilities, ref_dist, optimizer, mode; norm, verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes the variance of the profit allocation among the players
with respect to a pre-defined reward distribution function
for a game described by the utility function and the grand coalition of player_set.


Inputs
------
player_set : Vector
    Vector of the players of the coalition
utilities : AbstractDictionary or Function
    AbstractDictionary storing the utilities for every coalition.
    or Function to derive the utility by each coalition.
ref_dist : AbstrctDict
    Reference distribution by player
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
mode : AbstractCalcMode
    Calculation mode: enumerative technique
norm (optional, default nothing)
    Normalization denominator for every player
    Default value nothing, hence for every player the Normalization
    factor is 1.0
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function ref_least_core(
        player_set,
        utilities,
        ref_dist,
        optimizer,
        mode::EnumMode=EnumMode();
        norm=nothing,
        kwargs...
    )

    # create the objective function for the problem
    function var_objective(m, player_set)
        n_players = length(player_set)

        # initialize dictionary for the normalization denominator
        norm_den = isnothing(norm) ? Dict(zip(player_set, fill(1.0, n_players))) : norm

        obj = sum(
            ((m[:profit_dist][p] - ref_dist[p])/norm_den[p]).^2
            for p in player_set
        )
        @objective(m, Min, obj)
    end

    return specific_least_core(
        player_set,
        utilities,
        var_objective,
        optimizer,
        mode;
        kwargs...
    )
end