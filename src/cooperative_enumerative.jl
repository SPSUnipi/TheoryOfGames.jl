"""
    shapley_value(mode; verbose)

Function to calculte the shapley value for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
shapley_value : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function shapley_value(
        mode::EnumMode;
        verbose=true,
    )

    utilities = mode.utilities
    player_set = mode.player_set

    # initialize shapley dictionary
    type_output = eltype(values(utilities))
    sh_val = Dict{eltype(player_set), type_output}()

    n_pl = length(player_set)

    empty_coal = empty_set(player_set)

    for pl in player_set  # see https://en.wikipedia.org/wiki/Shapley_value
        sh_val[pl] = 1/n_pl*(
            sum(type_output[
                (utilities[union(Set(comb), [pl])] - utilities[Set(comb)])/ binomial(n_pl-1, length(comb))
                for comb in combinations(setdiff(player_set, [pl]))
            ])  # general case except empty set
            + 
            (utilities[union(empty_coal, [pl])] - utilities[empty_coal]) / binomial(n_pl-1, 0)  # case comb == empty_set
        )
    end

    return sh_val
end


"""
    least_core(mode, utilities, optimizer; verbose)

Function to calculte the least core profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
player_set : Vector
    Vector of the players of the coalition
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
        mode::EnumMode,
        optimizer;
        verbose=true
    )

    utilities = mode.utilities
    player_set = mode.player_set

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
    @variable(model_dist, minimum(values(utilities)) <= min_surplus <= maximum(values(utilities)))

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in comb_set; length(comb) > 0],
        sum([profit_dist[pl] for pl in comb]) >= utilities[Set(comb)] + min_surplus
    )

    # specify the objective maximize the benefit or remaining in the coalition
    @objective(model_dist, Max, min_surplus)

    # optimize
    optimize!(model_dist)

    lc_dist = Dict(zip(player_set, value.(profit_dist).data))

    return lc_dist
end



"""
    nucleolus(mode, utilities; verbose, tol)

Function to calculte the nucleolus profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
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
        mode::EnumMode,
        optimizer;
        verbose=true,
        tol=1e-5
    )

    utilities = mode.utilities
    player_set = mode.player_set

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
    @variable(model_dist, minimum(values(utilities)) <= min_surplus <= maximum(values(utilities)))

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == utilities[Set(player_set)])

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in comb_set; length(comb) > 0],
        sum([profit_dist[pl] for pl in comb]) >= utilities[Set(comb)] + min_surplus
    )

    # specify the objective
    @objective(model_dist, Max, min_surplus)

    # if option verbose setup progresss bar
    pbar = verbose ? ProgressBar(1:n_coalitions) : nothing

    fixed_constraints = 1  # initialize constraint to the null value

    # optimize the iteration
    optimize!(model_dist)

    # last visited combination, to avoid repetitions
    last_visited_combs = [empty_set(player_set)]
    visited_nodes = []

    while (fixed_constraints < n_coalitions)

        # calculate the value of the current min_surplus
        min_surplus_val = value(min_surplus)

        # get dual variable of con_subset_profit
        dual_val = dual.(con_subset_profit)

        # set of current visited nodes to avoid infinite loops
        visited_combs = []

        # when the variable is non-zero means that the constraint is binding
        # thus update the constraint to specify the corresponding value of min_surplus
        # as limit
        for comb in comb_set
            if dual_val[comb] > tol

                # update visited combs
                push!(visited_combs, comb)

                if comb ∉ last_visited_combs
                    # get current rhs of the constraint
                    pre_rhs = normalized_rhs(con_subset_profit[comb])

                    # disable min_surplus for the specific constraint
                    set_normalized_coefficient(con_subset_profit[comb], min_surplus, 0.0)

                    # update rhs
                    set_normalized_rhs(con_subset_profit[comb], pre_rhs + min_surplus_val)

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
    specific_in_core(mode, dist_objective, optimizer; verbose)

Function to calculte a stable profit distribution that belongs to the core
and maximizes a specific distribution objective specified by dist_objective
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
dist_objective : Function
    Function to build the objective function of the profit distribution
    It shall a function with two arguments, which are the JuMP model and the player set,
    and it shall build the objective functions using the JuMP @NLobjective or @objective
    commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function specific_in_core(
        mode::EnumMode,
        dist_objective::Function,
        optimizer;
        verbose=true,
    )

    utilities = mode.utilities
    player_set = mode.player_set

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
                                            # that min_surplus shall be zero to belong to the core
    )

    # The objective function is calculated using the dist_objective function
    dist_objective(model_dist, player_set)

    # optimize
    optimize!(model_dist)

    if termination_status(model_dist) in [JuMP.MOI.OPTIMAL, JuMP.MOI.LOCALLY_SOLVED]
        specific_in_core_dist = Dict(zip(player_set, value.(profit_dist).data))
    else
        specific_in_core_dist = nothing
    end

    return specific_in_core_dist
end


"""
    in_core(mode, optimizer; verbose, utilities)

Function to calculte a stable profit distribution that belongs to the core
for profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function in_core(mode::EnumMode, optimizer; kwargs...)
    
    # create feasibility problem
    in_core_objective(m, player_set) = @objective(m, Max, 0.0)

    return specific_in_core(
        mode,
        in_core_objective,  # the core is only a feasibility problem, thus specify constant obj. funciton
        optimizer,
        kwargs...
    )
end


"""
    var_core(mode, optimizer; verbose)

Function to calculte a stable profit distribution that belongs to the core
and minimizes the variance of the profit allocation among the plauers
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
var_core : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function var_core(mode::EnumMode, optimizer; kwargs...)

    # create the objective function for the problem
    function var_objective(m, player_set)
        obj = sum(m[:profit_dist].^2)/length(m[:profit_dist]) - sum(m[:profit_dist]/length(m[:profit_dist]))^2
        @objective(m, Min, obj)
    end

    return specific_in_core(
        mode,
        var_objective,
        optimizer;
        kwargs...
    )
end


"""
    ref_in_core(mode, ref_dist, optimizer; verbose)

Function to calculte a stable profit distribution that belongs to the core
and minimizes the variance of the profit allocation among the players
with respect to a pre-defined reward distribution function
for a game described by the utility function and the grand coalition of player_set.


Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
ref_dist : AbstrctDict
    Reference distribution by player
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
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
        mode::EnumMode,
        ref_dist,
        optimizer;
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

    return specific_in_core(
        mode,
        var_objective,
        optimizer; 
        kwargs...
    )
end


"""
    verify_in_core(profit_dist, mode, optimizer; verbose, utilities)

Function to calculte whether a given profit distribution belongs to the core.
The game shall be described by the utility function and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function verify_in_core(profit_dist, mode::EnumMode, optimizer; kwargs...)
    
    
    function verify_in_core_obj_and_constraint(m, player_set)
        # create feasibility problem
        @objective(m, Max, 0.0)

        # fix profit distribution
        for p in axes(m[:profit_dist])[1]
            fix(m[:profit_dist][p], profit_dist[p], force=true)
        end
    end

    # calculate return value
    ret_value = specific_in_core(
        mode,
        verify_in_core_obj_and_constraint,
        optimizer;
        kwargs...
    )

    # if return value is nothing, the problem is infeasible, hence the solution does not belong to the core
    return !isnothing(ret_value)
end


"""
    specific_least_core(mode, dist_objective, optimizer; verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes a specific objective for the profit allocation among the plauers,
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
dist_objective : Function
    Function to build the objective function of the profit distribution, after the least core
    is obtained
    It shall be a function with two arguments: the JuMP model of the problem and the player set,
    and the function shall build the desired objective functions using the 
    JuMP @NLobjective or @objective commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
utilities (optional, default nothing)

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function specific_least_core(
        mode::EnumMode,
        dist_objective::Function,
        optimizer; 
        verbose=true,
    )

    utilities = mode.utilities
    player_set = mode.player_set

    # total combinations
    comb_set = combinations(player_set)

    # initialize JuMP model
    model_dist = Model(optimizer)

    @variable(model_dist, minimum(values(utilities)) <= profit_dist[u in player_set] <= maximum(values(utilities)))

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == utilities[Set(player_set)])

    # the gain of the worst group of the current iteration
    @variable(model_dist, minimum(values(utilities)) <= min_surplus <= maximum(values(utilities)))

    # specify that the profit of each subset of the group is better off with the grand coalition
    @constraint(model_dist, con_subset_profit[comb in comb_set; length(comb) > 0],
        sum([profit_dist[pl] for pl in comb]) >= utilities[Set(comb)] + min_surplus
    )

    # specify the objective maximize the benefit or remaining in the coalition
    @objective(model_dist, Max, min_surplus)

    # optimize and get the value of min_surplus
    optimize!(model_dist)

    # fix the value of min_surplus
    fix(min_surplus, value(min_surplus), force=true)

    # Update objective function according to the desired dist_objective function
    dist_objective(model_dist, player_set)

    # re-optimize to get the result according with dist_objective
    optimize!(model_dist)

    specific_least_core_dist = Dict(zip(player_set, value.(profit_dist).data))

    return specific_least_core_dist
end


"""
    var_least_core(mode, optimizer; verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes the variance of the profit allocation among the players
for a game described by the utility function and the grand coalition of player_set.


Inputs
------
mode : EnumMode
    Calculation mode: enumerative technique
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
specific_least_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function var_least_core(
        mode::EnumMode,
        optimizer; 
        kwargs...
    )

    # create the objective function for the problem
    function var_objective(m, player_set)
        obj = sum(m[:profit_dist].^2)/length(m[:profit_dist]) - sum(m[:profit_dist]/length(m[:profit_dist]))^2
        @objective(m, Min, obj)
    end

    return specific_least_core(
        mode,
        var_objective,
        optimizer; 
        kwargs...
    )
end


"""
    ref_least_core(mode, ref_dist, optimizer; norm, verbose)

Function to calculte a stable profit distribution that belongs to the least core
and minimizes the variance of the profit allocation among the players
with respect to a pre-defined reward distribution function
for a game described by the utility function and the grand coalition of player_set.


Inputs
------
mode : AbstractCalcMode
    Calculation mode: enumerative technique
ref_dist : AbstrctDict
    Reference distribution by player
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
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
        mode::EnumMode,
        ref_dist,
        optimizer;
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
        mode,
        var_objective,
        optimizer;
        kwargs...
    )
end