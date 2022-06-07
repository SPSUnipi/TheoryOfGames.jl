"""
least_core(mode, player_set, utilities, optimizer; verbose)

Function to calculte the least core profit distribution for a game described by the utility function
and the grand coalition of player_set.
This function implements an iterative approach by using a column generating approach
to iteratively add constraints to the master problem.
To do so, the callback function stored in the mode category is exploited

Inputs
------
mode : IterMode
    Calculation mode that contains the reference to the callback function
    that given a profit distribution scheme, returns a tuple of
    the set of the coalition with the worst profit and its total benefit to be shared
    callback_worst_coalition accepts one argument (current profit sharing)
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
rtol : Number (optional, default 1e-6)
    Relative tolerance of the convergence process
atol : Number (optional, default 1e-6)
    Absolute tolerance of the convergence process
lower_bound : Number (optional, default nothing)
    Lower bound of the variables of the problem (benefit distribution and margin of the worst coalition)
upper_bound : Number (optional, default nothing)
    Upper bound of the variables of the problem (benefit distribution and margin of the worst coalition)
    When nothing, the value is automatically set to the benefit of the grand coalition
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
raw_outputs : Bool (optional, default false)
    When true, it returns all raw outputs
use_start_value : Bool (optional, default false)
    When true, in the iterative process the previous iteration value is used as initialization
    for the followin iteration
max_iter : Bool (optional, default 100)
    Maximum number of iterations of the process
preload_coalitions : Vector (optional, default empty)
    List of coalitions whose benefit shall be automatically
    included before the iterative procedure starts

Outputs
------
profit_distribution
    Benefit distribution by player
min_surplus
    Benefit of the coalition with the minimum surplus
history
"""
function least_core(
        mode::IterMode,
        optimizer;
        rtol=1e-6,
        atol=1e-6, 
        lower_bound=nothing,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
        use_start_value=false,
        max_iter=100,
        preload_coalitions=[],
        exclude_visited_coalitions=true,
    )

    player_set = mode.player_set
    callback_worst_coalition = mode.callback_worst_coalition
    callback_benefit_by_coalition = mode.callback_benefit_by_coalition

    # benefit of the grand coalition
    benefit_grand_coalition = callback_benefit_by_coalition(player_set)

    # upper_bound is nothing, set it to the benefit of the grand coalition
    if isnothing(upper_bound)
        upper_bound = benefit_grand_coalition
    end
    # lower_bound is nothing, set it to the benefit of the grand coalition
    if isnothing(lower_bound)
        lower_bound = -benefit_grand_coalition
    end

    # initialize JuMP model
    model_dist = Model(optimizer)

    # profit reward distribution by player
    @variable(model_dist, lower_bound <= profit_dist[u in player_set] <= upper_bound)

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == benefit_grand_coalition)

    # the gain of the worst group of the current iteration
    @variable(model_dist, -upper_bound <= min_surplus <= upper_bound)

    # specify the objective maximize the benefit or remaining in the coalition
    @objective(model_dist, Max, min_surplus)

    # Add constraints for preloaded coalitions
    @constraint(model_dist, preloaded_coalitions[precoal in preload_coalitions],
        sum(GenericAffExpr{Float64,VariableRef}[profit_dist[pl] for pl in precoal]) >= callback_benefit_by_coalition(precoal) + min_surplus
    )

    if !isempty(preload_coalitions)
        println("Preload $(length(preloaded_coalitions)) coalitions")
    end

    # initialization of the min_surplus of the inner problem
    lower_problem_min_surplus = lower_bound

    # list of constraints
    history = NamedTuple[]

    # initialization while condition
    continue_while = true
    iter = 0
    value_min_surplus = nothing

    # visited coalitions
    visited_coalitions = [Set.(preload_coalitions); Set(player_set); Set([])]

    # printing formats
    format_print_head = "{:<15s} {:<15s} {:<15s} {:<15s} {:<s}"  # for the header
    format_print_iter = "{:<15s} {:< 15.2e} {:< 15.2e} {:< 15.2f} {:<s}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Current value", "Tol. [%]", "Benefit distribution")
    end

    
    while continue_while && iter <= max_iter

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        value_min_surplus = value(min_surplus)

        # get current profit distribution
        current_profit_dist = value.(profit_dist)

        # get a vector in which each row contains a named tuple with:
        # (1) the coalition (worst_coal_set) with the least benefit [least_profitable_coalition],
        # (2) the total benefit of that coalition (worst_coal_benefit) [coalition_benefit], and
        # (3) the minimum surplus of that coalition [min_surplus]
        output_data = callback_worst_coalition(current_profit_dist)

        # get the minimum surplus of the iteration
        lower_problem_min_surplus = output_data[1].min_surplus

        # check if convergence has been reached
        if isapprox(value_min_surplus, lower_problem_min_surplus, rtol=rtol, atol=atol, norm=abs)
            # convergence reached
            continue_while = false
        else
            # convergence not reached: add new constraint

            if use_start_value
                # set initial start value
                value_start = value.(all_variables(model_dist))
            end

            for row in output_data

                least_profitable_coalition = Set(row.least_profitable_coalition)

                if !exclude_visited_coalitions || least_profitable_coalition ∉ visited_coalitions

                    # update visited_coalitions
                    push!(visited_coalitions, least_profitable_coalition)

                    # specify that the profit of each subset of the group is better off with the grand coalition
                    con_it = @constraint(
                        model_dist,
                        sum(GenericAffExpr{Float64,VariableRef}[profit_dist[pl] for pl in least_profitable_coalition]) >= row.coalition_benefit + min_surplus
                    )

                    # data of the iteration
                    iter_data = (
                        iteration=iter,
                        current_profit=current_profit_dist,
                        worst_coal=least_profitable_coalition,
                        benefit_coal=row.coalition_benefit,
                        value_min_surplus=value_min_surplus,
                        lower_problem_min_surplus=lower_problem_min_surplus,
                        constraint=con_it,
                    )
            
                    # add the iteration to the history
                    push!(history, iter_data)
                end
            end
        end

        # if verbose print log
        if verbose
            # update coalition
            benefit_distribution = string(current_profit_dist.data)

            rel_tol = compute_relative_tol(lower_problem_min_surplus, value_min_surplus)
            printfmtln(format_print_iter, iter, value_min_surplus, lower_problem_min_surplus, rel_tol*100, benefit_distribution)
        end
    end

    # result of the profit distribution
    profit_distribution = Dict(zip(player_set, value.(profit_dist).data))

    if raw_outputs
        return profit_distribution, value_min_surplus, history, visited_coalitions, model_dist
    else
        return profit_distribution
    end

end


"""
specific_least_core(mode, dist_objective, player_set, utilities, optimizer; verbose)

Function to calculte the least core profit distribution for a game described by the utility function
and the grand coalition of player_set.
This function implements an iterative approach by using a column generating approach
to iteratively add constraints to the master problem.
To do so, the callback function stored in the mode category is exploited

Inputs
------
mode : IterMode
    Calculation mode that contains the reference to the callback function
    that given a profit distribution scheme, returns a tuple of
    the set of the coalition with the worst profit and its total benefit to be shared
    callback_worst_coalition accepts one argument (current profit sharing)
dist_objective : Function
    Function to build the objective function of the profit distribution, after the least core
    is obtained
    It shall be a function with two arguments: the JuMP model of the problem and the player set,
    and the function shall build the desired objective functions using the 
    JuMP @NLobjective or @objective commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
rtol : Number (optional, default 1e-6)
    Relative tolerance of the convergence process
atol : Number (optional, default 1e-6)
    Absolute tolerance of the convergence process
lower_bound : Number (optional, default nothing)
    Lower bound of the variables of the problem (benefit distribution and margin of the worst coalition)
upper_bound : Number (optional, default nothing)
    Upper bound of the variables of the problem (benefit distribution and margin of the worst coalition)
    When nothing, the value is automatically set to the benefit of the grand coalition
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
raw_outputs : Bool (optional, default false)
    When true, it returns all raw outputs
use_start_value : Bool (optional, default false)
    When true, in the iterative process the previous iteration value is used as initialization
    for the followin iteration
max_iter : Bool (optional, default 100)
    Maximum number of iterations of the process

Outputs
------
profit_distribution
    Benefit distribution by player
min_surplus
    Benefit of the coalition with the minimum surplus
history
"""
function specific_least_core(
        mode::IterMode,
        dist_objective::Function,
        optimizer;
        rtol=1e-6,
        atol=1e-6, 
        lower_bound=nothing,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
        use_start_value=false,
        max_iter=100,
        exclude_visited_coalitions=true,
        kwargs...
    )

    if verbose
        println("PHASE 1: Start first least core analysis\n")
    end

    profit_distribution, min_surplus, history, visited_coalitions, model_dist = least_core(
        mode, optimizer;
        rtol=rtol,
        atol=atol,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        verbose=verbose,
        raw_outputs=true,
        max_iter=max_iter,
        exclude_visited_coalitions=exclude_visited_coalitions,
        kwargs...
    )
    
    if verbose
        println("PHASE 1: First least core analysis ended with value $min_surplus\n ----")
        println("PHASE 2: Start second phase to indentify the desired objective\n")
    end

    # get data of the IterMode
    player_set = mode.player_set
    callback_worst_coalition = mode.callback_worst_coalition
    callback_benefit_by_coalition = mode.callback_benefit_by_coalition

    # update objective function of the least_core model
    dist_objective(model_dist, player_set)

    # fix value of the min_surplus
    fix(model_dist[:min_surplus], min_surplus, force=true)

    # initialization while condition
    continue_while = true
    iter = history[end].iteration

    # printing formats
    format_print_head = "{:<15s} {:<15s} {:<15s} {:<15s} {:<s}"  # for the header
    format_print_iter = "{:<15s} {:< 15.2e} {:< 15.2e} {:< 15.2f} {:<s}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Current value", "Tol. [%]", "Benefit distribution")
    end

    # initialize variable
    current_profit_dist = nothing

    while continue_while && iter <= max_iter

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        # get current profit distribution
        current_profit_dist = value.(model_dist[:profit_dist])

        # get a vector in which each row contains a named tuple with:
        # (1) the coalition (worst_coal_set) with the least benefit [least_profitable_coalition],
        # (2) the total benefit of that coalition (worst_coal_benefit) [coalition_benefit], and
        # (3) the minimum surplus of that coalition [min_surplus]
        output_data = callback_worst_coalition(current_profit_dist)

        # get the minimum surplus of the iteration
        lower_problem_min_surplus = output_data[1].min_surplus

        # check if convergence has been reached
        if isapprox(lower_problem_min_surplus, min_surplus, rtol=rtol, atol=atol, norm=abs)
            # convergence reached
            continue_while = false
        else

            if use_start_value
                # set initial start value
                value_start = value.(all_variables(model_dist))
            end

            for row in output_data

                least_profitable_coalition = Set(row.least_profitable_coalition)

                if !exclude_visited_coalitions || least_profitable_coalition ∉ visited_coalitions

                    # update visited_coalitions
                    push!(visited_coalitions, least_profitable_coalition)

                    # specify that the profit of each subset of the group is better off with the grand coalition
                    con_it = @constraint(
                        model_dist,
                        sum(GenericAffExpr{Float64,VariableRef}[model_dist[:profit_dist][pl] for pl in least_profitable_coalition]) >= row.coalition_benefit + model_dist[:min_surplus]
                    )

                    # data of the iteration
                    iter_data = (
                        iteration=iter,
                        current_profit=current_profit_dist,
                        worst_coal=least_profitable_coalition,
                        benefit_coal=row.coalition_benefit,
                        value_min_surplus=min_surplus,
                        lower_problem_min_surplus=lower_problem_min_surplus,
                        constraint=con_it,
                    )
            
                    # add the iteration to the history
                    push!(history, iter_data)
                end
            end

            if use_start_value
                # set start value
                set_start_value.(all_variables(model_dist), value_start)
            end
        end

        # if verbose print log
        if verbose
            # update coalition
            benefit_distribution = string(current_profit_dist.data)

            rel_tol = compute_relative_tol(lower_problem_min_surplus, min_surplus)
            printfmtln(format_print_iter, iter, min_surplus, lower_problem_min_surplus, rel_tol*100, benefit_distribution)
        end
    end

    # result of the profit distribution
    profit_distribution = Dict(zip(player_set, current_profit_dist.data))
    

    if verbose
        println("PHASE 2: Completed second phase\n")
    end

    if raw_outputs
        return profit_distribution, min_surplus, history, model_dist
    else
        return profit_distribution
    end

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
        mode::IterMode,
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
        mode::IterMode,
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



"""
    specific_in_core(mode, player_set, utilities, optimizer; verbose)

Function to calculte a profit distribution that belongs to the core 
and maximizes a specific distribution objective specified by dist_objective 
for a game described by the utility function and the grand coalition of player_set.
This function implements an iterative approach by using a column generating approach
to iteratively add constraints to the master problem.
To do so, the callback function stored in the mode category is exploited

Inputs
------
mode : IterMode
    Calculation mode that contains the reference to the callback function
    that given a profit distribution scheme, returns a tuple of
    the set of the coalition with the worst profit and its total benefit to be shared
    callback_worst_coalition accepts one argument (current profit sharing)
dist_objective : Function
    Function to build the objective function of the profit distribution
    It shall a function with two arguments, which are the JuMP model and the player set,
    and it shall build the objective functions using the JuMP @NLobjective or @objective
    commands
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
rtol : Number (optional, default 1e-6)
    Relative tolerance of the convergence process
atol : Number (optional, default 1e-6)
    Absolute tolerance of the convergence process
lower_bound : Number (optional, default nothing)
    Lower bound of the variables of the problem (benefit distribution and margin of the worst coalition)
upper_bound : Number (optional, default nothing)
    Upper bound of the variables of the problem (benefit distribution and margin of the worst coalition)
    When nothing, the value is automatically set to the benefit of the grand coalition
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
raw_outputs : Bool (optional, default false)
    When true, it returns all raw outputs
use_start_value : Bool (optional, default false)
    When true, in the iterative process the previous iteration value is used as initialization
    for the followin iteration
max_iter : Bool (optional, default 100)
    Maximum number of iterations of the process
preload_coalitions : Vector (optional, default empty)
    List of coalitions whose benefit shall be automatically
    included before the iterative procedure starts

Outputs
------
profit_distribution
    Benefit distribution by player
min_surplus
    Benefit of the coalition with the minimum surplus
history
"""
function specific_in_core(
        mode::IterMode,
        dist_objective::Function,
        optimizer;
        rtol=1e-6,
        atol=1e-6, 
        lower_bound=nothing,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
        use_start_value=false,
        max_iter=100,
        preload_coalitions=[],
        exclude_visited_coalitions=true,
        kwargs...
    )

    player_set = mode.player_set
    callback_worst_coalition = mode.callback_worst_coalition
    callback_benefit_by_coalition = mode.callback_benefit_by_coalition

    # benefit of the grand coalition
    benefit_grand_coalition = callback_benefit_by_coalition(player_set)

    # lower_bound is nothing, set it to the benefit of the grand coalition
    if isnothing(lower_bound)
        lower_bound = -benefit_grand_coalition
    end
    # upper_bound is nothing, set it to the benefit of the grand coalition
    if isnothing(upper_bound)
        upper_bound = benefit_grand_coalition
    end

    # initialize JuMP model
    model_dist = Model(optimizer)

    # profit reward distribution by player
    @variable(model_dist, lower_bound <= profit_dist[u in player_set] <= upper_bound)

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == benefit_grand_coalition)

    # Add constraints for preloaded coalitions
    @constraint(model_dist, preloaded_coalitions[precoal in preload_coalitions],
        sum(GenericAffExpr{Float64,VariableRef}[profit_dist[pl] for pl in precoal]) >= callback_benefit_by_coalition(precoal) # + 0.0 # set 0.0 to make it belong to the core
    )

    if !isempty(preload_coalitions)
        println("Preload $(length(preloaded_coalitions)) coalitions")
    end

    # specify function
    dist_objective(model_dist, player_set)

    # initialization of the min_surplus of the inner problem
    lower_problem_min_surplus = lower_bound

    # list of constraints
    history = NamedTuple[]

    # initialization while condition
    continue_while = true
    iter = 0
    visited_coalitions = [Set.(preload_coalitions); Set(player_set); Set([])]

    # printing formats
    format_print_head = "{:<15s} {:<15s} {:<15s} {:<15s} {:<s}"  # for the header
    format_print_iter = "{:<15s} {:< 15.2e} {:< 15.2e} {:< 15.2f} {:<s}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Current value", "Tol. [%]", "Benefit allocation")
    end

    while continue_while && iter <= max_iter

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        # get current profit distribution
        current_profit_dist = value.(profit_dist)
        # get a vector in which each row contains a named tuple with:
        # (1) the coalition (worst_coal_set) with the least benefit [least_profitable_coalition],
        # (2) the total benefit of that coalition (worst_coal_benefit) [coalition_benefit], and
        # (3) the minimum surplus of that coalition [min_surplus]
        output_data = callback_worst_coalition(current_profit_dist)

        # get the minimum surplus of the iteration
        lower_problem_min_surplus = output_data[1].min_surplus

        # check if convergence has been reached: when the lower problem min surplus is non-negative
        if lower_problem_min_surplus * (1+rtol) + atol >= 0
            # convergence reached
            continue_while = false
        else
            # convergence not reached: add new constraint

            if use_start_value
                # set initial start value
                value_start = value.(all_variables(model_dist))
            end

            for row in output_data

                if !exclude_visited_coalitions || row.least_profitable_coalition ∉ visited_coalitions

                    least_profitable_coalition = Set(row.least_profitable_coalition)

                    # update visited_coalitions
                    push!(visited_coalitions, least_profitable_coalition)

                    # specify that the profit of each subset of the group is better off with the grand coalition
                    con_it = @constraint(
                        model_dist,
                        sum(GenericAffExpr{Float64,VariableRef}[profit_dist[pl] for pl in least_profitable_coalition]) >= row.coalition_benefit # + 0.0 # set 0.0 to make it belong to the core
                    )

                    # data of the iteration
                    iter_data = (
                        iteration=iter,
                        current_profit=current_profit_dist,
                        worst_coal=least_profitable_coalition,
                        benefit_coal=row.coalition_benefit,
                        value_min_surplus=0.0,
                        lower_problem_min_surplus=lower_problem_min_surplus,
                        constraint=con_it,
                    )
            
                    # add the iteration to the history
                    push!(history, iter_data)
                end
            end

            if use_start_value
                # set start value
                set_start_value.(all_variables(model_dist), value_start)
            end
        end

        # if verbose print log
        if verbose
            # update coalition
            benefit_distribution = string(current_profit_dist.data)

            rel_tol = compute_relative_tol(lower_problem_min_surplus, 0.0)
            printfmtln(format_print_iter, iter, lower_problem_min_surplus, 0.0, rel_tol*100, benefit_distribution)
        end
    end

    # result of the profit distribution
    profit_distribution = Dict(zip(player_set, value.(profit_dist).data))

    if raw_outputs
        return profit_distribution, 0.0, history, model_dist
    else
        return profit_distribution
    end

end


"""
    in_core(mode, optimizer; verbose, ...)

Function to calculte a stable profit distribution that belongs to the core
for profit distribution for a game described by the utility function
and the grand coalition of player_set.

Inputs
------
mode : IterMode
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
function in_core(mode::IterMode, optimizer; kwargs...)
    
    # create feasibility problem
    in_core_objective(m, player_set) = @objective(m, Max, 0.0)

    return specific_in_core(
        mode,
        in_core_objective,  # the core is only a feasibility problem, thus specify constant obj. funciton
        optimizer;
        kwargs...
    )
end


"""
    var_in_core(mode, optimizer; verbose)

Function to calculte a stable profit distribution that belongs to the core
and minimizes the variance of the profit allocation among the plauers
for a game described by the utility function and the grand coalition of player_set.

Inputs
------
mode : IterMode
    Calculation mode: enumerative technique
optimizer : Any
    Optimizer function for the JuMP model used for computation purposes
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
var_in_core : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function var_in_core(mode::IterMode, optimizer; kwargs...)

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
mode : IterMode
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
specific_ref_in_core_dist : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function ref_in_core(
        mode::IterMode,
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
function verify_in_core(
        profit_dist, 
        mode::IterMode,
        optimizer; 
        rtol=1e-6,
        atol=1e-6,
        kwargs...
    )
    
    
    out_data = mode.callback_worst_coalition(profit_dist)

    lower_problem_min_surplus = out_data[1].min_surplus

    # if return value is nothing, the problem is infeasible, hence the solution does not belong to the core
    return lower_problem_min_surplus*(1+rtol) + atol >= 0.0
end