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
start_time : (optional, default nothing)
    Specify the initial time of the method; if nothing the value is initialized by time()
rtol : Number (optional, default 1e-2)
    Relative tolerance of the convergence process
atol : Number (optional, default 1e-2)
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
best_objective_stop_option : String (optional, default nothing)
    Name of the option to stop the lower problem as it reaches a preset value.
    When this option is nothing, this feature is not used.
    When this option is non-nothing, in every iteration, a minimum convergence criterion is added
    so to stop the lower problem as soon as a minimum fesible objective function is reached.
    This minimum objective value is obtained with respect to the solution of the master problem
    multiplied by the factor "best_objective_stop_tolerance"
    If gurobi is used, this option is BestObjStop
best_objective_stop_tolerance : Number (optional, default 0.05)
    Tolerance used in the "best_objective_stop_option" approach
lower_relaxation_stop_option : String (optional, default nothing)
    Name of the option used to setup the stop criterion of the optimization as soon as
    the lowest bound reaches the tolerance specified by tolerance_lower_relaxation_stop option
tolerance_lower_relaxation_stop : double (optional, default 0.0)
    When lower_relaxation_stop_option is enabled, this option specifies the tolerance used
    to stop the loop

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
        start_time=nothing,
        rtol=1e-2,
        atol=1e-2, 
        lower_bound=nothing,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
        use_start_value=false,
        max_iter=100,
        preload_coalitions=[],
        exclude_visited_coalitions=true,
        best_objective_stop_option=nothing,
        best_objective_stop_tolerance=0.05,
        lower_relaxation_stop_option=nothing,
        tolerance_lower_relaxation_stop=0.0,
    )

    if !isnothing(best_objective_stop_option) && exclude_visited_coalitions
        println("WARNING: BestObjStop enabled alongside exclude_visited_coalitions. This may lead to infinite loops")
    end

    if isnothing(start_time)
        start_time = time()
    end

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
    precoal_constraint = @constraint(model_dist, [precoal in preload_coalitions],
        sum(GenericAffExpr{Float64,VariableRef}[profit_dist[pl] for pl in precoal]) >= callback_benefit_by_coalition(precoal) + min_surplus
    )

    if !isempty(preload_coalitions)
        println("Preload $(length(preload_coalitions)) coalitions")
    end

    # initialization of the min_surplus of the inner problem
    lower_problem_min_surplus = lower_bound

    # list of constraints
    history = NamedTuple[_create_history_row(
        0,
        time() - start_time,
        nothing,
        nothing,
        NaN,
        upper_bound,
        lower_bound,
        precoal_constraint,
    )]

    # initialization while condition
    continue_while = true
    iter = 0
    value_min_surplus = nothing

    # visited coalitions
    visited_coalitions = [Set.(preload_coalitions); Set(player_set); Set([])]
    looping_ongoing = false  # true when the iterative algorithm is looping on the same solution

    # printing formats
    format_print_head = "{:<15s} {:<15s} {:<15s} {:<15s} {:<s}"  # for the header
    format_print_iter = "{:<15s} {:< 15.2e} {:< 15.2e} {:< 15.2f} {:<s}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Lower bound", "Tol. [%]", "Benefit distribution")
    end

    
    while continue_while && iter <= max_iter

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        value_min_surplus = value(min_surplus)

        # get current profit distribution
        current_profit_dist = value.(profit_dist)

        # update best objective stop if applicable
        modify_solver_options = []
        if !isnothing(best_objective_stop_option)
            if looping_ongoing
                println("Looping identified on the same solution, skipped option $best_objective_stop_option")
            else
                best_obj_stop = value_min_surplus - max(abs(value_min_surplus) * best_objective_stop_tolerance, atol)
                push!(modify_solver_options, string(best_objective_stop_option)=>best_obj_stop)
            end
        end
        if !isnothing(lower_relaxation_stop_option)
            push!(modify_solver_options, string(lower_relaxation_stop_option)=>value_min_surplus - abs(value_min_surplus) * tolerance_lower_relaxation_stop)
        end

        # get a vector in which each row contains a named tuple with:
        # (1) a vector-like representing the fractional participance of each uer to the worst coalition
        #     when the value is 1, that user belongs to the worst coalition, 0 does not belong
        # (2) the coalition (worst_coal_set) with the least benefit [least_profitable_coalition],
        # (3) the total benefit of that coalition (worst_coal_benefit) [coalition_benefit], and
        # (4) the minimum surplus of that coalition [min_surplus]
        output_data = callback_worst_coalition(
            current_profit_dist;
            modify_solver_options=modify_solver_options,
            decompose_ANC_lower_obj_stop=value_min_surplus,
        )

        # get the minimum surplus of the iteration
        lower_problem_min_surplus = output_data[1].min_surplus

        # check if convergence has been reached
        if isapprox(value_min_surplus, lower_problem_min_surplus, rtol=rtol, atol=atol, norm=abs)
            # convergence reached
            continue_while = false

            # data of the iteration
            iter_data = _create_history_row(
                iter,
                time()-start_time,
                current_profit_dist,
                output_data[1].least_profitable_coalition_status,
                output_data[1].coalition_benefit,
                value_min_surplus,
                lower_problem_min_surplus,
                nothing,
            )
    
            # add the iteration to the history
            push!(history, iter_data)
        else
            # convergence not reached: add new constraint

            # initialize to true and if at least a solution is included, then set it to false
            looping_ongoing = true

            if use_start_value
                # set initial start value
                value_start = value.(all_variables(model_dist))
            end

            for row in output_data

                least_profitable_coalition_status = row.least_profitable_coalition_status
                least_profitable_coalition = Set(row.least_profitable_coalition)

                # turn looping ongoing to false if all solutions have already been visited
                looping_ongoing &= (least_profitable_coalition ∈ visited_coalitions)

                con_it = nothing

                if !exclude_visited_coalitions || least_profitable_coalition ∉ visited_coalitions

                    # update visited_coalitions
                    push!(visited_coalitions, least_profitable_coalition)

                    # specify that the profit of each subset of the group is better off with the grand coalition
                    con_it = @constraint(
                        model_dist,
                        sum(GenericAffExpr{Float64,VariableRef}[
                            least_profitable_coalition_status[pl] * profit_dist[pl] for pl in player_set
                        ]) >= row.coalition_benefit + min_surplus
                    )
                end

                # data of the iteration
                iter_data = _create_history_row(
                    iter,
                    time()-start_time,
                    current_profit_dist,
                    least_profitable_coalition_status,
                    row.coalition_benefit,
                    value_min_surplus,
                    lower_problem_min_surplus,
                    con_it,
                )
        
                # add the iteration to the history
                push!(history, iter_data)
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
start_time : (optional, default nothing)
    Specify the initial time of the method; if nothing the value is initialized by time()
rtol : Number (optional, default 1e-2)
    Relative tolerance of the convergence process
atol : Number (optional, default 1e-2)
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
best_objective_stop_option : String (optional, default nothing)
    Name of the option to stop the lower problem as it reaches a preset value.
    When this option is nothing, this feature is not used.
    When this option is non-nothing, in every iteration, a minimum convergence criterion is added
    so to stop the lower problem as soon as a minimum fesible objective function is reached.
    This minimum objective value is obtained with respect to the solution of the master problem
    multiplied by the factor "best_objective_stop_tolerance"
    If gurobi is used, this option is BestObjStop
best_objective_stop_tolerance : Number (optional, default 0.05)
    Tolerance used in the "best_objective_stop_option" approach
lower_relaxation_stop_option : String (optional, default nothing)
    Name of the option used to setup the stop criterion of the optimization as soon as
    the lowest bound reaches the tolerance specified by tolerance_lower_relaxation_stop option
tolerance_lower_relaxation_stop : double (optional, default 0.0)
    When lower_relaxation_stop_option is enabled, this option specifies the tolerance used
    to stop the loop

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
        start_time=nothing,
        rtol=1e-2,
        atol=1e-2, 
        lower_bound=nothing,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
        use_start_value=false,
        max_iter=100,
        exclude_visited_coalitions=true,
        best_objective_stop_option=nothing,
        best_objective_stop_tolerance=0.05,
        lower_relaxation_stop_option=nothing,
        tolerance_lower_relaxation_stop=0.0,
        kwargs...
    )

    if !isnothing(best_objective_stop_option) && exclude_visited_coalitions
        println("WARNING: BestObjStop enabled alongside exclude_visited_coalitions. This may lead to infinite loops")
    end

    if isnothing(start_time)
        start_time = time()
    end

    if verbose
        println("PHASE 1: Start first least core analysis\n")
    end

    profit_distribution, min_surplus, history, visited_coalitions, model_dist = least_core(
        mode, optimizer;
        start_time=start_time,
        rtol=rtol,
        atol=atol,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        verbose=verbose,
        raw_outputs=true,
        max_iter=max_iter,
        exclude_visited_coalitions=exclude_visited_coalitions,
        lower_relaxation_stop_option=lower_relaxation_stop_option,
        best_objective_stop_option=best_objective_stop_option,
        best_objective_stop_tolerance=best_objective_stop_tolerance,
        tolerance_lower_relaxation_stop=tolerance_lower_relaxation_stop,
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
        printfmtln(format_print_head, "Iteration", "Upper bound", "Lower bound", "Tol. [%]", "Benefit distribution")
    end

    # initialize variable
    current_profit_dist = nothing
    looping_ongoing = false  # true when the iterative algorithm is looping on the same solution

    while continue_while && iter <= max_iter

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        # get current profit distribution
        current_profit_dist = value.(model_dist[:profit_dist])

        # update best objective stop if applicable
        modify_solver_options = []
        if !isnothing(best_objective_stop_option)
            if looping_ongoing
                println("Looping identified on the same solution, skipped option $best_objective_stop_option")
            else
                best_obj_stop = min_surplus - max(abs(min_surplus) * best_objective_stop_tolerance, atol)
                push!(modify_solver_options, string(best_objective_stop_option)=>best_obj_stop)
            end
        end
        if !isnothing(lower_relaxation_stop_option)
            push!(modify_solver_options, string(lower_relaxation_stop_option)=>min_surplus - abs(min_surplus) * tolerance_lower_relaxation_stop)
        end

        # get a vector in which each row contains a named tuple with:
        # (1) a vector-like representing the fractional participance of each uer to the worst coalition
        #     when the value is 1, that user belongs to the worst coalition, 0 does not belong
        # (2) the coalition (worst_coal_set) with the least benefit [least_profitable_coalition],
        # (3) the total benefit of that coalition (worst_coal_benefit) [coalition_benefit], and
        # (4) the minimum surplus of that coalition [min_surplus]
        output_data = callback_worst_coalition(
            current_profit_dist;
            modify_solver_options=modify_solver_options,
            decompose_ANC_lower_obj_stop=min_surplus,
        )

        # get the minimum surplus of the iteration
        lower_problem_min_surplus = output_data[1].min_surplus

        # check if convergence has been reached
        if isapprox(lower_problem_min_surplus, min_surplus, rtol=rtol, atol=atol, norm=abs)
            # convergence reached
            continue_while = false

            # data of the iteration
            iter_data = _create_history_row(
                iter,
                time() - start_time,
                current_profit_dist,
                output_data[1].least_profitable_coalition_status,
                output_data[1].coalition_benefit,
                min_surplus,
                lower_problem_min_surplus,
                nothing,
            )
    
            # add the iteration to the history
            push!(history, iter_data)
        else

            # initialize to true and if at least a solution is included, then set it to false
            looping_ongoing = true

            if use_start_value
                # set initial start value
                value_start = value.(all_variables(model_dist))
            end

            for row in output_data

                least_profitable_coalition_status = row.least_profitable_coalition_status
                least_profitable_coalition = Set(row.least_profitable_coalition)

                # turn looping ongoing to false if all solutions have already been visited
                looping_ongoing &= (least_profitable_coalition ∈ visited_coalitions)

                con_it = nothing

                if !exclude_visited_coalitions || least_profitable_coalition ∉ visited_coalitions

                    # update visited_coalitions
                    push!(visited_coalitions, least_profitable_coalition)

                    # specify that the profit of each subset of the group is better off with the grand coalition
                    con_it = @constraint(
                        model_dist,
                        sum(GenericAffExpr{Float64,VariableRef}[
                            least_profitable_coalition_status[pl] * model_dist[:profit_dist][pl]
                            for pl in player_set
                        ]) >= row.coalition_benefit + model_dist[:min_surplus]
                    )
                end

                # data of the iteration
                iter_data = _create_history_row(
                    iter,
                    time() - start_time,
                    current_profit_dist,
                    least_profitable_coalition_status,
                    row.coalition_benefit,
                    min_surplus,
                    lower_problem_min_surplus,
                    con_it,
                )
        
                # add the iteration to the history
                push!(history, iter_data)
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
numerical_multiplier : Float (default 1e-3)
    Multiplier to adjust numerical issues
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
        numerical_multiplier=0.001,
        kwargs...
    )

    # create the objective function for the problem
    function var_objective(m, player_set)
        obj = numerical_multiplier * (sum(m[:profit_dist].^2)/length(m[:profit_dist]) - sum(m[:profit_dist]/length(m[:profit_dist]))^2)
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
start_time : (optional, default nothing)
    Specify the initial time of the method; if nothing the value is initialized by time()
rtol : Number (optional, default 1e-2)
    Relative tolerance of the convergence process
atol : Number (optional, default 1e-2)
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
best_objective_stop_option : String (optional, default nothing)
    Name of the option to stop the lower problem as it reaches a preset value.
    When this option is nothing, this feature is not used.
    When this option is non-nothing, in every iteration, a minimum convergence criterion is added
    so to stop the lower problem as soon as a minimum fesible objective function is reached.
    This minimum objective value after which the solver shall return the solution
    is provided as the option best_objective_stop_value
    If gurobi is used, this option is BestObjStop
best_objective_stop_value : Number (optional, default -1e75)
    Minimum objective function used for the solver to converge
    When the procedure starts looping
lower_relaxation_stop_option : String (optional, default nothing)
    Name of the option used to setup the stop criterion of the optimization as soon as
    the lowest bound reaches the tolerance specified by tolerance_lower_relaxation_stop option
tolerance_lower_relaxation_stop : double (optional, default 0.0)
    When lower_relaxation_stop_option is enabled, this option specifies the tolerance used
    to stop the loop

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
        start_time=nothing,
        rtol=1e-2,
        atol=1e-2, 
        lower_bound=nothing,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
        use_start_value=false,
        max_iter=100,
        preload_coalitions=[],
        exclude_visited_coalitions=true,
        best_objective_stop_option=nothing,
        best_objective_stop_value=-1e75,
        lower_relaxation_stop_option=nothing,
        tolerance_lower_relaxation_stop=0.0,
        kwargs...
    )

    if !isnothing(best_objective_stop_option) && exclude_visited_coalitions
        println("WARNING: BestObjStop enabled alongside exclude_visited_coalitions. This may lead to infinite loops")
    end

    if isnothing(start_time)
        start_time = time()
    end

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
    precoal_constraint = @constraint(model_dist, [precoal in preload_coalitions],
        sum(GenericAffExpr{Float64,VariableRef}[profit_dist[pl] for pl in precoal]) >= callback_benefit_by_coalition(precoal) # + 0.0 # set 0.0 to make it belong to the core
    )

    if !isempty(preload_coalitions)
        println("Preload $(length(preload_coalitions)) coalitions")
    end

    # specify function
    dist_objective(model_dist, player_set)

    # initialization of the min_surplus of the inner problem
    lower_problem_min_surplus = lower_bound

    # list of constraints
    history = NamedTuple[_create_history_row(
        0,
        time() - start_time,
        nothing,
        nothing,
        NaN,
        0.0,
        NaN,
        precoal_constraint,
    )]

    # initialization while condition
    continue_while = true
    iter = 0
    visited_coalitions = [Set.(preload_coalitions); Set(player_set); Set([])]
    looping_ongoing = false  # true when the iterative algorithm is looping on the same solution

    # printing formats
    format_print_head = "{:<15s} {:<15s} {:<15s} {:<15s} {:<s}"  # for the header
    format_print_iter = "{:<15s} {:< 15.2e} {:< 15.2e} {:< 15.2f} {:<s}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Lower bound", "Tol. [%]", "Benefit allocation")
    end

    while continue_while && iter <= max_iter

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        # get current profit distribution
        current_profit_dist = value.(profit_dist)

        # update best objective stop if applicable
        modify_solver_options = []
        if !isnothing(best_objective_stop_option)
            if looping_ongoing
                println("Looping identified on the same solution, skipped option $best_objective_stop_option")
            else
                push!(modify_solver_options, string(best_objective_stop_option)=>best_objective_stop_value)
            end
        end
        if !isnothing(lower_relaxation_stop_option)
            push!(modify_solver_options, string(lower_relaxation_stop_option)=>0.0)
        end
        
        # get a vector in which each row contains a named tuple with:
        # (1) a vector-like representing the fractional participance of each uer to the worst coalition
        #     when the value is 1, that user belongs to the worst coalition, 0 does not belong
        # (2) the coalition (worst_coal_set) with the least benefit [least_profitable_coalition],
        # (3) the total benefit of that coalition (worst_coal_benefit) [coalition_benefit], and
        # (4) the minimum surplus of that coalition [min_surplus]
        output_data = callback_worst_coalition(
            current_profit_dist;
            modify_solver_options=modify_solver_options,
            decompose_ANC_lower_obj_stop=0.0,
        )

        # get the minimum surplus of the iteration
        lower_problem_min_surplus = output_data[1].min_surplus
        last_coalition = Set(output_data[1].least_profitable_coalition)

        # check if convergence has been reached: when the lower problem min surplus is non-negative
        if lower_problem_min_surplus + max(abs(lower_problem_min_surplus) * rtol, atol) >= 0
            # convergence reached
            continue_while = false

            # data of the iteration
            iter_data = _create_history_row(
                iter,
                time() - start_time,
                current_profit_dist,
                output_data[1].least_profitable_coalition_status,
                output_data[1].coalition_benefit,
                0.0,
                lower_problem_min_surplus,
                nothing,
            )
    
            # add the iteration to the history
            push!(history, iter_data)
        else
            # convergence not reached: add new constraint

            # initialize to true and if at least a solution is included, then set it to false
            looping_ongoing = true

            if use_start_value
                # set initial start value
                value_start = value.(all_variables(model_dist))
            end

            for row in output_data

                least_profitable_coalition_status = row.least_profitable_coalition_status
                least_profitable_coalition = Set(row.least_profitable_coalition)

                # turn looping ongoing to false if all solutions have already been visited
                looping_ongoing &= (least_profitable_coalition ∈ visited_coalitions)

                con_it = nothing

                if !exclude_visited_coalitions || least_profitable_coalition ∉ visited_coalitions

                    # update visited_coalitions
                    push!(visited_coalitions, least_profitable_coalition)

                    # specify that the profit of each subset of the group is better off with the grand coalition
                    con_it = @constraint(
                        model_dist,
                        sum(GenericAffExpr{Float64,VariableRef}[
                            least_profitable_coalition_status[pl] * profit_dist[pl]
                            for pl in player_set
                        ]) >= row.coalition_benefit # + 0.0 # set 0.0 to make it belong to the core
                    )
                end

                # data of the iteration
                iter_data = _create_history_row(
                    iter,
                    time() - start_time,
                    current_profit_dist,
                    least_profitable_coalition_status,
                    row.coalition_benefit,
                    0.0,
                    lower_problem_min_surplus,
                    con_it,
                )
        
                # add the iteration to the history
                push!(history, iter_data)
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
            printfmtln(format_print_iter, iter, 0.0, lower_problem_min_surplus, rel_tol*100, benefit_distribution)
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
numerical_multiplier : Float (default 1e-3)
    Multiplier to adjust numerical issues
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status

Outputs
------
var_in_core : Dict
    Dictionary of the fair distributions of the profits among the players
"""
function var_in_core(mode::IterMode, optimizer; numerical_multiplier=0.001, kwargs...)

    # create the objective function for the problem
    function var_objective(m, player_set)
        obj = numerical_multiplier * (sum(m[:profit_dist].^2)/length(m[:profit_dist]) - sum(m[:profit_dist]/length(m[:profit_dist]))^2)
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
numerical_multiplier : Float (default 1e-3)
    Multiplier to adjust numerical issues
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
        numerical_multiplier=0.001,
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
        ) * numerical_multiplier
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
        rtol=1e-2,
        atol=1e-2,
        kwargs...
    )
    
    
    out_data = mode.callback_worst_coalition(profit_dist)

    lower_problem_min_surplus = out_data[1].min_surplus

    # if return value is nothing, the problem is infeasible, hence the solution does not belong to the core
    return lower_problem_min_surplus + max(abs(lower_problem_min_surplus)*rtol, atol) >= 0.0
end

