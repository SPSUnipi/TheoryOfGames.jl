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
lower_bound : Number (optional, default 0.0)
    Lower bound of the variables of the problem (benefit distribution and margin of the worst coalition)
upper_bound : Number (optional, default nothing)
    Upper bound of the variables of the problem (benefit distribution and margin of the worst coalition)
    When nothing, the value is automatically set to the benefit of the grand coalition
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
raw_outputs : Bool (optional, default false)
    When true, it returns all raw outputs

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
        lower_bound=0.0,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
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

    # initialize JuMP model
    model_dist = Model(optimizer)

    # profit reward distribution by player
    @variable(model_dist, lower_bound <= profit_dist[u in player_set] <= upper_bound)

    # specify that the total profit distribution cannot exceed the total benefit
    @constraint(model_dist, con_total_benefit, sum(profit_dist) == benefit_grand_coalition)

    # the gain of the worst group of the current iteration
    @variable(model_dist, lower_bound <= min_surplus <= upper_bound)

    # specify the objective maximize the benefit or remaining in the coalition
    @objective(model_dist, Max, min_surplus)

    # initialization of the min_surplus of the inner problem
    lower_problem_min_surplus = lower_bound

    # list of constraints
    history = NamedTuple[]

    # initialization while condition
    continue_while = true
    iter = 0
    added_coalition = "Init"

    # printing formats
    format_print_head = "{:<15s} {:^15s} {:^15s} {:^15s} {:<s}"  # for the header
    format_print_iter = "{:<15s} {:> 13.2e} {:> 13.2e} {:> 13.2e} {:<s}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Current value", "Tolerance", "Added coalition : benefit allocation")
    end

    while continue_while

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        value_min_surplus = value(min_surplus)

        # if verbose print log
        if verbose
            rel_tol = compute_relative_tol(lower_problem_min_surplus, value_min_surplus)
            printfmtln(format_print_iter, iter, lower_problem_min_surplus, value_min_surplus, rel_tol, added_coalition)
        end

        # check if convergence has been reached
        if isapprox(value_min_surplus, lower_problem_min_surplus, rtol=rtol, atol=atol, norm=abs)
            # convergence reached
            continue_while = false
        else
            # convergence not reached: add new constraint

            current_profit_dist = value.(profit_dist)

            # get (1) the coalition (worst_coal_set) with the least benefit,
            # (2) the total benefit of that coalition (worst_coal_benefit), and
            # (3) the minimum surplus of that coalition
            worst_coal_set, worst_coal_benefit, lower_problem_min_surplus = callback_worst_coalition(current_profit_dist)

            if worst_coal_set != Set(player_set)
                # specify that the profit of each subset of the group is better off with the grand coalition
                con_it = @constraint(
                    model_dist,
                    sum([profit_dist[pl] for pl in worst_coal_set]) >= worst_coal_benefit + min_surplus
                )
            else
                con_it = nothing
            end

            # data of the iteration
            iter_data = (
                current_profit=current_profit_dist,
                worst_coal=worst_coal_set,
                benefit_coal=worst_coal_benefit,
                value_min_surplus=value_min_surplus,
                lower_problem_min_surplus=lower_problem_min_surplus,
                constraint=con_it,
            )

            # update coalition
            added_coalition = string(worst_coal_set) * ": " * string(current_profit_dist.data)

            # add the iteration to the history
            push!(history, iter_data)
        end
    end

    # result of the profit distribution
    profit_distribution = Dict(zip(player_set, value.(profit_dist).data))
    min_surplus = value(min_surplus)

    if raw_outputs
        return profit_distribution, min_surplus, history, model_dist
    else
        return profit_distribution
    end

end


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
lower_bound : Number (optional, default 0.0)
    Lower bound of the variables of the problem (benefit distribution and margin of the worst coalition)
upper_bound : Number (optional, default nothing)
    Upper bound of the variables of the problem (benefit distribution and margin of the worst coalition)
    When nothing, the value is automatically set to the benefit of the grand coalition
verbose : Bool (optional, default true)
    When true, it shows a progress bar to describe the current execution status
raw_outputs : Bool (optional, default false)
    When true, it returns all raw outputs

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
        lower_bound=0.0,
        upper_bound=nothing,
        verbose=true,
        raw_outputs=false,
    )

    if verbose
        println("PHASE 1: Start first least core analysis\n")
    end

    profit_distribution, min_surplus, history, model_dist = least_core(
        mode, optimizer;
        rtol=rtol,
        atol=atol,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        verbose=verbose,
        raw_outputs=true
    )
    
    if verbose
        println("PHASE 1: First least core analysis ended with value $min_surplus\n ----")
        println("PHASE 2: Start second phase to indentify the desired objective\n")
    end

    # update objective function of the least_core model
    dist_objective(model_dist)

    # get JuMP variable of the profit distribution
    profit_dist = model_dist[:profit_dist]

    # initialization while condition
    continue_while = true
    iter = 0
    added_coalition = "Init"

    # printing formats
    format_print_head = "{:<15s} {:^15s} {:^15s} {:^15s} {:<s}"  # for the header
    format_print_iter = "{:<15s} {:> 13.2e} {:> 13.2e} {:> 13.2e} {:<s}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Current value", "Tolerance", "Added coalition : benefit allocation")
    end

    while continue_while

        # update counter
        iter += 1

        # optimize current model
        optimize!(model_dist)

        # get current profit distribution
        current_profit_dist = value.(profit_dist)

        # get (1) the coalition (worst_coal_set) with the least benefit,
        # (2) the total benefit of that coalition (worst_coal_benefit), and
        # (3) the minimum surplus of that coalition
        worst_coal_set, worst_coal_benefit, lower_problem_min_surplus = callback_worst_coalition(current_profit_dist)

        # if verbose print log
        if verbose
            rel_tol = compute_relative_tol(min_surplus, lower_problem_min_surplus)
            printfmtln(format_print_iter, iter, min_surplus, value_min_surplus, rel_tol, added_coalition)
        end
        
        con_it = nothing

        # check if convergence has been reached
        if isapprox(lower_problem_min_surplus, min_surplus, rtol=rtol, atol=atol, norm=abs)
            # convergence reached
            continue_while = false
        else
            # convergence not reached: add new constraint

            if worst_coal_set != Set(player_set)
                # specify that the profit of each subset of the group is better off with the grand coalition
                con_it = @constraint(
                    model_dist,
                    sum([profit_dist[pl] for pl in worst_coal_set]) >= worst_coal_benefit + min_surplus
                )
            end
        end
        

        # data of the iteration
        iter_data = (
            current_profit=current_profit_dist,
            worst_coal=worst_coal_set,
            benefit_coal=worst_coal_benefit,
            value_min_surplus=value_min_surplus,
            lower_problem_min_surplus=lower_problem_min_surplus,
            constraint=con_it,
        )

        # update coalition
        added_coalition = string(worst_coal_set) * ": " * string(current_profit_dist.data)

        # add the iteration to the history
        push!(history, iter_data)
    end

    # result of the profit distribution
    profit_distribution = Dict(zip(player_set, value.(profit_dist).data))
    min_surplus = value(min_surplus)

    if raw_outputs
        return profit_distribution, min_surplus, history, model_dist
    else
        return profit_distribution
    end

end