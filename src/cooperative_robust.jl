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

Outputs
------
profit_distribution
    Benefit distribution by player
min_surplus
    Benefit of the coalition with the minimum surplus
history
"""
function least_core(
        mode::RobustMode,
        optimizer;
        rtol=1e-6,
        atol=1e-6, 
        lower_bound=0.0,
        upper_bound=nothing,
        verbose=true,
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

    # printing formats
    format_print_head = "{:<12s} {:^12s} {:^12s} {:^12s}"  # for the header
    format_print_iter = "{:<12s} {:> 10.2e} {:> 10.2e} {:> 10.2e}"  # for the iterations

    # if verbose, print header
    if verbose
        printfmtln(format_print_head, "Iteration", "Upper bound", "Current value", "Tolerance")
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
            printfmtln(format_print_iter, iter, lower_problem_min_surplus, value_min_surplus, rel_tol)
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

            # specify that the profit of each subset of the group is better off with the grand coalition
            con_it = @constraint(
                model_dist,
                sum([profit_dist[pl] for pl in worst_coal_set]) >= worst_coal_benefit + min_surplus
            )

            # data of the iteration
            iter_data = (
                current_profit=current_profit_dist,
                worst_coal=worst_coal_set,
                benefit_coal=worst_coal_benefit,
                value_min_surplus=value_min_surplus,
                lower_problem_min_surplus=lower_problem_min_surplus,
                constraint=con_it,
            )

            # add the iteration to the history
            push!(history, iter_data)
        end
    end

    # result of the profit distribution
    profit_distribution = value.(profit_dist)
    min_surplus = value(min_surplus)

    return profit_distribution, min_surplus, history
end