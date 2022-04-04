"Return the empty set of the coalition"
empty_set(player_set) = Set(similar(player_set, 0))

"Get the total number of coalitions"
number_coalitions(player_set) = sum(binomial(length(player_set), k) for k=1:length(player_set))

"Function to get the types of the utility arguments and outputs"
function utility_io_types(player_set, utility::Function)
    # identify return type of utility for the example of an empty coalition
    empty_coalition = empty_set(player_set)
    empty_val = utility(empty_coalition)

    return Set{eltype(player_set)}, typeof(empty_val)
end


"""
    utility_combs(player_set, utility)

Function to calculate the utility for every combination of players.
This table may be used in Game Theory, such as in the calculation of the shapley value.
The function iterates all combinations of players in the player_set and execute the utility function
to identify the benefits to be shared.

Inputs
------
player_set : Vector
    Vector of the players
utility : Function
    Utility function that given any coalition returns the benefit of the coalition
    It shall be a function utility(::Vector)::T<:Number
verbose : Bool
    When true, it shows a progress bar to describe the current execution status

Outputs
-------
utilities : Dict
    Dictionary that specifies the utility of each combination of coalition in player_set
"""
function utility_combs(player_set, utility::Function; verbose=true)

    # get types of player_set and utility output
    empty_coalition = empty_set(player_set)
    empty_val = utility(empty_coalition)

    ptype = Set{eltype(empty_coalition)}
    utype = typeof(empty_val)

    # initialize return dictionary
    dict_ret = Dict{ptype, utype}(empty_coalition=>empty_val)

    # get players combinations
    combs = combinations(player_set)

    np = length(player_set)

    # number of combinations
    n_combs = number_coalitions(player_set)

    if verbose
        for comb in ProgressBar(Set.(combs), total=n_combs)
            dict_ret[comb] = utility(comb)
        end
    else
        for comb in Set.(combs)
            dict_ret[comb] = utility(comb)
        end
    end

    return dict_ret
end


"""

Auxiliary function to compute the relative tolerance between two quantities x and y

"""
compute_relative_tol(x, y, norm=abs, atol=1e-6) = norm(x-y) / (atol + max(x, y))