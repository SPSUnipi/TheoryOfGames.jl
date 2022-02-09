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
function shapley_value(player_set::Vector, utility::Function, mode::AbstractCalcMode=EnumMode(); verbose=true)
    
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