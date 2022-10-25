abstract type AbstractCalcMode end  # abstract type to calculate a quantity

"""
    EnumMode is a structure defining the enumerative representation of the
    utilities of a Game
"""
struct EnumMode <: AbstractCalcMode  # Enumerative technique
    player_set::Union{Set, Vector}  # AbstractVector of the players
    utilities  # AbstractDict of utility functions

    function EnumMode(player_set_::Union{Set, Vector}=[], utilities_::AbstractDict=Dict())
        return new(player_set_, utilities_)
    end

    function EnumMode(player_set::Union{Set, Vector}, utility::Function; verbose=true, kwargs...)
        return new(player_set, utility_combs(player_set, utility; verbose=verbose, kwargs...))
    end

    function EnumMode(input_file::String)
        # load YAML file
        data = YAML.load_file(file_name)
        player_set_ = get(data, "user_set", [])
        utilities_ = get(data, "user_set", Dict([]=>0.0))

        return new(player_set_, utilities_)
    end

    function EnumMode(::Any, ::Any)
        throw(ArgumentError("Invalid arguments for EnumMode"))
    end
end

"""
    IterMode is a structure defining the modality for interative identification
    of the benefit distribution mechanism

Fields
------
player_set
    Vector of players
callback_benefit_by_coalition : Function
    Callback function that is used to determine what is the benefit of a coalition
    and the total benefit of the coalition
callback_worst_coalition : Function
    Callback function that is used to determine what is the coalition with the worst benefit
    and the total benefit of the coalition.
    Arguments shall be:
    - (mandatory) a Dict-like container to describe the current reward distribution by user
    - kwargs (optional):
        - modify_solver_options: vector of pairs to iteratively change the optimization callback
            to speed up computation time. The implemented callbacks are:
            1. best_obj_stop: when the best_obj_stop callback is used in iterative methods,
                the corresponding option is updated here depending on the value
            2. lower_relaxation_stop_option: option to stop the execution as a lower bound is reached
    The function shall return a Vector of NamedTuple, where for every entry o shall contain:
    - least_profitable_coalition: members of the worst coalition, for result o
    - coalition_benefit: benefit of the coalition, for result o
    - min_surplus: minimum surplus of the coalition, for result o
"""
struct IterMode <: AbstractCalcMode  # Robust-Optimization technique
    # AbstractVector of the players
    player_set::Union{Set, Vector}
    # Callback function used to obtain the benefit of a coalition
    callback_benefit_by_coalition::Function
    # Callback function used in the iterative approaches
    callback_worst_coalition::Function

    function IterMode(player_set_::Union{Set, Vector}, callback_benefit_by_coalition_::Function, callback_worst_coalition_::Function)
        return new(player_set_, callback_benefit_by_coalition_, callback_worst_coalition_)
    end

    # function EnumMode(input_file::String)
    #     # load YAML file
    #     data = YAML.load_file(file_name)
    #     player_set_ = get(data, "user_set", [])
    #     utilities_ = get(data, "user_set", Dict([]=>0.0))

    #     return new(player_set_, utilities_)
    # end

    function IterMode(::Any, ::Any)
        throw(ArgumentError("Invalid arguments for EnumMode"))
    end
end