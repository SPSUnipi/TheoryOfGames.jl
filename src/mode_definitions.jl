abstract type AbstractCalcMode end  # abstract type to calculate a quantity

"""
    EnumMode is a structure defining the enumerative representation of the
    utilities of a Game
"""
struct EnumMode <: AbstractCalcMode  # Enumerative technique
    player_set  # AbstractVector of the players
    utilities  # AbstractDict of utility functions

    function EnumMode(player_set_=[], utilities_::AbstractDict=Dict())
        return new(player_set_, utilities_)
    end

    function EnumMode(player_set, utility::Function; verbose=true)
        return new(player_set, utility_combs(player_set, utility; verbose=verbose))
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
    RobustMode is a structure defining the modality for robust-like identification
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
    and the total benefit of the coalition
"""
struct RobustMode <: AbstractCalcMode  # Robust-Optimization technique
    # AbstractVector of the players
    player_set
    # Callback function used to obtain the benefit of a coalition
    callback_benefit_by_coalition::Function
    # Callback function used in the iterative approaches
    callback_worst_coalition::Function

    function RobustMode(player_set_, callback_benefit_by_coalition_::Function, callback_worst_coalition_::Function)
        return new(player_set_, callback_benefit_by_coalition_, callback_worst_coalition_)
    end

    # function EnumMode(input_file::String)
    #     # load YAML file
    #     data = YAML.load_file(file_name)
    #     player_set_ = get(data, "user_set", [])
    #     utilities_ = get(data, "user_set", Dict([]=>0.0))

    #     return new(player_set_, utilities_)
    # end

    function RobustMode(::Any, ::Any)
        throw(ArgumentError("Invalid arguments for EnumMode"))
    end
end