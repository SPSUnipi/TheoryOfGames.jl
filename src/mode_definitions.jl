abstract type AbstractCalcMode end  # abstract type to calculate a quantity

"""
    EnumMode is a structure defining the enumerative representation of the
    utilities of a Game
"""
struct EnumMode <: AbstractCalcMode  # Enumerative technique
    player_set
    utilities
    

    function EnumMode(player_set_::AbstractVector=[], utilities_::AbstractDict=Dict())
        return new(player_set_, utilities_)
    end

    function EnumMode(::Any, ::Any)
        throw(ArgumentError("Invalid arguments for EnumMode"))
    end
end

function EnumMode(player_set, utility::Function; verbose=true)
    return EnumMode(player_set, utility_combs(player_set, utility; verbose=verbose))
end


struct RobustMode <: AbstractCalcMode end  # Robust-Optimization technique