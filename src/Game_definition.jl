abstract type AbstractGame end

mutable struct Game{T} <: AbstractGame
    player_set::Vector{T}
    utility::Function
end