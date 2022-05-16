

# """
# load(output_file::AbstractString)
# Function to save the results and the model to the hard drive
# """
# function FileIO.load(output_file::AbstractString)
#     return load!(output_file, ModelEC())
# end

"""
    save(output_file::AbstractString, game_mode<:AbstractCalcMode)
Function to save the Mode element of Games.jl
"""
function FileIO.save(output_file::AbstractString, game_mode::AbstractCalcMode)
    f_names = fieldnames(typeof(game_mode))
    save_model = Dict([
        "mode"=>string(typeof(game_mode));
        string.(f_names) .=> [getfield(game_mode, field) for field in f_names]
    ])
    FileIO.save(output_file, save_model)
end

"""
    load(output_file::AbstractString)
Function to load the results and the model to the hard drive
"""
function FileIO.load(output_file::AbstractString, game_mode::AbstractCalcMode)

    mode_type = typeof(game_mode)

    data = FileIO.load(output_file)

    if data["mode"] != string(typeof(game_mode))
        return throw(ArgumentError("Input file mode does not fit the requested mode"))
    else
        f_names = fieldnames(typeof(game_mode))
        args = [data[string(f_name)] for f_name in f_names]
        ret_mode = mode_type(args...)
        return ret_mode
    end
end
