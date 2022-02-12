using Games
using Test
using YAML
using Gurobi
using JuMP

include("tests.jl")

"Function to test examples"
function test_example(example, testing_function, args...; kwargs...)
    # calculate simulations
    calc_solution = testing_function(example.player_set, example.utility, args...; kwargs...)

    path_solution = (string(@__DIR__) * "\\testcases\\" * string(testing_function) * "\\" * example.name * ".yml")
    
    if isfile(path_solution)
        # if the file exists run tests
        proven_solution = YAML.load_file(path_solution)
        
        @test Set(keys(calc_solution)) == Set(keys(proven_solution))
        @test all(isapprox(calc_solution[k] - proven_solution[k], 0.0, atol=1e-4) for k in keys(proven_solution))
    else
        # otherwise create the tests
        mkpath(dirname(path_solution))

        YAML.write_file(path_solution, calc_solution)
        @warn("Preloaded solution not found, then it has been created")
    end

end

example_list = [
    Examples.three_users_onlyone,
    Examples.three_users_atleasttwo,
    Examples.three_users_mapping
]

optimizer = optimizer_with_attributes(Gurobi.Optimizer)

functions_to_test = [least_core, nucleolus, in_core, var_core]


@testset "Game tests" begin

    @testset "shapley" begin
        for example in example_list
            test_example(example, shapley_value)
        end
    end

    for f in functions_to_test
        @testset :(string(f) begin
            for example in example_list
                test_example(example, f, optimizer)
            end
        end
    end
end