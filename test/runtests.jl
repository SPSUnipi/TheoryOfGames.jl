using Games
using Test
using YAML
using Gurobi
using Ipopt
using JuMP
using HiGHS

include("tests.jl")

"Function to test examples"
function test_example(example_name, testing_function, args...; kwargs...)
    # calculate simulations
    calc_solution = testing_function(args...; kwargs...)

    path_solution = (string(@__DIR__) * "\\testcases\\" * string(testing_function) * "\\" * example_name * ".yml")
    
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

enum_example_list = [
    (example.name, EnumMode(example.player_set, example.utility))
        for example in example_list
]

optimizer = optimizer_with_attributes(Gurobi.Optimizer)  # , "tol"=>1e-4)  #, "NonConvex"=>2)


@testset "Game tests" begin


    @testset "shapley" begin
        for (example_name, example_mode) in enum_example_list
            test_example(example_name, shapley_value, example_mode)
        end
    end

    @testset "least_core" begin
        for (example_name, example_mode) in enum_example_list
            test_example(example_name, least_core, example_mode, optimizer)
        end
    end

    @testset "nucleolus" begin
        for (example_name, example_mode) in enum_example_list
            test_example(example_name, nucleolus, example_mode, optimizer)
        end
    end

    @testset "in_core" begin
        for (example_name, example_mode) in enum_example_list
            test_example(example_name, in_core, example_mode, optimizer)
        end
    end

    @testset "var_core" begin
        for (example_name, example_mode) in enum_example_list
            test_example(example_name, var_core, example_mode, optimizer)
        end
    end

    @testset "var_least_core" begin
        for (example_name, example_mode) in enum_example_list
            test_example(example_name, var_least_core, example_mode, optimizer)
        end
    end

    @testset "ref_least_core" begin
        for (example_name, example_mode) in enum_example_list
            # reference distribution
            ref_dist = Dict(
                zip(example_mode.player_set, fill(0.0, length(example_mode.player_set)))
            )
            test_example(example_name, ref_least_core, example_mode, ref_dist, optimizer)
        end
    end

    @testset "ref_in_core" begin
        for (example_name, example_mode) in enum_example_list
            # reference distribution
            ref_dist = Dict(
                zip(example_mode.player_set, fill(0.0, length(example_mode.player_set)))
            )
            test_example(example_name, ref_least_core, example_mode, ref_dist, optimizer)
        end
    end

end