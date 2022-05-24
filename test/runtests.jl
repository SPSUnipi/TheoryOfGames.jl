using Games
using Test
using YAML
# using Gurobi
using Ipopt
using JuMP
using HiGHS
using FileIO
using JLD2

include("tests.jl")

# constant for the optimization and testing
const ATOL = 1e-4
const OPTIMIZER = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0, "tol"=>1e-6)
# optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=> 0)  #, "tol"=>1e-4)  #, "NonConvex"=>2)

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

"Function to create IterMode out of an example"
function to_IterMode(example)
    util_combs = utility_combs(example.player_set, example.utility)
    keys_no_grand_coalition = setdiff(keys(util_combs), [Set(), Set(example.player_set)])

    let util_combs=util_combs, keys_no_grand_coalition=keys_no_grand_coalition
        
        callback_benefit_by_coalition = (coal)->util_combs[Set(coal)]

        function callback_worst_coalition(profit_dist)
            min_surplus_combs = Dict(
                comb=>sum(Float64[profit_dist[c] for c in comb]) - util_combs[Set(comb)]
                for comb in keys_no_grand_coalition
            )
            min_surplus, least_benefit_coal = findmin(min_surplus_combs)
            return least_benefit_coal, util_combs[Set(least_benefit_coal)], min_surplus
        end

        return Games.IterMode(example.player_set, callback_benefit_by_coalition, callback_worst_coalition)
    end
end

"Function to create EnumMode out of an example"
to_EnumMode(example) = EnumMode(example.player_set, example.utility)

example_list = [
    Examples.three_users_onlyone,
    Examples.three_users_atleasttwo,
    Examples.three_users_atleasttwo_string,
    Examples.three_users_mapping,
]


@testset "Game tests - EnumMode" begin

    @testset "shapley" begin
        println("TEST SET - SHAPLEY")
        for example in example_list
            test_example("ENUM_" * example.name, shapley_value, to_EnumMode(example))
        end
    end

    @testset "least_core" begin
        println("TEST SET - LEAST CORE")
        for example in example_list
            test_example("ENUM_" * example.name, least_core, to_EnumMode(example), OPTIMIZER)
        end
    end

    @testset "nucleolus" begin
        println("TEST SET - NUCLEOLUS")
        for example in example_list
            test_example("ENUM_" * example.name, nucleolus, to_EnumMode(example), OPTIMIZER)
        end
    end

    @testset "in_core" begin
        println("TEST SET - IN CORE")
        for example in example_list
            test_example("ENUM_" * example.name, in_core, to_EnumMode(example), OPTIMIZER)
        end
    end

    @testset "verify_in_core" begin
        println("TEST SET - VERIFY IN CORE")
        for example in example_list

            player_set = example.player_set
            example_mode = to_EnumMode(example)

            # obtain an in-core solution
            val_dist = in_core(example_mode, OPTIMIZER)

            # test that the in-core solution is actually recognized in the core
            @test verify_in_core(val_dist, example_mode, OPTIMIZER) == true

            # create an artificial solution likely not to be in the core
            equal_vals = JuMP.Containers.DenseAxisArray(
                fill(sum(values(val_dist))/length(player_set), length(player_set)), player_set
            )

            # test that the artificial distribution does not belong to the core
            @test verify_in_core(equal_vals, example_mode, OPTIMIZER) == false
        end
    end

    @testset "var_core" begin
        println("TEST SET - VAR CORE")
        for example in example_list
            test_example("ENUM_" * example.name, var_core, to_EnumMode(example), OPTIMIZER)
        end
    end

    @testset "ref_in_core" begin
        println("TEST SET - REF IN CORE")
        for example in example_list
            # reference distribution
            ref_dist = Dict(
                zip(example.player_set, fill(0.0, length(example.player_set)))
            )
            test_example("ENUM_" * example.name, ref_in_core, to_EnumMode(example), ref_dist, OPTIMIZER)
        end
    end

    @testset "var_least_core" begin
        println("TEST SET - VAR LEAST CORE")
        for example in example_list
            test_example("ENUM_" * example.name, var_least_core, to_EnumMode(example), OPTIMIZER)
        end
    end

    @testset "ref_least_core" begin
        println("TEST SET - REF LEAST CORE")
        for example in example_list
            # reference distribution
            ref_dist = Dict(
                zip(example.player_set, fill(0.0, length(example.player_set)))
            )
            test_example("ENUM_" * example.name, ref_least_core, to_EnumMode(example), ref_dist, OPTIMIZER)
        end
    end

end

@testset "Games tests - IterMode" begin
    
    @testset "least_core" begin
        println("TEST SET - ROBUST LEAST CORE")
        for example in example_list

            # obtain output from robust least core
            result_value = least_core(to_IterMode(example), OPTIMIZER)

            # test that the solution belongs to the core
            @test verify_in_core(result_value, to_EnumMode(example), OPTIMIZER) == true
        end
    end

    @testset "var_least_core" begin
        println("TEST SET - VAR LEAST CORE")
        for example in example_list
            test_example("ITER_" * example.name, var_least_core, to_IterMode(example), OPTIMIZER)
        end
    end

    @testset "ref_least_core" begin
        println("TEST SET - REF LEAST CORE")
        for example in example_list
            # reference distribution
            ref_dist = Dict(
                zip(example.player_set, fill(0.0, length(example.player_set)))
            )
            test_example("ITER_" * example.name, ref_least_core, to_IterMode(example), ref_dist, OPTIMIZER)
        end
    end

    @testset "in_core" begin
        println("TEST SET - IN CORE")
        for example in example_list
            test_example("ITER_" * example.name, in_core, to_IterMode(example), OPTIMIZER)
        end
    end

    @testset "verify_in_core" begin
        println("TEST SET - VERIFY IN CORE")
        for example in example_list

            player_set = example.player_set
            example_mode = to_IterMode(example)

            # obtain an in-core solution
            val_dist = in_core(example_mode, OPTIMIZER)

            # test that the in-core solution is actually recognized in the core
            @test verify_in_core(val_dist, example_mode, OPTIMIZER) == true

            # create an artificial solution likely not to be in the core
            equal_vals = JuMP.Containers.DenseAxisArray(
                fill(sum(values(val_dist))/length(player_set), length(player_set)), player_set
            )

            # test that the artificial distribution does not belong to the core
            @test verify_in_core(equal_vals, example_mode, OPTIMIZER) == false
        end
    end

    @testset "var_core" begin
        println("TEST SET - VAR CORE")
        for example in example_list
            test_example("ITER_" * example.name, var_core, to_IterMode(example), OPTIMIZER)
        end
    end

    @testset "ref_in_core" begin
        println("TEST SET - REF IN CORE")
        for example in example_list
            # reference distribution
            ref_dist = Dict(
                zip(example.player_set, fill(0.0, length(example.player_set)))
            )
            test_example("ITER_" * example.name, ref_in_core, to_IterMode(example), ref_dist, OPTIMIZER)
        end
    end

end

@testset "IO tests" begin
    
    println("TEST IO")
    for example in example_list

        enum_mode = to_EnumMode(example)
    
        path_solution = (
            string(@__DIR__) * 
            "\\testcases\\io\\" * 
            string(typeof(enum_mode)) * "\\" 
            * example.name * ".jld2"
        )
    
        # obtain output from robust least core
        @test_nowarn save(path_solution, enum_mode)

        # test that the solution belongs to the core
        loaded_mode = load(path_solution, EnumMode())

        @test enum_mode.player_set == loaded_mode.player_set
        @test enum_mode.utilities == loaded_mode.utilities
    end

end