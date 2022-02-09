using Games
using Test


@testset "Game tests" begin

    @testset "shapley_0" begin

        utility(x) = 1 in x ? 1.0 : 0.0
        player_set = [1, 2, 3]
        
        sh_vals = shapley_value(player_set, utility)

        @test sh_vals == Dict(1=>1.0, 2=>0.0, 3=>0.0)
    end

    @testset "shapley_1" begin

        utility(x) = 1 in x && length(x) >=2 ? 1.0 : 0.0
        player_set = [1, 2, 3]
        
        sh_vals = shapley_value(player_set, utility)

        sh_vals_sol = Dict(1=>0.66666666, 2=>0.166666666, 3=>0.166666666)

        @test Set(keys(sh_vals)) == Set(keys(sh_vals_sol))
        @test all(isapprox(sh_vals[k] - sh_vals_sol[k], 0.0, atol=1e-6) for k in keys(sh_vals))

    end

end