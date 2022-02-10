module Examples

    struct Example
        name
        player_set
        utility
    end

    three_users_onlyone = Example(
        "three_users_onlyone",  # name of example
        [1, 2, 3],  # player_set
        x-> 1 in x ? 1.0 : 0.0,  # function
    )

    three_users_atleasttwo = Example(  # https://en.wikipedia.org/wiki/Shapley_value
        "three_users_atleasttwo",  # name of example
        [1, 2, 3],  # player_set
        x-> (1 in x) && length(x) > 1 ? 1.0 : 0.0,  # function
    )

    map_three_users_mapping = Dict( # https://math.stackexchange.com/questions/1108449/finding-the-nucleolus
        Set([])=>0.0,
        Set([1])=>2.,
        Set([2])=>5.,
        Set([3])=>4.,
        Set([1,2])=>14.,
        Set([1,3])=>18.,
        Set([2,3])=>9.,
        Set([1,2,3])=>24.,
    )

    three_users_mapping = Example(  # https://math.stackexchange.com/questions/1108449/finding-the-nucleolus
        "three_users_atleasttwo",  # name of example
        [1, 2, 3],  # player_set
        x->map_three_users_mapping[Set(x)],  # function
    )

end