module Examples

    struct Example
        name
        player_set
        utility
    end

    three_users_zeromargin = Example(
        "three_users_zeromargin",  # name of example
        [1, 2, 3],  # player_set
        x-> 1 in x ? 1.0 : 0.0,  # function
    )

    three_users = Example(
        "three_users",  # name of example
        [1, 2, 3],  # player_set
        x-> 1 in x ? 1.0 : 0.0,  # function
    )

end