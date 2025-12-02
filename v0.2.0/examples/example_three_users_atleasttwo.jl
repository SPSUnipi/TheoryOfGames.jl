# # Example 1: Key player but not alone
# This example aims to describe the simple use of the package for a 3-user case.
# The major reward allocation techniques supported by the package are considered,
# using the enumerative technique.

# The sample game considered in this example has the following characteristics:
# - there are 3 players: 1, 2 and 3
# - each user can join or not the coalition
# - the utility value is:
#   - always 0 if player 1 does not join the coalition
#   - +1 if player 1 joins the coalition but it is not alone

# In the following, the game is constructed.

# First, the packages are imported
using TheoryOfGames
using JuMP, Ipopt

# ## Initialization of the game

# Define the set of the users that can join or not the coalition
player_set = [1, 2, 3]

# Define the utility function that is a function that:
# - takes as input an iterable that represents a coalition
# - it returns a value in agreement to what discussed above and here summarized:
#   - it returns +1 if player 1 joins the coalition but it is not alone
#     ``(1 in x) && length(x)``
#   - it returns 0 if player 1 does not join the coalition (otherwise)
#     ``else condition`
utility = x->((1 in x) && (length(x) > 1)) ? 1.0 : 0.0

# Show the value of utility for the coalition where only user 1 joins
utility([1])

# Show the value if also user 2 joins
utility([1, 2])

# ## Calculation of selected reward allocation functions

# Define the Enumerative object
enum_obj = EnumMode(player_set, utility)

# Calculate shapley value
shapley_value(enum_obj)

# Define the optimizer needed for nucleolus and var least core techniques.
# We use the default Ipopt and disable the output for simplicity (``print_level=>0``).
OPTIMIZER = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# Calculate nucleolus
nucleolus(enum_obj, OPTIMIZER)

# Calculate var least core
var_least_core(enum_obj, OPTIMIZER)