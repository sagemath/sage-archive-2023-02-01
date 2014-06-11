r"""
2 Player normal form games

This module implements characteristic function cooperative games.
The main contribution is a class for a characteristic function game.
Methods to calculate the Shapley value (a fair way of sharing common
resources: https://www.youtube.com/watch?v=aThG4YAFErw) as well as
test properties of the game (monotonicity, super additivity) are also included.
This module implements 2 by 2 normal form (bi-matrix) games. A variety of operations on these games can be carried out:

- Identification of (weakly) dominated strategies;
- Identification of Best responses to a given strategy;
- Identification of Nash Equilibrium (this is done by interfacing with Gambit);
"""
from itertools import product
from sage.misc.package import is_package_installed


def normal_form_game(matrix1, matrix2):
    """
    Creates an instance of the Gambit Game class from
    2 matrices.
    """
    if is_package_installed('gambit') is not True:
        raise NotImplementedError("Optional package Gambit is not installed")
    else:
        from gambit.lib.libgambit import new_table

    #create a test for all matrices being the same shape.
    p1_strats = range(len(matrix1.rows()))
    p2_strats = range(len(matrix1.columns()))

    game = new_table([len(p1_strats), len(p2_strats)])

    for k in product(p1_strats, p2_strats):
            game[k][0] = int(matrix1[k])
            game[k][1] = int(matrix2[k])

    return game


def game_matrix(game):
    """
    Creates 2 Matrices based on a Gambit Game object.
    """
    # check that it is only two players.


def test_game():
    from sage.matrix.constructor import matrix
    a = matrix([[1, 2], [3, 4]])
    b = matrix([[3, 3], [1, 4]])
    return two_matrix_game(a, b)


