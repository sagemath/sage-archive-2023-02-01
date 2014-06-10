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

from sage.misc.package import is_package_installed
if is_package_installed('gambit') is not True:
            raise NotImplementedError("Optional package Gambit is not installed")

from gambit.lib.libgambit import new_table
from sage.matrix.constructor import matrix


def two_matrix_game(matrix1, matrix2):

    #create a test for all matrices being the same shape.
    p1_strats = len(matrix1.rows())
    p2_strats = len(matrix1.columns())

    game = new_table([p1_strats, p2_strats])

    for i in range(p1_strats):
        for j in range(p2_strats):
            game[i, j][0] = int(matrix1[i, j])
            game[i, j][1] = int(matrix2[i, j])

    return game


def test_game():
    a = matrix([[1, 2], [3, 4]])
    b = matrix([[3, 3], [1, 4]])
    return two_matrix_game(a, b)
