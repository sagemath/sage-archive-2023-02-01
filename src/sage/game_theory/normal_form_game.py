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

from gambit.lib.libgambit import *


def matrix_game(*args):
    num_players = len(args)

    #create a test for all matrices being the same shape.

    numrows = args[0].nrows()
    numcols = args[0].ncols()

    game_array =
