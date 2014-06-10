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

from gambit.lib.libgambit import Game


class NormalFormGame(Game):
    r"""
    An object representing a normal form game.

    INPUT:

    - Bi matrix: can be input as two matrices, an array of two dimensional lists OR a single matrix in which case it will be assumed that the game is zero sum.
    - Strategy names - default = ``False``, if a list of names is passed then all strategy vectors will be returned as dictionaries...

    EXAMPLES::
    Basic example of how to initiate a normal form game using two matrices. ::

        sage: A = matrix([[2,0],
        ....:             [5,1]])
        sage: B = matrix([[2,5],
        ....:             [0,1]])
        sage: game = NormalFormGame(A, B)
        sage: game.show()

        sage: D = matrix([[2, 5],
        ....:             [1, 0]])
        sage: game = NormalFormGame(D)
        sage: game.show()
    """

    def __init__(self):
        r"""
        """
        self.game = gambit.new_table([matrix1.nrows(), matrix1.ncols()])
        self.player1 = player1
        self.player2 = player2
        self.game.players[0].label = player1
        self.game.players[1].label = player2

        for row in range(matrix1.nrows()):
            for col in range(matrix1.ncols()):
                a = matrix1[row,col]
                b = matrix2[row,col]
                self.game[row, col][0] = 1
                self.game[row, col][1] = 2

    def _repr_(self):
        print self.game

    def best_responses(self):
        """
        Maybe return a dictionary for each player...
        """

    def dominated_strategies(self):
        """
        Identify any dominated strategies (dominated by pure but also by mixed, that's harder... I'm not entirely sure how to do it, but there will be a way somewhere...)
        """

    def nash_equilibrium(self):
        """
        Compute all equilibrium
        """
        solver = gambit.nash.ExternalEnumPureSolver()
        return solver.solver(self.game)
