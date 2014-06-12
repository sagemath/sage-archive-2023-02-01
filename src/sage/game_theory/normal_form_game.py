r"""
2 Player normal form games

This module implements 2 by 2 normal form (bi-matrix) games. A variety of operations on these games can be carried out:

- Identification of (weakly) dominated strategies;
- Identification of Best responses to a given strategy;
- Identification of Nash Equilibrium (this is done by interfacing with Gambit);
"""
from itertools import product
from sage.misc.package import is_package_installed
from gambit import Game
from gambit.nash import ExternalLCPSolver


class NormalFormGame(Game):
    def __new__(NormalFormGame, matrix1=False, matrix2=False, game=False):
        if matrix1 and matrix2:
            g = Game.new_table([len(matrix1.rows()), len(matrix2.rows())])
            g.__class__ = NormalFormGame
            return g
        if game:
            h = game
            h.__class__ = NormalFormGame
            return h

    def __init__(self, matrix1=False, matrix2=False, game=False ):

        self.matrix1 = matrix1
        self.matrix2 = matrix2
        if matrix1 and matrix2:
            p1_strats = range(len(matrix1.rows()))
            p2_strats = range(len(matrix1.columns()))

            for k in product(p1_strats, p2_strats):
                    self[k][0] = int(matrix1[k])
                    self[k][1] = int(matrix2[k])

    def obtain_Nash(self, algorithm="gambit"):
        r"""
        A function to return the Nash equilibrium for a game.
        Optional arguments can be used to specify the algorithm used.
        If no algorithm is passed then an attempt is made to use the most appropriate algorithm.

        INPUT:

        - ``game`` -- a gambit game object.

        - ``algorithm`` -- the following algorithms should be available through this function:

          - lrsnash (``algorithm="lrsnash"``). This algorithm is only suited for 2 player games. See the [insert website here] web site.

          - gambit (``algorithm="gambit"``). This algorithm is only suited for 2 player games. See the [insert website here] web site. NOTE THAT WE NEED TO GET THE ACTUAL NAME OF THE GAMBIT ALGORITHM

          - support enumeration (``algorithm="supportenum"``). This is a very inefficient algorithm (in essence a brute force approach).

          - If ``algorithm=False`` (default), a default solver dependent on the size of the game will be used.

        - ``maximization``

           - When set to ``True`` (default) it is assumed that players aim to maximise their utility.
           - When set to ``False`` (default) it is assumed that players aim to minimise their utility.
        """

        if algorithm == "gambit":
            solver = ExternalLCPSolver()
            return solver.solve(self)

        if not algorithm:
            if len(game.players) > 2:
                raise NotImplementedError("Nash equilibrium for games with more than 2 players have not been implemented yet. Please see the gambit website [LINK] that has a variety of available algorithms.")
            algorithm = "lrsnash"  # This will do for now: when we have a variety of algorithms we will test the best and set the best one.
        if algorithm == "lrsnash":
            if not is_package_installed('lrs'):  # This is different to how you've used above, this is better unless there are Sage common practices that point to your way.
                raise NotImplementedError("lrs is not installed")
            # Run lrs (Have something simple that takes opposite of matrices as input to lrs)
            return vector(1/3,1/3,1/3),vector(1/2,1/2)
        raise NotImplementedError("%s is not yet implemented" % algorithm)


    def game_matrix(game):
       """
       Creates 2 Matrices based on a Gambit Game object.
       """
       # check that it is only two players.


def test_game():
    from sage.matrix.constructor import matrix
    a = matrix([[1, 2], [3, 4]])
    b = matrix([[3, 3], [1, 4]])
    return NormalFormGame(matrix1=a, matrix2=b)



