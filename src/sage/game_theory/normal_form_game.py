r"""
2 Player normal form games

This module implements 2 by 2 normal form (bi-matrix) games. A variety of
operations on these games can be carried out:

- Identification of (weakly) dominated strategies;
- Identification of Best responses to a given strategy;
- Identification of Nash Equilibrium (this is done by interfacing with Gambit);
"""
from itertools import product
from sage.misc.package import is_package_installed
from gambit import Game
from gambit.nash import ExternalLCPSolver
from sage.matrix.constructor import matrix, copy
from subprocess import Popen, PIPE, call
import random
import string
from os import remove

class NormalFormGame(Game):
    r"""
    An object representing a Normal Form Game. Primarily used to compute the
    Nash Equilibrium.

    INPUT:

    If 2 matrices are passed in the list then the corresponding Normal
    Form Game is produced.
    If a single matrix is passed in the list then the corresponding zero
    sum game is produced.
    If a gambit game object is passed then the corresponding game is created.

    - ``payoff_matrices`` - a list of payoff matrices defining a
                            Normal Form game.
    - ``game`` - an instance of gambit.Game().

    EXAMPLES:

    A basic 2-player game constructed from matrices. ::

        sage: A = matrix([[1, 2], [3, 4]])
        sage: B = matrix([[3, 3], [1, 4]])
        sage: C = NormalFormGame([A, B])

    This can be given a title and the players can be named. ::

        sage: C.title = "Simple Game"
        sage: C.players[int(0)].label = "James"
        sage: C.players[int(1)].label = "Vince"
        sage: C
        NFG 1 R "Simple Game" { "James" "Vince" }
        <BLANKLINE>
        { { "1" "2" }
        { "1" "2" }
        }
        ""
        <BLANKLINE>
        {
        { "" 1, 3 }
        { "" 3, 1 }
        { "" 2, 3 }
        { "" 4, 4 }
        }
        1 2 3 4
        <BLANKLINE>

    We can also pass a Gambit game and create it manually.
    (Taken from [GAMBIT WEBSITE]) ::

        sage: gambitgame= Game.new_table([2, 2])
        sage: gambitgame[int(0), int(0)][int(0)] = int(8)
        sage: gambitgame[int(0), int(0)][int(1)] = int(8)
        sage: gambitgame[int(0), int(1)][int(0)] = int(2)
        sage: gambitgame[int(0), int(1)][int(1)] = int(10)
        sage: gambitgame[int(1), int(0)][int(0)] = int(10)
        sage: gambitgame[int(1), int(0)][int(1)] = int(2)
        sage: gambitgame[int(1), int(1)][int(0)] = int(5)
        sage: gambitgame[int(1), int(1)][int(1)] = int(5)
        sage: gambitgame.title = "A prisoner's dilemma game"
        sage: gambitgame.players[int(0)].label = "Alphonse"
        sage: gambitgame.players[int(1)].label = "Gaston"
        sage: g = NormalFormGame(game=gambitgame)
        sage: g
        NFG 1 R "A prisoner's dilemma game" { "Alphonse" "Gaston" }
        <BLANKLINE>
        { { "1" "2" }
        { "1" "2" }
        }
        ""
        <BLANKLINE>
        {
        { "" 8, 8 }
        { "" 10, 2 }
        { "" 2, 10 }
        { "" 5, 5 }
        }
        1 2 3 4
        <BLANKLINE>


    This can be solved using ``obtain_Nash``. ::

        sage: g.obtain_Nash()
        [[[0.0, 1.0], [0.0, 1.0]]]

    Here is an example of a 3 by 2 game ::

        sage: A = matrix([[3,3],
        ....:             [2,5],
        ....:             [0,6]])
        sage: B = matrix([[3,2],
        ....:             [2,6],
        ....:             [3,1]])
        sage: g = NormalFormGame(A,B)

    This particular game has 3 Nash equilibrium ::

        sage: g.obtain_Nash()
        [[[1.0, 0.0, 0.0], [1.0, 0.0]],
         [[0.8, 0.2, 0.0], [0.6666666667, 0.3333333333]],
         [[0.0, 0.3333333333, 0.6666666667], [0.3333333333, 0.6666666667]]]

    Here is a slightly larger game ::

        sage: A = matrix([[160, 205, 44],
        ....:       [175, 180, 45],
        ....:       [201, 204, 50],
        ....:       [120, 207, 49]])
        sage: B = matrix([[2, 2, 2],
        ....:             [1, 0, 0],
        ....:             [3, 4, 1],
        ....:             [4, 1, 2]])
        sage: g=NormalFormGame(A,B)
        sage: g.obtain_Nash()
        [[[0.0, 0.0, 0.75, 0.25], [0.0357142857, 0.9642857143, 0.0]]]

    One can also input a single matrix and then a zero sum game is constructed.
    Here is an instance of Rock-Paper-Scissors-Lizard-Spock ::

        sage: A = matrix([[0, -1, 1, 1, -1],
        ....:             [1, 0, -1, -1, 1],
        ....:             [-1, 1, 0, 1 , -1],
        ....:             [-1, 1, -1, 0, 1],
        ....:             [1, -1, 1, -1, 0]])
        sage: g = NormalFormGame(A)
        sage: g.obtain_Nash()
        [[[0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2]]]

    Here is a slightly longer game.
    Consider the following:

    An airline loses two suitcases belonging to two different travelers.
    Both suitcases happen to be identical and contain identical antiques.
    An airline manager tasked to settle the claims of both travelers explains
    that the airline is liable for a maximum of 10 per suitcase, and in
    order to determine an honest appraised value of the antiques the manager
    separates both travelers so they can’t confer, and asks them to write
    down the amount of their value at no less than 2 and no larger than 100.
    He also tells them that if both write down the same number, he will
    treat that number as the true dollar value of both suitcases and
    reimburse both travelers that amount.
    However, if one writes down a smaller number than the other,
    this smaller number will be taken as the true dollar value,
    and both travelers will receive that amount along with a bonus/malus:
    2 extra will be paid to the traveler who wrote down the lower value
    and a 2 deduction will be taken from the person who wrote down the higher amount.
    The challenge is: what strategy should both travelers follow to decide
    the value they should write down?

    In the following we create the game and solve it ::

    sage: K = 10  # Modifying this value lets us play with games of any size
    sage: A = matrix([[min(i,j) + 2 * sign(j-i)  for j in range(2, K+1)]  for i in range(2, K+1)])
    sage: B = matrix([[min(i,j) + 2 * sign(i-j)  for j in range(2, K+1)]  for i in range(2, K+1)])
    sage: g = NormalFormGame(A, B)
    sage: g.obtain_Nash()
    [[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]

    The equilibrium strategy is thus for both players to state that the value of their suitcase is 2.

    """

    def __new__(NormalFormGame, generator):
        r"""
        Creates an Instance of NormalFormGame.

        EXAMPLES:

        A simple 2x2 two player game. ::

            sage: g = NormalFormGame()

        TESTS:

        Raise an error if both matrix and game provided. ::

            sage: g = NormalFormGame(matrix1=4, game=5)
            Traceback (most recent call last):
            ...
            ValueError: Can't input both a matrix and a game

        """
        if type(generator) is list:
            g = Game.new_table([generator[0].dimensions()[0], generator[0].dimensions()[1]])
        if type(generator) is gambit.Game:
            g = generator
        else:
            g = Game.new_table([])

        g.__class__ = NormalFormGame
        return g

    def __init__(self, generator):
        r"""
        Initializes a Normal Form game and checks the inputs.
        """
        if type(generator) is list:
            if len(generator) == 1:
                generator.append(- generator[-1])
            self.matrix_to_game()
        if type(generator) is gambit.Game:
            self.payoff_matrices = self.game_to_matrix(generator)

    def game_to_matrix(game):
        r"""
        Sets ``self.payoff_matrices`` to be the payoff matrices
        associated with current game. This gets called at ``__init__``, but if
        ``self`` changes at any point, ``game_to_matrix`` can be called to
        update those changes.

        EXAMPLES:

        A simple 2 player game. ::

            sage: two_player = Game.new_table([2, 2])
            sage: simple = NormalFormGame(game=two_player)
            sage: simple[int(0), int(0)][int(0)] = int(8)
            sage: simple[int(0), int(0)][int(1)] = int(8)
            sage: simple[int(0), int(1)][int(0)] = int(2)
            sage: simple[int(0), int(1)][int(1)] = int(10)
            sage: simple[int(1), int(0)][int(0)] = int(10)
            sage: simple[int(1), int(0)][int(1)] = int(2)
            sage: simple[int(1), int(1)][int(0)] = int(5)
            sage: simple[int(1), int(1)][int(1)] = int(5)
            sage: simple.game_to_matrix()
            sage: simple.matrix1
            [ 8  2]
            [10  5]
            sage: simple.matrix2
            [ 8 10]
            [ 2  5]

        TESTS:

        Raise an error if a game with more than two players is used. ::

            sage: three_player = Game.new_table([1, 2, 3])
            sage: large_game = NormalFormGame(game=three_player)
            Traceback (most recent call last):
            ...
            ValueError: Only available for games with 2 players
        """
        payoff_matrices = [matrix(len(game.players[0].strategies), len(game.players[1].strategies)) for player in range(len(self.players))]

        if len(self.players) != 2:
            raise ValueError("Only available for games with 2 players")
        for k in list(game.contingencies):
            for player in range(len(self.players)):
                payoff_matrices[player][tuple(k)] = int(self[k][player])
        return payoff_matrices

    def matrix_to_game(self):
        r"""
        Builds a game based on ``self.matrix1`` and ``self.matrix2``. This
        gets called at ``__init__``, but if either matrix gets altered it can
        be called again to update ``self``.

        EXAMPLES:

        A zero-sum game. ::

            sage: single = matrix([[1, 0], [-2, 3]])
            sage: zero_game = NormalFormGame(matrix1=single)
            sage: zero_game.matrix1
            [ 1  0]
            [-2  3]
            sage: zero_game.matrix2
            [-1  0]
            [ 2 -3]

        TESTS:

        Raise error if matrices are not the same size. ::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4], [6, 6]])
            sage: error = NormalFormGame(matrix1=p1, matrix2=p2)
            Traceback (most recent call last):
            ...
            ValueError: Matrices must be the same size
        """
        if self.matrix1.dimensions() != self.matrix2.dimensions():
            raise ValueError("Matrices must be the same size")
        strategysizes = [range(self.payoff_matrices.dimensions()[0]), range(self.payoff_matrices.dimensions()[1])]
        for k in product(*strategysizes):
            for player in range(len(self.players)):
                self[k][player] = int(self.payoff_matrices[player][k])

    def obtain_Nash(self, algorithm="LCP", maximization=True):
        r"""
        A function to return the Nash equilibrium for a game.
        Optional arguments can be used to specify the algorithm used.
        If no algorithm is passed then an attempt is made to use the most
        appropriate algorithm.

        INPUT:

        - ``algorithm`` - the following algorithms should be available through
                          this function:
                * ``"lrs"`` - This algorithm is only suited for 2 player games.
                  See the [insert website here] web site.
                * ``"LCP"`` - This algorithm is only suited for 2 player games.
                  See the [insert website here] web site. NOTE THAT WE NEED TO
                  GET THE ACTUAL NAME OF THE GAMBIT ALGORITHM
                * ``"support enumeration"`` - This is a very inefficient
                  algorithm (in essence a brute force approach).

        - ``maximization``

           - When set to ``True`` (default) it is assumed that players aim to
             maximise their utility.
           - When set to ``False`` (default) it is assumed that players aim to
             minimise their utility.
        """
        if len(self.players) > 2:
            raise NotImplementedError("Nash equilibrium for games with "
                  "more than 2 players have not been implemented yet."
                  "Please see the gambit website [LINK] that has a variety"
                  " of available algorithms")

        if algorithm == "LCP":
            return [self._gambit_profile_to_list(profile) for profile in ExternalLCPSolver().solve(self)]

        if algorithm == "lrs":
            if not is_package_installed('lrs'):
                raise NotImplementedError("lrs is not installed")
            return self._solve_lrs()

        if algorithm == "support enumeration":
            raise NotImplementedError("Support enumeration is not implemented "
                                      "yet")
            return self._solve_enumeration()

    def _gambit_profile_to_list(self, gambitstrategy):
        gambitstrategy = eval(str(gambitstrategy)[str(gambitstrategy).index("["): str(gambitstrategy).index("]") + 1])
        profile = [gambitstrategy[:len(self.players[int(0)].strategies)]]
        for player in list(self.players)[1:]:
            previousplayerstrategylength = len(profile[-1])
            profile.append(gambitstrategy[previousplayerstrategylength: previousplayerstrategylength + len(player.strategies)])
        return profile


    def _solve_lrs(self):
        # file1 = id_generator(6)
        # game = open(file1, "w")
        # game.write("%s %s\n" % (len(self.players[0].strategies), len(self.players[1].strategies))

        # for i in list(self.matrix1):
        #     game.write("\n")
        #     game.write(" ".join([str(e) for e in i]))

        # game.write("\n")

        # for j in list(self.matrix2):
        #     game.write("\n")
        #     game.write(" ".join([str(e) for e in j]))
        # game.close()

        # # Write H representations for each player (really should automate this)
        # file2 = id_generator(6)
        # file3 = id_generator(6)
        # call(["setupnash", file1, file2, file3], stdout=PIPE)
        # # Solve game using lrs:
        # process = Popen(["nash", file2, file3], stdout=PIPE)
        # # Save output
        # lrs_output = [row for row in process.stdout]
        # # Delete lrs files, need to do this without writing hard files
        # for f in [file1, file2, file3]:
        #     remove(f)
        # return lrs_output
        pass

    def _solve_enumeration(self):
        pass


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    """
    A function to generate random file names for physical files needed to communicate with lrs.
    """
    return ''.join(random.choice(chars) for x in range(size))
