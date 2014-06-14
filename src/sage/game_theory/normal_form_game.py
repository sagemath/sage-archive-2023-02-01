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

    If only ``matrix1`` is provided, ``matrix2`` will be created as the
    negative of ``matrix1`` so that a zero-sum game is created. If a ``game``
    is provided, ``matrix1`` and ``matrix2`` will be generated automatically.

    - ``matrix1`` - a matrix representing the payoff for player1 in a 2 player
                    Normal Form game.
    - ``matrix2`` - a matrix representing the payoff for player2 in a 2 player
                    Normal Form game.
    - ``game`` - an instance of gambit.Game().

    EXAMPLES:

    A basic 2-player game constructed from matrices. ::

        sage: a = matrix([[1, 2], [3, 4]])
        sage: b = matrix([[3, 3], [1, 4]])
        sage: c = NormalFormGame(matrix1=a, matrix2=b)

    This can be given a title and the players can be named. ::

        sage: c.title = "Simple Game"
        sage: c.players[int(0)].label = "James"
        sage: c.players[int(1)].label = "Vince"
        sage: c
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

        sage: gam = Game.new_table([2, 2])
        sage: gam[int(0), int(0)][int(0)] = int(8)
        sage: gam[int(0), int(0)][int(1)] = int(8)
        sage: gam[int(0), int(1)][int(0)] = int(2)
        sage: gam[int(0), int(1)][int(1)] = int(10)
        sage: gam[int(1), int(0)][int(0)] = int(10)
        sage: gam[int(1), int(0)][int(1)] = int(2)
        sage: gam[int(1), int(1)][int(0)] = int(5)
        sage: gam[int(1), int(1)][int(1)] = int(5)
        sage: gam.title = "A prisoner's dilemma game"
        sage: gam.players[int(0)].label = "Alphonse"
        sage: gam.players[int(1)].label = "Gaston"
        sage: g = NormalFormGame(game=gam)
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

        sage: gam.obtain_Nash()
        [[[0.0, 1.0], [0.0, 1.0]]]

    Here is an example of a 3 by 2 game ::

        sage: A = matrix([[3,3],
        ....:             [2,5],
        ....:             [0,6]])
        sage: B = matrix([[3,2],
        ....:             [2,6],
        ....:             [3,1]])
        sage: game = NormalFormGame(A,B)


    This particular game has 3 Nash equilibrium::

        sage: game.obtain_Nash()
        [[[1.0, 0.0, 0.0], [1.0, 0.0]],
         [[0.8, 0.2, 0.0], [0.6666666667, 0.3333333333]],
         [[0.0, 0.3333333333, 0.6666666667], [0.3333333333, 0.6666666667]]]
    """

    def __new__(NormalFormGame, matrix1=False, matrix2=False, game=False):
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
        if matrix1 and game:
            raise ValueError("Can't input both a matrix and a game")
        if matrix1:
            g = Game.new_table([matrix1.dimensions()[0], matrix1.dimensions()[1]])
        elif game:
            g = game
        else:
            g = Game.new_table([])

        g.__class__ = NormalFormGame
        return g

    def __init__(self, matrix1=False, matrix2=False, game=False):
        r"""
        Initializes a Normal Form game and checks the inputs.
        """
        self.matrix1 = None
        self.matrix2 = None
        if not matrix1 and not game:
            pass
        elif not matrix1 and game:
            self.game_to_matrix()
        else:
            self.matrix1 = matrix1
            if not matrix2:
                self.matrix2 = - self.matrix1
            else:
                self.matrix2 = matrix2
            self.matrix_to_game()

    def game_to_matrix(self):
        r"""
        Sets ``self.matrix1`` and ``self.matrix2`` to be the payoff matrices
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
        self.matrix1 = matrix(len(self.players[0].strategies), len(self.players[1].strategies))
        self.matrix2 = copy(self.matrix1)
        if len(self.players) != 2:
            raise ValueError("Only available for games with 2 players")
        for k in list(self.contingencies):
            self.matrix1[tuple(k)] = int(self[k][0])
            self.matrix2[tuple(k)] = int(self[k][1])

    def matrix_to_game(self):
        r"""
        Builds a game based on ``self.matrix1`` and ``self.matrix2``. This
        gets called at ``__init__``, but it either matrix gets altered it can
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
        p1_strats = range(self.matrix1.dimensions()[0])
        p2_strats = range(self.matrix1.dimensions()[1])
        for k in product(p1_strats, p2_strats):
                self[k][0] = int(self.matrix1[k])
                self[k][1] = int(self.matrix2[k])

    def obtain_Nash(self, algorithm="LCP"):
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
            return [self._gambit_profile_to_vector(profile) for profile in ExternalLCPSolver().solve(self)]

        if algorithm == "lrs":
            if not is_package_installed('lrs'):
                raise NotImplementedError("lrs is not installed")
            return self._solve_lrs()

        if algorithm == "support enumeration":
            raise NotImplementedError("Support enumeration is not implemented "
                                      "yet")
            return self._solve_enumeration()

    def _gambit_profile_to_vector(self, gambitstrategy):
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
