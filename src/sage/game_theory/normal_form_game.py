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
from sage.matrix.constructor import matrix
from formatter import Formatter


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
        sage: g = NormalFormGame(gambitgame)
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
        sage: g = NormalFormGame([A, B])

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
        sage: g=NormalFormGame([A, B])
        sage: g.obtain_Nash()
        [[[0.0, 0.0, 0.75, 0.25], [0.0357142857, 0.9642857143, 0.0]]]

    One can also input a single matrix and then a zero sum game is constructed.
    Here is an instance of Rock-Paper-Scissors-Lizard-Spock ::

        sage: A = matrix([[0, -1, 1, 1, -1],
        ....:             [1, 0, -1, -1, 1],
        ....:             [-1, 1, 0, 1 , -1],
        ....:             [-1, 1, -1, 0, 1],
        ....:             [1, -1, 1, -1, 0]])
        sage: g = NormalFormGame([A])
        sage: g.obtain_Nash()
        [[[0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2]]]

    Here is a slightly longer game.
    Consider the following:

    An airline loses two suitcases belonging to two different travelers. Both
    suitcases happen to be identical and contain identical antiques. An
    airline manager tasked to settle the claims of both travelers explains
    that the airline is liable for a maximum of 10 per suitcase, and in order
    to determine an honest appraised value of the antiques the manager
    separates both travelers so they can’t confer, and asks them to write down
    the amount of their value at no less than 2 and no larger than 100. He
    also tells them that if both write down the same number, he will treat
    that number as the true dollar value of both suitcases and reimburse both
    travelers that amount.
    However, if one writes down a smaller number than the other, this smaller
    number will be taken as the true dollar value, and both travelers will
    receive that amount along with a bonus/malus: 2 extra will be paid to the
    traveler who wrote down the lower value and a 2 deduction will be taken
    from the person who wrote down the higher amount. The challenge is: what
    strategy should both travelers follow to decide the value they should
    write down?

    In the following we create the game and solve it ::

    sage: K = 10  # Modifying this value lets us play with games of any size
    sage: A = matrix([[min(i,j) + 2 * sign(j-i)  for j in range(2, K+1)]  for i in range(2, K+1)])
    sage: B = matrix([[min(i,j) + 2 * sign(i-j)  for j in range(2, K+1)]  for i in range(2, K+1)])
    sage: g = NormalFormGame([A, B])
    sage: g.obtain_Nash()
    [[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]

    The equilibrium strategy is thus for both players to state that the value
    of their suitcase is 2.

    """

    def __new__(NormalFormGame, generator):
        r"""
        Creates an Instance of NormalFormGame.
        """
        if type(generator) is list:
            g = Game.new_table([generator[0].dimensions()[0], generator[0].dimensions()[1]])
        elif type(generator) is Game:
            g = generator
        else:
            g = Game.new_table([])

        g.__class__ = NormalFormGame
        return g

    def __init__(self, generator):
        r"""
        Initializes a Normal Form game and checks the inputs.

        EXAMPLES:

        Check for type of inputs. ::

            sage: NormalFormGame(4)
            Traceback (most recent call last):
            ...
            TypeError: Generator function must be a list or Game

        Can have games with more than 2 players. ::

            sage: threeplayer = Game.new_table([2,2,2])
            sage: threegame = NormalFormGame(threeplayer)
            sage: threegame[int(0), int(0), int(0)][int(0)] = int(3)
            sage: threegame[int(0), int(0), int(0)][int(1)] = int(1)
            sage: threegame[int(0), int(0), int(0)][int(2)] = int(4)
            sage: threegame[int(0), int(0), int(1)][int(0)] = int(1)
            sage: threegame[int(0), int(0), int(1)][int(1)] = int(5)
            sage: threegame[int(0), int(0), int(1)][int(2)] = int(9)
            sage: threegame[int(0), int(1), int(0)][int(0)] = int(2)
            sage: threegame[int(0), int(1), int(0)][int(1)] = int(6)
            sage: threegame[int(0), int(1), int(0)][int(2)] = int(5)
            sage: threegame[int(0), int(1), int(1)][int(0)] = int(3)
            sage: threegame[int(0), int(1), int(1)][int(1)] = int(5)
            sage: threegame[int(0), int(1), int(1)][int(2)] = int(8)
            sage: threegame[int(1), int(0), int(0)][int(0)] = int(9)
            sage: threegame[int(1), int(0), int(0)][int(1)] = int(7)
            sage: threegame[int(1), int(0), int(0)][int(2)] = int(9)
            sage: threegame[int(1), int(0), int(1)][int(0)] = int(3)
            sage: threegame[int(1), int(0), int(1)][int(1)] = int(2)
            sage: threegame[int(1), int(0), int(1)][int(2)] = int(3)
            sage: threegame[int(1), int(1), int(0)][int(0)] = int(8)
            sage: threegame[int(1), int(1), int(0)][int(1)] = int(4)
            sage: threegame[int(1), int(1), int(0)][int(2)] = int(6)
            sage: threegame[int(1), int(1), int(1)][int(0)] = int(2)
            sage: threegame[int(1), int(1), int(1)][int(1)] = int(6)
            sage: threegame[int(1), int(1), int(1)][int(2)] = int(4)
            sage: threegame.obtain_Nash()
            Traceback (most recent call last):
            ...
            NotImplementedError: Nash equilibrium for games with more than 2 players have not been implemented yet. Please see the gambit website [LINK] that has a variety of available algorithms
        """
        if type(generator) is not list and type(generator) is not NormalFormGame:
            raise TypeError("Generator function must be a list or Game")

        if type(generator) is list:
            if len(generator) == 1:
                generator.append(-generator[-1])
            self.payoff_matrices = generator
            self.matrix_to_game()
        if type(generator) is Game and len(self.players) == 2:
            self.payoff_matrices = self.game_to_matrix()

    def game_to_matrix(self):
        r"""
        Sets ``self.payoff_matrices`` to be the payoff matrices associated
        with current game. This gets called at ``__init__``, but if ``self``
        changes at any point, ``game_to_matrix`` can be called to update those
        changes.

        EXAMPLES:

        A simple 2 player game. ::

            sage: two_player = Game.new_table([2, 2])
            sage: simple = NormalFormGame(two_player)
            sage: simple[int(0), int(0)][int(0)] = int(8)
            sage: simple[int(0), int(0)][int(1)] = int(8)
            sage: simple[int(0), int(1)][int(0)] = int(2)
            sage: simple[int(0), int(1)][int(1)] = int(10)
            sage: simple[int(1), int(0)][int(0)] = int(10)
            sage: simple[int(1), int(0)][int(1)] = int(2)
            sage: simple[int(1), int(1)][int(0)] = int(5)
            sage: simple[int(1), int(1)][int(1)] = int(5)
            sage: simple.payoff_matrices = simple.game_to_matrix()
            sage: simple.payoff_matrices
            [
            [ 8  2]  [ 8 10]
            [10  5], [ 2 5]
            ]

        TESTS:

        Raise an error if a game with more than two players is used. ::

            sage: three_player = Game.new_table([1, 2, 3])
            sage: large_game = NormalFormGame(three_player)
            sage: large_game.game_to_matrix()
            Traceback (most recent call last):
            ...
            ValueError: Only available for games with 2 players
        """

        if len(self.players) != 2:
            raise ValueError("Only available for games with 2 players")

        payoff_matrices = [matrix(len(self.players[0].strategies), len(self.players[1].strategies)) for player in range(len(self.players))]

        for k in list(self.contingencies):
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
            sage: zero_game = NormalFormGame([single])
            sage: zero_game.payoff_matrices
            [
            [ 1  0]  [-1  0]
            [-2  3], [ 2 -3]
            ]

        TESTS:

        Raise error if matrices are not the same size. ::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4], [6, 6]])
            sage: error = NormalFormGame([p1, p2])
            Traceback (most recent call last):
            ...
            ValueError: Matrices must be the same size
        """
        if len(set([matrix.dimensions() for matrix in self.payoff_matrices])) != 1:
            raise ValueError("Matrices must be the same size")
        strategysizes = [range(self.payoff_matrices[0].dimensions()[0]), range(self.payoff_matrices[0].dimensions()[1])]
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

        EXAMPLES:

        A game with 2 equilibria when ``maximization`` is ``True`` and 3 when
        ``maximization`` is ``False``. ::

            sage: A = matrix([[160, 205, 44],
            ....:       [175, 180, 45],
            ....:       [201, 204, 50],
            ....:       [120, 207, 49]])
            sage: B = matrix([[2, 2, 2],
            ....:             [1, 0, 0],
            ....:             [3, 4, 1],
            ....:             [4, 1, 2]])
            sage: g=NormalFormGame([A, B])
            sage: g.obtain_Nash(algorithm='lrs')
            [([0, 0, 3/4, 1/4], [1/28, 27/28, 0])]
            sage: g.obtain_Nash(algorithm='lrs', maximization=False)
            [([1, 0, 0, 0], [127/1212, 115/1212, 485/606]), ([0, 1, 0, 0], [0, 1/26, 25/26])]

        This particular game has 3 Nash equilibria. ::

            sage: A = matrix([[3,3],
            ....:             [2,5],
            ....:             [0,6]])
            sage: B = matrix([[3,2],
            ....:             [2,6],
            ....:             [3,1]])
            sage: g = NormalFormGame([A, B])
            sage: g.obtain_Nash(maximization=False)
            [[[1.0, 0.0, 0.0], [0.0, 1.0]]]

        2 random matrices. ::

            sage: player1 = matrix([[2, 8, -1, 1, 0],
            ....:                   [1, 1, 2, 1, 80],
            ....:                   [0, 2, 15, 0, -12],
            ....:                   [-2, -2, 1, -20, -1],
            ....:                   [1, -2, -1, -2, 1]])
            sage: player2 = matrix([[0, 8, 4, 2, -1],
            ....:                   [6, 14, -5, 1, 0],
            ....:                   [0, -2, -1, 8, -1],
            ....:                   [1, -1, 3, -3, 2],
            ....:                   [8, -4, 1, 1, -17]])
            sage: fivegame = NormalFormGame([player1, player2])
            sage: fivegame.obtain_Nash()
            [[[1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0, 0.0]]]
            sage: fivegame.obtain_Nash(algorithm='lrs')
            [([1, 0, 0, 0, 0], [0, 1, 0, 0, 0])]

        """
        if len(self.players) > 2:
            raise NotImplementedError("Nash equilibrium for games with more "
                                      "than 2 players have not been "
                                      "implemented yet. Please see the gambit "
                                      "website [LINK] that has a variety of "
                                      "available algorithms")

        if algorithm == "LCP":
            if maximization is False:
                for k in list(self.contingencies):
                    for player in range(len(self.players)):
                        self[k][player] *= -1

            output = ExternalLCPSolver().solve(self)
            nasheq = Formatter(output, self).format_gambit()

            if maximization is False:
                for k in list(self.contingencies):
                    for player in range(len(self.players)):
                        self[k][player] *= -1
            return nasheq

        if algorithm == "lrs":
            if not is_package_installed('lrs'):
                raise NotImplementedError("lrs is not installed")
            if maximization is False:
                min1 = - self.payoff_matrices[0]
                min2 = - self.payoff_matrices[1]
                nasheq = self._solve_lrs(min1, min2)
            else:
                nasheq = self._solve_lrs(self.payoff_matrices[0],
                                         self.payoff_matrices[1])
            return nasheq

        if algorithm == "support enumeration":
            raise NotImplementedError("Support enumeration is not implemented "
                                      "yet")
            return self._solve_enumeration()

    def _solve_lrs(self, m1, m2):
        r"""
        EXAMPLES:

        A simple game. ::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: C._solve_lrs(A, B)
            [([0, 1], [0, 1])]

        2 random matrices. ::

        sage: p1 = matrix([[-1, 4, 0, 2, 0],
        ....:              [-17, 246, -5, 1, -2],
        ....:              [0, 1, 1, -4, -4],
        ....:              [1, -3, 9, 6, -1],
        ....:              [2, 53, 0, -5, 0]])
        sage: p2 = matrix([[0, 1, 1, 3, 1],
        ....:              [3, 9, 44, -1, -1],
        ....:              [1, -4, -1, -3, 1],
        ....:              [1, 0, 1, 0, 0,],
        ....:              [1, -3, 1, 21, -2]])
        sage: biggame = NormalFormGame([p1, p2])
        sage: biggame._solve_lrs(p1, p2)
        [([0, 0, 0, 20/21, 1/21], [11/12, 0, 0, 1/12, 0]), ([0, 0, 0, 1, 0], [9/10, 0, 1/10, 0, 0])]
        """
        from sage.misc.temporary_file import tmp_filename
        from subprocess import Popen, PIPE
        # so that we don't call _Hrepresentation() twice.
        in_str = self._Hrepresentation(m1, m2)
        game1_str = in_str[0]
        game2_str = in_str[1]

        g1_name = tmp_filename()
        g2_name = tmp_filename()
        g1_file = file(g1_name, 'w')
        g2_file = file(g2_name, 'w')
        g1_file.write(game1_str)
        g1_file.close()
        g2_file.write(game2_str)
        g2_file.close()

        process = Popen(['nash', g1_name, g2_name], stdout=PIPE)
        lrs_output = [row for row in process.stdout]
        nasheq = Formatter(lrs_output).format_lrs()
        return nasheq

    def _solve_enumeration(self):
        pass

    def _scale_matrices(self):
        r"""
        Returns the ``payoff_matrices`` with integer coeficcients so that they
        can be solved using gambit.

        EXAMPLES:

        Simple example. ::

            sage: a = matrix([[1, 0], [1/3, 4]])
            sage: b = matrix([[2.5, 3], [-0.75, 4]])
            sage: c = NormalFormGame([a, b])
            sage: c._scale_matrices()
            sage: c.payoff_matrices[0]
            [24  0]
            [ 8 96]
            sage: c.payoff_matrices[1]
            [ 60.0000000000000  72.0000000000000]
            [-18.0000000000000  96.0000000000000]
        """
        from sage.rings.rational import Rational
        scalar = 1
        for player in range(len(self.players)):
            for k in list(self.contingencies):
                utility = Rational(self.payoff_matrices[player][tuple(k)])
                scalar *= utility.denom()

        self.payoff_matrices[0] *= scalar
        self.payoff_matrices[1] *= scalar
        self.matrix_to_game()

    def _Hrepresentation(self, m1, m2):
        r"""
        Creates the H-representation strings required to use lrs nash.

        EXAMPLES::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: print C._Hrepresentation(A, B)[0]
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 1 0 0
            0 0 1 0
            0 -3 -1 1
            0 -3 -4 1
            -1 1 1 0
            end
            <BLANKLINE>
            sage: print C._Hrepresentation(A, B)[1]
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 -1 -2 1
            0 -3 -4 1
            0 1 0 0
            0 0 1 0
            -1 1 1 0
            end
            <BLANKLINE>

        """
        from sage.geometry.polyhedron.misc import _to_space_separated_string
        m = len(self.players[0].strategies)
        n = len(self.players[1].strategies)
        midentity = list(matrix.identity(m))
        nidentity = list(matrix.identity(n))

        s = 'H-representation\n'
        s += 'linearity 1 ' + str(m + n + 1) + '\n'
        s += 'begin\n'
        s += str(m + n + 1) + ' ' + str(m + 2) + ' rational\n'
        for f in list(midentity):
            s += '0 ' + _to_space_separated_string(f) + ' 0 \n'
        for e in list(m2.transpose()):
            s += '0 ' + _to_space_separated_string(-e) + '  1 \n'
        s += '-1 '
        for g in range(m):
            s += '1 '
        s += '0 \n'
        s += 'end\n'

        t = 'H-representation\n'
        t += 'linearity 1 ' + str(m + n + 1) + '\n'
        t += 'begin\n'
        t += str(m + n + 1) + ' ' + str(n + 2) + ' rational\n'
        for e in list(m1):
            t += '0 ' + _to_space_separated_string(-e) + '  1 \n'
        for f in list(nidentity):
            t += '0 ' + _to_space_separated_string(f) + ' 0 \n'
        t += '-1 '
        for g in range(n):
            t += '1 '
        t += '0 \n'
        t += 'end\n'
        return s, t
