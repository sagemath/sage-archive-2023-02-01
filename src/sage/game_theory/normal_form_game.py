from collections import MutableMapping
from itertools import product
from sage.structure.sage_object import SageObject
from sage.misc.package import is_package_installed
from sage.matrix.constructor import matrix
from parser import Parser
from sage.rings.rational import Rational


class NormalFormGame(SageObject, MutableMapping):
    r"""
    An object representing a Normal Form Game. Primarily used to compute the
    Nash Equilibrium.

    INPUT:


    EXAMPLES:

    A basic 2-player game constructed from matrices. ::

        sage: A = matrix([[1, 2], [3, 4]])
        sage: B = matrix([[3, 3], [1, 4]])
        sage: C = NormalFormGame([A, B])
        sage: C.strategy_profiles
        {(0, 1): [2, 3], (1, 0): [3, 1], (0, 0): [1, 3], (1, 1): [4, 4]}

    The same game can be constructed manually. ::

        sage: f = NormalFormGame()
        sage: f.add_player(2)
        sage: f.add_player(2)
        sage: f[0,0][0] = 1
        sage: f[0,0][1] = 3
        sage: f[0,1][0] = 2
        sage: f[0,1][1] = 3
        sage: f[1,0][0] = 3
        sage: f[1,0][1] = 1
        sage: f[1,1][0] = 4
        sage: f[1,1][1] = 4
        sage: f.strategy_profiles
        {(0, 1): [2, 3], (1, 0): [3, 1], (0, 0): [1, 3], (1, 1): [4, 4]}

    We can add an extra strategy to the first player. ::

        sage: f.add_strategy(0)
        sage: f.strategy_profiles
        {(0, 1): [2, 3], (0, 0): [1, 3], (2, 1): [False, False], (2, 0): [False, False], (1, 0): [3, 1], (1, 1): [4, 4]}

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
        ....:             [175, 180, 45],
        ....:             [201, 204, 50],
        ....:             [120, 207, 49]])
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
        sage: g = NormalFormGame(A, 'zero-sum')
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
    def __delitem__(self, key):
        self.strategy_profiles.pop(key, None)

    def __getitem__(self, key):
        return self.strategy_profiles[key]

    def __iter__(self):
        return iter(self.strategy_profiles)

    def __len__(self):
        return len(self.strategy_profiles)

    def __setitem__(self, key, value):
        self.strategy_profiles[key] = value

    def __init__(self, arg1=None, arg2=None):
        r"""
        Initializes a Normal Form game and checks the inputs.

        EXAMPLES:

        Can have games with more than 2 players. ::

            sage: threegame = NormalFormGame()
            sage: threegame.add_player(2)
            sage: threegame.add_player(2)
            sage: threegame.add_player(2)
            sage: threegame[0, 0, 0][0] = 3
            sage: threegame[0, 0, 0][1] = 1
            sage: threegame[0, 0, 0][2] = 4
            sage: threegame[0, 0, 1][0] = 1
            sage: threegame[0, 0, 1][1] = 5
            sage: threegame[0, 0, 1][2] = 9
            sage: threegame[0, 1, 0][0] = 2
            sage: threegame[0, 1, 0][1] = 6
            sage: threegame[0, 1, 0][2] = 5
            sage: threegame[0, 1, 1][0] = 3
            sage: threegame[0, 1, 1][1] = 5
            sage: threegame[0, 1, 1][2] = 8
            sage: threegame[1, 0, 0][0] = 9
            sage: threegame[1, 0, 0][1] = 7
            sage: threegame[1, 0, 0][2] = 9
            sage: threegame[1, 0, 1][0] = 3
            sage: threegame[1, 0, 1][1] = 2
            sage: threegame[1, 0, 1][2] = 3
            sage: threegame[1, 1, 0][0] = 8
            sage: threegame[1, 1, 0][1] = 4
            sage: threegame[1, 1, 0][2] = 6
            sage: threegame[1, 1, 1][0] = 2
            sage: threegame[1, 1, 1][1] = 6
            sage: threegame[1, 1, 1][2] = 4
            sage: threegame.obtain_Nash()
            Traceback (most recent call last):
            ...
            NotImplementedError: Nash equilibrium for games with more than 2 players have not been implemented yet. Please see the gambit website [LINK] that has a variety of available algorithms

        TESTS:

        Raise error if matrices aren't the same size. ::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4], [6, 6]])
            sage: error = NormalFormGame([p1, p2])
            Traceback (most recent call last):
            ...
            ValueError: matrices must be the same size

        """
        self.players = []
        flag = 'n-player'
        if type(arg1) is list:
            matrices = arg1
            flag = 'two-matrix'
            if matrices[0].dimensions() != matrices[1].dimensions():
                raise ValueError("matrices must be the same size")
        elif arg2 == 'zero-sum':
            flag = 'two-matrix'
            matrices = [arg1]
            matrices.append(-arg1)

        if flag == 'two-matrix':
            self._two_matrix_game(matrices)

    def _two_matrix_game(self, matrices):
        r"""
        Populates ``self.strategy_profiles`` with the values from 2 matrices.

        EXAMPLES:

        A small example game. ::

            sage: s1 = matrix([[1, 0], [-2, 3]])
            sage: s2 = matrix([[3, 2], [-1, 0]])
            sage: two_game = NormalFormGame()
            sage: two_game._two_matrix_game([s1, s2])
        """
        self.players = []
        self.strategy_profiles = {}
        self.add_player(matrices[0].dimensions()[0])
        self.add_player(matrices[1].dimensions()[1])
        for key in self.strategy_profiles:
            self.strategy_profiles[key] = [matrices[0][key], matrices[1][key]]

    def add_player(self, num_strategies):
        r"""
        Adds a player to a NormalFormGame.

        INPUT:

        - ``num_strategies`` - the number of strategies the player should have.

        EXAMPLES:

        sage: g = NormalFormGame()
        sage: g.add_player(3)
        sage: g.add_player(4)
        sage: g.strategy_profiles
        {(0, 1): [False, False], (1, 2): [False, False], (0, 0): [False, False], (2, 1): [False, False], (0, 2): [False, False], (2, 0): [False, False], (1, 3): [False, False], (2, 3): [False, False], (2, 2): [False, False], (1, 0): [False, False], (0, 3): [False, False], (1, 1): [False, False]}
        """
        self.players.append(_Player(num_strategies))
        self._generate_strategy_profiles(True)

    def _generate_strategy_profiles(self, replacement):
        r"""
        Creates all the required keys for ``self.strategy_profiles``.
        """
        strategy_sizes = [range(p.num_strategies) for p in self.players]
        if replacement is True:
            self.strategy_profiles = {}
        for profile in product(*strategy_sizes):
            if profile not in self.strategy_profiles.keys():
                self.strategy_profiles[profile] = [False]*len(self.players)

    def add_strategy(self, player):
        r"""
        Adds a strategy to a player.

        INPUT:

        - ``player`` - the index of the player.

        EXAMPLES:


        """
        self.players[player].add_strategy()
        self._generate_strategy_profiles(False)

    def _is_complete(self):
        r"""
        Checks if ``strategy_profiles`` has been completed and returns a
        boolean.
        """
        for profile in self.strategy_profiles.values():
            for i in profile:
                if i == 0:
                    pass
                elif i is False:
                    return False
        return True

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
           - When set to ``False`` it is assumed that players aim to
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

        if not self._is_complete():
            raise ValueError("strategy_profiles hasn't been populated")

        if algorithm == "LCP":
            return self._solve_LCP(maximization)

        if algorithm == "lrs":
            if not is_package_installed('lrs'):
                raise NotImplementedError("lrs is not installed")
            m1, m2 = self._game_two_matrix()
            if maximization is False:
                min1 = - m1
                min2 = - m2
                nasheq = self._solve_lrs(min1, min2)
            else:
                nasheq = self._solve_lrs(m1, m2)
            return nasheq

    def _game_two_matrix(self):
        m1 = matrix(self.players[0].num_strategies, self.players[1].num_strategies)
        m2 = matrix(self.players[0].num_strategies, self.players[1].num_strategies)
        for key in self.strategy_profiles:
                m1[key] = self[key][0]
                m2[key] = self[key][1]
        return m1, m2

    def _solve_LCP(self, maximization):
        if not is_package_installed('gambit'):
                raise NotImplementedError("gambit is not installed")

        from gambit import Game
        from gambit.nash import ExternalLCPSolver

        strategy_sizes = [p.num_strategies for p in self.players]
        g = Game.new_table(strategy_sizes)

        scalar = 1
        for key in self.strategy_profiles:
            for player in range(len(self.players)):
                utility = Rational(self.strategy_profiles[key][player])
                scalar *= utility.denom()
                if maximization is False:
                    g[key][player] = - int(self.strategy_profiles[key][player])
                else:
                    g[key][player] = int(self.strategy_profiles[key][player])

        if scalar != 1:
            for key in self.strategy_profiles:
                for player in range(len(self.players)):
                    g[key][player] *= scalar

        output = ExternalLCPSolver().solve(g)
        nasheq = Parser(output, g).format_gambit()
        return nasheq

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
        nasheq = Parser(lrs_output).format_lrs()
        return nasheq

    def _solve_enumeration(self):
        # call _is_complete()
        # write algorithm at some point
        pass

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
        m = self.players[0].num_strategies
        n = self.players[1].num_strategies
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


class _Player():
    def __init__(self, num_strategies):
        self.num_strategies = num_strategies

    def add_strategy(self):
        self.num_strategies += 1