from itertools import product
from sage.structure.sage_object import SageObject
from sage.misc.package import is_package_installed
from formatter import Formatter
from player import _Player


class NormalFormGame(SageObject):
    r"""
    EXAMPLES:

    A basic 2-player game constructed from two matrices. ::

        sage: A = matrix([[1, 2], [3, 4]])
        sage: B = matrix([[3, 3], [1, 4]])
        sage: C = NormalFormGame([A, B])

    """
    def __init__(self, arg1=None, arg2=None):

        self.players = []
        if type(arg1) is list:
            matrices = arg1
            flag = 'two-matrix'
            if matrices[0].dimensions() != matrices[1].dimensions():
                raise ValueError("matrices must be the same size")
        # elif arg1 == 'bi-matrix':
        #     flag = arg1
        #     bimatrix = arg2
        elif arg1 == 'zero-sum':
            flag = 'two-matrix'
            matrices = [arg2]
            matrices.append(-arg2)

        if flag == 'two-matrix':
            self._two_matrix_game(matrices)
        # elif flag == 'bi-matrix':
        #     self._bimatrix_game(bimatrix)

    def _two_matrix_game(self, matrices):
        self.add_player(matrices[0].dimensions()[0])
        self.add_player(matrices[1].dimensions()[1])
        for key in self.strategy_profiles:
            self.strategy_profiles[key] = (matrices[0][key], matrices[1][key])

    # def _bimatrix_game(self, bimatrix):
    #     self.add_player(bimatrix.dimensions()[0])
    #     self.add_player(bimatrix.dimensions()[1])
    #     for key in self.strategy_profiles:
    #         self.strategy_profiles[key] = bimatrix[key]

    def add_player(self, num_strategies):
        self.players.append(_Player(num_strategies))
        self.generate_strategy_profiles()

    def generate_strategy_profiles(self):
        self.strategy_profiles = {}
        strategy_sizes = [range(p.num_strategies) for p in self.players]
        for profile in product(*strategy_sizes):
            self.strategy_profiles[profile] = False

    def add_strategy(self, player):
        self.players[player].add_strategy()
        self.generate_strategy_profiles()

    def define_utility_vector(self, strategy_profile, utility_vector):
        self.strategy_profiles[strategy_profile] = utility_vector

    def _is_complete(self):
        return all(self.strategy_profiles.values())

    def obtain_Nash(self, algorithm="LCP", maximization=True):
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

    def _solve_LCP(self, maximization):
        if not is_package_installed('gambit'):
                raise NotImplementedError("gambit is not installed")
        from gambit import Game
        from gambit.nash import ExternalLCPSolver
        strategy_sizes = [p.num_strategies for p in self.players]
        g = Game.new_table(strategy_sizes)
        for key in self.strategy_profiles:
            for player in range(len(self.players)):
                if maximization is False:
                    g[key][player] = - self.strategy_profiles[key][player]
                else:
                    g[key][player] = self.strategy_profiles[key][player]
        output = ExternalLCPSolver().solve(self)
        nasheq = Formatter(output, self).format_gambit()
        return nasheq

    def _solve_lrs(self):
        # call _is_complete()
        # create polytope??
        pass

    def _solve_enumeration(self):
        # call _is_complete()
        # write algorithm at some point
        pass
