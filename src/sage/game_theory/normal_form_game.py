from itertools import product
from sage.structure.sage_object import SageObject


class NormalFormGame(SageObject):

    def __init__(self):
        self.players = []
        self.strategy_profiles = {}

    def add_player(self, num_strategies):
        self.players.append(_Player(num_strategies))
        self.generate_strategy_profiles()

    def generate_strategy_profiles(self):
        strategy_sizes = [range(p.num_strategies) for p in self.players]
        for profile in product(strategy_sizes):
            self.strategy_profiles[profile] = False

    def add_strategy(self, player):
        self.players[player].add_strategy()
        self.generate_strategy_profiles()

    def define_utility_vector(self, strategy_profile, utility_vector):
        self.strategy_profiles[strategy_profile] = utility_vector

    def _is_complete(self):
        # will check that all strategy profiles have a legitimate payoff_vector
        # solve_nash etc can't be called until that's true.
        return all(self.strategy_profiles.values())

    def obtain_Nash(self, algorithm="LCP", maximization=True):
        # copy from old
        pass

    def _solve_gambit(self):
        # call _is_complete()
        # create a gambit.Game object using self.strategy_profiles
        # should be very quick
        pass

    def _solve_lrs(self):
        # call _is_complete()
        # create polytope??
        pass

    def _solve_enumeration(self):
        # call _is_complete()
        # write algorithm at some point
        pass


class _Player():
    def __init__(self, num_strategies):
        self.num_strategies = num_strategies

    def add_strategy(self):
        self.num_strategies += 1
