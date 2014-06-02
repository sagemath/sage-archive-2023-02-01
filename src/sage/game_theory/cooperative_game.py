from itertools import permutations

class CooperativeGame():
    def __init__(self, characteristic_function, player_list):
        self.char_fun = characteristic_function
        self.number_players = len(player_list)
        self.player_list = player_list

    def shapley_value(self):
        r"""
        Return the payoff vector for co-operative game.
        """
        payoff_vector = []
        for i in self.player_list:
            player_contribution = self.marginal_contributions(i)
            average  = sum(player_contribution) / len(player_contribution)
            payoff_vector.append(average)
        return payoff_vector

    def is_monotone():
        r"""
        Returns True if co-operative game is monotonic.
        """

    def is_superadditive():
        r"""
        Returns True if co-operative game is superadditive.
        """

    def marginal_contributions(self, player):
        contributions = []
        for i in permutations(self.player_list):
            contributions.append(self.marginal_of_pi(player, i))
        return contributions

    def marginal_of_pi(self, player, pi):
        player_and_pred = self.get_predecessors_and_player(player, pi)
        predecessors = [x for x in player_and_pred if x != player]
        if predecessors == None:
            predecessors = ()
        else:
            predecessors = tuple(predecessors)
        player_and_pred = tuple(player_and_pred)
        value = self.char_fun[player_and_pred] - self.char_fun[predecessors]
        return value

    def get_predecessors_and_player(self, player, permutation):
        predecessors = []
        for k in permutation:
            if k == player:
                predecessors.append(k)
                predecessors.sort()
                return predecessors
            else:
                predecessors.append(k)

"""
test_function = {(): 0,
                 ('A',): 6,
                 ('B',): 12,
                 ('C',): 42,
                 ('A', 'B',): 12,
                 ('A', 'C',): 42,
                 ('B', 'C',): 42,
                 ('A', 'B', 'C',): 42}

test_game = CooperativeGame(test_function, ['A', 'B', 'C'])
print test_game.shapley_value()
"""
