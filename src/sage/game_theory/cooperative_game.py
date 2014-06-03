r"""
Co-operative games with N players.

This module implements characteristic function cooperative games.
The main contribution is a class for a characteristic function game
that takes a characteristic function (as a dictionary) as an input.
Methods to calculate the Shapley value (a fair way of sharing common
resources: https://www.youtube.com/watch?v=aThG4YAFErw) as well as
test properties of the game (monotonicity, super additivity) are also included.

AUTHOR::
- James Campbell 06-2014: Original version

- Vince Knight 06-2014
"""

#*****************************************************************************
#       Copyright (C) 2014 James Campbell james.campbell@tanti.org.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import permutations
from sage.structure.sage_object import SageObject


class CooperativeGame(SageObject):
    r"""
    An object representing a co-operative game. Primarily used to compute the
    Shapley Value, but can also provide other information.

    INPUT:

    - characteristic_function - a dictionary containing all possible sets of players.
        * Key - Each set must be entered as a tuple, not a string. Apart from
                the empty set, they must all end with a comma.
        * Value - A real number representing each set of players contribution.

    - player_list - a list of all the players.

    - payoff_vector - default = ``False``, a list can be passed instead but
                      this will be overwritten if shapley_value is called.

    EXAMPLES::
    Basic example of how to implement a co-operative game. ::

        sage: integer_function = {(): 0,
        ....:                  (1,): 6,
        ....:                  (2,): 12,
        ....:                  (3,): 42,
        ....:                  (1, 2,): 12,
        ....:                  (1, 3,): 42,
        ....:                  (2, 3,): 42,
        ....:                  (1, 2, 3,): 42}
        sage: integer_game = CooperativeGame(integer_function, [1, 2, 3])

    We can also use strings instead of numbers. ::

        sage: letter_function = {(): 0,
        ....:                  ('A',): 6,
        ....:                  ('B',): 12,
        ....:                  ('C',): 42,
        ....:                  ('A', 'B',): 12,
        ....:                  ('A', 'C',): 42,
        ....:                  ('B', 'C',): 42,
        ....:                  ('A', 'B', 'C',): 42}
        sage: letter_game = CooperativeGame(letter_function, ['A', 'B', 'C'])
    """

    def __init__(self, characteristic_function, player_list, payoff_vector=False):
        self.char_fun = characteristic_function
        self.number_players = len(player_list)
        self.player_list = player_list
        self.payoff_vector = payoff_vector

    def shapley_value(self):
        r"""
        Return the payoff vector for co-operative game.

        EXAMPLES::
        A typical example of the use of shapley_value. ::

            sage: test_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: test_game = CooperativeGame(test_function, [1, 2, 3])
            sage: test_game.payoff_vector
            False
            sage: print test_game.shapley_value()
            [2, 5, 35]
            sage: test_game.payoff_vector
            [2, 5, 35]
        """
        payoff_vector = []
        for i in self.player_list:
            player_contribution = self.marginal_contributions(i)
            average = sum(player_contribution) / len(player_contribution)
            payoff_vector.append(average)
        self.payoff_vector = payoff_vector
        return payoff_vector

    def is_monotone(self):
        r"""
        Returns True if co-operative game is monotonic.
        """
        sets = list(self.char_fun.keys())
        for i, k in permutations(sets, 2):
            if set(i) <= set(k) and self.char_fun[i] > self.char_fun[k]:
                return False
            else:
                pass
        return True

    def is_superadditive(self):
        r"""
        Returns True if co-operative game is superadditive.
        """
        sets = list(self.char_fun.keys())
        for i, k in permutations(sets, 2):
            j = tuple(set(i) | set(k))
            if self.char_fun[j] < self.char_fun[i] + self.char_fun[k]:
                return False
            else:
                pass
        return True

    def marginal_contributions(self, player):
        contributions = []
        for i in permutations(self.player_list):
            contributions.append(self.marginal_of_pi(player, i))
        return contributions

    def marginal_of_pi(self, player, pi):
        predecessors, player_and_pred = self.get_predecessors(player, pi)
        if predecessors is None:
            predecessors = ()
        else:
            predecessors = tuple(predecessors)
        player_and_pred = tuple(player_and_pred)
        value = self.char_fun[player_and_pred] - self.char_fun[predecessors]
        return value

    def get_predecessors(self, player, permutation):
        pred = []
        play_and_pred = []
        for k in permutation:
            if k == player:
                play_and_pred.append(k)
                pred.sort()
                play_and_pred.sort()
                return pred, play_and_pred
            else:
                pred.append(k)
                play_and_pred.append(k)

    def show(self):
        r"""
        EXAMPLES::
        Typical use of the show function.::

            sage: letter_function = {(): 0,
            ....:                  ('A',): 6,
            ....:                  ('B',): 12,
            ....:                  ('C',): 42,
            ....:                  ('A', 'B',): 12,
            ....:                  ('A', 'C',): 42,
            ....:                  ('B', 'C',): 42,
            ....:                  ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function, ['A', 'B', 'C'], [14, 14, 14])
            sage: letter_game.show()
            A Co-operative Game with 3 players
            It's Characteristic Function is {('A',): 6, ('B', 'C'): 42, (): 0, ('C',): 42, ('A', 'B'): 12, ('B',): 12, ('A', 'C'): 42, ('A', 'B', 'C'): 42}
            And its Payoff Vector is [14, 14, 14]
            sage: letter_game.shapley_value()
            [2, 5, 35]
            sage: letter_game.show()
            A Co-operative Game with 3 players
            It's Characteristic Function is {('A',): 6, ('B', 'C'): 42, (): 0, ('C',): 42, ('A', 'B'): 12, ('B',): 12, ('A', 'C'): 42, ('A', 'B', 'C'): 42}
            And its Payoff Vector is [2, 5, 35]
        """
        np = self.number_players
        cf = self.char_fun
        pv = self.payoff_vector
        print "A Co-operative Game with %s players" % np
        print "It's Characteristic Function is %s" % cf
        if pv is False:
            print "And it has no Payoff Vector"
        else:
            print "And its Payoff Vector is %s" % pv
