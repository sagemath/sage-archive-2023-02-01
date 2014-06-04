r"""
Co-operative games with N players.

This module implements characteristic function cooperative games.
The main contribution is a class for a characteristic function game.
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
from itertools import permutations, combinations
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

    - payoff_vector - default = ``False``, a dictionary can be passed instead but
                      this will be overwritten if shapley_value is called.

    EXAMPLES::
    Basic example of how to implement a co-operative game. These functions will
    be used repeatedly in other examples. ::

        sage: integer_function = {(): 0,
        ....:                     (1,): 6,
        ....:                     (2,): 12,
        ....:                     (3,): 42,
        ....:                     (1, 2,): 12,
        ....:                     (1, 3,): 42,
        ....:                     (2, 3,): 42,
        ....:                     (1, 2, 3,): 42}
        sage: integer_game = CooperativeGame(integer_function)

    We can also use strings instead of numbers. ::

        sage: letter_function = {(): 0,
        ....:                    ('A',): 6,
        ....:                    ('B',): 12,
        ....:                    ('C',): 42,
        ....:                    ('A', 'B',): 12,
        ....:                    ('A', 'C',): 42,
        ....:                    ('B', 'C',): 42,
        ....:                    ('A', 'B', 'C',): 42}
        sage: letter_game = CooperativeGame(letter_function)

    From this we can now compute the Shapley Value. ::

        sage: letter_game.shapley_value()
        {'A': 2, 'B': 5, 'C': 35}

    We can test if it is Monotonic or Superadditive. ::

        sage: letter_game.is_monotone()
        True
        sage: letter_game.is_superadditive()
        True

    The show function will give us basic information about the game. ::

        sage: letter_game.show()
        A Co-operative Game with 3 players
        Characteristic Function is
             ('A',) : 6
             ('B', 'C') : 42
             () : 0
             ('C',) : 42
             ('A', 'B') : 12
             ('B',) : 12
             ('A', 'C') : 42
             ('A', 'B', 'C') : 42
        Payoff vector is {'A': 2, 'C': 35, 'B': 5}

    Finally, we can test 3 basic properties of the Payoff Vector. They are,
        * Is it is efficient?
        * Does it possess the nullplayer property?
        * Does it possesss the symmetry property? ::

        sage: letter_game.is_efficient()
        True
        sage: letter_game.nullplayer()
        True
        sage: letter_game.symmetry()
        True
    """

    def __init__(self, characteristic_function, payoff_vector=False):
        self.char_fun = characteristic_function
        self.player_list = list(characteristic_function.keys()[-1])
        self.number_players = len(self.player_list)
        self.payoff_vector = payoff_vector

    def shapley_value(self):
        r"""
        Return the payoff vector for co-operative game.

        EXAMPLES::
        A typical example of the use of shapley_value. ::

            sage: integer_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game.player_list
            [1, 2, 3]
            sage: integer_game.payoff_vector
            False
            sage: integer_game.shapley_value()
            {1: 2, 2: 5, 3: 35}
            sage: integer_game.payoff_vector
            {1: 2, 2: 5, 3: 35}
        """
        payoff_vector = {}
        for player in self.player_list:
            player_contribution = self.marginal_contributions(player)
            average = sum(player_contribution) / len(player_contribution)
            payoff_vector[player] = average
        self.payoff_vector = payoff_vector
        return payoff_vector

    def is_monotone(self):
        r"""
        Returns True if co-operative game is monotonic.

        EXAMPLES::
        Shows the use of is_monotone on a simple game. ::

            sage: integer_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game.is_monotone()
            True
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

        EXAMPLES::
        An example that returns False. ::

            sage: integer_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game.is_superadditive()
            False

        An example that returns True. ::

            sage: A_function = {(): 0,
            ....:               (1,): 6,
            ....:               (2,): 12,
            ....:               (3,): 42,
            ....:               (1, 2,): 18,
            ....:               (1, 3,): 48,
            ....:               (2, 3,): 55,
            ....:               (1, 2, 3,): 80}
            sage: A_game = CooperativeGame(A_function)
            sage: A_game.is_superadditive()
            True
        """
        sets = list(self.char_fun.keys())
        for i, k in permutations(sets, 2):
            if set(i) & set(k) != set():
                pass
            else:
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
            ....:                    ('A',): 6,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function, {'A': 14, 'B': 14, 'C': 14})
            sage: letter_game.show()
            A Co-operative Game with 3 players
            Characteristic Function is
                 ('A',) : 6
                 ('B', 'C') : 42
                 () : 0
                 ('C',) : 42
                 ('A', 'B') : 12
                 ('B',) : 12
                 ('A', 'C') : 42
                 ('A', 'B', 'C') : 42
            Payoff vector is {'A': 14, 'C': 14, 'B': 14}
            sage: letter_game.shapley_value()
            {'A': 2, 'C': 35, 'B': 5}
            sage: letter_game.show()
            A Co-operative Game with 3 players
            Characteristic Function is
                 ('A',) : 6
                 ('B', 'C') : 42
                 () : 0
                 ('C',) : 42
                 ('A', 'B') : 12
                 ('B',) : 12
                 ('A', 'C') : 42
                 ('A', 'B', 'C') : 42
            Payoff vector is {'A': 2, 'C': 35, 'B': 5}
        """
        np = self.number_players
        cf = self.char_fun
        pv = self.payoff_vector
        print "A Co-operative Game with %s players" % np
        print "Characteristic Function is"
        for key in cf:
            print "\t %s : %s" %(key, cf[key])
        if pv is False:
            print "And it has no payoff vector"
        else:
            print "Payoff vector is %s" % pv

    def is_efficient(self):
        r"""
        Returns True if the current payoff_vector is efficient.

        EXAMPLES::
        An efficient payoff_vector.::

            sage: letter_function = {(): 0,
            ....:                    ('A',): 6,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function, {'A': 14, 'B': 14, 'C': 14})
            sage: letter_game.is_efficient()
            True
        """
        if sum(self.payoff_vector.values()) == self.char_fun[tuple(self.player_list)]:
            return True
        else:
            return False

    def nullplayer(self):
        r"""
        Returns True if the current payoff_vector possesses the null player property.

        EXAMPLES::
        A payoff_vector that returns True. ::

            sage: letter_function = {(): 0,
            ....:                    ('A',): 0,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function, {'A': 0, 'B': 14, 'C': 14})
            sage: letter_game.nullplayer()
            True

        A payoff_vector that returns False. ::

            sage: A_function = {(): 0,
            ....:               (1,): 6,
            ....:               (2,): 12,
            ....:               (3,): 42,
            ....:               (1, 2,): 12,
            ....:               (1, 3,): 42,
            ....:               (2, 3,): 55,
            ....:               (1, 2, 3,): 55}
            sage: A_game = CooperativeGame(A_function, {1: 10, 2: 10, 3: 25})
            sage: A_game.nullplayer()
            True
        """
        sets = list(self.char_fun.keys())
        player = [(i,) for i in self.player_list if self.char_fun[(i,)] == 0]
        nulls = []
        for k in player:
            test = [set(m) for m in sets if k in m]
            results = []
            for j in test:
                results.append(self.char_fun[tuple(j.remove(k))] == self.char_fun[tuple(j)])
            if all(results):
                nulls.append(k)
            else:
                pass
        for i in nulls:
            if self.payoff_vector[i[0]] != 0:
                return False
            else:
                pass
        return True

    def symmetry(self):
        r"""
        Returns True if the current payoff_vector possesses the symmetry property.

        EXAMPLES::

            sage: integer_function = {(): 0,
            ....:                     (1,): 5,
            ....:                     (2,): 5,
            ....:                     (3,): 42,
            ....:                     (1, 2,): 12,
            ....:                     (1, 3,): 42,
            ....:                     (2, 3,): 42,
            ....:                     (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function, {1: 2, 2: 5, 3: 35})
            sage: integer_game.symmetry()
            False
        """
        sets = list(self.char_fun.keys())
        element = [i for i in sets if len(i) == 1]
        other = [i for i in sets]
        for j, k in combinations(element, 2):
            for m in other:
                junion = tuple(set(j) | set(m))
                kunion = tuple(set(k) | set(m))
                results = []
                results.append(self.char_fun[junion] == self.char_fun[kunion])
                if all(results) and self.payoff_vector[j[0]] == self.payoff_vector[k[0]]:
                    pass
                elif all(results):
                    return False
                else:
                    pass
        return True

    def is_additivity(self):
        r"""
        Returns True if the current payoff_vector possesses the additivity property.
        """
