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
        {'A': 2, 'C': 35, 'B': 5}

    We can test if it is Monotonic or Superadditive. ::

        sage: letter_game.is_monotone()
        True
        sage: letter_game.is_superadditive()
        False

    The show function will display basic information about the game. ::

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

    We can test 3 basic properties of the Payoff Vector. They are
        * Is it is efficient?
        * Does it possess the nullplayer property?
        * Does it possess the symmetry property?
    ::

        sage: letter_game.is_efficient()
        True
        sage: letter_game.nullplayer()
        True
        sage: letter_game.symmetry()
        True

    Any Payoff Vector can be passed to the game and these properties can once again be tested:

        sage: letter_game.payoff_vector = {'A': 0, 'C': 35, 'B': 3}
        sage: letter_game.is_efficient()
        False
        sage: letter_game.nullplayer()
        True
        sage: letter_game.symmetry()
        True
    """

    def __init__(self, characteristic_function, payoff_vector=False):
        r"""
        Initializes a co-operative game and checks the inputs.

        TESTS:

        An attempt to construct a game from an integer. ::

            sage: int_game = CooperativeGame(4)
            Traceback (most recent call last):
            ...
            TypeError: Characteristic function must be a dictionary

        This test checks that an error is raised when a key in the
        Characteristic Function is not a tuple. ::

            sage: tuple_function = {(): 0,
            ....:                  (1): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: tuple_game = CooperativeGame(tuple_function)
            Traceback (most recent call last):
            ...
            TypeError: Key must be a tuple

        The above test failed as `(1)` is not read as a tuple.
        """
        if type(characteristic_function) is not dict:
            raise TypeError("Characteristic function must be a dictionary")

        for key in list(characteristic_function.keys()):
            if type(key) is not tuple:
                raise TypeError("Key must be a tuple")

        self.ch_f = characteristic_function
        #self.player_list = list(characteristic_function.keys()[-1])
        self.player_list = max(characteristic_function.keys(), key = lambda key : len(key))
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
            (1, 2, 3)
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
        Shows the use of is_monotone on a simple game that is monotone. ::

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

        Shows the use of is_monotone for a game that is not monotone. ::

            sage: integer_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 10,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game.is_monotone()
            False
        """
        sets = list(self.ch_f.keys())
        return not any([set(p1) <= set(p2) and self.ch_f[p1] > self.ch_f[p2]
                                        for p1, p2 in permutations(sets, 2)])

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
        sets = list(self.ch_f.keys())
        for p1, p2 in permutations(sets, 2):
            if set(p1) & set(p2) == set():
                union = tuple(sorted(list(set(p1) | set(p2))))
                if self.ch_f[union] < self.ch_f[p1] + self.ch_f[p2]:
                    return False
        return True

    def marginal_contributions(self, player):
        contributions = []
        for pi in permutations(self.player_list):
            contributions.append(self.marginal_of_pi(player, pi))
        return contributions

    def marginal_of_pi(self, player, pi):
        predecessors, player_and_pred = self.get_predecessors(player, pi)
        if predecessors is None:
            predecessors = ()
        else:
            predecessors = tuple(predecessors)
        player_and_pred = tuple(player_and_pred)
        value = self.ch_f[player_and_pred] - self.ch_f[predecessors]
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
        Typical use of the show function with a given Payoff Vector.::

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

        Typical use of the show function with after calculating the Shapley value.::

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
        cf = self.ch_f
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

            sage: letter_function = {(): 0,
            ....:                    ('A',): 6,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function)
            sage: letter_game.payoff_vector = {'A': 10, 'B': 14, 'C': 14}
            sage: letter_game.is_efficient()
            False
        """
        return sum(self.payoff_vector.values()) == self.ch_f[self.player_list]

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
            ....:               (1,): 0,
            ....:               (2,): 12,
            ....:               (3,): 42,
            ....:               (1, 2,): 12,
            ....:               (1, 3,): 42,
            ....:               (2, 3,): 55,
            ....:               (1, 2, 3,): 55}
            sage: A_game = CooperativeGame(A_function, {1: 10, 2: 10, 3: 25})
            sage: A_game.nullplayer()
            False
        """
        allsets = list(self.ch_f.keys())
        players = [(i,) for i in self.player_list if self.ch_f[(i,)] == 0]
        nulls = []
        for player in players:
            test = [set(m) for m in allsets if player in m]
            results = []
            for sets in test:
                results.append(self.ch_f[tuple(sets.remove(player))] == self.ch_f[tuple(sets)])
            if all(results):
                nulls.append(player)

        for player in nulls:
            if self.payoff_vector[player[0]] != 0:
                return False
        return True

    def symmetry(self):
        r"""
        Returns True if the current payoff_vector possesses the symmetry property.

        EXAMPLES::

            sage: integer_function = {(): 0,
            ....:                     (1,): 12,
            ....:                     (2,): 12,
            ....:                     (3,): 42,
            ....:                     (1, 2,): 12,
            ....:                     (1, 3,): 42,
            ....:                     (2, 3,): 42,
            ....:                     (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function, {1: 2, 2: 5, 3: 35})
            sage: integer_game.symmetry()
            False
        """
        sets = list(self.ch_f.keys())
        element = [i for i in sets if len(i) == 1]
        for c1, c2 in combinations(element, 2):
            results = []
            for m in sets:
                junion = tuple(sorted(list(set(c1) | set(m))))
                kunion = tuple(sorted(list(set(c2) | set(m))))
                results.append(self.ch_f[junion] == self.ch_f[kunion])
            if all(results) and self.payoff_vector[c1[0]] == self.payoff_vector[c2[0]]:
                pass
            elif all(results):
                return False
            else:
                pass
        return True
