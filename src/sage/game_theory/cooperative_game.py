r"""
Co-operative games with N players.

This module implements a **basic** implementation of a characteristic function cooperative game.
The main contribution is a class for a characteristic function game.
Methods to calculate the Shapley value (a fair way of sharing common
resources: https://www.youtube.com/watch?v=aThG4YAFErw) as well as
test properties of the game (monotonicity, super additivity) are also included.

AUTHOR:

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
from sage.misc.misc import powerset
from sage.structure.sage_object import SageObject


class CooperativeGame(SageObject):
    r"""
    An object representing a co-operative game. Primarily used to compute the
    Shapley Value, but can also provide other information.

    INPUT:

    - characteristic_function - a dictionary containing all possible sets of players.
        * Key - Each set must be entered as a tuple, not a string. A single element
                tuple must end with a comma.
        * Value - A real number representing each set of players contribution.

    - payoff_vector - default = ``False``, a dictionary can be passed instead but
                      this will be overwritten if shapley_value is called.

    EXAMPLES:

    The type of game that is currently implemented is referred to as a Characteristic Function Game.
    This is a game on a set $\omega$ of players that is defined by a value function $v:C\to \mathbb{R}$
    where $C=2^{\Omega}$ is set of all coalitions of players.
    An example of such a game is shown below:

    \[
    v(c) = \begin{cases}
    0,&\text{ if }c=\emptyset\\
    6,&\text{ if }c=\{1\}\\
    12,&\text{ if }c=\{2\}\\
    42,&\text{ if }c=\{3\}\\
    12,&\text{ if }c=\{1,2\}\\
    42,&\text{ if }c=\{1,3\}\\
    42,&\text{ if }c=\{2,3\}\\
    42,&\text{ if }c=\{1, 2,3\}\\
    \end{cases}
    \]

    The function $v$ can be thought of as as a record of contribution of individuals
    and coalitions of individuals.
    Of interest, becomes how to fairly share the value of the grand coalition ($\omega$)?
    This class allows for such an answer to be formulated by calculating the Shapley value of the game.

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

    Characteristic function games can be of various types.

    A characteristic function game \(G=(N,v)\) is monotone if it satisfies \(v(C_2)\geq v(C_1) for all \(C_1\subseteq C_2\).
    A characteristic function game \(G=(N,v)\) is super-additive if it satisfies \(v(C_2)\geq v(C_1) for all \(C_1\subseteq C_2\) such that \(C_1\cap\C_2=\emptyset\).

    We can test if a game is Monotonic or Superadditive. ::

        sage: letter_game.is_monotone()
        True
        sage: letter_game.is_superadditive()
        False

    Instances have a basic representation that will display basic information about the game. ::

        sage: letter_game
        A 3 player Co-operative Game.


    It can be shown that the 'fair' payoff vector, referred to as the Shapley value is given by the following formula:

    \[
    \phi_i(G)=\frac{1}{N!}\sum_{\pi\in\Pi_n}\Delta_\pi^G(i)
    \]

    where the summation is over the permutations of the players and the marginal
    contributions of a player for a given permutation is given as:

    \[
    \Delta_\pi^G(i)=v(S_{\pi}(i)\cup i)-v(S_{\pi}(i))
    \]

    To compute the Shapley value in Sage is simple:

        sage: letter_game.shapley_value()
        {'A': 2, 'C': 35, 'B': 5}

    We can test 3 basic properties of a Payoff Vector $lambda$. They are
        * Efficiency:

        \[sum_{i=1}^N\lambda_i=v(\Omega)\]

        In other words: no value of the total coalition is lost.

        * The nullplayer property:

        If \(\exists\) \(i\) such that \(\(v(C\cup i)=v(C)\)\) for all \(C\in 2^{\Omega}\) then:

        \[\lambda_i=0\]

        In other words: if a player does not contribute to any coalition then that player should receive no payoff.

        * Does it possess the symmetry property?

        A payoff vector possesses the symmetry property if \(v(C\cup i)=v(C\cup j)\) for all \(C\in 2^{\Omega}\setminus\{i,j\}\) then:

        \[x_i=x_j\]

        If players contribute symmetrically then they should get the same payoff.

    ::

        sage: payoff_vector = letter_game.shapley_value()
        sage: letter_game.is_efficient(payoff_vector)
        True
        sage: letter_game.nullplayer(payoff_vector)
        True
        sage: letter_game.symmetry(payoff_vector)
        True

    Any Payoff Vector can be passed to the game and these properties can once again be tested:

        sage: payoff_vector = {'A': 0, 'C': 35, 'B': 3}
        sage: letter_game.is_efficient(payoff_vector)
        False
        sage: letter_game.nullplayer(payoff_vector)
        True
        sage: letter_game.symmetry(payoff_vector)
        True

    TESTS:

    Checks that the order within a key does not affect other functions. ::

        sage: letter_function = {(): 0,
        ....:                    ('A',): 6,
        ....:                    ('B',): 12,
        ....:                    ('C',): 42,
        ....:                    ('A', 'B',): 12,
        ....:                    ('C', 'A',): 42,
        ....:                    ('B', 'C',): 42,
        ....:                    ('B', 'A', 'C',): 42}
        sage: letter_game = CooperativeGame(letter_function)
        sage: letter_game.shapley_value()
        {'A': 2, 'C': 35, 'B': 5}
        sage: letter_game.is_monotone()
        True
        sage: letter_game.is_superadditive()
        False
        sage: letter_game.is_efficient({'A': 2, 'C': 35, 'B': 5})
        True
        sage: letter_game.nullplayer({'A': 2, 'C': 35, 'B': 5})
        True
        sage: letter_game.symmetry({'A': 2, 'C': 35, 'B': 5})
        True

    Any Payoff Vector can be passed to the game and these properties can once again be tested:

        sage: letter_game.is_efficient({'A': 0, 'C': 35, 'B': 3})
        False
        sage: letter_game.nullplayer({'A': 0, 'C': 35, 'B': 3})
        True
        sage: letter_game.symmetry({'A': 0, 'C': 35, 'B': 3})
        True
    """

    def __init__(self, characteristic_function):
        r"""
        Initializes a co-operative game and checks the inputs.

        TESTS:

        An attempt to construct a game from an integer. ::

            sage: int_game = CooperativeGame(4)
            Traceback (most recent call last):
            ...
            TypeError: Characteristic function must be a dictionary

        This test checks that an incorrectly entered singularly tuple will be
        changed into a tuple. In this case `(1)` becomes `(1,)`. ::

            sage: tuple_function = {(): 0,
            ....:                  (1): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: tuple_game = CooperativeGame(tuple_function)

        This test checks that if a key is not a tuple an error is raised. ::

            sage: error_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  12: 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: error_game = CooperativeGame(error_function)
            Traceback (most recent call last):
            ...
            TypeError: Key must be a tuple

        A test to ensure that the Characteristic Function is the power
        set of the grand coalition (ie all possible sub-coalitions). ::

            sage: incorrect_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: incorrect_game = CooperativeGame(incorrect_function)
            Traceback (most recent call last):
            ...
            ValueError: Characteristic Function must be the power set
        """
        if type(characteristic_function) is not dict:
            raise TypeError("Characteristic function must be a dictionary")

        self.ch_f = characteristic_function
        for key in self.ch_f:
            if len(str(key)) == 1 and type(key) is not tuple:
                self.ch_f[key, ] = self.ch_f.pop(key)
            elif type(key) is not tuple:
                raise TypeError("Key must be a tuple")
        for key in self.ch_f:
            sortedkey = tuple(sorted(list(key)))
            self.ch_f[sortedkey] = self.ch_f.pop(key)

        self.player_list = max(characteristic_function.keys(), key=lambda key: len(key))
        for coalition in powerset(self.player_list):
            if tuple(sorted(list(coalition))) not in sorted(self.ch_f.keys()):
                raise ValueError("Characteristic Function must be the power set")

        self.number_players = len(self.player_list)

    def shapley_value(self):
        r"""
        Return the payoff vector for co-operative game.

        EXAMPLES:

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
            sage: integer_game.shapley_value()
            {1: 2, 2: 5, 3: 35}
            sage: integer_game.payoff_vector
            {1: 2, 2: 5, 3: 35}

        A longer example of the shapley_value. ::

            sage: long_function = {(): 0,
            ....:                  (1,): 0,
            ....:                  (2,): 0,
            ....:                  (3,): 0,
            ....:                  (4,): 0,
            ....:                  (1, 2): 0,
            ....:                  (1, 3): 0,
            ....:                  (1, 4): 0,
            ....:                  (2, 3): 0,
            ....:                  (2, 4): 0,
            ....:                  (3, 4): 0,
            ....:                  (1, 2, 3): 0,
            ....:                  (1, 2, 4): 45,
            ....:                  (1, 3, 4): 40,
            ....:                  (2, 3, 4): 0,
            ....:                  (1, 2, 3, 4): 65}
            sage: long_game = CooperativeGame(long_function)
            sage: long_game.shapley_value()
            {1: 70/3, 2: 10, 3: 25/3, 4: 70/3}
        """
        payoff_vector = {}
        for player in self.player_list:
            player_contribution = self._marginal_contributions(player)
            average = sum(player_contribution) / len(player_contribution)
            payoff_vector[player] = average
        self.payoff_vector = payoff_vector
        return payoff_vector

    def is_monotone(self):
        r"""
        Returns True if co-operative game is monotonic.

        EXAMPLES:

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

        Shows the use of is_monotone for a longer game. ::

            sage: long_function = {(): 0,
            ....:                  (1,): 0,
            ....:                  (2,): 0,
            ....:                  (3,): 0,
            ....:                  (4,): 0,
            ....:                  (1, 2): 0,
            ....:                  (1, 3): 0,
            ....:                  (1, 4): 0,
            ....:                  (2, 3): 0,
            ....:                  (2, 4): 0,
            ....:                  (3, 4): 0,
            ....:                  (1, 2, 3): 0,
            ....:                  (1, 2, 4): 45,
            ....:                  (1, 3, 4): 40,
            ....:                  (2, 3, 4): 0,
            ....:                  (1, 2, 3, 4): 65}
            sage: long_game = CooperativeGame(long_function)
            sage: long_game.is_monotone()
            True
        """
        sets = list(self.ch_f.keys())
        return not any([set(p1) <= set(p2) and self.ch_f[p1] > self.ch_f[p2]
                        for p1, p2 in permutations(sets, 2)])

    def is_superadditive(self):
        r"""
        Returns True if co-operative game is superadditive.

        EXAMPLES:

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

        An example for is_superadditive with a longer game that returns True. ::

            sage: long_function = {(): 0,
            ....:                  (1,): 0,
            ....:                  (2,): 0,
            ....:                  (3,): 0,
            ....:                  (4,): 0,
            ....:                  (1, 2): 0,
            ....:                  (1, 3): 0,
            ....:                  (1, 4): 0,
            ....:                  (2, 3): 0,
            ....:                  (2, 4): 0,
            ....:                  (3, 4): 0,
            ....:                  (1, 2, 3): 0,
            ....:                  (1, 2, 4): 45,
            ....:                  (1, 3, 4): 40,
            ....:                  (2, 3, 4): 0,
            ....:                  (1, 2, 3, 4): 65}
            sage: long_game = CooperativeGame(long_function)
            sage: long_game.is_superadditive()
            True

        An example for is_superadditive with a longer game that returns False. ::

            sage: long_function = {(): 0,
            ....:                  (1,): 0,
            ....:                  (2,): 0,
            ....:                  (3,): 55,
            ....:                  (4,): 0,
            ....:                  (1, 2): 0,
            ....:                  (1, 3): 0,
            ....:                  (1, 4): 0,
            ....:                  (2, 3): 0,
            ....:                  (2, 4): 0,
            ....:                  (3, 4): 0,
            ....:                  (1, 2, 3): 0,
            ....:                  (1, 2, 4): 45,
            ....:                  (1, 3, 4): 40,
            ....:                  (2, 3, 4): 0,
            ....:                  (1, 2, 3, 4): 85}
            sage: long_game = CooperativeGame(long_function)
            sage: long_game.is_superadditive()
            False
        """
        sets = list(self.ch_f.keys())
        for p1, p2 in combinations(sets, 2):
            if set(p1) & set(p2) == set():
                union = tuple(sorted(list(set(p1) | set(p2))))
                if self.ch_f[union] < self.ch_f[p1] + self.ch_f[p2]:
                    return False
        return True

    def _marginal_contributions(self, player):
        r"""
        Returns a list of contributions specific to one player.

        INPUT:

        -player - A real number or string.

        EXAMPLES::

            sage: integer_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game._marginal_contributions(1)
            [6, 6, 0, 0, 0, 0]
        """
        contributions = []
        for pi in permutations(self.player_list):
            contributions.append(self.marginal_of_pi(player, pi))
        return contributions

    def marginal_of_pi(self, player, pi):
        r"""
        Returns a value for the players contribution in one permutation.

        INPUT:

        -player - A real number or string.

        -pi - A tuple which is the permutation that should be used.

        EXAMPLES::

            sage: integer_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game.marginal_of_pi(2, (2, 3, 1))
            12
        """
        predecessors, player_and_pred = self.get_predecessors(player, pi)
        if predecessors is None:
            predecessors = ()
        else:
            predecessors = tuple(predecessors)
        player_and_pred = tuple(player_and_pred)
        value = self.ch_f[player_and_pred] - self.ch_f[predecessors]
        return value

    def get_predecessors(self, player, pi):
        r"""
        Returns a list of all the predecessors of a player in a certain
        permutation and the same list including the original player
        (used elsewhere).

        INPUT:

        - player - A real number or string.

        - pi - A tuple which is the permutation that should be used.

        EXAMPLES::

            sage: integer_function = {(): 0,
            ....:                  (1,): 6,
            ....:                  (2,): 12,
            ....:                  (3,): 42,
            ....:                  (1, 2,): 12,
            ....:                  (1, 3,): 42,
            ....:                  (2, 3,): 42,
            ....:                  (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game.get_predecessors(1, (2, 3, 1))
            ([2, 3], [1, 2, 3])
            sage: integer_game.get_predecessors(2, (2, 3, 1))
            ([], [2])
            sage: integer_game.get_predecessors(3, (2, 3, 1))
            ([2], [2, 3])
        """
        pred = list(pi[:pi.index(player)])
        return sorted(pred), sorted(pred + [player])

    def _repr_(self):
        r"""

        Returns a concise description of the Game.

        EXAMPLES:

        Basic description of the game shown when calling the game instance. ::

            sage: letter_function = {(): 0,
            ....:                    ('A',): 6,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function)
            sage: letter_game
            A 3 player Co-operative Game.
        """
        np = self.number_players
        return "A %s player Co-operative Game." % np

    def _latex_(self):
        r"""

        Returns the LaTeX code representing the characteristic function.

        EXAMPLES::
        Basic description of the game shown when calling the game instance. ::

            sage: letter_function = {(): 0,
            ....:                    ('A',): 6,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function)
            sage: latex(letter_game)
            v(c) = \begin{cases}
            0,&\text{ if }c=\emptyset\\
            6,&\text{ if }c=\{A\}\\
            42,&\text{ if }c=\{C\}\\
            12,&\text{ if }c=\{B\}\\
            42,&\text{ if }c=\{B, C\}\\
            12,&\text{ if }c=\{A, B\}\\
            42,&\text{ if }c=\{A, C\}\\
            42,&\text{ if }c=\{A, B, C\}\\
            \end{cases}
        """
        np = self.number_players
        cf = self.ch_f
        output = "v(c) = \\begin{cases}\n"
        for key in sorted(cf.keys(), key=lambda key: len(key)) :
            if key == ():
                coalition = "\\emptyset"
            else:
                coalition = "\{"
                for player in key[:-1]:
                    coalition += "%s, " % player
                coalition += "%s\}" % key[-1]
            output += "%s,&\\text{ if }c=%s\\\\\n" % (cf[key], coalition)
        output += "\\end{cases}"
        return output

    def is_efficient(self, payoff_vector):
        r"""
        Returns True if the current payoff_vector is efficient.

        INPUT:

        - payoff_vector - a dictionary where the key is the player and the
                          value is their payoff.

        EXAMPLES:

        An efficient payoff_vector.::

            sage: letter_function = {(): 0,
            ....:                    ('A',): 6,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function)
            sage: letter_game.is_efficient({'A': 14, 'B': 14, 'C': 14})
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
            sage: letter_game.is_efficient({'A': 10, 'B': 14, 'C': 14})
            False

            A longer example for is_efficient. ::

            sage: long_function = {(): 0,
            ....:                  (1,): 0,
            ....:                  (2,): 0,
            ....:                  (3,): 0,
            ....:                  (4,): 0,
            ....:                  (1, 2): 0,
            ....:                  (1, 3): 0,
            ....:                  (1, 4): 0,
            ....:                  (2, 3): 0,
            ....:                  (2, 4): 0,
            ....:                  (3, 4): 0,
            ....:                  (1, 2, 3): 0,
            ....:                  (1, 2, 4): 45,
            ....:                  (1, 3, 4): 40,
            ....:                  (2, 3, 4): 0,
            ....:                  (1, 2, 3, 4): 65}
            sage: long_game = CooperativeGame(long_function)
            sage: long_game.is_efficient({1: 20, 2: 20, 3: 5, 4: 20})
            True
        """
        pl = tuple(sorted(list(self.player_list)))
        return sum(payoff_vector.values()) == self.ch_f[pl]

    def nullplayer(self, payoff_vector):
        r"""
        Returns True if the current payoff_vector possesses the null
        player property.

        INPUT:

        - payoff_vector - a dictionary where the key is the player and the
                          value is their payoff.

        EXAMPLES:

        A payoff_vector that returns True. ::

            sage: letter_function = {(): 0,
            ....:                    ('A',): 0,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function)
            sage: letter_game.nullplayer({'A': 0, 'B': 14, 'C': 14})
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
            sage: A_game = CooperativeGame(A_function)
            sage: A_game.nullplayer({1: 10, 2: 10, 3: 25})
            False

        A longer example for nullplayer. ::

            sage: long_function = {(): 0,
            ....:                  (1,): 0,
            ....:                  (2,): 0,
            ....:                  (3,): 0,
            ....:                  (4,): 0,
            ....:                  (1, 2): 0,
            ....:                  (1, 3): 0,
            ....:                  (1, 4): 0,
            ....:                  (2, 3): 0,
            ....:                  (2, 4): 0,
            ....:                  (3, 4): 0,
            ....:                  (1, 2, 3): 0,
            ....:                  (1, 2, 4): 45,
            ....:                  (1, 3, 4): 40,
            ....:                  (2, 3, 4): 0,
            ....:                  (1, 2, 3, 4): 65}
            sage: long_game = CooperativeGame(long_function)
            sage: long_game.nullplayer({1: 20, 2: 20, 3: 5, 4: 20})
            True

        TESTS:

        Checks that the function is going through all players. ::

            sage: A_function = {(): 0,
            ....:               (1,): 42,
            ....:               (2,): 12,
            ....:               (3,): 0,
            ....:               (1, 2,): 55,
            ....:               (1, 3,): 42,
            ....:               (2, 3,): 12,
            ....:               (1, 2, 3,): 55}
            sage: A_game = CooperativeGame(A_function)
            sage: A_game.nullplayer({1: 10, 2: 10, 3: 25})
            False
        """
        for player in self.player_list:
            coalitions = [coal for coal in self.ch_f if player in coal]
            results = []
            for coalit in coalitions:
                results.append(self.ch_f[coalit] == self.ch_f[tuple(
                               sorted(list(set(coalit) - {player})))])
            if all(results) and payoff_vector[player] != 0:
                return False
        return True

    def symmetry(self, payoff_vector):
        r"""
        Returns True if the current payoff_vector possesses the symmetry
        property.

        INPUT:

        - payoff_vector - a dictionary where the key is the player and the
                          value is their payoff.

        EXAMPLES:

        A Payoff Vector that returns True. ::

            sage: letter_function = {(): 0,
            ....:                    ('A',): 6,
            ....:                    ('B',): 12,
            ....:                    ('C',): 42,
            ....:                    ('A', 'B',): 12,
            ....:                    ('A', 'C',): 42,
            ....:                    ('B', 'C',): 42,
            ....:                    ('A', 'B', 'C',): 42}
            sage: letter_game = CooperativeGame(letter_function)
            sage: letter_game.symmetry({'A': 5, 'B': 14, 'C': 20})
            True

        A Payoff Vector that returns False. ::

            sage: integer_function = {(): 0,
            ....:                     (1,): 12,
            ....:                     (2,): 12,
            ....:                     (3,): 42,
            ....:                     (1, 2,): 12,
            ....:                     (1, 3,): 42,
            ....:                     (2, 3,): 42,
            ....:                     (1, 2, 3,): 42}
            sage: integer_game = CooperativeGame(integer_function)
            sage: integer_game.symmetry({1: 2, 2: 5, 3: 35})
            False

        A longer example for symmetry. ::

            sage: long_function = {(): 0,
            ....:                  (1,): 0,
            ....:                  (2,): 0,
            ....:                  (3,): 0,
            ....:                  (4,): 0,
            ....:                  (1, 2): 0,
            ....:                  (1, 3): 0,
            ....:                  (1, 4): 0,
            ....:                  (2, 3): 0,
            ....:                  (2, 4): 0,
            ....:                  (3, 4): 0,
            ....:                  (1, 2, 3): 0,
            ....:                  (1, 2, 4): 45,
            ....:                  (1, 3, 4): 40,
            ....:                  (2, 3, 4): 0,
            ....:                  (1, 2, 3, 4): 65}
            sage: long_game = CooperativeGame(long_function)
            sage: long_game.symmetry({1: 20, 2: 20, 3: 5, 4: 20})
            True
        """
        sets = list(self.ch_f.keys())
        element = [i for i in sets if len(i) == 1]
        for c1, c2 in combinations(element, 2):
            results = []
            for m in sets:
                junion = tuple(sorted(list(set(c1) | set(m))))
                kunion = tuple(sorted(list(set(c2) | set(m))))
                results.append(self.ch_f[junion] == self.ch_f[kunion])
            if (all(results) and payoff_vector[c1[0]] != payoff_vector[c2[0]]):
                return False
        return True
