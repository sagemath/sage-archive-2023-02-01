"""
Matching games.

This module implements a class for matching games (stable marriage problems)
[DI1989]_. At present the extended Gale Shapley algorithm is implemented
which can be used to obtain stable matchings.

AUTHOR:

    - James Campbell and Vince Knight 06-2014: Original version
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
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from copy import deepcopy
from sage.graphs.bipartite_graph import BipartiteGraph


class MatchingGame(SageObject):
    r"""
    An object representing a Matching Game. Includes an implementation of the
    extended Gale-Shapley algorithm.

    A matching game (also called a stable matching problem) models a situation
    of in a population of `N` suitors and `N` reviewers. Suitors and reviewers
    rank their preference and attempt to find a match.

    Formally, a matching game of size `N` is defined by two disjoint sets `S`
    and `R` of size `N`. Associated to each element of `S` and `R` is a
    preference list:

    .. MATH::

        f:S\to R^N
        \text{ and }
        g:R\to S^N

    Here is an example of matching game on 4 players:

    .. MATH::

        S = \{J, K, L, M\}\\
        R = \{A, B, C, D\}

    With preference functions:

    .. MATH::

        f(s) = \begin{cases}
        (A, D, C, B),& \text{ if } s=J\\
        (A, B, C, D),& \text{ if } s=K\\
        (B, D, C, A),& \text{ if } s=L\\
        (C, A, B, D),& \text{ if } s=M\\
        \end{cases}

        g(s) = \begin{cases}
        (L, J, K, M),& \text{ if } s=A\\
        (J, M, L, K),& \text{ if } s=B\\
        (K, M, L, J),& \text{ if } s=C\\
        (M, K, J, L),& \text{ if } s=D\\
        \end{cases}

    To implement the above game in Sage: ::

        sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
        ....:               'K': ('A', 'B', 'C', 'D'),
        ....:               'L': ('B', 'D', 'C', 'A'),
        ....:               'M': ('C', 'A', 'B', 'D')}
        sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
        ....:                 'B': ('J', 'M', 'L', 'K'),
        ....:                 'C': ('K', 'M', 'L', 'J'),
        ....:                 'D': ('M', 'K', 'J', 'L')}
        sage: m = MatchingGame([suitr_pref, reviewr_pref])
        sage: m
        A matching game with 8 players
        sage: m.suitors
        ['K', 'J', 'M', 'L']
        sage: m.reviewers
        ['A', 'C', 'B', 'D']

    A matching `M` is any bijection between `S` and `R`. If `s\in S` and
    `r\in R` are matched by `M` we denote:

    .. MATH::

        M(s)=r

    On any given matching game one intends to find a matching that is stable,
    in other words so that no one individual has an incentive to break their
    current match.

    Formally, a stable matchin is a matching that has no blocking pairs.
    A blocking pair is any pair `(s, r)` such that `M(s)\ne r` but `s` prefers
    r to `M(r)` and `r` prefers `s` to `M^{-1}(r)`.

    To obtain the stable matching in Sage we use the ``solve`` method which
    uses the extended Gale-Shapley algorithm [DI1989]_: ::

        sage: m.solve()
        {'A': ['J'],
         'C': ['K'],
         'B': ['M'],
         'D': ['L'],
         'K': ['C'],
         'J': ['A'],
         'M': ['B'],
         'L': ['D']}

    Matchings have a natural representations as bi-partite graph: ::

        sage: plot(m)

    The above plots the bi-partite graph associated with the matching.
    This plot can be accessed directly: ::

        sage: graph = m.bi_partite()
        sage: graph
        Bipartite graph on 8 vertices

    It is possible to initiate a matching game without having to name each
    suitor and reviewer: ::

        sage: n = 10
        sage: big_game = MatchingGame(n)
        sage: big_game.suitors
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: big_game.reviewers
        [-1, -2, -3, -4, -5, -6, -7, -8, -9, -10]

    If we attempt to obtain the stable matching for the above game,
    without defining the preference function we obtain an error: ::

        sage: big_game.solve()
        Traceback (most recent call last):
        ...
        ValueError: Suitor preferences are not complete

    To continue we have to populate the preference dictionary. Here
    is one example where the preferences are simply the corresponding
    element of the permutation group: ::

        sage: from itertools import permutations
        sage: suitr_preferences = list(permutations([-i-1 for i in range(n)]))
        sage: revr_preferences = list(permutations([i+1 for i in range(n)]))
        sage: for player in range(n):
        ....:     big_game.suitors[player].pref = suitr_preferences[player]
        ....:     big_game.reviewers[player].pref = revr_preferences[-player]
        sage: big_game.solve()
        {1: [-1],
         2: [-8],
         3: [-9],
         4: [-10],
         5: [-7],
         6: [-6],
         7: [-5],
         8: [-4],
         9: [-3],
        10: [-2],
        -2: [10],
        -10: [4],
        -9: [3],
        -8: [2],
        -7: [5],
        -6: [6],
        -5: [7],
        -4: [8],
        -3: [9],
        -1: [1]}


    It can be shown that the Gale Shapley algorithm will return the stable
    matching that is optimal from the point of view of the suitors and is in
    fact the worst possible matching from the point of view of the reviewers.
    To quickly obtain the matching that is optimal for the reviewers we use the
    ``solve`` method with the ``invert=True`` option: ::

        sage: left_dict = {'a': ('A', 'B', 'C'),
        ....:              'b': ('B', 'C', 'A'),
        ....:              'c': ('B', 'A', 'C')}
        sage: right_dict = {'A': ('b', 'c', 'a'),
        ....:               'B': ('a', 'c', 'b'),
        ....:               'C': ('a', 'b', 'c')}
        sage: quick_game = MatchingGame([left_dict, right_dict])
        sage: quick_game.solve()
        {'a': ['A'],
         'A': ['a'],
         'c': ['B'],
         'b': ['C'],
         'C': ['b'],
         'B': ['c']}
        sage: quick_game.solve(invert=True)
        {'a': ['B'],
         'A': ['c'],
         'c': ['A'],
         'b': ['C'],
         'C': ['b'],
         'B': ['a']}


    REFERENCES:

    .. [DI1989]  Gusfield, Dan, and Robert W. Irving.
       *The stable marriage problem: structure and algorithms.*
       Vol. 54. Cambridge: MIT press, 1989.
    """
    def __init__(self, generator):
        r"""
        Initializes a Matching Game and checks the inputs.

        TESTS:

        8 player letter game. ::

            sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
            ....:               'K': ('A', 'B', 'C', 'D'),
            ....:               'L': ('B', 'D', 'C', 'A'),
            ....:               'M': ('C', 'A', 'B', 'D')}
            sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
            ....:                 'B': ('J', 'M', 'L', 'K'),
            ....:                 'C': ('K', 'M', 'L', 'J'),
            ....:                 'D': ('M', 'K', 'J', 'L')}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.suitors
            ['K', 'J', 'M', 'L']
            sage: m.reviewers
            ['A', 'C', 'B', 'D']

        Also works for numbers. ::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])

        """
        self.suitors = []
        self.reviewers = []
        if type(generator) is Integer:
            for i in range(generator):
                self.add_suitor()
                self.add_reviewer()
        elif type(generator[0]) is dict and type(generator[1]) is dict:
            self._dict_game(generator[0], generator[1])
        else:
            raise TypeError("generator must be an integer or a list of 2 dictionaries.")

    def _dict_game(self, suitor_dict, reviwer_dict):
        r"""
        Populates the game from 2 dictionaries. One for reviewers and one for
        suitors.
        """
        for i in suitor_dict:
            self.add_suitor(i)
        for k in reviwer_dict:
            self.add_reviewer(k)

        for i in self.suitors:
            i.pref = suitor_dict[i.name]
        for k in self.reviewers:
            k.pref = reviwer_dict[k.name]

    def _repr_(self):
        r"""
        Returns a basic representation of the game stating how many players are in the game.

        EXAMPLES::

        Matching game with 2 reviewers and 2 suitors. ::

            sage: M = MatchingGame(2)
            sage: M
            A matching game with 4 players
        """
        n = len(self.reviewers) + len(self.suitors)
        return 'A matching game with %s players' % n

    def _latex_(self):
        r"""
        Creates the LaTeX representation of the dictionaries for suitors
        and reviewers.

        EXAMPLES::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: latex(g)
            Suitors
            0=\{(3, 4)\}
            1=\{(3, 4)\}
            Reviewers
            3=\{(0, 1)\}
            4=\{(1, 0)\}
        """
        suitor_dictionary, reviewer_dictionary = self.game_to_dict()
        output = "Suitors"
        for key in suitor_dictionary:
            output += "\n%s=\{%s\}" % (key, suitor_dictionary[key])
        output += "\nReviewers"
        for key in reviewer_dictionary:
            output += "\n%s=\{%s\}" % (key, reviewer_dictionary[key])
        return output

    def game_to_dict(self):
        r"""
        Obtains the preference dictionaries for a game.
        This is relevant if a game has been created by adding the preferences
        manually.

        EXAMPLES::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g.game_to_dict()
            ({0: (3, 4), 1: (3, 4)}, {3: (0, 1), 4: (1, 0)})
        """
        suitor_dictionary = {}
        reviewer_dictionary = {}
        for suitor in self.suitors:
            suitor_dictionary[suitor.name] = suitor.pref
        for reviewer in self.reviewers:
            reviewer_dictionary[reviewer.name] = reviewer.pref
        return suitor_dictionary, reviewer_dictionary

    def plot(self):
        r"""
        Creates the plot representing the stable matching for the game.
        Note that the game must be solved for this to work.

        EXAMPLES::

        An error is returned if the game is not solved:

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: plot(g)
            Traceback (most recent call last):
            ...
            ValueError: Game has not been solved yet

            sage: g.solve()
            {0: [3], 1: [4], 3: [0], 4: [1]}
            sage: plot(g)
        """
        pl = self.bi_partite()
        return pl.plot()

    def bi_partite(self):
        r"""
        Constructs a ``BipartiteGraph`` Object of the game.
        This method
        """
        self._is_solved()

        sol_dict = self._sol_dict()
        graph = BipartiteGraph(sol_dict)
        return graph

    def _is_solved(self):
        r"""
        Checks if the Game has been solved yet.
        """
        suitor_check = all(s.partner for s in self.suitors)
        reviewer_check = all(r.partner for r in self.reviewers)
        if not suitor_check or not reviewer_check:
            raise ValueError("Game has not been solved yet")

    def _is_complete(self):
        r"""
        Checks that all players have acceptable preferences.

        TESTS:

        Not enough reviewers. ::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Must have the same number of reviewers as suitors

        Suitors preferences make no sense. ::

            sage: suit = {0: (3, 8),
            ....:         1: (0, 0)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Suitor preferences are not complete

        """
        if len(self.suitors) != len(self.reviewers):
            raise ValueError("Must have the same number of reviewers as suitors")

        for suitor in self.suitors:
            if set(suitor.pref) != set(self.reviewers):
                raise ValueError("Suitor preferences are not complete")

        for reviewer in self.reviewers:
            if set(reviewer.pref) != set(self.suitors):
                raise ValueError("Reviewer preferences are not complete")

    def add_suitor(self, name=False):
        r"""
        Adds a suitor to the game.

        INPUTS:

        -``name`` - Can be a string or a number. If left blank will automatically
                    generate an integer.
        """
        if name is False:
            name = len(self.suitors) + 1
        new_suitor = _Player(name, 'suitor', len(self.reviewers))
        self.suitors.append(new_suitor)
        for r in self.reviewers:
            r.pref = [-1 for s in self.suitors]

    def add_reviewer(self, name=False):
        r"""
        Adds a reviewer to the game.

        INPUTS:

        -``name`` - Can be a string or number. If left blank will automatically
                    generate an integer.
        """
        if name is False:
            name = -len(self.reviewers) - 1
        new_reviewer = _Player(name, 'reviewer', len(self.suitors))
        self.reviewers.append(new_reviewer)
        for s in self.suitors:
            s.pref = [-1 for r in self.reviewers]

    def _sol_dict(self):
        r"""
        Creates a dictionary of the stable matching. Keys are the player,
        values are their partner as a single element list. This is to allow
        the creation of ``BipartiteGraph``.
        """
        self._is_solved()

        sol_dict = {}
        for s in self.suitors:
            sol_dict[s] = [s.partner]
        for r in self.reviewers:
            sol_dict[r] = [r.partner]
        return sol_dict

    def solve(self, invert=False):
        r"""
        Computes a stable matching for the game using the Gale-Shapley
        algorithm.
        """
        self._is_complete()

        for s in self.suitors:
            s.partner = False
        for r in self.reviewers:
            r.partner = False

        if invert:
            reviewers = deepcopy(self.suitors)
            suitors = deepcopy(self.reviewers)
        else:
            suitors = deepcopy(self.suitors)
            reviewers = deepcopy(self.reviewers)

        while len([s for s in suitors if s.partner is False]) != 0:
            s = [s for s in suitors if s.partner is False][0]
            r = next((x for x in reviewers if x.name == s.pref[0]), None)
            if r.partner is False:
                r.partner = s
                s.partner = r
            elif r.pref.index(s.name) < r.pref.index(r.partner.name):
                r.partner.partner = False
                r.partner = s
                s.partner = r
            else:
                s.pref = s.pref[1:]

        if invert:
            suitors, reviewers = reviewers, suitors

        for i, j in zip(self.suitors, suitors):
            i.partner = j.partner
        for i, j in zip(self.reviewers, reviewers):
            i.partner = j.partner

        return self._sol_dict()


class _Player():
    def __init__(self, name, player_type, len_pref):
        self.name = name
        self.type = player_type
        self.pref = [-1 for i in range(len_pref)]
        self.partner = False

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return repr(self.name)

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()
