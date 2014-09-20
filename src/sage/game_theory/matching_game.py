"""
Matching games.

This module implements a class for matching games (stable marriage problems)
[DI1989]_. At present the extended Gale-Shapley algorithm is implemented
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
    in a population of `N` suitors and `N` reviewers. Suitors and reviewers
    rank their preferences and attempt to find a match.

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
        A matching game with 4 suitors and 4 reviewers
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

    Formally, a stable matching is a matching that has no blocking pairs.
    A blocking pair is any pair `(s, r)` such that `M(s)\ne r` but `s` prefers
    r to `M(r)` and `r` prefers `s` to `M^{-1}(r)`.

    To obtain the stable matching in Sage we use the ``solve`` method which
    uses the extended Gale-Shapley algorithm [DI1989]_: ::

        sage: m.solve()
        {'K': ['C'],
         'J': ['A'],
         'M': ['B'],
         'L': ['D']}

    Matchings have a natural representations as bi-partite graph: ::

        sage: plot(m)

    The above plots the bi-partite graph associated with the matching.
    This plot can be accessed directly: ::

        sage: graph = m.bipartite()
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
        10: [-2]}


    It can be shown that the Gale-Shapley algorithm will return the stable
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
         'c': ['B'],
         'b': ['C']}
        sage: quick_game.solve(invert=True)
        {'A': ['c'],
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

        Can create a game from an integer which then requires
        a bespoke creation of preferences::

            sage: g = MatchingGame(3)
            sage: g
            A matching game with 3 suitors and 3 reviewers

        """
        self.suitors = []
        self.reviewers = []
        if type(generator) is Integer:
            for i in range(generator):
                self.add_suitor()
                self.add_reviewer()
        elif type(generator[0]) is dict and type(generator[1]) is dict:
            for i in generator[0]:
                self.add_suitor(i)
            for k in generator[1]:
                self.add_reviewer(k)

            for i in self.suitors:
                i.pref = generator[0][i.name]
            for k in self.reviewers:
                k.pref = generator[1][k.name]
        else:
            raise TypeError("generator must be an integer or a list of 2 dictionaries.")

        self.suitor_dictionary = {}
        self.reviewer_dictionary = {}
        for suitor in self.suitors:
            self.suitor_dictionary[suitor.name] = suitor.pref
        for reviewer in self.reviewers:
            self.reviewer_dictionary[reviewer.name] = reviewer.pref

    def _repr_(self):
        r"""
        Returns a basic representation of the game stating how many players are in the game.

        EXAMPLES:

        Matching game with 2 reviewers and 2 suitors. ::

            sage: M = MatchingGame(2)
            sage: M
            A matching game with 2 suitors and 2 reviewers
        """
        size = (len(self.reviewers), len(self.suitors))
        return 'A matching game with %s suitors and %s reviewers' % size

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
            0\to\{(3, 4)\}
            1\to\{(3, 4)\}
            Reviewers
            3\to\{(0, 1)\}
            4\to\{(1, 0)\}
        """
        output = "Suitors"
        for key in self.suitor_dictionary:
            output += "\n%s\\to\{%s\}" % (key, self.suitor_dictionary[key])
        output += "\nReviewers"
        for key in self.reviewer_dictionary:
            output += "\n%s\\to\{%s\}" % (key, self.reviewer_dictionary[key])
        return output

    def plot(self):
        r"""
        Creates the plot representing the stable matching for the game.
        Note that the game must be solved for this to work.

        EXAMPLES:

        An error is returned if the game is not solved::

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
            {0: [3], 1: [4]}
            sage: plot(g)
        """
        pl = self.bipartite()
        return pl.plot()

    def bipartite(self):
        r"""
        Constructs a ``BipartiteGraph`` Object of the game.
        This method is similar to the plot method.
        Note that the game must be solved for this to work.

        EXAMPLES:

        An error is returned if the game is not solved::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g.bipartite()
            Traceback (most recent call last):
            ...
            ValueError: Game has not been solved yet

            sage: g.solve()
            {0: [3], 1: [4]}
            sage: g.bipartite()
            Bipartite graph on 4 vertices
        """
        self._is_solved()

        graph = BipartiteGraph(self.sol_dict)
        return graph

    def _is_solved(self):
        r"""
        Raises an error if the Game has been solved yet.

        EXAMPLES:

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_solved()
            Traceback (most recent call last):
            ...
            ValueError: Game has not been solved yet
            sage: g.solve()
            {0: [3], 1: [4]}
            sage: g._is_solved()
        """
        suitor_check = all(s.partner for s in self.suitors)
        reviewer_check = all(r.partner for r in self.reviewers)
        if not suitor_check or not reviewer_check:
            raise ValueError("Game has not been solved yet")

    def _is_complete(self):
        r"""
        Raises an error if all players do not have acceptable preferences.

        EXAMPLES:

        Not enough reviewers. ::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Must have the same number of reviewers as suitors

        Not enough suitors. ::

            sage: suit = {0: (3, 4)}
            sage: revr = {1: (0, 2),
            ....:         3: (0, 1)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Must have the same number of reviewers as suitors

        Suitors preferences are incomplete. ::

            sage: suit = {0: (3, 8),
            ....:         1: (0, 0)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Suitor preferences are not complete

        Reviewer preferences are incomplete. ::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 2, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Reviewer preferences are not complete
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

        -``name`` -- Can be a string or a number; if left blank will automatically
         generate an integer.

        EXAMPLES:

        Creating a two player game::

            sage: g = MatchingGame(2)
            sage: g.suitors
            [1, 2]

        Adding a suitor without specifying a name::

            sage: g.add_suitor()
            sage: g.suitors
            [1, 2, 3]

        Adding a suitor while specifying a name::

            sage: g.add_suitor('D')
            sage: g.suitors
            [1, 2, 3, 'D']

        Note that now our game is no longer complete::

            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Must have the same number of reviewers as suitors

        Note that an error is raised if one tries to add a suitor
        with a name that already exists::

            sage: g.add_suitor('D')
            Traceback (most recent call last):
            ...
            ValueError: A suitor with name: D already exists
        """
        if name is False:
            name = len(self.suitors) + 1
        if name in [s.name for s in self.suitors]:
            raise ValueError("A suitor with name: %s already exists" % name)

        new_suitor = _Player(name, 'suitor', len(self.reviewers))
        self.suitors.append(new_suitor)
        for r in self.reviewers:
            r.pref = [-1 for s in self.suitors]

    def add_reviewer(self, name=False):
        r"""
        Adds a reviewer to the game.

        INPUTS:

        -``name`` -- Can be a string or number; if left blank will automatically
         generate an integer.

        EXAMPLES:

        Creating a two player game::

            sage: g = MatchingGame(2)
            sage: g.reviewers
            [-1, -2]

        Adding a suitor without specifying a name::

            sage: g.add_reviewer()
            sage: g.reviewers
            [-1, -2, -3]

        Adding a suitor while specifying a name::

            sage: g.add_reviewer(10)
            sage: g.reviewers
            [-1, -2, -3, 10]

        Note that now our game is no longer complete::

            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: Must have the same number of reviewers as suitors

        Note that an error is raised if one tries to add a reviewer
        with a name that already exists::

            sage: g.add_reviewer(10)
            Traceback (most recent call last):
            ...
            ValueError: A reviewer with name: 10 already exists
        """
        if name is False:
            name = -len(self.reviewers) - 1
        if name in [s.name for s in self.reviewers]:
            raise ValueError("A reviewer with name: %s already exists" % name)

        new_reviewer = _Player(name, 'reviewer', len(self.suitors))
        self.reviewers.append(new_reviewer)
        for s in self.suitors:
            s.pref = [-1 for r in self.reviewers]

    def solve(self, invert=False):
        r"""
        Computes a stable matching for the game using the Gale-Shapley
        algorithm.

        TESTS::

            sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
            ....:               'K': ('A', 'B', 'C', 'D'),
            ....:               'L': ('B', 'C', 'D', 'A'),
            ....:               'M': ('C', 'A', 'B', 'D')}
            sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
            ....:                 'B': ('J', 'M', 'L', 'K'),
            ....:                 'C': ('M', 'K', 'L', 'J'),
            ....:                 'D': ('M', 'K', 'J', 'L')}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve()
            {'K': ['D'], 'J': ['A'], 'M': ['C'], 'L': ['B']}

            sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
            ....:               'K': ('A', 'B', 'C', 'D'),
            ....:               'L': ('B', 'C', 'D', 'A'),
            ....:               'M': ('C', 'A', 'B', 'D')}
            sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
            ....:                 'B': ('J', 'M', 'L', 'K'),
            ....:                 'C': ('M', 'K', 'L', 'J'),
            ....:                 'D': ('M', 'K', 'J', 'L')}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve(invert=True)
            {'A': ['L'], 'C': ['M'], 'B': ['J'], 'D': ['K']}

            sage: suitr_pref = {1: (-1,)}
            sage: reviewr_pref = {-1: (1,)}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve()
            {1: [-1]}

            sage: suitr_pref = {}
            sage: reviewr_pref = {}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve()
            {}
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

        while not all(s.partner for s in suitors):
            s = [s for s in suitors if s.partner is False][0]
            r = next((x for x in reviewers if x == s.pref[0]), None)
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

        self.sol_dict = {}
        for s in self.suitors:
            self.sol_dict[s] = [s.partner]
        for r in self.reviewers:
            self.sol_dict[r] = [r.partner]

        if invert:
            return {key:self.sol_dict[key] for key in self.reviewers}
        return {key:self.sol_dict[key] for key in self.suitors}


class _Player(SageObject):
    r"""
    A class to act as a data holder for the players used of the matching games
    These instances are used when initiating players and to keep track of
    whether or not partners have a preference.
    """
    def __init__(self, name, player_type, len_pref):
        r"""
        TESTS::

            sage: from sage.game_theory.matching_game import _Player
            sage: p = _Player(10, 'suitor', 3)
            sage: p
            10
            sage: p.pref
            []
            sage: p.partner
            False
        """
        self.name = name
        self.type = player_type
        self.pref = []
        self.partner = False

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.game_theory.matching_game import _Player
            sage: p = _Player(10, 'suitor', 3)
            sage: d = {p : (1, 2, 3)}
            sage: d
            {10: (1, 2, 3)}
        """
        return hash(self.name)

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.game_theory.matching_game import _Player
            sage: p = _Player(10, 'suitor', 3)
            sage: p
            10

            sage: p = _Player('Karl', 'reviewer', 2)
            sage: p
            'Karl'
        """
        return repr(self.name)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from sage.game_theory.matching_game import _Player
            sage: p = _Player(10, 'suitor', 3)
            sage: q = _Player('Karl', 'reviewer', 2)
            sage: p == q
            False

            sage: from sage.game_theory.matching_game import _Player
            sage: p = _Player(10, 'suitor', 3)
            sage: q = _Player(10, 'reviewer', 2)
            sage: p == q
            True
        """
        return self.__repr__() == other.__repr__()
