# -*- coding: utf-8 -*-
r"""
Finite filtered complexes

AUTHORS:

- Guillaume Rousseau (2021-05)

This module implements the basic structures of finite filtered complexes.
A filtered complex is a simplicial complex, where each simplex is given
a weight, or "filtration value", such that the weight of a simplex is
greater than the weight of each of its faces.

The algorithm used in this module comes from [ZC2005]_.

EXAMPLES::

    sage: FilteredSimplicialComplex([([0], 0), ([1], 0), ([0, 1], 1)])
    Filtered complex on vertex set (0, 1) and with simplices ((0,) : 0), ((1,) : 0), ((0, 1) : 1)

Sage can compute persistent homology of simplicial complexes::

    sage: X = FilteredSimplicialComplex([([0], 0), ([1], 0), ([0, 1], 1)])
    sage: X.persistence_intervals(0)
    [(0, 1), (0, +Infinity)]

FilteredSimplicialComplex objects are mutable. Filtration values can be
set with the ``filtration`` method as follows::

    sage: X = FilteredSimplicialComplex() # returns an empty complex
    sage: X.persistence_intervals(1)
    []
    sage: X.filtration(Simplex([0, 2]), 0) # recursively adds faces
    sage: X.filtration(Simplex([0, 1]), 0)
    sage: X.filtration(Simplex([1, 2]), 0)
    sage: X.filtration(Simplex([0, 1, 2]), 1) # closes the circle
    sage: X.persistence_intervals(1)
    [(0, 1)]

The filtration value of a simplex can be accessed as well with the
``filtration`` method, by not specifying a filtration value in
the arguments. If the simplex is not in the complex, this returns
None::

    sage: X = FilteredSimplicialComplex([([0], 0), ([1], 0), ([0,1], 1)])
    sage: X.filtration(Simplex([0]))
    0
    sage: X.filtration(Simplex([1,2]))
    <BLANKLINE>

Filtration values can be accessed with function call and list
syntax as follows::

    sage: X = FilteredSimplicialComplex([([0], 0), ([1], 0), ([0,1], 1)])
    sage: s_1 = Simplex([0])
    sage: X[s_1]
    0
    sage: X(Simplex([0,1]))
    1
    sage: X(Simplex(['baba']))
    <BLANKLINE>

It is also possible to set the filtration value of a simplex with
the ``insert`` method, which takes as argument a list of vertices
rather than a ``Simplex``. This can make code more readable / clear::

    sage: X = FilteredSimplicialComplex()
    sage: X.insert(['a'], 0)
    sage: X.insert(['b', 'c'], 1)
    sage: X
    Filtered complex on vertex set ('a', 'b', 'c') and with simplices (('a',) : 0), (('c',) : 1), (('b',) : 1), (('b', 'c') : 1)

"""

# ****************************************************************************
#       Copyright (C) 2013 GUILLAUME ROUSSEAU <guillaume.rousseau@ens-lyon.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.homology.simplicial_complex import Simplex, SimplicialComplex
from sage.modules.free_module import FreeModule
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.infinity import infinity
from sage.misc.cachefunc import cached_method

class FilteredSimplicialComplex(SageObject):

    def __init__(self, simplices=[], warnings=False):
        r"""
        Define a filtered complex.

        :param simplices: list of simplices and filtration values
        :param warnings: set to True for more info on each insertion

        ``simplices`` should be a list of tuples ``(l, v)`` where
        ``l`` is a list of vertices and ``v`` is the corresponding filtration
        value.


        If ``warnings`` is set to True, the user will be warned when they try
        to change the filtration value of a simplex.


        EXAMPLES::

            sage: FilteredSimplicialComplex([([0], 0), ([1], 0), ([2], 1), ([0,1], 2.27)])
            Filtered complex on vertex set (0, 1, 2) and with simplices ((0,) : 0), ((1,) : 0), ((2,) : 1), ((0, 1) : 2.27000000000000)

        """

        # _vertices is the set of vertices on which the complex
        # is constructed
        self._vertices = set()

        # _filtration_dict has simplices as keys
        # and entries are corresponding filtration values
        self._filtration_dict = {}
        self._dimension = 0
        self._max_value = 0

        # when _warnings is set to True, insertion
        # will warn the user when something non-trivial
        # happens.
        self._warnings = warnings

        # Insert all simplices in the initial list
        for l, v in simplices:
            self.insert(l, v)

    def _get_value(self, s):
        r"""
        Return the filtration value of a simplex in the complex.

        :param s: simplex



        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 1), ([1], 2)])
            sage: X._get_value(Simplex([0]))
            1

        """
        if s in self._filtration_dict:
            return self._filtration_dict[s]
        else:
            return None

    def _insert(self, simplex, filtration_value):
        r"""
        Add a simplex to the complex.

        All faces of the simplex are added recursively if they are
        not already present, with the same value.
        If the simplex is already present, and the new value is lower
        than its current value in the complex, the value gets updated,
        otherwise it does not change. This propagates recursively to faces.
        If warnings have been enabled, this method will describe what it
        is doing during an insertion.

        :param simplex: simplex to be inserted
        :type simplex: Simplex
        :param filtration_value: value of the simplex

        EXAMPLES::

            sage: X = FilteredSimplicialComplex()
            sage: X._insert(Simplex([0]),3)
            sage: X
            Filtered complex on vertex set (0,) and with simplices ((0,) : 3)

        """
        # Keep track of whether the simplex is already in the complex
        # and if it should be updated or not
        update = False
        curr_value = self[simplex]
        if curr_value is not None:
            if self._warnings:
                print("Face {} is already in the complex.".format(simplex))
            if curr_value > filtration_value:
                if self._warnings:
                    warning_string = "However its value is {}".format(curr_value)
                    warning_string += ": updating it to {}".format(filtration_value)
                    print(warning_string)
                update = True
                self._filtration_dict.pop(simplex)
            else:
                if self._warnings:
                    print("Its value is {}: keeping it that way".format(curr_value))
                return

        # check that all faces are in the complex already.
        # If not, warn the user and add faces (recursively)
        faces = simplex.faces()
        if simplex.dimension() > 0:
            for f in faces:
                if self._warnings:
                    print("Also inserting face {} with value {}".format(f, filtration_value))
                self._insert(f, filtration_value)

        self._filtration_dict[simplex] = filtration_value
        self._dimension = max(self._dimension, simplex.dimension())
        self._max_value = max(self._max_value, filtration_value)
        self._vertices.update(simplex.set())
        self._persistent_homology.clear_cache()

    def insert(self, vertex_list, filtration_value):
        r"""
        Add a simplex to the complex.

        All faces of the simplex are added recursively if they are
        not already present, with the same value.
        If the simplex is already present, and the new value is lower
        than its current value in the complex, the value gets updated,
        otherwise it does not change. This propagates recursively to faces.
        If warnings have been enabled, this method will describe what it
        is doing during an insertion.

        This method calls self._insert, turning the list of vertices
        into a simplex.

        :param vertex_list: list of vertices
        :param filtration_value: desired value of the simplex to be
            added.

        EXAMPLES::

            sage: X = FilteredSimplicialComplex()
            sage: X.insert(Simplex([0]),3)
            sage: X
            Filtered complex on vertex set (0,) and with simplices ((0,) : 3)

        If the warning parameter was set to true, this method will print
        some info::

            sage: X = FilteredSimplicialComplex(warnings = True)
            sage: X.insert(Simplex([0, 1]), 2)
            Also inserting face (1,) with value 2
            Also inserting face (0,) with value 2
            sage: X.insert(Simplex([0]),1)
            Face (0,) is already in the complex.
            However its value is 2: updating it to 1
            sage: X.insert(Simplex([0]), 77)
            Face (0,) is already in the complex.
            Its value is 1: keeping it that way

        """
        self._insert(Simplex(vertex_list), filtration_value)

    def filtration(self, s, filtration_value=None):
        """
        Set filtration value of a simplex, or return value
        of existing simplex.

        :param s: Simplex for which to set or return value
        :param filtration_value: Optional, filtration value
            for the simplex.

        If no filtration value is specified, this function calls
        ``_get_value`` and hence returns the value of the simplex
        in the complex. If the simplex is not in the complex, this
        returns None.

        If ``filtration_value`` is set, this function inserts the
        simplex into the complex with the specified value.
        See documentation of ``insert`` for more details.

        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1], 1)])
            sage: X.filtration(Simplex([0, 1]))
            <BLANKLINE>
            sage: X.filtration(Simplex([0, 1]), 2)
            sage: X.filtration(Simplex([0, 1]))
            2
        """
        if filtration_value is None:
            return self._get_value(s)
        else:
            self._insert(s, filtration_value)

    def prune(self,threshold):
        """
        Return a copy of the filtered complex, where simplices above
        the threshold value have been removed.

        :param threshold: a real value, above which simplices are
            discarded.

        Simplices with filtration value exactly equal to ``threshold``
        are kept in the result.

        EXAMPLES::

            sage: a = FilteredSimplicialComplex()
            sage: a.insert([0], 0)
            sage: a.insert([0, 1], 1)
            sage: a.insert([0, 2], 2)
            sage: b = a.prune(1)
            sage: b
            Filtered complex on vertex set (0, 1) and with simplices ((0,) : 0), ((1,) : 1), ((0, 1) : 1)

        """
        result_complex = FilteredSimplicialComplex()
        for s in self._filtration_dict:
            if self[s] <= threshold:
                result_complex._insert(s, self[s])

        return result_complex

    @cached_method(key=lambda self,f,s,v:(f,s))
    def _persistent_homology(self, field=2, strict=True, verbose=False):
        """
        Compute the homology intervals of the complex.

        :param field: prime number modulo which homology is computed
            Default value: 2
        :param strict: if set to False, takes into account intervals
            of persistence 0. Default value: True
        :param verbose: if set to True, prints progress of computation.
            Default value: False.

        This method is called whenever betti numbers or intervals are
        computed, and the result is cached. It returns the list of
        intervals.

        ALGORITHM:

        The computation behind persistent homology is a matrix reduction.
        The algorithm implemented is described in [ZC2005]_.

        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1], 0), ([0,1], 2)])
            sage: X._persistent_homology()[0]
            [(0, 2), (0, +Infinity)]

        Some homology elements may have a lifespan or persistence of 0.
        They are usually discarded, but can be kept if necessary::

            sage: X = FilteredSimplicialComplex()
            sage: X.insert([0,1],1) # opens a hole and closes it instantly
            sage: X._persistent_homology(strict=False)[0]
            [(1, 1), (1, +Infinity)]

        REFERENCES:

        [ZC2005]_

        TESTS:

        This complex is used as a running example in [ZC2005]_:

            sage: l = [([0], 0), ([1], 0), ([2], 1), ([3], 1), ([0, 1], 1), ([1, 2], 1), ([0, 3], 2), ([2, 3], 2), ([0, 2], 3), ([0, 1, 2], 4), ([0, 2, 3], 5)]
            sage: X = FilteredSimplicialComplex(l)
            sage: X.persistence_intervals(0)
            [(0, 1), (1, 2), (0, +Infinity)]
            sage: X.persistence_intervals(1)
            [(3, 4), (2, 5)]
            sage: X.persistence_intervals(0, strict=False)
            [(0, 1), (1, 1), (1, 2), (0, +Infinity)]

        """

        # first, order the simplices in lexico order
        # on dimension, value and then arbitrary order
        # defined by the Simplex class.
        def key(s):
            d = self._get_value(s)
            return (s.dimension(), d, s)
        simplices = list(self._filtration_dict)
        simplices.sort(key=key)

        # remember the index of each simplex in a dict
        self._index_of_simplex = {}
        n = len(simplices)
        for i in range(n):
            self._index_of_simplex[simplices[i]] = i

        self._field_int = field
        self._field = GF(field)
        self._chaingroup = FreeModule(self._field, rank_or_basis_keys=simplices)

        # Initialize data structures for the algo
        self._marked = [False for i in range(n)]
        self._T = [None for i in range(n)]
        intervals = [[] for i in range(self._dimension+1)]
        self.pairs = []

        self._strict = strict
        self._verbose = verbose

        if self._verbose:
            print("Beginning first pass")

        for j in range(n):
            # if being verbose, print progress
            # every 1000 simplices.
            if self._verbose and j % 1000 == 0:
                print('{}/{}'.format(j, n))

            s = simplices[j]
            d = self._remove_pivot_rows(s, simplices)

            if d == 0:
                self._marked[j] = True
            else:
                max_index = self._max_index(d)
                t = simplices[max_index]
                self._T[max_index] = (s, d)
                self._add_interval(t, s, intervals)

        if self._verbose:
            print("First pass over, beginning second pass")

        for j in range(n):
            if j % 1000 == 0 and self._verbose:
                print('{}/{}'.format(j, n))

            s = simplices[j]
            if self._marked[j] and not self._T[j]:
                self._add_interval(s, None, intervals)

        if self._verbose:
            print("Second pass over")
        return intervals
        

    def _add_interval(self, s, t, intervals):
        r"""
        Add a new interval (i.e. homology element).
        This method should not be called by users, it is used in
        the ``_compute_persistence`` method. The simplex of
        death may be None, in which case the interval is infinite.

        :param s: birth simplex
        :param t: death simplex
        :param intervals: list of current intervals

        If ``t`` is not None, its dimension should be
        1 more than the dimension of ``s``.

        TESTS::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1, 2], 10)])
            sage: int_list = X._persistent_homology()
            sage: int_list[0]
            [(0, +Infinity), (10, +Infinity)]
            sage: X._add_interval(Simplex([0]), Simplex([1, 2]),int_list)
            sage: int_list[0]
            [(0, +Infinity), (10, +Infinity), (0, 10)]

        Infinite interval:
            
            sage: int_list2 = [[],[]]
            sage: X._add_interval(Simplex([1, 2]), None, int_list2)
            sage: int_list2[1]
            [(10, +Infinity)]

        """

        # figure out dimension of homology element
        # and indices of the two simplices. If the
        # closing simplex is None, then the interval
        # is infinite.
        k = s.dimension()
        i = self._filtration_dict[s]
        if not t:
            j = infinity
        else:
            j = self._filtration_dict[t]

        # Only add intervals of length 0 if
        # strict mode is not enabled.
        if i != j or (not self._strict):
            intervals[k].append((i, j))
            self.pairs.append((s, t))

    def _remove_pivot_rows(self, s, simplices):
        r"""
        Return the boundary chain of a simplex,
        from which pivot elements have been removed.

        This method implements the subroutine of the
        same name in [ZC2005]_. This method should not
        be called by users, it is used in the
        ``compute_persistence`` method.

        TESTS::

            sage: l = [([0], 0), ([1], 0), ([2], 1), ([3], 1), ([0, 1], 1), ([1, 2], 1), ([0, 3], 2), ([2, 3], 2), ([0, 2], 3), ([0, 1, 2], 4)]
            sage: X = FilteredSimplicialComplex(l)
            sage: X._persistent_homology()
            [[(0, 1), (1, 2), (0, +Infinity)], [(3, 4), (2, +Infinity)], []]
            sage: X._remove_pivot_rows(Simplex([0,1,2]), list(X._filtration_dict))
            0
            sage: X.insert([0,2,3],5)
            sage: X._remove_pivot_rows(Simplex([0,2,3]), list(X._filtration_dict))
            B[(2, 3)]

        """

        d = self._chaingroup()
        # Handle the case when the simplex is a vertex
        if s.dimension() == 0:
            return d

        # Initialize the boundary chain
        for (i, f) in enumerate(s.faces()):
            d += (-1)**i * self._chaingroup(f)

        # Remove all unmarked elements
        for (s, x_s) in d:
            j = self._index_of_simplex[s]
            if not self._marked[j]:
                d = d - x_s*self._chaingroup(s)

        # Reduce d until it is empty or until the simplex
        # with maximum index in the complex among all
        # non-zero terms is not in T.
        while d != 0:
            max_index = self._max_index(d)
            t = simplices[max_index]

            if not self._T[max_index]:
                break

            c = self._T[max_index][1]
            q = c[t]
            d = d - ((q**(-1))*c)

        return d

    def _max_index(self, d):
        r"""
        Return the maximal index of all simplices with nonzero
        coefficient in ``d``.

        This method is called in ``_remove_pivot_rows`` and
        ``compute_persistence``. It should not be called by users
        outside of those methods.

        TESTS::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1], 5), ([0, 1], 18), ([0, 2, 3], 32)])
            sage: X._persistent_homology()
            [[(5, 18), (0, +Infinity)], [], []]
            sage: a = X._chaingroup(Simplex([0, 1]))
            sage: b = X._chaingroup(Simplex([0, 3]))
            sage: d = a + b
            sage: X._max_index(d)
            6
        """
        currmax = -1
        for (s, x_s) in d:
            j = self._index_of_simplex[s]
            if j > currmax:
                currmax = j
        return currmax

    def persistence_intervals(self, dimension, field=2, strict=True):
        """
        Return the list of d-dimensional homology elements.

        :param dimension: integer, dimension for which to
            return intervals.
        :param field: prime number, modulo which persistent
            homology is computed. Default value 2.
        :param strict: if set to False, takes into account intervals
            of persistence 0. Default value: True

        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1], 1), ([0,1], 2)])
            sage: X.persistence_intervals(0)
            [(1, 2), (0, +Infinity)]

        """
        intervals = self._persistent_homology(field, strict, verbose=False)
        if dimension < len(intervals):
            return intervals[dimension][:]
        else:
            return []

    def betti_number(self, k, a, b, field=2, strict=True):
        r"""
        Return the ``k``-dimensional betti number from ``a`` to ``a+b``
   
        :param field: prime number, modulo which persistent
            homology is computed. Default value 2.
        :param strict: if set to False, takes into account intervals
            of persistence 0. Default value: True

        The betti number ``\beta_k^{a,a+b}`` counts the number of
        homology elements which are alive throughout the whole
        duration ``[a, a+b]``.

        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1], 0), ([0,1], 2)])
            sage: X.betti_number(0, 0.5, 1)
            2
            sage: X.betti_number(0, 1.5, 1)
            1

        If an element vanishes at time ``a+b`` exactly,
        it does not count towards the betti number::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1], 0), ([0,1], 2)])
            sage: X.betti_number(0, 1.5, 0.5)
            1

        """
        intervals = self._persistent_homology(field, strict, verbose=False)
        res = 0
        for (i, j) in intervals[k]:
            if (i <= a and a + b < j) and a >= 0:
                    res += 1
        return res

    def __call__(self, s):
        r"""
        Return the filtration value of a simplex

        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 1)])
            sage: X(Simplex([0]))
            1

        """
        return self._get_value(s)

    def __getitem__(self, s):
        r"""
        Return the filtration value of a simplex

        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 1)])
            sage: X[Simplex([0])]
            1

        """
        return self._get_value(s)

    def _repr_(self):
        """
        Print representation.

        If there are too many simplices or vertices, only prints the
        count for each.

        EXAMPLES::

            sage: X = FilteredSimplicialComplex([([0], 0), ([1], 0), ([0, 1], 1)])
            sage: X
            Filtered complex on vertex set (0, 1) and with simplices ((0,) : 0), ((1,) : 0), ((0, 1) : 1)
            sage: X.insert([0, 1, 2, 3, 4], 8)
            sage: X
            Filtered complex on 5 vertices and with 31 simplices

        """

        vert_count = len(self._vertices)
        simp_count = len(self._filtration_dict)
        if simp_count > 10 or vert_count > 10:
            vertex_string = "on {} vertices".format(vert_count)
            simplex_string = "with {} simplices".format(simp_count)
        else:
            vertex_string = "on vertex set {}".format(tuple(sorted(self._vertices)))
            simplex_string = "with simplices "
            simplex_list = ["({} : {})".format(s, self._filtration_dict[s]) for s in self._filtration_dict]
            simplex_string += ", ".join(simplex_list)

        return "Filtered complex " + vertex_string + " and " + simplex_string

    def _simplicial_(self):
        """
        Return the associated simplicial complex

        All simplices of the filtered simplicial complex are
        included in the resulting simplicial complex.

        EXAMPLES::

            sage: l = [([0],0), ([1],0), ([2],1), ([3],1), ([0, 1],1), ([1, 2],1), ([0, 3],2), ([2, 3],2), ([0, 2],3), ([0, 1, 2],4), ([0, 2, 3],5)]
            sage: a = FilteredSimplicialComplex(l)
            sage: b = SimplicialComplex(a)
            sage: b
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (0, 2, 3)}

        """
        return SimplicialComplex(self._filtration_dict)