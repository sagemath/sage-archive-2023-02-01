# distutils: language = c++
"""
Dancing Links internal pyx code

EXAMPLES::

    sage: from sage.combinat.matrices.dancing_links import dlx_solver
    sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
    sage: x = dlx_solver(rows)
    sage: x
    Dancing links solver for 6 columns and 6 rows

The number of solutions::

    sage: x.number_of_solutions()
    3

.. WARNING:: 
 
    The way it is coded, it can be iterated through only once. The
    next call to the above function gives wrong result::

        sage: x.number_of_solutions()
        0

We recreate the dancing links object and we find all solutions::

    sage: x = dlx_solver(rows)
    sage: sorted(x.solutions_iterator())
    [[0, 1], [2, 3], [4, 5]]

Return the first solution found when the computation is done in parallel::

    sage: sorted(x.first_solution_found_in_parallel(ncpus=8))
    [0, 1]

Find all solutions using some specific rows::

    sage: x_using_row_2 = x.restrict([2])
    sage: x_using_row_2
    Dancing links solver for 7 columns and 6 rows
    sage: list(x_using_row_2.solutions_iterator())
    [[2, 3]]

The two basic methods that are wrapped in this class are ``search`` which
returns ``1`` if a solution is found or ``0`` otherwise and ``get_solution``
which return the current solution::

    sage: x = dlx_solver(rows)
    sage: x.search()
    1
    sage: x.get_solution()
    [0, 1]
    sage: x.search()
    1
    sage: x.get_solution()
    [2, 3]
    sage: x.search()
    1
    sage: x.get_solution()
    [4, 5]
    sage: x.search()
    0
"""
#*****************************************************************************
#       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>
#       Copyright (C) 2015-2017 Sébastien Labbé <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from cpython.object cimport PyObject_RichCompare
from libcpp.vector cimport vector
from cysignals.signals cimport sig_on, sig_off

cdef extern from "dancing_links_c.h":
    cdef cppclass dancing_links:
        vector[int] solution
        int number_of_columns()
        void add_rows(vector[vector[int]] rows)
        int search()


cdef class dancing_linksWrapper:
    r"""
    A simple class that implements dancing links.

    The main methods to list the solutions are :meth:`search` and
    :meth:`get_solution`. You can also use :meth:`number_of_solutions` to count
    them.

    This class simply wraps a C++ implementation of Carlo Hamalainen.
    """
    cdef dancing_links _x
    cdef list _rows

    def __init__(self, rows):
        """
        Initialize our wrapper (self._x) as an actual C++ object.

        We must pass a list of rows at start up. There are no methods
        for resetting the list of rows, so this class acts as a one-time
        executor of the C++ code.

        TESTS::

            sage: rows = [[0,1,2], [1, 2]]
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: x = dlx_solver(rows)
            sage: x
            Dancing links solver for 3 columns and 2 rows
            sage: type(x)
            <type 'sage.combinat.matrices.dancing_links.dancing_linksWrapper'>
        """
        self._init_rows(rows)

    def __repr__(self):
        """
        The string representation of this wrapper is just the list of
        rows as supplied at startup.

        TESTS::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0]]
            sage: dlx_solver(rows)
            Dancing links solver for 3 columns and 3 rows
        """
        return "Dancing links solver for {} columns and {} rows".format(
                self.ncols(), self.nrows())

    def rows(self):
        r"""
        Return the list of rows.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0]]
            sage: x = dlx_solver(rows)
            sage: x.rows()
            [[0, 1, 2], [1, 2], [0]]
        """
        return self._rows

    def ncols(self):
        """
        Return the number of columns.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0]]
            sage: dlx = dlx_solver(rows)
            sage: dlx.ncols()
            3
        """
        return self._x.number_of_columns()

    def nrows(self):
        """
        Return the number of rows.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0]]
            sage: dlx = dlx_solver(rows)
            sage: dlx.nrows()
            3
        """
        return len(self._rows)

    def __reduce__(self):
        """
        This is used when pickling.

        TESTS::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: X = dlx_solver(rows)
            sage: X == loads(dumps(X))
            1
            sage: rows += [[2]]
            sage: Y = dlx_solver(rows)
            sage: Y == loads(dumps(X))
            0
        """
        return type(self), (self._rows,)

    def __richcmp__(dancing_linksWrapper left, dancing_linksWrapper right, int op):
        """
        Two dancing_linksWrapper objects are equal if they were
        initialised using the same row list.

        TESTS::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: X = dlx_solver(rows)
            sage: Z = dlx_solver(rows)
            sage: rows += [[2]]
            sage: Y = dlx_solver(rows)
            sage: X == Z
            1
            sage: X == Y
            0
        """
        return PyObject_RichCompare(left._rows, right._rows, op)

    def _init_rows(self, rows):
        """
        Initialize our instance of dancing_links with the given rows.

        This is for internal use by dlx_solver only.

        TESTS:

        This doctest tests ``_init_rows`` vicariously! ::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print(x.search())
            1

        The following example would crash in Sage's debug version
        from :trac:`13864` prior to the fix from :trac:`13882`::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: x = dlx_solver([])          # indirect doctest
            sage: x.get_solution()
            []

        """
        cdef vector[int] v
        cdef vector[vector[int]] vv

        self._rows = [row for row in rows]

        for row in self._rows:
            v.clear()

            for x in row:
                v.push_back(x)

            vv.push_back(v)

        sig_on()
        self._x.add_rows(vv)
        sig_off()

    def get_solution(self):
        """
        Return the current solution.

        After a new solution is found using the method :meth:`search` this
        method return the rows that make up the current solution.

        TESTS::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print(x.search())
            1
            sage: print(x.get_solution())
            [3, 0]
        """
        cdef size_t i
        cdef list s = []
        for i in range(self._x.solution.size()):
            s.append(self._x.solution.at(i))

        return s

    def search(self):
        """
        Search for a new solution.

        Return ``1`` if a new solution is found and ``0`` otherwise. To recover
        the solution, use the method :meth:`get_solution`.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print(x.search())
            1
            sage: print(x.get_solution())
            [3, 0]

        TESTS:

        Test that :trac:`11814` is fixed::

            sage: dlx_solver([]).search()
            0
            sage: dlx_solver([[]]).search()
            0

        If search is called once too often, it keeps returning 0::

            sage: x = dlx_solver([[0]])
            sage: x.search()
            1
            sage: x.search()
            0
            sage: x.search()
            0
        """
        sig_on()
        x = self._x.search()
        sig_off()
        return x

    def restrict(self, indices):
        r"""
        Return a dancing links solver solving the subcase which uses some
        given rows.

        For every row that is wanted in the solution, we add a new column
        to the row to make sure it is in the solution.

        INPUT:

        - ``indices`` -- list, row indices to be found in the solution

        OUTPUT:

            dancing links solver

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: d
            Dancing links solver for 6 columns and 6 rows
            sage: sorted(d.solutions_iterator())
            [[0, 1], [2, 3], [4, 5]]

        To impose that the 0th row is part of the solution, the rows of the new
        problem are::

            sage: d_using_0 = d.restrict([0])
            sage: d_using_0
            Dancing links solver for 7 columns and 6 rows
            sage: d_using_0.rows()
            [[0, 1, 2, 6], [3, 4, 5], [0, 1], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]

        After restriction the subproblem has one more columns and the same
        number of rows as the original one::

            sage: d.restrict([1]).rows()
            [[0, 1, 2], [3, 4, 5, 6], [0, 1], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]
            sage: d.restrict([2]).rows()
            [[0, 1, 2], [3, 4, 5], [0, 1, 6], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]

        This method allows to find solutions where the 0th row is part of a
        solution::

            sage: map(sorted, d.restrict([0]).solutions_iterator())
            [[0, 1]]

        Some other examples::

            sage: map(sorted, d.restrict([1]).solutions_iterator())
            [[0, 1]]
            sage: map(sorted, d.restrict([2]).solutions_iterator())
            [[2, 3]]
            sage: map(sorted, d.restrict([3]).solutions_iterator())
            [[2, 3]]

        Here there are no solution using both 0th and 3rd row::

            sage: list(d.restrict([0,3]).solutions_iterator())
            []

        TESTS::

            sage: d.restrict([]).rows()
            [[0, 1, 2], [3, 4, 5], [0, 1], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]
        """
        from copy import deepcopy
        rows = deepcopy(self._rows)
        ncols = self.ncols()
        for i,row_index in enumerate(indices):
            rows[row_index].append(ncols+i)
        return dlx_solver(rows)

    def split(self, column):
        r"""
        Return a dict of independent solvers.

        For each ``i``-th row containing a ``1`` in the ``column``, the
        dict associates the solver giving all solution using the ``i``-th
        row. 

        This is used for parallel computations.

        INPUT:

        - ``column`` -- integer, the column used to split the problem into
          independent subproblems

        OUTPUT:

            dict where keys are row numbers and values are dlx solvers

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: d
            Dancing links solver for 6 columns and 6 rows
            sage: sorted(d.solutions_iterator())
            [[0, 1], [2, 3], [4, 5]]

        After the split each subproblem has one more column and the same
        number of rows as the original problem::

            sage: D = d.split(0)
            sage: D
            {0: Dancing links solver for 7 columns and 6 rows,
             2: Dancing links solver for 7 columns and 6 rows,
             4: Dancing links solver for 7 columns and 6 rows}

        The (disjoint) union of the solutions of the subproblems is equal to the
        set of solutions shown above::

            sage: for x in D.values(): list(x.solutions_iterator())
            [[0, 1]]
            [[2, 3]]
            [[4, 5]]

        TESTS::

            sage: d.split(6)
            Traceback (most recent call last):
            ...
            ValueError: column(=6) must be in range(ncols) where ncols=6

        """
        if not 0 <= column < self.ncols():
            raise ValueError("column(={}) must be in range(ncols) "
                             "where ncols={}".format(column, self.ncols()))
        indices = [i for (i,row) in enumerate(self._rows) if column in row]
        return {i:self.restrict([i]) for i in indices}

    def solutions_iterator(self):
        r"""
        Return an iterator of the solutions.

        .. WARNING::

            This function can be used only once. To iterate through the
            solutions another time, one needs to recreate the dlx solver.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: list(d.solutions_iterator())
            [[0, 1], [2, 3], [4, 5]]

        As warned above, it can be used only once::

            sage: list(d.solutions_iterator())
            []
        """
        while self.search():
            yield self.get_solution()

    def first_solution_found_in_parallel(self, ncpus=1, column=None):
        r"""
        Return the first solution found after spliting the problem to
        allow parallel computation.

        Usefull when it is very hard just to find one solution to a given
        problem.

        INPUT:

        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time
        - ``column`` -- integer (default: ``None``), the column used to split
          the problem, if ``None`` a random column is chosen

        OUTPUT:

        list of rows or ``None`` if no solution is found

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: sorted(d.first_solution_found_in_parallel())
            [0, 1]
            sage: sorted(d.first_solution_found_in_parallel(ncpus=8))
            [0, 1]
            sage: sorted(d.first_solution_found_in_parallel(ncpus=8, column=4))
            [0, 1]

        When no solution is found::

            sage: rows = [[0,1,2], [2,3,4,5], [0,1,2,3]]
            sage: d = dlx_solver(rows)
            sage: [d.first_solution_found_in_parallel(column=i) for i in range(6)]
            [None, None, None, None, None, None]

        """
        if column is None:
            from random import randrange
            column = randrange(self.ncols())
        D = self.split(column)

        from sage.parallel.decorate import parallel
        @parallel(ncpus=ncpus)
        def first_solution(i):
            dlx = D[i]
            if dlx.search():
                return dlx.get_solution()
            else:
                return None

        K = sorted(D)
        for ((args, kwds), val) in first_solution(K):
            if not val is None:
                return val

    def _number_of_solutions_iterator(self, ncpus=1, column=0):
        r"""
        Return an iterator over the number of solutions using each row
        containing a ``1`` in the given ``column``.

        INPUT:

        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time
        - ``column`` -- integer (default: ``0``), the column used to split
          the problem

        OUTPUT:

            iterator of tuples (row number, number of solutions)

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: sorted(d._number_of_solutions_iterator(ncpus=2, column=3))
            [(1, 1), (3, 1), (5, 1)]

        ::

            sage: S = Subsets(range(5))
            sage: rows = [list(x) for x in S]
            sage: d = dlx_solver(rows)
            sage: d.number_of_solutions()
            52
            sage: sum(b for a,b in d._number_of_solutions_iterator(ncpus=2, column=3))
            52
        """
        D = self.split(column)
        from sage.parallel.decorate import parallel

        @parallel(ncpus=ncpus)
        def nb_sol(i):
            return D[i].number_of_solutions()
        K = sorted(D)
        for ((args, kwds), val) in nb_sol(K):
            yield args[0], val

    def number_of_solutions(self, ncpus=1, column=0):
        r"""
        Return the number of distinct solutions.

        INPUT:

        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time. If `ncpus>1` the dancing
          links problem is split into independent subproblems to
          allow parallel computation.
        - ``column`` -- integer (default: ``0``), the column used to split
          the problem (ignored if ``ncpus`` is ``1``)

        OUTPUT:

        integer

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows += [[0,2]]
            sage: rows += [[1]]
            sage: rows += [[3]]
            sage: x = dlx_solver(rows)
            sage: x.number_of_solutions()
            2

        ::

            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: x = dlx_solver(rows)
            sage: x.number_of_solutions(ncpus=2, column=3)
            3

        TESTS::

            sage: dlx_solver([]).number_of_solutions()
            0
        """
        cdef int N = 0
        if ncpus == 1:
            while self.search():
                N += 1
            return N
        else:
            it = self._number_of_solutions_iterator(ncpus, column)
            return sum(val for (k,val) in it)

def dlx_solver(rows):
    """
    Internal-use wrapper for the dancing links C++ code.

    EXAMPLES::

        sage: from sage.combinat.matrices.dancing_links import dlx_solver
        sage: rows = [[0,1,2]]
        sage: rows+= [[0,2]]
        sage: rows+= [[1]]
        sage: rows+= [[3]]
        sage: x = dlx_solver(rows)
        sage: print(x.search())
        1
        sage: print(x.get_solution())
        [3, 0]
        sage: print(x.search())
        1
        sage: print(x.get_solution())
        [3, 1, 2]
        sage: print(x.search())
        0
    """
    return dancing_linksWrapper(rows)


def make_dlxwrapper(s):
    """
    Create a dlx wrapper from a Python *string* s.

    This was historically used in unpickling and is kept for backwards
    compatibility. We expect s to be ``dumps(rows)`` where rows is the
    list of rows used to instantiate the object.

    TESTS::

        sage: from sage.combinat.matrices.dancing_links import make_dlxwrapper
        sage: rows = [[0,1,2]]
        sage: x = make_dlxwrapper(dumps(rows))
        sage: print(x.__str__())
        Dancing links solver for 3 columns and 1 rows
    """
    from sage.structure.sage_object import loads
    return dancing_linksWrapper(loads(s))
