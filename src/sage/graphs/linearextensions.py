"""
Linear Extensions of Directed Acyclic Graphs.

A linear extension of a directed acyclic graph is a total (linear) ordering on
the vertices that is compatible with the graph in the following sense:
if there is a path from x to y in the graph, the x appears before y in the
linear extension.

The algorithm implemented in this module is from "Generating Linear Extensions
Fast" by Preusse and Ruskey, which can be found at
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.52.3057 .  The algorithm
generates the extensions in constant amortized time (CAT) -- a constant amount
of time per extension generated, or linear in the number of extensions
generated.

EXAMPLES:

Here we generate the 5 linear extensions of the following directed
acyclic graph::

    sage: from sage.graphs.linearextensions import LinearExtensions
    sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
    sage: D.is_directed_acyclic()
    True
    sage: sorted(LinearExtensions(D))
    [[0, 1, 2, 3, 4],
     [0, 1, 2, 4, 3],
     [0, 2, 1, 3, 4],
     [0, 2, 1, 4, 3],
     [0, 2, 4, 1, 3]]

Notice how all of the total orders are compatible with the ordering
induced from the graph.

We can also get at the linear extensions directly from the graph.  From
the graph, the linear extensions are known as topological sorts ::

    sage: sorted(D.topological_sort_generator())
    [[0, 1, 2, 3, 4],
     [0, 1, 2, 4, 3],
     [0, 2, 1, 3, 4],
     [0, 2, 1, 4, 3],
     [0, 2, 4, 1, 3]]


"""
#*****************************************************************************
#      Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************
from copy import copy
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

class LinearExtensions(Parent):
    def __init__(self, dag):
        r"""
        Creates an object representing the class of all linear extensions
        of the directed acyclic graph \code{dag}.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: l == loads(dumps(l))
            True

        TESTS::

            sage: list(LinearExtensions(DiGraph({ })))
            [[]]

        """
        self._dag = copy(dag)
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: LinearExtensions(D)
            Linear extensions of Digraph on 5 vertices

        """
        return "Linear extensions of %s"%self._dag

    def _prepare(self):
        r"""
        The preprocessing routine in Figure 7 of "Generating Linear
        Extensions Fast" by Preusse and Ruskey.

        This is an in-place algorithm and the list self._le keeps track
        of the current linear extensions.  The boolean variable self.is_plus
        keeps track of the "sign".

        TESTS::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: list(LinearExtensions(DiGraph({ 0:[1], 1:[0] }))) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: The digraph must be acyclic to have linear extensions.

        """
        dag_copy = copy(self._dag)
        le = []
        a = []
        b = []

        # the preprocessing routine found in Figure 7 of
        # "Generating Linear Extensions Fast" by
        # Pruesse and Ruskey
        while dag_copy.num_verts() != 0:
            # find all the minimal elements of dag_copy
            minimal_elements = dag_copy.sources()
            if not minimal_elements:
                raise ValueError("The digraph must be acyclic to have linear extensions.")
            elif len(minimal_elements) == 1:
                le.append(minimal_elements[0])
                dag_copy.delete_vertex(minimal_elements[0])
            else:
                ap = minimal_elements[0]
                bp = minimal_elements[1]
                a.append(ap)
                b.append(bp)
                le.append(ap)
                le.append(bp)
                dag_copy.delete_vertex(ap)
                dag_copy.delete_vertex(bp)

        self._max_pair = len(a) - 1
        self._le = le
        self._a = a
        self._b = b
        self._is_plus = True

    def _switch(self, i):
        """
        This implements the Switch procedure described on page 7
        of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        If i == -1, then the sign is changed.  If i > 0, then self._a[i]
        and self._b[i] are transposed.

        Note that this meant to be called by the generate_linear_extensions
        method and is not meant to be used directly.

        TESTS::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: _ = l.list()
            sage: l._le = [0, 1, 2, 3, 4]
            sage: l._is_plus
            True
            sage: [e for e in l._switch(-1)]
            []
            sage: l._is_plus
            False
            sage: l._a
            [1, 4]
            sage: l._b
            [2, 3]
            sage: [e for e in l._switch(0)]
            []
            sage: l._le
            [0, 2, 1, 3, 4]
            sage: l._a
            [2, 4]
            sage: l._b
            [1, 3]

        """
        if i == -1:
            self._is_plus = not self._is_plus
        if i >= 0:
            a, b = self._a[i], self._b[i]
            a_index = self._le.index(a)
            b_index = self._le.index(b)
            self._le[a_index] = b
            self._le[b_index] = a

            self._b[i], self._a[i] = a, b

        if self._is_plus:
            yield self._le[:]

    def _move(self, element, direction):
        """
        This implements the Move procedure described on page 7
        of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        If direction is "left", then this transposes element with the
        element on its left.  If the direction is "right", then this
        transposes element with the element on its right.

        Note that this is meant to be called by the generate_linear_extensions
        method and is not meant to be used directly.

        TESTS::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: _ = l.list()
            sage: l._le = [0, 1, 2, 3, 4]
            sage: [e for e in l._move(1, "left")]
            [[1, 0, 2, 3, 4]]
            sage: l._le
            [1, 0, 2, 3, 4]
            sage: [e for e in l._move(1, "right")]
            [[0, 1, 2, 3, 4]]
            sage: l._le
            [0, 1, 2, 3, 4]

        """
        index = self._le.index(element)
        if direction == "right":
            self._le[index] = self._le[index+1]
            self._le[index+1] = element
        elif direction == "left":
            self._le[index] = self._le[index-1]
            self._le[index-1] = element
        else:
            raise ValueError("bad direction, it should be 'right'or 'left'")
        if self._is_plus:
            yield self._le[:]


    def _right(self, i, letter):
        """
        If letter =="b", then this returns True if and only if
        self._b[i] is incomparable with the elements to its right
        in self._le.  If letter == "a", then it returns True if
        and only if self._a[i] is incomparable with the element to its
        right in self._le and the element to the right is not
        self._b[i]

        This is the Right function described on page 8 of
        "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        Note that this is meant to be called by the generate_linear_extensions
        method and is not meant to be used directly.

        TESTS::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: _ = l.list()
            sage: l._le
            [0, 1, 2, 4, 3]
            sage: l._a
            [1, 4]
            sage: l._b
            [2, 3]
            sage: l._right(0, "a")
            False
            sage: l._right(1, "a")
            False
            sage: l._right(0, "b")
            False
            sage: l._right(1, "b")
            False

        """
        if letter == "a":
            x = self._a[i]
            yindex = self._le.index(x) + 1
            if yindex >= len(self._le):
                return False
            y = self._le[yindex]
            return y != self._b[i] and self._incomparable(x, y)

        elif letter == "b":
            x = self._b[i]
            yindex = self._le.index(x) + 1
            if yindex >= len(self._le):
                return False
            y = self._le[yindex]
            return self._incomparable(x, y)
        else:
            raise ValueError("bad letter, it should be 'a' or 'b'")

    def _generate_linear_extensions(self, i):
        """
        This a Python version of the GenLE routine found in Figure 8
        of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        Note that this is meant to be called by the list
        method and is not meant to be used directly.

        TESTS::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: l._prepare()
            sage: [e for e in l._generate_linear_extensions(l._max_pair)]
            [[0, 2, 1, 3, 4]]

        """
        if i >= 0:
            for e in self._generate_linear_extensions(i-1):
                yield e
            mrb = 0
            typical = False
            while self._right(i, "b"):
                mrb += 1
                for e in self._move(self._b[i], "right"):
                    yield e
                for e in self._generate_linear_extensions(i-1):
                    yield e
                mra = 0
                if self._right(i, "a"):
                    typical = True
                    cont = True
                    while cont:
                        mra += 1
                        for e in self._move(self._a[i], "right"):
                            yield e
                        for e in self._generate_linear_extensions(i-1):
                            yield e
                        cont = self._right(i, "a")
                if typical:
                    for e in self._switch(i-1):
                        yield e
                    for e in self._generate_linear_extensions(i-1):
                        yield e
                    if mrb % 2 == 1:
                        mla = mra -1
                    else:
                        mla = mra + 1
                    for _ in range(mla):
                        for e in self._move(self._a[i], "left"):
                            yield e
                        for e in self._generate_linear_extensions(i-1):
                            yield e

            if typical and (mrb % 2 == 1):
                for e in self._move(self._a[i], "left"):
                    yield e
            else:
                for e in self._switch(i-1):
                    yield e
            for e in self._generate_linear_extensions(i-1):
                yield e
            for _ in range(mrb):
                for e in self._move(self._b[i], "left"):
                    yield e
                for e in self._generate_linear_extensions(i-1):
                    yield e

    def __iter__(self):
        """
        Iterate over the linear extensions of the directed acyclic graph.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: sorted(LinearExtensions(D))
            [[0, 1, 2, 3, 4],
             [0, 1, 2, 4, 3],
             [0, 2, 1, 3, 4],
             [0, 2, 1, 4, 3],
             [0, 2, 4, 1, 3]]
        """
        self._prepare()
        yield self._le[:]
        for e in self._generate_linear_extensions(self._max_pair):
            yield e
        for e in self._switch(self._max_pair):
            yield e
        for e in self._generate_linear_extensions(self._max_pair):
            yield e

    def _incomparable(self, x, y):
        """
        Returns True if vertices x and y are incomparable in the directed
        acyclic graph when thought of as a poset.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D); l._prepare()
            sage: l._incomparable(0,1)
            False
            sage: l._incomparable(1,2)
            True
        """
        return (not self._dag.shortest_path(x, y)) and (not self._dag.shortest_path(y, x))
