r"""
Module of trace monoids (free partially commutative monoids).

EXAMPLES:

The following example demonstrates a monoid creation::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a'))); M
    Trace monoid on 3 generators ([a], [b], [c]) over Free monoid on 3 generators (a, b, c) with independence relation {(a, c), (c, a)}

Different monoid elements can be equal because of partially commutative multiplication::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
    sage: c*a*b == a*c*b
    True

Lets ensure that it is a monoid::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
    sage: M in Monoids()
    True

AUTHORS:

- Pavlo Tokariev (2019-05-31): initial version

"""
# ****************************************************************************
#       Copyright (C) 2019 Pavlo Tokariev <pavlo.tokariev@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import print_function

import operator
from collections import OrderedDict
from itertools import repeat, chain, product

from sage.misc.cachefunc import cached_method

from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.monoid import Monoid_class
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.sets.set import Set
from sage.structure.element import MonoidElement
from sage.structure.parent import Set_generic
from sage.structure.unique_representation import UniqueRepresentation


class TraceMonoidElement(MonoidElement):
    """
    Element of a trace monoid. Also known as a trace.

    Elements of trace monoid is actually a equivalence classes
    of related free monoid over some equivalence relation
    that in the case is presented as independence relation.

    EXAMPLES:

    Basic TraceMonoid elements usage::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
        sage: M.<a,b,c,d> = TraceMonoid(I=I)
        sage: x = b*a*d*a*c*b
        sage: x^3
        [b*a^2*d*b^2*c*a^2*d*b^2*c*a^2*d*b*c]
        sage: x^0
        1
        sage: x.lex_normal_form()
        b*a^2*d*b*c
        sage: x.foata_normal_form()
        (b, a*d, a, b*c)
    """

    def __init__(self, M, x):
        MonoidElement.__init__(self, M)
        self._repr = x

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: hash(a)
            1914282862589934403  # 64-bit
            139098947            # 32-bit
            sage: hash(b)
            2996819001369607946  # 64-bit
            13025034             # 32-bit
            sage: hash(a*b)
            7114093379175463612  # 64-bit
            2092317372           # 32-bit
        """
        return hash(self._repr)

    def _repr_(self):
        """
        Textual representation of a trace.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: a*b
            [a*b]
            sage: b*a
            [b*a]
            sage: d*a
            [a*d]
        """
        if self == self.parent().one():
            return "1"
        return "[{}]".format(self._repr)

    def lex_normal_form(self):
        """
        Returns lexicographic normal form of the trace.

        OUTPUT: free monoid element

        SEEALSO::

            :meth:`sage.monoids.trace_monoid.TraceMonoid._compute_lex_normal_form`

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: (a*b).lex_normal_form()
            a*b
            sage: (b*a).lex_normal_form()
            b*a
            sage: (d*a).lex_normal_form()
            a*d
        """
        return self._repr

    def foata_normal_form(self):
        """
        Returns Foata normal form of the trace.

        OUTPUT: tuple of free monoid elements

        SEEALSO::

            :meth:`sage.monoids.trace_monoid.TraceMonoid._compute_foata_normal_form`

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.foata_normal_form()
            (b, a*d, a, b*c)
        """
        return self.parent()._compute_foata_normal_form(self._repr)

    def _mul_(self, other):
        r"""
        Concatenates one equivalence class with another.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: a*b*c == a*c*b
            True
        """
        return self.parent(self._repr * other._repr)

    def _flat_elements(self):
        r"""
        Return flatten list of generator numbers representing the trace.

        OUTPUT: list of generator indexes

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a^3*d*a*c*b^2
            sage: x._flat_elements()
            [b, a, a, a, a, d, b, b, c]
        """
        return [g for g, times in self._repr for _ in range(times)]

    @cached_method
    def dependence_graph(self):
        r"""
        Return dependence graph of the trace.

        It is a directed graph where all dependent (non-commutative)
        generators are connected by edges which
        direction depend on the generator position in the trace.

        OUTPUT: directed graph of generator indexes

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.dependence_graph()
            Digraph on 6 vertices
        """
        elements = self._flat_elements()
        independence = self.parent()._independence
        graph = {}

        for i, e in enumerate(elements):
            edges = []
            for v in graph:
                if (e, elements[v]) not in independence:
                    edges.append((v, i))
            graph[i] = []
            for v1, v2 in edges:
                graph[v1].append(v2)

        return DiGraph(graph)

    @cached_method
    def hasse_diagram(self, alg="naive"):
        r"""
        Return Hasse diagram of the trace.

        Hasse diagram is a dependence graph without transitive edges.

        INPUT:

        - ``alg`` -- string (default: ``'naive'``); defines algorithm that will be used
          to compute Hasse diagram; there are two variants: ``'naive'`` and ``'min'``.

        OUTPUT: directed graph of generator indexes

        SEEALSO::

            :meth:`sage.monoids.trace_monoid.TraceMonoidElement.naive_hasse_digram`,
            :meth:`sage.monoids.trace_monoid.TraceMonoidElement.min_hasse_diagram`.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.hasse_diagram()
            Digraph on 6 vertices

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.hasse_diagram(alg='naive') == x.hasse_diagram(alg='min')
            True
            sage: y = b*a^3*d*a*c*b^2
            sage: y.hasse_diagram(alg='naive') == y.hasse_diagram(alg='min')
            True
        """
        if alg == "naive":
            return self.naive_hasse_diagram()
        elif alg == "min":
            return self.min_hasse_diagram()
        else:
            raise ValueError("`alg` option must be `naive` "
                             "or `min`, got `{}`.".format(alg))

    def min_hasse_diagram(self):
        r"""
        Return Hasse diagram of the trace.

        OUTPUT: directed graph of generator indexes

        SEEALSO::

            :meth:`sage.monoids.trace_monoid.TraceMonoidElement.hasse_digram`,
            :meth:`sage.monoids.trace_monoid.TraceMonoidElement.naive_hasse_diagram`.
        """
        elements = self._flat_elements()
        elements.reverse()
        independence = self.parent()._independence
        reachable = dict()
        min = set()
        graph = DiGraph({})

        for i, x in enumerate(elements):
            reachable[i] = set()
            front = min.copy()
            while front:
                used = set()
                for j in list(front):
                    y = elements[j]
                    if (x, y) not in independence:
                        graph.add_edge(i, j)
                        reachable[i].add(j)
                        reachable[i].update(reachable[j])
                        if j in min:
                            min.remove(j)
                        used.add(j)
                forbidden = set(chain.from_iterable(reachable[v] for v in used))
                front = set(dest for _, dest in graph.outgoing_edges(front, labels=False))
                front = front - forbidden

            min.add(i)

        length = len(elements)
        graph.relabel(length - 1 - i for i in range(length))
        return graph

    def naive_hasse_diagram(self):
        r"""
        Return Hasse diagram of the trace.

        ALGORITHM:

        In loop check for every two pair of edges
        if they have common vertex, remove their transitive edge.

        OUTPUT: directed graph of generator indexes

        SEEALSO::

            :meth:`sage.monoids.trace_monoid.TraceMonoidElement.hasse_digram`,
            :meth:`sage.monoids.trace_monoid.TraceMonoidElement.min_hasse_diagram`.
        """
        d = self.dependence_graph()
        h = d.copy()

        for e1 in d.edges():
            for e2 in d.edges():
                if e1[1] == e2[0]:
                    h.delete_edge((e1[0], e2[1]))

        return h

    def _richcmp_(self, other, op):
        r"""
        Compare two traces by their lexicographic normal forms.

        ALGORITHM:

        Transform each trace to lexicographic form and then compare them
        lexicographically in terms of free monoid.

        OUTPUT: 0 if self == other, -1 if self < other, 1 if self > other

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: a^2 > a
            True
            sage: a*b < b*a
            True
            sage: a*c*b == a*b*c
            True
        """
        return self._repr._richcmp_(other._repr, op)

    def alphabet(self):
        r"""
        Return alphabet of a trace.

        OUTPUT: set of free monoid generators

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.alphabet()
            {b, a, d, c}
        """
        return Set([g for g, _ in self._repr])

    def projection(self, letters):
        r"""
        Return a trace that formed from ``self`` by erasing ``letters``

        INPUT:

        - ``letters`` -- set of generators; defines set of letters that will be
          used to filter the trace.

        OUTPUT: a trace

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = M(b*a*d*a*c*b)
            sage: x.projection({a,b})
            [b*a^2*b]
            sage: x.projection({b,d,c})
            [b*d*b*c]
        """
        return self.parent(reduce(
            operator.mul,
            [x for x in self._flat_elements() if x in letters],
            self.parent()._base.one()
        ))


class TraceMonoid(Monoid_class, UniqueRepresentation):
    r"""
    Return a free partially commuting monoid (trace monoid) on `n` generators
    over independence relation `I`.

    We construct a trace monoid by specifing:

    - a free monoid and independence relation
    - or generator names and independence relation,
      FreeMonoid is constructed automatically then

    INPUT:

    - ``M`` -- a free monoid

    - ``I`` -- commutation relation between generators
      (or their names if the ``names`` are given)

    - ``names`` -- names of generators

    OUTPUT: A trace monoid.

    EXAMPLES::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: F = TraceMonoid(names=('a', 'b', 'c'), I=Set({('a','c'), ('c','a')})); F
    Trace monoid on 3 generators ([a], [b], [c]) over Free monoid on 3 generators (a, b, c) with independence relation {(a, c), (c, a)}
    sage: x = F.gens()
    sage: x[0]*x[1]**5 * (x[0]*x[2])
    [a*b^5*a*c]

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
    sage: latex(M)
    \langle a, b, c \mid ac=ca \rangle

    TESTS::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
        sage: M.number_of_words(3) == len(M.words(3))
        True
    """

    Element = TraceMonoidElement

    def __init__(self, M=None, I=None, names=None):
        r"""
        Initializes TraceMonoid.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: TestSuite(M).run()
        """
        if M:
            pass
        elif names:
            M = FreeMonoid(len(names), names=names)
        else:
            raise ValueError()

        Monoid_class.__init__(self, names=names or [str(g) for g in M.gens()])
        self._base = M

        if not I:
            raise ValueError()
        if not isinstance(I, Set_generic):
            I = Set(I)

        base_gens = set(self._base.gens())
        el = next(iter(I))[0]
        if isinstance(el, str):
            f = self._base.monoid_generators()
            reversed_family = {str(f[k]): k for k in f.keys()}
            I = Set([(self._base.gen(reversed_family[e1]), self._base.gen(reversed_family[e2])) for e1, e2 in I])
        else:
            for (e1, e2) in I:
                if e1 not in base_gens or e2 not in base_gens:
                    raise ValueError()

        self._independence = I

    def ngens(self):
        """
        The number of generators of the monoid.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.ngens()
            3
        """
        return self._base.ngens()

    def one(self):
        """
        Neutral element of the monoid.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.one()
            1
        """
        return self(1)

    def gen(self, i=0):
        """
        The `i`-th generator of the monoid.

        INPUT:

        - ``i`` -- integer (default: 0)

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.gen(1)
            [b]
            sage: M.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        return self.element_class(self, self._base.gen(i))

    def cardinality(self):
        return self._base.cardinality()

    def _compute_dependence_stack(self, x):
        r"""
        Return generator stacks formed from trace
        subelements with respect to non-commutativity.

        OUTPUT: used generators and list of stacks as tuple

        ALGORITHM:

        Let `x` be a word of monoid; we scan `x` from right to left;
        when processing a letter `a` it is pushed on its stack and a
        marker is pushed on the stack of all the letters `b` ( `b \neq a` )
        which do not commute with `a`.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = b*a*d*a*c*b
            sage: M._compute_dependence_stack(x)
            ({a, b, c, d},
             OrderedDict([(a, [False, False, True, True, False]), (b, [True, False, False, False, True]), (c, [True, False, False, False]), (d, [False, False, True, False])]))
        """
        independence = self._independence
        generators_set = set(e for e, _ in x)
        stacks = OrderedDict(sorted((g, []) for g in generators_set))
        for generator, times in reversed(list(x)):
            stacks[generator].extend(repeat(True, times))
            for other_gen in generators_set:
                if other_gen == generator:
                    continue
                if (generator, other_gen) not in independence:
                    stacks[other_gen].extend(repeat(False, times))
        return generators_set, stacks

    @cached_method
    def _compute_lex_normal_form(self, x):
        r"""
        Return lexicographic normal form of the free monoid
        element in free monoid terms.

        OUTPUT: trace monoid element

        ALGORITHM:

        Take among the letters being on the top of some stack that
        letter `a` being minimal with respect to the given lexicographic
        ordering. We pop a marker from each stack corresponding to a
        letter `b` ( `b \neq a` ) which does not commute with `a`. We repeat
        this loop until all stacks are empty.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: M._compute_lex_normal_form(c*a*c*b*a^2)
            c*a*b*c*a^2
        """
        if not x._element_list:
            return x
        generators_set, stacks = self._compute_dependence_stack(x)
        independence = self._independence

        elements = []
        while any(stacks.values()):
            for generator, g_stack in stacks.items():
                if g_stack and g_stack[-1]:
                    g_stack.pop()
                    elements.append(generator)
                    for other_gen in generators_set:
                        if other_gen != generator and \
                                (generator, other_gen) not in independence:
                            stacks[other_gen].pop()
                    break

        return reduce(operator.mul, elements)

    @cached_method
    def _compute_foata_normal_form(self, x):
        r"""
        Return Foata normal form of the monoid element.

        OUTPUT: tuple of steps

        ALGORITHM:

        Within a loop we form the set using letters being
        on the top of stacks; arranging the letters in the lexicographic
        order yields a step of the Foata normal form;
        This loop is repeated until all stacks are empty.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = b*a*d*a*c*b
            sage: M._compute_foata_normal_form(x)
            (b, a*d, a, b*c)
            sage: y = b*a*a*d*b*a*b*c^2*a
            sage: M._compute_foata_normal_form(y)
            (b, a*d, a, b, a, b*c, c, a)
        """
        if not x._element_list:
            return tuple()

        generators_set, stacks = self._compute_dependence_stack(x)
        independence = self._independence

        steps = []
        while any(stacks.values()):
            step = []
            for generator, g_stack in stacks.items():
                if g_stack and g_stack[-1]:
                    g_stack.pop()
                    step.append(generator)

            for g in step:
                for other_gen in generators_set:
                    if other_gen != g and (g, other_gen) not in independence:
                        stacks[other_gen].pop()

            steps.append(step)

        return tuple(reduce(operator.mul, step) for step in steps)

    def _element_constructor_(self, x):
        """
        Return ``x`` coerced into this trace monoid.

        One can create a free monoid element from the integer 1,
        free monoid elements of the same generators as internal one,
        and coerce everything that can coerce free monoid.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c,d> = FreeMonoid()
            sage: I = ((a,d), (d,a), (b,c), (c,b))
            sage: M.<ac,bc,cc,dc> = TraceMonoid(F, I=I)
            sage: x = b*a*d*a*c*b
            sage: M(x)
            [b*a^2*d*b*c]
        """
        if isinstance(x, TraceMonoidElement) and x.parent() is self:
            return x
        if isinstance(x, TraceMonoidElement) and x.parent() == self:
            return self.element_class(self, x._base)

        x = self._compute_lex_normal_form(self._base(x))
        return self.element_class(self, x)

    @cached_method
    def independence(self):
        r"""
        Return independence relation over the monoid.

        OUTPUT: set of commuting generator pairs.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c> = FreeMonoid()
            sage: I=Set(((a,c), (c,a)))
            sage: M.<ac,bc,cc> = TraceMonoid(F, I=I)
            sage: M.independence() == {(a,c), (c,a)}
            True
        """
        return Set(self._independence)

    @cached_method
    def dependence(self):
        r"""
        Return dependence relation over the monoid.

        OUTPUT: set of non-commuting generator pairs.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.dependence()
            {(a, a), (b, b), (b, a), (a, b), (c, b), (b, c), (c, c)}
        """
        return Set(pair for pair in product(self._base.gens(), repeat=2)
                   if pair not in self._independence)

    @cached_method
    def dependence_graph(self):
        r"""
        Return graph of dependence relation.

        OUTPUT: dependence graph with generators as vertices

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c> = FreeMonoid()
            sage: M.<ai,bi,ci> = TraceMonoid(F, I=((a,c), (c,a)))
            sage: M.dependence_graph() == Graph({a:[a,b], b:[b], c:[c,b]})
            True
        """
        return Graph(list(set(frozenset((e1, e2)) if e1 != e2 else (e1, e2)
                              for e1, e2 in self.dependence())), loops=True)

    @cached_method
    def independence_graph(self):
        r"""
        Return graph of dependence relation.

        OUTPUT: dependence graph with generators as vertices

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: F.<a,b,c> = FreeMonoid()
            sage: M.<ai,bi,ci> = TraceMonoid(F, I=((a,c), (c,a)))
            sage: M.independence_graph() == Graph({a:[c], b:[], c:[]})
            True
        """
        g = Graph(list(set(map(frozenset, self._independence))))
        g.add_vertices(self._base.gens())
        return g

    @cached_method
    def dependence_polynomial(self, t=None):
        r"""
        Return dependence polynomial.

        The polynomial is defined as follows: `\\sum{i}{(-1)^i c_i t^i}`,
        where `c_i` equals to number of full subgraphs
        of size `i` in the independence graph.

        OUTPUT: polynomial over integer ring

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); M
            Trace monoid on 4 generators ([a], [b], [c], [d]) over Free monoid on 4 generators (a, b, c, d) with independence relation {(a, d), (b, c), (c, b), (d, a)}
            sage: M.dependence_polynomial()
            1/(2*t^2 - 4*t + 1)
        """
        if t is None:
            R = PolynomialRing(ZZ, 't')
            t = R.gen()
        clique_seq = self.independence_graph().clique_polynomial().coefficients()
        return Integer(1) / sum(((-1) ** i) * coeff * (t ** i)
                                for i, coeff in enumerate(clique_seq))

    @cached_method
    def number_of_words(self, length):
        r"""
        Return number of unique words of defined length.

        INPUT:

        - ``length`` -- integer; defines size of words what number should be computed

        OUTPUT: words number as integer

        EXAMPLES:

        Get number of words of size 3 ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); M
            Trace monoid on 4 generators ([a], [b], [c], [d]) over Free monoid on 4 generators (a, b, c, d) with independence relation {(a, d), (b, c), (c, b), (d, a)}
            sage: M.number_of_words(3)
            48
        """
        psr = PowerSeriesRing(ZZ, default_prec=length + 1)
        return psr(self.dependence_polynomial()).coefficients()[length]

    @cached_method
    def words(self, length):
        r"""
        Return all lexicographic forms of defined length.

        INPUT:

        - ``length`` -- integer; defines size of words

        OUTPUT: set of traces of size ``length``

        EXAMPLES:

        All words of size 2 ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: M.words(2)
            {[a^2], [a*b], [a*c], [a*d], [b*a], [b^2], [b*c], [b*d], [c*a], [c^2], [c*d], [d*b], [d*c], [d^2]}

        Get number of words of size 3 ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: len(M.words(3))
            48

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','b'), ('b','a'), ('b', 'c'), ('c', 'b')))
            sage: for i in range(10):
            ....:    assert len(M.words(i)) == M.number_of_words(i)
            sage: True
            True
        """
        if length < 0:
            raise ValueError("Bad length of words. Expected zero or positive number.")
        if length == 0:
            return {self(1)}
        if length == 1:
            return set(self.gens())

        return {
            word * suffix for word in self.words(length - 1)
            for suffix in self.gens()
            if not ((list(word._repr)[-1][0], suffix._repr) in self._independence
                    and list(word._repr)[-1][0] > suffix._repr)
        }

    def _repr_(self):
        r"""
        Textual representation of trace monoids.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); M
            Trace monoid on 4 generators ([a], [b], [c], [d]) over Free monoid on 4 generators (a, b, c, d) with independence relation {(a, d), (b, c), (c, b), (d, a)}
        """
        return "Trace monoid on {!s} generators {!s} " \
               "over {} with independence relation {{{}}}" \
            .format(self.ngens(), self.gens(), self._base, repr(sorted(self.independence()))[1:-1])

    def _latex_(self):
        r"""
        LaTeX representation of trace monoids.

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); latex(M)
            \langle a, b, c, d \mid ad=da,bc=cb \rangle
        """
        return "\\langle {} \\mid {} \\rangle".format(
            repr(self._base.gens())[1:-1],
            ",".join(
                "{0!r}{1!r}={1!r}{0!r}".format(v1, v2)
                for v1, v2 in sorted(set(map(frozenset, self._independence)))
            )
        )
