r"""
Module of trace monoid related structures.

Contains TraceMonoid, TraceMonoidElement and FoataForm classes.

EXAMPLES::

The following example demonstrates a monoid creation ::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a'))); M
    Trace monoid on 3 generators (a, b, c) over independence relation {(c, a), (a, c)}

Different monoid elements can be equal because of partially commutative multiplication ::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
    sage: c*a*b == a*c*b
    True

Lets ensure that it is a monoid ::

    sage: from sage.monoids.trace_monoid import TraceMonoid
    sage: TraceMonoid() in Monoids()
    True

AUTHORS::

- Pavlo Tokariev (2019-05-26): initial version

"""

# ******************************************************************************
#  Copyright (C) 2019      Pavlo Tokariev (KhNU) <pavlo.tokariev@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************

from __future__ import print_function

import operator
from collections import OrderedDict
from itertools import repeat, chain, product, combinations_with_replacement

from sage.misc.cachefunc import cached_method

from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.sets.set import Set
from sage.structure.parent import Set_generic


class TraceMonoidElement(FreeMonoidElement):
    """
    Element of a trace monoid.

    EXAMPLES::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
        sage: M.<a,b,c,d> = TraceMonoid(I=I)
        sage: x = b*a*d*a*c*b
        sage: x**3
        b*a*d*a*c*b^2*a*d*a*c*b^2*a*d*a*c*b
        sage: x**0
        1
        sage: x.lexic_norm_form()
        b*a^2*d*b*c
        sage: x.foata_norm_form()
        (b)(a*d)(a)(b*c)

    """

    def _dependence_stack(self):
        r"""
        Return generator stacks formed from trace
        subelements with respect to non-commutativity.

        OUTPUT: used generators and list of stacks as tuple

        ALGORITHM:

        Let `x` be a word of monoid; we scan `x` from right to left;
        when processing a letter `a` it is pushed on its stack and a
        marker is pushed on the stack of all the letters `b` ( `b \neq a` )
        which do not commute with `a`.
        """
        independence = self.parent()._independence
        generators_set = sorted(e[0] for e in self._element_list)
        stacks = OrderedDict(sorted((g, []) for g in generators_set))
        for generator, amount in reversed(self._element_list):
            stacks[generator].extend(repeat(True, amount))
            for other_gen in generators_set:
                if other_gen == generator:
                    continue
                if (generator, other_gen) not in independence:
                    stacks[other_gen].extend(repeat(False, amount))
        return generators_set, stacks

    @cached_method
    def lexic_norm_form(self):
        r"""
        Return lexicographic normal form of the monoid element.

        OUTPUT: trace monoid element

        ALGORITHM:

        Take among the letters being on the top of some stack that
        letter `a` being minimal with respect to the given lexicographic
        ordering. We pop a marker from each stack corresponding to a
        letter `b` ( `b \neq a` ) which does not commute with `a`. We repeat
        this loop until all stacks are empty.

        EXAMPLES::

            Get the normal form ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: (c*a*c*b*a^2).lexic_norm_form()
            a*c^2*b*a^2
        """
        monoid = self.parent()
        if not self._element_list:
            return self
        generators_set, stacks = self._dependence_stack()
        independence = monoid._independence

        elements = []
        while True:
            empty_stacks = []
            for generator, g_stack in stacks.items():
                if g_stack:
                    empty_stacks.append(False)
                    if g_stack[-1]:
                        g_stack.pop()
                        elements.append(generator)
                        for other_gen in generators_set:
                            if other_gen != generator and \
                                    (generator, other_gen) not in independence:
                                stacks[other_gen].pop()
                        break
                else:
                    empty_stacks.append(True)

            if all(empty_stacks):
                break

        return monoid([(e, 1) for e in elements])

    @cached_method
    def foata_norm_form(self):
        r"""
        Return Foata normal form of the monoid element.

        OUTPUT: tuple of steps

        ALGORITHM:

        Within a loop we form the set using letters being
        on the top of stacks; arranging the letters in the lexicographic
        order yields a step of the Foata normal form;
        This loop is repeated until all stacks are empty.

        EXAMPLES::

            Get the normal form ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); M
            Trace monoid on 4 generators (a, b, c, d) over independence relation {(a, d), (b, c), (c, b), (d, a)}
            sage: x = b*a*d*a*c*b
            sage: x.foata_norm_form()
            (b)(a*d)(a)(b*c)
        """
        monoid = self.parent()
        if not self._element_list:
            return tuple()

        generators_set, stacks = self._dependence_stack()
        independence = monoid._independence

        steps = []
        while True:
            empty_stacks = []
            step = []
            for generator, g_stack in stacks.items():
                if g_stack:
                    if g_stack[-1]:
                        g_stack.pop()
                        step.append(generator)
                    empty_stacks.append(False)
                else:
                    empty_stacks.append(True)

            if all(empty_stacks):
                break

            for g in step:
                for other_gen in generators_set:
                    if other_gen != g and (g, other_gen) not in independence:
                        stacks[other_gen].pop()

            if step:
                steps.append(step)

        return tuple(monoid(list((v, 1) for v in step)) for step in steps)

    def _plain_elements(self):
        repeated_el_iter = (repeat(e, times) for e, times in self._element_list)
        return list(chain.from_iterable(repeated_el_iter))

    @cached_method
    def dependence_graph(self):
        r"""
        Return dependence graph of the trace.

        It is a directed graph where all dependent (non-commutative)
        generators are connected by edges which
        direction depend on the generator position in the trace.

        OUTPUT: directed graph of generator indexes
        """
        elements = self._plain_elements()
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
          to compute Hasse diagram; there are two variants: `naive` and `min`.

        OUTPUT: directed graph of generator indexes
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
        """
        elements = self._plain_elements()
        elements.reverse()
        independence = self.parent()._independence
        min = set()
        graph = DiGraph({})

        def next_front(front, removed):
            removed = set(dest for _, dest in graph.outgoing_edges(removed, labels=False))
            front = set(dest for _, dest in graph.outgoing_edges(front, labels=False))
            return front - removed

        for i, x in enumerate(elements):
            front = min.copy()
            while front:
                removed = set()
                for j in list(front):
                    y = elements[j]
                    if (x, y) not in independence:
                        if j in min:
                            min.remove(j)
                        removed.add(j)
                        graph.add_edge(i, j)
                front = next_front(front, removed)

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

        Transform each trace to lexicographic form and then compare them.
        Equality means that they have the same "meaning"
        in terms of trace theory

        OUTPUT: 0 if self == other, -1 if self < other, 1 if self > other
        """
        other = other.lexic_norm_form()
        return super(TraceMonoidElement, self.lexic_norm_form())._richcmp_(other, op)

    def alphabet(self):
        r"""
        Return alphabet of a trace.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.alphabet()
            {b, a, d, c}

        OUTPUT: set of generators
        """
        return Set(self.to_list())

    def projection(self, letters):
        r"""
        Return a trace that formed from `self` using filtering by `letters`

        INPUT:

        - ``letters`` -- set of generators; defines set of letters that will be
          used to filter the trace.

        EXAMPLES::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: x = b*a*d*a*c*b
            sage: x.projection({a,b})
            b*a^2*b
            sage: x.projection({b,d,c})
            b*d*c*b

        OUTPUT: a trace
        """
        return reduce(operator.mul, [x for x in self.to_list() if x in letters])


class TraceMonoid(FreeMonoid):
    r"""
    Return a free partially commuting monoid (trace monoid) on `n` generators
    over independence relation `I`.

    We construct a trace monoid by specifing:

    - the number of generators and a relation
    - the names of the generators and a relation

    INPUT:

    - ``n`` -- number of generators

    -  ``names`` -- names of generators

    - ``I`` -- commutation relation between generators
      (or their names if the ``names`` are given)

    OUTPUT:

    A trace monoid.

    EXAMPLES::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: F = TraceMonoid(names=('a', 'b', 'c'), I=Set({('a','c'), ('c','a')})); F
        Trace monoid on 3 generators (a, b, c) over independence relation {(c, a), (a, c)}
        sage: x = F.gens()
        sage: x[0]*x[1]**5 * (x[0]*x[2])
        a*b^5*a*c
        sage: F = TraceMonoid(3, 'a')
        sage: F
        Trace monoid on 3 generators (a0, a1, a2) over independence relation {}

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a'))); M
        Trace monoid on 3 generators (a, b, c) over independence relation {(c, a), (a, c)}
        sage: latex(M)
        \langle a, b, c \mid ac=ca \rangle

    TESTS::

        sage: from sage.monoids.trace_monoid import TraceMonoid
        sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
        sage: M.number_of_words(3) == len(M.words(3))
        True
    """

    Element = TraceMonoidElement

    def __init__(self, n=None, names=None, I=None):
        if n and not names:
            names = 'x'
        if names and not n:
            n = len(names)
        if not n and not names:
            n = 0
        elif n and names:
            if isinstance(names, str):
                names = [names + str(i) for i in range(n)]
            elif len(names) != n:
                raise ValueError(
                    "Argument `n` should be equal to size "
                    "of `names`, got n={}, names={}.".format(n, names))

        FreeMonoid.__init__(self, n, names)

        if I is None:
            I = Set()
        elif n and len(I) > 0:
            el = next(iter(I))[0]
            if isinstance(el, str):
                f = self.monoid_generators()
                reversed_family = {str(f[k]): k for k in f.keys()}
                I = ((reversed_family[e1], reversed_family[e2]) for e1, e2 in I)
        elif not n and I:
            raise ValueError(
                "The monoid contains zero elements, no relation is allowed.")
        if not isinstance(I, Set_generic):
            I = Set(I)

        self._independence = I

    @cached_method
    def independence(self):
        r"""
        Return independence relation over the monoid.

        OUTPUT: set of commuting generator pairs.

        EXAMPLES::

            Print the relation ::

                sage: from sage.monoids.trace_monoid import TraceMonoid
                sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
                sage: sorted(M.independence)
                [(a, c), (c, a)]

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.independence() == I
            True
        """
        f = self.monoid_generators()
        return Set((f[v1], f[v2]) for v1, v2 in self._independence)

    @cached_method
    def _named_set_without_duplicates(self):
        r"""
        Return set of commuting generators without pair permutations.

        OUTPUT: list of generator pairs

        TESTS::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M._named_set_without_duplicates()
            [(a, c)]
        """
        f = self.monoid_generators()
        return [(f[v1], f[v2]) for v1, v2 in set(map(frozenset, self._independence))]

    @cached_method
    def _dependence(self):
        r"""
        Return set of commuting generator indexes without pair permutations.

        OUTPUT: list of generator index pairs
        """
        return Set(
            pair for pair in product(range(self.ngens()), repeat=2)
            if pair not in self._independence
        )

    @cached_method
    def dependence(self):
        r"""
        Return dependence relation over the monoid.

        OUTPUT: set of non-commuting generator pairs.

        EXAMPLES:

        Print the relation ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: M.<a,b,c> = TraceMonoid(I=(('a','c'), ('c','a')))
            sage: M.dependence
            {(c, c), (b, b), (b, a), (a, b), (c, b), (b, c), (a, a)}

        """
        f = self.monoid_generators()
        return Set((f[v1], f[v2]) for v1, v2 in self._dependence())

    @cached_method
    def dependence_graph(self):
        r"""
        Return graph of dependence relation.

        OUTPUT: dependence graph with generators as vertices
        """
        f = self.monoid_generators()
        return Graph([
            (f[v1], f[v2]) for v1, v2 in combinations_with_replacement(range(self.ngens()), 2)
            if (v1, v2) not in self._independence
        ], loops=True)

    @cached_method
    def independence_graph(self):
        r"""
        Return graph of dependence relation.

        OUTPUT: dependence graph with generators as vertices
        """
        g = Graph(self._named_set_without_duplicates())
        g.add_vertices(self.gens())
        return g

    @cached_method
    def dependence_polynomial(self, t=None):
        r"""
        Return dependence polynomial.

        The polynomial is defined as follows: `\\sum{i}{(-1)^i c_i t^i}`,
        where `c_i` equals to number of full subgraphs
        of size `i` in the independence graph.

        OUTPUT: polynomial over integer ring

        EXAMPLES:

        Get the polynomial ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I); M
            Trace monoid on 4 generators (a, b, c, d) over independence relation {(a, d), (b, c), (c, b), (d, a)}
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
            Trace monoid on 4 generators (a, b, c, d) over independence relation {(a, d), (b, c), (c, b), (d, a)}
            sage: M.number_of_words(3)
            48
            """
        psr = PowerSeriesRing(ZZ, default_prec=length + 1)
        return psr(self.dependence_polynomial()).coefficients()[length]

    @cached_method
    def words(self, length):
        r"""
        Return all words in lexicographic form of defined length.

        INPUT:

        - ``length`` -- integer; defines size of words
            what number should be computed

        OUTPUT: words number as integer

        EXAMPLES:

        Get number of words of size 3 ::

            sage: from sage.monoids.trace_monoid import TraceMonoid
            sage: I = (('a','d'), ('d','a'), ('b','c'), ('c','b'))
            sage: M.<a,b,c,d> = TraceMonoid(I=I)
            sage: len(M.words(3))
            48
            """
        if length < 0:
            raise ValueError("Bad length of words. Expected zero or positive number.")
        if length == 0:
            return [self(1)]
        if length == 1:
            return list(self.gens())

        return [
            word * self.gen(suffix) for word in self.words(length - 1)
            for suffix in range(self.ngens())
            if not ((word._element_list[-1][0], suffix) in self._independence
                    and word._element_list[-1][0] > suffix)
        ]

    def _repr_(self):
        return "Trace monoid on {!s} generators {!s} " \
               "over independence relation {!r}" \
            .format(self.ngens(), self.gens(), self.independence())

    def _latex_(self):
        return "\\langle {} \\mid {} \\rangle".format(
            repr(self.gens())[1:-1],
            ",".join(
                "{0}{1}={1}{0}".format(v1, v2)
                for v1, v2 in self._named_set_without_duplicates()
            )
        )
