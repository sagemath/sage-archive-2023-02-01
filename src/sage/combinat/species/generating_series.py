r"""
Generating Series

This file makes a number of extensions ot lazy power series by
endowing them with some semantic content for how they're to be
interpreted.

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevent section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatse10.html.
One notable difference is that we use power-sum symmetric functions
as the coefficients of our cycle index series.

TESTS::

    sage: from sage.combinat.species.stream import Stream, _integers_from
    sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
    sage: p = SFAPower(QQ)
    sage: CIS = CycleIndexSeriesRing(QQ)

::

    sage: geo1 = CIS((p([1])^i  for i in _integers_from(0)))
    sage: geo2 = CIS((p([2])^i  for i in _integers_from(0)))
    sage: s = geo1 * geo2
    sage: s[0]
    p[]
    sage: s[1]
    p[1] + p[2]
    sage: s[2]
    p[1, 1] + p[2, 1] + p[2, 2]
    sage: s[3]
    p[1, 1, 1] + p[2, 1, 1] + p[2, 2, 1] + p[2, 2, 2]

Whereas the coefficients of the above test are homogeneous with
respect to total degree, the following test groups with respect to
weighted degree where each variable x_i has weight i.

::

    sage: def g():
    ...       for i in _integers_from(0):
    ...           yield p([2])^i
    ...           yield p(0)
    sage: geo1 = CIS((p([1])^i  for i in _integers_from(0)))
    sage: geo2 = CIS(g())
    sage: s = geo1 * geo2
    sage: s[0]
    p[]
    sage: s[1]
    p[1]
    sage: s[2]
    p[1, 1] + p[2]
    sage: s[3]
    p[1, 1, 1] + p[2, 1]
    sage: s[4]
    p[1, 1, 1, 1] + p[2, 1, 1] + p[2, 2]
"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from series import LazyPowerSeriesRing, LazyPowerSeries
from stream import Stream, _integers_from
from sage.rings.all import Integer, moebius, lcm, divisors
from sage.combinat.partition import Partition, Partitions
from functools import partial
from sage.combinat.sf.all import SFAPower
from sage.misc.cachefunc import cached_method, cached_function

@cached_function
def OrdinaryGeneratingSeriesRing(R):
    """
    Returns the ring of ordinary generating series. Note that is is
    just a LazyPowerSeriesRing whose elements have some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
        sage: R = OrdinaryGeneratingSeriesRing(QQ); R
        Lazy Power Series Ring over Rational Field
        sage: R([1]).coefficients(4)
        [1, 1, 1, 1]
        sage: R([1]).counts(4)
        [1, 1, 1, 1]

    TESTS: We test to make sure that caching works.

    ::

        sage: R is OrdinaryGeneratingSeriesRing(QQ)
        True
    """
    return OrdinaryGeneratingSeriesRing_class(R)

class OrdinaryGeneratingSeriesRing_class(LazyPowerSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: R == loads(dumps(R))
            True
        """
        LazyPowerSeriesRing.__init__(self, R, OrdinaryGeneratingSeries)

class OrdinaryGeneratingSeries(LazyPowerSeries):
    def count(self, n):
        """
        Returns the number of structures on a set of size n.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.count(10)
            10
        """
        return self.coefficient(n)

    def counts(self, n):
        """
        Returns the number of structures on a set for size i for i in
        range(n).

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.counts(10)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return [self.count(i) for i in range(n)]


@cached_function
def ExponentialGeneratingSeriesRing(R):
    """
    Returns the ring of ordinary generating series. Note that is is
    just a LazyPowerSeriesRing whose elements have some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
        sage: R = ExponentialGeneratingSeriesRing(QQ); R
        Lazy Power Series Ring over Rational Field
        sage: R([1]).coefficients(4)
        [1, 1, 1, 1]
        sage: R([1]).counts(4)
        [1, 1, 2, 6]

    TESTS: We test to make sure that caching works.

    ::

        sage: R is ExponentialGeneratingSeriesRing(QQ)
        True
    """
    return ExponentialGeneratingSeriesRing_class(R)

class ExponentialGeneratingSeriesRing_class(LazyPowerSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: R == loads(dumps(R))
            True
        """
        LazyPowerSeriesRing.__init__(self, R, ExponentialGeneratingSeries)

class ExponentialGeneratingSeries(LazyPowerSeries):
    def count(self, n):
        """
        Returns the number of structures of size n.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: f = R([1])
            sage: [f.count(i) for i in range(7)]
            [1, 1, 2, 6, 24, 120, 720]
        """
        return factorial_stream[n] * self.coefficient(n)

    def counts(self, n):
        """
        Returns the number of structures on a set for size i for i in
        range(n).

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.counts(5)
            [0, 1, 4, 18, 96]
        """
        return [self.count(i) for i in range(n)]

    def functorial_composition(self, y):
        r"""
        Returns the exponential generating series which is the functorial
        composition of self with y.

        If `f = \sum_{n=0}^{\infty} f_n \frac{x^n}{n!}` and
        `g = \sum_{n=0}^{\infty} f_n \frac{x^n}{n!}`, then
        functorial composition `f \Box g` is defined as


        .. math::

             f \Box g = \sum_{n=0}^{\infty} f_{g_n} \frac{x^n}{n!}



        REFERENCES:

        - Section 2.2 of BLL.

        EXAMPLES::

            sage: G = species.SimpleGraphSpecies()
            sage: g = G.generating_series()
            sage: g.coefficients(10)
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]
        """
        return self._new(partial(self._functorial_compose_gen, y), lambda a,b: 0, self, y)

    def _functorial_compose_gen(self, y, ao):
        """
        Returns a generator for the coefficients of the functorial
        composition of self with y.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: E2 = E.restricted(min=2, max=3)
            sage: WP = species.SubsetSpecies()
            sage: P2 = E2*E
            sage: g1 = WP.generating_series()
            sage: g2 = P2.generating_series()
            sage: g = g1._functorial_compose_gen(g2, 0)
            sage: [g.next() for i in range(10)]
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]
        """
        n = 0
        while True:
            yield self.count(y.count(n))/factorial_stream[n]
            n += 1

def factorial_gen():
    """
    A generator for the factorials starting at 0.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import factorial_gen
        sage: g = factorial_gen()
        sage: [g.next() for i in range(5)]
        [1, 1, 2, 6, 24]
    """
    z = Integer(1)
    yield z
    yield z
    n = Integer(2)
    while True:
        z *= n
        yield z
        n += 1

factorial_stream = Stream(factorial_gen())



@cached_function
def CycleIndexSeriesRing(R):
    """
    Returns the ring of cycle index series. Note that is is just a
    LazyPowerSeriesRing whose elements have some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
        sage: R = CycleIndexSeriesRing(QQ); R
        Cycle Index Series Ring over Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
        sage: R([1]).coefficients(4)
        [1, 1, 1, 1]

    TESTS: We test to make sure that caching works.

    ::

        sage: R is CycleIndexSeriesRing(QQ)
        True
    """
    return CycleIndexSeriesRing_class(R)

class CycleIndexSeriesRing_class(LazyPowerSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: R = CycleIndexSeriesRing(QQ); R
            Cycle Index Series Ring over Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
            sage: R == loads(dumps(R))
            True
        """
        R = SFAPower(R)
        LazyPowerSeriesRing.__init__(self, R, CycleIndexSeries)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: CycleIndexSeriesRing(QQ)
            Cycle Index Series Ring over Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
        """
        return "Cycle Index Series Ring over %s"%self.base_ring()


class CycleIndexSeries(LazyPowerSeries):
    def count(self, t):
        """
        Returns the number of structures corresponding to a certain cycle
        type.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SFAPower(QQ)
            sage: CIS = CycleIndexSeriesRing(p)
            sage: f = CIS([0, p([1]), 2*p([1,1]),3*p([2,1])])
            sage: f.count([1])
            1
            sage: f.count([1,1])
            4
            sage: f.count([2,1])
            6
        """
        t = Partition(t)
        return t.aut() * self.coefficient_cycle_type(t)

    def coefficient_cycle_type(self, t):
        """
        Returns the coefficient of a cycle type t.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SFAPower(QQ)
            sage: CIS = CycleIndexSeriesRing(p)
            sage: f = CIS([0, p([1]), 2*p([1,1]),3*p([2,1])])
            sage: f.coefficient_cycle_type([1])
            1
            sage: f.coefficient_cycle_type([1,1])
            2
            sage: f.coefficient_cycle_type([2,1])
            3
        """
        t = Partition(t)
        p = self.coefficient(t.size())
        return p.coefficient(t)


    def stretch(self, k):
        r"""
        Returns the stretch of a cycle index series by a positive integer
        `k`.

        If

        .. math::

           f = \sum_{n=0}^{\infty} f_n(x_1, x_2, \ldots, x_n),

        then the stretch `g` of `f` by `k` is

        .. math::

           g = \sum_{n=0}^{\infty} f_n(x_k, x_{2k}, \ldots, x_{nk}).

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SFAPower(QQ)
            sage: CIS = CycleIndexSeriesRing(p)
            sage: f = CIS([p([1])])
            sage: f.stretch(3).coefficients(10)
            [p[3], 0, 0, p[3], 0, 0, p[3], 0, 0, p[3]]
        """
        return self._new(partial(self._stretch_gen, k), lambda ao: k*ao, self)

    def _stretch_gen(self, k, ao):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SFAPower(QQ)
            sage: CIS = CycleIndexSeriesRing(p)
            sage: f = CIS([p([1])])
            sage: g = f._stretch_gen(2,0)
            sage: [g.next() for i in range(10)]
            [p[2], 0, p[2], 0, p[2], 0, p[2], 0, p[2], 0]
        """
        from sage.combinat.partition import Partition
        BR = self.base_ring()
        zero = BR(0)

        stretch_k = lambda p: Partition([k*i for i in p])

        yield self.coefficient(0).map_support(stretch_k)

        n = 1
        while True:
            for i in range(k-1):
                yield BR(0)
            yield self.coefficient(n).map_support(stretch_k)
            n += 1

    def isotype_generating_series(self):
        """
        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()
            sage: f = cis.isotype_generating_series()
            sage: f.coefficients(10)
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        R = self.base_ring().base_ring()
        OGS = OrdinaryGeneratingSeriesRing(R)()
        return OGS._new(self._ogs_gen, lambda ao: ao, self)

    def _ogs_gen(self, ao):
        """
        Returns a generator for the coefficients of the ordinary generating
        series obtained from a cycle index series.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()
            sage: g = cis._ogs_gen(0)
            sage: [g.next() for i in range(10)]
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        for i in range(ao):
            yield 0
        for i in _integers_from(ao):
            yield sum( self.coefficient(i).coefficients() )

    def generating_series(self):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: cis = P.cycle_index_series()
            sage: f = cis.generating_series()
            sage: f.coefficients(5)
            [1, 1, 1, 5/6, 5/8]
        """
        R = self.base_ring().base_ring()
        EGS = ExponentialGeneratingSeriesRing(R)()
        return EGS._new(self._egs_gen, lambda ao: ao, self)

    def _egs_gen(self, ao):
        """
        Returns a generator for the coefficients of the exponential
        generating series obtained from a cycle index series.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()
            sage: g = cis._egs_gen(0)
            sage: [g.next() for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        for i in range(ao):
            yield 0
        for i in _integers_from(ao):
            yield self.coefficient(i).coefficient([1]*i)


    def functorial_composition(self, g):
        r"""
        Returns the functorial composition of self and g.

        If `F`, `G`, and `H` are combinatorial
        species such that ` H = F \Box G `, then the cycle index
        series `Z_H = Z_F \Box Z_G` is defined as

        EXAMPLES::

            sage: S = species.SimpleGraphSpecies()
            sage: S.cycle_index_series().coefficients(5)
            [p[],
             p[1],
             p[1, 1] + p[2],
             4/3*p[1, 1, 1] + 2*p[2, 1] + 2/3*p[3],
             8/3*p[1, 1, 1, 1] + 4*p[2, 1, 1] + 2*p[2, 2] + 4/3*p[3, 1] + p[4]]
        """
        return self._new(partial(self._functorial_compose_gen, g), lambda a,b: 0, self, g)

    def _functorial_compose_gen(self, g, ao):
        """
        Returns s generator for the coefficients of the functorial
        composition of self with g.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: E2 = species.SetSpecies(size=2)
            sage: WP = species.SubsetSpecies()
            sage: P2 = E2*E
            sage: P2_cis = P2.cycle_index_series()
            sage: WP_cis = WP.cycle_index_series()
            sage: g = WP_cis._functorial_compose_gen(P2_cis,0)
            sage: [g.next() for i in range(5)]
            [p[],
             p[1],
             p[1, 1] + p[2],
             4/3*p[1, 1, 1] + 2*p[2, 1] + 2/3*p[3],
             8/3*p[1, 1, 1, 1] + 4*p[2, 1, 1] + 2*p[2, 2] + 4/3*p[3, 1] + p[4]]
        """
        p = self.parent().base_ring()
        n = 0
        while True:
            res = p(0)
            for s in Partitions(n):
                t = g._cycle_type(s)
                q = self.count(t) / s.aut()
                res += q*p(s)
            yield res
            n += 1

    def _cycle_type(self, s):
        """
        EXAMPLES::

            sage: cis = species.PartitionSpecies().cycle_index_series()
            sage: [cis._cycle_type(p) for p in Partitions(3)]
            [[3, 1, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: cis = species.PermutationSpecies().cycle_index_series()
            sage: [cis._cycle_type(p) for p in Partitions(3)]
            [[3, 1, 1, 1], [2, 2, 1, 1], [1, 1, 1, 1, 1, 1]]
            sage: cis = species.SetSpecies().cycle_index_series()
            sage: [cis._cycle_type(p) for p in Partitions(3)]
            [[1], [1], [1]]
        """
        if s == []:
            return self._card(0)
        res = []
        for k in range(1, self._upper_bound_for_longest_cycle(s)+1):
            e = 0
            for d in divisors(k):
                m = moebius(d)
                if m == 0:
                    continue
                u = s.power(k/d)
                e += m*self.count(u)
            res.extend([k]*int(e/k))
        res.reverse()
        return Partition(res)


    def _upper_bound_for_longest_cycle(self, s):
        """
        EXAMPLES::

            sage: cis = species.PartitionSpecies().cycle_index_series()
            sage: cis._upper_bound_for_longest_cycle([4])
            4
            sage: cis._upper_bound_for_longest_cycle([3,1])
            3
            sage: cis._upper_bound_for_longest_cycle([2,2])
            2
            sage: cis._upper_bound_for_longest_cycle([2,1,1])
            2
            sage: cis._upper_bound_for_longest_cycle([1,1,1,1])
            1
        """
        if s == []:
            return 1
        return min(self._card(sum(s)), lcm(list(s)))

    def _card(self, n):
        """
        Returns the number of structures on an underlying set of size n for
        the species associated with self. This is just n! times the
        coefficient of p[1]n in self.

        EXAMPLES::

            sage: cis = species.PartitionSpecies().cycle_index_series()
            sage: cis._card(4)
            15
        """
        p = self.coefficient(n)
        return factorial_stream[n]*p.coefficient([1]*n)


    def _compose_gen(self, y, ao):
        """
        Returns a generator for the coefficients of the composition of this
        cycle index series and the cycle index series y. This overrides the
        the method defined in LazyPowerSeries.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: g = E_cis._compose_gen(C.cycle_index_series(),0)
            sage: [g.next() for i in range(4)]
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
        """
        assert y.coefficient(0) == 0
        y_powers = Stream(y._power_gen())

        parent = self.parent()
        res =  parent.sum_generator(self._compose_term(self.coefficient(i), y_powers)
                                    for i in _integers_from(0))

        for i in _integers_from(0):
            yield res.coefficient(i)

    def _compose_term(self, p, y_powers):
        """
        Returns the composition of one term in self with y.

        INPUT:


        -  ``p`` - a term in self

        -  ``y_powers`` - a stream for the powers of y
           starting with y


        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: C_cis = C.cycle_index_series()
            sage: c_powers = Stream(C_cis._power_gen())
            sage: p2 = E_cis.coefficient(2); p2
            1/2*p[1, 1] + 1/2*p[2]
            sage: E_cis._compose_term(p2, c_powers).coefficients(4)
            [0, 0, 1/2*p[1, 1] + 1/2*p[2], 1/2*p[1, 1, 1] + 1/2*p[2, 1]]
        """
        parent = self.parent()
        if p == 0:
            return parent(0)

        res = []
        #Go through all the partition, coefficient pairs in the term p
        for m, c in p:
            res_t = parent.term(c, 0)

            for e,v in enumerate(m.to_exp()):
                if v == 0:
                    continue
                res_t = res_t * y_powers[v-1].stretch(e+1)
            res.append(res_t)

        return parent.sum(res)

    def weighted_composition(self, y_species):
        """
        Returns the composition of this cycle index series with the cycle
        index series of y_species where y_species is a weighted species.
        Note that this is basically the same algorithm as composition
        except we can't use the optimization that the powering of cycle
        index series commutes with 'stretching'.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: E_cis.weighted_composition(C).coefficients(4)
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
            sage: E(C).cycle_index_series().coefficients(4)
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
        """
        base_ring = self.base_ring()
        y = y_species.cycle_index_series(base_ring)
        assert y.coefficient(0) == 0
        return self._new(partial(self._weighted_compose_gen, y_species), lambda a,b:a*b, self, y)


    def _weighted_compose_gen(self, y_species, ao):
        """
        Returns an iterator for the composition of this cycle index series
        and the cycle index series of the weighted species y_species.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: g = E_cis._weighted_compose_gen(C,0)
            sage: [g.next() for i in range(4)]
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
        """
        parent = self.parent()
        res =  parent.sum_generator(self._weighted_compose_term(self.coefficient(i), y_species)
                                    for i in _integers_from(0))

        for i in _integers_from(0):
            yield res.coefficient(i)

    def _weighted_compose_term(self, p, y_species):
        """
        Returns the weighted composition of one term in self with y.

        INPUT:


        -  ``p`` - a term in self

        -  ``y_species`` - a species


        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: p2 = E_cis.coefficient(2); p2
            1/2*p[1, 1] + 1/2*p[2]
            sage: E_cis._weighted_compose_term(p2, C).coefficients(4)
            [0, 0, 1/2*p[1, 1] + 1/2*p[2], 1/2*p[1, 1, 1] + 1/2*p[2, 1]]
        """
        parent = self.parent()
        if p == 0:
            return parent(0)

        base_ring = self.base_ring().base_ring()

        res = []
        #Go through all the partition, coefficient pairs in the term p
        for m, c in p:
            res_t = parent.term(c, 0)

            for e,v in enumerate(m.to_exp()):
                if v == 0:
                    continue
                res_t = res_t * (y_species.weighted(y_species._weight**(e+1)).cycle_index_series(base_ring)**v).stretch(e+1)
            res.append(res_t)

        return parent.sum(res)



