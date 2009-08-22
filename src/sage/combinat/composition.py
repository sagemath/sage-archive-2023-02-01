"""
Integer compositions

A composition `c` of a nonnegative integer `n` is a list of positive integers
(the *parts* of the compositions) with total sum `n`.

This module provides tools for manipulating compositions and enumerated
sets of compositions.

EXAMPLES::

    sage: Composition([5, 3, 1, 3])
    [5, 3, 1, 3]
    sage: list(Compositions(4))
    [[1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3], [2, 1, 1], [2, 2], [3, 1], [4]]


AUTHORS:

- Mike Hansen, Nicolas M. Thiery
- MuPAD-Combinat developers (algorithms and design inspiration)
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen       <mhansen@gmail.com>
#                     2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#              http://www.gnu.org/licenses/
#*****************************************************************************

import sage.combinat.skew_partition
from combinat import CombinatorialClass, CombinatorialObject, InfiniteAbstractCombinatorialClass
from cartesian_product import CartesianProduct
from integer_list import IntegerListsLex
import __builtin__
from sage.rings.integer import Integer
from sage.rings.arith import binomial
import misc

def Composition(co=None, descents=None, code=None):
    """
    Integer compositions

    A composition of a nonnegative integer `n` is a list
    `(i_1,\dots,i_k)` of positive integers with total sum `n`.

    TODO: mathematical definition of descents, code, ...

    EXAMPLES:

    The simplest way to create a composition is by specifying its
    entries as a list, tuple (or other iterable)::

        sage: Composition([3,1,2])
        [3, 1, 2]
        sage: Composition((3,1,2))
        [3, 1, 2]
        sage: Composition(i for i in range(2,5))
        [2, 3, 4]

    You can alternatively create a composition from the list of its
    descents::

        sage: Composition([1, 1, 3, 4, 3]).descents()
        [0, 1, 4, 8, 11]
        sage: Composition(descents=[1,0,4,8,11])
        [1, 1, 3, 4, 3]

    You can also create a composition from its code::

        sage: Composition([4,1,2,3,5]).to_code()
        [1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0]
        sage: Composition(code=_)
        [4, 1, 2, 3, 5]
        sage: Composition([3,1,2,3,5]).to_code()
        [1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0]
        sage: Composition(code=_)
        [3, 1, 2, 3, 5]
    """

    if descents is not None:
        if isinstance(descents, tuple):
            return from_descents(descents[0], nps=descents[1])
        else:
            return from_descents(descents)
    elif code is not None:
        return from_code(code)
    elif isinstance(co, Composition_class):
        return co
    else:
        co = list(co)
        assert(co in Compositions())
        #raise ValueError, "invalid composition"
        return Composition_class(co)


class Composition_class(CombinatorialObject):
    __doc__ = Composition.__doc__

    def conjugate(self):
        r"""
        Returns the conjugate of the composition comp.

        Algorithm from mupad-combinat.

        EXAMPLES::

            sage: Composition([1, 1, 3, 1, 2, 1, 3]).conjugate()
            [1, 1, 3, 3, 1, 3]
        """
        comp = self
        if comp == []:
            return Composition([])
        n = len(comp)
        coofcp = [sum(comp[:j])-j+1 for j in range(1,n+1)]

        cocjg = []
        for i in range(n-1):
            cocjg += [i+1 for _ in range(0, (coofcp[n-i-1]-coofcp[n-i-2]))]
        cocjg += [n for j in range(coofcp[0])]

        return Composition([cocjg[0]] + [cocjg[i]-cocjg[i-1]+1 for i in range(1,len(cocjg)) ])


    def complement(self):
        """
        Returns the complement of the composition co. The complement is the
        reverse of co's conjugate composition.

        EXAMPLES::

            sage: Composition([1, 1, 3, 1, 2, 1, 3]).conjugate()
            [1, 1, 3, 3, 1, 3]
            sage: Composition([1, 1, 3, 1, 2, 1, 3]).complement()
            [3, 1, 3, 3, 1, 1]
        """
        return Composition([element for element in reversed(self.conjugate())])


    def __add__(self, other):
        """
        Returns the concatenation of two compositions.

        EXAMPLES::

            sage: Composition([1, 1, 3]) + Composition([4, 1, 2])
            [1, 1, 3, 4, 1, 2]

        TESTS::

            sage: Composition([]) + Composition([]) == Composition([])
            True
        """
        return Composition(list(self)+list(other))

    def size(self):
        """
        Returns the size of the composition, that is the sum of its parts.

        EXAMPLES::

            sage: Composition([7,1,3]).size()
            11
        """
        return sum(self)

    @staticmethod
    def sum(compositions):
        """
        Returns the concatenation of the given compositions.

        INPUT:

        - ``compositions`` -- a list (or iterable) of compositions

        EXAMPLES::

            sage: sage.combinat.composition.Composition_class.sum([Composition([1, 1, 3]), Composition([4, 1, 2]), Composition([3,1])])
            [1, 1, 3, 4, 1, 2, 3, 1]

        Any iterable can be provided as input::

            sage: sage.combinat.composition.Composition_class.sum([Composition([i,i]) for i in [4,1,3]])
            [4, 4, 1, 1, 3, 3]

        Empty inputs are handled gracefully::

            sage: sage.combinat.composition.Composition_class.sum([]) == Composition([])
            True
        """
        return sum(compositions, Composition([]))

    def finer(self):
        """
        Returns the set of compositions which are finer than self.

        EXAMPLES::

            sage: C = Composition([3,2]).finer()
            sage: C.cardinality()
            8
            sage: list(C)
            [[1, 1, 1, 1, 1], [1, 1, 1, 2], [1, 2, 1, 1], [1, 2, 2], [2, 1, 1, 1], [2, 1, 2], [3, 1, 1], [3, 2]]
        """
        return CartesianProduct(*[Compositions(i) for i in self]).map(Composition_class.sum)

    def is_finer(self, co2):
        """
        Returns True if the composition self is finer than the composition
        co2; otherwise, it returns False.

        EXAMPLES::

            sage: Composition([4,1,2]).is_finer([3,1,3])
            False
            sage: Composition([3,1,3]).is_finer([4,1,2])
            False
            sage: Composition([1,2,2,1,1,2]).is_finer([5,1,3])
            True
            sage: Composition([2,2,2]).is_finer([4,2])
            True
        """
        co1 = self
        if sum(co1) != sum(co2):
            #Error: compositions are not of the same size
            raise ValueError, "compositions self (= %s) and co2 (= %s) must be of the same size"%(self, co2)


        sum1 = 0
        sum2 = 0
        i1 = 0
        for i2 in range(len(co2)):
            sum2 += co2[i2]
            while sum1 < sum2:
                sum1 += co1[i1]
                i1 += 1
            if sum1 > sum2:
                return False

        return True

    def fatten(self, grouping):
        """
        Returns the composition fatter than self, obtained by grouping
        together consecutive parts according to grouping.

        INPUT:

        - ``grouping`` -- a composition whose sum is the length of self

        EXAMPLES:

        Let us start with the composition::

            sage: c = Composition([4,5,2,7,1])

        With `grouping = (1,\dots,1)`, `c` is left unchanged::

            sage: c.fatten(Composition([1,1,1,1,1]))
            [4, 5, 2, 7, 1]

        With `grouping = (5)`, this yields the coarser composition above `c`::

            sage: c.fatten(Composition([5]))
            [19]

        Other values for ``grouping`` yield (all the) other compositions
        coarser to `c`::

            sage: c.fatten(Composition([2,1,2]))
            [9, 2, 8]
            sage: c.fatten(Composition([3,1,1]))
            [11, 7, 1]

        TESTS::

            sage: Composition([]).fatten(Composition([]))
            []
            sage: c.fatten(Composition([3,1,1])).__class__ == c.__class__
            True
        """
        result = [None] * len(grouping)
        j = 0
        for i in range(len(grouping)):
            result[i] = sum(self[j:j+grouping[i]])
            j += grouping[i]
        return Composition_class(result)

    def fatter(self):
        """
        Returns the set of compositions which are fatter than self.

        Complexity for generation: O(size(c)) memory, O(size(result)) time

        EXAMPLES::

            sage: C = Composition([4,5,2]).fatter()
            sage: C.cardinality()
            4
            sage: list(C)
            [[4, 5, 2], [4, 7], [9, 2], [11]]

        Some extreme cases::

            sage: list(Composition([5]).fatter())
            [[5]]
            sage: list(Composition([]).fatter())
            [[]]
            sage: list(Composition([1,1,1,1]).fatter()) == list(Compositions(4))
            True
        """
        return Compositions(len(self)).map(self.fatten)

    def refinement(self, co2):
        """
        Returns the refinement composition of self and co2.

        EXAMPLES::

            sage: Composition([1,2,2,1,1,2]).refinement([5,1,3])
            [3, 1, 2]
        """
        co1 = self
        if sum(co1) != sum(co2):
            #Error: compositions are not of the same size
            raise ValueError, "compositions self (= %s) and co2 (= %s) must be of the same size"%(self, co2)

        sum1 = 0
        sum2 = 0
        i1 = -1
        result = []
        for i2 in range(len(co2)):
            sum2 += co2[i2]
            i_res = 0
            while sum1 < sum2:
                i1 += 1
                sum1 += co1[i1]
                i_res += 1

            if sum1 > sum2:
                return None

            result.append(i_res)

        return Composition(result)

    def major_index(self):
        """
        Returns the major index of the composition co. The major index is
        defined as the sum of the descents.

        EXAMPLES::

            sage: Composition([1, 1, 3, 1, 2, 1, 3]).major_index()
            31
        """
        co = self
        lv = len(co)
        if lv == 1:
            return 0
        else:
            return sum([(lv-(i+1))*co[i] for i in range(lv)])


    def to_code(self):
        """
        Returns the code of the composition self. The code of a composition
        is a list of length self.size() of 1s and 0s such that there is a 1
        wherever a new part starts.

        EXAMPLES::

            sage: Composition([4,1,2,3,5]).to_code()
            [1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0]
        """

        if self == []:
            return [0]

        code = []
        for i in range(len(self)):
            code += [1] + [0]*(self[i]-1)

        return code


    def descents(self, final_descent=False):
        """
        Returns the list of descents of the composition co.

        EXAMPLES::

            sage: Composition([1, 1, 3, 1, 2, 1, 3]).descents()
            [0, 1, 4, 5, 7, 8, 11]
        """
        s = -1
        d = []
        for i in range(len(self)):
            s += self[i]
            d += [s]
        if len(self) != 0 and final_descent:
            d += [len(self)-1]
        return d


    def peaks(self):
        """
        Returns a list of the peaks of the composition self.  The
        peaks of a composition are the descents which do not
        immediately follow another descent.

        EXAMPLES::

            sage: Composition([1, 1, 3, 1, 2, 1, 3]).peaks()
            [4, 7]
        """
        descents = dict((d,True) for d in self.descents())
        return [i+1 for i in range(len(self))
                if i not in descents and i+1 in descents]

    def to_skew_partition(self, overlap=1):
        """
        Returns the skew partition obtained from the composition co. The
        parameter overlap indicates the number of boxes that are covered by
        boxes of the previous line.

        EXAMPLES::

            sage: Composition([3,4,1]).to_skew_partition()
            [[6, 6, 3], [5, 2]]
            sage: Composition([3,4,1]).to_skew_partition(overlap=0)
            [[8, 7, 3], [7, 3]]
        """
        outer = []
        inner = []
        sum_outer = -1*overlap

        for k in range(len(self)-1):
            outer += [ self[k]+sum_outer+overlap  ]
            sum_outer += self[k]-overlap
            inner += [ sum_outer + overlap ]

        if self != []:
            outer += [self[-1]+sum_outer+overlap]
        else:
            return [[],[]]

        return sage.combinat.skew_partition.SkewPartition(
            [ filter(lambda x: x != 0, [l for l in reversed(outer)]),
              filter(lambda x: x != 0, [l for l in reversed(inner)])])


##############################################################


def Compositions(n=None, **kwargs):
    r"""
    Sets of integer Compositions.

    A composition `c` of a nonnegative integer `n` is a list of
    positive integers with total sum `n`.

    See also: ``Composition``, ``Partitions``, ``IntegerVectors``

    EXAMPLES:

    There are 8 compositions of 4::

        sage: Compositions(4).cardinality()
        8

    Here is the list of them::

        sage: list(Compositions(4))
        [[1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3], [2, 1, 1], [2, 2], [3, 1], [4]]

    You can use the ``.first()`` method to get the 'first' composition of
    a number::

        sage: Compositions(4).first()
        [1, 1, 1, 1]

    You can also calculate the 'next' composition given the current
    one::

        sage: Compositions(4).next([1,1,2])
        [1, 2, 1]



    If `n` is not specified, this returns the combinatorial class of
    all (non-negative) integer compositions::

        sage: Compositions()
        Compositions of non-negative integers
        sage: [] in Compositions()
        True
        sage: [2,3,1] in Compositions()
        True
        sage: [-2,3,1] in Compositions()
        False

    If n is specified, it returns the class of compositions of n::

        sage: Compositions(3)
        Compositions of 3
        sage: list(Compositions(3))
        [[1, 1, 1], [1, 2], [2, 1], [3]]
        sage: Compositions(3).cardinality()
        4

    The following examples show how to test whether or not an object
    is a composition::

        sage: [3,4] in Compositions()
        True
        sage: [3,4] in Compositions(7)
        True
        sage: [3,4] in Compositions(5)
        False

    Similarly, one can check whether or not an object is a composition
    which satisfies further constraints::

        sage: [4,2] in Compositions(6, inner=[2,2])
        True
        sage: [4,2] in Compositions(6, inner=[2,3])
        False
        sage: [4,1] in Compositions(5, inner=[2,1], max_slope = 0)
        True

    Note that the given constraints should be compatible::

        sage: [4,2] in Compositions(6, inner=[2,2], min_part=3) #
        True

    The options length, min_length, and max_length can be used to set
    length constraints on the compositions. For example, the
    compositions of 4 of length equal to, at least, and at most 2 are
    given by::

        sage: Compositions(4, length=2).list()
        [[3, 1], [2, 2], [1, 3]]
        sage: Compositions(4, min_length=2).list()
        [[3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: Compositions(4, max_length=2).list()
        [[4], [3, 1], [2, 2], [1, 3]]

    Setting both min_length and max_length to the same value is
    equivalent to setting length to this value::

        sage: Compositions(4, min_length=2, max_length=2).list()
        [[3, 1], [2, 2], [1, 3]]

    The options inner and outer can be used to set part-by-part
    containment constraints. The list of compositions of 4 bounded
    above by [3,1,2] is given by::

        sage: list(Compositions(4, outer=[3,1,2]))
        [[3, 1], [2, 1, 1], [1, 1, 2]]

    Outer sets max_length to the length of its argument. Moreover, the
    parts of outer may be infinite to clear the constraint on specific
    parts. This is the list of compositions of 4 of length at most 3
    such that the first and third parts are at most 1::

        sage: list(Compositions(4, outer=[1,oo,1]))
        [[1, 3], [1, 2, 1]]

    This is the list of compositions of 4 bounded below by [1,1,1]::

        sage: list(Compositions(4, inner=[1,1,1]))
        [[2, 1, 1], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    The options min_slope and max_slope can be used to set constraints
    on the slope, that is the difference `p[i+1]-p[i]` of two
    consecutive parts. The following is the list of weakly increasing
    compositions of 4::

        sage: Compositions(4, min_slope=0).list()
        [[4], [2, 2], [1, 3], [1, 1, 2], [1, 1, 1, 1]]

    Here are the weakly decreasing ones::

        sage: Compositions(4, max_slope=0).list()
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]


    The following is the list of compositions of 4 such that two
    consecutive parts differ by at most one::

        sage: Compositions(4, min_slope=-1, max_slope=1).list()
        [[4], [2, 2], [2, 1, 1], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    The constraints can be combined together in all reasonable ways.
    This is the list of compositions of 5 of length between 2 and 4
    such that the difference between consecutive parts is between -2
    and 1::

        sage: Compositions(5, max_slope=1, min_slope=-2, min_length=2, max_length=4).list()
        [[3, 2], [3, 1, 1], [2, 3], [2, 2, 1], [2, 1, 2], [2, 1, 1, 1], [1, 2, 2], [1, 2, 1, 1], [1, 1, 2, 1], [1, 1, 1, 2]]

    We can do the same thing with an outer constraint::

        sage: Compositions(5, max_slope=1, min_slope=-2, min_length=2, max_length=4, outer=[2,5,2]).list()
        [[2, 3], [2, 2, 1], [2, 1, 2], [1, 2, 2]]

    However, providing incoherent constraints may yield strange
    results. It is up to the user to ensure that the inner and outer
    compositions themselves satisfy the parts and slope constraints.

    Note that if you specify min_part=0, then the objects produced may
    have parts equal to zero. This violates the internal assumptions
    that the Composition class makes. Use at your own risk, or
    preferably consider using ``IntegerVectors`` instead::

        sage: list(Compositions(2, length=3, min_part=0))
        doctest:... RuntimeWarning: Currently, setting min_part=0 produces Composition objects which violate internal assumptions.  Calling methods on these objects may produce errors or WRONG results!
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

        sage: list(IntegerVectors(2, 3))
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

    The generation algorithm is constant amortized time, and handled
    by the generic tool ``IntegerListsLex``.

    TESTS::

        sage: C = Compositions(4, length=2)
        sage: C == loads(dumps(C))
        True

        sage: Compositions(6, min_part=2, length=3)
        Compositions of the integer 6 satisfying constraints length=3, min_part=2


        sage: [2, 1] in Compositions(3, length=2)
        True
        sage: [2,1,2] in Compositions(5, min_part=1)
        True
        sage: [2,1,2] in Compositions(5, min_part=2)
        False

        sage: Compositions(4, length=2).cardinality()
        3
        sage: Compositions(4, min_length=2).cardinality()
        7
        sage: Compositions(4, max_length=2).cardinality()
        4
        sage: Compositions(4, max_part=2).cardinality()
        5
        sage: Compositions(4, min_part=2).cardinality()
        2
        sage: Compositions(4, outer=[3,1,2]).cardinality()
        3

        sage: Compositions(4, length=2).list()
        [[3, 1], [2, 2], [1, 3]]
        sage: Compositions(4, min_length=2).list()
        [[3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: Compositions(4, max_length=2).list()
        [[4], [3, 1], [2, 2], [1, 3]]
        sage: Compositions(4, max_part=2).list()
        [[2, 2], [2, 1, 1], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: Compositions(4, min_part=2).list()
        [[4], [2, 2]]
        sage: Compositions(4, outer=[3,1,2]).list()
        [[3, 1], [2, 1, 1], [1, 1, 2]]
        sage: Compositions(4, outer=[1,oo,1]).list()
        [[1, 3], [1, 2, 1]]
        sage: Compositions(4, inner=[1,1,1]).list()
        [[2, 1, 1], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: Compositions(4, min_slope=0).list()
        [[4], [2, 2], [1, 3], [1, 1, 2], [1, 1, 1, 1]]
        sage: Compositions(4, min_slope=-1, max_slope=1).list()
        [[4], [2, 2], [2, 1, 1], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: Compositions(5, max_slope=1, min_slope=-2, min_length=2, max_length=4).list()
        [[3, 2], [3, 1, 1], [2, 3], [2, 2, 1], [2, 1, 2], [2, 1, 1, 1], [1, 2, 2], [1, 2, 1, 1], [1, 1, 2, 1], [1, 1, 1, 2]]
        sage: Compositions(5, max_slope=1, min_slope=-2, min_length=2, max_length=4, outer=[2,5,2]).list()
        [[2, 3], [2, 2, 1], [2, 1, 2], [1, 2, 2]]
    """
    if n is None:
        assert(len(kwargs) == 0)
        return Compositions_all()
    else:
        if len(kwargs) == 0:
            if isinstance(n, (int,Integer)):
                return Compositions_n(n)
            else:
                raise ValueError, "n must be an integer"
        else:
            # FIXME: should inherit from IntegerListLex, and implement repr, or _name as a lazy attribute
            kwargs['name'] = "Compositions of the integer %s satisfying constraints %s"%(n, ", ".join( ["%s=%s"%(key, kwargs[key]) for key in sorted(kwargs.keys())] ))
            kwargs['element_constructor'] = Composition_class
            if 'min_part' not in kwargs:
                kwargs['min_part'] = 1
            elif kwargs['min_part'] == 0:
                from warnings import warn
                warn("Currently, setting min_part=0 produces Composition objects which violate internal assumptions.  Calling methods on these objects may produce errors or WRONG results!", RuntimeWarning)

            if 'outer' in kwargs:
                kwargs['ceiling'] = kwargs['outer']
                if 'max_length' in kwargs:
                    kwargs['max_length'] = min( len(kwargs['outer']), kwargs['max_length'])
                else:
                    kwargs['max_length'] = len(kwargs['outer'])
                del kwargs['outer']

            if 'inner' in kwargs:
                inner = kwargs['inner']
                kwargs['floor'] = inner
                del kwargs['inner']
                # Should this be handled by integer lists lex?
                if 'min_length' in kwargs:
                    kwargs['min_length'] = max( len(inner), kwargs['min_length'])
                else:
                    kwargs['min_length'] = len(inner)
            return IntegerListsLex(n, **kwargs)


# Allows to unpickle old constrained Compositions_constraints objects.
class Compositions_constraints(IntegerListsLex):
    def __setstate__(self, data):
        """
        TESTS::

            # This is the unpickling sequence for Compositions(4, max_part=2) in sage <= 4.1.1
            sage: pg_Compositions_constraints = unpickle_global('sage.combinat.composition', 'Compositions_constraints')
            sage: si = unpickle_newobj(pg_Compositions_constraints, ())
            sage: pg_make_integer = unpickle_global('sage.rings.integer', 'make_integer')
            sage: unpickle_build(si, {'constraints':{'max_part':pg_make_integer('2')}, 'n':pg_make_integer('4')})
            sage: si
            Integer lists of sum 4 satisfying certain constraints
            sage: si.list()
            [[2, 2], [2, 1, 1], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        """
        n = data['n']
        self.__class__ = IntegerListsLex
        constraints = {
                       'min_part' : 1,
                       'element_constructor' : Composition_class}
        constraints.update(data['constraints'])
        self.__init__(n, **constraints)

class Compositions_all(InfiniteAbstractCombinatorialClass):
    def __init__(self):
        """
        TESTS::

            sage: C = Compositions()
            sage: C == loads(dumps(C))
            True
        """
        pass
    def __repr__(self):
        """
        TESTS::

            sage: repr(Compositions())
            'Compositions of non-negative integers'
        """
        return "Compositions of non-negative integers"

    object_class = Composition_class

    def __contains__(self, x):
        """
        TESTS::

            sage: [2,1,3] in Compositions()
            True
            sage: [] in Compositions()
            True
            sage: [-2,-1] in Compositions()
            False
            sage: [0,0] in Compositions()
            True
        """
        if isinstance(x, Composition_class):
            return True
        elif isinstance(x, __builtin__.list):
            for i in range(len(x)):
                if not isinstance(x[i], (int, Integer)):
                    return False
                if x[i] < 0:
                    return False
            return True
        else:
            return False

    def _infinite_cclass_slice(self, n):
        """
        Needed by InfiniteAbstractCombinatorialClass to build __iter__.

        TESTS::

            sage: Compositions()._infinite_cclass_slice(4) == Compositions(4)
            True
            sage: it = iter(Compositions())    # indirect doctest
            sage: [it.next() for i in range(10)]
            [[], [1], [1, 1], [2], [1, 1, 1], [1, 2], [2, 1], [3], [1, 1, 1, 1], [1, 1, 2]]
        """
        return Compositions_n(n)


class Compositions_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS::

            sage: C = Compositions(3)
            sage: C == loads(dumps(C))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS::

            sage: repr(Compositions(3))
            'Compositions of 3'
        """
        return "Compositions of %s"%self.n

    def __contains__(self, x):
        """
        TESTS::

            sage: [2,1,3] in Compositions(6)
            True
            sage: [2,1,2] in Compositions(6)
            False
            sage: [] in Compositions(0)
            True
            sage: [0] in Compositions(0)
            True
        """
        return x in Compositions() and sum(x) == self.n

    def cardinality(self):
        """
        TESTS::

            sage: Compositions(3).cardinality()
            4
            sage: Compositions(0).cardinality()
            1
        """
        if self.n >= 1:
            return 2**(self.n-1)
        elif self.n == 0:
            return 1
        else:
            return 0

    def list(self):
        """
        TESTS::

            sage: Compositions(4).list()
            [[1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3], [2, 1, 1], [2, 2], [3, 1], [4]]
            sage: Compositions(0).list()
            [[]]
        """
        if self.n == 0:
            return [Composition_class([])]

        result = []
        for i in range(1,self.n+1):
            result += map(lambda x: [i]+x[:], Compositions_n(self.n-i).list())

        return [Composition_class(r) for r in result]




# Those belong to the Composition class

def from_descents(descents, nps=None):
    """
    Returns a composition from the list of descents.

    EXAMPLES::

        sage: Composition([1, 1, 3, 4, 3]).descents()
        [0, 1, 4, 8, 11]
        sage: sage.combinat.composition.from_descents([1,0,4,8],12)
        [1, 1, 3, 4, 3]
        sage: sage.combinat.composition.from_descents([1,0,4,8,11])
        [1, 1, 3, 4, 3]
    """

    d = [x+1 for x in descents]
    d.sort()

    if d == []:
        if nps == 0:
            return []
        else:
            return [nps]

    if nps is not None:
        if nps < max(d):
            #Error: d is not included in [1,...,nps-1]
            return None
        elif nps > max(d):
            d.append(nps)

    co = [d[0]]
    for i in range(len(d)-1):
        co += [ d[i+1]-d[i] ]

    return Composition(co)

def from_code(code):
    """
    Return the composition from its code.The code of a composition is a
    list of length self.size() of 1s and 0s such that there is a 1
    wherever a new part starts.

    EXAMPLES::

        sage: import sage.combinat.composition as composition
        sage: Composition([4,1,2,3,5]).to_code()
        [1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0]
        sage: composition.from_code(_)
        [4, 1, 2, 3, 5]
        sage: Composition([3,1,2,3,5]).to_code()
        [1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0]
        sage: composition.from_code(_)
        [3, 1, 2, 3, 5]
    """
    if code == [0]:
        return []

    L = filter(lambda x: code[x]==1, range(len(code))) #the positions of the letter 1
    return Composition([L[i]-L[i-1] for i in range(1, len(L))] + [len(code)-L[-1]])

