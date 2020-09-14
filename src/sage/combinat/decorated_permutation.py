# -*- coding: utf-8 -*-
r"""
Decorated permutations

AUTHORS:

- Martin Rubey (2020): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2020 Martin Rubey <martin.rubey at tuwien.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.arith.all import factorial
from sage.rings.integer import Integer
from .permutation import Permutations
from .subset import Subsets
from .colored_permutations import SignedPermutations

class DecoratedPermutation(Element,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A decorated permutation.

    A decorated permutation is a signed permutation where all
    non-fixed points have positive sign.
    """
    @staticmethod
    def __classcall_private__(cls, pi, check=True):
        """
        Create a decorated permutation.

        EXAMPLES::

            sage: DecoratedPermutation([2, 1, 3])
            [2, 1, 3]

            sage: DecoratedPermutation([2, 1, -3])
            [2, 1, -3]

        """
        pi = list(pi)
        return DecoratedPermutations(len(pi))(pi, check=check)

    def __init__(self, parent, pi):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: elt = S([2, 1, -3])
            sage: TestSuite(elt).run()
        """
        self._pi = pi
        Element.__init__(self, parent)

    def __hash__(self):
        r"""
        TESTS::

            sage: S = DecoratedPermutations(3)
            sage: elt = S([2, 1, -3])
            sage: hash(elt) # random
            915443076393556996
        """
        return hash(self._pi)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: DecoratedPermutation([2, 1, -3])
            [2, 1, -3]
        """
        return repr(list(self._pi))

    def _richcmp_(self, other, op):
        """
        Do the comparison.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: elt1 = S([2, 1, -3])
            sage: elt2 = DecoratedPermutation([2, 1, -3])
            sage: elt1 == elt2
            True
        """
        return richcmp(self._pi, other._pi, op)

    def size(self):
        """
        Return the size of the decorated permutation.

        EXAMPLES::

            sage: DecoratedPermutation([2, 1, -3]).size()
            3
        """
        return len(self._pi)

    def __iter__(self):
        """
        Iterate over the values of the decorated permutation.

        EXAMPLES::

            sage: list(DecoratedPermutation([2, 1, -3]))
            [2, 1, -3]
        """
        return iter(self._pi)                 
    
    def to_signed_permutation(self):
        """
        Return ``self`` as a signed permutation.

        EXAMPLES::

            sage: DecoratedPermutation([2, 1, -3]).to_signed_permutation()
            [2, 1, -3]
        """
        return SignedPermutations(len(self._pi))(list(self._pi))

class DecoratedPermutations(UniqueRepresentation, Parent):
    r"""
    Class of all decorated permutations of `n`.

    A decorated permutation is a signed permutation where all
    non-fixed points have positive sign.

    INPUT:

    - `n` -- an integer, the size of the decorated permutations.

    EXAMPLES:

    This will create an instance to manipulate the decorated
    permutations of size 3::

        sage: S = DecoratedPermutations(3); S
        Decorated permutations of size 3
        sage: S.cardinality()
        16

    """
    def __init__(self, n):
        r"""
        Initialize ``self``.

        TESTS::

            sage: S = DecoratedPermutations(4)
            sage: TestSuite(S).run()
        """
        self._n = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: DecoratedPermutations(4)
            Decorated permutations of size 4
        """
        return "Decorated permutations of size %s" % self._n

    def __contains__(self, pi):
        """
        Check if ``pi`` is in ``self``.

        TESTS::

            sage: S = DecoratedPermutations(3)
            sage: [2, 1, -3] in S
            True
            sage: [2, -1, 3] in S
            False
        """
        if isinstance(pi, DecoratedPermutation):
            return len(pi._pi) == self._n

        values = [v for v in pi]
        if len(values) != self._n:
            return False
        abs_values = [abs(v) for v in values]
        for i, (v, abs_v) in enumerate(zip(values, abs_values), 1):
            if i != abs_v and v < 0:
                return False
        return sorted(abs_values) == list(range(1, self._n+1))

    def _element_constructor_(self, pi, check=True):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: elt = S([2, 1, -3])
            sage: elt.parent() is S
            True
        """
        if isinstance(pi, DecoratedPermutation):
            if pi.parent() is self:
                return pi
            raise ValueError("Cannot convert between decorated permutations of different sizes")

        pi = tuple(pi)
        if check and pi not in self:
            raise ValueError('invalid decorated permutation')
        return self.element_class(self, pi)

    Element = DecoratedPermutation

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: S = DecoratedPermutations(3)
            sage: S._an_element_()
            [1, 2, 3]

        """
        return self.element_class(self, list(range(1, self._n + 1)))

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of decorated permutations of size `n` is equal to

        .. MATH::

            \sum_{k=0^n} \frac{n!}{k!}

        EXAMPLES::

            sage: [DecoratedPermutations(n).cardinality() for n in range(11)]
            [1, 2, 5, 16, 65, 326, 1957, 13700, 109601, 986410, 9864101]
        """
        return Integer(sum(factorial(self._n)//factorial(k) for k in range(self._n+1)))

    def __iter__(self):
        r"""
        Iterator on the decorated permutations of size `n`.

        TESTS::

            sage: S = DecoratedPermutations(2); S.list()
            [[1, 2], [-1, 2], [1, -2], [-1, -2], [2, 1]]
            sage: sum(1 for a in S)
            5
        """
        for sigma in Permutations(self._n):
            F = sigma.fixed_points()
            for X in Subsets(F):
                tau = list(sigma)
                for i in X:
                    tau[i-1] = -tau[i-1]
                yield DecoratedPermutations(self._n)(tau)
