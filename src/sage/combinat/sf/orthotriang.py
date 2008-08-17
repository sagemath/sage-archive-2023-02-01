"""
Symmetric functions defined by orthogonality and triangularity.


One characterization of Schur functions is that they are upper triangularly
related to the monomial symmetric functions and orthogonal with respect to
the Hall scalar product.  We can use the class
SymmetricFunctionAlgebra_orthotriang to obtain the Schur functions from
this definition.

sage: from sage.combinat.sf.sfa import zee
sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
sage: m = SFAMonomial(QQ)
sage: s =  SymmetricFunctionAlgebra_orthotriang(QQ, m, zee, 's', 'Schur functions')
sage: s([2,1])^2
s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

sage: s2 = SFASchur(QQ)
sage: s2([2,1])^2
s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.combinat.combinat import CombinatorialClass
from sage.libs.symmetrica.all import hall_littlewood
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
import sfa
import sage.combinat.partition
from sage.matrix.all import matrix, MatrixSpace
from sage.rings.all import ZZ, QQ

class SymmetricFunctionAlgebra_orthotriang(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, base, scalar, prefix, name, leading_coeff=None):
        """
        EXAMPLES:
            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: m = SFAMonomial(QQ)
            sage: s =  SymmetricFunctionAlgebra_orthotriang(QQ, m, zee, 's', 'Schur functions')
            sage: s == loads(dumps(s))
            True
        """
        self._combinatorial_class = sage.combinat.partition.Partitions()
        self._one = sage.combinat.partition.Partition_class([])
        self._sf_base = base
        self._element_class = SymmetricFunctionAlgebraElement_orthotriang
        self._scalar = scalar
        self._prefix = prefix
        self._name = name
        self._leading_coeff = leading_coeff
        CombinatorialAlgebra.__init__(self, R)

        self._self_to_base_cache = {}
        self._base_to_self_cache = {}

    def _coerce_start(self, x):
        """
        Coerce other symmetric functions into self through the base.

        EXAMPLES:
            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: m = SFAMonomial(QQ)
            sage: s =  SymmetricFunctionAlgebra_orthotriang(QQ, m, zee, 's', 'Schur functions')
            sage: s._coerce_start(m([2,1]))
            -2*s[1, 1, 1] + s[2, 1]
        """
        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._sf_base(x)
            return self._from_cache(x, self._base_cache, self._base_to_self_cache)
        else:
            raise TypeError

    def _base_cache(self, n):
        """
        Computes the change of basis between self and base for the
        homogenous component of size n.

        EXAMPLES:
            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: m = SFAMonomial(QQ)
            sage: s =  SymmetricFunctionAlgebra_orthotriang(QQ, m, zee, 's', 'Schur functions')
            sage: s._base_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(s._base_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -1), ([2], 1)])]
            sage: l(s._self_to_base_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], 1), ([2], 1)])]
        """
        if n in self._self_to_base_cache:
            return
        else:
            self._self_to_base_cache[n] = {}

        self._gram_schmidt(n, self._sf_base, self._scalar, self._self_to_base_cache,\
                           leading_coeff=self._leading_coeff, upper_triangular=True)
        self._invert_morphism(n, self.base_ring(), self._self_to_base_cache, \
                              self._base_to_self_cache, to_other_function = self._to_base)

    def _to_base(self, part):
        """
        Returns a function which takes in a partition part2 and returns
        the coefficient of part2 in the expansion of self(part) in
        base.

        NOTE:
            We assume that self._gram_schmidt has been called before
            self._to_base is called.

        EXAMPLES:
            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: m = SFAMonomial(QQ)
            sage: s =  SymmetricFunctionAlgebra_orthotriang(QQ, m, zee, 's', 'Schur functions')
            sage: m(s([2,1]))
            2*m[1, 1, 1] + m[2, 1]
            sage: f = s._to_base(Partition([2,1]))
            sage: [f(p) for p in Partitions(3)]
            [0, 1, 2]
        """
        f = lambda part2: self._self_to_base_cache[part].get(part2, 0)
        return f

    def _multiply(self, left, right):
        """
        Returns left*right by converting both to the base and then
        converting back to self.

        EXAMPLES:
            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: m = SFAMonomial(QQ)
            sage: s =  SymmetricFunctionAlgebra_orthotriang(QQ, m, zee, 's', 'Schur functions')
            sage: s([1])*s([2,1]) #indirect doctest
            s[2, 1, 1] + s[2, 2] + s[3, 1]
        """
        return self( self._sf_base(left)*self._sf_base(right) )


class SymmetricFunctionAlgebraElement_orthotriang(sfa.SymmetricFunctionAlgebraElement_generic):
    pass
