"""
Orthogonal Symmetric Functions

AUTHORS:

- Travis Scrimshaw (2013-11-10): Initial version
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
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
import sfa
import sage.combinat.partition as partition
from sage.misc.cachefunc import cached_method
from sage.rings.all import ZZ, QQ, Integer
from sage.matrix.all import matrix

class SymmetricFunctionAlgebra_orthogonal(sfa.SymmetricFunctionAlgebra_generic):
    r"""
    The orthogonal symmetric function basis (or orthogonal basis, to be short).

    The orthogonal basis `\{ o_{\lambda} \}` where `\lambda` is taken over
    all partitions is defined by the following change of basis with the
    Schur functions:

    .. MATH::

        s_{\lambda} = \sum_{\mu} \left( \sum_{\nu \in H} c^{\lambda}_{\mu\nu}
        \right) o_{\mu}

    where `H` is the set of all partitions with even-width rows and
    `c^{\lambda}_{\mu\nu}` is the usual Littlewood-Richardson (LR)
    coefficients. By the properties of LR coefficients, this can be shown to
    be a upper unitriangular change of basis.

    .. NOTE::

        This is only a filtered basis, not a `\ZZ`-graded basis. However this
        does respect the induced `(\ZZ/2\ZZ)`-grading.

    INPUT:

    - ``Sym`` -- an instance of the ring of the symmetric functions

    REFERENCES:

    - [ChariKleber2000]_

    EXAMPLES:

    Here are the first few orthogonal symmetric functions, in various bases::

        sage: Sym = SymmetricFunctions(QQ)
        sage: o = Sym.o()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: p = Sym.p()
        sage: s = Sym.s()
        sage: m = Sym.m()

        sage: p(o([1]))
        p[1]
        sage: m(o([1]))
        m[1]
        sage: e(o([1]))
        e[1]
        sage: h(o([1]))
        h[1]
        sage: s(o([1]))
        s[1]

        sage: p(o([2]))
        -p[] + 1/2*p[1, 1] + 1/2*p[2]
        sage: m(o([2]))
        -m[] + m[1, 1] + m[2]
        sage: e(o([2]))
        -e[] + e[1, 1] - e[2]
        sage: h(o([2]))
        -h[] + h[2]
        sage: s(o([2]))
        -s[] + s[2]

        sage: p(o([3]))
        -p[1] + 1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]
        sage: m(o([3]))
        -m[1] + m[1, 1, 1] + m[2, 1] + m[3]
        sage: e(o([3]))
        -e[1] + e[1, 1, 1] - 2*e[2, 1] + e[3]
        sage: h(o([3]))
        -h[1] + h[3]
        sage: s(o([3]))
        -s[1] + s[3]

        sage: Sym = SymmetricFunctions(ZZ)
        sage: o = Sym.o()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: s = Sym.s()
        sage: m = Sym.m()
        sage: p = Sym.p()
        sage: m(o([4]))
        -m[1, 1] + m[1, 1, 1, 1] - m[2] + m[2, 1, 1] + m[2, 2] + m[3, 1] + m[4]
        sage: e(o([4]))
        -e[1, 1] + e[1, 1, 1, 1] + e[2] - 3*e[2, 1, 1] + e[2, 2] + 2*e[3, 1] - e[4]
        sage: h(o([4]))
        -h[2] + h[4]
        sage: s(o([4]))
        -s[2] + s[4]

    Some examples of conversions the other way::

        sage: o(h[3])
        o[1] + o[3]
        sage: o(e[3])
        o[1, 1, 1]
        sage: o(m[2,1])
        o[1] - 2*o[1, 1, 1] + o[2, 1]
        sage: o(p[3])
        o[1, 1, 1] - o[2, 1] + o[3]

    Some multiplication::

        sage: o([2,1,1]) * o([2])
        o[1, 1] + o[1, 1, 1, 1] + 2*o[2, 1, 1] + o[2, 2] + o[2, 2, 1, 1]
         + o[3, 1] + o[3, 1, 1, 1] + o[3, 2, 1] + o[4, 1, 1]
        sage: o([1,1]) * o([2,1])
        o[1] + o[1, 1, 1] + 2*o[2, 1] + o[2, 1, 1, 1] + o[2, 2, 1]
         + o[3] + o[3, 1, 1] + o[3, 2]

    Examples of the Hopf algebra structure::

        sage: o([1]).antipode()
        -o[1]
        sage: o([2]).antipode()
        -o[] + o[1, 1]
        sage: o([1]).coproduct()
        o[] # o[1] + o[1] # o[]
        sage: o([2]).coproduct()
        o[] # o[] + o[] # o[2] + o[1] # o[1] + o[2] # o[]
        sage: o([1]).counit()
        0
        sage: o.one().counit()
        1
    """
    def __init__(self, Sym):
        """
        Initialize ``self``.

        TESTS::

            sage: o = SymmetricFunctions(QQ).o()
            sage: TestSuite(o).run()
        """
        sfa.SymmetricFunctionAlgebra_generic.__init__(self, Sym, "orthogonal",
                                                      'o', graded=False)

        # Setup the coercions
        s = Sym.schur()
        M = s.module_morphism(self._s_to_o_on_basis, codomain=self,
                              triangular='upper', unitriangular=True)
        M.register_as_coercion()
        (~M).register_as_coercion()

    def _multiply(self, r, l):
        r"""
        Return the product of ``r`` and ``l`` by coercing to the Schur
        functions and then pulling back.

        EXAMPLES::

            sage: o = SymmetricFunctions(QQ).o()
            sage: o([2]) * o([1,1]) # indirect doctest
            o[1, 1] + o[2] + o[2, 1, 1] + o[3, 1]
        """
        s = self.realization_of().schur()
        return self(s(r) * s(l))

    def counit(self, x):
        r"""
        Return the counit of ``x`` in ``self``.

        EXAMPLES::

            sage: o = SymmetricFunctions(QQ).o()
            sage: o.an_element()
            2*o[] + 2*o[1] + 3*o[2]
            sage: o.counit(o.an_element())
            -1
        """
        s = self.realization_of().schur()
        return s.counit(s(x))

    @cached_method
    def _s_to_o_on_basis(self, lam):
        r"""
        Return the Schur symmetric function ``s[lam]`` expanded in
        the orthogonal basis, where ``lam`` is a partition.

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - the expansion of ``s[lam]`` in the orthogonal basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.schur()
            sage: o = Sym.orthogonal()
            sage: o._s_to_o_on_basis(Partition([]))
            o[]
            sage: o._s_to_o_on_basis(Partition([4,2,1]))
            o[1] + 2*o[2, 1] + o[2, 2, 1] + o[3] + o[3, 1, 1] + o[3, 2] + o[4, 1] + o[4, 2, 1]
            sage: s(o._s_to_o_on_basis(Partition([3,1]))) == s[3,1]
            True
        """
        import sage.libs.lrcalc.lrcalc as lrcalc
        from sage.combinat.partition import Partitions
        R = self.base_ring()
        # TODO: This is a brute force check and could likely be simplified
        return self._from_dict({ mu: R.sum( lrcalc.lrcoef_unsafe(lam, mu, nu)
                                            for j in range(sum(lam)-sum(mu)+1)
                                            for nu in Partitions(j)
                                                if all(x % 2 == 0 for x in nu) )
                                 for i in range(sum(lam)+1) for mu in Partitions(i) })

