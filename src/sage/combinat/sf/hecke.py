r"""
Hecke Character Basis

The basis of symmetric functions given by characters of the
Hecke algebra (of type `A`).

AUTHORS:

- Travis Scrimshaw (2017-08): Initial version
"""
#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
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

from __future__ import absolute_import

from sage.combinat.partition import _Partitions
from sage.combinat.sf.multiplicative import SymmetricFunctionAlgebra_multiplicative
from sage.matrix.all import matrix
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.rings.integer_ring import ZZ

class HeckeCharacter(SymmetricFunctionAlgebra_multiplicative):
    r"""
    Basis of the symmetric functions that gives the characters of the
    Hecke aglebra in analogy to the Frobenious formula for the
    symmetric group.

    Consider the Hecke algebra `H_n(q)` with quadratic relations::

        T_i^2 = (q - 1) T_i + q.

    Let `\mu` be a partition of `n` with length `\ell`. The character
    `\chi` of a `H_n(q)`-representation is completely determined by
    the elements `T_{\gamma_{\mu}}`, where::

        \gamma_{\mu} = (\mu_1, \ldots, 1) (\mu_2 + \mu_1, \ldots, 1 + \mu_1)
        \cdots (n, \ldots, 1 + \sum_{i < \ell} \mu_i),

    (written in cycle notation). We define a basis of the symmetric
    functions by

        \bar{q}_{\mu} = \sum_{\lambda \vdash n}
        \chi^{\lambda}(T_{\gamma_{\mu}}) s_{\lambda}.

    INPUT:

    - ``sym`` -- the ring of symmetric functions
    - ``q`` -- (default: ``'q'``) the parameter `q`

    EXAMPLES::

        sage: q = ZZ['q'].fraction_field().gen()
        sage: Sym = SymmetricFunctions(q.parent())
        sage: qbar = Sym.hecke_character(q)
        sage: qbar[2] * qbar[3] * qbar[3,1]
        qbar[3, 3, 2, 1]

        sage: s = Sym.s()
        sage: s(qbar[2,1])
        -s[1, 1, 1] + (q-1)*s[2, 1] + q*s[3]
        sage: qbar(s[2,1])
        (q/(q^2+q+1))*qbar[1, 1, 1] + ((q-1)/(q^2+q+1))*qbar[2, 1]
         - (1/(q^2+q+1))*qbar[3]

    We compute character tables for Hecke algebras, which correspond
    to the transition matrix from the `\bar{q}` basis to the Schur
    basis::

        sage: qbar.transition_matrix(s, 1)
        [1]
        sage: qbar.transition_matrix(s, 2)
        [ q -1]
        [ 1  1]
        sage: qbar.transition_matrix(s, 3)
        [  q^2    -q     1]
        [    q q - 1    -1]
        [    1     2     1]
        sage: qbar.transition_matrix(s, 4)
        [      q^3      -q^2         0         q        -1]
        [      q^2   q^2 - q        -q    -q + 1         1]
        [      q^2 q^2 - 2*q   q^2 + 1  -2*q + 1         1]
        [        q   2*q - 1     q - 1     q - 2        -1]
        [        1         3         2         3         1]

    REFERENCES:

    - [Ram1991]_
    - [RR1997]_
    """
    def __init__(self, sym, q='q'):
        r"""
        Initialize ``self``.

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(ZZ['q']))
            sage: qbar = Sym.qbar()
            sage: TestSuite(qbar).run()

            sage: Sym = SymmetricFunctions(QQ)
            sage: qbar = Sym.qbar(q=2)
            sage: TestSuite(qbar).run()
        """
        self.q = sym.base_ring()(q)
        SymmetricFunctionAlgebra_multiplicative.__init__(self, sym,
            basis_name="Hecke character with q={}".format(self.q),
            prefix="qbar")
        self._s = sym.schur()
        self._self_to_s_cache = {}
        self._s_to_self_cache = {}

        # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
        category = ModulesWithBasis(self._sym.base_ring())
        self   .register_coercion(SetMorphism(Hom(self._s, self, category), self._s_to_self))
        self._s.register_coercion(SetMorphism(Hom(self, self._s, category), self._self_to_s))

    def _s_to_self(self, x):
        r"""
        Isomorphism from the Schur basis into ``self``

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis
        - ``x`` -- an element of the Schur basis

        OUTPUT:

        - an element of ``self`` equivalent to ``x``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: qbar = Sym.qbar(q=2)
            sage: s = Sym.schur()
            sage: qbar._s_to_self(s[2,1])
            2/7*qbar[1, 1, 1] + 1/7*qbar[2, 1] - 1/7*qbar[3]

        This is for internal use only. Instead use::

            sage: qbar(s[2,1])
            2/7*qbar[1, 1, 1] + 1/7*qbar[2, 1] - 1/7*qbar[3]
        """
        return self._from_cache(x, self._s_cache, self._s_to_self_cache, q=self.q)

    def _self_to_s(self, x):
        r"""
        Isomorphism from ``self`` to the Schur basis

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis
        - ``x`` -- an element of the basis ``self``

        OUTPUT:

        - an element of the Schur basis equivalent to ``x``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: qbar = Sym.qbar(q=2)
            sage: s = Sym.schur()
            sage: qbar._self_to_s(qbar[2,1])
            -s[1, 1, 1] + s[2, 1] + 2*s[3]

        This is for internal use only. Instead use::

            sage: s(qbar[2,1])
            -s[1, 1, 1] + s[2, 1] + 2*s[3]
        """
        return self._s._from_cache(x, self._s_cache, self._self_to_s_cache, q=self.q)

    def _qbar_to_s(self, mu):
        r"""
        Return `\bar{q}_{\mu}` in the Schur basis.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: Sym = SymmetricFunctions(q.parent())
            sage: qbar = Sym.hecke_character()
            sage: qbar._qbar_to_s(Partition([2]))
            -s[1, 1] + q*s[2]
            sage: qbar._qbar_to_s(Partition([4]))
            -s[1, 1, 1, 1] + q*s[2, 1, 1] - q^2*s[3, 1] + q^3*s[4]
            sage: qbar._qbar_to_s(Partition([2,1]))
            -s[1, 1, 1] + (q-1)*s[2, 1] + q*s[3]
        """
        if not mu:
            return self._s.one()
        ret = self._s.one()
        q = self.q
        mone = -self.base_ring().one()
        return self._s.prod(sum(mone**(r-m) * q**(m-1)
                                * self._s[_Partitions([m] + [1]*(r-m))]
                                for m in range(1, r+1))
                            for r in mu)

    def _s_cache(self, n):
        r"""
        Compute the change of basis between the `\bar{q}` polynomials
        and the Schur functions for partitions of size ``n``.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``n`` -- positive integer

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(ZZ['q']))
            sage: qbar = Sym.qbar()
            sage: qbar._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(qbar._s_to_self_cache[2])
            [([1, 1], [([1, 1], q/(q + 1)), ([2], -1/(q + 1))]),
             ([2], [([1, 1], 1/(q + 1)), ([2], 1/(q + 1))])]
            sage: l(qbar._self_to_s_cache[2])
            [([1, 1], [([1, 1], 1), ([2], 1)]), ([2], [([1, 1], -1), ([2], q)])]
            sage: qbar = Sym.qbar()
            sage: qbar._s_cache(3)
            sage: l(qbar._s_to_self_cache[2])
            [([1, 1], [([1, 1], q/(q + 1)), ([2], -1/(q + 1))]),
             ([2], [([1, 1], 1/(q + 1)), ([2], 1/(q + 1))])]
        """
        def coeff(part):
            return lambda mu: self._qbar_to_s(part)[mu]
        self._s._invert_morphism(n, self.base_ring(),
                                 self._s_to_self_cache,
                                 self._self_to_s_cache,
                                 to_self_function=coeff)

    def coproduct_on_generators(self, r):
        r"""
        Return the coproduct on the generator `\bar{q}_r` of ``self``.

        Define the coproduct on `\bar{q}_r` by

        .. MATH::

            \Delta(\bar{q}_r) = \bar{q}_0 \otimes \bar{q}_r
            + (q - 1) \sum_{j=1}^{r-1} \bar{q}_j \otimes \bar{q}_{r-j}
            + \bar{q}_r \otimes \bar{q}_0

        .. WARNING::

            The formula given here is not been formally proven.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: Sym = SymmetricFunctions(q.parent())
            sage: qbar = Sym.hecke_character()
            sage: s = Sym.s()
            sage: qbar[2].coproduct()
            qbar[] # qbar[2] + (q-1)*qbar[1] # qbar[1] + qbar[2] # qbar[]
        """
        def P(i): return _Partitions([i]) if i else _Partitions([])
        T = self.tensor_square()
        one = self.base_ring().one()
        q = self.q
        return T.sum_of_terms(((P(j), P(r-j)), one if j in [0,r] else q-one)
                              for j in range(r+1))

