"""
Hochschild Complexes
"""

########################################################################
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.category_types import ChainComplexes
from sage.categories.tensor import tensor
from sage.combinat.free_module import CombinatorialFreeModule
#from sage.rings.all import ZZ
#from sage.matrix.constructor import matrix
#from sage.homology.chain_complex import ChainComplex_class

import itertools

class HochschildComplex(UniqueRepresentation, Parent):
    r"""
    The Hochschild complex.

    Let `A` be an algebra over a commutative ring `R` which is a
    projective `R`-module, and `M` an `A`-bimodule. The *Hochschild complex*
    is the chain complex given by

    .. MATH::

        C_n(A, M) := M \otimes A^{\otimes n}

    with the boundary operator

    .. MATH::

        d_i(m \otimes a_1 \otimes \cdots \otimes a_n) = \begin{cases}
        m a_1 \otimes \cdots \otimes a_n & \text{if } i = 0, \\
        (-1)^n a_n m \otimes a_1 \otimes \cdots \otimes a_{n-1}
        & \text{if } i = n, \\
        (-1)^i m \otimes a_1 \otimes \cdots \otimes a_i a_{i+1} \otimes
        \cdots \otimes a_n & \text{otherwise.}
        \end{cases}

    *Hochschild homology* is the homology of this complex, and
    alternatively, the Hochschild homology can be described by
    `HH_n(A, M) = \operatorname{Tor}_n^A(A, M)`.

    *Hochschild cohomology* is the homology of the dual complex and
    can be described by `HH^n(A, M) = \operatorname{Ext}^n_A(A, M)`.

    REFERENCES:

    - :wikipedia:`Hochschild_homology`
    """
    def __init__(self, A, M):
        """
        Initialize ``self``.

        EXAMPLES::
        """
        self._A = A
        self._M = M
        Parent.__init__(self, category=ChainComplexes(A.base_ring()))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        """
        return "Hochschild complex of {} with coefficients in {}".format(self._A, self._M)

    def free_module(self, d):
        """
        Return the free module in degree ``d``
        """
        return tensor([self._M] + [self._A]*d)

    @cached_method
    def trivial_module(self):
        """
        Return the trivial module of ``self``.
        """
        return CombinatorialFreeModule(self._A.base_ring(), [])

    def boundary(self, d):
        """
        Return the boundary operator in degree ``d``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: H = HochschildComplex(SGA, SGA)
            sage: phi = H.boundary(2)
            sage: psi = H.boundary(3)
            sage: comp = phi * psi
            sage: for b in H.free_module(3).basis():
            ....:     assert comp(b) == 0
        """
        R = self._A.base_ring()
        one = R.one()
        if d == 0:
            t = self.trivial_module()
            zero = t.zero()
            return self.free_module(0).module_morphism(lambda x: zero, codomain=t)
        Fd = self.free_module(d-1)
        Fd1 = self.free_module(d)
        mone = -one
        def on_basis(k):
            p = self._M.monomial(k[0]) * self._A.monomial(k[1])
            ret = Fd._from_dict({(m,) + k[2:]: c for m,c in p}, remove_zeros=False)
            for i in range(1, d):
                p = self._A.monomial(k[i]) * self._A.monomial(k[i+1])
                ret += mone**i * Fd._from_dict({k[:i] + (m,) + k[i+2:]: c
                                                   for m,c in p}, remove_zeros=False)
            p = self._A.monomial(k[-1]) * self._M.monomial(k[0])
            ret += mone**d * Fd._from_dict({(m,) + k[1:-1]: c for m,c in p},
                                           remove_zeros=False)
            return ret
        return Fd1.module_morphism(on_basis, codomain=Fd)

    def coboundary(self, d):
        """
        Return the coboundary morphism of degree ``d``.
        """
        if self._A.category() is not self._A.category().FiniteDimensional():
            raise NotImplementedError("the algebra must be finite dimensional")
        bdry = self.boundary(d)
        dom = bdry.domain()
        cod = bdry.codomain()
        return cod.module_morphism(matrix=bdry.matrix().transpose(), codomain=dom)

    def homology(self, d):
        """
        Return the ``d``-th homology group.
        """
        if self._A.category() is not self._A.category().FiniteDimensional():
            raise NotImplementedError("the algebra must be finite dimensional")
        bdry = self.boundary(d)
        bdry1 = self.boundary(d+1)
        ker = bdry.kernel()
        im_retract = ker.submodule([ker.retract(b) for b in bdry1.image_basis()])
        return ker.quotient_module(im_retract)

    def cohomology(self, d):
        """
        Return the ``d``-th cohomology group.
        """
        cb = self.coboundary(d)
        cb1 = self.coboundary(d+1)
        ker = cb1.kernel()
        im_retract = ker.submodule([ker.retract(b) for b in cb.image_basis()])
        return ker.quotient_module(im_retract)

# This is superceded by #19613
class TrivialRepresentation(CombinatorialFreeModule):
    def __init__(self, GA):
        self._group_algebra = GA
        CombinatorialFreeModule.__init__(self, GA.base_ring(), ['v'])

    def _repr_(self):
        return "Trivial representation of {}".format(self._group_algebra)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(QQ, 3)
                sage: V = TrivialRepresentation(SGA)
                sage: x = V.an_element()
                sage: all(x * b == x for b in SGA.basis())
                True
                sage: all(b * x == x for b in SGA.basis())
                True
            """
            if (isinstance(scalar, Element)
                    and scalar.parent() is self.parent()._group_algebra):
                d = self.monomial_coefficients(copy=True)
                d['v'] *= sum(scalar.coefficients())
                return self.parent()._from_dict(d)
            return CombinatorialFreeModule.Element._acted_upon_(scalar, self_on_left)

        _rmul_ = _lmul_ = _acted_upon_

