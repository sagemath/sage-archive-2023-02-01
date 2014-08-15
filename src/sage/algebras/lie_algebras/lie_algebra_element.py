"""
Lie Algebra Elements

AUTHORS:

- Travis Scrimshaw (2005-05-04): Initial implementation
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.misc import repr_lincomb
from copy import copy
from functools import total_ordering
from sage.structure.element import ModuleElement, RingElement, coerce_binop
from sage.structure.sage_object import SageObject
from sage.combinat.free_module import CombinatorialFreeModuleElement
from sage.structure.element_wrapper import ElementWrapper

# TODO: Have the other classes inherit from this?
# TODO: Should this be a mixin class (or moved to the category)?
class LieAlgebraElement_generic(ModuleElement):
    """
    Generic methods for all Lie algebra elements.
    """
    def _mul_(self, other):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::
        """
        if self == 0 or other == 0:
            return self.parent().zero()
        # Otherwise we lift to the UEA
        return self.lift() * other.lift()

# TODO: factor out parts of CombinatorialFreeModuleElement into a SparseFreeModuleElement?
# TODO: Call this class LieAlgebraElement_sparse?
# TODO: Inherit from LieAlgebraElement_generic?
class LieAlgebraElement(CombinatorialFreeModuleElement):
    """
    A Lie algebra element.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        """
        if not self._monomial_coefficients:
            return '0'
        return repr_lincomb(self.list())

    # Default implementation
    def _latex_monomial(self, m):
        """
        Return a `\LaTeX` representation of the monomial ``m``.

        EXAMPLES::
        """
        from sage.misc.latex import latex
        return latex(m)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::
        """
        if not self._monomial_coefficients:
            return '0'
        return repr_lincomb(self.list(), repr_monomial=self._latex_monomial, is_latex=True)

    def _mul_(self, y):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::
        """
        if self == 0 or y == 0:
            return self.parent().zero()
        # Otherwise we lift to the UEA
        return self.lift() * y.lift()

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of ``self`` in ``codomain`` under the map that sends
        the images of the generators of the parent of ``self`` to the
        tuple of elements of ``im_gens``.

        EXAMPLES::
        """
        s = codomain.zero()
        if not self: # If we are 0
            return s
        names = self.parent().variable_names()
        return codomain.sum(c * t._im_gens_(codomain, im_gens, names)
                            for t, c in self._monomial_coefficients.iteritems())

    # TODO: Move to the category/lift morphism
    def lift(self):
        """
        Lift ``self`` to the universal enveloping algebra.

        EXAMPLES::
        """
        UEA = self.parent().universal_enveloping_algebra()
        gen_dict = UEA.gens_dict()
        s = UEA.zero()
        if not self:
            return s
        for t, c in self._monomial_coefficients.iteritems():
            if isinstance(t, LieBracket):
                s += c * t.lift(gen_dict)
            else:
                s += c * gen_dict[t._name]
        return s

    def is_constant(self):
        """
        Check if ``self`` is a constant (i.e. zero).

        EXAMPLES::
        """
        return notself._monomial_coefficients

    def dict(self):
        """
        Return ``self`` as a dictionary mapping monomials to coefficients.

        EXAMPLES::
        """
        return copy(self._monomial_coefficients)

    def list(self):
        """
        Return ``self`` as a list of pairs ``(m, c)`` where ``m`` is a
        monomial and ``c`` is the coefficient.

        EXAMPLES::
        """
        return sorted(self._monomial_coefficients.items())

class LieAlgebraElementWrapper(ElementWrapper):
    """
    Wrap an element as a Lie algebra element.
    """
    def __eq__(self, rhs):
        """
        Check equality.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',3], representation='matrix')
            sage: L.bracket(L.e(2), L.e(1)) == -L.bracket(L.e(1), L.e(2))
            True
        """
        if not isinstance(rhs, LieAlgebraElementWrapper):
            return self.value == 0 and rhs == 0
        return self.parent() == rhs.parent() and self.value == rhs.value

    def __nonzero__(self):
        """
        Check non-zero.
        """
        return bool(self.value)

    def _add_(self, rhs):
        """
        Add ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(QQ, R)
            sage: x + y
            x + y
        """
        return self.__class__(self.parent(), self.value + rhs.value)

    def _sub_(self, rhs):
        """
        Subtract ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(QQ, R)
            sage: x - y
            x - y
        """
        return self.__class__(self.parent(), self.value - rhs.value)

    # This seems to work with 2 underscores and I don't understand why...
    def _mul_(self, x):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L.<x,y> = LieAlgebra(QQ, S)
            sage: x*y - y*x
            (2,3) - (1,3)
        """
        if self.value == 0 or x == 0:
            return self.parent().zero()
        # Otherwise we lift to the UEA
        return self.lift() * x.lift()

    def _acted_upon_(self, scalar, self_on_left=False):
        """
        Return the action of a scalar on ``self``.

        EXAMPLES::
        """
        # This was copied and IDK if it still applies:
        # With the current design, the coercion model does not have
        # enough information to detect apriori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        if hasattr( scalar, 'parent' ) and scalar.parent() != self.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self.base_ring().has_coerce_map_from(scalar.parent()):
                scalar = self.base_ring()( scalar )
            else:
                return None
        if self_on_left:
            return self.__class__(self.parent(), self.value * scalar)
        return self.__class__(self.parent(), scalar * self.value)

    def __neg__(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(QQ, R)
            sage: -x
            -x
        """
        return self.__class__(self.parent(), -self.value)

    def __getitem__(self, i):
        """
        Redirect the ``__getitem__()`` to the wrapped element.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',2], representation='matrix')
            sage: m = L.e(0)
            sage: m[0,0]
            0
            sage: m[0][1]
            1
        """
        return self.value.__getitem__(i)

