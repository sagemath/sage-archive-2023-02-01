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

#from sage.misc.abstract_method import abstract_method
#from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
#from sage.misc.misc import repr_lincomb
from copy import copy
#from functools import total_ordering
#from sage.structure.element import ModuleElement, RingElement, coerce_binop
#from sage.structure.sage_object import SageObject
from sage.combinat.free_module import CombinatorialFreeModuleElement
from sage.structure.element_wrapper import ElementWrapper

# TODO: Refactor this with Element wrapper to a generic SageObject wrapper?
@total_ordering
class LieGenerator(SageObject):
    """
    A wrapper around an object so it can ducktype with and do
    comparison operations with :class:`LieBracket`.
    """
    def __init__(self, name):
        """
        Initalize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: x = LieGenerator('x')
            sage: TestSuite(x).run()
        """
        self._name = name

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: LieGenerator('x')
            x
        """
        return self._name

    _latex_ = _repr_

    def __eq__(self, rhs):
        """
        Compare equals.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: x == y
            False
            sage: z = LieGenerator('x')
            sage: x == z
            True
        """
        return isinstance(rhs, LieGenerator) and self._name == rhs._name

    def __lt__(self, rhs):
        """
        Compare less than.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: x < y
            True
            sage: y < x
            False
            sage: z = LieBracket(x, y)
            sage: x < z
            True
        """
        if isinstance(rhs, LieGenerator):
            return self._name < rhs._name
        if isinstance(rhs, LieBracket):
            return not rhs.__lt__(self) # Clearly self != rhs
        return False

    def _im_gens_(self, codomain, im_gens, names):
        """
        Return the image of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: im = H(Lyn.lie_algebra_generators()['x']); im # indirect doctest
            x
            sage: im.parent() is H
            True
        """
        x = im_gens[names.index(self._name)]
        return im_gens[names.index(self._name)]

    def to_word(self):
        """
        Return ``self`` as a word in the variable names.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
            sage: x = LieGenerator('x')
            sage: x.to_word()
            ['x']
        """
        return [self._name]

@total_ordering
class LieBracket(SageObject): # Does this need to be SageObject?
    """
    A Lie bracket. This is the building blocks for Lie algebra elements.
    """
    def __init__(self, l, r):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: z = LieBracket(x, y)
            sage: TestSuite(z).run()
        """
        self._left = l
        self._right = r

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: LieBracket(x, y)
            [x, y]
        """
        return "[{!s}, {!s}]".format(self._left, self._right)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: z = LieBracket(x, y)
            sage: latex(z)
            \left[ x , y \right]
        """
        from sage.misc.latex import latex
        return "\\left[" + latex(self._left) + "," + latex(self._right) + "\\right]"

    def __getitem__(self, i):
        r"""
        Return the `i`-th item of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: z = LieBracket(x, y)
            sage: z[0]
            x
            sage: z[1] is y
            True
            sage: z[2]
            Traceback (most recent call last):
            ...
            IndexError: must be either 0 or 1
        """
        if i == 0:
            return self._left
        if i == 1:
            return self._right
        raise IndexError("must be either 0 or 1")

    def __eq__(self, rhs):
        """
        Check equality.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: b = LieBracket(x, y)
            sage: c = LieBracket(y, x)
            sage: b == c
            False
            sage: b == x
            False
            sage: a = LieBracket(x, y)
            sage: a == b
            True
            sage: a == [x, y]
            True
        """
        if isinstance(rhs, list):
            if len(rhs) != 2:
                return False
            return self._left == rhs[0] and self._right == rhs[1]

        if not isinstance(rhs, LieBracket):
            return False
        return self._left == rhs._left and self._right == rhs._right

    def __lt__(self, rhs):
        """
        Check less than.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: z = LieGenerator('z')
            sage: b = LieBracket(x, y)
            sage: b < x
            False
            sage: b < z
            False
            sage: c = LieBracket(x, z)
            sage: b < c
            True
        """
        if not isinstance(rhs, LieBracket):
            return False
        if self._left < rhs._left:
            return True
        elif self._left == rhs._left:
            return self._right < rhs._right
        return False

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: z = LieGenerator('z')
            sage: b = LieBracket(x, y)
            sage: hash(b) == hash(b)
            True
        """
        return hash((self._left, self._right))

    def _im_gens_(self, codomain, im_gens, names):
        """
        Return the image of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: x,y,z = Lyn.gens()
            sage: im = H(Lyn([z, [[x, y], x]])); im # indirect doctest
            -[z, [x, [x, y]]]
            sage: im.parent() is H
            True
        """
        return codomain.bracket(self._left._im_gens_(codomain, im_gens, names),
                                self._right._im_gens_(codomain, im_gens, names))

    def lift(self, UEA_gens_dict):
        """
        Lift ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: x,y,z = Lyn.gens()
            sage: a = Lyn([z, [[x, y], x]]); a
            [x, [x, [y, z]]] + [x, [[x, z], y]] - [[x, y], [x, z]]
            sage: a.lift() # indirect doctest
            x^2*y*z - 2*x*y*x*z + y*x^2*z - z*x^2*y + 2*z*x*y*x - z*y*x^2
        """
        if isinstance(self._left, LieBracket):
            l = self._left.lift(UEA_gens_dict)
        else:
            l = UEA_gens_dict[self._left._name]

        if isinstance(self._right, LieBracket):
            r = self._right.lift(UEA_gens_dict)
        else:
            r = UEA_gens_dict[self._right._name]

        return l*r - r*l

    def to_word(self):
        """
        Return ``self`` as a word expressed in the variable names.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: b = LieBracket(x, y)
            sage: c = LieBracket(b, x)
            sage: c.to_word()
            ['x', 'y', 'x']
        """
        return self._left.to_word() + self._right.to_word()

class GradedLieBracket(LieBracket):
    """
    A Lie bracket in a graded Lie algebra.
    """
    def __init__(self, l, r, grade):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, GradedLieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: b = GradedLieBracket(x, y, 2)
            sage: TestSuite(b).run()
        """
        self._grade = grade
        LieBracket.__init__(self, l, r)

    def __lt__(self, rhs):
        """
        Check less than.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, GradedLieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: z = LieGenerator('z')
            sage: b = GradedLieBracket(x, y, 2)
            sage: b < x
            False
            sage: b < z
            False
            sage: c = GradedLieBracket(x, z, 2)
            sage: b < c
            True
            sage: c = GradedLieBracket(x, z, 1)
            sage: b < c
            False
        """
        if isinstance(rhs, GradedLieBracket) and self._grade != rhs._grade:
            return self._grade < rhs._grade
        if isinstance(rhs, LieGenerator):
            return False # Lie generators have grade 1 and our grade > 1
        return LieBracket.__lt__(self, rhs)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, GradedLieBracket
            sage: x = LieGenerator('x')
            sage: y = LieGenerator('y')
            sage: z = LieGenerator('z')
            sage: b = GradedLieBracket(x, y, 2)
            sage: hash(b) == hash(b)
            True
        """
        return hash((self._grade, self._left, self._right))

# TODO: Have the other classes inherit from this?
# TODO: Should this be a mixin class (or moved to the category)?
#class LieAlgebraElement_generic(ModuleElement):
#    """
#    Generic methods for all Lie algebra elements.
#    """
#    def __mul__(self, other):
#        """
#        If we are multiplying two non-zero elements, automatically
#        lift up to the universal enveloping algebra.
#
#        EXAMPLES::
#        """
#        if self == 0 or other == 0:
#            return self.parent().zero()
#        # Otherwise we lift to the UEA
#        return self.lift() * other

# TODO: Factor out parts of CombinatorialFreeModuleElement into a SparseFreeModuleElement?
# TODO: Do we want a dense version?
class LieAlgebraElement(CombinatorialFreeModuleElement):
    """
    A Lie algebra element.
    """
    # Need to bypass the coercion model
    def __mul__(self, y):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: y*x
            x*y - z
        """
        if self == 0 or y == 0:
            return self.parent().zero()
        if y in self.base_ring():
            return y * self
        # Otherwise we lift to the UEA
        return self.lift() * y

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of ``self`` in ``codomain`` under the map that sends
        the images of the generators of the parent of ``self`` to the
        tuple of elements of ``im_gens``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: x,y,z = Lyn.gens()
            sage: im = H(Lyn([z, [[x, y], x]])); im # indirect doctest
            -[z, [x, [x, y]]]
            sage: im.parent() is H
            True
        """
        s = codomain.zero()
        if not self: # If we are 0
            return s
        names = self.parent().variable_names()
        return codomain.sum(c * t._im_gens_(codomain, im_gens, names)
                            for t, c in self._monomial_coefficients.iteritems())

    # TODO: Move to the category/lift morphism?
    def lift(self):
        """
        Lift ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: x.lift().parent() == L.universal_enveloping_algebra()
            True
        """
        UEA = self.parent().universal_enveloping_algebra()
        gen_dict = UEA.gens_dict()
        s = UEA.zero()
        if not self:
            return s
        for t, c in self._monomial_coefficients.iteritems():
            s += c * gen_dict[t]
        return s

    def is_constant(self):
        """
        Check if ``self`` is a constant (i.e. zero).

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: a = x + y
            sage: a.is_constant()
            False
            sage: L.zero().is_constant()
            True
        """
        return not self._monomial_coefficients

    def dict(self):
        """
        Return ``self`` as a dictionary mapping monomials to coefficients.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: a = 3*x - 1/2*z
            sage: a.dict()
            {'x': 3, 'z': -1/2}
        """
        return copy(self._monomial_coefficients)

    def list(self):
        """
        Return ``self`` as a list of pairs ``(m, c)`` where ``m`` is a
        monomial and ``c`` is the coefficient.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: a = 3*x - 1/2*z
            sage: a.list()
            [('x', 3), ('z', -1/2)]
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

            sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
            sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
            True

        The next doctests show similar behavior, although on elements of
        other classes::

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 3)
            sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
            True

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
            sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
            True
        """
        if not isinstance(rhs, LieAlgebraElementWrapper):
            return self.value == 0 and rhs == 0
        return self.parent() == rhs.parent() and self.value == rhs.value

    def __ne__(self, rhs):
        """
        Check non-equality.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
            sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
            False
            sage: L.zero() == 0
            True
            sage: L.zero() != 0
            False

        The next doctests show similar behavior, although on elements of
        other classes::

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 3)
            sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
            False
            sage: L.zero() == 0
            True
            sage: L.zero() != 0
            False

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
            sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
            False
            sage: L.zero() == 0
            True
            sage: L.zero() != 0
            False
        """
        return not (self == rhs)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x + y
            x + y
        """
        return repr(self.value)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x')
            sage: L.<x0,x1,x2> = LieAlgebra(associative=R.gens())
            sage: latex(x0 + x1)
            x_{0} + x_{1}
        """
        from sage.misc.latex import latex
        return latex(self.value)

    def _add_(self, rhs):
        """
        Add ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x + y
            x + y
        """
        return self.__class__(self.parent(), self.value + rhs.value)

    def _sub_(self, rhs):
        """
        Subtract ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x - y
            x - y
        """
        return self.__class__(self.parent(), self.value - rhs.value)

    # We need to bypass the coercion framework
    # We let the universal enveloping algebra handle the rest if both
    #   arguments are non-zero
    def __mul__(self, x):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        .. TODO::

            Write tests for this method once :trac:`16822` is
            implemented.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L.<x,y> = LieAlgebra(associative=S.gens())
            sage: u = x*3; u
            3*(1,2,3)
            sage: parent(u) == L
            True
            sage: u = x*(3/2); u
            3/2*(1,2,3)
            sage: parent(u) == L
            True
            sage: elt = x*y - y*x; elt  # not tested: needs #16822
            sage: S(elt)  # not tested: needs #16822
            (2,3) - (1,3)
        """
        if self.value == 0 or x == 0:
            return self.parent().zero()
        if x in self.base_ring():
            return x * self
        # Otherwise we lift to the UEA
        return self.lift() * x

    def __div__(self, x, self_on_left=False ):
        """
        Division by coefficients.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 3)
            sage: x = L.an_element(); x
            p1 + p2 + p3 + q1 + q2 + q3 + z
            sage: x / 2
            1/2*p1 + 1/2*p2 + 1/2*p3 + 1/2*q1 + 1/2*q2 + 1/2*q3 + 1/2*z
        """
        if self_on_left:
            return self * (~x)
        return (~x) * self

    def _acted_upon_(self, scalar, self_on_left=False):
        """
        Return the action of a scalar on ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: 3*x
            3*x
            sage: parent(3*x) == parent(x)
            True
            sage: x / 2
            1/2*x
            sage: y * (1/2)
            1/2*y
            sage: y * 1/2
            1/2*y
            sage: 1/2 * y
            1/2*y
            sage: QQ(1/2) * y
            1/2*y
        """
        # This was copied and IDK if it still applies (TCS):
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
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: -x
            -x
        """
        return self.__class__(self.parent(), -self.value)

    def __getitem__(self, i):
        """
        Redirect the ``__getitem__()`` to the wrapped element.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
            sage: m = L.gen(0)
            sage: m[0,0]
            0
            sage: m[0][1]
            1
        """
        return self.value.__getitem__(i)

class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    def __init__(self, parent, value):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: z = L.z()
            sage: z.value.is_immutable()
            True
        """
        value.set_immutable() # Make the matrix immutable for hashing
        LieAlgebraElementWrapper.__init__(self, parent, value)

