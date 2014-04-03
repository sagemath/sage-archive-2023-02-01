r"""
Jordan Algebras

AUTHORS:

- Travis Scrimshaw (2014-04-02): initial version
"""

#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import AlgebraElement
from sage.categories.algebras import Algebras
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.misc.cachefunc import cached_method
#from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.all import QQ

class JordanAlgebra(Parent, UniqueRepresentation):
    r"""
    A (special) Jordan algebra `A^+` from an associative algebra `A`.

    A *Jordan algebra* is an algebra over a ring `R` whose multiplication
    satisfies the following axioms:

    - `xy = yx`, and
    - `(xy)(xx) = x(y(xx))` (the Jordan identity).

    These axiom imply that a Jordan algebra is power-associative and the
    following generalization of Jordan's identity: `(x^m y) x^n = x^m (yx^n)`
    for all `m, n \in \ZZ_{>0}`.

    Let `A` be an associative algebra over a ring `R` of characteristic
    not equal to 2. We construct a Jordan algebra `A^+` by considering
    the multiplication as

    .. MATH::

        x \circ y = \frac{xy + yx}{2}.

    Often the multiplication is written as `x \circ y` to avoid confusion
    with the product in the associative algebra `A`. We note that if `A` is
    commutative then this reduces to the usual multiplication in `A`.

    Jordan algebras constructed in this fashion, or their subalgebras,
    are called *special*. All other Jordan algebras are called *exceptional*.

    INPUT:

    - ``A`` -- an associative algebra

    EXAMPLES:

    We consider the base algebra `A` as the free algebra on 3 generators::

        sage: F.<x,y,z> = FreeAlgebra(QQ)
        sage: J = JordanAlgebra(F); J
        Jordan algebra of Free Algebra on 3 generators (x, y, z) over Rational Field
        sage: a,b,c = J.gens()
        sage: a*b
        1/2*x*y + 1/2*y*x
        sage: b*a
        1/2*x*y + 1/2*y*x

    Jordan algebras are typically non-associative::

        sage: (a*b)*c
        1/4*x*y*z + 1/4*y*x*z + 1/4*z*x*y + 1/4*z*y*x
        sage: a*(b*c)
        1/4*x*y*z + 1/4*x*z*y + 1/4*y*z*x + 1/4*z*y*x

    We check the Jordan identity::

        sage: (a*b)*(a*a) == a*(b*(a*a))
        True
        sage: x = a + c
        sage: y = b - 2*a
        sage: (x*y)*(x*x) == x*(y*(x*x))
        True

    REFERENCES:

    - :wikipedia:`Jordan_algebra`
    """
    def __init__(self, A, names=None):
        """
        Initialize ``self``.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: TestSuite(J).run()
        """
        self._A = A
        R = A.base_ring()
        if R.characteristic() == 2:
            raise ValueError("the base ring cannot have characteristic 2")
        # This isn't right, but we need non-assoc algebras. Therefore we need #10963.
        # Similarly, this algebra is only unital if ``A`` is unital
        if A in AlgebrasWithBasis(R):
            cat = (AlgebrasWithBasis(R), CommutativeAlgebras(R))
        elif A in Algebras(R):
            cat = CommutativeAlgebras(R)
        else:
            raise ValueError("A is not an algebra")

        Parent.__init__(self, base=R, names=names, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: JordanAlgebra(F)
            Jordan algebra of Free Algebra on 3 generators (x, y, z) over Rational Field
        """
        return "Jordan algebra of {}".format(self._A)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J(5)
            5
            sage: elt = J(x + 2*x*y); elt
            x + 2*x*y
            sage: elt.parent() is J
            True
        """
        if isinstance(x, JordanAlgebra.Element):
            if x.parent() is self:
                return x
            return self.element_class(self, self._A(x._x))

        if x in A:
            return self.element_class(self, x)

        return self.element_class(self, self._A(x))

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        .. WARNING::

            These do not generate the Jordan algebra `A^+` as an algebra
            since `A^+` is isomorphic to `A` as `R`-modules.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.gens()
            (x, y, z)
        """
        return tuple(map(lambda x: self.element_class(self, x), self._A.gens()))

    @cached_method
    def zero(self):
        """
        Return the element 0.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.zero()
            0
        """
        return self.element_class(self, self._A.zero())

    @cached_method
    def one(self):
        """
        Return the element 1 if it exists.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: J = JordanAlgebra(F)
            sage: J.one()
            1
        """
        return self.element_class(self, self._A.one())

    class Element(AlgebraElement):
        """
        An element of a Jordan algebra.
        """
        def __init__(self, parent, x):
            """
            Initialize ``self``.
            """
            self._x = x
            AlgebraElement.__init__(self, parent)

        def _repr_(self):
            """
            Return a string representation of ``self``.
            """
            return repr(self._x)

        def _latex_(self):
            """
            Return a latex representation of ``self``.
            """
            from sage.misc.latex import latex
            return latex(self._x)

        def __nonzero__(self):
            """
            Return if ``self`` is non-zero.
            """
            return self._x.__nonzero__()

        def __eq__(self, other):
            """
            Check equality.
            """
            if not isinstance(other, JordanAlgebra.Element):
                return False
            if other.parent() != self.parent():
                return False
            return self._x == other._x

        def __ne__(self, other):
            """
            Check inequality.
            """
            return not self.__eq__(other)

        def _add_(self, other):
            """
            Add ``self`` and ``other``.
            """
            return self.__class__(self.parent(), self._x + other._x)

        def _neg_(self):
            """
            Negate ``self``.
            """
            return self.__class__(self.parent(), -self._x)

        def _sub_(self, other):
            """
            Subtract ``self`` and ``other``.
            """
            return self.__class__(self.parent(), self._x - other._x)

        def _mul_(self, other):
            """
            Multiply ``self`` and ``other``.
            """
            x = self._x
            y = other._x
            # This is safer than dividing by 2
            return self.__class__(self.parent(), (x*y + y*x) * ~QQ(2))

        def _lmul_(self, other):
            """
            Multiply ``self`` and ``other`` by the left action.
            """
            return self.__class__(self.parent(), self._x * other)

        def _rmul_(self, other):
            """
            Multiply ``self`` and ``other`` by the right action.
            """
            return self.__class__(self.parent(), other * self._x)

