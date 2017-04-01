# -*- coding: utf-8 -*-
"""
Lie Algebra Elements

AUTHORS:

- Travis Scrimshaw (2013-05-04): Initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from functools import total_ordering

#from sage.misc.abstract_method import abstract_method
#from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.misc import repr_lincomb
from sage.combinat.free_module import CombinatorialFreeModule
#from sage.structure.element import ModuleElement, RingElement, coerce_binop
from sage.structure.element cimport have_same_parent, coercion_model
from sage.structure.element_wrapper cimport ElementWrapper
from sage.structure.sage_object cimport SageObject

# TODO: Inherit from IndexedFreeModuleElement and make cdef once #22632 is merged
# TODO: Do we want a dense version?
class LieAlgebraElement(CombinatorialFreeModule.Element):
    """
    A Lie algebra element.
    """
    # Need to bypass the coercion model
    def __mul__(self, y):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: y*x
            x*y - z
        """
        if self.is_zero() or y.is_zero():
            return self.parent().zero()
        if y in self.base_ring():
            return y * self
        # Otherwise we lift to the UEA
        return self.lift() * y

    #def _im_gens_(self, codomain, im_gens):
    #    """
    #    Return the image of ``self`` in ``codomain`` under the map that sends
    #    the images of the generators of the parent of ``self`` to the
    #    tuple of elements of ``im_gens``.
    #
    #    EXAMPLES::
    #    """
    #    s = codomain.zero()
    #    if not self: # If we are 0
    #        return s
    #    names = self.parent().variable_names()
    #    return codomain.sum(c * t._im_gens_(codomain, im_gens, names)
    #                        for t, c in self._monomial_coefficients.iteritems())

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

cdef class LieAlgebraElementWrapper(ElementWrapper):
    """
    Wrap an element as a Lie algebra element.

    TESTS:

    We check comparisons::

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

    Check inequality::

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
        sage: L.an_element()
        X + Y + Z
        sage: L.an_element() == 0
        False
        sage: L.an_element() != 0
        True

        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
        sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
        False
        sage: L.zero() == 0
        True
        sage: L.zero() != 0
        False
        sage: L.zero() >= 0
        True
        sage: L.zero() < 0
        False
    """

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

    def _ascii_art_(self):
        """
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: L = LieAlgebra(associative=s)
            sage: P = Partition([4,2,2,1])
            sage: x = L.basis()[P]
            sage: ascii_art(x)
            s
             ****
             **
             **
             *
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(self.value)

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: L = LieAlgebra(associative=s)
            sage: P = Partition([4,2,2,1])
            sage: x = L.basis()[P]
            sage: unicode_art(x)
            s
             ┌┬┬┬┐
             ├┼┼┴┘
             ├┼┤
             ├┼┘
             └┘
        """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art(self.value)

    def __nonzero__(self):
        """
        Return if ``self`` is non-zero.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: bool(L.zero())
            False
            sage: bool(x + y)
            True
        """
        return bool(self.value)

    cpdef _add_(self, right):
        """
        Add ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x + y
            x + y
        """
        return type(self)(self._parent, self.value + right.value)

    cpdef _sub_(self, right):
        """
        Subtract ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x - y
            x - y
        """
        return type(self)(self._parent, self.value - right.value)

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
            return self._parent.zero()
        if x in self.base_ring():
            return self._acted_upon_(x, True)
        # Otherwise we lift to the UEA
        return self.lift() * x

    def __div__(self, x):
        """
        Division by coefficients.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 3)
            sage: x = L.an_element(); x
            p1
            sage: x / 2
            1/2*p1
        """
        return self * (~x)

    cpdef _acted_upon_(self, scalar, bint self_on_left):
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
            return type(self)(self._parent, self.value * scalar)
        return type(self)(self._parent, scalar * self.value)

    def __neg__(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: -x
            -x
        """
        return type(self)(self._parent, -self.value)

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

# TODO: Also used for vectors, find a better name
cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    """
    Lie algebra element wrapper around a matrix.
    """
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

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    """
    An element of a Lie algebra given by structure coefficients.
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: x - 3/2 * y
            x - 3/2*y
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self._parent._print_options['scalar_mult'],
                            repr_monomial=self._parent._repr_generator,
                            strip_one=True)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: latex(elt)
            x - \frac{3}{2}y
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self._parent._print_options['scalar_mult'],
                            latex_scalar_mult=self._parent._print_options['latex_scalar_mult'],
                            repr_monomial=self._parent._latex_term,
                            is_latex=True, strip_one=True)

    cpdef bracket(self, right):
        """
        Return the Lie bracket ``[self, right]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}})
            sage: x.bracket(y)
            z
            sage: y.bracket(x)
            -z
            sage: (x + y - z).bracket(x - y + z)
            -2*y - 2*z
        """
        if not have_same_parent(self, right):
            self, right = coercion_model.canonical_coercion(self, right)
        return self._bracket_(right)

    # We need this method because the LieAlgebra.bracket method (from the
    #   category) calls this, where we are guaranteed to have the same parent.
    cpdef _bracket_(self, right):
        """
        Return the Lie bracket ``[self, right]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}})
            sage: x._bracket_(y)
            z
            sage: y._bracket_(x)
            -z
        """
        P = self._parent
        cdef dict s_coeff = P._s_coeff
        d = P.dimension()
        cdef list ret = [P.base_ring().zero()]*d
        cdef int i1, i2, i3
        cdef StructureCoefficientsElement rt = <StructureCoefficientsElement> right
        for i1 in range(d):
            c1 = self.value[i1]
            if not c1:
                continue
            for i2 in range(d):
                c2 = rt.value[i2]
                if not c2:
                    continue
                prod_c1_c2 = c1 * c2
                if (i1, i2) in s_coeff:
                    v = s_coeff[i1, i2]
                    for i3 in range(d):
                        ret[i3] += prod_c1_c2 * v[i3]
                elif (i2, i1) in s_coeff:
                    v = s_coeff[i2, i1]
                    for i3 in range(d):
                        ret[i3] -= prod_c1_c2 * v[i3]
        return type(self)(P, P._M(ret))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: list(elt)
            [('x', 1), ('y', -3/2)]
        """
        zero = self.base_ring().zero()
        I = self.parent()._indices
        cdef int i
        for i,v in enumerate(self.value):
            if v != zero:
                yield (I[i], v)

    cpdef to_vector(self):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = x + 3*y - z/2
            sage: a.to_vector()
            (1, 3, -1/2)
        """
        return self.value

    def lift(self):
        """
        Return the lift of ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: l = elt.lift(); l
            x - 3/2*y
            sage: l.parent()
            Noncommutative Multivariate Polynomial Ring in x, y
             over Rational Field, nc-relations: {y*x: x*y - x}
        """
        UEA = self.parent().universal_enveloping_algebra()
        gens = UEA.gens()
        return UEA.sum(c * gens[i] for i, c in self.value.iteritems())

    cpdef dict monomial_coefficients(self, bint copy=True):
        """
        Return the monomial coefficients of ``self`` as a dictionary.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = 2*x - 3/2*y + z
            sage: a.monomial_coefficients()
            {'x': 2, 'y': -3/2, 'z': 1}
            sage: a = 2*x - 3/2*z
            sage: a.monomial_coefficients()
            {'x': 2, 'z': -3/2}
        """
        I = self._parent._indices
        return {I[i]: v for i,v in self.value.iteritems()}

    def __getitem__(self, i):
        """
        Return the coefficient of the basis element indexed by ``i``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
            sage: elt = x - 3/2 * y
            sage: elt['y']
            -3/2
        """
        return self.value[self._parent._indices.index(i)]

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

