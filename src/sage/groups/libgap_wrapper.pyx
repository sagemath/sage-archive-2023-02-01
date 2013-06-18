"""
LibGAP-based Groups

This module provides helper class for wrapping GAP groups via
:mod:`~sage.libs.gap.libgap`. See :mod:`~sage.groups.free_group` for an
example how they are used.

The parent class keeps track of the libGAP element object, to use it
in your Python parent you have to derive both from the suitable group
parent and :class:`ParentLibGAP` ::

    sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
    sage: from sage.groups.group import Group
    sage: class FooElement(ElementLibGAP):
    ...       pass
    sage: class FooGroup(Group, ParentLibGAP):
    ...       Element = FooElement
    ...       def __init__(self):
    ...           lg = libgap(libgap.CyclicGroup(3))    # dummy
    ...           ParentLibGAP.__init__(self, lg)
    ...           Group.__init__(self)

Note how we call the constructor of both superclasses to initialize
``Group`` and ``ParentLibGAP`` separately. The parent class implements
its output via LibGAP::

    sage: FooGroup()
    <pc group of size 3 with 1 generators>
    sage: type(FooGroup().gap())
    <type 'sage.libs.gap.element.GapElement'>

The element class is a subclass of
:class:`~sage.structure.element.MultiplicativeGroupElement`. To use
it, you just inherit from :class:`ElementLibGAP` ::

    sage: element = FooGroup().an_element()
    sage: element
    f1

The element class implements group operations and printing via LibGAP::

    sage: element._repr_()
    'f1'
    sage: element * element
    f1^2

AUTHORS:

- Volker Braun
"""

##############################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.libs.gap.libgap import libgap
from sage.libs.gap.element cimport GapElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.structure.element cimport Element


class ParentLibGAP(object):
    """
    A class for parents to keep track of the GAP parent.

    INPUT:

    - ``libgap_parent`` -- the libgap element that is the parent in
      GAP.

    EXAMPLES::

        sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
        sage: from sage.groups.group import Group
        sage: class FooElement(ElementLibGAP):
        ...       pass
        sage: class FooGroup(Group, ParentLibGAP):
        ...       Element = FooElement
        ...       def __init__(self):
        ...           lg = libgap(libgap.CyclicGroup(3))    # dummy
        ...           ParentLibGAP.__init__(self, lg)
        ...           Group.__init__(self)
        sage: FooGroup()
        <pc group of size 3 with 1 generators>
    """

    def __init__(self, libgap_parent):
        """
        The Python constructor.

        TESTS::

            sage: G = FreeGroup(3)
            sage: TestSuite(G).run()
        """
        assert isinstance(libgap_parent, GapElement)
        self._libgap = libgap_parent

    def gap(self):
        """
        Returns the gap representation of self

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES::

            sage: G = FreeGroup(3);  G
            Free Group on generators {x0, x1, x2}
            sage: G.gap()
            <free group on the generators [ x0, x1, x2 ]>
            sage: G.gap().parent()
            C library interface to GAP
            sage: type(G.gap())
            <type 'sage.libs.gap.element.GapElement'>

        This can be useful, for example, to call GAP functions that
        are not wrapped in Sage::

            sage: G = FreeGroup(3)
            sage: H = G.gap()
            sage: H.DirectProduct(H)
            <fp group on the generators [ f1, f2, f3, f4, f5, f6 ]>
            sage: H.DirectProduct(H).RelatorsOfFpGroup()
            [ f1^-1*f4^-1*f1*f4, f1^-1*f5^-1*f1*f5, f1^-1*f6^-1*f1*f6, f2^-1*f4^-1*f2*f4,
              f2^-1*f5^-1*f2*f5, f2^-1*f6^-1*f2*f6, f3^-1*f4^-1*f3*f4, f3^-1*f5^-1*f3*f5,
              f3^-1*f6^-1*f3*f6 ]
        """
        return self._libgap

    _gap_ = gap

    @cached_method
    def _gap_gens(self):
        """
        Return the generators as a LibGAP object

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES:

            sage: G = FreeGroup(2)
            sage: G._gap_gens()
            [ x0, x1 ]
            sage: type(_)
            <type 'sage.libs.gap.element.GapElement_List'>
        """
        return self._libgap.GeneratorsOfGroup()

    @cached_method
    def ngens(self):
        """
        Return the number of generators of self.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: G = FreeGroup(2)
            sage: G.ngens()
            2

        TESTS::

            sage: type(G.ngens())
            <type 'sage.rings.integer.Integer'>
        """
        return self._gap_gens().Length().sage()

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        TESTS::

            sage: G.<a,b> =FreeGroup()
            sage: G     # indirect doctest
            Free Group on generators {a, b}
        """
        return self._libgap._repr_()

    def gen(self, i):
        """
        Return the `i`-th generator of self.

        .. warning::

            Indexing starts at `0` as usual in Sage/Python. Not as in
            GAP, where indexing starts at `1`.

        INPUT:

        - ``i`` -- integer between `0` (inclusive) and :meth:`ngens`
          (exclusive). The index of the generator.

        OUTPUT:

        The `i`-th generator of the group.

        EXAMPLES::

            sage: G = FreeGroup('a, b')
            sage: G.gen(0)
            a
            sage: G.gen(1)
            b
        """
        if not (0 <= i < self.ngens()):
            raise ValueError('i must be in range(ngens)')
        gap = self._gap_gens()[i]
        return self.element_class(gap, parent=self)

    @cached_method
    def gens(self):
        """
        Returns the generators of the group.

        EXAMPLES::

            sage: G = FreeGroup(2)
            sage: G.gens()
            (x0, x1)
            sage: H = FreeGroup('a, b, c')
            sage: H.gens()
            (a, b, c)

        :meth:`generators` is an alias for :meth:`gens` ::

            sage: G = FreeGroup('a, b')
            sage: G.generators()
            (a, b)
            sage: H = FreeGroup(3, 'x')
            sage: H.generators()
            (x0, x1, x2)
        """
        return tuple( self.gen(i) for i in range(self.ngens()) )

    generators = gens

    @cached_method
    def one(self):
        """
        Returns the identity element of self

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: G.one()
            1
            sage: G.one() == G([])
            True
            sage: G.one().Tietze()
            ()
        """
        return self.element_class(self.gap().Identity(), parent=self)

    def _an_element_(self):
        """
        Returns an element of self.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: G._an_element_()
            a*b
        """
        from sage.misc.all import prod
        return prod(self.gens())



cdef class ElementLibGAP(MultiplicativeGroupElement):
    """
    A class for LibGAP-based Sage group elements

    INPUT:

    - ``libgap_element`` -- the libgap element that is being wrapped.

    - ``parent`` -- the Sage parent.

    EXAMPLES::

        sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
        sage: from sage.groups.group import Group
        sage: class FooElement(ElementLibGAP):
        ...       pass
        sage: class FooGroup(Group, ParentLibGAP):
        ...       Element = FooElement
        ...       def __init__(self):
        ...           lg = libgap(libgap.CyclicGroup(3))    # dummy
        ...           ParentLibGAP.__init__(self, lg)
        ...           Group.__init__(self)
        sage: FooGroup()
        <pc group of size 3 with 1 generators>
        sage: FooGroup().gens()
        (f1,)
    """

    def __init__(self, libgap_element, parent):
        """
        The Python constructor

        TESTS::

            sage: G = FreeGroup(2)
            sage: g = G.an_element()
            sage: TestSuite(g).run()
        """
        MultiplicativeGroupElement.__init__(self, parent)
        assert isinstance(libgap_element, GapElement)
        assert isinstance(parent, ParentLibGAP)
        self._libgap = libgap_element

    cpdef GapElement gap(self):
        """
        Returns a LibGAP representation of the element

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: x
            a*b*a^-1*b^-1
            sage: xg = x.gap()
            sage: xg
            a*b*a^-1*b^-1
            sage: type(xg)
            <type 'sage.libs.gap.element.GapElement'>
        """
        return self._libgap

    _gap_ = gap

    def is_one(self):
        """
        Test whether the group element is the trivial element.

        OUTPUT:

        Boolean.

        EXAMPLES:

            sage: G.<a,b> = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: x.is_one()
            False
            sage: (x * ~x).is_one()
            True
        """
        return self == self.parent().one()

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: a._repr_()
            'a'
            sage: type(a)
            <class 'sage.groups.free_group.FreeGroup_class_with_category.element_class'>

            sage: x = G([1, 2, -1, -2])
            sage: x._repr_()
            'a*b*a^-1*b^-1'
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: y._repr_()
            'b^3*a*b^-3'

            sage: G.one()
            1
        """
        if self.is_one():
            return '1'
        else:
            return self._libgap._repr_()

    def _latex_(self):
        """
        Return a LaTeX representation

        OUTPUT:

        String. A valid LaTeX math command sequence.

        EXAMPLES::

            sage: F.<a,b,c> = FreeGroup()
            sage: f = F([1, 2, 2, -3, -1]) * c^15 * a^(-23)
            sage: f._latex_()
            'a\\cdot b^{2}\\cdot c^{-1}\\cdot a^{-1}\\cdot c^{15}\\cdot a^{-23}'

            sage: F = FreeGroup(3)
            sage: f = F([1, 2, 2, -3, -1]) * F.gen(2)^11 * F.gen(0)^(-12)
            sage: f._latex_()
            'x_{0}\\cdot x_{1}^{2}\\cdot x_{2}^{-1}\\cdot x_{0}^{-1}\\cdot x_{2}^{11}\\cdot x_{0}^{-12}'

            sage: F.<a,b,c> = FreeGroup()
            sage: G = F /  (F([1, 2, 1, -3, 2, -1]), F([2, -1]))
            sage: f = G([1, 2, 2, -3, -1]) * G.gen(2)^15 * G.gen(0)^(-23)
            sage: f._latex_()
            'a\\cdot b^{2}\\cdot c^{-1}\\cdot a^{-1}\\cdot c^{15}\\cdot a^{-23}'

            sage: F = FreeGroup(4)
            sage: G = F.quotient((F([1, 2, 4, -3, 2, -1]), F([2, -1])))
            sage: f = G([1, 2, 2, -3, -1]) * G.gen(3)^11 * G.gen(0)^(-12)
            sage: f._latex_()
            'x_{0}\\cdot x_{1}^{2}\\cdot x_{2}^{-1}\\cdot x_{0}^{-1}\\cdot x_{3}^{11}\\cdot x_{0}^{-12}'
        """
        import re
        s = self._repr_()
        s = re.sub('([a-z]|[A-Z])([0-9]+)', '\g<1>_{\g<2>}', s)
        s = re.sub('(\^)(-)([0-9]+)', '\g<1>{\g<2>\g<3>}', s)
        s = re.sub('(\^)([0-9]+)', '\g<1>{\g<2>}', s)
        s = s.replace('*', '\cdot ')
        return s

    cpdef MonoidElement _mul_(left, MonoidElement right):
        """
        Multiplication of group elements

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x*y    # indirect doctest
            a*b*a^-1*b^2*a*b^-3
            sage: y*x    # indirect doctest
            b^3*a*b^-3*a*b*a^-1*b^-1
            sage: x*y == x._mul_(y)
            True
            sage: y*x == y._mul_(x)
            True
        """
        P = left.parent()
        return P.element_class(left.gap() * right.gap(), parent=P)

    cdef int _cmp_c_impl(left, Element right):
        """
        This method implements comparison.

        TESTS::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x == x*y*y^(-1)     # indirect doctest
            True
            sage: cmp(x,y)
            -1
            sage: x < y
            True
        """
        return cmp((<ElementLibGAP>left)._libgap,
                   (<ElementLibGAP>right)._libgap)

    def __richcmp__(left, right, int op):
        """
        Boilerplate for Cython elements

        See :mod:`~sage.structure.element` for details.
        """
        return (<Element>left)._richcmp(right, op)

    def __cmp__(left, right):
        """
        Boilerplate for Cython elements

        See :mod:`~sage.structure.element` for details.
        """
        return (<Element>left)._cmp(right)

    cpdef MultiplicativeGroupElement _div_(left, MultiplicativeGroupElement right):
        """
        Division of group elements.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x/y # indirect doctest
            a*b*a^-1*b^2*a^-1*b^-3
            sage: y/x # indirect doctest
            b^3*a*b^-2*a*b^-1*a^-1
            sage: x/y == x.__div__(y)
            True
            sage: x/y == y.__div__(x)
            False
        """
        P = left.parent()
        return P.element_class(left.gap() / right.gap(), parent=P)

    def __pow__(self, n, dummy):
        """
        Implement exponentiation.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: y^(2) # indirect doctest
            b^3*a^2*b^-3
            sage: x^(-3) # indirect doctest
            b*a*b^-1*a^-1*b*a*b^-1*a^-1*b*a*b^-1*a^-1
            sage: y^3  ==  y.__pow__(3)
            True
        """
        if n not in IntegerRing():
            raise TypeError("exponent must be an integer")
        P = self.parent()
        return P.element_class(self.gap().__pow__(n), parent=P)

    def __invert__(self):
        """
        Return the inverse of self.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x.__invert__()
            b*a*b^-1*a^-1
            sage: y.__invert__()
            b^3*a^-1*b^-3
            sage: ~x
            b*a*b^-1*a^-1
            sage: x.inverse()
            b*a*b^-1*a^-1
        """
        P = self.parent()
        return P.element_class(self.gap().Inverse(), parent=P)

    inverse = __invert__


