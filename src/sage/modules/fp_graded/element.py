r"""
Elements of finitely presented graded modules

This class implements construction and basic manipulation of elements of the
Sage parent :class:`sage.modules.fp_graded.module.FPModule`, which models
finitely presented modules over connected graded algebras.

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

#*****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement


class FPElement(IndexedFreeModuleElement):
    r"""
    A module element of a finitely presented graded module over
    a connected graded algebra.

    TESTS::

        sage: from sage.modules.fp_graded.module import FPModule
        sage: FPModule(SteenrodAlgebra(2), [0])([Sq(2)])
        <Sq(2)>
    """
    def lift_to_free(self):
        r"""
        A lift of this element to the free module F,
        if this element is in a quotient of F.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: x = M([Sq(1), 1])
            sage: x
            <Sq(1), 1>
            sage: x.parent()
            Finitely presented left module on 2 generators and 1 relation over mod 2 Steenrod algebra, milnor basis
            sage: x.lift_to_free()
            <Sq(1), 1>
            sage: x.lift_to_free().parent()
            Finitely presented free left module on 2 generators over mod 2 Steenrod algebra, milnor basis
        """
        C = self.parent().j.codomain()
        return C(self.coefficients())


    @cached_method
    def degree(self):
        r"""
        The degree of ``self``.

        OUTPUT:

        The integer degree of ``self`` or raise an error if the zero element.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: x = M.an_element(7)

            sage: x
            <Sq(0,0,1), Sq(3,1)>
            sage: x.degree()
            7

        The zero element has no degree::

            sage: (x-x).degree()
            Traceback (most recent call last):
            ...
            ValueError: the zero element does not have a well-defined degree

        TESTS::

            sage: N = FPModule(SteenrodAlgebra(2), [0], [[Sq(2)]])
            sage: y = Sq(2)*N.generator(0)
            sage: y == 0
            True
            sage: y.degree()
            Traceback (most recent call last):
            ...
            ValueError: the zero element does not have a well-defined degree
        """
        if self.is_zero():
            raise ValueError("the zero element does not have a well-defined degree")
        return self.lift_to_free().degree()


    def coefficients(self):
        """
        Return a list of all coefficients of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra()
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: x = M([Sq(1), 1])
            sage: x.coefficients()
            [Sq(1), 1]
            sage: y = Sq(2) * M.generator(1)
            sage: y.coefficients()
            [0, Sq(2)]
        """
        return [self[i] for i in sorted(self.parent().basis().keys())]


    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: [M.an_element(n) for n in range(1,10)]
            [<Sq(1), 1>,
             <Sq(2), Sq(1)>,
             <Sq(0,1), Sq(2)>,
             <Sq(1,1), Sq(3)>,
             <Sq(2,1), Sq(4)>,
             <Sq(0,2), Sq(5)>,
             <Sq(0,0,1), Sq(3,1)>,
             <Sq(1,0,1), Sq(1,2)>,
             <Sq(2,0,1), Sq(2,2)>]
        """
        return self.lift_to_free()._repr_()


    def _lmul_(self, a):
        r"""
        Act by left multiplication on this element by ``a``.

        INPUT:

        - ``a`` -- an element of the algebra the parent module is defined over

        OUTPUT:

        The module element `a \cdot x` where `x` is ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FPModule(A2, [0,3], [[Sq(2)*Sq(4), Sq(3)]])
            sage: A2.Sq(2)*M.generator(1)
            <0, Sq(2)>
            sage: A2.Sq(2)*(A2.Sq(1)*A2.Sq(2)*M.generator(0) + M.generator(1))
            <Sq(2,1), Sq(2)>

        TESTS::

            sage: elements = [M.an_element(n) for n in range(1,10)]
            sage: a = A2.Sq(3)
            sage: [a*x for x in elements]
            [<Sq(1,1), 0>,
             0,
             <Sq(3,1), Sq(3)>,
             <0, Sq(1,1)>,
             0,
             <Sq(3,2), Sq(3,1)>,
             <Sq(3,0,1), Sq(7)>,
             <Sq(1,1,1), Sq(5,1)>,
             <0, Sq(3,2)>]
        """
        return self.parent()(a*self.lift_to_free())


    def vector_presentation(self):
        r"""
        A coordinate vector representing ``self`` when it is non-zero.

        These are coordinates with respect to the basis chosen by
        :meth:`sage.modules.fp_graded.module.FPModule.basis_elements`.
        When the element is zero, it has no well defined degree, and this
        function returns ``None``.

        OUTPUT:

        A vector of elements in the ground field of the algebra for
        this module when this element is non-zero.  Otherwise, the
        value ``None``.

        .. SEEALSO::

            :meth:`sage.modules.fp_graded.module.FPModule.vector_presentation`
            :meth:`sage.modules.fp_graded.module.FPModule.basis_elements`
            :meth:`sage.modules.fp_graded.module.FPModule.element_from_coordinates`

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FPModule(A2, (0,1))

            sage: x = M.an_element(7)
            sage: v = x.vector_presentation(); v
            (1, 0, 0, 0, 0, 1, 0)
            sage: type(v)
            <class 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>

            sage: V = M.vector_presentation(7)
            sage: v in V
            True

            sage: M.element_from_coordinates(v, 7) == x
            True

        We can use the basis for the module elements in the degree of `x`,
        together with the coefficients `v` to recreate the element `x`::

            sage: basis = M.basis_elements(7)
            sage: x_ = sum( [c*b for (c,b) in zip(v, basis)] ); x_
            <Sq(0,0,1), Sq(3,1)>
            sage: x == x_
            True

        TESTS::

            sage: M.zero().vector_presentation() is None
            True
        """
        # We cannot represent the zero element since it does not have a degree,
        # and we therefore do not know which vector space it belongs to.
        #
        # In this case, we could return the integer value 0 since coercion would
        # place it inside any vector space.  However, this will not work for
        # homomorphisms, so we return None to be consistent.
        try:
            degree = self.lift_to_free().degree()
        except ValueError:
            return None

        F_n = self.parent().vector_presentation(degree)
        return F_n.quotient_map()(self.lift_to_free().vector_presentation())


    def __bool__(self):
        r"""
        Determine if this element is non-zero.

        OUTPUT:

        The boolean value ``True`` if this element is non-zero
        and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,2,4], [[Sq(4),Sq(2),0]])
            sage: M(0) != 0
            False
            sage: M((Sq(6), 0, Sq(2))) == 0
            False
            sage: a = M((Sq(1)*Sq(2)*Sq(1)*Sq(4), 0, 0))
            sage: b = M((0, Sq(2)*Sq(2)*Sq(2), 0))
            sage: a != 0
            True
            sage: bool(b)
            True
            sage: bool(a + b)
            False
        """
        pres = self.vector_presentation()
        if pres is None:
            return False
        return bool(pres)


    def __eq__(self, other):
        r"""
        True iff ``self`` and ``other`` are equal.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: x = M([Sq(1), 1])
            sage: x
            <Sq(1), 1>
            sage: x == x
            True
            sage: x == M.zero()
            False
            sage: x-x == M.zero()
            True
        """
        try:
            return (self - other).is_zero()
        except TypeError:
            return False


    def normalize(self):
        r"""
        A normalized form of ``self``.

        OUTPUT:

        An instance representing the same module element as ``self`` in
        normalized form.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,2,4], [[Sq(4),Sq(2),0]])

            sage: m = M((Sq(6), 0, Sq(2))); m
            <Sq(6), 0, Sq(2)>
            sage: m.normalize()
            <Sq(6), 0, Sq(2)>
            sage: m == m.normalize()
            True

            sage: n = M((Sq(4), Sq(2), 0)); n
            <Sq(4), Sq(2), 0>
            sage: n.normalize()
            0
            sage: n == n.normalize()
            True
        """
        if self.is_zero():
            return self.parent().zero()

        v = self.vector_presentation()
        return self.parent().element_from_coordinates(v, self.degree())

