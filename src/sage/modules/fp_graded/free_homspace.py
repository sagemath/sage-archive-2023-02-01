r"""
The set of homomorphisms of finitely generated free graded left modules

For an overview of the free module API, see :doc:`free_module`.

EXAMPLES::

    sage: from sage.modules.fp_graded.free_module import FreeGradedModule
    sage: A = SteenrodAlgebra(2)
    sage: F1 = FreeGradedModule(A, (1,3))
    sage: F2 = FreeGradedModule(A, (2,3))
    sage: homset = Hom(F1, F2)
    sage: homset
    Set of Morphisms from Finitely presented free left module on 2 generators ...
    sage: homset([F2((Sq(1), 1)), F2((0, Sq(2)))])
    Module homomorphism of degree 2 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<Sq(1), 1>, <0, Sq(2)>]
    sage: TestSuite(homset).run()

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

#*****************************************************************************
#       Copyright (C) 2021 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import Homset
from sage.misc.cachefunc import cached_method


class FreeGradedModuleHomspace(Homset):
    """
    Homspace between two free graded modules.
    """
    # In the category framework, Elements of the class FPModule are of the
    # class FPElement, see
    # http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html#implementing-the-category-framework-for-the-elements
    from .free_morphism import FreeGradedModuleMorphism

    Element = FreeGradedModuleMorphism

    def _element_constructor_(self, values):
        r"""
        Construct an element of ``self``.

        This function is used internally by the call method when creating
        homomorphisms.

        INPUT:

        - ``values`` -- a tuple of values (i.e. elements of the
          codomain for this homset) corresponding bijectively to the generators
          of the domain of this homset, or the zero integer constant

        OUTPUT:

        An instance of the morphism class.  The returned morphism is defined
        by mapping the module generators in the domain to the given values.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FreeGradedModule(A2, (1,3))
            sage: L = FreeGradedModule(A2, (2,5))
            sage: H = Hom(F, L)

            sage: values = (A2.Sq(4)*L.generator(0), A2.Sq(3)*L.generator(1))
            sage: f = H(values); f
            Module homomorphism of degree 5 defined by sending the generators
              [<1, 0>, <0, 1>]
            to
              [<Sq(4), 0>, <0, Sq(3)>]

            sage: H(0)
            The trivial homomorphism
        """
        from .free_morphism import FreeGradedModuleMorphism
        if isinstance(values, FreeGradedModuleMorphism):
            return values
        elif values == 0 or all(v.is_zero() for v in values):
            return self.zero()
        else:
            return self.element_class(self, values)


    def _an_element_(self):
        r"""
        Return a morphism belonging to ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FreeGradedModule(A2, (1,3))
            sage: L = FreeGradedModule(A2, (2,3))
            sage: H = Hom(F, L)
            sage: H._an_element_()
            The trivial homomorphism
        """
        return self.zero()


    @cached_method
    def zero(self):
        r"""
        Return the trivial morphism of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FreeGradedModule(A2, (1,3))
            sage: L = FreeGradedModule(A2, (2,3))
            sage: H = Hom(F, L)
            sage: H.zero()
            The trivial homomorphism
        """
        return self.element_class(self, self.codomain().zero())


    def identity(self):
        r"""
        Return the identity morphism, if ``self`` is an endomorphism set.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: L = FreeGradedModule(A2, (2,3))
            sage: H = Hom(L, L)
            sage: H.identity()
            The identity homomorphism

        TESTS::

            sage: F = FreeGradedModule(A2, (1,3))
            sage: H = Hom(F, L)
            sage: H.identity()
            Traceback (most recent call last):
            ...
            TypeError: this homspace does not consist of endomorphisms
        """
        if self.is_endomorphism_set():
            return self.element_class(self, self.codomain().generators())
        else:
            raise TypeError('this homspace does not consist of endomorphisms')

