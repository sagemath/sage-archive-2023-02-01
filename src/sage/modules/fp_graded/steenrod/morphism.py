r"""
Homomorphisms of finitely presented modules over the Steenrod algebra

This class implements construction and basic manipulation of homomorphisms
between finitely presented graded modules over the mod `p`
Steenrod algebra.

For an overview of the API, see :doc:`module`.

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.

"""

#*****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu>
#             and          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import Hom

from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra_generic
from sage.modules.fp_graded.morphism import FPModuleMorphism
from .profile import enveloping_profile_elements


class FPA_ModuleMorphism(FPModuleMorphism):
    r"""
    Create a homomorphism between finitely presented graded modules over
    the mod `p` Steenrod algebra.

    INPUT:

    - ``parent`` -- A homspace object.

    - ``values`` -- A list of elements in the codomain.  Each element
      corresponds to a module generator in the domain.

    - ``check`` -- boolean (default: ``True``); if ``True``, check
      that the morphism is well-defined.

    OUTPUT: A module homomorphism defined by sending the generator with
    index `i` to the corresponding element in ``values``.

    .. NOTE:: Never use this constructor explicitly, but rather the parent's
        call method, or this class' __call__ method.  The reason for this
        is that the dynamic type of the element class changes as a
        consequence of the category system.

    TESTS:

        sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
        sage: # Trying to map the generators of a non-free module into a
        sage: # free module:
        sage: A = SteenrodAlgebra(2)
        sage: F.<a2,a3> = FPA_Module(A, [2,3])
        sage: Q.<x2,x3> = FPA_Module(A, [2,3], relations=[[Sq(6), Sq(5)]])
        sage: Hom(F, Q)((Sq(1)*x2, x3))
        Traceback (most recent call last):
         ...
        ValueError: ill-defined homomorphism: degrees do not match

    Trying to map the generators of a non-free module into a free module::

        sage: w = Hom(Q, F)((a2, a3))
        Traceback (most recent call last):
         ...
        ValueError: relation Sq(6)*x2 + Sq(5)*x3 is not sent to zero

    """
    def __init__(self, parent, values, check=True):
        r"""
        TESTS::

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M = FPA_Module(A2, [0], relations=[[Sq(1)]])
            sage: N = FPA_Module(A2, [0], relations=[[Sq(4)],[Sq(1)]])
            sage: f = Hom(M,N)([A2.Sq(3)*N.generator(0)])
            sage: TestSuite(f).run()
        """
        FPModuleMorphism.__init__(self, parent, values, check)


    def profile(self):
        r"""
        A finite profile over which this homomorphism can be defined.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: one = Hom(M,M).identity()
            sage: one.profile()
            (2, 1)
            sage: zero = Hom(M,M).zero()
            sage: zero.profile()
            (2, 1)
            sage: A_fin = SteenrodAlgebra(2, profile=(2,1))
            sage: M_fin = M.change_ring(A_fin)

        Change the ring of the module M::

            sage: M_fin.change_ring(A) is M
            True

        We can change rings to the finite sub-Hopf algebra defined by
        the profile we just computed::

            sage: one_fin = one.change_ring(A_fin)
            sage: one_fin.domain()
            Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]

        And if we change back to the full Steenrod algebra, we are back were
        we started::

            sage: one_fin.change_ring(A) == one
            True

        """
        def _flatten(f):
            return [coeffifient for value in f.values()\
                for coeffifient in value.dense_coefficient_list()]

        elements = _flatten(self.domain()._j) +\
            _flatten(self.codomain()._j) +\
            _flatten(self)

        elements = [a for a in elements if a not in (0, 1)]

        profile = enveloping_profile_elements(elements)

        # Avoid returning the zero profile because it triggers a corner case
        # in FPModule.resolution().
        #
        # XXX: Fix FPModule.resolution().
        #
        return (1,) if profile == (0,) else profile


    def is_injective(self, verbose=False):
        r"""
        Determine if this homomorphism is injective.

        INPUT:

        - ``verbose`` -- A boolean to control if log messages should be emitted.
          (optional, default: ``False``)

        OUTPUT: The boolean value ``True`` if this homomorphism has a trivial
        kernel, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = FPA_Module(A, [0], [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: g = Hom(S, M).zero()
            sage: g.is_injective()
            False
            sage: z = Hom(FPA_Module(A, []), M).zero()
            sage: z.is_injective()
            True
            sage: z.is_zero()
            True

        """
        algebra = self.base_ring()

        finite_algebra = SteenrodAlgebra_generic(algebra.prime(), profile=self.profile())

        return FPModuleMorphism.is_injective(
            self.change_ring(finite_algebra),
            verbose=verbose)


    def kernel_inclusion(self, top_dim=None, verbose=False):
        r"""
        The kernel of this homomorphism.

        INPUT:

        - ``verbose`` -- A boolean to control if log messages should be emitted.
          (optional, default: ``False``)

        OUTPUT: An injective homomorphism into the domain ``self`` which is
        onto the kernel of this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = FPA_Module(A, [0], [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: k = f.kernel_inclusion()
            sage: k == 0
            True

        Since k is both trivial and injective, its domain should
        be the zero module::

            sage: k.domain().is_trivial()
            True

            sage: g = Hom(S, M)([M([Sq(3),Sq(2)])])
            sage: h = g.kernel_inclusion()
            sage: h.is_identity()
            True
            sage: ker = h.domain();
            sage: ker is S
            True

        So `g` had to be trivial::

            sage: g.is_zero()
            True

        """
        return self._action(FPModuleMorphism.kernel_inclusion, verbose)


    def cokernel_projection(self, verbose=False):
        r"""
        Compute the map to the cokernel of ``self``.

        OUTPUT:

        The natural projection from the codomain of this homomorphism
        to its cokernel.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A1 = SteenrodAlgebra(2, profile=(2,1))
            sage: M = FPA_Module(A1, [0], [[Sq(2)]])
            sage: F = FPA_Module(A1, [0])

            sage: r = Hom(F, M)([A1.Sq(1)*M.generator(0)])
            sage: co = r.cokernel_projection(); co
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]
              To:   Finitely presented left module on 1 generator and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]
              Defn: g[0] |--> g[0]

            sage: co.domain().is_trivial()
            False
        """
        from .module import FPA_Module
        new_relations = ([x.dense_coefficient_list()
                          for x in self.codomain().relations()] +
                         [x.dense_coefficient_list() for x in self._values]

        coker = FPA_Module(self.base_ring(),
                    self.codomain().generator_degrees(),
                    relations=tuple(new_relations))

        projection = Hom(self.codomain(), coker)(coker.generators())

        return projection


    def image(self, verbose=False):
        r"""
        Compute the image of this homomorphism.

        INPUT:

        - ``verbose`` -- A boolean to control if log messages should be emitted.
          (optional, default: ``False``)

        OUTPUT: An injective homomorphism into the codomain of ``self`` which is
        onto the image of this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = FPA_Module(A, [0], [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: i = f.image(); i
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[1] |--> g[1]
            sage: i.codomain() is M
            True

        Lift the map ``f`` over the inclusion ``i``::

            sage: f_ = f.lift(i)
            sage: f_.is_injective()
            True
            sage: f_.is_surjective()
            True

            sage: g = Hom(S, M)([M([Sq(3),Sq(2)])])
            sage: j = g.image(); j
            Module morphism:
              From: Finitely presented left module on 0 generators and 0 relations over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 2 relations over mod 2 Steenrod algebra, milnor basis

        So ``g`` had to be trivial::

            sage: g.is_zero()
            True

        """
        return self._action(FPModuleMorphism.image, verbose)


    def _resolve_kernel(self, top_dim=None, verbose=False):
        r"""
        Resolve the kernel of this homomorphism by a free module.

        INPUT:

        - ``verbose`` -- A boolean to enable progress messages. (optional,
          default: ``False``)

        OUTPUT: A homomorphism `j: F \rightarrow D` where `D` is the domain of
        this homomorphism, `F` is free and such that `\ker(self) = \operatorname{im}(j)`.

        TESTS:

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: F = FPA_Module(A, [0,0])
            sage: L = FPA_Module(A, [0,0], [[Sq(3),Sq(0,1)], [0,Sq(2)]])
            sage: f = Hom(F, L)([L([Sq(2),0]), L([0, Sq(2)])])
            sage: f._resolve_kernel()
            Module morphism:
              From: Finitely presented left module on 3 generators and 0 relations over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 0 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0, 0] |--> g[0, 1]
                    g[3, 0] |--> Sq(0,1)*g[0, 0]
                    g[3, 1] |--> Sq(3)*g[0, 0]
        """
        return self._action(FPModuleMorphism._resolve_kernel, verbose)


    def _resolve_image(self, top_dim=None, verbose=False):
        r"""
        Resolve the image of this homomorphism by a free module.

        INPUT:

        - ``verbose`` -- A boolean to enable progress messages. (optional,
          default: ``False``)

        OUTPUT: A homomorphism `j: F \rightarrow C` where `C` is the codomain
        of this homomorphism, `F` is free, and
        `\operatorname{im}(self) = \operatorname{im}(j)`.

        TESTS:

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: F = FPA_Module(A, [0,0])
            sage: L = FPA_Module(A, [0,0], [[Sq(3),Sq(0,1)], [0,Sq(2)]])
            sage: f = Hom(F, L)([L([Sq(2),0]), L([0, Sq(2)])])
            sage: f._resolve_image()
            Module morphism:
              From: Finitely presented left module on 1 generator and 0 relations over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[2] |--> Sq(2)*g[0, 0]
        """
        return self._action(FPModuleMorphism._resolve_image, verbose)


    def _action(self, method, verbose=False):
        r"""
        Changes the ground ring to a finite algebra, acts by the given method
        and changes back into the original ground ring before returning.

        TESTS:

            sage: from sage.modules.fp_graded.steenrod.module import FPA_Module
            sage: from sage.modules.fp_graded.morphism import FPModuleMorphism
            sage: A = SteenrodAlgebra(2)
            sage: F = FPA_Module(A, [0])
            sage: L = FPA_Module(A, [0], [[Sq(3)]])
            sage: f = Hom(F, L)([L([Sq(2)])])
            sage: f._action(FPModuleMorphism._resolve_image, verbose=True)
            Computing using the profile:
            (2, 1)
            Resolving the image in the range of dimensions [0, 8]: 0 1 2 3 4 5 6 7 8.
            Module morphism:
              From: Finitely presented left module on 1 generator and 0 relations over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[2] |--> Sq(2)*g[0]
        """
        small_profile = self.profile()

        if verbose:
            print('Computing using the profile:')
            print(small_profile)

        algebra = self.base_ring()

        # Choose a finite sub Hopf-algebra of the original algebra.
        finite_algebra = SteenrodAlgebra_generic(algebra.prime(), profile=small_profile)

        # Perform the chosen action on the module after having changed rings
        # to the finite algebra.
        fp_result = method(
            self.change_ring(finite_algebra),
            verbose=verbose)

        # Change back to the original algebra and also from FPModule
        # to FPA_Module, and return the result.
        #
        # This is very clunky. Clean it up!
        f = fp_result.change_ring(self.base_ring())
        from .module import FPA_Module
        try:
            M = FPA_Module(f.domain()._j)
        except AttributeError: # f.domain() is free
            M = FPA_Module(f.domain())
        try:
            N = FPA_Module(f.codomain()._j)
        except AttributeError:
            N = FPA_Module(f.codomain())
        new_values = [N.linear_combination(zip(N.generators(),
                                               v.dense_coefficient_list()))
                      for v in f.values()]
        return Hom(M, N)(new_values)
