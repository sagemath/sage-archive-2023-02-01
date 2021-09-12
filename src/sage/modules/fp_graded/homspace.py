r"""
The set of homomorphisms of finitely presented graded modules

This class implements methods for construction and basic
manipulation of homsets of finitely presented graded modules over a connected
graded `k`-algebra, where `k` is a field.

.. NOTE:: This class is intended for private use by
    :class:`sage.modules.fp_steenrod.fpa_homspace.FPA_ModuleHomspace`.

TESTS:

    sage: from sage.modules.fp_graded.module import FP_Module
    sage: from sage.misc.sage_unittest import TestSuite
    sage: A = SteenrodAlgebra(2, profile=(3,2,1))
    sage: F = FP_Module(A, [1,3])
    sage: L = FP_Module(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])
    sage: homset = Hom(F, L); homset
    Set of Morphisms from Finitely presented left module on 2 generators ...
    sage: homset.an_element()
    Module homomorphism of degree 0 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<0, 0>, <Sq(1), 0>]
    sage: homset([L((A.Sq(1), 1)), L((0, A.Sq(2)))])
    Module homomorphism of degree 2 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<Sq(1), 1>, <0, Sq(2)>]
    sage: Hom(F, L) ([L((A.Sq(1), 1)), L((0, A.Sq(2)))]).kernel()
    Module homomorphism of degree 0 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      (<0, 1>, <Sq(0,1), 0>)

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

from __future__ import absolute_import

from sage.categories.homset import Homset

from sage.categories.homset import Hom

def is_FP_ModuleHomspace(x):
    r"""
    Check if the given object is of type FP_ModuleHomspace.

    OUTPUT: A boolean which is True if ``x`` is of type FP_ModuleHomspace.

    EXAMPLES::

        sage: from sage.modules.fp_graded.module import FP_Module
        sage: from sage.modules.fp_graded.homspace import is_FP_ModuleHomspace
        sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
        sage: F = FP_Module(A2, [1,3])
        sage: L = FP_Module(A2, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])
        sage: is_FP_ModuleHomspace(Hom(F, L))
        True
        sage: is_FP_ModuleHomspace(0)
        False
    """
    return isinstance(x, FP_ModuleHomspace)


class FP_ModuleHomspace(Homset):
    # FP_ModuleMorphism contains reference to is_FP_ModuleHomspace, so this import
    # statement must not appear before that function.
    from .morphism import FP_ModuleMorphism

    # In the category framework, Elements of the class FP_ModuleHomspace are of the
    # class FP_ModuleMorphism, see
    # http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html#implementing-the-category-framework-for-the-elements
    Element = FP_ModuleMorphism

    def _element_constructor_(self, values):
        r"""
        Constructs a morphism contained in this homset.

        This function is not part of the public API, but is used by :meth:Hom
        method to create morphisms.

        INPUT:

        - ``values`` -- An iterable of FP_Elements of the codomain.

        OUTPUT: A module homomorphism in this homspace sending the generators
        of the domain module to the given values.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FP_Module(A2, [1,3])
            sage: L = FP_Module(A2, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: homset = Hom(F, L)
            sage: v1 = L([A2.Sq(1), 1])
            sage: v2 = L([0, A2.Sq(2)])
            sage: f = homset._element_constructor_([v1, v2])

        Rather than calling ``_element_constructor_`` explicitly, one
        can call it implicitly, using the ``call`` syntax::

            sage: f = homset([v1, v2]); f
              Module homomorphism of degree 2 defined by sending the generators
                [<1, 0>, <0, 1>]
              to
                [<Sq(1), 1>, <0, Sq(2)>]

        One can construct a homomorphism from another homomorhism::

            sage: g = homset(f)
            sage: f == g
            True

        And there is a convenient way of making the trivial homomorphism::

            sage: z = homset(0); z
            The trivial homomorphism.
        """
        if isinstance(values, self.element_class):
            return values
        elif values == 0:
            return self.zero()
        else:
            return self.element_class(self, values)


    def an_element(self, n=0):
        r"""
        Create a homomorphism belonging to this homset.

        INPUT:

        - ``n`` -- an integer degree.  (optional, default: 0)

        OUTPUT: A module homomorphism of degree ``n``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FP_Module
            sage: A = SteenrodAlgebra(2)
            sage: HZ = FP_Module(A, [0], relations=[[Sq(1)]])

            sage: Hom(HZ, HZ).an_element(3)
            Module homomorphism of degree 3 defined by sending the generators
              [<1>]
            to
              [<Sq(0,1)>]

        TESTS:

            sage: K = FP_Module(A, [0, 0], [[Sq(2), 0]]) # Using a zero coefficient in the relations.
            sage: Hom(K, K).an_element(4)
            Module homomorphism of degree 4 defined by sending the generators
              [<1, 0>, <0, 1>]
            to
              [<0, 0>, <Sq(4), 0>]

            sage: K = FP_Module(A, [0, 0], [[Sq(2), 0], [0,0], [Sq(4), Sq(2)*Sq(2)]])
            sage: Hom(K, K).an_element(n=3)
            Module homomorphism of degree 3 defined by sending the generators
              [<1, 0>, <0, 1>]
            to
              [<0, 0>, <Sq(0,1), 0>]
        """
        return self._basis_elements(n, basis=False)


    def basis_elements(self, n):
        r"""
        Compute a basis for the vector space of degree ``n`` morphisms.

        INPUT:

        - ``n`` -- an integer degree.

        OUTPUT: A basis for the set of all module homomorphisms of degree ``n``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FP_Module
            sage: A = SteenrodAlgebra(2)
            sage: Hko = FP_Module(A, [0], relations=[[Sq(2)], [Sq(1)]])

            sage: Hom(Hko, Hko).basis_elements(21)
            [Module homomorphism of degree 21 defined by sending the generators
               [<1>]
             to
               [<Sq(0,0,3) + Sq(0,2,0,1)>],
             Module homomorphism of degree 21 defined by sending the generators
               [<1>]
             to
               [<Sq(8,2,1)>]]
        """
        return self._basis_elements(n, basis=True)


    def zero(self):
        r"""
        Create the trivial homomorphism in this homset.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FP_Module(A2, [1,3])
            sage: L = FP_Module(A2, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: z = Hom(F, L).zero(); z
            The trivial homomorphism.

            sage: z(F.an_element(5))
            <0, 0>

            sage: z(F.an_element(23))
            <0, 0>
        """
        return self.element_class(self, [self.codomain().zero() for g in self.domain().generator_degrees()])


    def identity(self):
        r"""
        Create the identity homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: L = FP_Module(A2, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: id = Hom(L, L).identity(); id
            The identity homomorphism.

            sage: e = L.an_element(5)
            sage: e == id(e)
            True

        It is an error to call this function when the homset is not a
        set of endomorphisms::

            sage: F = FP_Module(A2, [1,3])
            sage: Hom(F,L).identity()
            Traceback (most recent call last):
            ...
            TypeError: this homspace does not consist of endomorphisms
        """
        if self.is_endomorphism_set():
            return self.element_class(self, self.codomain().generators())
        else:
            raise TypeError('this homspace does not consist of endomorphisms')


    def _basis_elements(self, n, basis):
        r"""
        Compute a basis for the vector space of degree ``n`` homomorphisms.

        This function is private and used by :meth:`basis_elements` and
        :meth:`an_element`.

        INPUT:

        - ``n`` -- an integer degree.
        - ``basis`` -- boolean to decide if a basis should be returned, or just
          a single homomorphism.

        OUTPUT: A basis for the set of all module homomorphisms of degree ``n``
        if ``basis`` is True.  Otherwise a single element is returned.  In the
        latter case, this homomorphism is non-trivial if the vector space of all
        homomorphisms is non-trivial.

        TESTS:

            sage: from sage.modules.fp_graded.module import FP_Module
            sage: A = SteenrodAlgebra(2)
            sage: Hko = FP_Module(A, [0], relations=[[Sq(2)], [Sq(1)]])
            sage: Hom(Hko, Hko)._basis_elements(21, basis=True)
            [Module homomorphism of degree 21 defined by sending the generators
               [<1>]
             to
               [<Sq(0,0,3) + Sq(0,2,0,1)>],
             Module homomorphism of degree 21 defined by sending the generators
               [<1>]
             to
               [<Sq(8,2,1)>]]

            sage: Hom(Hko, Hko)._basis_elements(21, basis=False)
            Module homomorphism of degree 21 defined by sending the generators
              [<1>]
            to
              [<Sq(0,0,3) + Sq(0,2,0,1)>]

            sage: F = FP_Module(A, [0])
            sage: Hom(F, Hko)._basis_elements(21, basis=False)
            Module homomorphism of degree 21 defined by sending the generators
              [<1>]
            to
              [<Sq(0,2,0,1)>]

            sage: Hom(F, Hko)._basis_elements(21, basis=False)
            Module homomorphism of degree 21 defined by sending the generators
              [<1>]
            to
              [<Sq(0,2,0,1)>]

            Hom(FPA_Module([0], A, [[Sq(1)]]), FPA_Module([-2], A, [[Sq(1)]])).an_element(0)
            The trivial homomorphism.

        Test corner cases involving trivial modules:

            sage: F = FP_Module(A, [0]) # A module without relations.
            sage: Z0 = FP_Module(A, []) # A trivial module.
            sage: Z1 = FP_Module(A, [0], [[1]]) # A trivial module with a redundant generator and relation.

            Hom(FPA_Module([-1], A), F)._basis_elements(0, basis=True)
            []
            Hom(FPA_Module([-1], A), F)._basis_elements(0, basis=False)
            The trivial homomorphism.

            sage: from itertools import product
            sage: for D,C in product([(F, 'Free'), (Hko, 'Hko'), (Z0, 'Trivial'), (Z1, 'Trivial with redundant generator')], repeat=2):
            ....:     print('Hom(%s, %s):' % (D[1], C[1]))
            ....:     print('  basis==False:\n  %s' % Hom(D[0], C[0])._basis_elements(n=7, basis=False))
            ....:     print('  basis==True:\n  %s' % Hom(D[0], C[0])._basis_elements(n=7, basis=True))
            Hom(Free, Free):
              basis==False:
              Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(0,0,1)>]
              basis==True:
              [Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(0,0,1)>], Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(1,2)>], Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(4,1)>], Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(7)>]]
            Hom(Free, Hko):
              basis==False:
              Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(0,0,1)>]
              basis==True:
              [Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(0,0,1)>]]
            Hom(Free, Trivial):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Free, Trivial with redundant generator):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Hko, Free):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Hko, Hko):
              basis==False:
              Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(0,0,1)>]
              basis==True:
              [Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(0,0,1)>]]
            Hom(Hko, Trivial):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Hko, Trivial with redundant generator):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial, Free):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial, Hko):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial, Trivial):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial, Trivial with redundant generator):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial with redundant generator, Free):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial with redundant generator, Hko):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial with redundant generator, Trivial):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
            Hom(Trivial with redundant generator, Trivial with redundant generator):
              basis==False:
              The trivial homomorphism.
              basis==True:
              []
        """
        from .morphism import _CreateRelationsMatrix

        M = self.domain()
        N = self.codomain()

        def _trivial_case():
            '''
            The return value if there are no non-trivial homomorphisms.
            '''
            if basis:
                # Since the vector space of homomorphisms is trivial, the basis
                # is the empty set.
                return []
            else:
                # Since the vector space of homomorphisms is trivial, it contains
                # only the trivial homomorphism.
                return self.zero()

        # Deal with the trivial cases first.  Note that this covers the case
        # where the domain or codomain have no generators.
        if N.is_trivial() or M.is_trivial():
            return _trivial_case()

        # Then deal with the case where the domain has no relations.
        elif not M.has_relations():
            res = []
            num_generators = len(M.generators())
            for i, g in enumerate(M.generators()):
                # The i'th generator can go to any of these basis elements:
                base = N[(g.degree() + n)]
                for value in base:
                    values = [N.zero() if i != j else value for j in range(num_generators)]
                    res.append(Hom(M,N)(values))
                    if not basis:
                        return res[0]

        else:
            # Note that this list is non-empty since we dealt with the trivial
            # case above.
            source_degs = [g.degree() + n for g in M.generators()]

            # Note that this list is non-empty since we dealt with the free
            # case above.
            target_degs = [r.degree() + n for r in M.relations()]

            block_matrix, R = _CreateRelationsMatrix(
                N, [r.dense_coefficient_list() for r in M.relations()], source_degs, target_degs)

            ker = R.right_kernel()

            res = []
            for b in ker.basis():
                n = 0

                xs = []
                for j,X in enumerate(block_matrix[0]):
                    k = X.domain().dimension()
                    xs.append(N.element_from_coordinates(b[n:n+k], source_degs[j]))
                    n += k

                res.append(Hom(M, N)(xs))
                if not basis:
                    return res[0]

        # If the code above found a non-trivial homomorphism and ``basis==False``,
        # it will have terminated by now.
        if len(res) == 0:
            return _trivial_case()
        else:
            return res

