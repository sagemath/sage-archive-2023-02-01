r"""
Homomorphisms of finitely generated free graded left modules

For an overview of the free module API, see :doc:`free_module`.

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

#*****************************************************************************
#       Copyright (C) 2019 Robert R. Bruner <rrb@math.wayne.edu>
#                     and  Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from inspect import isfunction

from sage.categories.homset import Hom
from sage.categories.morphism import Morphism


class FreeGradedModuleMorphism(Morphism):
    r"""
    Create a homomorphism between finitely generated free graded modules.

    INPUT:

    - ``parent`` -- a homspace in the category of finitely generated free
      modules

    - ``values`` -- a list of elements in the codomain; each element
      corresponds (by their ordering) to a module generator in the domain

    EXAMPLES::

        sage: from sage.modules.fp_graded.free_module import FreeGradedModule
        sage: A = SteenrodAlgebra(2)
        sage: F1 = FreeGradedModule(A, (4,5), names='b')
        sage: F2 = FreeGradedModule(A, (3,4), names='c')
        sage: F3 = FreeGradedModule(A, (2,3), names='d')
        sage: H1 = Hom(F1, F2)
        sage: H2 = Hom(F2, F3)
        sage: f = H1( ( F2((Sq(4), 0)), F2((0, Sq(4))) ) )
        sage: g = H2( ( F3((Sq(2), 0)), F3((Sq(3), Sq(2))) ) )
        sage: g*f
        Module homomorphism of degree 4 defined by sending the generators
          [b_{4}, b_{5}]
        to
          [(Sq(0,2)+Sq(3,1)+Sq(6))*d_{2}, (Sq(1,2)+Sq(7))*d_{2} + (Sq(0,2)+Sq(3,1)+Sq(6))*d_{3}]

    TESTS::

    A non-example because the degree is not well-defined::

        sage: M = FreeGradedModule(A, (0, 0))
        sage: N = FreeGradedModule(A, (0,))
        sage: H = Hom(M, N)
        sage: g = N.generator(0)
        sage: H([Sq(1)*g, Sq(2)*g])
        Traceback (most recent call last):
        ...
        ValueError: ill-defined homomorphism: degrees do not match
    """

    def __init__(self, parent, values):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0, 0))
            sage: N = FreeGradedModule(A, (0,))
            sage: H = Hom(M, N)
            sage: g = N.generator(0)
            sage: TestSuite(H).run()
            sage: TestSuite(g).run()
        """
        from .free_homspace import FreeGradedModuleHomspace
        if not isinstance(parent, FreeGradedModuleHomspace):
            raise TypeError('the parent (%s) must be a f.p. free module homset' % parent)

        # Get the values.
        C = parent.codomain()
        D = parent.domain()
        if isfunction(values):
            _values = [C(values(g)) for g in D.generators()]
        elif values == 0:
            _values = len(D.generator_degrees())*[C(0)]
        else:
            _values = [C(a) for a in values]

        # Check the homomorphism is well defined.
        if len(D.generator_degrees()) != len(_values):
            raise ValueError('the number of values must equal the number of '
                'generators in the domain; invalid argument: %s' % _values)

        # Compute the degree.
        if all(v.is_zero() for v in _values):
            # The zero homomorphism does not get a degree.
            _degree = None
        else:
            degrees = []
            for i, value in enumerate(_values):
                if not value.is_zero():
                    x = value.degree()
                    xx = D.generator_degrees()[i]
                    degrees.append(x-xx)

            _degree = min(degrees)
            if _degree != max(degrees):
                raise ValueError('ill-defined homomorphism: degrees do not match')

        self._degree = _degree
        self._values = _values

        Morphism.__init__(self, parent)


    def degree(self):
        r"""
        The degree of ``self``.

        OUTPUT:

        The degree of this homomorphism. Raise an error if this is
        the trivial homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeGradedModule(A, (0,1)), FreeGradedModule(A, (0,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f.degree()
            5

        The zero homomorphism has no degree::

            sage: homspace.zero().degree()
            Traceback (most recent call last):
            ...
            ValueError: the zero morphism does not have a well-defined degree
        """
        if self._degree is None:
            # The zero morphism has no degree.
            raise ValueError("the zero morphism does not have a well-defined degree")
        return self._degree


    def values(self):
        r"""
        The values under ``self`` corresponding to the generators of
        the domain module.

        OUTPUT:

        A sequence of elements of the codomain module.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeGradedModule(A, (0,1)), FreeGradedModule(A, (2,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f.values()
            [Sq(5)*g_{2}, Sq(3,1)*g_{2}]
            sage: homspace.zero().values()
            [0, 0]
        """
        return self._values


    def _richcmp_(self, other, op):
        r"""
        Compare this homomorphism to the given homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeGradedModule(A, (0,1)), FreeGradedModule(A, (2,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f._richcmp_(f, op=2)
            True
            sage: f._richcmp_(f, op=3)
            False
        """
        try:
            same = (self - other).is_zero()
        except ValueError:
            return False

        # Equality
        if op == 2:
            return same

        # Non-equality
        if op == 3:
            return not same

        return False


    def _add_(self, g):
        r"""
        The pointwise sum of `self`` and ``g``.

        Pointwise addition of two homomorphisms `f` and `g` with the same
        domain and codomain is given by the formula `(f+g)(x) = f(x) + g(x)`
        for every `x` in the domain of `f`.

        INPUT:

        - ``g`` -- a homomorphism with the same domain and codomain as ``self``

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeGradedModule(A, (0,1)), FreeGradedModule(A, (2,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: ff = f.__add__(f)
            sage: ff.is_zero()
            True
            sage: ff.__add__(f) == f
            True
            sage: ff = f + f
            sage: ff.is_zero()
            True
        """
        if self.domain() != g.domain():
            raise ValueError('morphisms do not have the same domain')
        elif self.codomain() != g.codomain():
            raise ValueError('morphisms do not have the same codomain')
        elif self.is_zero():
            return g
        elif g.is_zero():
            return self
        elif self.degree() and g.degree() and self.degree() != g.degree():
            raise ValueError('morphisms do not have the same degree')

        v = [self(x) + g(x) for x in self.domain().generators()]
        return self.parent()(v)


    def _neg_(self):
        r"""
        The additive inverse of ``self`` with respect to the group
        structure given by pointwise sum.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeGradedModule(A, (0,1)), FreeGradedModule(A, (2,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f_inverse = -f; f_inverse
            Module homomorphism of degree 7 defined by sending the generators
              [g_{0}, g_{1}]
            to
              [Sq(5)*g_{2}, Sq(3,1)*g_{2}]
            sage: (f + f_inverse).is_zero()
            True
        """
        return self.parent()([-x for x in self.values()])


    def _sub_(self, g):
        r"""
        The pointwise difference between ``self`` and ``g``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeGradedModule(A, (0,1)), FreeGradedModule(A, (2,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: values2 = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: g = homspace(values2)
            sage: f - g
            The trivial homomorphism
        """
        return self + (-g)


    # Define __mul__ rather than _mul_, since we want to allow
    # "multiplication" by morphisms from different homsets.
    def __mul__(self, g):
        r"""
        The composition of ``g`` followed by ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: N = FreeGradedModule(A, (2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: values2 = [Sq(2)*M.generator(0)]
            sage: g = Hom(N, M)(values2)
            sage: fg = f * g; fg
            Module homomorphism of degree 7 defined by sending the generators
              [g_{2}]
            to
              [(Sq(4,1)+Sq(7))*g_{2}]
            sage: fg.is_endomorphism()
            True

        TESTS::

            sage: fg == f.__mul__(g)
            True
            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f * f
            Traceback (most recent call last):
            ...
            ValueError: morphisms are not composable
        """
        if self.parent().domain() != g.parent().codomain():
            raise ValueError('morphisms are not composable')
        homset = Hom(g.parent().domain(), self.parent().codomain())
        return homset([self(g(x)) for x in g.domain().generators()])


    def is_zero(self):
        r"""
        Decide if ``self`` is trivial.

        OUTPUT:

        The boolean value ``True`` if this homomorphism is trivial, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: N = FreeGradedModule(A, (2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.is_zero()
            False
            sage: (f-f).is_zero()
            True
        """
        return all(v.is_zero() for v in self.values())

    __bool__ = is_zero


    def is_identity(self):
        r"""
        Return if ``self`` is the identity endomorphism.

        OUTPUT:

        The boolean value ``True`` if this homomorphism is the identity
        and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: N = FreeGradedModule(A, (2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.is_identity()
            False
            sage: id = Hom(M, M)(M.generators()); id
            The identity homomorphism
            sage: id.is_identity()
            True
        """
        return (self.parent().is_endomorphism_set() and
                self.parent().identity() == self)

    def __call__(self, x):
        r"""
        Evaluate the homomorphism at the given domain element ``x``.

        INPUT:

        - ``x`` -- an element of the domain of this morphism

        OUTPUT:

        The module element of the codomain which is the value of ``x``
        under this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: N = FreeGradedModule(A, (2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.__call__(M.generator(0))
            Sq(5)*g_{2}
            sage: f.__call__(M.generator(1))
            Sq(3,1)*g_{2}
        """
        if x.parent() != self.domain():
            raise ValueError('cannot evaluate morphism on element not in the domain')

        value = sum((c * v for c, v in zip(x.dense_coefficient_list(), self.values())),
                    self.codomain().zero())

        return value


    def _repr_(self):
        r"""
        A string representation of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: N = FreeGradedModule(A, (2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]

            sage: Hom(M, N)(values)._repr_()
            'Module homomorphism of degree 7 defined by sending the generators\n
              [g_{0}, g_{1}]\nto\n  [Sq(5)*g_{2}, Sq(3,1)*g_{2}]'

            sage: Hom(M, N).zero()._repr_()
            'The trivial homomorphism'

            sage: Hom(M, M).identity()._repr_()
            'The identity homomorphism'
        """
        if self.is_zero():
            return "The trivial homomorphism"
        elif self.is_identity():
            return "The identity homomorphism"
        else:
            r = "Module homomorphism of degree {} defined by sending the generators\n  {}\nto\n  {}"
            return r.format(self.degree(), self.domain().generators(), self.values())


    def vector_presentation(self, n):
        r"""
        The restriction of ``self`` to the domain module elements of degree ``n``.

        The restriction of a non-zero module homomorphism to the vector space of
        module elements of degree `n` is a linear function into the vector space
        of elements of degree `n+d` belonging to the codomain.  Here `d` is the
        degree of this homomorphism.

        When this homomorphism is zero, it has no well defined degree so the
        function cannot be presented since we do not know the degree of its
        codomain.  In this case, an error is raised.

        INPUT:

        - ``n`` -- an integer degree

        OUTPUT:

        A linear function of finite dimensional vector spaces over the
        ground field of the algebra for this module.  The domain is isomorphic
        to the vector space of domain elements of degree ``n`` of this free
        module, via the choice of basis given by
        :meth:`sage.modules.fp_graded.free_module.FreeGradedModule.basis_elements`.
        If the morphism is zero, an error is raised.

        .. SEEALSO::

            :meth:`sage.modules.fp_graded.free_module.FreeGradedModule.vector_presentation`,
            :meth:`sage.modules.fp_graded.free_module.FreeGradedModule.basis_elements`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: N = FreeGradedModule(A, (2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.vector_presentation(0)
            Vector space morphism represented by the matrix:
            [0 1]
            Domain: Vector space of dimension 1 over Finite Field of size 2
            Codomain: Vector space of dimension 2 over Finite Field of size 2
            sage: f.vector_presentation(1)
            Vector space morphism represented by the matrix:
            [0 0 0]
            [0 1 0]
            Domain: Vector space of dimension 2 over Finite Field of size 2
            Codomain: Vector space of dimension 3 over Finite Field of size 2
            sage: f.vector_presentation(2)
            Vector space morphism represented by the matrix:
            [0 0 1 1]
            [0 0 0 0]
            Domain: Vector space of dimension 2 over Finite Field of size 2
            Codomain: Vector space of dimension 4 over Finite Field of size 2

        TESTS::

            sage: F = FreeGradedModule(A, (0,))
            sage: z = Hom(F, F)([0])
            sage: z.is_zero()
            True
            sage: z.vector_presentation(0)
            Traceback (most recent call last):
            ...
            ValueError: the zero map has no vector presentation
        """
        # The trivial map has no degree, so we can not create the codomain
        # of the linear transformation.
        if self.is_zero():
            raise ValueError("the zero map has no vector presentation")

        D_n = self.domain().vector_presentation(n)
        C_n = self.codomain().vector_presentation(n + self.degree())

        values = [self(e) for e in self.domain().basis_elements(n)]
        return Hom(D_n, C_n)([
            C_n.zero() if e.is_zero() else e.vector_presentation() for e in values])

    def to_fp_module(self):
        r"""
        Create a finitely presented module from ``self``.

        OUTPUT:

        The finitely presented module with presentation equal to ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: F1 = FreeGradedModule(A, (2,))
            sage: F2 = FreeGradedModule(A, (0,))
            sage: v = F2([Sq(2)])
            sage: pres = Hom(F1, F2)([v])
            sage: M = pres.to_fp_module(); M
            Finitely presented left module on 1 generator and 1 relation over
             mod 2 Steenrod algebra, milnor basis
            sage: M.generator_degrees()
            (0,)
            sage: M.relations()
            [Sq(2)*g_{0}]
        """
        from .module import FPModule
        return FPModule(algebra=self.base_ring(),
                         generator_degrees=self.codomain().generator_degrees(),
                         relations=tuple([r.dense_coefficient_list() for r in self.values()]))

