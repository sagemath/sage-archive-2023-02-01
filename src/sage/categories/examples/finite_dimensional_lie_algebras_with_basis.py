r"""
Examples of a finite dimensional Lie algebra with basis
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import LieAlgebras
from sage.modules.free_module import FreeModule
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.examples.lie_algebras import LieAlgebraFromAssociative as BaseExample

class AbelianLieAlgebra(Parent, UniqueRepresentation):
    r"""
    An example of a finite dimensional Lie algebra with basis:
    the abelian Lie algebra.

    Let `R` be a commutative ring, and `M` an `R`-module. The
    *abelian Lie algebra* on `M` is the `R`-Lie algebra
    obtained by endowing `M` with the trivial Lie bracket
    (`[a, b] = 0` for all `a, b \in M`).

    This class illustrates a minimal implementation of a finite dimensional
    Lie algebra with basis.

    INPUT:

    - ``R`` -- base ring

    - ``names`` -- list of strings, to be used as identifiers for
      the basis elements of `M`

    - ``M`` -- an `R`-module (default: the free `R`-module of
      rank ``len(names)``) to serve as the ground space for the
      Lie algebra

    - ``ambient`` -- (optional) a Lie algebra; if this is set,
      then the resulting Lie algebra is declared a Lie subalgebra
      of ``ambient``

    OUTPUT:

    The abelian Lie algebra on `M`.

    .. TODO::

       Have I correctly explained these?

    .. TODO::

       Why am I required to specify ``names`` if the actual
       elements of the Lie algebra end up being anonymous
       vectors anyway?

           sage: from sage.categories.examples.finite_dimensional_lie_algebras_with_basis import *
           sage: AbelianLieAlgebra(QQ, ['x','y'])
           An example of a finite dimensional Lie algebra with
            basis: the abelian Lie algebra with generators
            ('x', 'y') over Rational Field
           sage: _.gens()
           ((1, 0), (0, 1))

    """
    @staticmethod
    def __classcall_private__(cls, R, names, M=None, ambient=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.categories.examples.finite_dimensional_lie_algebras_with_basis import AbelianLieAlgebra
            sage: A1 = AbelianLieAlgebra(QQ, 'x,y,z')
            sage: A2 = AbelianLieAlgebra(QQ, ['x','y','z'])
            sage: A3 = AbelianLieAlgebra(QQ, ['x','y','z'], FreeModule(QQ, 3))
            sage: A1 is A2 and A2 is A3
            True
        """
        if isinstance(names, str):
            names = names.split(',')
        if M is None:
            M = FreeModule(R, len(names))
        elif len(names) != M.dimension():
            raise ValueError("number of generators is not correct")
        else:
            M = M.change_ring(R)
        return super(AbelianLieAlgebra, cls).__classcall__(cls, R, tuple(names), M, ambient)

    def __init__(self, R, names, M, ambient):
        """
        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: TestSuite(L).run()
        """
        self._ordered_indices = names
        self._M = M
        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        if ambient is None:
            ambient = self
        else:
            cat = cat.Subobjects()
        self._ambient = ambient
        Parent.__init__(self, base=R, names=names, category=cat)

    def _repr_(self):
        """
        EXAMPLES::

            sage: LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            An example of a finite dimensional Lie algebra with basis:
             the abelian Lie algebra with generators ('a', 'b', 'c')
             over Rational Field
        """
        ret = "An example of a finite dimensional Lie algebra with basis:" \
              " the abelian Lie algebra with generators {!r} over {}".format(
                         self.variable_names(), self.base_ring())
        B = self._M.basis_matrix()
        if not B.is_one():
            ret += " with basis matrix:\n{!r}".format(B)
        return ret

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L(0)
            (0, 0, 0)
            sage: M = FreeModule(ZZ, 3)
            sage: L(M([1, -2, 2]))
            (1, -2, 2)
            sage: a,b,c = L.lie_algebra_generators()
            sage: X = L.subalgebra([a+b, 2*a+c])
            sage: x,y = X.basis()
            sage: L(x)
            (1, 0, 1)
            sage: L(x+y)
            (1, 1, -1)
        """
        if isinstance(x, AbelianLieAlgebra.Element):
            x = x.value
        return self.element_class(self, self._M(x))

    @cached_method
    def zero(self):
        """
        Return the zero element.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.zero()
            (0, 0, 0)
        """
        return self.element_class(self, self._M.zero())

    def basis_matrix(self):
        """
        Return the basis matrix of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.basis_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return self._M.basis_matrix()

    def ambient(self):
        """
        Return the ambient Lie algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.inject_variables()
            Defining a, b, c
            sage: S = L.subalgebra([2*a+b, b + c], 'x,y')
            sage: S.ambient() == L
            True
        """
        return self._ambient

    def subalgebra(self, gens, names='x'):
        """
        Return the Lie subalgebra of ``self`` generated by the
        elements of the iterable ``gens``.

        This currently requires the ground ring `R` to be a field.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.inject_variables()
            Defining a, b, c
            sage: L.subalgebra([2*a+b, b + c], 'x,y')
            An example of a finite dimensional Lie algebra with basis:
             the abelian Lie algebra with generators ('x', 'y')
             over Rational Field with basis matrix:
            [   1    0 -1/2]
            [   0    1    1]
        """
        if isinstance(names, str):
            names = names.split(',')
        if len(names) == 1 and len(gens) != 1:
            names = tuple( names[0] + str(i) for i in range(len(gens)) )
        N = self._M.subspace([g.value for g in gens])
        return AbelianLieAlgebra(self.base_ring(), names, N, self._ambient)

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.basis()
            Finite family {'a': (1, 0, 0), 'c': (0, 0, 1), 'b': (0, 1, 0)}
        """
        names = self.variable_names()
        d = {names[i]: self.element_class(self, b)
             for i,b in enumerate(self._M.basis())}
        return Family(d)

    lie_algebra_generators = basis

    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.gens()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        """
        G = self.lie_algebra_generators()
        return tuple(G[i] for i in self.variable_names())

    def free_module(self):
        """
        Return ``self`` as a free module.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.free_module()
            Vector space of dimension 3 over Rational Field

            sage: L.inject_variables()
            Defining a, b, c
            sage: S = L.subalgebra([2*a+b, b + c], 'x,y')
            sage: S.free_module()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 -1/2]
            [   0    1    1]
        """
        return self._M

    class Element(BaseExample.Element):
        def __iter__(self):
            """
            Iterate over ``self`` by returning pairs ``(i, c)`` where ``i``
            is the index of the basis element and ``c`` is the corresponding
            coefficient.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.inject_variables()
                Defining a, b, c
                sage: elt = 2*a - c
                sage: list(elt)
                [('a', 2), ('c', -1)]
            """
            I = self.parent()._ordered_indices
            zero = self.parent().base_ring().zero()
            for i, c in self.value.iteritems():
                if c != zero:
                    yield (I[i], c)

        def __getitem__(self, i):
            """
            Redirect the ``__getitem__()`` to the wrapped element unless
            ``i`` is a basis index.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.inject_variables()
                Defining a, b, c
                sage: elt = 2*a + b - c
                sage: elt[0]
                2
                sage: elt['a']
                2
                sage: elt['c']
                -1
            """
            if i in self.parent()._ordered_indices:
                i = self.parent()._ordered_indices.index(i)
            return self.value.__getitem__(i)

        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.inject_variables()
                Defining a, b, c
                sage: a.bracket(c)
                (0, 0, 0)
                sage: a.bracket(b).bracket(c)
                (0, 0, 0)
            """
            return self.parent().zero()

        def lift(self):
            """
            Return the lift of ``self`` to the universal enveloping algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.inject_variables()
                Defining a, b, c
                sage: elt = 2*a + 2*b + 3*c
                sage: elt.lift()
                2*a + 2*b + 3*c
            """
            UEA = self.parent().universal_enveloping_algebra()
            gens = UEA.gens()
            return UEA.sum(c * gens[i] for i, c in self.value.iteritems())

        def to_vector(self):
            """
            Return ``self`` as a vector.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.inject_variables()
                Defining a, b, c
                sage: elt = 2*a + 2*b + 3*c
                sage: elt.to_vector()
                (2, 2, 3)
            """
            return self.value

Example = AbelianLieAlgebra

