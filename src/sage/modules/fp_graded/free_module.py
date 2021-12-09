r"""
Finitely generated free graded left modules over connected graded algebras.

This class implements methods for construction and basic manipulation
of finitely generated free graded modules over connected graded
algebras with a given graded basis.

==========
User guide
==========

Let `p` be a prime number. The mod `p` Steenrod algebra `A_p`
is a connected algebra over the finite field of `p` elements.  All modules
presented here will be defined over `A_p`, or one of its sub-Hopf algebras.
E.g.::

    sage: A = SteenrodAlgebra(p=2)

The constructor of the module class takes as arguments the algebra
over which the module is defined and an ordered tuple of degrees for
the generators::

    sage: from sage.modules.fp_graded.free_module import FreeGradedModule
    sage: M = FreeGradedModule(algebra=A, generator_degrees=(0,1)); M
    Finitely presented free left module on 2 generators over mod 2 Steenrod algebra, milnor basis

The resulting free module will have generators in the degrees given to its
constructor::

    sage: M.generator_degrees()
    (0, 1)

The connectivity of a module over a connected graded algebra is the minimum
degree of all its module generators.  Thus, if the module is non-trivial, the
connectivity is an integer::

    sage: M.connectivity()
    0

---------------
Module elements
---------------

For an `A`-module with generators `\{g_i\}_{i=1}^N`, any homogeneous element
of degree `n` has the form

.. MATH::

    x = \sum_{i=1}^N a_i\cdot g_i\,,

where `a_i\in A_{n-\deg(g_i)}` for all `i`.  The ordered set `\{a_i\}` 
is referred to as the coefficients of `x`.

Module elements are displayed by their algebra coefficients::

    sage: M.an_element(n=5)
    <Sq(2,1), Sq(4)>

    sage: M.an_element(n=15)
    <Sq(0,0,0,1), Sq(1,2,1)>

The generators are themselves elements of the module::

    sage: M.generators()
    [<1, 0>, <0, 1>]

Producing elements from a given set of coefficients is possible using the
module class ()-method::

    sage: coeffs=[Sq(5), Sq(1,1)]
    sage: x = M(coeffs); x
    <Sq(5), Sq(1,1)>

The module action produces new elements::

    sage: Sq(2)*x
    <Sq(4,1) + Sq(7), Sq(3,1)>

Each non-zero element has a well-defined degree::

    sage: x.degree()
    5

But the zero element has not::

    sage: zero = M.zero(); zero
    0
    sage: zero.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero element does not have a well-defined degree

Any two elements can be added as long as they are in the same degree::

    sage: y = M.an_element(5); y
    <Sq(2,1), Sq(4)>
    sage: x + y
    <Sq(2,1) + Sq(5), Sq(1,1) + Sq(4)>

or when at least one of them is zero::

    sage: x + zero == x
    True

Finally, additive inverses exist::

    sage: x - x
    0

For every integer `n`, the set of module elements of degree `n` form a
vector space over the ground field `k`.  A basis for this vector space can be
computed::

    sage: M.basis_elements(5)
    [<Sq(2,1), 0>, <Sq(5), 0>, <0, Sq(1,1)>, <0, Sq(4)>]

together with a corresponding vector space presentation::

    sage: M.vector_presentation(5)
    Vector space of dimension 4 over Finite Field of size 2

Given any element, its coordinates with resepct to this basis can be computed::

    sage: v = x.vector_presentation(); v
    (0, 1, 1, 0)

Going the other way, any element can be constructed by specifying its
coordinates::

    sage: x_ = M.element_from_coordinates((0,1,1,0), 5)
    sage: x_
    <Sq(5), Sq(1,1)>
    sage: x_ == x
    True

--------------------
Module homomorphisms
--------------------

Homomorphisms of free graded `A`-modules `M\to N` are linear maps of their
underlying `k`-vector spaces which commute with the `A`-module structure.

To create a homomorphism, first create the object modelling the set of all
such homomorphisms using the free function ``Hom``::

    sage: M = FreeGradedModule(A, (0,1))
    sage: N = FreeGradedModule(A, (2,))
    sage: homspace = Hom(M, N); homspace
    Set of Morphisms from Finitely presented free left module on 2 generators over mod 2 Steenrod algebra, milnor basis to Finitely presented free left module on 1 generator over mod 2 Steenrod algebra, milnor basis in Category of finite dimensional graded modules with basis over mod 2 Steenrod algebra, milnor basis

Just as module elements, homomorphisms are created using the ()-method
of the homspace object. The only argument is a list of module elements in the
codomain, corresponding to the module generators of the domain::

    sage: g = N([1])  # the generator of the codomain module.
    sage: values = [Sq(2)*g, Sq(2)*Sq(1)*g]
    sage: f = homspace(values)

The resulting homomorphism is the one sending the `i`-th generator of the
domain to the `i`-th codomain value given::

    sage: f
    Module homomorphism of degree 4 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<Sq(2)>, <Sq(0,1) + Sq(3)>]

Convenience methods exist for creating the trivial morphism::

    sage: homspace.zero()
    The trivial homomorphism.

as well as the identity endomorphism::

    sage: Hom(M, M).identity()
    The identity homomorphism.

Homomorphisms can be evaluated on elements of the domain module::

    sage: v1 = f(Sq(7)*M.generator(0)); v1
    <Sq(3,2)>

    sage: v2 = f(Sq(17)*M.generator(1)); v2
    <Sq(11,3) + Sq(13,0,1) + Sq(17,1)>

and they respect the module action::

    sage: v1 == Sq(7)*f(M.generator(0))
    True

    sage: v2 == Sq(17)*f(M.generator(1))
    True

Any non-trivial homomorphism has a well-defined degree::

    sage: f.degree()
    4

but just as module elements, the trivial homomorphism does not::

    sage: zero_map = homspace.zero()
    sage: zero_map.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero morphism does not have a well-defined degree

Any two homomorphisms can be added as long as they are of the same degree::

    sage: f2 = homspace([Sq(2)*g, Sq(3)*g])
    sage: f + f2
    Module homomorphism of degree 4 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [0, <Sq(0,1)>]

or when at least one of them is zero::

    sage: f + zero_map == f
    True

Finally, additive inverses exist::

    sage: f - f
    The trivial homomorphism.

The restriction of a homomorphism to the vector space of `n`-dimensional module
elements is a linear transformation::

    sage: f_4 = f.vector_presentation(4); f_4
    Vector space morphism represented by the matrix:
    [0 1 0]
    [1 1 1]
    [0 1 0]
    [0 0 0]
    Domain: Vector space of dimension 4 over Finite Field of size 2
    Codomain: Vector space of dimension 3 over Finite Field of size 2

This is compatible with the vector presentations of its domain and codomain
modules::

    sage: f.domain() is M
    True
    sage: f.codomain() is N
    True
    sage: f_4.domain() is M.vector_presentation(4)
    True
    sage: f_4.codomain() is N.vector_presentation(4 + f.degree())
    True


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

from sage.misc.cachefunc import cached_method
from sage.modules.free_module import VectorSpace
from sage.rings.infinity import PlusInfinity
from sage.categories.graded_modules import GradedModules
from sage.combinat.free_module import CombinatorialFreeModule

from .free_element import FreeGradedModuleElement

class FreeGradedModule(CombinatorialFreeModule):
    r"""
    Create a finitely generated free graded module over a connected
    graded algebra, with generators in specified degrees.

    INPUT:

    - ``algebra`` -- the graded connected algebra over which the module is
      defined. This algebra must be equipped with a graded basis.

    - ``generator_degrees`` -- tuple of integers defining the number
      of generators of the module, and their degrees.

    TESTS::

        sage: from sage.modules.fp_graded.free_module import FreeGradedModule
        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: M = FreeGradedModule(E, (-1,3))
        sage: M
        Finitely presented free left module on 2 generators over The exterior algebra of rank 3 over Rational Field
        sage: M.generator_degrees()
        (-1, 3)
        sage: a, b = M.generators()
        sage: (x*y*b).degree()
        5

        sage: A = SteenrodAlgebra(2)
        sage: FreeGradedModule(A, (-2,2,4))
        Finitely presented free left module on 3 generators over mod 2 Steenrod algebra, milnor basis
    """
    def __init__(self, algebra, generator_degrees):
        r"""
        Create a finitely generated free graded module over a connected graded
        algebra.

        OUTPUT: The finitely generated free graded module on generators with
        degrees given by ``generator_degrees``.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: TestSuite(FreeGradedModule(A, (-2,2,4))).run()
        """
        # If generator_degrees is [d_0, d_1, ...], then
        # the generators are indexed by (0,d_0), (1,d_1), ...
        keys = list(enumerate(generator_degrees))
        self._generator_keys = keys

        if not algebra.base_ring().is_field():
            raise ValueError('the ground ring of the algebra must be a field')

        # Call the base class constructor.
        CombinatorialFreeModule.__init__(self, algebra,
                                         basis_keys=keys,
                                         element_class=FreeGradedModuleElement,
                                         category=GradedModules(algebra).WithBasis().FiniteDimensional())


    def generator_degrees(self):
        r"""
        The degrees of the module generators.

        OUTPUT: A tuple containing the degrees of the generators for this
        module, in the order that the generators were given when this module
        was constructed.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (-2,2,4))
            sage: M.generator_degrees()
            (-2, 2, 4)
        """
        return tuple(a[1] for a in self._generator_keys)


    def is_trivial(self):
        r"""
        Decide if this module is trivial or not.

        OUTPUT: The boolean value ``True`` if the module is trivial, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: FreeGradedModule(A, (-2,2,4)).is_trivial()
            False
            sage: FreeGradedModule(A, ()).is_trivial()
            True
        """
        return not len(self.generator_degrees())


    def connectivity(self):
        r"""
        The connectivity of this module.

        OUTPUT: An integer equal to the minimal degree of all the generators, if
        this module is non-trivial.  Otherwise, `+\infty`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (-2,2,4))
            sage: M.connectivity()
            -2

        TESTS::

            sage: M = FreeGradedModule(A, ())
            sage: M.is_trivial()
            True
            sage: M.connectivity()
            +Infinity
        """
        return min(self.generator_degrees() + (PlusInfinity(),))


    def _element_constructor_(self, coefficients):
        r"""
        Construct any element of the module.

        This function is used internally by the ()-method when creating
        module elements, and should not be called by the user explicitly.

        INPUT:

        - ``coefficients`` -- A tuple of coefficient (i.e. elements of the
        algebra for this module), an element of FreeGradedModule, or the zero integer
        constant.

        OUTPUT: An instance of the element class with coefficients from
        ``coefficients``, the element ``coefficients`` if it already was an
        element, or the zero module element.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))

            sage: zero = M(0); zero
            0

            sage: e = M((Sq(4), Sq(2), 1)); e
            <Sq(4), Sq(2), 1>

            sage: e is M(e)
            True
        """
        if isinstance(coefficients, self.element_class):
            return coefficients
        elif not coefficients:
            return self.zero()
        else:
            B = self.basis()
            return sum(c*B[b] for (c,b) in zip(coefficients, self._generator_keys))


    def an_element(self, n=None):
        r"""
        Return an element of the module.

        This function chooses deterministically an element of the module in the
        given degree.

        INPUT:

        - ``n`` --  the degree of the element to construct.  If the default
          value ``None`` is given, a degree will be chosen by the function.

        OUTPUT: An element of the given degree.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M.an_element(172)
            <Sq(0,0,2,0,1,0,1), Sq(0,4,0,0,1,0,1), Sq(7,1,0,0,1,0,1)>

        Zero is the only element in the trivial module::

            sage: FreeGradedModule(A, ()).an_element()
            0
        """
        if not len(self.generator_degrees()):
            return self.zero()

        if n is None:
            n = max(self.generator_degrees()) + 7

        coefficients = []
        for g in self.generator_degrees():
            basis = self.base_ring().basis(n - g) if n >= g else ()
            # All of the algebra generators in basis will bring the
            # module generator in dimension g to dimension
            # g + (topDimension - g) = topDimension.  Picking any one of them
            # will do, so we pick the one with index (g (mod l)).
            l = len(basis)
            if l:
                coefficients.append(basis[g % l])
            else:
                coefficients.append(self.base_ring().zero())

        return self(coefficients)


    def _repr_(self):
        r"""
        Construct a string representation of the module.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M._repr_()
            'Finitely presented free left module on 3 generators over mod 2 Steenrod algebra, milnor basis'
        """
        return "Finitely presented free left module on %s generator%s over %s"\
            %(len(self.generator_degrees()), "" if len(self.generator_degrees()) == 1 else "s",
              self.base_ring())


    @cached_method
    def basis_elements(self, n):
        r"""
        A basis for the vector space of degree ``n`` module elements.
        This returns a basis as a vector space over the base field,
        not a basis as a free module over the algebra.

        INPUT:

        - ``n`` -- an integer.

        OUTPUT: A sequence of homogeneous module elements of degree ``n``
        which is a basis for the vector space of all degree ``n`` module
        elements.

        .. SEEALSO::
            :meth:`vector_presentation`, :meth:`element_from_coordinates`

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M.basis_elements(8)
            [<Sq(1,0,1), 0, 0>,
             <Sq(2,2), 0, 0>,
             <Sq(5,1), 0, 0>,
             <Sq(8), 0, 0>,
             <0, Sq(0,2), 0>,
             <0, Sq(3,1), 0>,
             <0, Sq(6), 0>,
             <0, 0, Sq(1,1)>,
             <0, 0, Sq(4)>]
        """
        basis_n = []
        for i, generator_degree in enumerate(self.generator_degrees()):
            l = n - generator_degree
            basis_n += [a*self.generator(i) for a in self.base_ring().basis(l)]

        return basis_n


    @cached_method
    def element_from_coordinates(self, coordinates, n):
        r"""
        The module element of degree ``n`` having the given coordinates
        with respect to the basis of module elements given by
        :meth:`basis_elements`.

        INPUT:

        - ``coordinates`` -- a sequence of elements of the ground field.
        - ``n`` -- an integer.

        OUTPUT: A module element of degree ``n``.

        .. SEEALSO::
            :meth:`vector_presentation`, and :meth:`basis_elements`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: x = M.element_from_coordinates((0,1,0,1), 5); x
            <Sq(5), Sq(4)>
            sage: basis = M.basis_elements(5)
            sage: y = 0*basis[0] + 1*basis[1] + 0*basis[2] + 1*basis[3]
            sage: x == y
            True

            sage: M.element_from_coordinates((0,0,0,0), 5)
            0
        """
        basis_elements = self.basis_elements(n)

        if len(coordinates) != len(basis_elements):
            raise ValueError('the given coordinate vector has incorrect length: %d.  '
                  'It should have length %d' % (len(coordinates), len(basis_elements)))

        # Adding the condition `if c != 0` improved performance dramatically in this
        # real life example:
        #
        # sage: rels = [ [Sq(1),0,0,0], [Sq(2),0,0,0], [Sq(4),0,0,0], [Sq(8),0,0,0], [0,Sq(1),0,
        # ....: 0], [0,Sq(2),0,0], [0,Sq(4),0,0], [Sq(31),Sq(14),0,0], [0,Sq(20),0,0], [0,0,Sq(1
        # ....: ),0], [0,0,Sq(2),0], [0,Sq(31),Sq(6),0], [0,0,Sq(8),0], [0,0,0,Sq(1)], [0,0,Sq(3
        # ....: 1),Sq(2)], [0,0,0,Sq(4)], [0,0,0,Sq(8)] ]
        # ....:
        # ....: M = FPA_Module([0, 17, 42, 71], A, relations=rels)
        # sage: res = M.resolution(2, top_dim=30, verbose=True)
        #
        # This function was called a total of 2897 times during the computation,
        # and the total running time of the entire computation dropped from
        # 57 to 21 seconds by adding the optimization.
        #
        element = sum([c*element for c, element in zip(coordinates, basis_elements) if c != 0])
        if element == 0:
            # The previous sum was over the empty list, yielding the integer
            # 0 as a result, rather than a module element.
            # Fix this by returning the zero element.
            return self.zero()
        else:
            # The sum defining element is of the correct type, so return it.
            return element


    def __getitem__(self, n):
        r"""
        A vector space isomorphic to the vector space of module elements of
        degree ``n``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: V = M[4]; V
            Vector space of dimension 4 over Finite Field of size 2
            sage: V.dimension()
            4

        .. SEEALSO::
            This function is an alias for :meth:`vector_presentation`.
        """
        return self.vector_presentation(n)


    @cached_method
    def vector_presentation(self, n):
        r"""
        A vector space over the ground field of the module algebra,
        isomorphic to the degree ``n`` elements of this module.

        Let `\mathcal{k}` be the ground field of the algebra over this module is defined,
        and let `M_n` be the vector space of module elements of degree ``n``.

        The return value of this function is the vector space
        `\mathcal{k}^{r}` where `r = dim(M_n)`.

        The isomorphism between `k^{r}` and `M_n` is given by the
        bijection taking the standard basis element `e_i` to the `i`-th
        element of the array returned by :meth:`basis_elements`.

        INPUT:

        - ``n`` -- an integer degree.

        OUTPUT: A vector space over the ground field of the algebra over which
        this module is defined, isomorphic to the vector space of module
        elements of degree ``n``.

        .. SEEALSO::
            :meth:`basis_elements`, :meth:`element_from_coordinates`

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A1 = SteenrodAlgebra(2, profile=[2,1])
            sage: M = FreeGradedModule(A1, (0,))
            sage: M.vector_presentation(3)
            Vector space of dimension 2 over Finite Field of size 2
            sage: M.basis_elements(3)
            [<Sq(0,1)>, <Sq(3)>]
            sage: [M.vector_presentation(i).dimension() for i in range(-2, 9)]
            [0, 0, 1, 1, 1, 2, 1, 1, 1, 0, 0]
        """
        return VectorSpace(self.base_ring().base_ring(), len(self.basis_elements(n)))


    @cached_method
    def generator(self, index):
        r"""
        Return the module generator with the given index.

        OUTPUT: An instance of the element class of this parent.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M.generator(0)
            <1, 0, 0>
            sage: M.generator(1)
            <0, 1, 0>
            sage: M.generator(2)
            <0, 0, 1>
        """
        try:
            return self.monomial(self._generator_keys[index])
        except IndexError:
            raise ValueError('the parent module has generators in the index '\
                'range [0, %s]; generator %s does not exist' %\
                (len(self.generator_degrees()) - 1, index))


    def generators(self):
        r"""
        Return all the module generators.

        OUTPUT: A list consisting instances of the element class of this
        parent.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: M.generators()
            [<1, 0>, <0, 1>]
        """
        return [self.generator(i) for i in range(len(self.generator_degrees()))]


    def _Hom_(self, Y, category):
        r"""
        The internal hook used by the free function
        :meth:`sage.categories.homset.hom.Hom` to create homsets involving this
        parent.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: M._Hom_(M, category=None)
            Set of Morphisms from Finitely presented free left module on 2 generators over mod 2 Steenrod algebra, milnor basis to Finitely presented free left module on 2 generators over mod 2 Steenrod algebra, milnor basis in Category of finite dimensional graded modules with basis over mod 2 Steenrod algebra, milnor basis
        """
        from .free_homspace import FreeGradedModuleHomspace
        return FreeGradedModuleHomspace(self, Y, category)


    def suspension(self, t):
        r"""
        Suspend the module by the given integer degree.

        INPUT:

        - ``t`` -- An integer.

        OUTPUT: A module which is isomorphic to this module by a shift
        of degrees by the integer ``t``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M.suspension(4).generator_degrees()
            (4, 6, 8)
            sage: M.suspension(-4).generator_degrees()
            (-4, -2, 0)
        """
        return FreeGradedModule(algebra=self.base_ring(),
            generator_degrees=tuple([g + t for g in self.generator_degrees()]))


    def to_fp_module(self):
        """
        Create a finitely presented module from this free module.

        OUTPUT: the finitely presented module having same set of generators
        as this module, no relations.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: F = FreeGradedModule(A, (-2,2,4))
            sage: F.to_fp_module()
            Finitely presented left module on 3 generators and 0 relations over mod 2 Steenrod algebra, milnor basis
        """
        from .module import FPModule
        return FPModule(algebra=self.base_ring(),
                   generator_degrees=self.generator_degrees(),
                   relations=())
        
