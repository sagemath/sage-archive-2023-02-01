# -*- coding: utf-8 -*-
"""
Homology and cohomology with a basis

This module provides homology and cohomology vector spaces suitable
for computing cup products and cohomology operations.

REFERENCES:

.. [G-DR03] \R. González-Díaz and P. Réal, *Computation of cohomology
   operations on finite simplicial complexes* in Homology,
   Homotopy and Applications 5 (2003), 83-93.

.. [G-DR99] \R. González-Díaz and P. Réal, *A combinatorial method for
   computing Steenrod squares* in J. Pure Appl. Algebra 139 (1999), 89-108.

AUTHORS:

- John H. Palmieri, Travis Scrimshaw (2015-09)
"""

########################################################################
#       Copyright (C) 2015 John H. Palmieri <palmieri@math.washington.edu>
#                          Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.misc.cachefunc import cached_method
from sage.categories.algebras import Algebras
from sage.categories.modules import Modules
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.sets.family import Family
from simplicial_complex import SimplicialComplex

class HomologyVectorSpaceWithBasis(CombinatorialFreeModule):
    r"""
    Homology (or cohomology) vector space.

    This provides enough structure to allow the computation of cup
    products and cohomology operations. See the class
    :class:`CohomologyRing` (which derives from this) for examples.

    It also requires field coefficients (hence the "VectorSpace" in
    the name of the class).

    .. NOTE::

        This is not intended to be created directly by the user, but
        instead via the methods
        :meth:`~sage.homology.cell_complex.GenericCellComplex.homology_with_basis` and
        :meth:`~sage.homology.cell_complex.GenericCellComplex.cohomology_ring`
        for the class of :class:`cell
        complexes<sage.homology.cell_complex.GenericCellComplex>`.

    INPUT:

    - ``base_ring`` -- must be a field
    - ``cell_complex`` -- the cell complex whose homology we are
      computing
    - ``cohomology`` -- (default: ``False``) if ``True``, return
      the cohomology as a module
    - ``category`` -- (optional) a subcategory of modules with basis

    EXAMPLES:

    Homology classes are denoted by ``h_{d,i}`` where ``d`` is the
    degree of the homology class and ``i`` is their index in the list
    of basis elements in that degree. Cohomology classes are denoted
    ``h^{1,0}``::

        sage: RP2 = cubical_complexes.RealProjectivePlane()
        sage: RP2.homology_with_basis(GF(2))
        Homology module of Cubical complex with 21 vertices and 81 cubes
         over Finite Field of size 2
        sage: RP2.cohomology_ring(GF(2))
        Cohomology ring of Cubical complex with 21 vertices and 81 cubes
         over Finite Field of size 2
        sage: simplicial_complexes.Torus().homology_with_basis(QQ)
        Homology module of Minimal triangulation of the torus
         over Rational Field

    To access a basis element, use its degree and index (0 or 1 in the 1st
    cohomology group of a torus)::

        sage: H = simplicial_complexes.Torus().cohomology_ring(QQ)
        sage: H.basis(1)
        Finite family {(1, 0): h^{1,0}, (1, 1): h^{1,1}}
        sage: x = H.basis()[1,0]; x
        h^{1,0}
        sage: y = H.basis()[1,1]; y
        h^{1,1}
        sage: 2*x-3*y
        2*h^{1,0} - 3*h^{1,1}

    You can compute cup products of cohomology classes::

        sage: x.cup_product(y)
        h^{2,0}
        sage: y.cup_product(x)
        -h^{2,0}
        sage: x.cup_product(x)
        0

    This works with simplicial, cubical, and `\Delta`-complexes::

        sage: Klein_c = cubical_complexes.KleinBottle()
        sage: H = Klein_c.cohomology_ring(GF(2))
        sage: x,y = H.basis(1)
        sage: x.cup_product(x)
        h^{2,0}
        sage: x.cup_product(y)
        0
        sage: y.cup_product(y)
        h^{2,0}

        sage: Klein_d = delta_complexes.KleinBottle()
        sage: H = Klein_d.cohomology_ring(GF(2))
        sage: u,v = H.basis(1)
        sage: u.cup_product(u)
        h^{2,0}
        sage: u.cup_product(v)
        0
        sage: v.cup_product(v)
        h^{2,0}

    The basis elements in the simplicial complex case have been chosen
    differently; apply the change of basis `x \mapsto a + b`, `y \mapsto
    b` to see the same product structure. ::

        sage: Klein_s = simplicial_complexes.KleinBottle()
        sage: H = Klein_s.cohomology_ring(GF(2))
        sage: a,b = H.basis(1)
        sage: a.cup_product(a)
        0
        sage: a.cup_product(b)
        h^{2,0}
        sage: (a+b).cup_product(a+b)
        h^{2,0}
        sage: b.cup_product(b)
        h^{2,0}
    """
    def __init__(self, base_ring, cell_complex, cohomology=False, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: H = RP2.homology_with_basis(QQ)
            sage: TestSuite(H).run()
            sage: H = RP2.homology_with_basis(GF(2))
            sage: TestSuite(H).run()
            sage: H = RP2.cohomology_ring(GF(2))
            sage: TestSuite(H).run()
            sage: H = RP2.cohomology_ring(GF(5))
            sage: TestSuite(H).run()
            sage: H = simplicial_complexes.ComplexProjectivePlane().cohomology_ring()
            sage: TestSuite(H).run()
        """
        # phi is the associated chain contraction.
        # M is the homology chain complex.
        phi, M = cell_complex.algebraic_topological_model(base_ring)
        if cohomology:
            phi = phi.dual()
            # We only need the rank of M in each degree, and since
            # we're working over a field, we don't need to dualize M
            # if working with cohomology.
        category = Modules(base_ring).WithBasis().Graded().FiniteDimensional().or_subcategory(category)
        self._contraction = phi
        self._complex = cell_complex
        self._cohomology = cohomology
        self._graded_indices = {deg: range(M.free_module_rank(deg))
                                for deg in range(cell_complex.dimension()+1)}
        indices = [(deg, i) for deg in self._graded_indices
                   for i in self._graded_indices[deg]]
        CombinatorialFreeModule.__init__(self, base_ring, indices, category=category)

    def basis(self, d=None):
        """
        Return (the degree ``d`` homogeneous component of) the basis
        of this graded vector space.

        INPUT:

        - ``d`` -- (optional) the degree

        EXAMPLES::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: H = RP2.homology_with_basis(QQ)
            sage: H.basis()
            Finite family {(0, 0): h_{0,0}}
            sage: H.basis(0)
            Finite family {(0, 0): h_{0,0}}
            sage: H.basis(1)
            Finite family {}
            sage: H.basis(2)
            Finite family {}
        """
        if d is None:
            return Family(self._indices, self.monomial)
        else:
            indices = [(d, i) for i in self._graded_indices.get(d, [])]
            return Family(indices, self.monomial)

    def degree_on_basis(self, i):
        r"""
        Return the degree of the basis element indexed by ``i``.

        EXAMPLES::

            sage: H = simplicial_complexes.Torus().homology_with_basis(GF(7))
            sage: H.degree_on_basis((2,0))
            2
        """
        return i[0]

    def contraction(self):
        r"""
        The chain contraction associated to this homology computation.

        That is, to work with chain representatives of homology
        classes, we need the chain complex `C` associated to the cell
        complex, the chain complex `H` of its homology (with trivial
        differential), chain maps `\pi: C \to H` and `\iota: H \to C`,
        and a chain contraction `\phi` giving a chain homotopy between
        `1_C` and `\iota \circ \pi`.

        OUTPUT: `\phi`

        See :class:`~sage.homology.chain_homotopy.ChainContraction` for information
        about chain contractions, and see
        :func:`~sage.homology.algebraic_topological_model.algebraic_topological_model`
        for the construction of this particular chain contraction `\phi`.

        EXAMPLES::

            sage: H = simplicial_complexes.Simplex(2).homology_with_basis(QQ)
            sage: H.contraction()
            Chain homotopy between:
              Chain complex endomorphism of Chain complex with at most 3 nonzero terms over Rational Field
              and Chain complex endomorphism of Chain complex with at most 3 nonzero terms over Rational Field

        From the chain contraction, one can also recover the maps `\pi`
        and `\iota`::

            sage: phi = H.contraction()
            sage: phi.pi()
            Chain complex morphism:
              From: Chain complex with at most 3 nonzero terms over Rational Field
              To: Chain complex with at most 1 nonzero terms over Rational Field
            sage: phi.iota()
            Chain complex morphism:
              From: Chain complex with at most 1 nonzero terms over Rational Field
              To: Chain complex with at most 3 nonzero terms over Rational Field
        """
        return self._contraction

    def complex(self):
        """
        The cell complex whose homology is being computed.

        EXAMPLES::

            sage: H = simplicial_complexes.Simplex(2).homology_with_basis(QQ)
            sage: H.complex()
            The 2-simplex
        """
        return self._complex

    def _repr_(self):
        """
        EXAMPLES::

            sage: simplicial_complexes.Torus().homology_with_basis(QQ)
            Homology module of Minimal triangulation of the torus
             over Rational Field
        """
        if self._cohomology:
            base = "Cohomology"
        else:
            base = "Homology"
        return base + " module of {} over {}".format(self._complex, self.base_ring())

    def _repr_term(self, i):
        """
        Return ``'h_{i[0],i[1]}'`` for homology, ``'h^{i[0],i[1]}'`` for
        cohomology, for the basis element indexed by ``i``.

        EXAMPLES::

            sage: H = simplicial_complexes.Torus().homology_with_basis(QQ)
            sage: H.basis()[1,0] # indirect doctest
            h_{1,0}
            sage: latex(H.basis()[1,1]) # indirect doctest
            h_{1,1}
            sage: co = simplicial_complexes.KleinBottle().cohomology_ring(GF(2))
            sage: co.basis()[1,0] # indirect doctest
            h^{1,0}

        """
        sym = '^' if self._cohomology else '_'
        return 'h{}{{{},{}}}'.format(sym, i[0], i[1])

    _latex_term = _repr_term

    @cached_method
    def _to_cycle_on_basis(self, i):
        """
        Return the (co)cycle representative of the basis element
        indexed by ``i``.

        .. SEEALSO::

            :meth:`HomologyVectorSpaceWithBasis.Element.to_cocycle`

        EXAMPLES::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: H = S2.homology_with_basis(QQ)
            sage: H._to_cycle_on_basis((2,0))
            -(0, 1, 2) + (0, 1, 3) - (0, 2, 3) + (1, 2, 3)

            sage: S2.cohomology_ring(QQ)._to_cycle_on_basis((2,0))
            \chi_(0, 1, 3)
            sage: S2.cohomology_ring(QQ)._to_cycle_on_basis((0,0))
            \chi_(0,) + \chi_(1,) + \chi_(2,) + \chi_(3,)

            sage: RP3 = simplicial_complexes.RealProjectiveSpace(3)
            sage: H = RP3.cohomology_ring(GF(2))
            sage: H._to_cycle_on_basis((0,0))
            \chi_(1,) + \chi_(2,) + \chi_(3,) + \chi_(4,) + \chi_(5,) + \chi_(6,)
             + \chi_(7,) + \chi_(8,) + \chi_(9,) + \chi_(10,) + \chi_(11,)
            sage: H._to_cycle_on_basis((1,0))
            \chi_(1, 2) + \chi_(1, 3) + \chi_(1, 4) + \chi_(1, 7)
             + \chi_(1, 10) + \chi_(2, 4) + \chi_(2, 6) + \chi_(2, 9)
             + \chi_(2, 10) + \chi_(2, 11) + \chi_(3, 4) + \chi_(3, 5)
             + \chi_(3, 11) + \chi_(4, 8) + \chi_(4, 9) + \chi_(5, 9)
             + \chi_(5, 10) + \chi_(7, 9) + \chi_(8, 10)
            sage: H._to_cycle_on_basis((2,0))
            \chi_(2, 3, 8) + \chi_(2, 7, 8) + \chi_(3, 4, 8) + \chi_(3, 5, 9)
             + \chi_(3, 6, 7) + \chi_(3, 6, 8) + \chi_(3, 6, 10)
             + \chi_(3, 8, 9) + \chi_(3, 9, 10) + \chi_(4, 5, 7)
             + \chi_(4, 5, 9) + \chi_(5, 6, 7) + \chi_(5, 7, 8)
            sage: H._to_cycle_on_basis((3,0))
            \chi_(3, 4, 5, 9)
        """
        vec = self.contraction().iota().in_degree(i[0]).column(i[1])
        chains = self.complex().n_chains(i[0], self.base_ring(),
                                         cochains=self._cohomology)
        return chains.from_vector(vec)

    class Element(CombinatorialFreeModuleElement):
        def to_cycle(self):
            r"""
            (Co)cycle representative of this homogeneous (co)homology class.

            EXAMPLES::

                sage: S2 = simplicial_complexes.Sphere(2)
                sage: H = S2.homology_with_basis(QQ)
                sage: h20 = H.basis()[2,0]; h20
                h_{2,0}
                sage: h20.to_cycle()
                -(0, 1, 2) + (0, 1, 3) - (0, 2, 3) + (1, 2, 3)

            Chains are written as linear combinations of simplices
            `\sigma`. Cochains are written as linear combinations of
            characteristic functions `\chi_{\sigma}` for those
            simplices::

                sage: S2.cohomology_ring(QQ).basis()[2,0].to_cycle()
                \chi_(0, 1, 3)
                sage: S2.cohomology_ring(QQ).basis()[0,0].to_cycle()
                \chi_(0,) + \chi_(1,) + \chi_(2,) + \chi_(3,)
            """
            if not self.is_homogeneous():
                raise ValueError("only defined for homogeneous elements")
            return sum(c * self.parent()._to_cycle_on_basis(i) for i,c in self)

class CohomologyRing(HomologyVectorSpaceWithBasis):
    """
    The cohomology ring.

    .. NOTE::

        This is not intended to be created directly by the user, but
        instead via the
        :meth:`cohomology ring<sage.homology.cell_complex.GenericCellComplex.cohomology_ring>`
        of a :class:`cell
        complex<sage.homology.cell_complex.GenericCellComplex>`.

    INPUT:

    - ``base_ring`` -- must be a field
    - ``cell_complex`` -- the cell complex whose homology we are
      computing

    EXAMPLES::

        sage: CP2 = simplicial_complexes.ComplexProjectivePlane()
        sage: H = CP2.cohomology_ring(QQ)
        sage: H.basis(2)
        Finite family {(2, 0): h^{2,0}}
        sage: x = H.basis(2)[2,0]

    The product structure is the cup product::

        sage: x.cup_product(x)
        h^{4,0}
        sage: x * x
        h^{4,0}

    There are mod 2 cohomology operations defined, also::

        sage: Hmod2 = CP2.cohomology_ring(GF(2))
        sage: y = Hmod2.basis(2)[2,0]
        sage: y.Sq(2)
        h^{4,0}
    """
    def __init__(self, base_ring, cell_complex):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: H = RP2.cohomology_ring(GF(2))
            sage: TestSuite(H).run()
            sage: H = RP2.cohomology_ring(GF(5))
            sage: TestSuite(H).run()
        """
        cat = Algebras(base_ring).WithBasis().Graded().FiniteDimensional()
        HomologyVectorSpaceWithBasis.__init__(self, base_ring, cell_complex, True, cat)

    def _repr_(self):
        """
        EXAMPLES::

            sage: simplicial_complexes.Torus().cohomology_ring(QQ)
            Cohomology ring of Minimal triangulation of the torus
             over Rational Field
        """
        return "Cohomology ring of {} over {}".format(self._complex, self.base_ring())

    @cached_method
    def one(self):
        """
        The multiplicative identity element.

        EXAMPLES::

            sage: H = simplicial_complexes.Torus().cohomology_ring(QQ)
            sage: H.one()
            h^{0,0}
            sage: all(H.one() * x == x == x * H.one() for x in H.basis())
            True
        """
        one = self.base_ring().one()
        d = {(0,i): one for i in self._graded_indices[0]}
        return self._from_dict(d, remove_zeros=False)

    @cached_method
    def product_on_basis(self, li, ri):
        r"""
        The cup product of the basis elements indexed by ``li`` and ``ri``
        in this cohomology ring.

        INPUT:

        - ``li``, ``ri`` -- index of a cohomology class

        .. SEEALSO::

            :meth:`CohomologyRing.Element.cup_product` -- the
            documentation for this method describes the algorithm.

        EXAMPLES::

            sage: RP3 = simplicial_complexes.RealProjectiveSpace(3)
            sage: H = RP3.cohomology_ring(GF(2))
            sage: c = H.basis()[1,0]
            sage: c.cup_product(c).cup_product(c) # indirect doctest
            h^{3,0}

            sage: T = simplicial_complexes.Torus()
            sage: x,y = T.cohomology_ring(QQ).basis(1)
            sage: x.cup_product(y)
            h^{2,0}
            sage: x.cup_product(x)
            0

            sage: one = T.cohomology_ring(QQ).basis()[0,0]
            sage: x.cup_product(one)
            h^{1,0}
            sage: one.cup_product(y) == y
            True
            sage: one.cup_product(one)
            h^{0,0}
            sage: x.cup_product(y) + y.cup_product(x)
            0

        This also works with cubical complexes::

            sage: T = cubical_complexes.Torus()
            sage: x,y = T.cohomology_ring(QQ).basis(1)
            sage: x.cup_product(y)
            -h^{2,0}
            sage: x.cup_product(x)
            0

        and `\Delta`-complexes::

            sage: T_d = delta_complexes.Torus()
            sage: a,b = T_d.cohomology_ring(QQ).basis(1)
            sage: a.cup_product(b)
            h^{2,0}
            sage: b.cup_product(a)
            -h^{2,0}
            sage: RP2 = delta_complexes.RealProjectivePlane()
            sage: w = RP2.cohomology_ring(GF(2)).basis()[1,0]
            sage: w.cup_product(w)
            h^{2,0}

        A non-connected example::

            sage: K = cubical_complexes.Torus().disjoint_union(cubical_complexes.Torus())
            sage: a,b,c,d = K.cohomology_ring(QQ).basis(1)
            sage: x,y = K.cohomology_ring(QQ).basis(0)
            sage: a.cup_product(x) == a
            True
            sage: a.cup_product(y)
            0
        """
        B = self.basis()
        scomplex = self.complex()
        base_ring = self.base_ring()
        deg_left = li[0]
        deg_right = ri[0]
        deg_tot = deg_left + deg_right
        left_cycle = self._to_cycle_on_basis(li)
        right_cycle = self._to_cycle_on_basis(ri)
        n_chains_left = scomplex.n_chains(deg_left, base_ring)
        n_chains_right = scomplex.n_chains(deg_right, base_ring)

        result = {}
        H = scomplex.homology_with_basis(base_ring)
        for gamma_index in H._graded_indices.get(deg_tot, []):
            gamma_coeff = base_ring.zero()
            for cell, coeff in H._to_cycle_on_basis((deg_tot, gamma_index)):
                for (c, left_cell, right_cell) in scomplex.alexander_whitney(cell, deg_left):
                    if c:
                        left = n_chains_left(left_cell)
                        right = n_chains_right(right_cell)
                        gamma_coeff += c * coeff * left_cycle.eval(left) * right_cycle.eval(right)
            if gamma_coeff != base_ring.zero():
                result[(deg_tot, gamma_index)] = gamma_coeff
        return self._from_dict(result, remove_zeros=False)

    class Element(HomologyVectorSpaceWithBasis.Element):
        def cup_product(self, other):
            r"""
            Return the cup product of this element and ``other``.

            Algorithm: see González-Díaz and Réal [G-DR03]_, p. 88.
            Given two cohomology classes, lift them to cocycle
            representatives via the chain contraction for this
            complex, using
            :meth:`~HomologyVectorSpaceWithBasis.Element.to_cycle`. In
            the sum of their dimensions, look at all of the homology
            classes `\gamma`: lift each of those to a cycle
            representative, apply the Alexander-Whitney diagonal map
            to each cell in the cycle, evaluate the two cocycles on
            these factors, and multiply. The result is the value of
            the cup product cocycle on this homology class. After this
            has been done for all homology classes, since homology and
            cohomology are dual, one can tell which cohomology class
            corresponds to the cup product.

            .. SEEALSO::

                :meth:`CohomologyRing.product_on_basis`

            EXAMPLES::

                sage: RP3 = simplicial_complexes.RealProjectiveSpace(3)
                sage: H = RP3.cohomology_ring(GF(2))
                sage: c = H.basis()[1,0]
                sage: c.cup_product(c)
                h^{2,0}
                sage: c * c * c
                h^{3,0}

            We can also take powers::

                sage: RP2 = simplicial_complexes.RealProjectivePlane()
                sage: a = RP2.cohomology_ring(GF(2)).basis()[1,0]
                sage: a**0
                h^{0,0}
                sage: a**1
                h^{1,0}
                sage: a**2
                h^{2,0}
                sage: a**3
                0

            A non-connected example::

                sage: K = cubical_complexes.Torus().disjoint_union(cubical_complexes.Sphere(2))
                sage: a,b = K.cohomology_ring(QQ).basis(2)
                sage: a**0
                h^{0,0} + h^{0,1}

            """
            return self * other

        def Sq(self, i):
            r"""
            Return the result of applying `Sq^i` to this element.

            INPUT:

            - ``i`` -- nonnegative integer

            .. WARNING::

               This is only implemented for simplicial complexes.

            This cohomology operation is only defined in
            characteristic 2.

            Algorithm: see González-Díaz and Réal [G-DR99]_,
            Corollary 3.2.

            EXAMPLES::

                sage: RP2 = simplicial_complexes.RealProjectiveSpace(2)
                sage: x = RP2.cohomology_ring(GF(2)).basis()[1,0]
                sage: x.Sq(1)
                h^{2,0}

                sage: K = RP2.suspension()
                sage: K.set_immutable()
                sage: y = K.cohomology_ring(GF(2)).basis()[2,0]
                sage: y.Sq(1)
                h^{3,0}

                sage: RP4 = simplicial_complexes.RealProjectiveSpace(4)
                sage: H = RP4.cohomology_ring(GF(2))
                sage: x = H.basis()[1,0]
                sage: y = H.basis()[2,0]
                sage: z = H.basis()[3,0]
                sage: x.Sq(1) == y
                True
                sage: z.Sq(1)  # long time
                h^{4,0}

            TESTS::

                sage: T = cubical_complexes.Torus()
                sage: x = T.cohomology_ring(GF(2)).basis()[1,0]
                sage: x.Sq(1)
                Traceback (most recent call last):
                ...
                NotImplementedError: Steenrod squares are only implemented for simplicial complexes
                sage: S2 = simplicial_complexes.Sphere(2)
                sage: x = S2.cohomology_ring(GF(7)).basis()[2,0]
                sage: x.Sq(1)
                Traceback (most recent call last):
                ...
                ValueError: Steenrod squares are only defined in characteristic 2
            """
            P = self.parent()
            scomplex = P.complex()
            if not isinstance(scomplex, SimplicialComplex):
                raise NotImplementedError('Steenrod squares are only implemented for simplicial complexes')
            base_ring = P.base_ring()
            if base_ring.characteristic() != 2:
                raise ValueError('Steenrod squares are only defined in characteristic 2')
            # We keep the same notation as in [G-DR99].
            # The trivial cases:
            if i == 0:
                # Sq^0 is the identity.
                return self

            # Construct each graded component of ``self``
            ret = P.zero()
            H = scomplex.homology_with_basis(base_ring)
            deg_comp = {}
            for index,coeff in self:
                d = deg_comp.get(index[0], {})
                d[index] = coeff
                deg_comp[index[0]] = d

            # Do the square on each graded componenet of ``self``.
            for j in deg_comp:
                # Make it into an actual element
                m = j + i
                if not P._graded_indices.get(m, []) or i > j:
                    continue
                elt = P._from_dict(deg_comp[j], remove_zeros=False)
                if i == j:
                    ret += elt.cup_product(elt)
                    continue

                n = j - i
                # Now assemble the indices over which the sums take place.
                # S(n) is defined to be floor((m+1)/2) + floor(n/2).
                S_n = (m+1) // 2 + n // 2
                if n == 0:
                    sums = [[S_n]]
                else:
                    sums = [[i_n] + l for i_n in range(S_n, m+1)
                            for l in sum_indices(n-1, i_n, S_n)]
                # At this point, 'sums' is a list of lists of the form
                # [i_n, i_{n-1}, ..., i_0]. (It is reversed from the
                # obvious order because this is closer to the order in
                # which the face maps will be applied.)  Now we sum over
                # these, according to the formula in [G-DR99], Corollary 3.2.
                result = {}
                cycle = elt.to_cycle()
                n_chains = scomplex.n_chains(j, base_ring)
                for gamma_index in H._graded_indices.get(m, []):
                    gamma_coeff = base_ring.zero()
                    for cell, coeff in H._to_cycle_on_basis((m, gamma_index)):
                        for indices in sums:
                            indices = list(indices)
                            left = cell
                            right = cell
                            # Since we are working with a simplicial complex, 'cell' is a simplex.
                            if not m % 2:
                                left_endpoint = m
                                while indices:
                                    right_endpoint = indices[0] - 1
                                    for k in range(left_endpoint, indices.pop(0), -1):
                                        left = scomplex.face(left, k)
                                    try:
                                        left_endpoint = indices[0] - 1
                                        for k in range(right_endpoint, indices.pop(0), -1):
                                            right = scomplex.face(right, k)
                                    except IndexError:
                                        pass
                                for k in range(right_endpoint, -1, -1):
                                    right = scomplex.face(right, k)
                            else:
                                right_endpoint = m
                                while indices:
                                    left_endpoint = indices[0] - 1
                                    try:
                                        for k in range(right_endpoint, indices.pop(0), -1):
                                            right = scomplex.face(right, k)
                                        right_endpoint = indices[0] - 1
                                    except IndexError:
                                        pass
                                    for k in range(left_endpoint, indices.pop(0), -1):
                                        left = scomplex.face(left, k)
                                for k in range(right_endpoint, -1, -1):
                                    right = scomplex.face(right, k)

                            left = n_chains(left)
                            right = n_chains(right)
                            gamma_coeff += coeff * cycle.eval(left) * cycle.eval(right)
                    if gamma_coeff != base_ring.zero():
                        result[(m, gamma_index)] = gamma_coeff
                ret += P._from_dict(result, remove_zeros=False)
            return ret

def sum_indices(k, i_k_plus_one, S_k_plus_one):
    r"""
    This is a recursive function for computing the indices for the
    nested sums in González-Díaz and Réal [G-DR99]_, Corollary 3.2.

    In the paper, given indices `i_n`, `i_{n-1}`, ..., `i_{k+1}`,
    given `k`, and given `S(k+1)`, the number `S(k)` is defined to be

    .. MATH::

        S(k) = -S(k+1) + floor(k/2) + floor((k+1)/2) + i_{k+1},

    and `i_k` ranges from `S(k)` to `i_{k+1}-1`. There are two special
    cases: if `k=0`, then `i_0 = S(0)`. Also, the initial case of
    `S(k)` is `S(n)`, which is set in the method :meth:`Sq` before
    calling this function. For this function, given `k`, `i_{k+1}`,
    and `S(k+1)`, return a list consisting of the allowable possible
    indices `[i_k, i_{k-1}, ..., i_1, i_0]` given by the above
    formula.

    INPUT:

    - ``k`` -- non-negative integer
    - ``i_k_plus_one`` -- the positive integer `i_{k+1}`
    - ``S_k_plus_one`` -- the integer `S(k+1)`

    EXAMPLES::

        sage: from sage.homology.homology_vector_space_with_basis import sum_indices
        sage: sum_indices(1, 3, 3)
        [[1, 0], [2, 1]]
        sage: sum_indices(0, 4, 2)
        [[2]]
    """
    S_k = -S_k_plus_one + k//2 + (k+1)//2 + i_k_plus_one
    if k == 0:
        return [[S_k]]
    return [[i_k] + l for i_k in range(S_k, i_k_plus_one)
            for l in sum_indices(k-1, i_k, S_k)]

