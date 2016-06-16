"""
Simplicial Sets
"""
#*****************************************************************************
#  Copyright (C) 2015 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.sets_cat import Sets

class SimplicialSets(Category_singleton):
    r"""
    The category of simplicial sets.

    A simplicial set `X` is a collection of sets `X_i`, indexed by
    the non-negative integers, together with maps

    .. math::

        d_i: X_n \to X_{n-1}, \quad 0 \leq i \leq n \quad \text{(face maps)} \\
        s_j: X_n \to X_{n+1}, \quad 0 \leq j \leq n \quad \text{(degeneracy maps)}

    satisfying the *simplicial identities*:

    .. math::

        d_i d_j &= d_{j-1} d_i \quad \text{if } i<j \\
        d_i s_j &= s_{j-1} d_i \quad \text{if } i<j \\
        d_j s_j &= 1 = d_{j+1} s_j \\
        d_i s_j &= s_{j} d_{i-1} \quad \text{if } i>j+1 \\
        s_i s_j &= s_{j+1} s_{i} \quad \text{if } i \leq j

    Morphisms are sequences of maps `f_i : X_i \to Y_i` which commute
    with the face and degeneracy maps.

    EXAMPLES::

        sage: from sage.categories.simplicial_sets import SimplicialSets
        sage: C = SimplicialSets(); C
        Category of simplicial sets

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.simplicial_sets import SimplicialSets
            sage: SimplicialSets().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class ParentMethods:
        def is_finite(self):
            """
            True if this simplicial set is finite, i.e., has a finite number
            of nondegenerate simplices.

            As of this writing, all simplicial sets implemented in
            Sage are finite.

            EXAMPLES::

                sage: simplicial_sets.Torus().is_finite()
                True
            """
            return SimplicialSets.Finite() in self.categories()

        def is_pointed(self):
            """
            True if this simplicial set is pointed, i.e., has a base point.

            EXAMPLES::

                sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
                sage: v = AbstractSimplex(0)
                sage: w = AbstractSimplex(0)
                sage: e = AbstractSimplex(1)
                sage: X = SimplicialSet({e: (v, w)})
                sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
                sage: X.is_pointed()
                False
                sage: Y.is_pointed()
                True
            """
            return SimplicialSets.Pointed() in self.categories()

        def set_base_point(self, point):
            """
            Return a copy of this simplicial set in which the base point is
            set to ``point``.

            INPUT:

            - ``point`` -- a 0-simplex in this simplicial set

            EXAMPLES::

                sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
                sage: v = AbstractSimplex(0, name='v_0')
                sage: w = AbstractSimplex(0, name='w_0')
                sage: e = AbstractSimplex(1)
                sage: X = SimplicialSet({e: (v, w)})
                sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
                sage: Y.base_point()
                w_0
                sage: X_star = X.set_base_point(w)
                sage: X_star.base_point()
                w_0
                sage: Y_star = Y.set_base_point(v)
                sage: Y_star.base_point()
                v_0

            TESTS::

                sage: X.set_base_point(e)
                Traceback (most recent call last):
                ...
                ValueError: the "point" is not a zero-simplex
                sage: pt = AbstractSimplex(0)
                sage: X.set_base_point(pt)
                Traceback (most recent call last):
                ...
                ValueError: the point is not a simplex in this simplicial set
            """
            from sage.homology.simplicial_set import SimplicialSet
            if point.dimension() != 0:
                raise ValueError('the "point" is not a zero-simplex')
            if point not in self._simplices:
                raise ValueError('the point is not a simplex in this '
                                 'simplicial set')
            return SimplicialSet(self.face_data(), base_point=point)

    class ElementMethods:
        pass

    class MorphismMethods:
        pass

    class Finite(CategoryWithAxiom):
        """
        Category of finite simplicial sets.

        The objects are simplicial sets with finitely many
        non-degenerate simplices.
        """
        class ParentMethods:
            pass

        class ElementMethods:
            pass

        class MorphismMethods:
            pass

    class SubcategoryMethods:
        def Pointed(self):
            """
            A simplicial set is *pointed* if it has a distinguished base
            point.

            EXAMPLES::

                sage: from sage.categories.simplicial_sets import SimplicialSets
                sage: SimplicialSets().Pointed().Finite()
                Category of finite pointed simplicial sets
                sage: SimplicialSets().Finite().Pointed()
                Category of finite pointed simplicial sets
            """
            return self._with_axiom("Pointed")

    class Pointed(CategoryWithAxiom):
        class ParentMethods:
            def base_point(self):
                """
                Return this simplicial set's base point

                EXAMPLES::

                    sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
                    sage: v = AbstractSimplex(0, name='*')
                    sage: e = AbstractSimplex(1)
                    sage: S1 = SimplicialSet({e: (v, v)}, base_point=v)
                    sage: S1.is_pointed()
                    True
                    sage: S1.base_point()
                    *
                """
                return self._basepoint

            def base_point_map(self, domain=None):
                """
                Return a map from a one-point space to this one, with image the
                base point.

                This raises an error if this simplicial set does not have a
                base point.

                INPUT:

                - ``domain`` -- optional, default ``None``. Use this
                  to specify a particular one-point space as the
                  domain. The default behavior is to use the
                  :func:`sage.homology.simplicial_set.Point` function
                  to use a standard one-point space.

                EXAMPLES::

                    sage: T = simplicial_sets.Torus()
                    sage: f = T.base_point_map(); f
                    Simplicial set morphism:
                      From: Point
                      To:   Torus
                      Defn: [*] --> [(v_0, v_0)]
                    sage: S3 = simplicial_sets.Sphere(3)
                    sage: g = S3.base_point_map()
                    sage: f.domain() == g.domain()
                    True
                    sage: RP3 = simplicial_sets.RealProjectiveSpace(3)
                    sage: temp = simplicial_sets.Simplex(0)
                    sage: pt = temp.set_base_point(temp.n_cells(0)[0])
                    sage: h = RP3.base_point_map(domain=pt)
                    sage: f.domain() == h.domain()
                    False
                """
                from sage.homology.simplicial_set import Point
                if domain is None:
                    domain = Point()
                else:
                    if len(domain._simplices) > 1:
                        raise ValueError('domain has more than one nondegenerate simplex')
                src = domain.base_point()
                target = self.base_point()
                return domain.Hom(self)({src: target})

            def unset_base_point(self):
                """
                Return a copy of this simplicial set in which the base point has
                been forgotten.

                EXAMPLES::

                    sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
                    sage: v = AbstractSimplex(0, name='v_0')
                    sage: w = AbstractSimplex(0, name='w_0')
                    sage: e = AbstractSimplex(1)
                    sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
                    sage: Y.is_pointed()
                    True
                    sage: Y.base_point()
                    w_0
                    sage: Z = Y.unset_base_point()
                    sage: Z.is_pointed()
                    False
                """
                from sage.homology.simplicial_set import SimplicialSet
                return SimplicialSet(self.face_data())

            def fat_wedge(self, n):
                """
                The $n$-th fat wedge of this pointed simplicial set.

                This is the subcomplex of the $n$-fold product `X^n`
                consisting of those points in which at least one
                factor is the base point. Thus when $n=2$, this is the
                wedge of the simplicial set with itself, but when $n$
                is larger, the fat wedge is larger than the $n$-fold
                wedge.

                EXAMPLES::

                    sage: S1 = simplicial_sets.Sphere(1)
                    sage: S1.fat_wedge(0)
                    Point
                    sage: S1.fat_wedge(1)
                    S^1
                    sage: S1.fat_wedge(2).fundamental_group()
                    Finitely presented group < e0, e1 |  >
                    sage: S1.fat_wedge(4).homology()
                    {0: 0, 1: Z x Z x Z x Z, 2: Z^6, 3: Z x Z x Z x Z}
                """
                from sage.homology.simplicial_set import Point
                if n == 0:
                    return Point()
                if n == 1:
                    return self
                return self.product(*[self]*(n-1)).fat_wedge()

            def smash_product(self, other):
                """
                The smash product of this simplicial set with ``other``.

                EXAMPLES::

                    sage: S1 = simplicial_sets.Sphere(1)
                    sage: RP2 = simplicial_sets.RealProjectiveSpace(2)
                    sage: X = S1.smash_product(RP2)
                    sage: X.homology(base_ring=GF(2))
                    {0: Vector space of dimension 0 over Finite Field of size 2,
                     1: Vector space of dimension 0 over Finite Field of size 2,
                     2: Vector space of dimension 1 over Finite Field of size 2,
                     3: Vector space of dimension 1 over Finite Field of size 2}

                    sage: T = S1.product(S1)
                    sage: X = T.smash_product(S1)
                    sage: X.homology(reduced=False)
                    {0: Z, 1: 0, 2: Z x Z, 3: Z}
                """
                prod = self.product(other)
                wedge = prod.wedge_as_subset()
                return prod.quotient(wedge)




        class ElementMethods:
            pass

        class MorphismMethods:
            pass

        
