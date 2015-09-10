r"""
Induced morphisms on homology

This module implements morphisms on homology induced by morphisms of
simplicial complexes. It requires working with field coefficients.

See :class:`InducedHomologyMorphism` for documentation.

AUTHORS:

- John H. Palmieri (2015.09)
"""

########################################################################
#       Copyright (C) 2015 John H. Palmieri <palmieri@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

# To do: implement morphisms of cubical complexes, with methods
#   - domain
#   - codomain
#   - associated_chain_complex_morphism
# Once this is done, the code here ought to work without modification.

from sage.structure.sage_object import SageObject
from sage.rings.rational_field import QQ

class InducedHomologyMorphism(SageObject):
    """
    An element of this class is a morphism of (co)homology groups
    induced by a map of simplicial complexes. It requires working
    with field coefficients.

    INPUTS:

    - ``map`` -- the map of simplicial complexes
    - ``base_ring`` -- a field (optional, default ``QQ``)
    - ``cohomology`` -- boolean (optional, default ``False``). If
      ``True``, return the induced map in cohomology rather than
      homology.

    .. note::

        This is not intended to be used directly by the user, but instead
        via the method
        :meth:`simplicial_complex_morphism.SimplicialComplexMorphism.induced_homology_morphism`.

    EXAMPLES::

        sage: S1 = simplicial_complexes.Sphere(1)
        sage: H = Hom(S1, S1)
        sage: f = H({0:0, 1:2, 2:1})  # f switches two vertices
        sage: f_star = f.induced_homology_morphism(QQ, cohomology=True)
        sage: f_star.to_matrix(1)
        [-1]

        sage: T = simplicial_complexes.Torus()
        sage: y = T.homology_with_basis(1, QQ).basis()[1]
        sage: y.to_cycle()
        (0, 3) - (0, 6) + (3, 6)

    Since `(0,3) - (0,6) + (3,6)` is a cycle representing a homology
    class in the torus, we can define a map `S^1 \to T` inducing an
    inclusion on `H_1`::

        sage: Hom(S1, T)({0:0, 1:3, 2: 6})
        Simplicial complex morphism {0: 0, 1: 3, 2: 6} from Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)} to Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 14 facets
        sage: g = Hom(S1, T)({0:0, 1:3, 2: 6})
        sage: g_star = g.induced_homology_morphism(QQ)
        sage: g_star.to_matrix(0)
        [1]
        sage: g_star.to_matrix(1)
        [0]
        [1]
    """
    def __init__(self, map, base_ring=None, cohomology=False):
        """
        INPUTS:

        - ``map`` -- the map of simplicial complexes
        - ``base_ring`` -- a field (optional, default ``QQ``)
        - ``cohomology`` -- boolean (optional, default ``False``). If
          ``True``, return the induced map in cohomology rather than
          homology.

        EXAMPLES::

            sage: from sage.homology.homology_morphism import InducedHomologyMorphism
            sage: K = simplicial_complexes.RandomComplex(8, 3)
            sage: H = Hom(K,K)
            sage: id = H.identity()
            sage: f = InducedHomologyMorphism(id, QQ)
            sage: f.to_matrix(0) == 1  and  f.to_matrix(1) == 1  and  f.to_matrix(2) == 1
            True
            sage: f = InducedHomologyMorphism(id, ZZ)
            Traceback (most recent call last):
            ...
            ValueError: the coefficient ring must be a field
            sage: S1 = simplicial_complexes.Sphere(1).barycentric_subdivision()
            sage: S1.is_mutable()
            True
            sage: g = Hom(S1, S1).identity()
            sage: h = g.induced_homology_morphism(QQ)
            Traceback (most recent call last):
            ...
            ValueError: the domain and codomain complexes must be immutable
            sage: S1.set_immutable()
            sage: g = Hom(S1, S1).identity()
            sage: h = g.induced_homology_morphism(QQ)
        """
        if map.domain().is_mutable() or map.codomain().is_mutable():
            raise ValueError('the domain and codomain complexes must be immutable')
        if base_ring is None:
            base_ring = QQ
        if not base_ring.is_field():
            raise ValueError('the coefficient ring must be a field')

        self._cohomology = cohomology
        self._map = map
        self._base_ring = base_ring

    def base_ring(self):
        """
        The base ring for this map

        EXAMPLES::

            sage: K = simplicial_complexes.Simplex(2)
            sage: H = Hom(K,K)
            sage: id = H.identity()
            sage: id.induced_homology_morphism(QQ).base_ring()
            Rational Field
            sage: id.induced_homology_morphism(GF(13)).base_ring()
            Finite Field of size 13
        """
        return self._base_ring

    def to_matrix(self, deg):
        """
        The matrix for this map in degree ``deg``

        INPUTS:

        - ``deg`` -- integer

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S1_b = S1.barycentric_subdivision()
            sage: S1_b.set_immutable()
            sage: d = {(0,): 0, (0,1): 1, (1,): 2, (1,2): 0, (2,): 1, (0,2): 2}
            sage: f = Hom(S1_b, S1)(d)
            sage: h = f.induced_homology_morphism(QQ)
            sage: h.to_matrix(1)
            [2]
        """
        base_ring = self.base_ring()
        if self._cohomology:
            domain = self._map.codomain()
            codomain = self._map.domain()
        else:
            domain = self._map.domain()
            codomain = self._map.codomain()
        phi_codomain, _ = codomain.algebraic_topological_model(base_ring)
        phi_domain, _ = domain.algebraic_topological_model(base_ring)
        return phi_codomain.pi().in_degree(deg) * self._map.associated_chain_complex_morphism(self.base_ring(), cochain=self._cohomology).in_degree(deg) * phi_domain.iota().in_degree(deg)

    def __call__(self, elt):
        """
        Evaluate this map on ``elt``, an element of (co)homology.

        INPUT:

        - ``elt`` -- informally, an element of the domain of this
          map. More formally, an element of
          :class:`homology_vector_space_with_basis.HomologyVectorSpaceWithBasis`.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: f = {0:0, 1:2, 2:1}
            sage: H = Hom(S1,S1)
            sage: g = H(f)
            sage: h = g.induced_homology_morphism(QQ)
            sage: x = S1.homology_with_basis(1).basis()[0]
            sage: x
            h_{1,0}
            sage: h(x)  # indirect doctest
            -h_{1,0}
        """
        deg = elt.parent().degree()
        base_ring = self.base_ring()
        if self._cohomology:
            codomain = self._map.domain().cohomology_with_basis(deg, base_ring)         
            if elt.parent().complex() != self._map.codomain():
                raise ValueError('element is not a cohomology class for the correct complex')
        else:
            codomain = self._map.codomain().homology_with_basis(deg, base_ring)
            if elt.parent().complex() != self._map.domain():
                raise ValueError('element is not a homology class for the correct complex')

        return codomain.from_vector(self.to_matrix(deg) * elt.to_vector())

    def __eq__(self, other):
        """
        Return ``True`` if and only if this map agrees with ``other``.

        INPUTS:

        - ``other`` -- another induced homology morphism

        This automatically returns ``False`` if the morphisms have
        different domains, codomains, base rings, or values for their
        cohomology flags

        Otherwise, determine this by computing the matrices for this
        map and ``other`` using the (same) basis for the homology
        vector spaces.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: K = simplicial_complexes.Simplex(2)
            sage: f = Hom(S1, K)({0: 0, 1:1, 2:2})
            sage: g = Hom(S1, K)({0: 0, 1:0, 2:0})
            sage: f.induced_homology_morphism(QQ) == g.induced_homology_morphism(QQ)
            True
            sage: f.induced_homology_morphism(QQ) == g.induced_homology_morphism(GF(2))
            False
            sage: id = Hom(K, K).identity()   # different domain
            sage: f.induced_homology_morphism(QQ) == id.induced_homology_morphism(QQ)
            False
        """
        if (self._map.domain() != other._map.domain()
            or self._map.codomain() != other._map.codomain()
            or self.base_ring() != other.base_ring()
            or self._cohomology != other._cohomology):
            return False
        dim = min(self._map.domain().dimension(), self._map.codomain().dimension())
        return all(self.to_matrix(d) == other.to_matrix(d) for d in range(dim+1))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: K = simplicial_complexes.Simplex(1)
            sage: f = Hom(K,K).identity()
            sage: f.induced_homology_morphism(QQ)
            Homology morphism induced by Simplicial complex morphism {0: 0, 1: 1} from Simplicial complex with vertex set (0, 1) and facets {(0, 1)} to Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
            sage: f.induced_homology_morphism(QQ, cohomology=True)
            Cohomology morphism induced by Simplicial complex morphism {0: 0, 1: 1} from Simplicial complex with vertex set (0, 1) and facets {(0, 1)} to Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
        """
        name = "Homology" if not self._cohomology else "Cohomology"
        return "{} morphism induced by {}".format(name, self._map)
