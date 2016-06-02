r"""
Morphisms and homsets for simplicial sets

AUTHORS:

- John H. Palmieri (2016-03)

This module implements morphisms and homsets of simplicial sets.
"""

#*****************************************************************************
#  Copyright (C) 2016 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

from sage.homology.simplicial_set import AbstractSimplex, NonDegenerateSimplex, SimplicialSet
from sage.matrix.constructor import matrix, zero_matrix
from sage.rings.integer_ring import ZZ
from sage.homology.chain_complex_morphism import ChainComplexMorphism
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom, Homset
from sage.categories.simplicial_sets import SimplicialSets

class SimplicialSetHomset(Homset):
    r"""
    Set of morphisms between simplicial sets.

    Given a hom set, one can use it to construct a morphism `f` by
    specifying a dictionary, the keys of which are the nondegenerate
    simplices in the domain, and the value corresponding to `\sigma`
    is the simplex `f(\sigma)` in the codomain.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
        sage: v = NonDegenerateSimplex(0, name='v')
        sage: w = NonDegenerateSimplex(0, name='w')
        sage: e = NonDegenerateSimplex(1, name='e')
        sage: f = NonDegenerateSimplex(1, name='f')
        sage: X = SimplicialSet({e: (v, w), f: (w, v)})
        sage: Y = SimplicialSet({e: (v, v)})

    Define the hom set::

        sage: H = Hom(X, Y)

    Now define a morphism by specifying a dictionary::

        sage: H({v: v, w: v, e: e, f: e})
        Simplicial set morphism:
          From: Simplicial set with 4 non-degenerate simplices
          To:   Simplicial set with 2 non-degenerate simplices
          Defn: [v, w, e, f] --> [v, v, e, e]
    """
    def __call__(self, f):
        r"""
        INPUT:

        - ``f`` -- a dictionary with keys the simplices of the domain
          and values simplices of the codomain

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = S1.n_cells(0)[0]
            sage: e = S1.n_cells(1)[0]
            sage: f = {v0: v0, e: v0.apply_degeneracies(0)} # constant map
            sage: Hom(S1, S1)(f)
            Simplicial set endomorphism of S^1
              Defn: [v_0, sigma_1] --> [v_0, Simplex obtained by applying degeneracy s_0 to v_0]
        """
        return SimplicialSetMorphism(f, self.domain(), self.codomain())

    def diagonal_morphism(self):
        r"""
        Return the diagonal morphism in `Hom(S, S \times S)`.

        EXAMPLES::

            sage: RP2 = simplicial_sets.RealProjectiveSpace(2)
            sage: Hom(RP2, RP2.product(RP2)).diagonal_morphism()
            Simplicial set morphism:
              From: RP^2
              To:   Simplicial set with 31 non-degenerate simplices
              Defn: [1, f, f * f] --> [(1, 1), (f, f), (f * f, f * f)]
        """
        domain = self.domain()
        codomain = self.codomain()
        if not hasattr(codomain, 'factors'):
            raise ValueError('diagonal morphism is only defined for Hom(X, XxX)')
        factors = codomain.factors()
        if len(factors) != 2 or factors[0] != domain or factors[1] != domain:
            raise ValueError('diagonal morphism is only defined for Hom(X, XxX)')
        f = {}
        for i in range(domain.dimension()+1):
            for s in domain.n_cells(i):
                f[s] = dict(codomain._translation)[((s, ()), (s, ()))]
        return self(f)

    def identity(self):
        """
        Return the identity morphism in `Hom(S, S)`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: Hom(S1, S1).identity()
            Simplicial set endomorphism of S^1
              Defn: [v_0, sigma_1] --> [v_0, sigma_1]
            sage: T = simplicial_sets.Torus()
            sage: Hom(S1, T).identity()
            Traceback (most recent call last):
            ...
            TypeError: identity map is only defined for endomorphism sets
        """
        if not self.is_endomorphism_set():
            raise TypeError("identity map is only defined for endomorphism sets")
        f = {}
        domain = self.domain()
        for i in range(domain.dimension()+1):
            for s in domain.n_cells(i):
                f[s] = s
        return self(f)

    def an_element(self):
        """
        Return an element of ``self``: a constant map.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: Hom(S2, S1).an_element()
            Simplicial set morphism:
              From: S^2
              To:   S^1
              Defn: [v_0, sigma_2] --> [v_0, Simplex obtained by applying degeneracies s_1 s_0 to v_0]
        """
        domain = self.domain()
        codomain = self.codomain()
        target = codomain.n_cells(0)[0]
        const = {}
        for i in range(domain.dimension() + 1):
            for s in domain.n_cells(i):
                const[s] = target
            target = target.apply_degeneracies(i)
        return self(const)


class SimplicialSetMorphism(Morphism):
    def __init__(self, f, domain, codomain):
        """
        A morphism of simplicial sets.

        INPUTS:

        - ``f`` -- dictionary defining the map
        - ``domain`` -- simplicial set
        - ``codomain`` -- simplicial set

        The keys of the dictionary are the nondegenerate simplices of
        the domain, the corresponding values are simplices in the
        codomain.

        EXAMPLES::

            sage: from sage.homology.simplicial_set_morphism import SimplicialSetMorphism
            sage: K = simplicial_sets.Simplex(1)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = K.n_cells(0)[0]
            sage: v1 = K.n_cells(0)[1]
            sage: e01 = K.n_cells(1)[0]
            sage: w = S1.n_cells(0)[0]
            sage: sigma = S1.n_cells(1)[0]

            sage: f = {v0: w, v1: w, e01: sigma}
            sage: SimplicialSetMorphism(f, K, S1)
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: [(0,), (1,), (0, 1)] --> [v_0, v_0, sigma_1]

        The same map can be defined as follows::

            sage: H = Hom(K, S1)
            sage: H(f)
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: [(0,), (1,), (0, 1)] --> [v_0, v_0, sigma_1]

        A constant map::

            sage: g = {v0: w, v1: w, e01: w.apply_degeneracies(0)}
            sage: SimplicialSetMorphism(g, K, S1)
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: [(0,), (1,), (0, 1)] --> [v_0, v_0, Simplex obtained by applying degeneracy s_0 to v_0]

        A non-map::

            sage: h = {w: v0, sigma: e01}
            sage: SimplicialSetMorphism(h, S1, K)
            Traceback (most recent call last):
            ...
            ValueError: the dictionary does not define a map of simplicial sets
        """
        if not isinstance(domain, SimplicialSet) or not isinstance(codomain, SimplicialSet):
            raise ValueError('the domain and codomain must be simplicial sets')

        # We have to check that the proposed map commutes with the
        # face maps. (The degeneracy maps should work automatically,
        # since our simplicial sets have "free" degeneracies.)
        for simplex in f:
            # Compare f[d_i (simplex)] to d_i f[simplex]. Since
            # d_i(simplex) may be degenerate, we have to be careful
            # when applying f to it. We can skip vertices and start
            # with 1-simplices.
            bad = False
            for i in range(simplex.dimension()+1):
                face_f = codomain.face(f[simplex], i)
                face = domain.face(simplex, i)
                if face is None:
                    f_face = None
                elif face.is_nondegenerate():
                    f_face = f[face]
                else:
                    f_face = f[face.nondegenerate()].apply_degeneracies(*face.degeneracies())
                if face_f != f_face:
                    bad = True
                    break
            if bad:
                raise ValueError('the dictionary does not define a map of simplicial sets')
        self._dictionary = f
        Morphism.__init__(self, Hom(domain, codomain, SimplicialSets()))

    def __call__(self, x):
        """
        INPUT: a simplex of the domain.

        Return its image under this morphism.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(1)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = K.n_cells(0)[0]
            sage: v1 = K.n_cells(0)[1]
            sage: e01 = K.n_cells(1)[0]
            sage: w = S1.n_cells(0)[0]
            sage: sigma = S1.n_cells(1)[0]
            sage: d = {v0: w, v1: w, e01: sigma}
            sage: f = Hom(K, S1)(d)
            sage: f(e01) # indirect doctest
            sigma_1
        """
        try:
            return self._dictionary[x]
        except KeyError:
            raise ValueError('element is not a simplex in the domain')

    def _repr_type(self):
        """
        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: f = Hom(S1,S1).identity()
            sage: f._repr_type()
            'Simplicial set'
        """
        return "Simplicial set"

    def _repr_defn(self):
        """
        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: f = Hom(S1,S1).identity()
            sage: f._repr_defn()
            '[v_0, sigma_1] --> [v_0, sigma_1]'
        """
        d = self._dictionary
        keys = sorted(d.keys())
        return "{} --> {}".format(keys, [d[x] for x in keys])

    def is_pointed(self):
        """
        Return ``True`` if this is a pointed map.

        That is, return ``True`` if the domain and codomain are
        pointed and this morphism preserves the base point.

        EXAMPLES::

        



        """
        return (self.domain().is_pointed() and self.codomain().is_pointed() 
                and self(self.domain().base_point()) == self.codomain().base_point())

    def associated_chain_complex_morphism(self, base_ring=ZZ,
                                          augmented=False, cochain=False):
        """
        Return the associated chain complex morphism of ``self``.

        INPUT:

        - ``base_ring`` -- default ``ZZ``
        - ``augmented`` -- boolean, default ``False``. If ``True``,
          return the augmented complex.
        - ``cochain`` -- boolean, default ``False``. If ``True``,
          return the cochain complex.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = S1.n_cells(0)[0]
            sage: e = S1.n_cells(1)[0]
            sage: f = {v0: v0, e: v0.apply_degeneracies(0)} # constant map
            sage: g = Hom(S1, S1)(f)
            sage: g.associated_chain_complex_morphism().to_matrix()
            [1|0]
            [-+-]
            [0|0]
        """
        # One or the other chain complex is trivial between these
        # dimensions:
        max_dim = max(self.domain().dimension(), self.codomain().dimension())
        min_dim = min(self.domain().dimension(), self.codomain().dimension())
        matrices = {}
        if augmented is True:
            m = matrix(base_ring,1,1,1)
            if not cochain:
                matrices[-1] = m
            else:
                matrices[-1] = m.transpose()
        for dim in range(min_dim+1):
            X_faces = list(self.domain().n_cells(dim))
            Y_faces = list(self.codomain().n_cells(dim))
            num_faces_X = len(X_faces)
            num_faces_Y = len(Y_faces)
            mval = [0 for i in range(num_faces_X * num_faces_Y)]
            for idx,x in enumerate(X_faces):
                y = self(x)
                if y.is_nondegenerate():
                    mval[idx + (Y_faces.index(y) * num_faces_X)] = 1
            m = matrix(base_ring, num_faces_Y, num_faces_X, mval, sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        for dim in range(min_dim+1,max_dim+1):
            try:
                l1 = len(self.codomain().n_cells(dim))
            except KeyError:
                l1 = 0
            try:
                l2 = len(self.domain().n_cells(dim))
            except KeyError:
                l2 = 0
            m = zero_matrix(base_ring,l1,l2,sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        if not cochain:
            return ChainComplexMorphism(matrices,
                    self.domain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=False),
                    self.codomain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=False))
        else:
            return ChainComplexMorphism(matrices,
                    self.codomain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=True),
                    self.domain().chain_complex(base_ring=base_ring, augmented=augmented, cochain=True))

    def induced_homology_morphism(self, base_ring=None, cohomology=False):
        """
        The map in (co)homology induced by this map

        INPUT:

        - ``base_ring`` -- must be a field (optional, default ``QQ``)

        - ``cohomology`` -- boolean (optional, default ``False``). If
          ``True``, the map induced in cohomology rather than homology.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: e = NonDegenerateSimplex(1, name='e')
            sage: f = NonDegenerateSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, w), f: (w, v)})
            sage: Y = SimplicialSet({e: (v, v)})
            sage: H = Hom(X, Y)
            sage: f = H({v: v, w: v, e: e, f: e})
            sage: g = f.induced_homology_morphism()
            sage: g.to_matrix()
            [1|0]
            [-+-]
            [0|2]
            sage: g3 = f.induced_homology_morphism(base_ring=GF(3), cohomology=True)
            sage: g3.to_matrix()
            [2|0]
            [-+-]
            [0|1]
        """
        from homology_morphism import InducedHomologyMorphism
        return InducedHomologyMorphism(self, base_ring, cohomology)
