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

from sage.homology.simplicial_set import SimplicialSet_infinite, \
    PushoutOfSimplicialSets, PullbackOfSimplicialSets
from sage.matrix.constructor import matrix, zero_matrix
from sage.rings.integer_ring import ZZ
from sage.homology.chain_complex_morphism import ChainComplexMorphism
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom, Homset
from sage.categories.simplicial_sets import SimplicialSets
from sage.homology.homology_morphism import InducedHomologyMorphism

class SimplicialSetHomset(Homset):
    r"""
    Set of morphisms between simplicial sets.

    Given a hom set, one can use it to construct a morphism `f` by
    specifying a dictionary, the keys of which are the nondegenerate
    simplices in the domain, and the value corresponding to `\sigma`
    is the simplex `f(\sigma)` in the codomain.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: v = AbstractSimplex(0, name='v')
        sage: w = AbstractSimplex(0, name='w')
        sage: e = AbstractSimplex(1, name='e')
        sage: f = AbstractSimplex(1, name='f')
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
        r"""
        A morphism of simplicial sets.

        INPUTS:

        - ``f`` -- dictionary defining the map
        - ``domain`` -- simplicial set
        - ``codomain`` -- simplicial set

        The keys of the dictionary are the nondegenerate simplices of
        the domain, the corresponding values are simplices in the
        codomain.

        Actually, the keys do not need to include all of the
        nondegenerate simplices: if `\sigma` is a face of `\tau`, then
        only the image of `\tau` needs to be specified.

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

        Also, this map can be defined by specifying where the
        1-simplex goes; the vertices then go where they have to, to
        satisfy the condition `d_i \circ f = f \circ d_i`::

            sage: H = Hom(K, S1)
            sage: H({e01: sigma})
            Simplicial set morphism:
              From: 1-simplex
              To:   S^1
              Defn: [(0,), (1,), (0, 1)] --> [v_0, v_0, sigma_1]

        A constant map::

            sage: g = {e01: w.apply_degeneracies(0)}
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

        Another non-map::

            sage: h = {w: v0, v0: w, sigma: e01}
            sage: SimplicialSetMorphism(h, S1, K)
            Traceback (most recent call last):
            ...
            ValueError: at least one simplex in the defining dictionary is not in the domain

        An improperly partially defined map::

            sage: h = {w: v0}
            sage: SimplicialSetMorphism(h, S1, K)
            Traceback (most recent call last):
            ...
            ValueError: the image of at least one simplex in the domain is not defined

        A (good) partially defined map::

            sage: S5 = simplicial_sets.Sphere(5)
            sage: s = S5.n_cells(5)[0]
            sage: one = S5.Hom(S5)({s: s})
            sage: one
            Simplicial set endomorphism of S^5
              Defn: [v_0, sigma_5] --> [v_0, sigma_5]
            sage: one._dictionary
            {v_0: v_0, sigma_5: sigma_5}
        """
        if (not isinstance(domain, SimplicialSet_infinite)
            or not isinstance(codomain, SimplicialSet_infinite)):
            raise ValueError('the domain and codomain must be simplicial sets')
        if any(x.nondegenerate() not in
               domain.nondegenerate_simplices() for x in f.keys()):
            raise ValueError('at least one simplex in the defining '
                             'dictionary is not in the domain')
        # Remove degenerate simplices from the domain specification.
        d = {sigma:f[sigma] for sigma in f if sigma.is_nondegenerate()}
        # For each simplex in d.keys(), add its faces, and the faces
        # of its faces, etc., to d.
        for simplex in d.keys():
            faces = domain.faces(simplex)
            add = []
            if faces:
                for (i,sigma) in enumerate(faces):
                    nondegen = sigma.nondegenerate()
                    if nondegen not in d:
                        add.append((sigma,i,simplex))
            while add:
                (sigma,i,tau) = add.pop()
                # sigma is the ith face of tau.
                face_f = codomain.face(d[tau], i)
                degens = sigma.degeneracies()
                x = face_f
                for j in degens:
                    x = codomain.face(x, j)
                d[sigma.nondegenerate()] = x
                faces = domain.faces(sigma.nondegenerate())
                if faces:
                    for (i,rho) in enumerate(faces):
                        nondegen = rho.nondegenerate()
                        if nondegen not in d:
                            add.append((rho,i,sigma))
        # Now check that the proposed map commutes with the face
        # maps. (The degeneracy maps should work automatically.)
        for simplex in d:
            # Compare d[d_i (simplex)] to d_i d[simplex]. Since
            # d_i(simplex) may be degenerate, we have to be careful
            # when applying f to it. We can skip vertices and start
            # with 1-simplices.
            bad = False
            for i in range(simplex.dimension()+1):
                face_f = codomain.face(d[simplex], i)
                face = domain.face(simplex, i)
                if face is None:
                    f_face = None
                elif face.is_nondegenerate():
                    f_face = d[face]
                else:
                    nondegen = face.nondegenerate()
                    f_face = d[nondegen].apply_degeneracies(*face.degeneracies())
                if face_f != f_face:
                    bad = True
                    break
            if bad:
                raise ValueError('the dictionary does not define a map of simplicial sets')
        if any(x not in d.keys() for x in domain.nondegenerate_simplices()):
            raise ValueError('the image of at least one simplex in '
                             'the domain is not defined')
        self._dictionary = d
        Morphism.__init__(self, Hom(domain, codomain, SimplicialSets()))

    def __eq__(self, other):
        """
        Two morphisms are equal iff their domains are the same, their
        codomains are the same, and their defining dictionaries are
        the same.

        EXAMPLES::

            sage: S = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: T_c = T.constant_map() * T.base_point_map()
            sage: S_c = S.constant_map() * S.base_point_map()
            sage: T_c == S_c
            True
            sage: T.constant_map() == S.constant_map()
            False
            sage: K = simplicial_sets.Sphere(1)
            sage: K.constant_map() == S.constant_map()
            False

            sage: Point = simplicial_sets.Point()
            sage: f = Point._map_from_empty_set()
            sage: Empty = f.domain()
            sage: g = Empty.constant_map()
            sage: f == g
            True
        """
        return (self.domain() == other.domain() and self.codomain() == other.codomain()
                and self._dictionary == other._dictionary)

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: v,w = S0.n_cells(0)
            sage: H = Hom(S0, S0)
            sage: H({v:v, w:w}) != H({v:w, w:v})
            True
            sage: H({v:v, w:w}) != H({w:w, v:v})
            False
        """
        return not self == other

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
            return self._dictionary[x.nondegenerate()].apply_degeneracies(*x.degeneracies())
        except KeyError:
            raise ValueError('element is not a simplex in the domain')

    def _composition_(self, right, homset):
        """
        INPUT:

        - ``self``, ``right`` -- maps
        - ``homset`` -- a homset

        ASSUMPTION:

        The codomain of ``right`` is contained in the domain of
        ``self``.  This assumption should be verified by the
        ``Map.__mul__`` method in ``categories/map.pyx``, so we don't
        need to check it here.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: f = S1.Hom(S1).identity()
            sage: f * f # indirect doctest
            Simplicial set endomorphism of S^1
              Defn: [v_0, sigma_1] --> [v_0, sigma_1]
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: g = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: g
            Simplicial set morphism:
              From: S^1
              To:   Simplicial set with 6 non-degenerate simplices
              Defn: [v_0, sigma_1] --> [(v_0, v_0), (sigma_1, Simplex obtained by applying degeneracy s_0 to v_0)]
            sage: (g*f).image()
            Simplicial set with 2 non-degenerate simplices
            sage: f.image().homology()
            {0: 0, 1: Z}
        """
        if self.is_identity():
            return right
        if right.is_identity():
            return self
        d = {}
        for sigma in right._dictionary:
            d[sigma] = self(right(sigma))
        return homset(d)

    def image(self):
        """
        The image of this morphism as a subsimplicial set of the codomain.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: f = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: f
            Simplicial set morphism:
              From: S^1
              To:   Simplicial set with 6 non-degenerate simplices
              Defn: [v_0, sigma_1] --> [(v_0, v_0), (sigma_1, Simplex obtained by applying degeneracy s_0 to v_0)]
            sage: f.image()
            Simplicial set with 2 non-degenerate simplices
            sage: f.image().homology()
            {0: 0, 1: Z}
        """
        return self.codomain().subsimplicial_set(self._dictionary.values())

    def is_identity(self):
        """
        True if this morphism is an identity map.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(1)
            sage: v0 = K.n_cells(0)[0]
            sage: v1 = K.n_cells(0)[1]
            sage: e01 = K.n_cells(1)[0]
            sage: L = simplicial_sets.Simplex(2).n_skeleton(1)
            sage: w0 = L.n_cells(0)[0]
            sage: w1 = L.n_cells(0)[1]
            sage: w2 = L.n_cells(0)[2]
            sage: f01 = L.n_cells(1)[0]
            sage: f02 = L.n_cells(1)[1]
            sage: f12 = L.n_cells(1)[2]

            sage: d = {v0:w0, v1:w1, e01:f01}
            sage: f = K.Hom(L)(d)
            sage: f.is_identity()
            False
            sage: d = {w0:v0, w1:v1, w2:v1, f01:e01, f02:e01, f12: v1.apply_degeneracies(0,)}
            sage: g = L.Hom(K)(d)
            sage: (g*f).is_identity()
            True
            sage: (f*g).is_identity()
            False
            sage: (f*g).induced_homology_morphism().to_matrix(1)
            [0]

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP5.n_skeleton(2).inclusion_map().is_identity()
            False
            sage: RP5.n_skeleton(5).inclusion_map().is_identity()
            True
        """
        return (self.domain() == self.codomain()
                and all(a == b for a,b in self._dictionary.items()))

    def is_surjective(self):
        """
        Return ``True`` if this map is surjective.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP2.inclusion_map().is_surjective()
            False

            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.quotient_map().is_surjective()
            True

            sage: K = RP5_2.pullback(RP5_2.quotient_map(), RP5_2.base_point_map())
            sage: f = K.universal_property(RP2.inclusion_map(), RP2.constant_map())
            sage: f.is_surjective()
            True
        """
        return self.image() == self.codomain()

    def is_injective(self):
        """
        Return ``True`` if this map is surjective.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP2.inclusion_map().is_injective()
            True

            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.quotient_map().is_injective()
            False

            sage: K = RP5_2.pullback(RP5_2.quotient_map(), RP5_2.base_point_map())
            sage: f = K.universal_property(RP2.inclusion_map(), RP2.constant_map())
            sage: f.is_injective()
            True
        """
        domain = self.domain()
        for n in range(domain.dimension()+1):
            input = domain.n_cells(n)
            output = set([self(sigma) for sigma in input if self(sigma).is_nondegenerate()])
            if len(input) > len(output):
                return False
        return True

    def is_bijective(self):
        """
        Return ``True`` if this map is surjective.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP2.inclusion_map().is_bijective()
            False

            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.quotient_map().is_bijective()
            False

            sage: K = RP5_2.pullback(RP5_2.quotient_map(), RP5_2.base_point_map())
            sage: f = K.universal_property(RP2.inclusion_map(), RP2.constant_map())
            sage: f.is_bijective()
            True
        """
        return self.is_injective() and self.is_surjective()

    def is_pointed(self):
        """
        Return ``True`` if this is a pointed map.

        That is, return ``True`` if the domain and codomain are
        pointed and this morphism preserves the base point.

        EXAMPLES::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: f = Hom(S0,S0).identity()
            sage: f.is_pointed()
            True
            sage: v = S0.n_cells(0)[0]
            sage: w = S0.n_cells(0)[1]
            sage: g = Hom(S0,S0)({v:v, w:v})
            sage: g.is_pointed()
            True
            sage: t = Hom(S0,S0)({v:w, w:v})
            sage: t.is_pointed()
            False
        """
        return (self.domain().is_pointed() and self.codomain().is_pointed()
                and self(self.domain().base_point()) == self.codomain().base_point())

    def pushout(self, *others):
        """
        The pushout of this morphism along with ``others``.

        INPUT:

        - ``others`` -- morphisms of simplicial sets, the domains of
          which must all equal that of ``self``.

        This returns the pushout as a simplicial set. See
        :class:`sage.homology.simplicial_set.PushoutOfSimplicialComplexes`
        for more documentation and examples.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: init_T = T._map_from_empty_set()
            sage: init_K = K._map_from_empty_set()
            sage: D = init_T.pushout(init_K) # the disjoint union as a pushout
            sage: D
            Simplicial set with 12 non-degenerate simplices
        """
        domain = self.domain()
        if any(domain != f.domain() for f in others):
            raise ValueError('the domains of the maps must be equal')
        return PushoutOfSimplicialSets((self,) + others)

    def pullback(self, *others):
        """
        The pullback of this morphism along with ``others``.

        INPUT:

        - ``others`` -- morphisms of simplicial sets, the codomains of
          which must all equal that of ``self``.

        This returns the pullback as a simplicial set. See
        :class:`sage.homology.simplicial_set.PullbackOfSimplicialComplexes`
        for more documentation and examples.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: term_T = T.constant_map()
            sage: term_K = K.constant_map()
            sage: P = term_T.pullback(term_K) # the product as a pullback
            sage: P
            Simplicial set with 150 non-degenerate simplices
        """
        codomain = self.codomain()
        if any(codomain != f.codomain() for f in others):
            raise ValueError('the codomains of the maps must be equal')
        return PullbackOfSimplicialSets((self,) + others)

    def equalizer(self, other):
        r"""
        The equalizer of this map with ``other``.

        INPUT:

        - ``other`` -- a morphism with the same domain and codomain as this map

        If the two maps are `f, g: X \to Y`, then the equalizer `P` is
        constructed as the pullback ::

            P ----> X
            |       |
            V       V
            X --> X x Y

        where the two maps `X \to X \times Y` are `(1,f)` and `(1,g)`.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: x = AbstractSimplex(0, name='x')
            sage: evw = AbstractSimplex(1, name='vw')
            sage: evx = AbstractSimplex(1, name='vx')
            sage: ewx = AbstractSimplex(1, name='wx')
            sage: X = SimplicialSet({evw: (w, v), evx: (x, v)})
            sage: Y = SimplicialSet({evw: (w, v), evx: (x, v), ewx: (x, w)})

        Here `X` is a wedge of two 1-simplices (a horn, that is), and
        `Y` is the boundary of a 2-simplex. The map `f` includes the
        two 1-simplices into `Y`, while the map `g` maps both
        1-simplices to the same edge in `Y`. ::

            sage: f = Hom(X, Y)({v:v, w:w, x:x, evw:evw, evx:evx})
            sage: g = Hom(X, Y)({v:v, w:x, x:x, evw:evx, evx:evx})
            sage: P = f.equalizer(g)
            sage: P
            Simplicial set with 3 non-degenerate simplices
        """
        domain = self.domain()
        codomain = self.codomain()
        if domain != other.domain() or codomain != other.codomain():
            raise ValueError('the maps must have the same domain and the same codomain')
        prod = domain.product(codomain)
        one = domain.Hom(domain).identity()
        f = prod.universal_property(one, self)
        g = prod.universal_property(one, other)
        return PullbackOfSimplicialSets([f, g])

    def coequalizer(self, other):
        r"""
        The coequalizer of this map with ``other``.

        INPUT:

        - ``other`` -- a morphism with the same domain and codomain as this map

        If the two maps are `f, g: X \to Y`, then the coequalizer `P` is
        constructed as the pushout ::

            X v Y --> Y
              |       |
              V       V
              Y ----> P

        where the upper left corner is the coproduct of `X` and `Y`
        (the wedge if they are pointed, the disjoint union otherwise),
        and the two maps `X \amalg Y \to Y` are `f \amalg 1` and `g
        \amalg 1`.

        EXAMPLES::

            sage: L = simplicial_sets.Simplex(2)
            sage: pt = L.n_cells(0)[0]
            sage: e = L.n_cells(1)[0]
            sage: K = L.subsimplicial_set([e])
            sage: f = K.inclusion_map()
            sage: v,w = K.n_cells(0)
            sage: g = Hom(K,L)({v:pt, w:pt, e:pt.apply_degeneracies(0)})
            sage: P = f.coequalizer(g)
            sage: P
            Simplicial set with 5 non-degenerate simplices
        """
        domain = self.domain()
        codomain = self.codomain()
        if domain != other.domain() or codomain != other.codomain():
            raise ValueError('the maps must have the same domain and the same codomain')
        coprod = domain.coproduct(codomain)
        one = codomain.Hom(codomain).identity()
        f = coprod.universal_property(self, one)
        g = coprod.universal_property(other, one)
        return PushoutOfSimplicialSets([f, g])

    def mapping_cone(self):
        r"""
        The mapping cone defined by this map.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: v_0, sigma_1 = S1.nondegenerate_simplices()
            sage: K = simplicial_sets.Simplex(2).n_skeleton(1)

        The mapping cone will be a little smaller if we use only
        pointed simplicial sets. `S^1` is already pointed, but not
        `K`. ::

            sage: L = K.set_base_point(K.n_cells(0)[0])
            sage: u,v,w = L.n_cells(0)
            sage: e,f,g = L.n_cells(1)
            sage: h = L.Hom(S1)({u:v_0, v:v_0, w:v_0, e:sigma_1, f:v_0.apply_degeneracies(0), g:sigma_1})
            sage: h
            Simplicial set morphism:
              From: Simplicial set with 6 non-degenerate simplices
              To:   S^1
              Defn: [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2)] --> [v_0, v_0, v_0, sigma_1, Simplex obtained by applying degeneracy s_0 to v_0, sigma_1]
            sage: h.induced_homology_morphism().to_matrix()
            [1|0]
            [-+-]
            [0|2]
            sage: X = h.mapping_cone()
            sage: X.homology() == simplicial_sets.RealProjectiveSpace(2).homology()
            True
        """
        dom = self.domain()
        cone = dom.cone()
        i = cone.map_from_X()
        return self.pushout(i)

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

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')
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
        return InducedHomologyMorphism(self, base_ring, cohomology)

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
