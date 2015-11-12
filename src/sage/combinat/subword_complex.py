r"""
Subword complex

Fix a Coxeter system `(W,S)`. The subword complex `\mathcal{SC}(Q,w)`
associated to a word `Q \in S^*` and an element `w \in W` the
simplicial whose ground set is the set of positions in `Q` and whose
facets are complements of sets of positions defining a reduced
expression for `w`.

A subword complex is a shellable sphere if and only if the Demazure
product of `Q` equals `w`, otherwise it is a shellable ball.

AUTHORS:

- Christian Stump: initial version
- Vincent Pilaud: greedy flip algorithm, minor improvements, documentation

REFERENCES:

.. [KnuMil] Knutson and Miller. Subword complexes in Coxeter groups. Adv. Math., 184(1):161-176, 2004.
.. [PilStu] Pilaud and Stump. Brick polytopes of spherical subword complexes and generalized associahedra. Adv. Math. 276:1-61, 2015.
"""
#*****************************************************************************
#       Copyright (C) 2012      Christian Stump <christian.stump@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.homology.simplicial_complex import SimplicialComplex, Simplex
from sage.geometry.cone import Cone
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.matrix.all import Matrix
from copy import copy
from sage.combinat.subword_complex_c import _flip_c, _construct_facets_c
from sage.geometry.polyhedron.constructor import Polyhedron


class SubwordComplex(SimplicialComplex, Parent):
    r"""
    Fix a Coxeter system `(W,S)`. The subword complex
    `\mathcal{SC}(Q,w)` associated to a word `Q \in S^*` and an
    element `w \in W` the simplicial whose ground set is the set of
    positions in `Q` and whose facets are complements of sets of
    positions defining a reduced expression for `w`.

    A subword complex is a shellable sphere if and only if the
    Demazure product of `Q` equals `w`, otherwise it is a shellable
    ball.

    EXAMPLES:

    As an example, dual associahedra are subword complexes in type
    `A_{n-1}` given by the word `[1, \dots, n, 1, \dots, n, 1, \dots,
    n-1, \dots, 1, 2, 1]` and the permutation `w_0`.

    ::

        sage: W = CoxeterGroup(['A',2], index_set=[1,2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w); SC
        Subword complex of type ['A', 2] for Q = [1, 2, 1, 2, 1] and pi = [1, 2, 1]
        sage: SC.facets()
        [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

    REFERENCES: [KnuMil]_, [PilStu]_
    """

    # standard functions

    def __init__(self, Q, w, inversion_set_indices=None, algorithm="inductive"):
        r"""
        Initialize the subword complex `\mathcal{SC}(Q,w)`.

        INPUT:

        - ``Q`` -- word on the simple generators of the Coxeter group.
        - ``w`` -- element of the Coxeter group.
        - ``invertion_set_indices`` -- (default: None)
        - ``algorithm`` -- (default: inductive) choice of the
          algorithm to generate the subword complex. Options are
          ``inductive`` or ``greedy``. The second option is
          recommended when `|Q|` is closed to `\ell(w) +
          \mathrm{rank}(W)`.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], index_set=[1,2,3])
            sage: w = W.from_reduced_word([1,2,3,1,2,1])
            sage: SC = SubwordComplex([1,2,3,1,2,3,1,2,1], w); SC
            Subword complex of type ['A', 3] for Q = [1, 2, 3, 1, 2, 3, 1, 2, 1] and pi = [1, 2, 3, 1, 2, 1]
            sage: len(SC)
            14
        """
        W = w.parent()
        I = W.index_set()
        if not all(i in I for i in Q):
            raise ValueError("All elements in Q = %s must be contained in the index set %s" % (Q, W.index_set()))
        self._Q = Q
        self._pi = w
        if algorithm == "inductive":
            Fs = _construct_facets_c(Q, w)
        elif algorithm == "greedy":
            Fs, Rs = _greedy_flip_algorithm(Q, w)
        else:
            raise ValueError("The optional argument algorithm can be either inductive or greedy")
        if Fs == []:
            raise ValueError("The word %s does not contain a reduced expression for %s" % (Q, w.reduced_word()))
        SimplicialComplex.__init__(self, maximal_faces=Fs,
                                   maximality_check=False)
        self.__custom_name = 'Subword complex'
        self._W = W
        if hasattr(W, '_type'):
            self._cartan_type = W._type
        else:
            self._cartan_type = None
        self._facets_dict = None
        if algorithm == "greedy":
            _facets_dict = {}
            for i in range(len(Fs)):
                X = self(Fs[i], facet_test=False)
                X._extended_root_conf_indices = Rs[i]
                _facets_dict[tuple(sorted(Fs[i]))] = X
            self._facets_dict = _facets_dict
        else:
            self._facets_dict = {}

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SubwordComplex([1,2,1,2,1], w)
            Subword complex of type ['A', 2] for Q = [1, 2, 1, 2, 1] and pi = [1, 2, 1]
        """
        if self._cartan_type is None:
            return 'Subword complex of unknown type for Q = %s and pi = %s' % (self._Q, self._pi.reduced_word())
        else:
            return 'Subword complex of type %s for Q = %s and pi = %s' % (self.cartan_type(), self._Q, self._pi.reduced_word())

    def __eq__(self, other):
        r"""
        Compare the subword complexes ``self`` and ``other``.

        INPUT:

        - ``other`` -- another subword complex.

        EXAMPLE::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC1 = SubwordComplex([1,2,1,2,1], w); SC2 = SubwordComplex([1,2,1,2,1], w); SC1 == SC2
            True
        """
        return self.word() == other.word() and self.pi() == other.pi()

    def __call__(self, F, facet_test=True):
        r"""
        Create a facet of ``self``.

        INPUT:

        - ``F`` -- an iterable of positions.
        - ``facet_test`` -- boolean (default: True) tells whether or not the facet ``F``
          should be tested before creation.

        OUTPUT: the facet of ``self`` at positions given by ``F``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: F = SC([1,2]); F
            (1, 2)

        """
        F = tuple(F)
        if self._facets_dict is not None and self._facets_dict != dict() and F in self._facets_dict:
            return self._facets_dict[F]
        else:
            return SubwordComplexFacet(self, F, facet_test=facet_test)

    def __contains__(self, F):
        r"""
        Tests if ``self`` contains a given iterable ``F``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w  = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
            sage: [0,1] in SC
            True
            sage: [0,2] in SC
            False
            sage: [0,1,5] in SC
            False
            sage: [0] in SC
            False
            sage: ['a','b'] in SC
            False
        """
        W = self.group()
        Q = self.word()
        if not all(i in range(len(Q)) for i in F):
            return False
        return W.from_reduced_word(Q[i] for i in range(len(Q)) if i not in F) == self.pi()

    def list(self):
        r"""
        Return the list of facets of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: list(SC)
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
        """
        return [F for F in self]

    # getting the stored properties

    def group(self):
        r"""
        Return the group associated to ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.group()
            Irreducible finite Coxeter group of rank 2 and type A2
        """
        return self._W

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.cartan_type()
            ['A', 2]
        """
        from sage.combinat.root_system.cartan_type import CartanType
        if self._cartan_type is None:
            raise ValueError("No Cartan type defined for %s" % (self._W))
        elif len(self._cartan_type) == 1:
            return CartanType([self._cartan_type[0]['series'],
                               self._cartan_type[0]['rank']])
        else:
            return CartanType(self._cartan_type)

    def word(self):
        r"""
        Return the word in the simple generators associated to ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.word()
            [1, 2, 1, 2, 1]
        """
        return copy(self._Q)

    def pi(self):
        r"""
        Return the element in the Coxeter group associated to ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.pi().reduced_word()
            [1, 2, 1]
        """
        return self._pi

    def facets(self):
        r"""
        Return all facets of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
        """
        if self._facets_dict:
            return [self._facets_dict[tuple(F)] for F in self._facets]
        else:
            return [self(F, facet_test=False) for F in self._facets]

    def __iter__(self):
        r"""
        Return an iterator on the facets of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: for X in SC:
            ...       print X
            (0, 1)
            (0, 4)
            (1, 2)
            (2, 3)
            (3, 4)
        """
        return iter(self.facets())

    def greedy_facet(self, side="positive"):
        r"""
        Return the negative (or positive) greedy facet of ``self``.

        This is the lexicographically last (or first) facet of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.greedy_facet(side="positive")
            (0, 1)
            sage: SC.greedy_facet(side="negative")
            (3, 4)
        """
        return SubwordComplexFacet(self, _greedy_facet(self.word(),
                                                       self.pi(), side=side))

    # topological properties

    def is_sphere(self):
        r"""
        Return ``True`` if the subword complex ``self`` is a sphere.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], index_set=[1,2,3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_sphere()
            True

            sage: SC = SubwordComplex([3,2,1,3,2,3], w)
            sage: SC.is_sphere()
            False
        """
        W = self._pi.parent()
        w = W.demazure_product(self._Q)
        return w == self._pi

    def is_ball(self):
        r"""
        Return ``True`` if the subword complex ``self`` is a ball.

        This is the case if and only if it is not a sphere.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], index_set=[1,2,3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_ball()
            False

            sage: SC = SubwordComplex([3,2,1,3,2,3], w)
            sage: SC.is_ball()
            True
        """
        return not self.is_sphere()

    def is_pure(self):
        r"""
        Return ``True`` since all subword complexes are pure.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], index_set=[1,2,3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_pure()
            True
        """
        return True

    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.dimension()
            1
        """
        return self._facets[0].dimension()

    # root and weight

    @cached_method
    def is_root_independent(self):
        r"""
        Return ``True`` if ``self`` is root-independent, that is if the root configuration
        of any (or equivalently all) facets is linearly independent.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.is_root_independent()
            True

            sage: SC = SubwordComplex([1,2,1,2,1,2], W.w0)
            sage: SC.is_root_independent()
            False
        """
        from sage.matrix.all import matrix
        M = matrix(self.greedy_facet(side="negative").root_configuration())
        return M.rank() == max(M.ncols(), M.nrows())

    @cached_method
    def is_double_root_free(self):
        r"""
        Return ``True`` if ``self`` is double-root-free, that is if the root configurations
        of all facets do not contain a root twice.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.is_double_root_free()
            True

            sage: SC = SubwordComplex([1,1,2,2,1,1], w)
            sage: SC.is_double_root_free()
            True

            sage: SC = SubwordComplex([1,2,1,2,1,2], w)
            sage: SC.is_double_root_free()
            False
        """
        if not self.is_root_independent():
            size = self.dimension() + 1
            for F in self:
                conf = F._root_configuration_indices()
                if len(set(conf)) < size:
                    return False
        return True

    def kappa_preimages(self):
        """
        Return a dictionary containing facets of ``self`` as keys,
        and list of elements of ``self.group()`` as values.

        .. SEEALSO::

            :meth:`kappa_preimage <sage.combinat.subword_complex.SubwordComplexFacet.kappa_preimage>`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: kappa_pres = SC.kappa_preimages()
            sage: for F in SC: print F, [ w.reduced_word() for w in kappa_pres[F] ]
            (0, 1) [[]]
            (0, 4) [[2, 1], [2]]
            (1, 2) [[1]]
            (2, 3) [[1, 2]]
            (3, 4) [[1, 2, 1]]
        """
        return {F: F.kappa_preimage() for F in self}

    def compatibility_fan(self, F=None):
        r"""
        Return the root fan of ``self`` with respect to the initial facet ``F``.

        INPUT:

        - ``F`` -- (default: ``greedy_facet(self)``) an iterable of positions.

        .. SEEALSO::

            :func:`root_fan <sage.combinat.subword_complex.SubwordComplexFacet.root_fan>`

        EXAMPLES::
        """
        if F is None:
            F = self.greedy_facet()
        return self(F).compatibility_fan()

    def brick_fan(self):
        r"""
        Return the brick fan of ``self``.

        It is the normal fan of the brick polytope of ``self``. It is
        formed by the cones generated by the weight configurations of
        the facets of ``self``.

        .. SEEALSO::

            :func:`weight_cone <sage.combinat.subword_complex.SubwordComplexFacet.weight_cone>`

        EXAMPLES::
        """
        from sage.geometry.fan import Fan
        return Fan([F.weight_cone() for F in self])

    # brick polytope

    def brick_vectors(self, coefficients=None):
        r"""
        Return the list of all brick vectors of facets of ``self``.

        .. SEEALSO::

            :func:`brick_vector <sage.combinat.subword_complex.SubwordComplexFacet.brick_vector>`

        EXAMPLES::
        """
        return [F.brick_vector(coefficients=coefficients) for F in self]

    def minkowski_summand(self, i):
        r"""
        Return the ``i``th Minkowski summand of ``self``.

        EXAMPLES::
        """
        return Polyhedron([F.extended_weight_configuration()[i] for F in self])

    def brick_polytope(self, coefficients=None):
        r"""
        Return the brick polytope of ``self``.

        This polytope is the convex hull of the brick vectors of ``self``.
        The coefficients

        .. SEEALSO::

            :meth:`brick_vectors`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)

            sage: X = SC.brick_polytope(); X
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices

            sage: Y = SC.brick_polytope(coefficients=[1,2]); Y
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices

            sage: X == Y
            False
        """
        BV = self.brick_vectors(coefficients=coefficients)
        if not self.group().is_crystallographic():
            from sage.rings.all import CC, QQ
            print "Caution: the polytope is build with rational vertices."
            BV = [[QQ(CC(v).real()) for v in V] for V in BV]
        return Polyhedron(BV)

    def word_orbit_partition(self):
        r"""
        Return the partition of the positions of `Q` into orbits of letter rotation.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: SC.word_orbit_partition()
            [[0, 1, 2, 3, 4]]

            sage: W = CoxeterGroup(['A',3],index_set=[1,2,3])
            sage: w = W.from_reduced_word([1,2,3,1,2,1])
            sage: SC = SubwordComplex([1,2,3,1,2,3,1,2,1],w)
            sage: SC.word_orbit_partition()
            [[0, 2, 3, 5, 6, 8], [1, 4, 7]]
        """
        Q = self.word()
        positions = range(len(Q))
        W = self.group()
        N = W.long_element().length()
        w = W.w0
        I = W.index_set()
        I_decomp = set([frozenset([i, I[w(W.index_set()[i] + 1) - 1 - N]])
                        for i in I])
        return [[i for i in positions if Q[i] in I_part]
                for I_part in I_decomp]

    def barycenter(self):
        """
        Return the barycenter of the brick polytope of ``self``.

        .. SEEALSO::

            :meth:`brick_polytope`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.barycenter()
            (2/3, 4/3)
        """
        facets = self.facets()
        if not self.is_root_independent():
            facets = [F for F in facets if F.is_vertex()]
        return sum(F.brick_vector() for F in facets) / len(facets)

    # cambrian constructions

    def cover_relations(self, label=False):
        N = self.group().long_element().length()
        F = self.greedy_facet(side="positive")
        Fs = set([F])
        seen = set([F])
        covers = []
        while Fs != set():
            F = Fs.pop()
            seen.add(F)
            conf = F._extended_root_configuration_indices()
            for i in F:
                if conf[i] < N:
                    G = F.flip(i)
                    if label:
                        covers.append((F, G, i))
                    else:
                        covers.append((F, G))
                    if G not in seen:
                        Fs.add(G)
        return covers

    def increasing_flip_graph(self, label=True):
        from sage.graphs.digraph import DiGraph
        return DiGraph(self.cover_relations(label=label))

    def interval(self, I, J):
        G = self.increasing_flip_graph()
        paths = G.all_paths(I, J)
        return set(K for path in paths for K in path)

    def increasing_flip_poset(self):
        from sage.combinat.posets.posets import Poset
        cov = self.cover_relations()
        if not self.is_root_independent():
            Fs = [F for F in self if F.is_vertex()]
            cov = [(a, b) for a, b in cov if a in Fs and b in Fs]
        return Poset(((), cov), facade=True)


class SubwordComplexFacet(Simplex, Element):
    r"""
    A facet of a subword complex.

    Facets of the subword complex `\mathcal{SC}(Q,w)` are complements of sets of positions
    in `Q` defining a reduced expression for `w`.

    EXAMPLES::

        sage: W = CoxeterGroup(['A',2], index_set=[1,2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w)
        sage: F = SC[0]
        sage: F
        (0, 1)
        sage: type(F)
        <class 'sage.combinat.subword_complex.SubwordComplexFacet'>

    """

    # standard functions

    def __init__(self, parent, positions, facet_test=True):
        r"""
        Initializes a facet of the subword complex ``parent``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: F = SC([1,2]); F
            (1, 2)
        """
        if facet_test and positions not in parent:
            raise ValueError("The given iterable %s is not a facet of the subword complex %s" % (positions, parent))
        Element.__init__(self, parent)
        Simplex.__init__(self, sorted(positions))
        self._extended_root_conf_indices = None
        self._extended_weight_conf = None

    def __eq__(self, other):
        r"""
        Compare the subword complexes facets ``self`` and ``other``.

        INPUT:

        - ``other`` -- another subword complex facet.

        EXAMPLE::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: F1 = SC([0,1]); F2 = SC([0,1]); F1 == F2
            True
        """
        return self.parent() == other.parent() and self.tuple() == other.tuple()

    # roots

    def _extended_root_configuration_indices(self):
        r"""
        Return the indices of the roots in ``self.group().roots()`` of
        the extended root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The extended
        root configuration of a facet `I` of `\mathcal{SC}(Q,w)` is
        the sequence `\mathsf{r}(I, 1), \dots, \mathsf{r}(I, m)` of
        roots defined by `\mathsf{r}(I, k) = \Pi Q_{[k-1]
        \smallsetminus I} (\alpha_{q_k})`, where `\Pi Q_{[k-1]
        \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        .. SEEALSO::

            :meth:`extended_root_configuration`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F._extended_root_configuration_indices()
            [0, 2, 3, 2, 1]
        """
        if self._extended_root_conf_indices is None:
            self._extended_root_conf_indices = _extended_root_configuration_indices(self.parent().group(), self.parent().word(), self)
        return self._extended_root_conf_indices

    def _root_configuration_indices(self):
        r"""
        Return the indices of the roots in ``self.group().roots()`` of the root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The root configuration
        of a facet `I = [i_1, \dots, i_n]` of `\mathcal{SC}(Q,w)` is the sequence `\mathsf{r}(I, i_1), \dots, \mathsf{r}(I, i_n)`
        of roots defined by `\mathsf{r}(I, k) = \Pi Q_{[k-1] \smallsetminus I} (\alpha_{q_k})`, where
        `\Pi Q_{[k-1] \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        .. SEEALSO::

            :meth:`root_configuration`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F._root_configuration_indices()
            [2, 3]
        """
        indices = self._extended_root_configuration_indices()
        return [indices[i] for i in self]

    def extended_root_configuration(self):
        r"""
        Return the extended root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The extended root configuration
        of a facet `I` of `\mathcal{SC}(Q,w)` is the sequence `\mathsf{r}(I, 1), \dots, \mathsf{r}(I, m)`
        of roots defined by `\mathsf{r}(I, k) = \Pi Q_{[k-1] \smallsetminus I} (\alpha_{q_k})`, where
        `\Pi Q_{[k-1] \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        The extended root configuration is used to perform flips efficiently.

        .. SEEALSO::

            :meth:`flip`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.extended_root_configuration()
            [(1, 0), (1, 1), (-1, 0), (1, 1), (0, 1)]
        """
        Phi = self.parent().group().roots()
        return [Phi[i] for i in self._extended_root_configuration_indices()]

    def root_configuration(self):
        r"""
        Return the root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The root configuration
        of a facet `I = [i_1, \dots, i_n]` of `\mathcal{SC}(Q,w)` is the sequence `\mathsf{r}(I, i_1), \dots, \mathsf{r}(I, i_n)`
        of roots defined by `\mathsf{r}(I, k) = \Pi Q_{[k-1] \smallsetminus I} (\alpha_{q_k})`, where
        `\Pi Q_{[k-1] \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.root_configuration()
            [(1, 1), (-1, 0)]
        """
        Phi = self.parent().group().roots()
        return [Phi[i] for i in self._root_configuration_indices()]

    def kappa_preimage(self):
        r"""
        Return the fiber of ``self`` under the `\kappa` map.

        The `\kappa` map sends an element `w \in W` to the unique facet of `I \in \mathcal{SC}(Q,w)`
        such that the cone generated by `w(\Phi^+)` is contained in the cone
        generated by the root configuration of `I`.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)

            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.kappa_preimage()
            [(1,4)(2,3)(5,6)]

            sage: F = SC([0,4]); F
            (0, 4)
            sage: F.kappa_preimage()
            [(1,3)(2,5)(4,6), (1,2,6)(3,4,5)]
        """
        W = self.parent().group()
        N = W.long_element().length()
        root_conf = self._root_configuration_indices()
        return [w for w in W if all(w(i + 1) <= N for i in root_conf)]

    def is_vertex(self):
        r"""
        Return ``True`` if ``self`` is a vertex of the brick polytope
        of ``self.parent``.

        A facet is a vertex of the brick polytope if its root cone is
        pointed. Note that this property is always satisfied for
        root-independent subword complexes.

        .. SEEALSO::

            :meth:`root_cone`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',1], index_set=[1])
            sage: w = W.from_reduced_word([1])
            sage: SC = SubwordComplex([1,1,1],w)
            sage: F = SC([0,1]); F.is_vertex()
            True
            sage: F = SC([0,2]); F.is_vertex()
            False

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1,2,1],w)
            sage: F = SC([0,1,2,3]); F.is_vertex()
            True
            sage: F = SC([0,1,2,6]); F.is_vertex()
            False
        """
        S = self.parent()
        if S.is_root_independent():
            return True
        return self.root_cone().is_strictly_convex()

    @cached_method
    def root_cone(self):
        r"""
        Return the polyhedral cone generated by the root configuration
        of ``self``.

        .. SEEALSO::

            :meth:`root_configuration`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',1], index_set=[1])
            sage: w = W.from_reduced_word([1])
            sage: SC = SubwordComplex([1,1,1],w)
            sage: F = SC([0,2]); F.root_cone()
            1-d cone in 1-d lattice N
        """
        return Cone(self.root_configuration())

    def compatibility_vectors(self):
        r"""
        Return the rays of the compatibility fan of the subword
        complex ``self.parent()`` with respect to the initial facet
        ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], index_set=[1,2,3])
            sage: w = W.from_reduced_word([1,2,3,1,2,1])
            sage: SC = SubwordComplex([1,2,3,1,2,3,1,2,1],w)
            sage: F = SC([0,1,2]); F.compatibility_vectors()
            [(-1, 0, 0),
             (0, -1, 0),
             (0, 0, -1),
             (1, 0, 0),
             (1, 1, 0),
             (1, 1, 1),
             (0, 1, 0),
             (0, 1, 1),
             (0, 0, 1)]
            sage: F = SC([0,2,6]); F.compatibility_vectors()
            [(-1, 0, 0),
             (0, 0, 1),
             (0, -1, 0),
             (1, 0, 1),
             (1, 0, 0),
             (1, 1, 0),
             (0, 0, -1),
             (0, 1, 0),
             (0, 1, 1)]
        """
        if not self.parent().is_root_independent():
            raise ValueError("compatibility fan only defined for "
                             "root-independent subword complexes")
        ext_root_conf = self.extended_root_configuration()
        root_conf = self.root_configuration()
        M = Matrix(root_conf)
        M_inverse = M.inverse()
        rays = []
        for j in range(len(ext_root_conf)):
            root = ext_root_conf[j]
            if j in self:
                rays.append(-root * M_inverse)
            else:
                new_root = root * M_inverse
                rays.append(vector([(1 if i < j else -1) * new_root[self.tuple().index(i)] for i in self]))
        return rays

    def compatibility_fan(self):
        r"""
        Return the compatibility fan of the subword complex ``self.parent()`` with respect to the initial facet ``self``.

        The compatibility fan of `\mathcal{SC}(Q,w)` with respect to an initial facet `I`
        is the fan whose rays are the compatibility vectors of the positions in `Q`
        with respect to `I` and whose maximal cones are the cones generated by the
        compatibility vectors corresponding to the facets of `\mathcal{SC}(Q,w)`.

        .. SEEALSO::

            :meth:`compatibility_vectors`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], index_set=[1,2,3])
            sage: w = W.from_reduced_word([1,2,3,1,2,1])
            sage: SC = SubwordComplex([1,2,3,1,2,3,1,2,1],w)
            sage: F = SC([0,1,2])
            sage: CF = F.compatibility_fan(); CF
            Rational polyhedral fan in 3-d lattice N
            sage: CF.rays()
            N(-1,  0,  0),
            N( 0, -1,  0),
            N( 0,  0, -1),
            N( 1,  0,  0),
            N( 1,  1,  0),
            N( 1,  1,  1),
            N( 0,  1,  0),
            N( 0,  1,  1),
            N( 0,  0,  1)
            in 3-d lattice N

            sage: F = SC([0,2,6])
            sage: CF = F.compatibility_fan(); CF
            Rational polyhedral fan in 3-d lattice N
            sage: CF.rays()
            N(-1,  0,  0),
            N( 0,  0,  1),
            N( 0, -1,  0),
            N( 1,  0,  1),
            N( 1,  0,  0),
            N( 1,  1,  0),
            N( 0,  0, -1),
            N( 0,  1,  0),
            N( 0,  1,  1)
            in 3-d lattice N
        """
        from sage.geometry.fan import Fan
        return Fan(self.parent().facets(), self.compatibility_vectors())

    def root_polytope(self, orbits):
        """
        """
        rays = self.compatibility_vectors()
        rhocheck = sum(rays[i] for i in range(len(rays)) if i not in self)
        inequalities = []
        for x in range(len(rays)):
            for orbit in orbits:
                if x in orbit:
                    i = min([j for j in self if j in orbit])
                    c = rhocheck[self.tuple().index(i)]
                    inequalities.append([c] + list(rays[x]))
        return Polyhedron(ieqs=inequalities)

    def upper_root_configuration(self):
        r"""
        Return the positive roots of the root configuration of ``self``.

        EXAMPLE::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.root_configuration()
            [(1, 1), (-1, 0)]
            sage: F.upper_root_configuration()
            [(1, 0)]
        """
        conf = self._root_configuration_indices()
        W = self.parent().group()
        Phi = W.roots()
        N = len(Phi) / 2
        return [Phi[i - N] for i in conf if i >= N]

    def product_of_upper_roots(self):
        r"""
        Return the product of the upper roots of ``self``.

        FIX: Connect to non-crossing partitions.

        .. SEEALSO::

            :meth:`upper_root_configuration`

        EXAMPLES::
        """
        W = self.parent().group()
        return W.prod(W.root_to_reflection(beta)
                      for beta in self.upper_root_configuration())

    # weights

    def extended_weight_configuration(self, coefficients=None):
        r"""
        Return the extended weight configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The extended
        weight configuration of a facet `I` of `\mathcal{SC}(Q,w)` is
        the sequence `\mathsf{w}(I, 1), \dots, \mathsf{w}(I, m)` of
        weights defined by `\mathsf{w}(I, k) = \Pi Q_{[k-1]
        \smallsetminus I} (\omega_{q_k})`, where `\Pi Q_{[k-1]
        \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        The extended weight configuration is used to compute the brick vector.

        .. SEEALSO::

            :meth:`brick_vector`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.extended_weight_configuration()
            [(2/3, 1/3), (1/3, 2/3), (-1/3, 1/3), (1/3, 2/3), (-1/3, 1/3)]
        """
        if coefficients is not None or self._extended_weight_conf is None:
            W = self.parent().group()
            Lambda = W.fundamental_weights()
            if coefficients is not None:
                Lambda = [coefficients[i] * Lambda[i]
                          for i in range(len(Lambda))]
            Q = self.parent().word()

            Phi = W.roots()

            V_weights = []

            pi = W.one()
            for i in range(len(Q)):
                wi = Q[i]
                fund_weight = Lambda[W.index_set()[wi]]
                V_weights.append(sum(fund_weight[j] * Phi[(~pi)(j + 1) - 1]
                                     for j in range(len(fund_weight))))
                if i not in self:
                    pi = pi.apply_simple_reflection(wi)
            if coefficients is None:
                self._extended_weight_conf = V_weights
            return V_weights
        else:
            return self._extended_weight_conf

    def weight_configuration(self):
        r"""
        Return the weight configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The weight
        configuration of a facet `I = [i_1, \dots, i_n]` of
        `\mathcal{SC}(Q,w)` is the sequence `\mathsf{w}(I, i_1),
        \dots, \mathsf{w}(I, i_n)` of weights defined by
        `\mathsf{w}(I, k) = \Pi Q_{[k-1] \smallsetminus I}
        (\omega_{q_k})`, where `\Pi Q_{[k-1] \smallsetminus I}` is the
        product of the simple reflections `q_i` for `i \in [k-1]
        \smallsetminus I` in this order.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.weight_configuration()
            [(1/3, 2/3), (-1/3, 1/3)]
        """
        extended_configuration = self.extended_weight_configuration()
        return [extended_configuration[i] for i in self]

    @cached_method
    def weight_cone(self):
        r"""
        Return the polyhedral cone generated by the weight configuration of ``self``.

        .. SEEALSO::

            :meth:`weight_configuration`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: WC = F.weight_cone(); WC
            2-d cone in 2-d lattice N
            sage: WC.rays()
            N( 1, 2),
            N(-1, 1)
            in 2-d lattice N
        """
        return Cone(self.weight_configuration())

    def brick_vector(self, coefficients=None):
        r"""
        Return the brick vector of ``self``, that is, the sum of the weight vectors
        in the extended weight configuration.

        .. SEEALSO::

            :meth:`extended_weight_configuration`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.extended_weight_configuration()
            [(2/3, 1/3), (1/3, 2/3), (-1/3, 1/3), (1/3, 2/3), (-1/3, 1/3)]
            sage: F.brick_vector()
            (2/3, 7/3)
        """
        extended_weight_conf = self.extended_weight_configuration(coefficients=coefficients)
        return sum(weight for weight in extended_weight_conf)

    # flip

    def flip(self, i, return_position=False):
        r"""
        Return the facet obtained after flipping position ``i`` in ``self``.

        INPUT:

        - ``i`` -- position in the word `Q` (integer).
        - ``return_position`` -- boolean (default: False) tells
          whether the new position should be returned as well.

        OUTPUT:

        - The new subword complex facet.
        - The new position if ``return_position`` is ``True``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.flip(1)
            (2, 3)
            sage: F.flip(1, return_position=True)
            ((2, 3), 3)
        """
        F = set([j for j in self])
        R = [j for j in self._extended_root_configuration_indices()]
        j = _flip_c(self.parent().group(), F, R, i)
        new_facet = SubwordComplexFacet(self.parent(), F)
        new_facet._extended_root_conf_indices = tuple(R)
        if return_position:
            return new_facet, j
        else:
            return new_facet

    # plot and show

    def plot(self, list_colors=[], labels=[], thickness=3, fontsize=14,
             shift=(0, 0), compact=False, roots=True, **args):
        r"""
        In type `A` or `B`, plot a pseudoline arrangement representing
        the facet ``self``.

        Pseudoline arrangements are graphical representations of
        facets of types A or B subword complexes.

        INPUT:

        - ``list_colors`` -- list (default: ``[]``) to personnalize the colors of the pseudolines.
        - ``labels`` -- list (default: ``[]``) to personnalize the labels of the pseudolines.
        - ``thickness`` -- integer (default: ``3``) for the tnickness of the pseudolines.
        - ``fontsize`` -- integer (default: ``14``) for the size of the font used for labels.
        - ``shift`` -- couple of coordinates (default: ``(0,0)``) to change the origin.
        - ``compact`` -- boolean (default: ``False``) to require a more compact representation.
        - ``roots`` -- boolean (default: ``True``) to print the extended root configuration.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F.plot()
            Graphics object consisting of 26 graphics primitives

            sage: W = CoxeterGroup(['B',3])
            sage: c = W.from_reduced_word([1,2,3])
            sage: Q = c.reduced_word()*2 + W.w0.coxeter_sorting_word(c)
            sage: SC = SubwordComplex(Q, W.w0)
            sage: F = SC[15]; F.plot()
            Graphics object consisting of 53 graphics primitives

        REFERENCES: [PilStu]_
        """
        # check that the type is A or B
        # TODO

        # if any(x is 'A' for x in self.parent().group().series()):
        #     type = 'A'
        # elif any(x is 'B' for x in self.parent().group().series()):
        #     type = 'B'
        # else:
        #     raise ValueError("Plotting is currently only implemented in types A or B")

        # import plot facilities
        from sage.plot.line import line
        from sage.plot.text import text
        from sage.plot.colors import colors
        from sage.combinat.permutation import Permutation

        # get properties
        x = 1
        W = self.parent().group()
        Q = self.parent().word()
        if type == 'A':
            last = W.rank()
        else:
            last = W.rank() - 1
        permutation = Permutation(range(1, last + 2))
        x_max = .5

        # list the pseudolines to be drawn
        pseudolines = [[(shift[0], shift[1] + i), .5] for i in range(last + 1)]
        pseudolines_type_B = [[] for i in range(last + 1)]
        contact_points = []
        root_labels = []
        pseudoline_labels = []
        if labels is not False:
            pseudoline_labels += [(pseudoline,
                                   (shift[0] - .1, shift[1] + pseudoline),
                                   "center") for pseudoline in range(last + 1)]
        if roots:
            extended_root_conf = self.extended_root_configuration()
        for position in range(len(Q)):
            y = W.index_set()[Q[position]]
            if type == 'B' and y == 0:
                pseudoline = permutation(1) - 1
                x = pseudolines[pseudoline].pop()
                if compact:
                    x_max = max(x + 1, x_max)
                else:
                    x = x_max
                    x_max += 1
                if position in self:
                    pseudolines[pseudoline] += [(shift[0] + x + 1,
                                                 shift[1]), x + 1]
                    contact_points += [[(shift[0] + x + .5, shift[1] - .2),
                                        (shift[0] + x + .5, shift[1])]]
                else:
                    pseudolines_type_B[pseudoline] = pseudolines[pseudoline] + [(shift[0] + x + .5, shift[1]), (shift[0] + x + .5, shift[1] - .2)]
                    pseudolines[pseudoline] = [(shift[0] + x + .6, shift[1] - .2), (shift[0] + x + .6, shift[1]), .5]
                if roots:
                    root_labels.append((extended_root_conf[position],
                                        (shift[0] + x + .25, shift[1] - .2)))
            else:
                if type == 'B':
                    y -= 1
                pseudoline1 = permutation(y + 1) - 1
                pseudoline2 = permutation(y + 2) - 1
                x = max(pseudolines[pseudoline1].pop(),
                        pseudolines[pseudoline2].pop())
                if compact:
                    x_max = max(x + 1, x_max)
                else:
                    x = x_max
                    x_max += 1
                if position in self:
                    pseudolines[pseudoline1] += [(shift[0] + x + 1,
                                                  shift[1] + y), x + 1]
                    pseudolines[pseudoline2] += [(shift[0] + x + 1,
                                                  shift[1] + y + 1), x + 1]
                    contact_points += [[(shift[0] + x + .5, shift[1] + y),
                                        (shift[0] + x + .5, shift[1] + y + 1)]]
                else:
                    pseudolines[pseudoline1] += [(shift[0] + x + .6,
                                                  shift[1] + y),
                                                 (shift[0] + x + .6,
                                                  shift[1] + y + 1), x + 1]
                    pseudolines[pseudoline2] += [(shift[0] + x + .5,
                                                  shift[1] + y + 1),
                                                 (shift[0] + x + .5,
                                                  shift[1] + y), x + 1]
                    permutation = permutation._left_to_right_multiply_on_left(Permutation((y + 1, y + 2)))
                if roots:
                    root_labels.append((extended_root_conf[position],
                                        (shift[0] + x + .35,
                                         shift[1] + y + .5)))
                if labels is not False:
                    pseudoline_labels += [(pseudoline1, (shift[0] + x + .35,
                                                         shift[1] + y + .05),
                                           "bottom"),
                                          (pseudoline2, (shift[0] + x + .35,
                                                         shift[1] + y + .95),
                                           "top")]

        # transform list to real lines
        list_colors += ['red', 'blue', 'green', 'orange', 'yellow', 'purple']
        list_colors += colors.keys()
        thickness = max(thickness, 2)
        L = line([(1, 1)])
        for contact_point in contact_points:
            L += line(contact_point, rgbcolor=[0, 0, 0],
                      thickness=thickness - 1)
        for pseudoline in range(last + 1):
            pseudolines[pseudoline].pop()
            pseudolines[pseudoline].append((shift[0] + x_max,
                                            shift[1] + permutation.inverse()(pseudoline + 1) - 1))
            L += line(pseudolines[pseudoline], color=list_colors[pseudoline],
                      thickness=thickness)
            if type == 'B':
                L += line(pseudolines_type_B[pseudoline],
                          color=list_colors[pseudoline],
                          thickness=thickness, linestyle="--")
        for root_label in root_labels:
            L += text(root_label[0], root_label[1], rgbcolor=[0, 0, 0],
                      fontsize=fontsize, vertical_alignment="center",
                      horizontal_alignment="right")
        if len(labels) < last + 1:
            labels = range(1, last + 2)
        for pseudoline_label in pseudoline_labels:
            L += text(labels[pseudoline_label[0]], pseudoline_label[1],
                      color=list_colors[pseudoline_label[0]],
                      fontsize=fontsize,
                      vertical_alignment=pseudoline_label[2],
                      horizontal_alignment="right")
        if labels is not False:
            for pseudoline in range(last):
                L += text(labels[pseudoline],
                          (shift[0] + x_max + .1,
                           shift[1] + permutation.inverse()(pseudoline + 1) - 1),
                          color=list_colors[pseudoline], fontsize=fontsize,
                          vertical_alignment="center",
                          horizontal_alignment="left")
        return L

    def show(self, *kwds, **args):
        """
        Show the facet ``self``.

        .. SEEALSO::

            :meth:`plot`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F.show()
        """
        return self.plot().show(axes=False, *kwds, **args)


def _greedy_facet(Q, w, side="negative", n=None, pos=0, l=None, elems=[]):
    r"""
    Return the (positive or negative) *greedy facet* of ``self``.
    """
    if side == "negative":
        pass
    elif side == "positive":
        Q.reverse()
        w = w.inverse()
    else:
        raise ValueError("The optional argument side is not positive "
                         "or negative")

    if n is None:
        n = len(Q)
    if l is None:
        l = w.length()

    if l == 0:
        return elems + range(pos, n)
    elif n < l:
        return []

    s = Q[pos]

    if w.has_left_descent(s):
        X = _greedy_facet(Q, w.apply_simple_reflection_left(s),
                          n=n, pos=pos + 1, l=l - 1, elems=elems)
    else:
        X = []

    if X == []:
        X = _greedy_facet(Q, w, n=n, pos=pos + 1, l=l, elems=elems + [pos])

    if side == "positive":
        X = [n - 1 - i for i in X]
        Q.reverse()
        w = w.inverse()

    return set(X)


def _extended_root_configuration_indices(W, Q, F):
        V_roots = []
        pi = W.one()
        for i in range(len(Q)):
            wi = Q[i]
            V_roots.append((~pi)(W.index_set()[wi] + 1) - 1)
            if i not in F:
                pi = pi.apply_simple_reflection(wi)
        return V_roots


def _greedy_flip_algorithm(Q, w):
    W = w.parent()
    F = _greedy_facet(Q, w, side="positive")
    R = _extended_root_configuration_indices(W, Q, F)
    facet_list = [F]
    extended_root_conf_indices_list = [R]
    flip_to_ancestors = [-1]
    next_index = 0
    while flip_to_ancestors != []:
        has_new_child = False
        for i in sorted(F):
            if (not has_new_child) and (i >= next_index):
                j = _flip_c(W, F, R, i, side="positive")
                if j != i:
                    flip_to_ancestors.append(j)
                    next_index = i + 1
                    has_new_child = True
                    facet_list.append([x for x in F])
                    extended_root_conf_indices_list.append([x for x in R])
        if not has_new_child:
            i = flip_to_ancestors.pop()
            if i != -1:
                j = _flip_c(W, F, R, i, side="negative")
                next_index = j + 1
    return facet_list, extended_root_conf_indices_list
