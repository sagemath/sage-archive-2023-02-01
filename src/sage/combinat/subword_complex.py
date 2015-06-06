r"""
Subword complex

AUTHORS:

- Christian Stump
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
from sage.combinat.combination import Combinations
from sage.geometry.cone import Cone
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.matrix.all import Matrix
from copy import copy
from sage.misc.classcall_metaclass import ClasscallMetaclass

class SubwordComplex(SimplicialComplex,Parent):
    r"""
    Fix a Coxeter system (W,S). The subword complex Delta(Q,w) associated to a word Q in S and an element w in W
    is defined to be the simplicial complex with facets being complements Q/Q' for which Q' is a reduced expression for w

    A subword complex is a shellable sphere if and only if the Demazure product of Q equals w,
    otherwise it is a shellable ball.

    EXAMPLES::

        As an example, dual associahedra are subword complexes in type A_{n-1} given by the
        word [n-1,...,1,n-1,...,1,n-1,...,2,...,n-1,n-2,n-1] and the permutation w_0

        sage: W = CoxeterGroup(['A',2],index_set=[1,2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: C = SubwordComplex([2,1,2,1,2],w); C
        Subword complex of type ['A', 2] for Q = [2, 1, 2, 1, 2] and pi = [1, 2, 1]

        sage: C.facets()
        [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
    """
    def __init__(self, Q, w, g=0, inversion_set_indices=None, algorithm="inductive"):
        W = w.parent()
        I = W.index_set()
        if not all( i in I for i in Q ):
            raise ValueError, "All elements in Q = %s must be contained in the index set %s"%(Q,W.index_set())
        if algorithm != "inductive" and not g == 0:
            raise ValueError("The genus can only be nonzero if the inductive algorithm is used.")
        Vs = range(len(Q))
        self._Q = Q
        self._pi = w
        self._g = g
        if algorithm == "inductive":
            Fs = _construct_facets(Q,w,g=g)
            if g > 0 and inversion_set_indices is not None:
                Fs_new = []
                for F in Fs:
                    ext_root_indices = _extended_root_configuration_indices(W, Q, F)
                    F_inv = [ ext_root_indices[i] for i in range(len(Q)) if i not in F ]
                    print sorted(F_inv)
                    if sorted(F_inv) == sorted(inversion_set_indices):
                        Fs_new.append(F)
                Fs = Fs_new
        elif algorithm == "greedy":
            Fs, Rs = _greedy_flip_algorithm(Q,w)
        else:
            raise ValueError, "The optional argument algorithm can be either inductive or greedy"
        if Fs == []:
            raise ValueError, "The word %s does not contain a reduced expression for %s"%(Q,w.reduced_word())
        SimplicialComplex.__init__(self, vertex_set=Vs, maximal_faces=Fs, vertex_check=False, maximality_check=False)
        self.__custom_name = 'Subword complex'
        self._W = W
        self._cartan_type = W._type
        self._facets_dict = None
        if algorithm == "greedy":
            _facets_dict = {}
            for i in range(len(Fs)):
                X = self(Fs[i], facet_test=False)
                X._extended_root_conf_indices = Rs[i]
                _facets_dict[ tuple(sorted(Fs[i])) ] = X
            self._facets_dict = _facets_dict
        else:
            self._facets_dict = {}

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SubwordComplex([2,1,2,1,2],w)
            Subword complex of type ['A', 2] for Q = [2, 1, 2, 1, 2] and pi = [1, 2, 1]
        """
        repr_str = 'Subword complex of type %s for Q = %s and pi = %s'%(self.cartan_type(),self._Q,self._pi.reduced_word())
        if self._g > 0:
            repr_str += " of genus %s"%self._g
        return repr_str

    def __cmp__(self,other):
        return self.word() == other.word() and self.pi() == other.pi()

    def __call__(self, F, facet_test=True):
        F = tuple(F)
        if self._facets_dict is not None and self._facets_dict != dict() and F in self._facets_dict:
            return self._facets_dict[F]
        else:
            return SubwordComplexFacet(self, F, facet_test=facet_test)

    def __contains__(self,F):
        r"""
        Tests if ``self`` contains a given iterable ``F``.

        EXAMPLES::
        
            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w  = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: S.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
            sage: [0,1] in S
            True
            sage: [0,2] in S
            False
            sage: [0,1,5] in S
            False
            sage: [0] in S
            False
            sage: ['a','b'] in S
            False
        """
        W = self.group()
        Q = self.word()
        if not all( i in range(len(Q)) for i in F ):
            return False
        else:
            return W.from_word( Q[i] for i in range(len(Q)) if i not in F ) == self.pi()

    def list(self):
        return [F for F in self]

    # getting the stored properties

    def group(self):
        r"""
        Returns the group associated to ``self``.

        EXAMPLES::
        
            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: S.group()
            Irreducible finite Coxeter group of rank 2 and type A2
        """
        return self._W

    def cartan_type(self):
        r"""
        Returns the Cartan type of self.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: C = SubwordComplex([2,1,2,1,2],w)
            sage: C.cartan_type()
            ['A', 2]
        """
        from sage.combinat.root_system.cartan_type import CartanType
        if len(self._cartan_type) == 1:
            return CartanType([self._cartan_type[0]['series'],self._cartan_type[0]['rank']])
        else:
            return CartanType(self._cartan_type)

    def word(self):
        r"""
        Returns the word in the simple generators associated to self.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: C = SubwordComplex([2,1,2,1,2],w)
            sage: C.word()
            [2, 1, 2, 1, 2]
        """
        return copy(self._Q)

    def pi(self):
        r"""
        Returns the element in the Coxeter group associated to ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: C = SubwordComplex([2,1,2,1,2],w)
            sage: C.pi().reduced_word()
            [1, 2, 1]
        """
        return self._pi

    def genus(self):
        return self._g

    def facets(self):
        r"""
        Returns all facets of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: S.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
        """
        if self._facets_dict:
            return [ self._facets_dict[tuple(F)] for F in self._facets ]
        else:
            return [ self(F, facet_test=False) for F in self._facets ]

    def __iter__(self):
        r"""
        Returns an iterator on the facets of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: for X in S:
            ...       print X
            (0, 1)
            (0, 4)
            (1, 2)
            (2, 3)
            (3, 4)
        """
        return iter(self.facets())

    def greedy_facet(self,side="positive"):
        r"""
        Returns the negative (or positive) greedy facet of ``self``.

        This is the lexicographically last (or first) facet of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: S.greedy_facet(side="positive")
            (0, 1)
            sage: S.greedy_facet(side="negative")
            (3, 4)
        """
        return SubwordComplexFacet(self,_greedy_facet(self.word(),self.pi(),g=self.genus(),side=side))

    # topological properties

    def is_sphere(self):
        """
        Returns True if the subword complex is a sphere.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3],index_set=[1,2,3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: C = SubwordComplex([3,2,3,2,3],w)
            sage: C.is_sphere()
            True

            sage: C = SubwordComplex([3,2,1,3,2,3],w)
            sage: C.is_sphere()
            False
        """
        W = self._pi.parent()
        w = W.demazure_product(self._Q)
        return w == self._pi

    def is_ball(self):
        """
        Returns True if the subword complex is a ball. This is the case if and only if it is not a sphere.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3],index_set=[1,2,3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: C = SubwordComplex([3,2,3,2,3],w)
            sage: C.is_ball()
            False

            sage: C = SubwordComplex([3,2,1,3,2,3],w)
            sage: C.is_ball()
            True
        """
        return not self.is_sphere()

    def is_pure(self):
        """
        Returns True since subword complexes are pure in general.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3],index_set=[1,2,3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: C = SubwordComplex([3,2,3,2,3],w)
            sage: C.is_pure()
            True
        """
        return True

    def dimension(self):
        """
        Returns the dimension of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: C = SubwordComplex([1,2,1,2,1],W.w0)
            sage: C.dimension()
            1
        """
        return self._facets[0].dimension()

    @cached_method
    def is_root_independent(self):
        """
        Returns True if ``self`` is root independent. This is, if the root configuration of any (or equivalently all) facets is linearly independent.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: C = SubwordComplex([1,2,1,2,1],W.w0)
            sage: C.is_root_independent()
            True

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: C = SubwordComplex([1,2,1,2,1,2],W.w0)
            sage: C.is_root_independent()
            False
        """
        from sage.matrix.all import matrix
        M = matrix(self.greedy_facet(side="negative").root_configuration())
        return M.rank() == max( M.ncols(), M.nrows() )

    @cached_method
    def is_double_root_free(self):
        """
        Returns True if ``self`` is double root free. This is, if the root configurations of all facets do not contain a root twice.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: C = SubwordComplex([1,2,1,2,1],w)
            sage: C.is_double_root_free()
            True

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: C = SubwordComplex([1,1,2,2,1,1],w)
            sage: C.is_double_root_free()
            True
            
            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: C = SubwordComplex([1,2,1,2,1,2],w)
            sage: C.is_double_root_free()
            False
        """
        if not self.is_root_independent():
            size = self.dimension() + 1
            for F in self:
                conf = F._root_configuration_indices()
                if len( set( conf ) ) < size:
                    return False
        return True

    # root and weight properties

    @cached_method
    def word_inversions(self):
        """
        FIX: New name???
        Returns the list of roots `w(\alpha_{q_i})` where `w` is the prefix of the ``self.word()``
        until position `i-1` applied on the `i`-th letter `q_i` of ``self.word()``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: S = SubwordComplex([1,2,1,2,1],W.w0)
            sage: S.word_inversions()
            [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1)]
        """
        W = self._W
        Q = self._Q
        Phi = W.roots()
        pi = W.identity()
        roots = []
        for i in range(len(Q)):
            wi = Q[i]
            roots.append( Phi[(~pi)(W._index_set[wi]+1)-1] )
            pi = pi.apply_simple_reflection(wi)
        return roots

    def kappa_preimages(self):
        """
        Returns a dictionary containing facets of ``self`` as keys,
        and list of elements of ``self.group()`` as values FIX: tba.

        FIX: Add seealso

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([1,2,1,2,1],w)
            sage: kappa_pres = S.kappa_preimages()
            sage: for F in S: print F, [ w.reduced_word() for w in kappa_pres[F] ]
            (0, 1) [[]]
            (0, 4) [[2, 1], [2]]
            (1, 2) [[1]]
            (2, 3) [[1, 2]]
            (3, 4) [[1, 2, 1]]
        """
        return dict( ( F, F.kappa_preimage() ) for F in self )

    def brick_vectors(self,coefficients=None):
        """
        Returns the list of all brick_vectors of facets of ``self``.

        FIX: Add seealso

        EXAMPLES::

            sage: tba
        """
        return [ F.brick_vector(coefficients=coefficients) for F in self ]

    def minkowski_summand(self,i):
        """
        Returns the ``i``th Minkowski summand of ``self``.

        EXAMPLES::

            sage: tba
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        if not self.group().is_crystallographic():
            from sage.rings.all import CC,QQ
            print "Caution: the polytope is build with rational vertices."
            BV = [ [ QQ(CC(v).real()) for v in V ] for V in BV ]
        return Polyhedron([F.extended_weight_configuration()[i] for F in self])

    def brick_polytope(self,coefficients=None):
        """
        Returns the `brick polytope of self.
        FIX: explain coefficients

        EXAMPLES::
        
            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: S = SubwordComplex([1,2,1,2,1],W.w0)

            sage: X = S.brick_polytope(); X
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices

            sage: Y = S.brick_polytope(coefficients=[1,2]); Y
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices

            sage: X == Y
            False
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        BV = self.brick_vectors(coefficients=coefficients)
        if not self.group().is_crystallographic():
            from sage.rings.all import CC,QQ
            print "Caution: the polytope is build with rational vertices."
            BV = [ [ QQ(CC(v).real()) for v in V ] for V in BV ]
        return Polyhedron( BV )

    def word_orbit_decomposition(self):
        """
        FIX: what was that???
        """
        Q = self.word()
        positions = range(len(Q))
        W = self.group()
        N = W.nr_reflections()
        w = W.w0
        I = W.index_set()
        I_decomp = set( [ frozenset( [ i, I[ w( W._index_set[i]+1 ) - 1 - N ] ] ) for i in I ] )
        return [ [ i for i in positions if Q[i] in I_part ] for I_part in I_decomp ]

    def baricenter(self):
        """
        Returns the baricenter of the brick polytope of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: S = SubwordComplex([1,2,1,2,1],W.w0)
            sage: S.baricenter()
            (2/3, 4/3)
        """
        facets = self.facets()
        if not self.is_root_independent():
            facets = [ F for F in facets if F.is_vertex() ]
        return sum( F.brick_vector() for F in facets ) / len(facets)

    # cambrian constructions

    def cover_relations(self,label=False):
        N = self.group().nr_reflections()
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
                        covers.append((F,G,i))
                    else:
                        covers.append((F,G))
                    if G not in seen:
                        Fs.add(G)
        return covers

    def increasing_flip_graph(self,label=True):
        from sage.graphs.digraph import DiGraph
        return DiGraph( self.cover_relations(label=label) )

    def interval(self,I,J):
        G = self.increasing_flip_graph()
        paths = G.all_paths(I,J)
        return set( K for path in paths for K in path )

    def increasing_flip_poset(self):
        from sage.combinat.posets.posets import Poset
        cov = self.cover_relations() 
        if not self.is_root_independent():
            Fs = [ F for F in self if F.is_vertex() ]
            cov = [ (a,b) for a,b in cov if a in Fs and b in Fs ]
        return Poset( ((),cov), facade=True )

class SubwordComplexFacet(Simplex,Element):
    r"""
    A facet of a subword complex.
    """
    
    def __init__(self, parent, positions, facet_test=True):
        r"""
        Initializes a facet of the subword complex ``parent``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: S = SubwordComplex([1,2,1,2,1],W.w0)
            sage: F = S([1,2]); F
            (1, 2)
        """
        if facet_test and positions not in parent:
            raise ValueError, "The given iterable %s is not a facet of the subword complex %s"%(positions,parent)
        Element.__init__(self, parent)
        Simplex.__init__(self, sorted(positions))
        self._extended_root_conf_indices = None
        self._extended_weight_conf = None

    def __eq__(self,other):
        return self.tuple() == other.tuple()

    def _extended_root_configuration_indices(self):
        r"""
        Returns the indices of the roots in ``self.group().roots()`` of the extended root configuration.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: F = S([1,2]); F
            (1, 2)
            sage: F._extended_root_configuration_indices()
            [1, 2, 4, 2, 0]
        """
        if self._extended_root_conf_indices is None:
            self._extended_root_conf_indices = _extended_root_configuration_indices(self.parent().group(), self.parent().word(), self)
        return self._extended_root_conf_indices

    def _root_configuration_indices(self):
        r"""
        Returns the indices of the roots in ``self.group().roots()`` of the root configuration.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: F = S([1,2]); F
            (1, 2)
            sage: F._root_configuration_indices()
            [2, 4]
        """
        indices = self._extended_root_configuration_indices()
        return [ indices[i] for i in self ]

    def extended_root_configuration(self):
        r"""
        Returns the extended root configuration.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: F = S([1,2]); F
            (1, 2)
            sage: F.extended_root_configuration()
            [(0, 1), (1, 1), (0, -1), (1, 1), (1, 0)]
        """
        Phi = self.parent().group().roots()
        return [ Phi[i] for i in self._extended_root_configuration_indices() ]

    def root_configuration(self):
        r"""
        Returns the root configuration.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2],index_set=[1,2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: S = SubwordComplex([2,1,2,1,2],w)
            sage: F = S([1,2]); F
            (1, 2)
            sage: F.root_configuration()
            [(1, 1), (0, -1)]
        """
        Phi = self.parent().group().roots()
        return [ Phi[i] for i in self._root_configuration_indices() ]

    def upper_root_configuration(self):
        conf = self._root_configuration_indices()
        W = self.parent().group()
        Phi = W.roots()
        N = len(Phi)/2
        return [ Phi[i-N] for i in conf if i >= N ]

    def extended_weight_configuration(self,coefficients=None):
        """
        Returns the multiset

        EXAMPLES::

            tba
        """
        if coefficients is not None or self._extended_weight_conf is None:
            W = self.parent().group()
            I = W.index_set()
            Lambda = W.fundamental_weights()
            if coefficients is not None:
                Lambda = [ coefficients[i] * Lambda[i] for i in range(len(Lambda)) ]
            Q = self.parent().word()

            Phi = W.roots()

            V_weights = []

            pi = W.identity()
            for i in range(len(Q)):
                wi = Q[i]
                fund_weight = Lambda[W._index_set[wi]]
                V_weights.append( sum(fund_weight[j]*Phi[ (~pi)(j+1)-1 ] for j in range(len(fund_weight)) ) )
                if i not in self:
                    pi = pi.apply_simple_reflection(wi)
            if coefficients is None:
                self._extended_weight_conf = V_weights
            return V_weights
        else:
            return self._extended_weight_conf

    def weight_configuration(self):
        extended_configuration = self.extended_weight_configuration()
        return [ extended_configuration[i] for i in self ]

    def brick_vector(self,coefficients = None):
        weight_conf = self.extended_weight_configuration(coefficients=coefficients)
        return sum(weight for weight in weight_conf)

    def rays_of_root_fan(self):
        ext_root_conf = self.extended_root_configuration()
        root_conf = self.root_configuration()
        M = Matrix( root_conf )
        M_inverse = M.inverse()
        rays = []
        for j in range(len(ext_root_conf)):
            root = ext_root_conf[j]
            if j in self:
#                rays.append(-root)
                rays.append(-root*M_inverse)
            else:
                new_root = root*M_inverse
#                rays.append( vector([( 1 if i < j else -1 ) * new_root[self.tuple().index(i)] for i in self ])*M )
                rays.append( vector([( 1 if i < j else -1 ) * new_root[self.tuple().index(i)] for i in self ]) )
        return rays

    def root_fan(self):
        from sage.geometry.fan import Fan
        root_conf = self.root_configuration()
        return Fan( self.parent().facets(), self.rays_of_root_fan() )

    def root_polytope(self, orbits):
        from sage.geometry.polyhedron.constructor import Polyhedron
        M_inverse = Matrix( self.root_configuration() ).inverse()
        rays = self.rays_of_root_fan()
#        rhocheck = sum( rays[i] for i in range(len(rays)) if i not in self )*M_inverse
        rhocheck = sum( rays[i] for i in range(len(rays)) if i not in self )
#        print rhocheck
        inequalities = []
#        orbits = self.parent().word_orbit_decomposition()
        for x in range(len(rays)):
            for orbit in orbits:
                if x in orbit:
                    i = min([j for j in self if j in orbit])
                    c = rhocheck[self.tuple().index(i)]
                    inequalities.append( [c] + list(rays[x]) )
        return Polyhedron(ieqs=inequalities)

    def kappa_preimage(self):
        W = self.parent().group()
        N = W.nr_reflections()
        root_conf = self._root_configuration_indices()
        return [ w for w in W if all( w(i+1) <= N for i in root_conf ) ]

    def product_of_upper_roots(self):
        W = self.parent().group()
        return W.prod( W.root_to_reflection(beta) for beta in self.upper_root_configuration() )

    @cached_method
    def local_cone(self):
        return Cone(self.root_configuration())

    def is_vertex(self):
        S = self.parent()
        if not S.is_root_independent():
            W = S.group()
            N = W.nr_reflections()
            root_conf = self._root_configuration_indices()
            for w in W:
                if all( w(i+1) <= N for i in root_conf ):
                    return False
        return True

    def flip(self,i,return_position=False):
        F = set([j for j in self])
        R = [j for j in self._extended_root_configuration_indices()]
        j = _flip(self.parent().group(), F, R, i)
        new_facet = SubwordComplexFacet(self.parent(), F)
        new_facet._extended_root_conf_indices = tuple(R)
        if return_position:
            return new_facet, j
        else:
            return new_facet

    def plot(self, compact=False, list_colors=[], thickness=3, roots=True, labels=[], fontsize=14, shift=(0,0), **args):
        if any( x is not 'A' for x in self.parent().group().series() ):
            raise ValueError, "Plotting is currently only implemented in type A"
        from sage.plot.line import line 
        from sage.plot.text import text 
        from sage.plot.colors import colors
        from sage.combinat.permutation import Permutation
        x = 1
        W = self.parent().group()
        Q = self.parent().word()
        n = W.rank()
        permutation = Permutation(range(1,n+2))
        pseudolines = [[(shift[0]+0,shift[1]+i),.5] for i in range(n+1)]
        x_max = .5
        contact_points = []
        root_labels = []
        pseudoline_labels = []
        if labels is not False:
            pseudoline_labels += [(pseudoline, (shift[0]-.1, shift[1]+pseudoline), "center") for pseudoline in range(n+1)]
        if roots:
            extended_root_conf = self.extended_root_configuration()
        for position in range(len(Q)):
            y = W._index_set[Q[position]]
            pseudoline1 = permutation(y+1)-1
            pseudoline2 = permutation(y+2)-1
            if compact:
                x = max(pseudolines[pseudoline1].pop(), pseudolines[pseudoline2].pop())
                x_max = max(x+1, x_max)
            else:
                pseudolines[pseudoline1].pop()
                pseudolines[pseudoline2].pop()
                x = x_max
                x_max += 1
            if position in self:
                pseudolines[pseudoline1] += [(shift[0]+x+1,shift[1]+y), x+1]
                pseudolines[pseudoline2] += [(shift[0]+x+1,shift[1]+y+1), x+1]
                contact_points += [[(shift[0]+x+.5,shift[1]+y), (shift[0]+x+.5, shift[1]+y+1)]]                
            else:
                pseudolines[pseudoline1] += [(shift[0]+x+.6,shift[1]+y), (shift[0]+x+.6,shift[1]+y+1), x+1]
                pseudolines[pseudoline2] += [(shift[0]+x+.5,shift[1]+y+1), (shift[0]+x+.5,shift[1]+y), x+1]
                permutation = permutation._left_to_right_multiply_on_left(Permutation((y+1, y+2)))
            if roots:
                root_labels.append((extended_root_conf[position],(shift[0]+x+.35,shift[1]+y+.5)))
            if labels is not False:
                pseudoline_labels += [(pseudoline1, (shift[0]+x+.35,shift[1]+y+.05), "bottom"), (pseudoline2, (shift[0]+x+.35,shift[1]+y+.95), "top")]
        list_colors += ['red', 'blue', 'green', 'orange', 'yellow', 'purple'] + colors.keys()
        thickness = max(thickness, 2)
        L = line([(1,1)])
        for contact_point in contact_points:
            L += line(contact_point, rgbcolor=[0,0,0], thickness=thickness-1)
        for pseudoline in range(n+1):
            pseudolines[pseudoline].pop()
            pseudolines[pseudoline].append((shift[0]+x_max, shift[1]+permutation.inverse()(pseudoline+1)-1))
            L += line(pseudolines[pseudoline], color=list_colors[pseudoline], thickness=thickness)
        for root_label in root_labels:
            L += text(root_label[0], root_label[1], rgbcolor=[0,0,0], fontsize=fontsize, vertical_alignment="center", horizontal_alignment="right")
        if len(labels) < n+1:
            labels = range(1,n+2)
        for pseudoline_label in pseudoline_labels:
            L += text(labels[pseudoline_label[0]], pseudoline_label[1], color=list_colors[pseudoline_label[0]], fontsize=fontsize, vertical_alignment=pseudoline_label[2], horizontal_alignment="right")
        if labels is not False:
            for pseudoline in range(n+1):
                L += text(labels[pseudoline], (shift[0]+x_max+.1, shift[1]+permutation.inverse()(pseudoline+1)-1), color=list_colors[pseudoline], fontsize=fontsize, vertical_alignment="center", horizontal_alignment="left")
        return L

    def show(self, *kwds, **args):
        return self.plot().show(axes=False,*kwds,**args)

    def show(self, **kwds):
        """
        
        EXAMPLES::
        
        """
        self.plot(**kwds).show(axes=False)

def _construct_facets(Q, w, g=0, n=None, pos=0, l=None):
    r"""
    Returns the list of facets of the subword complex associated to the word Q and the element w in a Coxeter group W.
    """
    if n is None:
        n = len(Q)
    if l is None:
        first = True
        l = w.length()
    else:
        first = False

    if l == 0 and g == 0:
        return [range(pos,n)]
    elif n < l+pos+2*g:
        return []

    s = Q[pos]
    X = _construct_facets(Q,w,g=g,n=n,pos=pos+1,l=l)
    for f in X:
        f.append(pos)

    if w.has_left_descent(s):
        Y = _construct_facets(Q,w.apply_simple_reflection_left(s),g=g,  n=n,pos=pos+1,l=l-1)
        Y = X+Y
    elif g > 0:
        Y = _construct_facets(Q,w.apply_simple_reflection_left(s),g=g-1,n=n,pos=pos+1,l=l+1)
        Y = X+Y
    else:
        Y = X
    if first:
        return sorted([ sorted(x) for x in Y ])
    else:
        return Y

def _greedy_facet(Q,w,g=0,side="negative",n=None, pos=0,l=None,elems=[]):
    r"""
    Returns the (positive or negative) *greedy facet* of ``self``.
    """
    if side == "negative":
        pass
    elif side == "positive":
        Q.reverse()
        w = w.inverse()
    else:
        raise ValueError, "The optional argument side is not positive or negative"

    if n is None:
        n = len(Q)
    if l is None:
        l = w.length()

    if l == 0 and g == 0:
        return elems + range(pos,n)
    elif n < l+pos+2*g:
        return []

    s = Q[pos]

    if w.has_left_descent(s):
        X = _greedy_facet(Q,w.apply_simple_reflection_left(s),g=g,n=n,pos=pos+1,l=l-1,elems=elems)
    elif g > 0:
        X = _greedy_facet(Q,w.apply_simple_reflection_left(s),g=g-1,n=n,pos=pos+1,l=l+1,elems=elems)
    else:
        X = []

    if X == []:
        X = _greedy_facet(Q,w,g=g,n=n,pos=pos+1,l=l,elems=elems+[pos])
        
    if side == "positive":
        X = [ n - 1 - i for i in X ]
        Q.reverse()
        w = w.inverse()

    return set(X)

def _extended_root_configuration_indices(W, Q, F):
        V_roots = []
        pi = W.identity()
        for i in range(len(Q)):
            wi = Q[i]
            V_roots.append((~pi)(W._index_set[wi]+1)-1)
            if i not in F:
                pi = pi.apply_simple_reflection(wi)
        return V_roots

def _flip(W, positions, extended_root_conf_indices, i, side="both"):
    r = extended_root_conf_indices[i]
    nr_ref = W.nr_reflections()
    r_minus = (r + nr_ref) % (2*nr_ref)
    positions.remove(i)
    j = i
    for k in range(len(extended_root_conf_indices)):
        if j == i and k < i and k not in positions and extended_root_conf_indices[k] == r_minus and side != "positive":
            j = k
        elif j == i and k > i and k not in positions and extended_root_conf_indices[k] == r and side != "negative":
            j = k
    positions.add(j)
    if j != i:
        t = W.reflections()[ min(r,r_minus) ]
        for k in range(min(i,j)+1,max(i,j)+1):
            extended_root_conf_indices[k] = t( extended_root_conf_indices[k]+1 ) - 1
    return j

def _greedy_flip_algorithm(Q, w):
    W = w.parent()
    F = _greedy_facet(Q,w,side="positive")
    R = _extended_root_configuration_indices(W, Q, F)
    facet_list = [F]
    extended_root_conf_indices_list = [R]
    flip_to_ancestors = [-1]
    next_index = 0
    while flip_to_ancestors != []:
        has_new_child = False
        for i in sorted(F):
            if (not has_new_child) and (i >= next_index):
                j = _flip(W, F, R, i, side="positive")
                if j != i:
                    flip_to_ancestors.append(j)
                    next_index = i+1
                    has_new_child = True                    
                    facet_list.append([x for x in F])
                    extended_root_conf_indices_list.append([x for x in R])
        if not has_new_child:
            i = flip_to_ancestors.pop()
            if i != -1:
                j = _flip(W, F, R, i, side="negative")
                next_index = j+1
    return facet_list, extended_root_conf_indices_list
