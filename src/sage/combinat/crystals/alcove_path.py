r"""
Alcove paths

AUTHORS:

- Brant Jones (2008): initial version
- Arthur Lubovsky (2013-03-07): rewritten to implement affine type

Special thanks to: Nicolas Borie, Anne Schilling, Travis Scrimshaw, and
Nicolas Thiery.
"""
#*****************************************************************************
# Copyright (C) 2008 Brant Jones <brant at math.ucdavis.edu>
# Copyright (C) 2013 Arthur Lubovsky <alubovsky at albany.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_crystals import FiniteCrystals
from sage.graphs.all import DiGraph
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.all import vector
from sage.rings.integer import Integer
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.misc.misc_c import prod
from sage.categories.sets_cat import Sets
from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from copy import copy
from sage.misc.latex import latex

class CrystalOfAlcovePaths(UniqueRepresentation, Parent):
    r"""
    Crystal of alcove paths generated from a "straight-line" path to the
    negative of a given dominant weight.

    INPUT:

    - ``cartan_type`` -- Cartan type of a finite or affine untwisted root
      system.

    - ``weight`` -- Dominant weight as a list of (integral) coefficients of
      the fundamental weights.

    - ``highest_weight_crystal`` -- (Default: ``True``) If ``True``
      returns the highest weight crystal.  If ``False`` returns an
      object which is close to being isomorphic to the tensor product
      of Kirillov-Reshetikhin crystals of column shape in the
      following sense: We get all the vertices, but only some of the
      edges.  We'll call the included edges pseudo-Demazure.  They are
      all non-zero edges and the 0-edges not at the end of a 0-string
      of edges, i.e.  not those with `f_{0}(b) = b'` with
      `\phi_0(b) =1`.  (Whereas Demazure 0-edges are those that
      are not at the beginning of a zero string) In this case the
      weight `[c_1, c_2, \ldots, c_k]` represents
      `\sum_{i=1}^k c_i \omega_i`.

      .. NOTE::

          If ``highest_weight_crystal`` = ``False``, since we do not
          get the full crystal, ``TestSuite`` will fail on the
          Stembridge axioms.

    .. SEEALSO::

        - :class:`Crystals`

    EXAMPLES:

    The following example appears in Figure 2 of [LP2008]_::

        sage: C = crystals.AlcovePaths(['G',2],[0,1])
        sage: G = C.digraph()
        sage: GG = DiGraph({
        ....:     ()        : {(0)         : 2 },
        ....:     (0)       : {(0,8)       : 1 },
        ....:     (0,1)     : {(0,1,7)     : 2 },
        ....:     (0,1,2)   : {(0,1,2,9)   : 1 },
        ....:     (0,1,2,3) : {(0,1,2,3,4) : 2 },
        ....:     (0,1,2,6) : {(0,1,2,3)   : 1 },
        ....:     (0,1,2,9) : {(0,1,2,6)   : 1 },
        ....:     (0,1,7)   : {(0,1,2)     : 2 },
        ....:     (0,1,7,9) : {(0,1,2,9)   : 2 },
        ....:     (0,5)     : {(0,1)       : 1, (0,5,7) : 2 },
        ....:     (0,5,7)   : {(0,5,7,9)   : 1 },
        ....:     (0,5,7,9) : {(0,1,7,9)   : 1 },
        ....:     (0,8)     : {(0,5)       : 1 },
        ....:     })
        sage: G.is_isomorphic(GG)
        True
        sage: for (u,v,i) in G.edges(): print (u.integer_sequence() , v.integer_sequence(), i)
        ([], [0], 2)
        ([0], [0, 8], 1)
        ([0, 1], [0, 1, 7], 2)
        ([0, 1, 2], [0, 1, 2, 9], 1)
        ([0, 1, 2, 3], [0, 1, 2, 3, 4], 2)
        ([0, 1, 2, 6], [0, 1, 2, 3], 1)
        ([0, 1, 2, 9], [0, 1, 2, 6], 1)
        ([0, 1, 7], [0, 1, 2], 2)
        ([0, 1, 7, 9], [0, 1, 2, 9], 2)
        ([0, 5], [0, 1], 1)
        ([0, 5], [0, 5, 7], 2)
        ([0, 5, 7], [0, 5, 7, 9], 1)
        ([0, 5, 7, 9], [0, 1, 7, 9], 1)
        ([0, 8], [0, 5], 1)

    Alcove path crystals are a discrete version of Littelmann paths.
    We verify that the alcove path crystal is isomorphic to the LS
    path crystal::

        sage: C1 = crystals.AlcovePaths(['C',3],[2,1,0])
        sage: g1 = C1.digraph() #long time
        sage: C2 = crystals.LSPaths(['C',3],[2,1,0])
        sage: g2 = C2.digraph() #long time
        sage: g1.is_isomorphic(g2, edge_labels=True) #long time
        True

    The preferred initialization method is via explicit weights rather than a Cartan type
    and the coefficients of the fundamental weights::

        sage: R = RootSystem(['C',3])
        sage: P = R.weight_lattice()
        sage: La = P.fundamental_weights()
        sage: C = crystals.AlcovePaths(2*La[1]+La[2]); C
        Highest weight crystal of alcove paths of type ['C', 3] and weight 2*Lambda[1] + Lambda[2]
        sage: C1==C
        True

    We now explain the data structure::

        sage: C = crystals.AlcovePaths(['A',2],[2,0]) ; C
        Highest weight crystal of alcove paths of type ['A', 2] and weight 2*Lambda[1]
        sage: C._R.lambda_chain()
        [(alpha[1], 0), (alpha[1] + alpha[2], 0), (alpha[1], 1), (alpha[1] + alpha[2], 1)]

    The previous list gives the initial "straight line" path from the
    fundamental alcove `A_o` to its translation  `A_o - \lambda` where
    `\lambda = 2\omega_1` in this example. The initial path for weight
    `\lambda` is called the `\lambda`-chain. This path is constructed from
    the ordered pairs `(\beta, k)`, by crossing the hyperplane orthogonal to
    `\beta` at height `-k`. We can view a plot of this path as follows::

        sage: x=C( () )
        sage: x.plot() # not tested - outputs a pdf

    An element of the crystal is given by a subset of the `\lambda`-chain.
    This subset indicates the hyperplanes where the initial path should be
    folded. The highest weight element is given by the empty subset. ::

        sage: x
        ()
        sage: x.f(1).f(2)
        ((alpha[1], 1), (alpha[1] + alpha[2], 1))
        sage: x.f(1).f(2).integer_sequence()
        [2, 3]
        sage: C([2,3])
        ((alpha[1], 1), (alpha[1] + alpha[2], 1))
        sage: C([2,3]).is_admissible() #check if a valid vertex
        True
        sage: C([1,3]).is_admissible() #check if a valid vertex
        False

    Alcove path crystals now works in affine type (:trac:`14143`)::

        sage: C = crystals.AlcovePaths(['A',2,1],[1,0,0]) ; C
        Highest weight crystal of alcove paths of type ['A', 2, 1] and weight Lambda[0]
        sage: x=C(  () )
        sage: x.f(0)
        ((alpha[0], 0),)
        sage: C.R
        Root system of type ['A', 2, 1]
        sage: C.weight
        Lambda[0]

    Test that the tensor products of Kirillov-Reshetikhin crystals
    minus non-pseudo-Demazure arrows is in bijection with alcove path
    construction::

        sage: K = crystals.KirillovReshetikhin(['B',3,1],2,1)
        sage: T = crystals.TensorProduct(K,K)
        sage: g = T.digraph() #long time
        sage: for e in g.edges(): #long time
        ....:     if e[0].phi(0) == 1 and e[2] == 0: #long time
        ....:         g.delete_edge(e)  #long time

        sage: C = crystals.AlcovePaths(['B',3,1],[0,2,0], highest_weight_crystal=False)
        sage: g2 = C.digraph_fast() #long time
        sage: g.is_isomorphic(g2, edge_labels = True) #long time
        True

    .. NOTE::

        In type `C_n^{(1)}`, the Kirillov-Reshetikhin crystal is not connected
        when restricted to pseudo-Demazure arrows, hence the previous example will
        fail for type `C_n^{(1)}` crystals.

    ::

        sage: R = RootSystem(['B',3])
        sage: P = R.weight_lattice()
        sage: La = P.fundamental_weights()
        sage: D = crystals.AlcovePaths(2*La[2], highest_weight_crystal=False)
        sage: C == D
        True

    .. WARNING:: Weights from finite root systems index non-highest weight crystals.

    REFERENCES:

    .. [LP2008]  C. Lenart and A. Postnikov. *A combinatorial model for
       crystals of Kac-Moody algebras*. Trans. Amer. Math. Soc. 360 (2008),
       4349-4381.
    """

    @staticmethod
    def __classcall_private__(cls, starting_weight, cartan_type=None,
                              highest_weight_crystal=None):
        """
        Classcall to mend the input.

        Internally, the
        :class:`~sage.combinat.crystals.alcove_path.CrystalOfAlcovePaths`
        code works with a ``starting_weight`` that is in the weight space
        associated to the crystal. The user can, however, also input a
        ``cartan_type`` and the coefficients of the fundamental weights as
        ``starting_weight``. This code transforms the input into the right
        format (also necessary for :class:`UniqueRepresentation`).

        TESTS::

            sage: C = crystals.AlcovePaths(['A',2,1], [1,0,0])
            sage: C2 = crystals.AlcovePaths(CartanType(['A',2,1]), (1,0,0))
            sage: C is C2
            True
            sage: R = RootSystem(['B',2,1])
            sage: La = R.weight_space().basis()
            sage: B1 = crystals.AlcovePaths(['B',2,1],[0,0,1])
            sage: B2 = crystals.AlcovePaths(La[2])
            sage: B1 is B2
            True
        """
        if isinstance(cartan_type, bool): # new style signature, optional arguments leak over
            highest_weight_crystal = cartan_type

        elif isinstance(cartan_type, list) or isinstance(cartan_type, tuple): #old style signature
            #switch positional arguments
            cartan_type, starting_weight = CartanType(starting_weight), cartan_type

            if highest_weight_crystal == False:
                cartan_type = cartan_type.classical()

            R = RootSystem(cartan_type)
            P = R.weight_space()
            Lambda = P.basis()
            offset = R.index_set()[Integer(0)]
            starting_weight = P.sum(starting_weight[j-offset]*Lambda[j] for j in R.index_set())

        #set defaults
        if highest_weight_crystal is None:
            highest_weight_crystal = True

        if not starting_weight.is_dominant():
            raise ValueError("{0} is not a dominant weight".format(starting_weight))


        return super(CrystalOfAlcovePaths, cls).__classcall__(cls,
                starting_weight, highest_weight_crystal)


    def __init__(self, starting_weight, highest_weight_crystal):
        r"""
        Initialize ``self``.

        TESTS::

            sage: C = crystals.AlcovePaths(['G',2],[0,1])
            sage: TestSuite(C).run()

            sage: C = crystals.AlcovePaths(['A',2,1],[1,0,0])
            sage: TestSuite(C).run() #long time

            sage: C = crystals.AlcovePaths(['A',2,1],[1,0],False)
            sage: TestSuite(C).run(skip="_test_stembridge_local_axioms") #long time
        """
        ##########################################################################
        # NOTE:
        # If cartan_type.is_affine() == True and highest weight crystal == False,
        # since we only use the positive roots of the *finite* root system
        # to get the crystal we set self._finite_cartan_type is true
        #
        # We want the indexing set to include 0 so use the affine type notation
        # for the Cartan type.
        ##########################################################################
        cartan_type = starting_weight.parent().cartan_type()

        self.weight = starting_weight
        self.R = RootSystem(cartan_type)
        self._highest_weight_crystal = highest_weight_crystal
        self._cartan_type = cartan_type


        if cartan_type.is_finite() and highest_weight_crystal == True:
            Parent.__init__(self, category = FiniteCrystals() )
            self._R = RootsWithHeight(starting_weight)
            self._finite_cartan_type = True
        elif cartan_type.is_finite() and highest_weight_crystal == False:
            Parent.__init__(self, category = FiniteCrystals() )
            self._R = RootsWithHeight(starting_weight)
            self._finite_cartan_type = True
            self._cartan_type = cartan_type.affine()
        else:
            Parent.__init__(self, category = HighestWeightCrystals())
            self._R = RootsWithHeight(starting_weight)
            self._finite_cartan_type = False


        self.module_generators = ( self.element_class(self, ()), )

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2,1], [1,0,0])
            sage: C
            Highest weight crystal of alcove paths of type ['A', 2, 1] and weight Lambda[0]
            sage: C = crystals.AlcovePaths(['A',2,1], [1,0], False)
            sage: C
            Crystal of alcove paths of type ['A', 2, 1] and weight Lambda[1]
        """
        if self._highest_weight_crystal:
            return "Highest weight crystal of alcove paths of type %s and weight %s"%(self._cartan_type, self.weight)
        return "Crystal of alcove paths of type %s and weight %s"%(self._cartan_type, self.weight)

    def _element_constructor_(self, data):
        """
        Construct an element of ``self`` from ``data``.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[3,2])
            sage: C([8,9])
            ((alpha[1], 2), (alpha[1] + alpha[2], 4))
        """
        if isinstance(data, tuple):
            return self.element_class(self, data)
        elif isinstance(data, list):
            lambda_chain = self._R.lambda_chain()
            #data starts indexing at 0
            return self.element_class(self, tuple(sorted([lambda_chain[i] for i in data])))

    def vertices(self):
        """
        Return a list of all the vertices of the crystal.

        The vertices are represented as lists of integers recording the folding
        positions.

        One can compute all vertices of the crystal by finding all the
        admissible subsets of the `\lambda`-chain  (see method
        is_admissible, for definition).  We use the breath first
        search algorithm.

        .. WARNING::

            This method is (currently) only useful for the case when
            ``highest_weight_crystal = False``, where you cannot always
            reach all vertices of the crystal using crystal operators,
            starting from the highest weight vertex.  This method is
            typically slower than generating the crystal graph using
            crystal operators.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['C',2],[1,0])
            sage: C.vertices()
            [[], [0], [0, 1], [0, 1, 2]]
            sage: C = crystals.AlcovePaths(['C',2,1],[2,1],False)
            sage: len(C.vertices())
            80

        The number of elements reachable using the crystal operators from the
        module generator::

            sage: len(list(C))
            55
        """
        lambda_chain = self._R.lambda_chain()
        len_lambda_chain = len(lambda_chain)
        W = WeylGroup(self._R._cartan_type, prefix='s')
        s = W.simple_reflections()
        highest_weight_crystal = self._highest_weight_crystal

        if highest_weight_crystal == True:
            successors = 'bruhat_upper_covers'
        else:
            successors = 'quantum_bruhat_successors'

        # lst contains ordered pairs (w,l) l= list of positions that get
        # you to the word, it needs to be refreshed

        #initialization
        lst=[]
        for i in range(len_lambda_chain):
            associated_reflection = lambda_chain[i].root.associated_reflection()
            if len(associated_reflection) == 1:
                lst.append( (prod([ s[j] for j in associated_reflection ]), [i]) )

        l=copy(lst)

        while True:
            lst2 = []
            for x in lst:
                suc = getattr(x[0], successors)()
                for j in range(x[1][-1]+1, len_lambda_chain):
                    temp = x[0] * prod(
                            [ s[k] for k in lambda_chain[j].root.associated_reflection() ])
                    if temp in suc:
                        lst2.append((temp,x[1]+[j]))
                        l.append((temp,x[1]+[j]))
            if lst2 == []:
                break
            else :
                lst = lst2

        return [ [] ] + [i[1] for i in l]

    def digraph_fast(self, depth=None):
        r"""
        Return the crystal :class:`graph <DiGraph>` with maximum depth
        ``depth`` deep starting at the module generator. Significant speed up
        for highest_weight_crystals of affine type.

        EXAMPLES::

            sage: crystals.AlcovePaths(['A',2], [1,1]).digraph_fast(depth=3)
            Digraph on 7 vertices

        TESTS:

        The following example demonstrates the speed improvement.
        The speedup in non-affine types is small however::

            sage: cartan_type = ['A',2,1] #long time
            sage: weight = [1,1,0] #long time
            sage: depth = 5 #long time
            sage: C = crystals.AlcovePaths(cartan_type, weight) #long time
            sage: %timeit C.digraph_fast(depth) # not tested
            10 loops, best of 3: 171 ms per loop
            sage: %timeit C.digraph(subset=C.subcrystal(max_depth=depth, direction='lower')) #not tested
            1 loops, best of 3: 19.7 s per loop
            sage: G1 = C.digraph_fast(depth) #long time
            sage: G2 = C.digraph(subset=C.subcrystal(max_depth=depth, direction='lower')) #long time
            sage: G1.is_isomorphic(G2, edge_labels=True) #long time
            True

        """
        if not self._highest_weight_crystal:
            return super(CrystalOfAlcovePaths, self).digraph()

        if self._cartan_type.is_affine() and depth is None:
            depth = 10
        I = self.index_set()

        rank = 0
        G = { self.module_generators[0]: {} }
        visited = { self.module_generators[0] }

        while depth is None or rank < depth:
            recently_visited = set()
            for x in visited:
                G.setdefault(x, {}) # does nothing if there's a default
                for i in I:
                    xfi = x.f(i)
                    if xfi is not None:
                        G[x][xfi] = i
                        recently_visited.add(xfi)
            if len(recently_visited) == 0: # No new nodes, nothing more to do
                break
            rank += 1
            visited = recently_visited

        return DiGraph(G)

class CrystalOfAlcovePathsElement(ElementWrapper):
    """
    Crystal of alcove paths element.

    INPUT:

    - ``data`` -- a list of folding positions in the lambda chain (indexing
      starts at 0) or a tuple of :class:`RootsWithHeight` giving folding
      positions in the lambda chain.

    EXAMPLES::

        sage: C = crystals.AlcovePaths(['A',2],[3,2])
        sage: x = C ( () )
        sage: x.f(1).f(2)
        ((alpha[1], 2), (alpha[1] + alpha[2], 4))
        sage: x.f(1).f(2).integer_sequence()
        [8, 9]
        sage: C([8,9])
        ((alpha[1], 2), (alpha[1] + alpha[2], 4))
    """
    def __iter__(self):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[1,0])
            sage: lst = list(C)
            sage: for i in lst[2]: i
            (alpha[1], 0)
            (alpha[1] + alpha[2], 0)
        """
        return iter(self.value)

    def is_admissible(self):
        r"""
        Diagnostic test to check if ``self`` is a valid element of the crystal.

        If ``self.value`` is given by

        .. MATH::

            (\beta_1, i_1), (\beta_2, i_2), \ldots, (\beta_k, i_k),

        for highest weight crystals this checks if the sequence

        .. MATH::

            1 \rightarrow s_{\beta_1} \rightarrow
            s_{\beta_1}s_{\beta_2} \rightarrow \cdots \rightarrow
            s_{\beta_1}s_{\beta_2} \ldots s_{\beta_k}

        is a path in the Bruhat graph. If ``highest_weight_crystal=False``,
        then the method checks if the above sequence is a path in the quantum
        Bruhat graph.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[1,1]); C
            Highest weight crystal of alcove paths of type ['A', 2] and weight Lambda[1] + Lambda[2]
            sage: roots = sorted(list(C._R._root_lattice.positive_roots())); roots
            [alpha[1], alpha[1] + alpha[2], alpha[2]]
            sage: r1 = C._R(roots[0],0); r1
            (alpha[1], 0)
            sage: r2 = C._R(roots[2],0); r2
            (alpha[2], 0)
            sage: r3 = C._R(roots[1],1); r3
            (alpha[1] + alpha[2], 1)
            sage: x = C( ( r1,r2) )
            sage: x.is_admissible()
            True
            sage: x = C( (r3,) ); x
            ((alpha[1] + alpha[2], 1),)
            sage: x.is_admissible()
            False
            sage: C = crystals.AlcovePaths(['C',2,1],[2,1],False)
            sage: C([7,8]).is_admissible()
            True
            sage: C = crystals.AlcovePaths(['A',2],[3,2])
            sage: C([2,3]).is_admissible()
            True

        .. TODO:: Better doctest
        """
        W = WeylGroup(self.parent()._R._cartan_type, prefix='s')
        s = W.simple_reflections()
        highest_weight_crystal = self.parent()._highest_weight_crystal

        if highest_weight_crystal == True:
            successors = 'bruhat_upper_covers'
        else:
            successors = 'quantum_bruhat_successors'

        #start at the identity
        w = W.one()
        for i in self:
            t = prod( [ s[j] for j in  i.root.associated_reflection() ] )
            successor = w * t
            if successor not in getattr(w, successors)():
               return False
            w = successor
        return True

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[1,1])
            sage: C([1,2])._latex_()
            [(\alpha_{1} + \alpha_{2}, 0), (\alpha_{1}, 0)]
        """
        return [ (latex(i.root),i.height) for i in self.value ]

    @cached_in_parent_method
    def integer_sequence(self):
        r"""
        Return a list of integers corresponding to positions in
        the `\lambda`-chain where it is folded.

        .. TODO::

            Incorporate this method into the ``_repr_`` for finite Cartan type.

        .. NOTE::

            Only works for finite Cartan types and indexing starts at 0.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[3,2])
            sage: x = C( () )
            sage: x.f(1).f(2).integer_sequence()
            [8, 9]
        """
        lambda_chain = self.parent()._R.lambda_chain()
        return [lambda_chain.index(j) for j in self.value]

    def phi(self, i):
        r"""
        Return the distance to the end of the `i`-string.

        This method overrides the generic implementation in the category of
        crystals since this computation is more efficient.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[1,1])
            sage: [c.phi(1) for c in C]
            [1, 2, 0, 1, 0, 0, 1, 0]
            sage: [c.phi(2) for c in C]
            [1, 0, 2, 0, 1, 1, 0, 0]
        """
        highest_weight_crystal = self.parent()._highest_weight_crystal
        positions, gi = self._gi(i)

        m=max(gi)

        if not highest_weight_crystal and i == 0:
            raise NotImplementedError
            # I think the M below should still work in this case

        M = Integer(m)/2 - Integer(1)/2
        return M

    def epsilon(self, i):
        r"""
        Return the distance to the start of the `i`-string.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[1,1])
            sage: [c.epsilon(1) for c in C]
            [0, 0, 1, 1, 0, 2, 0, 1]
            sage: [c.epsilon(2) for c in C]
            [0, 1, 0, 0, 1, 0, 2, 1]
        """
        #crude but functional
        j = 0
        temp = self
        temp = temp.e(i)
        while temp is not None:
            j+=1
            temp = temp.e(i)

        return j

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[2,0])
            sage: for i in C: i.weight()
            2*Lambda[1]
            Lambda[2]
            Lambda[1] - Lambda[2]
            -2*Lambda[1] + 2*Lambda[2]
            -Lambda[1]
            -2*Lambda[2]


        """
        root_space = self.parent().R.root_space()
        weight = -self.parent().weight
        for i in self.value[::-1]:
            root = root_space(i.root)
            weight = -i.height*root + weight.reflection(root)

        return -weight

    #def __repr__(self):
        #return str(self.integer_sequence())

    def plot(self):
        r"""
        Return a plot ``self``.

        .. NOTE::

            Currently only implemented for types `A_2`, `B_2`, and `C_2`.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[2,0])
            sage: x = C( () ).f(1).f(2)
            sage: x.plot() # Not tested - creates a pdf
        """
        ct = self.parent()._R._cartan_type.dual()
        word = self.parent()._R.word()
        integer_sequence = self.integer_sequence()
        foldings = [False for i in word]
        for i in integer_sequence:
            foldings[i] = True
        affine_ambient_space = RootSystem(ct.affine()).ambient_space()
        return affine_ambient_space.plot() + affine_ambient_space.plot_alcove_walk( word, foldings=foldings, labels=False)

    def __eq__(self, other):
        r"""
        Test equality of ``self.value`` and ``other.value``.

        EXAMPLES::

            sage: C=crystals.AlcovePaths(['B',2],[1,0])
            sage: lst=list(C)
            sage: lst[2] == lst[2]
            True
            sage: lst[2] == lst[1]
            False
        """
        #note: may want to use _eq_ for coercion
        try:
            return self.value == other.value
        except (NameError, AttributeError):
            return False

    def __lt__(self, other):
        r"""
        Test if ``self.value`` is less than ``other.value`` in dictionary order.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[2,0])
            sage: x = C( () )
            sage: x.__lt__(x.f(1))
            True
            sage: a=x.f(1) ; b = x.f(1).f(1).f(2)
            sage: a.__lt__(b)
            False
        """
        return self.value < other.value

    def __gt__(self, other):
        r"""
        Test if ``self.value`` is greater than ``other.value`` in dictionary
        order.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[2,0])
            sage: x = C( () )
            sage: x.__gt__(x.f(1))
            False
            sage: a=x.f(1) ; b = x.f(1).f(1).f(2)
            sage: a.__gt__(b)
            True
        """
        return self.value > other.value

    def _folding_data(self, i):
        r"""
        Compute information needed to build the graph `g_{\alpha_i}`.
        Results of this method are sent to _gi for further processing.

        INPUT:

        - ``i`` -- element of the index_set of the underlying root_system.

        OUTPUT:

        A dictionary where the keys are of type RootsWithHeight which record
        positions where `\pm \alpha_i` shows up in the folded `\lambda` chain.
        The values are `1` if `\alpha_i` is in the corresponding position in
        the folded `\lambda`-chain, `-1` if `-\alpha_i` is in the corresponding
        position in the folded `\lambda`-chain.

        .. NOTE::

            *infinity* is a special key that records the "sign at infinity".

        ::

            sage: C=crystals.AlcovePaths(['A',2],[1,1])
            sage: x=C( () ).f(1)
            sage: x._folding_data(2)
            {(alpha[1] + alpha[2], 1): 1, 'infinity': 1, (alpha[2], 0): 1}
        """
        Parent = self.parent()

        #self.value contains the admissible sequence as a tuple of Element

        finite_cartan_type = Parent._finite_cartan_type  # bool
        J = list(self.value)

        #NOTE: R is a RootsWithHeight object and NOT a RootSystem object
        R = Parent._R
        weight = Parent.weight

        signs = {}

        # 0 arrows in the case of finite Cartan type
        # always allow 0 arrows
        if finite_cartan_type and i == 0:
            Beta = R._root_lattice.highest_root()
        elif i in self.index_set():
            Beta = R._root_lattice.simple_root(i)

        max_height_Beta = weight.scalar(Beta.associated_coroot())

        if len(J) == 0:
            for k in range( max_height_Beta ) :
                x = R(Beta, k)
                signs[x]=self._sign(Beta)
            signs['infinity'] = self._sign(Beta)

        elif len(J) > 0 :
            #NOTE: we assume J is sorted by order on Element of RootsWithHeight

            for k in  range( max_height_Beta ):
                x = R(Beta, k)
                if x <= J[0]:
                    signs[x] = self._sign(Beta)

            for j in range( len(J)  ):

                Beta = Beta.reflection(J[j].root)
                sign_Beta = self._sign(Beta)
                max_height_Beta = weight.scalar(
                    (sign_Beta * Beta).associated_coroot())


                # some optimization so we don't initialize too many objects
                # range(c1,c2) can be replaced by range(max_height_Beta) but it
                # checks unnecessary extra things

                c1 = J[j]._cmp_v[0] * max_height_Beta
                if j == len(J) - 1:
                    c2 = max_height_Beta
                else:
                    c2 = min (max_height_Beta, J[j+1]._cmp_v[0]*max_height_Beta + 1)

                for k in range(c1,c2):

                    x=R( sign_Beta * Beta , k)

                    if (
                        ( j < len(J) - 1 and J[j] < x <= J[j+1] ) or
                        ( j == len(J) - 1 and J[j] < x)
                    ):
                        signs[x] = sign_Beta

            signs['infinity'] = sign_Beta # tail sign tells something about last step
                                          # in g_alpha

        if finite_cartan_type and i==0:
            signs = { x : -signs[x] for x in signs.keys() }

        return signs

    def e(self, i):
        r"""
        Return the `i`-th crystal raising operator on ``self``.

        INPUT:

        - ``i`` -- element of the index set of the underlying root system.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['A',2],[2,0]); C
            Highest weight crystal of alcove paths of type ['A', 2] and weight 2*Lambda[1]
            sage: x = C( () )
            sage: x.e(1)
            sage: x.f(1) == x.f(1).f(2).e(2)
            True
        """
        Parent = self.parent()
        finite_cartan_type = Parent._finite_cartan_type

        J = list(self.value)
        positions, gi = self._gi(i)

        m=max(gi)
        m_index = len(gi)-1-list(reversed(gi)).index(m) # last max in gi


        if finite_cartan_type and i==0 :
            M = Integer(m)/2 + Integer(1)/2
        else:
            M = Integer(m)/2 - Integer(1)/2


        KR_test = finite_cartan_type and i==0 and m_index < len(gi) - 1
        KR_test = KR_test and M >= 1

        ######################################################################
        # NOTE:
        # In the KR_case we want to insure that positions[m_index] is in J
        # If m_index > 0 then it's always true
        # If m_index == 0 then M >=1 guarantees this
        ######################################################################

        if ( (not finite_cartan_type or i!=0) and m_index < len(gi)-1  # alpha_i is a simple root
            ) or KR_test:


            J.remove(positions[m_index])
            if m_index+1 < len(positions): # if m_index+1 != 'infinity'
                                           # i.e. positions[m_index+1] makes sense
                J.append(positions[m_index+1])
            return_value = Parent ( tuple( sorted(J) ) )

            # we attach to each admissible sequence a list
            # which encodes a path (via root operators) from the () generator
            # to the admissible sequence
            # this is useful for investing the crystal

            try:
                return_value.i_string = self.i_string + [['e',i]]
            except AttributeError:
                return_value.i_string = [['e',i]]

            return return_value
        else:
            return None

    @cached_method
    def _gi(self, i):
        r"""
        Compute information needed to build the graph `g_{\alpha_i}`.
        This graph is used to apply the `i`-th crystal operator.

        INPUT:

        - ``i`` - element of the index_set of the underlying root_system.

        OUTPUT:

        A tuple ``(positions, gi)``:

        - ``positions`` -- is a list of RootsWithHeight. These appear sorted in
          their natural order, and record where  `\pm \alpha_i` shows up in
          the folded `\lambda`-chain.

        - ``gi`` -- is a list of integers recording the height
          (up to affine transformation)  of `\pm \alpha_i`
          in the folded `\lambda`-chain whose location is recorded by
          ``positions``.

        .. NOTE::

            - ``positions`` has length one less than ``gi`` since it does not
              contain the position 'infinity'.

            - To get the real `g_{\alpha_i}` one has to divide by 2 and add 1/2
              or divide by 2 and subtract 1/2 depending on if
              ``self._finite_cartan_type==True and i == 0``
              or not. This is done in crystal operator methods.

        EXAMPLES::

            sage: C=crystals.AlcovePaths(['A',2],[1,1])
            sage: x=C( () ).f(1)
            sage: x._gi(2)
            ([(alpha[2], 0), (alpha[1] + alpha[2], 1)], [1, 3, 5])
        """
        signs = self._folding_data(i)
        positions = sorted( [ x for x in signs.keys() if x != 'infinity' ] )

        if len(positions)==0 :
            return ( positions, [ signs['infinity'] ] )

        gi = [ signs[ positions[0] ] ]
        for j in range(1,len(positions)):
            gi.append(
                    gi[j-1] +
                    signs[positions[j-1]]*self._eps(positions[j-1]) + signs[positions[j]] )
        gi.append(  gi[-1] +
            signs[positions[-1]]*self._eps(positions[-1]) + signs['infinity'] )

        return (positions, gi)

    def f(self, i):
        r"""
        Returns the `i`-th crystal lowering operator on ``self``.

        INPUT:

        - ``i`` -- element of the index_set of the underlying root_system.

        EXAMPLES::

            sage: C=crystals.AlcovePaths(['B',2],[1,1])
            sage: x=C(  () )
            sage: x.f(1)
            ((alpha[1], 0),)
            sage: x.f(1).f(2)
            ((alpha[1], 0), (alpha[1] + alpha[2], 2))

        """
        Parent = self.parent()
        finite_cartan_type = Parent._finite_cartan_type

        # get a copy in a form of a list of self.value
        J = list(self.value)
        positions, gi = self._gi(i)

        m=max(gi)
        m_index=gi.index(m)


        if finite_cartan_type and i==0 :

            # python doesn't handle fractions natively
            M = Integer(m)/2 + Integer(1)/2
        else:
            M = Integer(m)/2 - Integer(1)/2


        # boolian determining when to move a folding in KR case
        KR_test = finite_cartan_type and i==0
        KR_test = KR_test and M > 1

        # In the KR case, we return a value other than None when
        # `\alpha_i` is in position m_index - 1
        # (The following relies on a technical condition (C2) )
        # note insert reference
        #
        # if m_index - 1 == 0 then M > 1 and (C2) forces
        # `\alhpa_i` in positions[m_index - 1]
        #
        # otherwise if m_index - 1 > 0 then (C2) is enough

        if ( (not finite_cartan_type or i!=0) and M > 0  # alpha_i is a simple root
           ) or KR_test :# KR case

            J.append(positions[m_index-1])
            if m_index < len(positions): # if m_index != 'infinity'
                                         # thus positions[m_index] makes sense
                J.remove(positions[m_index])
            return_value = Parent ( tuple( sorted(J) ) )

            # we attach to each admissible sequence a list
            # which encodes a path (via root operators) from the generator ()

            try:
                return_value.i_string = self.i_string + [['f',i]]
            except AttributeError:
                return_value.i_string = [['f',i]]

            return return_value
        else:
            return None

    @staticmethod
    def _sign(root):
        r"""
        Return `1` if root is a positive root, and `-1` if root is a negative
        root.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import CrystalOfAlcovePathsElement
            sage: rl = RootSystem(['A',2]).root_lattice()
            sage: x = rl.from_vector(vector([0,1]))
            sage: CrystalOfAlcovePathsElement._sign(x)
            1
        """
        if root.is_positive_root():
            return 1
        else:
            return -1

    def _eps(self, root):
        r"""
        Return `-1` if root is in ``self.value``, otherwise return `1`.

        EXAMPLES::

            sage: C = crystals.AlcovePaths(['C',2],[3,2])
            sage: x = C( () ).f(1).f(2); x
            ((alpha[1], 2), (2*alpha[1] + alpha[2], 4))
            sage: x._eps(x.value[0])
            -1
            sage: R = C._R
            sage: y = R ( x.value[0].root, 1 ); y
            (alpha[1], 1)
            sage: x._eps(y)
            1

        """
        if root in self.value:
            return -1
        else:
            return 1

CrystalOfAlcovePaths.Element = CrystalOfAlcovePathsElement
#deprecate the old name
from sage.misc.superseded import deprecated_function_alias
ClassicalCrystalOfAlcovePaths = deprecated_function_alias(14143, CrystalOfAlcovePaths)

class RootsWithHeight(UniqueRepresentation, Parent):
    r"""
    Data structure of the ordered pairs `(\beta,k)`,
    where `\beta` is a positive root and `k` is a non-negative integer. A total
    order is implemented on this set, and depends on the weight.

    INPUT:

    - ``cartan_type`` -- Cartan type of a finite or affine untwisted root
      system

    - ``weight`` -- dominant weight as a list of (integral) coefficients of
      the fundamental weights

    EXAMPLES::

        sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
        sage: R = RootsWithHeight(['A',2],[1,1]); R
        Roots with height of Cartan type ['A', 2] and dominant weight Lambda[1] + Lambda[2]

        sage: r1 = R._root_lattice.from_vector(vector([1,0])); r1
        alpha[1]
        sage: r2 = R._root_lattice.from_vector(vector([1,1])); r2
        alpha[1] + alpha[2]

        sage: x = R(r1,0); x
        (alpha[1], 0)
        sage: y = R(r2,1); y
        (alpha[1] + alpha[2], 1)
        sage: x < y
        True
    """

    @staticmethod
    def __classcall_private__(cls, starting_weight, cartan_type = None):
        """
        Classcall to mend the input.

        Internally, the RootsWithHeight code works with a ``starting_weight`` that
        is in the ``weight_space`` associated to the crystal. The user can, however,
        also input a ``cartan_type`` and the coefficients of the fundamental weights
        as ``starting_weight``. This code transforms the input into the right
        format (also necessary for UniqueRepresentation).

        TESTS::
            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: S = RootsWithHeight(CartanType(['A',2]), (3,2))
            sage: R is S
            True

            sage: R = RootSystem(['B',2,1])
            sage: La = R.weight_space().basis()
            sage: C = RootsWithHeight(['B',2,1],[0,0,1])
            sage: B = RootsWithHeight(La[2])
            sage: B is C
            True

        """
        if cartan_type is not None:
            cartan_type, starting_weight = CartanType(starting_weight), cartan_type

            R = RootSystem(cartan_type)
            P = R.weight_space()
            Lambda = P.basis()
            offset = R.index_set()[Integer(0)]
            starting_weight = P.sum(starting_weight[j-offset]*Lambda[j] for j in R.index_set())

        return super(RootsWithHeight, cls).__classcall__(cls, starting_weight)


    def __init__(self, weight):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: TestSuite(R).run()
        """
        Parent.__init__(self, category = Sets() )

        cartan_type = weight.parent().cartan_type()
        self._cartan_type = cartan_type
        self._root_system = RootSystem(cartan_type)
        self._root_lattice = self._root_system.root_lattice()
        self._weight_lattice = self._root_system.weight_lattice()
        self.weight = weight

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: RootsWithHeight(['A',2],[3,2])
            Roots with height of Cartan type ['A', 2] and dominant weight 3*Lambda[1] + 2*Lambda[2]
        """
        return "Roots with height of Cartan type %s and dominant weight %s"%(
            self._root_system.cartan_type(), self.weight)

    def _max_height(self, root):
        r"""
        If root is `\beta`, return `k = \langle \lambda, \beta^{\vee} \rangle`.

        Only ordered pairs of the form `(\beta, l)` for `0 \leq l < k` are
        allowed.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: C = RootsWithHeight(['A',3],[3,2,0])
            sage: x = C._root_lattice.from_vector(vector([1,1])); x
            alpha[1] + alpha[2]
            sage: C._max_height(x)
            5
        """
        return self.weight.scalar(root.associated_coroot())

    @cached_method
    def word(self):
        r"""
        Gives the initial alcove path (`\lambda`-chain) in terms of simple
        roots. Used for plotting the path.

        .. NOTE::

            Currently only implemented for finite Cartan types.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: R.word()
            [2, 1, 2, 0, 1, 2, 1, 0, 1, 2]
        """
        cartan_type = self._root_system.cartan_type()
        if not cartan_type.is_finite():
            raise NotImplementedError
        lambda_chain = [ x.root for x in self.lambda_chain() ]

        coroot_lattice = RootSystem(cartan_type).coroot_lattice()
        cohighest_root = coroot_lattice.highest_root()

        word = []
        for i in range(len(lambda_chain)):
            beta = lambda_chain[i]
            for j in reversed(range(i)):
                beta = beta.reflection(lambda_chain[j])
            #beta is now a simple root or the highest root

            coroot = beta.associated_coroot()
            support = coroot.support() # the path is in dual affine space
            if len(support) == 1: # beta is a simple root
                word.append(support[0])
            elif coroot == -cohighest_root:
                word.append(0)
            else:
                assert False, 'should never get here'

        return word

    @cached_method
    def lambda_chain(self):
        r"""
        Return the unfolded `\lambda`-chain.

        .. NOTE:: Only works in root systems of finite type.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[1,1]); R
            Roots with height of Cartan type ['A', 2] and dominant weight Lambda[1] + Lambda[2]
            sage: R.lambda_chain()
            [(alpha[2], 0), (alpha[1] + alpha[2], 0), (alpha[1], 0), (alpha[1] + alpha[2], 1)]
        """
        if not self._root_lattice.cartan_type().is_finite():
            raise ValueError("Cartan type {0} is not finite".format(self._root_lattice.cartan_type()))

        l=[]
        for i in self._root_lattice.positive_roots():
            for j in range(self._max_height(i)):
                l.append(self(i,j))

        return sorted(l)

    def _element_constructor_(self, root, height):
        r"""
        Construct a :class:`RootsWithHeightElement` with ``self`` as the parent.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: rl = RootSystem(['A',2]).root_lattice()
            sage: x = rl.from_vector(vector([1,1])); x
            alpha[1] + alpha[2]
            sage: R = RootsWithHeight(['A',2],[1,1]); R
            Roots with height of Cartan type ['A', 2] and dominant weight Lambda[1] + Lambda[2]
            sage: y = R(x,1); y
            (alpha[1] + alpha[2], 1)
        """
        root = self._root_lattice.from_vector(vector(root))
        return self.element_class(self, root, height)

    def _an_element_(self):
        r"""

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: R._an_element_()
            (alpha[1], 0)
        """
        return self( self._root_lattice.from_vector(vector([1])), 0 )

class RootsWithHeightElement(Element):
    r"""
    Element of :class:`RootsWithHeight`.

    INPUT:

    - ``root`` -- A positive root `\beta` in our root system
    - ``height`` -- Is an integer, such that
      `0 \leq l \leq \langle \lambda, \beta^{\vee} \rangle`

    EXAMPLES::

        sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
        sage: rl = RootSystem(['A',2]).root_lattice()
        sage: x = rl.from_vector(vector([1,1])); x
        alpha[1] + alpha[2]
        sage: R = RootsWithHeight(['A',2],[1,1]); R
        Roots with height of Cartan type ['A', 2] and dominant weight Lambda[1] + Lambda[2]
        sage: y = R(x, 1); y
        (alpha[1] + alpha[2], 1)
    """
    def __init__(self, parent, root, height):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: rl = RootSystem(['A',2]).root_lattice()
            sage: x = rl.from_vector(vector([1,1]))
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: y = R(x, 1); y
            (alpha[1] + alpha[2], 1)
            sage: TestSuite(x).run()
        """
        Element.__init__(self, parent)
        max_height = parent._max_height(root)

        # make sure the height is in the right range, this also catches negative
        # roots

        if not 0 <= height < max_height:
            raise ValueError("%d out of allowed range [%d,%d)"%(height, 0, max_height))

        v = [height/max_height]
        v.extend( [ x/max_height for x in root.associated_coroot().to_vector() ] )
        #v.insert(0, height/max_height)

        # the map from (root, height) --> _cmp_v is injective

        self._cmp_v = tuple(v)
        self.root = root
        self.height = height

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: rl = RootSystem(['A',2]).root_lattice()
            sage: vec = rl.from_vector(vector([1,1])); vec
            alpha[1] + alpha[2]
            sage: R(vec,1)
            (alpha[1] + alpha[2], 1)
        """
        return "(%s, %s)" % (self.root, self.height)

    def __hash__(self):
        r"""

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: rl = RootSystem(['A',2]).root_lattice()
            sage: root = rl.from_vector(vector([1,1]))
            sage: vec = R(root,0)
            sage: hash(vec) == hash(vec)
            True
        """
        return hash(self._cmp_v)

    def __eq__(self, other):
        r"""

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: rl = RootSystem(['A',2]).root_lattice()
            sage: v1 = rl.from_vector(vector([1,1]))
            sage: v2 = rl.from_vector(vector([1]))
            sage: x1 = R(v1,1) ; x2 = R(v1,0) ; x3 = R(v2,1)
            sage: x1.__eq__(x1)
            True
            sage: x1.__eq__(x2)
            False
            sage: x1.__eq__(x3)
            False
        """
        try:
            return self._cmp_v == other._cmp_v
        except (NameError, AttributeError):
            return False

    def __cmp__(self, other):
        r"""
        Define a total order on :class:`RootsWithHeightElement`. This defines
        the initial `\lambda`-chain.

        EXAMPLES::

            sage: from sage.combinat.crystals.alcove_path import RootsWithHeight
            sage: R = RootsWithHeight(['A',2],[3,2])
            sage: rl = RootSystem(['A',2]).root_lattice()
            sage: v1 = rl.from_vector(vector([1,1]))
            sage: v2 = rl.from_vector(vector([1]))
            sage: x1 = R(v1,1) ; x2 = R(v1,0) ; x3 = R(v2,1)
            sage: x1.__cmp__(x2)
            1
            sage: x1.__cmp__(x3)
            -1

        """
        # I suspect that if you redefine this method to produce a
        # different (valid)  `\lambda`-chain the rest of the
        # code should still work.
        #todo: check if self and other have the same parent ?
        #assert self.parent() is other.parent(), "elements have different parents"
        return cmp(self._cmp_v, other._cmp_v)

RootsWithHeight.Element = RootsWithHeightElement

#####################################################################
# Test code, by comparing with existing crystal implementations.
#####################################################################

def _test_some_specific_examples(clss=CrystalOfAlcovePaths):
    r"""
    Test against some specific (finite type) examples.

    EXAMPLES::

        sage: from sage.combinat.crystals.alcove_path import _test_some_specific_examples
        sage: _test_some_specific_examples(crystals.AlcovePaths)
        G2 example passed.
        C3 example passed.
        B3 example 1 passed.
        B3 example 2 passed.
        True
    """
    # This appears in Lenart.
    C = clss(['G',2],[0,1])
    G = C.digraph()

    GT = DiGraph({
        ()        : {(0)         : 2 },
        (0)       : {(0,8)       : 1 },
        (0,1)     : {(0,1,7)     : 2 },
        (0,1,2)   : {(0,1,2,9)   : 1 },
        (0,1,2,3) : {(0,1,2,3,4) : 2 },
        (0,1,2,6) : {(0,1,2,3)   : 1 },
        (0,1,2,9) : {(0,1,2,6)   : 1 },
        (0,1,7)   : {(0,1,2)     : 2 },
        (0,1,7,9) : {(0,1,2,9)   : 2 },
        (0,5)     : {(0,1)       : 1, (0,5,7) : 2 },
        (0,5,7)   : {(0,5,7,9)   : 1 },
        (0,5,7,9) : {(0,1,7,9)   : 1 },
        (0,8)     : {(0,5)       : 1 }
        })

    if (G.is_isomorphic(GT) != True):
        return False
    else:
        print "G2 example passed."

    # Some examples from Hong--Kang:

    # type C, ex. 8.3.5, pg. 189
    C = clss(['C',3],[0,0,1])
    G = C.digraph()
    GT = DiGraph({
        ():{ (0): 3},
        (0):{ (0, 6): 2},
        (0, 1):{ (0, 1, 3): 3, (0, 1, 7): 1},
        (0, 1, 2):{ (0, 1, 2, 3): 3},
        (0, 1, 2, 3):{ (0, 1, 2, 3, 8): 2},
        (0, 1, 2, 3, 4):{ (0, 1, 2, 3, 4, 5): 3},
        (0, 1, 2, 3, 8):{ (0, 1, 2, 3, 4): 2},
        (0, 1, 3):{ (0, 1, 3, 7): 1},
        (0, 1, 3, 7):{ (0, 1, 2, 3): 1, (0, 1, 3, 7, 8): 2},
        (0, 1, 3, 7, 8):{ (0, 1, 2, 3, 8): 1},
        (0, 1, 7):{ (0, 1, 2): 1, (0, 1, 3, 7): 3},
        (0, 6):{ (0, 1): 2, (0, 6, 7): 1},
        (0, 6, 7):{ (0, 1, 7): 2}
        })

    if (G.is_isomorphic(GT) != True):
        return False
    else:
        print "C3 example passed."

    # type B, fig. 8.1 pg. 172
    C = clss(['B',3],[2,0,0])
    G = C.digraph()

    GT = DiGraph({
        ():{ (6): 1},
        (0):{ (0, 7): 2},
        (0, 1):{ (0, 1, 11): 3},
        (0, 1, 2):{ (0, 1, 2, 9): 2},
        (0, 1, 2, 3):{ (0, 1, 2, 3, 10): 1},
        (0, 1, 2, 3, 10):{ (0, 1, 2, 3, 4): 1},
        (0, 1, 2, 9):{ (0, 1, 2, 3): 2, (0, 1, 2, 9, 10): 1},
        (0, 1, 2, 9, 10):{ (0, 1, 2, 3, 10): 2},
        (0, 1, 5):{ (0, 1, 2): 3, (0, 1, 5, 9): 2},
        (0, 1, 5, 9):{ (0, 1, 2, 9): 3, (0, 1, 5, 9, 10): 1},
        (0, 1, 5, 9, 10):{ (0, 1, 2, 9, 10): 3},
        (0, 1, 8):{ (0, 1, 5): 3},
        (0, 1, 8, 9):{ (0, 1, 5, 9): 3, (0, 1, 8, 9, 10): 1},
        (0, 1, 8, 9, 10):{ (0, 1, 5, 9, 10): 3},
        (0, 1, 11):{ (0, 1, 8): 3},
        (0, 7):{ (0, 1): 2, (0, 7, 11): 3},
        (0, 7, 8):{ (0, 7, 8, 9): 2},
        (0, 7, 8, 9):{ (0, 1, 8, 9): 2},
        (0, 7, 8, 9, 10):{ (0, 1, 8, 9, 10): 2},
        (0, 7, 11):{ (0, 1, 11): 2, (0, 7, 8): 3},
        (6):{ (0): 1, (6, 7): 2},
        (6, 7):{ (0, 7): 1, (6, 7, 11): 3},
        (6, 7, 8):{ (0, 7, 8): 1, (6, 7, 8, 9): 2},
        (6, 7, 8, 9):{ (6, 7, 8, 9, 10): 1},
        (6, 7, 8, 9, 10):{ (0, 7, 8, 9, 10): 1},
        (6, 7, 11):{ (0, 7, 11): 1, (6, 7, 8): 3}
        })

    if (G.is_isomorphic(GT) != True):
        return False
    else:
        print "B3 example 1 passed."

    C = clss(['B',3],[0,1,0])
    G = C.digraph()

    GT = DiGraph({
        ():{ (0): 2},
        (0):{ (0, 1): 1, (0, 7): 3},
        (0, 1):{ (0, 1, 7): 3},
        (0, 1, 2):{ (0, 1, 2, 8): 2},
        (0, 1, 2, 3):{ (0, 1, 2, 3, 5): 1, (0, 1, 2, 3, 9): 3},
        (0, 1, 2, 3, 4):{ (0, 1, 2, 3, 4, 5): 1},
        (0, 1, 2, 3, 4, 5):{ (0, 1, 2, 3, 4, 5, 6): 2},
        (0, 1, 2, 3, 5):{ (0, 1, 2, 3, 5, 9): 3},
        (0, 1, 2, 3, 5, 9):{ (0, 1, 2, 3, 4, 5): 3},
        (0, 1, 2, 3, 9):{ (0, 1, 2, 3, 4): 3, (0, 1, 2, 3, 5, 9): 1},
        (0, 1, 2, 5):{ (0, 1, 2, 3, 5): 2},
        (0, 1, 2, 8):{ (0, 1, 2, 3): 2},
        (0, 1, 2, 8, 9):{ (0, 1, 2, 3, 9): 2},
        (0, 1, 7):{ (0, 1, 2): 3, (0, 1, 7, 8): 2},
        (0, 1, 7, 8):{ (0, 1, 7, 8, 9): 3},
        (0, 1, 7, 8, 9):{ (0, 1, 2, 8, 9): 3},
        (0, 2):{ (0, 1, 2): 1, (0, 2, 5): 2},
        (0, 2, 5):{ (0, 2, 5, 8): 1},
        (0, 2, 5, 8):{ (0, 1, 2, 5): 1},
        (0, 7):{ (0, 1, 7): 1, (0, 2): 3}
        })

    if (G.is_isomorphic(GT) != True):
        return False
    else:
        print "B3 example 2 passed."

    # type B, fig. 8.3 pg. 174

    return True

def compare_graphs(g1, g2, node1, node2):
    r"""
    Compare two edge-labeled :class:`graphs <DiGraph>` obtained from
    ``Crystal.digraph()``, starting from the root nodes of each graph.

    - ``g1`` -- :class:`graphs <DiGraph>`, first digraph
    - ``g2`` -- :class:`graphs <DiGraph>`, second digraph
    - ``node1`` -- element of ``g1``
    - ``node2`` -- element of ``g2``

    Traverse ``g1`` starting at ``node1`` and compare this graph with
    the one obtained by traversing ``g2`` starting with ``node2``.
    If the graphs match (including labels) then return ``True``.
    Return ``False`` otherwise.

    EXAMPLES::

        sage: from sage.combinat.crystals.alcove_path import compare_graphs
        sage: G1 = crystals.Tableaux(['A',3], shape=[1,1]).digraph()
        sage: C = crystals.AlcovePaths(['A',3],[0,1,0])
        sage: G2 = C.digraph()
        sage: compare_graphs(G1, G2, C( () ), G2.vertices()[0])
        True
    """
    for out_edge in g1.outgoing_edges( node1 ):
        matched = False
        for o2 in g2.outgoing_edges( node2 ):
            if o2[2] == out_edge[2]:
                if matched == True:
                    print "ERROR:  Two edges with the same label for ", out_edge, " exist."
                    return False
                matched = True
                result = compare_graphs(g1, g2, out_edge[1], o2[1])
                if result == False:
                    return False
        if matched == False:
            print "ERROR:  No matching edge for ", out_edge, "."
            return False
    return True

def _test_against_tableaux(R, N, k, clss=CrystalOfAlcovePaths):
    r"""
    Test :class:`~sage.combinat.crystals.alcove_path.CrystalOfAlcovePaths`
    against all of the tableaux crystals of type `R` in rank `N` with
    highest weight given by a partition of `k`.

    EXAMPLES::

        sage: from sage.combinat.crystals.alcove_path import _test_against_tableaux
        sage: _test_against_tableaux(['A',3], 3, 2)
        ** Shape  [2]
          T has  10  nodes.
          C weight  [2, 0, 0]
          C has  10  nodes.
          Compare graphs:  True
        ** Shape  [1, 1]
          T has  6  nodes.
          C weight  [0, 1, 0]
          C has  6  nodes.
          Compare graphs:  True
    """
    from sage.combinat.partition import Partitions
    from sage.combinat.crystals.tensor_product import CrystalOfTableaux
    shapes = Partitions(k).list()
    for shape in shapes:
        print "** Shape ", shape
        T = CrystalOfTableaux(R, shape = shape)
        ct = len(T.list())
        print "  T has ", ct, " nodes."
        #T.digraph().show(edge_labels=True)
        H = T.digraph()
        weight = T.module_generators[0].weight()
        w = [ weight.scalar(RootSystem(R).ambient_space().simple_coroot(i)) for i in range(1,N+1) ]
        print "  C weight ", w

        C = clss(R , w)

        cc = len(C.list())
        #C.digraph().show(edge_labels=True)
        G = C.digraph()
        print "  C has ", cc, " nodes."
        if cc != ct:
            print "FAIL: number of nodes differ.", cc, ct
            return
        print "  Compare graphs: ", compare_graphs(G, H, C(()), H.vertices()[0])

def _test_with_lspaths_crystal(cartan_type, weight, depth=10):
    r"""
    Test if the digraphs generated are isomorphic to the ones generated by
    LS-path model.

    INPUT:

    - ``cartan_type`` -- Cartan type of a finite or affine untwisted root
      system
    - ``weight`` -- dominant weight as a list of (integral) coefficients of the
      fundamental weights
    - ``depth`` -- starting at the module generator how deep do you want to
      generate the crystal, useful for affine types

    EXAMPLES::

        sage: from sage.combinat.crystals.alcove_path import _test_with_lspaths_crystal
        sage: _test_with_lspaths_crystal(['A',3,1],[1,0,0,0],10) #long time
        True
        sage: _test_with_lspaths_crystal(['G',2,1],[1,0,0,0,0],10) #long time
        True
    """
    from sage.combinat.crystals.littelmann_path import CrystalOfLSPaths
    G1 = CrystalOfAlcovePaths(cartan_type, weight).digraph_fast(depth)
    C = CrystalOfLSPaths(cartan_type, weight)
    G2 = C.digraph(subset=C.subcrystal(max_depth=depth, direction='lower'))

    return G1.is_isomorphic(G2, edge_labels=True)

