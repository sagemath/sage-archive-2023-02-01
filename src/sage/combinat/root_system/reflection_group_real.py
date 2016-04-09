r"""
Finite real reflection groups
-------------------------------

AUTHORS:

- Christian Stump (initial version 2011--2015)

.. NOTE::

    - For definitions and classification types of finite complex reflection groups, see :wikipedia:`Complex_reflection_group`.
    - Uses the GAP3 package *chevie* available at `Jean Michel's website <http://webusers.imj-prg.fr/~jean.michel/gap3/>`_.

.. WARNING:: works only if the GAP3 package Chevie is available.

.. TODO::

    - Properly provide root systems for real reflection groups
    - Element class should be unique to be able to work with large groups without creating elements multiple times.
"""
#*****************************************************************************
#       Copyright (C) 2011-2015 Christian Stump <christian.stump at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from copy import copy
from sage.misc.all import prod
from sage.misc.cachefunc import cached_function, cached_method, cached_in_parent_method
from sage.categories.category import Category
from sage.categories.finite_permutation_groups import FinitePermutationGroups
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.combinat.root_system.cartan_type import CartanType
from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.rings.all import ZZ, QQ
from sage.matrix.all import Matrix, identity_matrix
from sage.matrix.matrix import is_Matrix
from sage.interfaces.gap3 import GAP3Record, gap3
from sage.interfaces.gap import gap
from sage.combinat.words.word import Word
from sage.rings.arith import gcd, lcm
from sage.combinat.root_system.reflection_group_complex import ComplexReflectionGroup, IrreducibleComplexReflectionGroup
from sage.categories.coxeter_groups import CoxeterGroups
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.root_system.coxeter_group import is_chevie_available

from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField, E

def ReflectionGroup(*args,**kwds):
    r"""
    Construct a finite (complex or real) reflection group as a Sage
    permutation group by fetching the permutation representation of the
    generators from chevie's database.

    INPUT:

    can be one or multiple of the following:

    - triple `(r,p,n)` with `p` divides `r`, which denotes the group
      `G(r,p,n)`

    - integer between `4` and `37`, which denotes an exceptional
      irreducible complex reflection group

    - finite Cartan-Killing type

    EXAMPLES:

    Finite reflection groups can be constructed from

    Cartan-Killing classification types::

        sage: W = ReflectionGroup(['A',3]); W                           # optional - chevie
         Irreducible real reflection group of rank 3 and type A3        # optional - chevie

        sage: W = ReflectionGroup(['H',4]); W                           # optional - chevie
         Irreducible real reflection group of rank 4 and type H4

        sage: W = ReflectionGroup(['I',5]); W                           # optional - chevie
         Irreducible real reflection group of rank 2 and type I2(5)

    the complex infinite family `G(r,p,n)` with `p` divides `r`::

        sage: W = ReflectionGroup((1,1,4)); W                           # optional - chevie
        Irreducible complex reflection group of rank 3 and type A3

        sage: W = ReflectionGroup((2,1,3)); W                           # optional - chevie
        Irreducible complex reflection group of rank 3 and type B3

    Chevalley-Shepard-Todd exceptional classification types::

        sage: W = ReflectionGroup(23); W                                # optional - chevie
         Irreducible complex reflection group of rank 3 and type H3

    multiples of the above::

        sage: W = ReflectionGroup(['A',2],['B',2]); W
        Reducible real reflection group of rank 4 and type A2 x B2

        sage: W = ReflectionGroup(['A',2],4); W
        Reducible complex reflection group of rank 4 and type A2 x ST4

        sage: W = ReflectionGroup((4,2,2),4); W
        Reducible complex reflection group of rank 4 and type G(4,2,2) x ST4
    """
    if not is_chevie_available():
        raise ImportError("The GAP3 package 'chevie' is needed to work with (complex) reflection groups")
    gap3.load_package("chevie")

    W_types     = []
    is_complex  = False
    for arg in args:
        # preparsing
        if type(arg) is list:
            X = tuple(arg)
        else:
            X = arg

        # precheck for valid input data
        if not ( is_Matrix(X) or isinstance(X,CartanMatrix) or isinstance(X,tuple) or ( X in ZZ and 4 <= X <= 37 ) ):
            raise ValueError("The input data (%s) is not valid for reflection groups."%X)

        # check for real vs complex
        elif X in ZZ or ( isinstance(X,tuple) and len(X) == 3 ):
            is_complex = True

        # transforming two reducible types
        if X == (2,2,2):
            W_types.extend([(1,1,2),(1,1,2)])
        elif X == ('I',2):
            W_types.extend([('A',1),('A',1)])
        else:
            W_types.append(X)

    for index_set_kwd in ['index_set','hyperplane_index_set','reflection_index_set']:
        index_set = kwds.get(index_set_kwd, None)
        if index_set is not None:
            from sage.sets.family import Family
            if type(index_set) in [list,tuple]:
                kwds[index_set_kwd] = Family(index_set, lambda x: index_set.index(x))
            elif type(index_set) is dict:
                kwds[index_set_kwd] = Family(index_set)
            else:
                raise ValueError('The keyword %s must be a list, tuple, or dict'%index_set_kwd)

    if len(W_types) == 1:
        if is_complex is True:
            cls = IrreducibleComplexReflectionGroup
        else:
            cls = IrreducibleRealReflectionGroup
    else:
        if is_complex is True:
            cls = ComplexReflectionGroup
        else:
            cls = RealReflectionGroup
    return cls(tuple(W_types), index_set=kwds.get('index_set', None),
                               hyperplane_index_set=kwds.get('hyperplane_index_set', None),
                               reflection_index_set=kwds.get('reflection_index_set', None) )

class RealReflectionGroup(ComplexReflectionGroup):
    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None):
        r"""

        TESTS::

            sage: W = ReflectionGroup(['A',3])
            sage: TestSuite(W).run()
        """
        W_types = tuple( tuple( W_type ) if isinstance(W_type,(list,tuple)) else W_type for W_type in W_types )
        cartan_types = []
        for W_type in W_types:
            W_type = CartanType(W_type)
            if not W_type.is_finite() or not W_type.is_irreducible():
                raise ValueError("The given Cartan type of a component is not irreducible and finite")
            cartan_types.append( W_type )
        if len(W_types) == 1:
            cls = IrreducibleComplexReflectionGroup
        else:
            cls = ComplexReflectionGroup
        cls.__init__(self, W_types, index_set               = index_set,
                                    hyperplane_index_set    = hyperplane_index_set,
                                    reflection_index_set    = reflection_index_set)
        N = self.nr_reflections()
        self._is_positive_root = [None] + [True]*N + [False]*N

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2],['I',5],['I',6]); W
            Reducible real reflection group of rank 9 and type A3 x B2 x I2(5) x G2
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible real reflection group of rank %s and type %s'%(self._rank,type_str)

    def __iter__(self):
        from sage.combinat.root_system.reflection_group_c import Iterator
        return iter(Iterator(self))

    def _iterator_tracking_words(self):
        r"""
        Return an iterator through the elements of ``self`` together
        with the words in the simple generators.

        The iterator is a breadth first search through the right weak
        order of ``self``.

        .. REMARK::

            In order to save space, the fact that the right weak order
            is graded is used.

        .. TODO::

            This algorithm could be still much optimized:

            - the right weak order is self-dual under the action of the
              longest element, so one only needs to search through the
              first half.

            - the coset decomposition used in chevie is much quicker.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])
            sage: for w in W._iterator_tracking_words(): print w
            ((), ())
            ((1,4)(2,3)(5,6), (0,))
            ((1,3)(2,5)(4,6), (1,))
            ((1,6,2)(3,5,4), (0, 1))
            ((1,2,6)(3,4,5), (1, 0))
            ((1,5)(2,4)(3,6), (0, 1, 0))
        """
        I = self.gens()
        index_list = range(len(I))

        level_set_old   = set()
        level_set_cur   = [ (self.one(), tuple()) ]
        while level_set_cur:
            level_set_new = []
            for x, word in level_set_cur:
                yield x, word
                for i in index_list:
                    y = x._mul_(I[i])
                    if y not in level_set_old:
                        level_set_old.add(y)
                        level_set_new.append((y, word+(i,)))
            level_set_old = set( elt[0] for elt in level_set_cur )
            level_set_cur = level_set_new

    @cached_method
    def bipartite_index_set(self):
        r"""
        Return the bipartite index set of a real reflection group.

        EXAMPLES::

            sage: W = ReflectionGroup(["A",5])
            sage: W.bipartite_index_set()
            [[0, 2, 4], [1, 3]]

            sage: W = ReflectionGroup(["A",5],index_set=['a','b','c','d','e'])
            sage: W.bipartite_index_set()
            [['a', 'c', 'e'], ['b', 'd']]
        """
        index_family = self._index_set
        keys = index_family.keys()
        L,R = self._gap_group.BipartiteDecomposition().sage()
        L = [ i for i in keys if index_family[i]+1 in L ]
        R = [ i for i in keys if index_family[i]+1 in R ]
        return [L,R]

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2])
            sage: W.cartan_type()
            [['A', 3], ['B', 2]]

            sage: W = ReflectionGroup(['A',3])
            sage: W.cartan_type()
            ['A', 3]
        """
        if len(self._type) == 1:
            ct = self._type[0]
            return CartanType([ct['series'],ct['rank']])
        else:
            return [ W.cartan_type() for W in self.irreducible_components() ]

    def invariant_form(self):
        r"""
        Returns the form that is invariant under the action of ``self``.

        EXAMPLES::

            sage: tba
        """
        C       = self.cartan_matrix()
        n       = self.rank()

        if self.is_crystallographic():
            ring = QQ
        else:
            ring = UniversalCyclotomicField()

        from sage.matrix.constructor import zero_matrix
        form = zero_matrix(ring,n,n)

        for j in range(n):
            for i in range(j):
                if C[i,j] != 0:
                    form[j,j] = form[i,i]*C[i,j]/C[j,i]
            if form[j,j] == 0:
                form[j,j] = 1
        for j in range(n):
            for i in range(j):
                form[i,j] = C[i,j] * form[i,i] / 2
                form[j,i] = C[j,i] * form[j,j] / 2

        form.set_immutable()
        return form

    def simple_root(self,i):
        r"""
        Return the simple root with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])
            sage: W.simple_root(0)
            (1, 0, 0)
        """
        return self.simple_roots()[self._index_set[i]]

    def positive_roots(self):
        r"""
        Return the simple root with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2])
            sage: W.positive_roots()
                [(1, 0, 0, 0, 0),
                 (0, 1, 0, 0, 0),
                 (0, 0, 1, 0, 0),
                 (0, 0, 0, 1, 0),
                 (0, 0, 0, 0, 1),
                 (1, 1, 0, 0, 0),
                 (0, 1, 1, 0, 0),
                 (0, 0, 0, 1, 1),
                 (1, 1, 1, 0, 0),
                 (0, 0, 0, 2, 1)]

            sage: W = ReflectionGroup(['A',3])
            sage: W.positive_roots()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 1, 1)]
        """
        return self.roots()[:self.nr_reflections()]

    def almost_positive_roots(self):
        r"""
        Return the almost positive roots of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2])
            sage: W.almost_positive_roots()
                [(-1, 0, 0, 0, 0),
                 (0, -1, 0, 0, 0),
                 (0, 0, -1, 0, 0),
                 (0, 0, 0, -1, 0),
                 (0, 0, 0, 0, -1),
                 (1, 0, 0, 0, 0),
                 (0, 1, 0, 0, 0),
                 (0, 0, 1, 0, 0),
                 (0, 0, 0, 1, 0),
                 (0, 0, 0, 0, 1),
                 (1, 1, 0, 0, 0),
                 (0, 1, 1, 0, 0),
                 (0, 0, 0, 1, 1),
                 (1, 1, 1, 0, 0),
                 (0, 0, 0, 2, 1)]


            sage: W = ReflectionGroup(['A',3])
            sage: W.almost_positive_roots()
                [(-1, 0, 0),
                 (0, -1, 0),
                 (0, 0, -1),
                 (1, 0, 0),
                 (0, 1, 0),
                 (0, 0, 1),
                 (1, 1, 0),
                 (0, 1, 1),
                 (1, 1, 1)]
        """
        return [ -beta for beta in self.simple_roots() ] + self.positive_roots()

    def root_to_reflection(self,root):
        r"""
        Return the reflection along the given ``root``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])
            sage: for beta in W.roots(): print W.root_to_reflection(beta)
            (1,4)(2,3)(5,6)
            (1,3)(2,5)(4,6)
            (1,5)(2,4)(3,6)
            (1,4)(2,3)(5,6)
            (1,3)(2,5)(4,6)
            (1,5)(2,4)(3,6)
        """
        Phi = self.roots()
        R = self.reflections()
        i = Phi.index(root)+1
        j = Phi.index(-root)+1
        for r in R:
            if r(i) == j:
                return r
        raise ValueError("There is a bug in root_to_reflection!")

    def reflection_to_positive_root(self,r):
        r"""
        Return the positive root orthogonal to the given reflection.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])
            sage: for r in W.reflections(): print W.reflection_to_positive_root(r)
            (1, 0)
            (0, 1)
            (1, 1)
        """
        Phi = self.roots()
        N = len(Phi)/2
        for i in range(1,N+1):
            if r(i) == i+N:
                return Phi[i-1]
        raise ValueError("There is a bug in reflection_to_positive_root!")

    @cached_method
    def fundamental_weights(self):
        r"""
        Return the fundamental weights of ``self`` in terms of the simple roots.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2])
            sage: W.fundamental_weights()
            [(3/4, 1/2, 1/4, 0, 0),
             (1/2, 1, 1/2, 0, 0),
             (1/4, 1/2, 3/4, 0, 0),
             (0, 0, 0, 1, 1/2),
             (0, 0, 0, 1, 1)]

            sage: W = ReflectionGroup(['A',3])
            sage: W.fundamental_weights()
            [(3/4, 1/2, 1/4), (1/2, 1, 1/2), (1/4, 1/2, 3/4)]
        """
        m = self.cartan_matrix().transpose().inverse()
        S = self.simple_roots()
        zero = S[0] - S[0]
        weights = [ sum( [ m[i,j] * S[j] for j in range(len(S)) ], zero ) for i in range(len(S)) ]
        for weight in weights:
            weight.set_immutable()
        return weights

    def fundamental_weight(self,i):
        r"""
        Return the fundamental weight with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])
            sage: [ W.fundamental_weight(i) for i in W.index_set() ]
            [(3/4, 1/2, 1/4), (1/2, 1, 1/2), (1/4, 1/2, 3/4)]
        """
        return self.fundamental_weights()[self._index_set[i]]

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        .. TODO::

            Move this method to the CoxeterGroups category. The issue
            with this is that the indexing of a Coxeter group is not
            handled in the category, so that ``self._index_set`` is
            not required to do what is expected here.

        EXAMPLES::

            sage: G = ReflectionGroup(['A',3])
            sage: G.coxeter_matrix()
            [1 3 2]
            [3 1 3]
            [2 3 1]
        """
        from sage.rings.integer_ring import ZZ
        from sage.matrix.all import MatrixSpace

        S = self.simple_reflections()
        I = self.index_set()
        I_inv = self._index_set
        n = self.rank()
        MS = MatrixSpace(ZZ, n)
        m = MS(0)
        for i in I:
            for j in I:
                m[I_inv[i],I_inv[j]] = (S[i]*S[j]).order()
        return m

    def permutahedron(self,point=None):
        r"""
        Return the permutahedron of ``self``.

        This is the convex hull of the point ``point`` in the weight
        basis under the action of ``self`` on the underlying vector
        space `V`.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2])
            sage: W.permutahedron()
            A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 192 vertices

            sage: W = ReflectionGroup(['A',3])
            sage: W.permutahedron()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices
        """
        n = self.rank()
        weights = self.fundamental_weights()
        if point is None:
            point = [1]*n
        v = sum( point[i] * weights[i] for i in range(n) )
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron( vertices=[ v*(~w).as_matrix() for w in self] )

    @cached_method
    def right_coset_representatives(self,J):
        r"""
        Return the right coset representatives of ``self`` for the
        parabolic subgroup generated by the simple reflections in ``J``.

        EXAMPLES::

            sage: W = ReflectionGroup(["A",3])
            sage: for J in Subsets([0,1,2]): print W.right_coset_representatives(J)
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4)(2,8)(3,5)(7,10)(9,11), (1,7)(2,4)(5,6)(8,10)(11,12), (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,4,6)(2,3,11)(5,8,9)(7,10,12), (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,7)(2,6)(3,9)(4,5)(8,12)(10,11), (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,2,3,12)(4,5,10,11)(6,7,8,9), (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,6)(2,9)(3,8)(5,11)(7,12), (1,8)(2,7)(3,6)(4,10)(9,12), (1,10,9,5)(2,6,8,12)(3,11,7,4), (1,12,3,2)(4,11,10,5)(6,9,8,7), (1,3)(2,12)(4,10)(5,11)(6,8)(7,9), (1,5,12)(2,9,4)(3,10,8)(6,7,11), (1,8,11)(2,5,7)(3,12,4)(6,10,9), (1,11,8)(2,7,5)(3,4,12)(6,9,10), (1,12,5)(2,4,9)(3,8,10)(6,11,7), (1,3,7,9)(2,11,6,10)(4,8,5,12), (1,9,7,3)(2,10,6,11)(4,12,5,8), (1,11)(3,10)(4,9)(5,7)(6,12), (1,9)(2,8)(3,7)(4,11)(5,10)(6,12)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4)(2,8)(3,5)(7,10)(9,11), (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,4,6)(2,3,11)(5,8,9)(7,10,12), (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,2,3,12)(4,5,10,11)(6,7,8,9), (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,6)(2,9)(3,8)(5,11)(7,12), (1,3)(2,12)(4,10)(5,11)(6,8)(7,9), (1,5,12)(2,9,4)(3,10,8)(6,7,11), (1,3,7,9)(2,11,6,10)(4,8,5,12)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,7)(2,4)(5,6)(8,10)(11,12), (1,4,6)(2,3,11)(5,8,9)(7,10,12), (1,7)(2,6)(3,9)(4,5)(8,12)(10,11), (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,2,3,12)(4,5,10,11)(6,7,8,9), (1,10,9,5)(2,6,8,12)(3,11,7,4), (1,12,3,2)(4,11,10,5)(6,9,8,7), (1,8,11)(2,5,7)(3,12,4)(6,10,9), (1,12,5)(2,4,9)(3,8,10)(6,11,7), (1,11)(3,10)(4,9)(5,7)(6,12)]
            [(), (1,4)(2,8)(3,5)(7,10)(9,11), (1,7)(2,4)(5,6)(8,10)(11,12), (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,8)(2,7)(3,6)(4,10)(9,12), (1,12,3,2)(4,11,10,5)(6,9,8,7), (1,3)(2,12)(4,10)(5,11)(6,8)(7,9), (1,11,8)(2,7,5)(3,4,12)(6,9,10), (1,9,7,3)(2,10,6,11)(4,12,5,8)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4,6)(2,3,11)(5,8,9)(7,10,12), (1,2,3,12)(4,5,10,11)(6,7,8,9)]
            [(), (1,4)(2,8)(3,5)(7,10)(9,11), (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,3)(2,12)(4,10)(5,11)(6,8)(7,9)]
            [(), (1,7)(2,4)(5,6)(8,10)(11,12), (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,12,3,2)(4,11,10,5)(6,9,8,7)]
            [()]
        """
        from sage.combinat.root_system.reflection_group_complex import _gap_return
        J_inv = [ self._index_set[j]+1 for j in J ]
        S = str(gap3('ReducedRightCosetRepresentatives(%s,ReflectionSubgroup(%s,%s))'%(self._gap_group._name,self._gap_group._name,J_inv)))
        exec('L = ' + _gap_return(S))
        return L

    class Element(ComplexReflectionGroup.Element):

        def _compute_reduced_word(self):
            r"""
            Computes a reduced word and stores it into ``self._reduced_word``.

            TESTS::

                sage: W = ReflectionGroup(['A',2])
                sage: for w in W: w._reduced_word = None; w._compute_reduced_word()
                sage: [ w._reduced_word for w in W ]
                [[], [0], [1], [0, 1], [1, 0], [0, 1, 0]]
            """
            self._reduced_word = CoxeterGroups.ElementMethods.reduced_word.__func__(self)

        def length(self):
            r"""
            Return the length of ``self`` in generating reflections.

            This is the minimal numbers of generating reflections needed
            to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])
                sage: for w in W:
                ....:   print w.reduced_word(), w.length()
                 0
                0 1
                1 1
                01 2
                10 2
                010 3
            """
            if not self._reduced_word is None:
                return len(self._reduced_word)
            else:
                N = self.parent().nr_reflections()
                return sum( 1 for i in range(N) if not self.parent()._is_positive_root[self(i+1)] )

        def has_left_descent(self, i):
            r"""
            Return whether ``i`` is a left descent of ``self``.

            This is done by testing whether ``i`` is mapped by ``self``
            to a negative root.

            EXAMPLES::

                sage: W = ReflectionGroup(["A",3])
                sage: s = W.simple_reflections()
                sage: (s[1]*s[2]).has_left_descent(1)
                True
                sage: (s[1]*s[2]).has_left_descent(2)
                False
            """
            W = self.parent()
            return not W._is_positive_root[self(W._index_set[i]+1)]

        def has_descent(self, i, side='left', positive=False):
            r"""
            Return whether ``i`` is a descent (or ascent) of ``self``.

            This is done by testing whether ``i`` is mapped by ``self``
            to a negative root.

            INPUT:

            - ``i`` - an index of a simple reflection
            - ``side`` (default: 'right') - 'left' or 'right'
            - ``positive`` (default: ``False``) - a boolean

            EXAMPLES::

                sage: W = ReflectionGroup(["A",3])
                sage: s = W.simple_reflections()
                sage: (s[1]*s[2]).has_descent(1)
                True
                sage: (s[1]*s[2]).has_descent(2)
                False
            """
            if not isinstance(positive, bool):
                raise TypeError("%s is not a boolean"%(bool))

            if i not in self.parent().hyperplane_index_set():
                raise ValueError("The given index %s is not in the index set"%i)

            negative = not positive

            if side == 'left':
                return self.has_left_descent(i) is negative
            elif side == 'right':
                return self.has_right_descent(i) is negative
            else:
                raise ValueError("The method 'has_descent' needs the input 'side' to be either 'left' or 'right'.")

        def act_on_root(self,root):
            r"""
            Return the root obtained by applying ``self`` on ``root``.

            EXAMPLES::

                sage: W = CoxeterGroup(['A',2],implementation='chevie')
                sage: for w in W:
                ....:     print w.reduced_word(), [ w.act_on_root(beta) for beta in W.roots() ]
                     [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]
                0 [(-1, 0), (1, 1), (0, 1), (1, 0), (-1, -1), (0, -1)]
                1 [(1, 1), (0, -1), (1, 0), (-1, -1), (0, 1), (-1, 0)]
                01 [(0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1), (1, 0)]
                10 [(-1, -1), (1, 0), (0, -1), (1, 1), (-1, 0), (0, 1)]
                010 [(0, -1), (-1, 0), (-1, -1), (0, 1), (1, 0), (1, 1)]
            """
            Phi = self.parent().roots()
            return Phi[ (~self)(Phi.index(root)+1)-1 ]

        def inversion_set(self):
            r"""
            Return the inversion set of ``self``.

            This is the set `\{ \beta \in \Phi^+ : ``self``(\beta) \in \Phi^- \}`.

            EXAMPLES::

                sage: W = CoxeterGroup(['A',2],implementation='chevie')
                sage: for w in W:
                ....:     print w.reduced_word(), w.inversion_set()
                 []
                0 [(1, 0)]
                1 [(0, 1)]
                01 [(0, 1), (1, 1)]
                10 [(1, 0), (1, 1)]
                010 [(0, 1), (1, 0), (1, 1)]
            """
            Phi_plus = set(self.parent().positive_roots())
            return [ root for root in Phi_plus if self.act_on_root(root) not in Phi_plus ]

        @cached_in_parent_method
        def right_coset_representatives(self):
            r"""
            Return the right coset representatives of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])
                sage: for w in W: print w.reduced_word(), [ v.reduced_word() for v in w.right_coset_representatives() ]
                 [word: , word: 1, word: 0, word: 10, word: 01, word: 010]
                0 [word: , word: 1, word: 10]
                1 [word: , word: 0, word: 01]
                01 [word: ]
                10 [word: ]
                010 [word: , word: 1, word: 0]
            """
            from sage.combinat.root_system.reflection_group_complex import _gap_return
            W = self.parent()
            T = W.reflections()
            T_fix = [ i+1 for i in T.keys() if self.fix_space().is_subspace(T[i].fix_space()) ]
            S = str(gap3('ReducedRightCosetRepresentatives(%s,ReflectionSubgroup(%s,%s))'%(W._gap_group._name,W._gap_group._name,T_fix)))
            exec('L = ' + _gap_return(S,coerce_obj='W'))
            return L

        def left_coset_representatives(self):
            r"""
            Return the left coset representatives of ``self``.

            .. SEEALSO:: :meth:`right_coset_representatives`

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])
                sage: for w in W: print w.reduced_word(), [ v.reduced_word() for v in w.left_coset_representatives() ]
                 [word: , word: 1, word: 0, word: 01, word: 10, word: 010]
                0 [word: , word: 1, word: 01]
                1 [word: , word: 0, word: 10]
                01 [word: ]
                10 [word: ]
                010 [word: , word: 1, word: 0]
            """
            return [ (~w) for w in self.right_coset_representatives() ]

class IrreducibleRealReflectionGroup(RealReflectionGroup, IrreducibleComplexReflectionGroup):

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: for i in [2..7]: print ReflectionGroup(["I",i])
            Reducible real reflection group of rank 2 and type A1 x A1
            Irreducible real reflection group of rank 2 and type A2
            Irreducible real reflection group of rank 2 and type B2
            Irreducible real reflection group of rank 2 and type I2(5)
            Irreducible real reflection group of rank 2 and type G2
            Irreducible real reflection group of rank 2 and type I2(7)
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible real reflection group of rank %s and type %s'%(self._rank,type_str)

    class Element(RealReflectionGroup.Element,IrreducibleComplexReflectionGroup.Element):
        pass

