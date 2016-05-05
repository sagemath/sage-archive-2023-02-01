r"""
Finite real reflection groups
-------------------------------

Let `V` be a finite-dimensional real vector space. A reflection of
`V` is an operator `r \in \operatorname{GL}(V)` that has order `2`
and fixes pointwise a hyperplane in `V`.
In the present implementation, finite real reflection groups are
tied with a root system.

Finite real reflection groups with root systems have been classified
according to finite Cartan-Killing types.
For more definitions and classification types of finite complex
reflection groups, see :wikipedia:`Complex_reflection_group`.

The point of entry to work with reflection groups is :func:`~sage.combinat.root_system.reflection_group_real.ReflectionGroup`
which can be used with finite Cartan-Killing types::

    sage: ReflectionGroup(['A',2])                                      # optional - gap3
    Irreducible real reflection group of rank 2 and type A2
    sage: ReflectionGroup(['F',4])                                      # optional - gap3
    Irreducible real reflection group of rank 4 and type F4
    sage: ReflectionGroup(['H',3])                                      # optional - gap3
    Irreducible real reflection group of rank 3 and type H3

AUTHORS:

- Christian Stump (initial version 2011--2015)

.. WARNING::

    Uses the GAP3 package *Chevie* which is available as an
    experimental package (installed by ``sage -i gap3``) or to
    download by hand from `Jean Michel's website
    <http://webusers.imj-prg.fr/~jean.michel/gap3/>`_.

.. TODO::

    - Implement descents, left/right descents, ``has_descent``,
      ``first_descent`` directly in this class, since the generic
      implementation is much slower.
"""
#*****************************************************************************
#       Copyright (C) 2011-2016 Christian Stump <christian.stump at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.root_system.cartan_type import CartanType, CartanType_abstract
from sage.rings.all import ZZ, QQ
from sage.matrix.matrix import is_Matrix
from sage.interfaces.gap3 import gap3
from sage.combinat.root_system.reflection_group_complex import ComplexReflectionGroup, IrreducibleComplexReflectionGroup
from sage.categories.coxeter_groups import CoxeterGroups
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.root_system.coxeter_group import is_chevie_available
from sage.misc.sage_eval import sage_eval
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
from sage.combinat.root_system.reflection_group_c import reduced_word_c

def ReflectionGroup(*args,**kwds):
    r"""
    Construct a finite (complex or real) reflection group as a Sage
    permutation group by fetching the permutation representation of the
    generators from chevie's database.

    INPUT:

    can be one or multiple of the following:

    - a triple `(r, p, n)` with `p` divides `r`, which denotes the group
      `G(r, p, n)`

    - an integer between `4` and `37`, which denotes an exceptional
      irreducible complex reflection group

    - a finite Cartan-Killing type

    EXAMPLES:

    Finite reflection groups can be constructed from

    Cartan-Killing classification types::

        sage: W = ReflectionGroup(['A',3]); W                           # optional - gap3
         Irreducible real reflection group of rank 3 and type A3

        sage: W = ReflectionGroup(['H',4]); W                           # optional - gap3
         Irreducible real reflection group of rank 4 and type H4

        sage: W = ReflectionGroup(['I',5]); W                           # optional - gap3
         Irreducible real reflection group of rank 2 and type I2(5)

    the complex infinite family `G(r,p,n)` with `p` divides `r`::

        sage: W = ReflectionGroup((1,1,4)); W                           # optional - gap3
        Irreducible real reflection group of rank 3 and type A3

        sage: W = ReflectionGroup((2,1,3)); W                           # optional - gap3
        Irreducible real reflection group of rank 3 and type B3

    Chevalley-Shepard-Todd exceptional classification types::

        sage: W = ReflectionGroup(23); W                                # optional - gap3
        Irreducible real reflection group of rank 3 and type H3

    Cartan types and matrices::

        sage: ReflectionGroup(CartanType(['A',2]))                      # optional - gap3
        Irreducible real reflection group of rank 2 and type A2

        sage: ReflectionGroup(CartanType((['A',2],['A',2])))            # optional - gap3
        Reducible real reflection group of rank 4 and type A2 x A2

        sage: C = CartanMatrix(['A',2])                                 # optional - gap3
        sage: ReflectionGroup(C)                                        # optional - gap3
        Irreducible real reflection group of rank 2 and type A2

    multiples of the above::

        sage: W = ReflectionGroup(['A',2],['B',2]); W                   # optional - gap3
        Reducible real reflection group of rank 4 and type A2 x B2

        sage: W = ReflectionGroup(['A',2],4); W                         # optional - gap3
        Reducible complex reflection group of rank 4 and type A2 x ST4

        sage: W = ReflectionGroup((4,2,2),4); W                         # optional - gap3
        Reducible complex reflection group of rank 4 and type G(4,2,2) x ST4
    """
    if not is_chevie_available():
        raise ImportError("the GAP3 package 'chevie' is needed to work with (complex) reflection groups")
    gap3.load_package("chevie")

    error_msg = "the input data (%s) is not valid for reflection groups"

    W_types = []
    is_complex = False
    for arg in args:
        # preparsing
        if isinstance(arg, list):
            X = tuple(arg)
        else:
            X = arg

        # precheck for valid input data
        if not (isinstance(X, (CartanType_abstract,tuple)) or (X in ZZ and 4 <= X <= 37)):
            raise ValueError(error_msg%X)

        # transforming two reducible types and an irreducible type
        if isinstance(X, CartanType_abstract):
            if not X.is_finite():
                raise ValueError(error_msg%X)
            if hasattr(X,"cartan_type"):
                X = X.cartan_type()
            if X.is_irreducible():
                W_types.extend([(X.letter, X.n)])
            else:
                W_types.extend([(x.letter, x.n) for x in X.component_types()])

        elif X == (2,2,2) or X == ('I',2):
            W_types.extend([('A',1), ('A',1)])

        elif X == (2,2,3):
            W_types.extend([('A', 3)])

        else:
            W_types.append(X)

    # converting the real types given as complex types
    # and then checking for real vs complex
    for i,W_type in enumerate(W_types):
        if W_type in ZZ:
            if W_type == 23:
                W_types[i] = ('H', 3)
            elif W_type == 28:
                W_types[i] = ('F', 4)
            elif W_type == 30:
                W_types[i] = ('H', 4)
            elif W_type == 35:
                W_types[i] = ('E', 6)
            elif W_type == 36:
                W_types[i] = ('E', 7)
            elif W_type == 37:
                W_types[i] = ('E', 8)
        if isinstance(W_type,tuple) and len(W_type) == 3:
            if W_type[0] == W_type[1] == 1:
                W_types[i] = ('A', W_type[2]-1)
            elif W_type[0] == 2 and W_type[1] == 1:
                W_types[i] = ('B', W_type[2])
            elif W_type[0] == W_type[1] == 2:
                W_types[i] = ('D', W_type[2])
            elif W_type[0] == W_type[1] and W_type[2] == 2:
                W_types[i] = ('I', W_type[0])

        W_type = W_types[i]
        # check for real vs complex
        if W_type in ZZ or (isinstance(W_type, tuple) and len(W_type) == 3):
            is_complex = True

    for index_set_kwd in ['index_set', 'hyperplane_index_set', 'reflection_index_set']:
        index_set = kwds.get(index_set_kwd, None)
        if index_set is not None:
            if isinstance(index_set, (list, tuple)):
                kwds[index_set_kwd] = tuple(index_set)
            else:
                raise ValueError('the keyword %s must be a list or tuple'%index_set_kwd)

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
    return cls(tuple(W_types),
               index_set=kwds.get('index_set', None),
               hyperplane_index_set=kwds.get('hyperplane_index_set', None),
               reflection_index_set=kwds.get('reflection_index_set', None))

class RealReflectionGroup(ComplexReflectionGroup):
    """
    A real reflection group given as a permutation group.

    .. SEEALSO::

        :func:`ReflectionGroup`
    """
    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: TestSuite(W).run()                                    # optional - gap3
        """
        W_types = tuple([tuple(W_type) if isinstance(W_type, (list,tuple)) else W_type
                         for W_type in W_types])
        cartan_types = []
        for W_type in W_types:
            W_type = CartanType(W_type)
            if not W_type.is_finite() or not W_type.is_irreducible():
                raise ValueError("the given Cartan type of a component is not irreducible and finite")
            cartan_types.append( W_type )
        if len(W_types) == 1:
            cls = IrreducibleComplexReflectionGroup
        else:
            cls = ComplexReflectionGroup
        cls.__init__(self, W_types, index_set            = index_set,
                                    hyperplane_index_set = hyperplane_index_set,
                                    reflection_index_set = reflection_index_set)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2],['I',5],['I',6])  # optional - gap3
            sage: W._repr_()                                            # optional - gap3
            'Reducible real reflection group of rank 9 and type A3 x B2 x I2(5) x G2'
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible real reflection group of rank %s and type %s'%(self._rank,type_str)

    def iteration(self, algorithm="breadth", tracking_words=True):
        r"""
        Return an iterator going through all elements in ``self``.

        INPUT:

        - ``algorithm`` (default: ``'breadth'``) -- must be one of
          the following:

          * ``'breadth'`` - iterate over in a linear extension of the
            weak order
          * ``'depth'`` - iterate by a depth-first-search

        - ``tracking_words`` (default: ``True``) -- whether or not to keep
          track of the reduced words and store them in ``_reduced_word``

        .. NOTE::

            The fastest iteration is the depth first algorithm without
            tracking words. In particular, ``'depth'`` is ~1.5x faster.

        EXAMPLES::

            sage: W = ReflectionGroup(["B",2])                          # optional - gap3

            sage: for w in W.iteration("breadth",True):                 # optional - gap3
            ....:     print("%s %s"%(w, w._reduced_word))               # optional - gap3
            () []
            (1,3)(2,6)(5,7) [1]
            (1,5)(2,4)(6,8) [0]
            (1,7,5,3)(2,4,6,8) [0, 1]
            (1,3,5,7)(2,8,6,4) [1, 0]
            (2,8)(3,7)(4,6) [1, 0, 1]
            (1,7)(3,5)(4,8) [0, 1, 0]
            (1,5)(2,6)(3,7)(4,8) [0, 1, 0, 1]

            sage: for w in W.iteration("depth", False): w               # optional - gap3
            ()
            (1,3)(2,6)(5,7)
            (1,5)(2,4)(6,8)
            (1,3,5,7)(2,8,6,4)
            (1,7)(3,5)(4,8)
            (1,7,5,3)(2,4,6,8)
            (2,8)(3,7)(4,6)
            (1,5)(2,6)(3,7)(4,8)
        """
        from sage.combinat.root_system.reflection_group_c import Iterator
        return iter(Iterator(self, N=self._number_of_reflections,
                             algorithm=algorithm, tracking_words=tracking_words))

    def __iter__(self):
        r"""
        Return an iterator going through all elements in ``self``.

        For options and faster iteration see :meth:`iteration`.

        EXAMPLES::

            sage: W = ReflectionGroup(["B",2])                          # optional - gap3

            sage: for w in W: print("%s %s"%(w, w._reduced_word))       # optional - gap3
            () []
            (1,3)(2,6)(5,7) [1]
            (1,5)(2,4)(6,8) [0]
            (1,7,5,3)(2,4,6,8) [0, 1]
            (1,3,5,7)(2,8,6,4) [1, 0]
            (2,8)(3,7)(4,6) [1, 0, 1]
            (1,7)(3,5)(4,8) [0, 1, 0]
            (1,5)(2,6)(3,7)(4,8) [0, 1, 0, 1]
        """
        return self.iteration(algorithm="breadth", tracking_words=True)

    @cached_method
    def bipartite_index_set(self):
        r"""
        Return the bipartite index set of a real reflection group.

        EXAMPLES::

            sage: W = ReflectionGroup(["A",5])                          # optional - gap3
            sage: W.bipartite_index_set()                               # optional - gap3
            [[1, 3, 5], [2, 4]]

            sage: W = ReflectionGroup(["A",5],index_set=['a','b','c','d','e'])  # optional - gap3
            sage: W.bipartite_index_set()                               # optional - gap3
            [['a', 'c', 'e'], ['b', 'd']]
        """
        L, R = self._gap_group.BipartiteDecomposition().sage()
        inv = self._index_set_inverse
        L = [i for i in self._index_set if inv[i] + 1 in L]
        R = [i for i in self._index_set if inv[i] + 1 in R]
        return [L, R]

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.cartan_type()                                       # optional - gap3
            ['A', 3]

            sage: W = ReflectionGroup(['A',3], ['B',2])                 # optional - gap3
            sage: W.cartan_type()                                       # optional - gap3
            A3xB2
        """
        if len(self._type) == 1:
            ct = self._type[0]
            return CartanType([ct['series'], ct['rank']])
        else:
            return CartanType([W.cartan_type() for W in self.irreducible_components()])

    def simple_root(self, i):
        r"""
        Return the simple root with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.simple_root(1)                                      # optional - gap3
            (1, 0, 0)
        """
        return self.simple_roots()[i]

    def positive_roots(self):
        r"""
        Return the positive roots of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3], ['B',2])                 # optional - gap3
            sage: W.positive_roots()                                    # optional - gap3
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

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.positive_roots()                                    # optional - gap3
            [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 1, 1)]
        """
        return self.roots()[:self._number_of_reflections]

    def almost_positive_roots(self):
        r"""
        Return the almost positive roots of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3], ['B',2])                 # optional - gap3
            sage: W.almost_positive_roots()                             # optional - gap3
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

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.almost_positive_roots()                             # optional - gap3
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
        return [-beta for beta in self.simple_roots()] + self.positive_roots()

    def root_to_reflection(self, root):
        r"""
        Return the reflection along the given ``root``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: for beta in W.roots(): W.root_to_reflection(beta)     # optional - gap3
            (1,4)(2,3)(5,6)
            (1,3)(2,5)(4,6)
            (1,5)(2,4)(3,6)
            (1,4)(2,3)(5,6)
            (1,3)(2,5)(4,6)
            (1,5)(2,4)(3,6)
        """
        Phi = self.roots()
        R = self.reflections()
        i = Phi.index(root) + 1
        j = Phi.index(-root) + 1
        for r in R:
            if r(i) == j:
                return r
        raise AssertionError("there is a bug in root_to_reflection")

    def reflection_to_positive_root(self, r):
        r"""
        Return the positive root orthogonal to the given reflection.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: for r in W.reflections(): print W.reflection_to_positive_root(r)  # optional - gap3
            (1, 0)
            (0, 1)
            (1, 1)
        """
        Phi = self.roots()
        N = len(Phi) / 2
        for i in range(1, N+1):
            if r(i) == i + N:
                return Phi[i-1]
        raise AssertionError("there is a bug in reflection_to_positive_root")

    @cached_method
    def fundamental_weights(self):
        r"""
        Return the fundamental weights of ``self`` in terms of the simple roots.

        The fundamental weights are defined by
        `s_j(\omega_i) = \omega_i - \delta_{i=j}\alpha_j`
        for the simple reflection `s_j` with corresponding simple
        roots `\alpha_j`.

        In other words, the transpose Cartan matrix sends the weight
        basis to the root basis. Observe again that the action here is
        defined as a right action, see the example below.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3], ['B',2])                 # optional - gap3
            sage: W.fundamental_weights()                               # optional - gap3
            Finite family {1: (3/4, 1/2, 1/4, 0, 0), 2: (1/2, 1, 1/2, 0, 0), 3: (1/4, 1/2, 3/4, 0, 0), 4: (0, 0, 0, 1, 1/2), 5: (0, 0, 0, 1, 1)}

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.fundamental_weights()                               # optional - gap3
            Finite family {1: (3/4, 1/2, 1/4), 2: (1/2, 1, 1/2), 3: (1/4, 1/2, 3/4)}

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: S = W.simple_reflections()                            # optional - gap3
            sage: N = W.fundamental_weights()                           # optional - gap3
            sage: for i in W.index_set():                               # optional - gap3
            ....:     for j in W.index_set():                           # optional - gap3
            ....:         print i, j, N[i], N[i]*S[j].to_matrix()       # optional - gap3
            1 1 (3/4, 1/2, 1/4) (-1/4, 1/2, 1/4)
            1 2 (3/4, 1/2, 1/4) (3/4, 1/2, 1/4)
            1 3 (3/4, 1/2, 1/4) (3/4, 1/2, 1/4)
            2 1 (1/2, 1, 1/2) (1/2, 1, 1/2)
            2 2 (1/2, 1, 1/2) (1/2, 0, 1/2)
            2 3 (1/2, 1, 1/2) (1/2, 1, 1/2)
            3 1 (1/4, 1/2, 3/4) (1/4, 1/2, 3/4)
            3 2 (1/4, 1/2, 3/4) (1/4, 1/2, 3/4)
            3 3 (1/4, 1/2, 3/4) (1/4, 1/2, -1/4)
        """
        from sage.sets.family import Family
        m = self.cartan_matrix().transpose().inverse()
        Delta = tuple(self.simple_roots())
        zero = Delta[0].parent().zero()
        weights = [sum([m[i,j] * sj for j,sj in enumerate(Delta)], zero)
                   for i in range(len(Delta))]
        for weight in weights:
            weight.set_immutable()
        return Family({ind:weights[i] for i,ind in enumerate(self._index_set)})

    def fundamental_weight(self, i):
        r"""
        Return the fundamental weight with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: [ W.fundamental_weight(i) for i in W.index_set() ]    # optional - gap3
            [(3/4, 1/2, 1/4), (1/2, 1, 1/2), (1/4, 1/2, 3/4)]
        """
        return self.fundamental_weights()[i]

    @cached_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::

            sage: G = ReflectionGroup(['A',3])                          # optional - gap3
            sage: G.coxeter_matrix()                                    # optional - gap3
            [1 3 2]
            [3 1 3]
            [2 3 1]
        """
        return self.cartan_type().coxeter_matrix()

    def permutahedron(self, point=None):
        r"""
        Return the permutahedron of ``self``.

        This is the convex hull of the point ``point`` in the weight
        basis under the action of ``self`` on the underlying vector
        space `V`.

        INPUT:

        - ``point`` -- optional, a point given by its coordinates in
          the weight basis (default is `(1, 1, 1, \ldots)`)

        .. NOTE::

            The result is expressed in the root basis coordinates.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.permutahedron()                                     # optional - gap3
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull
            of 24 vertices

            sage: W = ReflectionGroup(['A',3],['B',2])                  # optional - gap3
            sage: W.permutahedron()                                     # optional - gap3
            A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 192 vertices

        TESTS::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.permutahedron([3,5,8])                              # optional - gap3
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull
            of 24 vertices
        """
        n = self.rank()
        weights = self.fundamental_weights()
        if point is None:
            point = [1] * n
        v = sum(point[i] * wt for i, wt in enumerate(weights))
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(vertices=[v*w.to_matrix() for w in self])

    @cached_method
    def right_coset_representatives(self, J):
        r"""
        Return the right coset representatives of ``self`` for the
        parabolic subgroup generated by the simple reflections in ``J``.

        EXAMPLES::

            sage: W = ReflectionGroup(["A",3])                          # optional - gap3
            sage: for J in Subsets([1,2,3]): W.right_coset_representatives(J)   # optional - gap3
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4)(2,8)(3,5)(7,10)(9,11),
             (1,7)(2,4)(5,6)(8,10)(11,12), (1,2,10)(3,6,5)(4,7,8)(9,12,11),
             (1,4,6)(2,3,11)(5,8,9)(7,10,12), (1,6,4)(2,11,3)(5,9,8)(7,12,10),
             (1,7)(2,6)(3,9)(4,5)(8,12)(10,11),
             (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,2,3,12)(4,5,10,11)(6,7,8,9),
             (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,6)(2,9)(3,8)(5,11)(7,12),
             (1,8)(2,7)(3,6)(4,10)(9,12), (1,10,9,5)(2,6,8,12)(3,11,7,4),
             (1,12,3,2)(4,11,10,5)(6,9,8,7), (1,3)(2,12)(4,10)(5,11)(6,8)(7,9),
             (1,5,12)(2,9,4)(3,10,8)(6,7,11), (1,8,11)(2,5,7)(3,12,4)(6,10,9),
             (1,11,8)(2,7,5)(3,4,12)(6,9,10), (1,12,5)(2,4,9)(3,8,10)(6,11,7),
             (1,3,7,9)(2,11,6,10)(4,8,5,12), (1,9,7,3)(2,10,6,11)(4,12,5,8),
             (1,11)(3,10)(4,9)(5,7)(6,12), (1,9)(2,8)(3,7)(4,11)(5,10)(6,12)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4)(2,8)(3,5)(7,10)(9,11),
             (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,4,6)(2,3,11)(5,8,9)(7,10,12),
             (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,2,3,12)(4,5,10,11)(6,7,8,9),
             (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,6)(2,9)(3,8)(5,11)(7,12),
             (1,3)(2,12)(4,10)(5,11)(6,8)(7,9),
             (1,5,12)(2,9,4)(3,10,8)(6,7,11), (1,3,7,9)(2,11,6,10)(4,8,5,12)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,7)(2,4)(5,6)(8,10)(11,12),
             (1,4,6)(2,3,11)(5,8,9)(7,10,12),
             (1,7)(2,6)(3,9)(4,5)(8,12)(10,11),
             (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,2,3,12)(4,5,10,11)(6,7,8,9),
             (1,10,9,5)(2,6,8,12)(3,11,7,4), (1,12,3,2)(4,11,10,5)(6,9,8,7),
             (1,8,11)(2,5,7)(3,12,4)(6,10,9), (1,12,5)(2,4,9)(3,8,10)(6,11,7),
             (1,11)(3,10)(4,9)(5,7)(6,12)]
            [(), (1,4)(2,8)(3,5)(7,10)(9,11), (1,7)(2,4)(5,6)(8,10)(11,12),
             (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,6,4)(2,11,3)(5,9,8)(7,12,10),
             (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,5,9,10)(2,12,8,6)(3,4,7,11),
             (1,8)(2,7)(3,6)(4,10)(9,12), (1,12,3,2)(4,11,10,5)(6,9,8,7),
             (1,3)(2,12)(4,10)(5,11)(6,8)(7,9),
             (1,11,8)(2,7,5)(3,4,12)(6,9,10), (1,9,7,3)(2,10,6,11)(4,12,5,8)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4,6)(2,3,11)(5,8,9)(7,10,12),
             (1,2,3,12)(4,5,10,11)(6,7,8,9)]
            [(), (1,4)(2,8)(3,5)(7,10)(9,11), (1,2,10)(3,6,5)(4,7,8)(9,12,11),
             (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,5,9,10)(2,12,8,6)(3,4,7,11),
             (1,3)(2,12)(4,10)(5,11)(6,8)(7,9)]
            [(), (1,7)(2,4)(5,6)(8,10)(11,12), (1,10,2)(3,5,6)(4,8,7)(9,11,12),
             (1,12,3,2)(4,11,10,5)(6,9,8,7)]
            [()]
        """
        from sage.combinat.root_system.reflection_group_complex import _gap_return
        J_inv = [self._index_set_inverse[j] + 1 for j in J]
        S = str(gap3('ReducedRightCosetRepresentatives(%s,ReflectionSubgroup(%s,%s))' % (self._gap_group._name, self._gap_group._name, J_inv)))
        return sage_eval(_gap_return(S), locals={'self': self})

    class Element(ComplexReflectionGroup.Element):

        @lazy_attribute
        def _reduced_word(self):
            r"""
            Computes a reduced word and stores it into ``self._reduced_word``.
            The words are in ``range(n)`` and not in the index set.

            TESTS::

                sage: W = ReflectionGroup(['A',2])                      # optional - gap3
                sage: [w._reduced_word for w in W]                      # optional - gap3
                [[], [1], [0], [0, 1], [1, 0], [0, 1, 0]]
            """
            return reduced_word_c(self.parent(),self)

        def reduced_word_in_reflections(self):
            r"""
            Return a word in the reflections to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2], index_set=['a','b'], reflection_index_set=['A','B','C']) # optional - gap3
                sage: [(w.reduced_word(), w.reduced_word_in_reflections()) for w in W]  # optional - gap3
                [([], []),
                 (['b'], ['B']),
                 (['a'], ['A']),
                 (['a', 'b'], ['A', 'B']),
                 (['b', 'a'], ['A', 'C']),
                 (['a', 'b', 'a'], ['C'])]

            .. SEEALSO:: :meth:`reduced_word`
            """
            if self.is_one():
                return []

            W = self.parent()
            r = self.reflection_length()
            R = W.reflections()
            I = W.reflection_index_set()
            word = []
            while r > 0:
                for i in I:
                    w = R[i]._mul_(self)
                    if w.reflection_length() < r:
                        word += [i]
                        r -= 1
                        self = w
                        break
            return word

        def length(self):
            r"""
            Return the length of ``self`` in generating reflections.

            This is the minimal numbers of generating reflections needed
            to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print("%s %s"%(w.reduced_word(), w.length())) # optional - gap3
                [] 0
                [2] 1
                [1] 1
                [1, 2] 2
                [2, 1] 2
                [1, 2, 1] 3
            """
            return len(self._reduced_word)

        def has_left_descent(self, i):
            r"""
            Return whether ``i`` is a left descent of ``self``.

            This is done by testing whether ``i`` is mapped by ``self``
            to a negative root.

            EXAMPLES::

                sage: W = ReflectionGroup(["A",3])                      # optional - gap3
                sage: s = W.simple_reflections()                        # optional - gap3
                sage: (s[1]*s[2]).has_left_descent(1)                   # optional - gap3
                True
                sage: (s[1]*s[2]).has_left_descent(2)                   # optional - gap3
                False
            """
            W = self.parent()
            return self(W._index_set_inverse[i]+1) > W._number_of_reflections

        def has_descent(self, i, side='left', positive=False):
            r"""
            Return whether ``i`` is a descent (or ascent) of ``self``.

            This is done by testing whether ``i`` is mapped by ``self``
            to a negative root.

            INPUT:

            - ``i`` -- an index of a simple reflection
            - ``side`` (default: ``'right'``) -- ``'left'`` or ``'right'``
            - ``positive`` (default: ``False``) -- a boolean

            EXAMPLES::

                sage: W = ReflectionGroup(["A",3])                      # optional - gap3
                sage: s = W.simple_reflections()                        # optional - gap3
                sage: (s[1]*s[2]).has_descent(1)                        # optional - gap3
                True
                sage: (s[1]*s[2]).has_descent(2)                        # optional - gap3
                False
            """
            if not isinstance(positive, bool):
                raise TypeError("%s is not a boolean"%(bool))

            if i not in self.parent().index_set():
                raise ValueError("the given index %s is not in the index set"%i)

            negative = not positive

            if side == 'left':
                return self.has_left_descent(i) is negative
            elif side == 'right':
                return self.has_right_descent(i) is negative
            else:
                raise ValueError("the method 'has_descent' needs the input 'side' to be either 'left' or 'right'")

        def act_on_root(self, root, side="right"):
            r"""
            Return the root obtained by applying ``self`` on ``root``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print("%s %s"%(w.reduced_word(),              # optional - gap3
                ....:           [w.act_on_root(beta) for beta in W.positive_roots()]))  # optional - gap3
                [] [(1, 0), (0, 1), (1, 1)]
                [2] [(1, 1), (0, -1), (1, 0)]
                [1] [(-1, 0), (1, 1), (0, 1)]
                [1, 2] [(-1, -1), (1, 0), (0, -1)]
                [2, 1] [(0, 1), (-1, -1), (-1, 0)]
                [1, 2, 1] [(0, -1), (-1, 0), (-1, -1)]

                sage: elt = W.from_reduced_word([1,2])                  # optional - gap3
                sage: [ elt.act_on_root(beta, side="left") for beta in W.positive_roots() ]  # optional - gap3
                [(0, 1), (-1, -1), (-1, 0)]
            """
            Phi = self.parent().roots()
            if side == "left":
                w = ~self
            elif side == "right":
                w = self
            else:
                raise ValueError("the action on roots must be on the left or on the right")

            return Phi[w(Phi.index(root)+1) - 1]

        def inversion_set(self, side="right"):
            r"""
            Return the inversion set of ``self``.

            This is the set `\{\beta \in \Phi^+ : s(\beta) \in \Phi^-\}`,
            where `s` is ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print("%s %s"%(w.reduced_word(), w.inversion_set()))  # optional - gap3
                [] []
                [2] [(0, 1)]
                [1] [(1, 0)]
                [1, 2] [(1, 0), (1, 1)]
                [2, 1] [(0, 1), (1, 1)]
                [1, 2, 1] [(0, 1), (1, 0), (1, 1)]

                sage: W.from_reduced_word([1,2]).inversion_set(side="left") # optional - gap3
                [(0, 1), (1, 1)]
            """
            Phi_plus = set(self.parent().positive_roots())
            if side == "left":
                w = ~self
            elif side == "right":
                w = self
            else:
                raise ValueError("the action on roots must be on the left or on the right")

            return [root for root in Phi_plus if w.act_on_root(root) not in Phi_plus]

        @cached_in_parent_method
        def right_coset_representatives(self):
            r"""
            Return the right coset representatives of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     rcr = w.right_coset_representatives()         # optional - gap3
                ....:     print("%s %s"%(w.reduced_word(),              # optional - gap3
                ....:                    [v.reduced_word() for v in rcr]))  # optional - gap3
                [] [[], [2], [1], [2, 1], [1, 2], [1, 2, 1]]
                [2] [[], [2], [1]]
                [1] [[], [1], [1, 2]]
                [1, 2] [[]]
                [2, 1] [[]]
                [1, 2, 1] [[], [2], [2, 1]]
            """
            from sage.combinat.root_system.reflection_group_complex import _gap_return
            W = self.parent()
            T = W.reflections()
            T_fix = [i + 1 for i in T.keys()
                     if self.fix_space().is_subspace(T[i].fix_space())]
            S = str(gap3('ReducedRightCosetRepresentatives(%s,ReflectionSubgroup(%s,%s))' % (W._gap_group._name, W._gap_group._name, T_fix)))
            return sage_eval(_gap_return(S, coerce_obj='W'),
                             locals={'self': self, 'W': W})

        def left_coset_representatives(self):
            r"""
            Return the left coset representatives of ``self``.

            .. SEEALSO:: :meth:`right_coset_representatives`

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     lcr = w.left_coset_representatives()          # optional - gap3
                ....:     print("%s %s"%(w.reduced_word(),              # optional - gap3
                ....:                    [v.reduced_word() for v in lcr]))  # optional - gap3
                [] [[], [2], [1], [1, 2], [2, 1], [1, 2, 1]]
                [2] [[], [2], [1]]
                [1] [[], [1], [2, 1]]
                [1, 2] [[]]
                [2, 1] [[]]
                [1, 2, 1] [[], [2], [1, 2]]
            """
            return [ (~w) for w in self.right_coset_representatives() ]

class IrreducibleRealReflectionGroup(RealReflectionGroup, IrreducibleComplexReflectionGroup):

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: for i in [2..7]: ReflectionGroup(["I", i])             # optional - gap3
            Reducible real reflection group of rank 2 and type A1 x A1
            Irreducible real reflection group of rank 2 and type A2
            Irreducible real reflection group of rank 2 and type C2
            Irreducible real reflection group of rank 2 and type I2(5)
            Irreducible real reflection group of rank 2 and type G2
            Irreducible real reflection group of rank 2 and type I2(7)
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible real reflection group of rank %s and type %s'%(self._rank,type_str)

    class Element(RealReflectionGroup.Element, IrreducibleComplexReflectionGroup.Element):
        pass

