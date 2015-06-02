r"""
Coxeter groups

AUTHORS:

- Christian Stump

.. note::

    - For definitions and classification types of finite complex reflection groups, see http://en.wikipedia.org/wiki/Complex_reflection_group.
    - Uses the GAP3 package *chevie*.

Version: 2011-04-26

EXAMPLES::


"""
#*****************************************************************************
#       Copyright (C) 2011 Christian Stump <christian.stump at lacim.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
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
from sage.combinat.root_system.complex_reflection_group import FiniteComplexReflectionGroup, IrreducibleFiniteComplexReflectionGroup, is_chevie_available
from sage.categories.coxeter_groups import CoxeterGroups
from sage.combinat.root_system.cartan_matrix import CartanMatrix

from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField

UCF = UniversalCyclotomicField()
E = UCF.gen

class FiniteCoxeterGroupChevie(FiniteComplexReflectionGroup):
    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None):
        r"""

        TESTS::

            sage: W = CoxeterGroups().example()
            sage: TestSuite(W).run()
        """
        W_types = tuple( tuple( W_type ) if isinstance(W_type,(list,tuple)) else W_type for W_type in W_types )
        cartan_types = []
        for W_type in W_types:
            W_type = CartanType(W_type)
            assert W_type.is_finite()
            assert W_type.is_irreducible()
            cartan_types.append( W_type )
        if len(W_types) == 1:
            cls = IrreducibleFiniteComplexReflectionGroup
        else:
            cls = FiniteComplexReflectionGroup
        cls.__init__(self, W_types, index_set=index_set,
                                    hyperplane_index_set=hyperplane_index_set,
                                    reflection_index_set=reflection_index_set,
                                    is_coxeter_group = True)
        N = self.nr_reflections()
        self._is_positive_root = [None] + [ True ] * N + [False]*N 

    __iter__ = CoxeterGroups.ParentMethods.__iter__.__func__

    def _repr_(self):
        r"""
        Returns the string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.coxeter_group_chevie import CoxeterGroupChevie
            sage: W = CoxeterGroupChevie(['A',3],['B',2],['I',5],['I',6]); W
            Reducible finite Coxeter group of rank 9 and type A3 x B2 x I2(5) x G2
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible finite Coxeter group of rank %s and type %s'%(self._rank,type_str)

    @cached_method
    def bipartite_index_set(self):
        r"""
        Returns the bipartite index set of a finite real reflection group.

        EXAMPLES::

            sage: from sage.combinat.root_system.coxeter_group_chevie import CoxeterGroupChevie
            sage: W = CoxeterGroupChevie(["A",5])
            sage: W.bipartite_index_set()
            [[0, 2, 4], [1, 3]]
 
            sage: W = CoxeterGroupChevie(["A",5],index_set=['a','b','c','d','e'])
            sage: W.bipartite_index_set()
            [['a', 'c', 'e'], ['b', 'd']]
        """
        index_family = self._index_set
        keys = index_family.keys()
        L,R = self._gap_group.BipartiteDecomposition().sage()
        L = [ i for i in keys if index_family[i]+1 in L ]
        R = [ i for i in keys if index_family[i]+1 in R ]
        return [L,R]

    def irreducible_components(self):
        r"""
        Returns a list containing the irreducible components of ``self`` as finite reflection groups.

        EXAMPLES::

            tba
        """
        if self.nr_irreducible_components() == 1:
            irr_comps = [self]
        else:
            irr_comps = []
            for W_type in self._type:
                W_str = [ W_type["series"], W_type["rank"] ]
                if W_type["series"] == "I":
                    W_str[1] = W_type["bond"]
                irr_comps.append( CoxeterGroupChevie(W_str) )
        return irr_comps

    def cartan_type(self):
        if len(self._type) == 1:
            ct = self._type[0]
            return CartanType([ct['series'],ct['rank']])
        else:
            return [ W.cartan_type() for W in self.irreducible_components() ]

    @cached_method
    def cartan_matrix(self):
        from sage.matrix.constructor import matrix
        return matrix(self._gap_group.CartanMat().sage())

    def simple_root(self,i):
        return self.simple_roots()[self._index_set[i]]

    def almost_positive_roots(self):
        return [ -beta for beta in self.simple_roots() ] + self.positive_roots()

    def root_to_reflection(self,root):
        Phi = self.roots()
        R = self.reflections()
        i = Phi.index(root)+1
        j = Phi.index(-root)+1
        for r in R:
            if r(i) == j:
                return r
        raise ValueError, "There is a bug in root_to_reflection!"

    def reflection_to_positive_root(self,r):
        Phi = self.roots()
        N = len(Phi)/2
        for i in range(1,N+1):
            if r(i) == i+N:
                return Phi[i-1]
        raise ValueError, "There is a bug in reflection_to_positive_root!"

    @cached_method
    def fundamental_weights(self):
        m = self.cartan_matrix().transpose().inverse()
        S = self.simple_roots()
        zero = S[0] - S[0]
        weights = [ sum( [ m[i,j] * S[j] for j in range(len(S)) ], zero ) for i in range(len(S)) ]
        for weight in weights:
            weight.set_immutable()
        return weights

    def fundamental_weight(self,i):
        return self.fundamental_weights()[self._index_set[i]]

    def permutahedron(self,coefficients=None):
        n = self.rank()
        weights = self.fundamental_weights()
        if coefficients is None:
            coefficients = [1]*n
        v = sum( coefficients[i] * weights[i] for i in range(n) )
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron( vertices=[ v*(~w).as_matrix() for w in self] )

    class Element(FiniteComplexReflectionGroup.Element):

        has_descent = CoxeterGroups.ElementMethods.has_descent.__func__
        reduced_word = cached_in_parent_method(CoxeterGroups.ElementMethods.reduced_word.__func__)

        def has_left_descent(self, i):
            r"""
            Returns whether ``i`` is a descent of ``self`` by testing
            whether ``i`` is mapped to a negative root.

            EXAMPLES::

                sage: from sage.combinat.root_system.coxeter_group_chevie import CoxeterGroupChevie
                sage: W = CoxeterGroupChevie(["A",3])
                sage: s = W.simple_reflections()
                sage: (s[1]*s[2]).has_left_descent(1)
                True
                sage: (s[1]*s[2]).has_left_descent(2)
                False
            """
            W = self.parent()
            assert i in W.index_set()
            return not W._is_positive_root[self(W._index_set[i]+1)]

        def act_on_root(self,root):
            Phi = self.parent().roots()
            return Phi[ (~self)(Phi.index(root)+1)-1 ]

        def inversion_set(self):
            Phi_plus = set(self.parent().positive_roots())
            return [ root for root in Phi_plus if self.act_on_root(root) not in Phi_plus ]

class IrreducibleFiniteCoxeterGroupChevie(FiniteCoxeterGroupChevie, IrreducibleFiniteComplexReflectionGroup):

    def _repr_(self):
        r"""
        Returns the string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.coxeter_group_chevie import CoxeterGroupChevie
            sage: for i in [2..7]: print CoxeterGroupChevie(["I",i])
            Reducible finite Coxeter group of rank 2 and type A1 x A1
            Irreducible finite Coxeter group of rank 2 and type A2
            Irreducible finite Coxeter group of rank 2 and type B2
            Irreducible finite Coxeter group of rank 2 and type I2(5)
            Irreducible finite Coxeter group of rank 2 and type G2
            Irreducible finite Coxeter group of rank 2 and type I2(7)
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible finite Coxeter group of rank %s and type %s'%(self._rank,type_str)

    class Element(FiniteCoxeterGroupChevie.Element,IrreducibleFiniteComplexReflectionGroup.Element):
        pass

def CoxeterGroupChevie(*args,**kwds):
    """
    INPUT:

     - every argument should be a finite cartan type, or coercible into; see :class:`CartanType`)

    OUTPUT:

    Returns the Coxeter group as a finite reflection group, see :func:`ComplexReflectionGroup`.

    .. warning:: works only if the GAP3 package Chevie is available, use :func:`WeylGroup` otherwise.

    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_group_chevie import CoxeterGroupChevie
        sage: W = CoxeterGroupChevie(["A",2]); W                  # optional (requires chevie)
        Permutation Group with generators [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6)]

        sage: W.category()                                  # optional (requires chevie)
        Join of Category of finite permutation groups and Category of finite coxeter groups
    """
    assert is_chevie_available()
    gap3.load_package("chevie")

    W_types = []
    for arg in args:
        if type(arg) is list:
            X = tuple(arg)
        else:
            X = arg
        assert is_Matrix(X) or isinstance(X,CartanMatrix) or isinstance(X,tuple), "The input is not valid."
        if X == ('I',2):
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
                raise ValueError, 'The keyword %s must be a list, tuple, or dict'%index_set_kwd

    if len(W_types) == 1:
        cls = IrreducibleFiniteCoxeterGroupChevie
    else:
        cls = FiniteCoxeterGroupChevie
    return cls(tuple(W_types), index_set=kwds.get('index_set', None),
                               hyperplane_index_set=kwds.get('hyperplane_index_set', None),
                               reflection_index_set=kwds.get('reflection_index_set', None) )
