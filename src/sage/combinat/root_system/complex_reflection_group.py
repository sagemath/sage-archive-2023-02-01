r"""
Finite (complex) reflection groups

AUTHORS:

- Christian Stump

.. note::

    - For definitions and classification types of finite complex reflection groups, see http://en.wikipedia.org/wiki/Complex_reflection_group.
    - Uses the GAP3 package *chevie*.

Version: 2011-04-26

TODO:

- Element class should be unique to be able to work with large groups without creating elements multiple times
"""
#*****************************************************************************
#       Copyright (C) 2015 Christian Stump <christian.stump at lacim.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.all import prod, add
from sage.misc.cachefunc import cached_function, cached_method, cached_in_parent_method
from sage.categories.category import Category
from sage.categories.permutation_groups import PermutationGroups
from sage.categories.complex_reflection_groups import ComplexReflectionGroups, WellGeneratedComplexReflectionGroups
from sage.categories.coxeter_groups import CoxeterGroups
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
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
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
from sage.rings.arith import gcd, lcm
from sage.modules.free_module_element import vector
from sage.combinat.root_system.cartan_matrix import CartanMatrix

UCF = UniversalCyclotomicField()
E = UCF.gen

class FiniteComplexReflectionGroup(UniqueRepresentation, PermutationGroup_generic):

    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None, is_coxeter_group=False):
        r"""
        TESTS::

            sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
            sage: W = ComplexReflectionGroups().example()
            sage: TestSuite(W).run()
        """
        W_components = []
        reflection_type = []
        for W_type in W_types:
            if W_type == (1,1,1):
                raise ValueError, "The one element group is not considered a reflection group."
            elif W_type in ZZ:
                call_str = 'ComplexReflectionGroup(%s)'%W_type
            elif isinstance(W_type,CartanMatrix):
                call_str = 'PermRootGroup(IdentityMat(%s),%s)'%(W_type._rank,str(W_type._M._gap_()))
            elif is_Matrix(W_type):
                call_str = 'PermRootGroup(IdentityMat(%s),%s)'%(W_type._rank,str(W_type._gap_()))
            elif is_coxeter_group:
                if W_type[0] == "I":
                    call_str = 'CoxeterGroup("I",2,%s)'%W_type[1]
                else:
                    call_str = 'CoxeterGroup("%s",%s)'%W_type
            else:
                call_str = 'ComplexReflectionGroup%s'%str(W_type)
            W_components.append(gap3(call_str))
            X = list(W_components[-1].ReflectionType())
            if len(X) > 1:
                raise ValueError, "Your input data %s is not valid."%W_type
            X = X[0]
            type_dict = dict()
            type_dict["series"] = X.series.sage()
            type_dict["rank"] = X.rank.sage()
            type_dict["indices"] = X.indices.sage()
            if hasattr(X.ST,"sage"):
                type_dict["ST"] = X.ST.sage()
            elif hasattr(X.p,"sage") and hasattr(X.q,"sage"):
                type_dict["ST"] = ( X.p.sage(), X.q.sage(), X.rank.sage() )
            elif hasattr(X.bond,"sage"):
                type_dict["bond"] = X.bond.sage()
            if type_dict["series"] == "B" and X.cartanType.sage() == 1:
                type_dict["series"] = "C"
            reflection_type.append( type_dict )

        self._type = reflection_type
        self._gap_group = prod(W_components)
        generators = [str(x) for x in self._gap_group.generators]
        self._index_set = index_set
        self._hyperplane_index_set = hyperplane_index_set
        self._reflection_index_set = reflection_index_set

        self._elements = None
        self._conjugacy_classes = {}
        self._conjugacy_classes_representatives = None
        self._reflection_representation = None

        self._rank = self._gap_group.rank.sage()
        if len(generators) == self._rank:
            category = WellGeneratedComplexReflectionGroups()
            if all(str(W_comp).find('CoxeterGroup') >= 0 for W_comp in W_components):
                category = Category.join([category,CoxeterGroups()])
        else:
            category = ComplexReflectionGroups()

        category = Category.join([category,PermutationGroups()]).Finite()
        #if len(self._type) == 1:
            #category = category.Irreducible()

        PermutationGroup_generic.__init__(self, gens = generators, canonicalize=False, category = category)

        from sage.sets.family import Family
        l_set = range(len(self.gens()))
        if self._index_set is None:
            self._index_set = Family(dict( (i,i) for i in l_set))
        else:
            assert sorted(self._index_set.values()) == l_set
        self._index_set_inverse = self._index_set.inverse_family()
        Nstar_set = range(self.nr_reflecting_hyperplanes())
        if self._hyperplane_index_set is None:
            self._hyperplane_index_set = Family(dict( (i,i) for i in Nstar_set))
        else:
            assert sorted(self._hyperplane_index_set.values()) == Nstar_set
        self._hyperplane_index_set_inverse = self._hyperplane_index_set.inverse_family()
        N_set = range(self.nr_reflections())
        if self._reflection_index_set is None:
            self._reflection_index_set = Family(dict( (i,i) for i in N_set))
        else:
            assert sorted(self._reflection_index_set.values()) == N_set
        self._reflection_index_set_inverse = self._reflection_index_set.inverse_family()

    def _irrcomp_repr_(self,W_type):
        type_str = ''
        if "ST" in W_type:
            if W_type["ST"] in ZZ:
                type_str += "ST" + str(W_type["ST"])
            else:
                type_str += 'G' + str(W_type["ST"]).replace(' ','')
        else:
            type_str += str(W_type["series"])
            if W_type["series"] == "I":
                type_str += '2(' + str(W_type["bond"]) + ')'
            else:
                type_str += str(W_type["rank"])
        return type_str

    def _repr_(self):
        r"""
        Returns the string representation of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup(25,[4,1,4],[1,1,4],[5,5,2]); W
            Reducible finite complex reflection group of rank 12 and type ST25 x G(4,1,4) x A3 x I2(5)
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible finite complex reflection group of rank %s and type %s'%(self._rank,type_str)

    def __iter__(self):
        r"""
        Returns an iterator going through all elements in ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: for w in W.__iter__(): print w
            ()
            (1,4)(2,3)(5,6)
            (1,3)(2,5)(4,6)
            (1,6,2)(3,5,4)
            (1,2,6)(3,4,5)
            (1,5)(2,4)(3,6)
        """
        if self._elements is not None and len(self._elements) == self.cardinality():
            for w in self._elements:
                yield w
        else:
            self._elements = []
            inv_dict = dict( (self._index_set[i],i) for i in self._index_set.keys() )
            for w,word in self.magma_closure_iter(I=self.gens(),return_word=True):
                if w._reduced_word is None:
                    w._reduced_word = Word( inv_dict[i] for i in word )
                self._elements.append(w)
                yield w

    __len__ = ComplexReflectionGroups.Finite.ParentMethods.cardinality.__func__

    @cached_method
    def index_set(self):
        r"""
        Returns the index set of the simple reflections of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: W.index_set()
            [0, 1, 2]
            sage: W = ComplexReflectionGroup((1,1,4),index_set=[1,3,'asdf'])
            sage: W.index_set()
            [1, 3, 'asdf']
            sage: W = ComplexReflectionGroup((1,1,4),index_set={'a':0,'b':1,'c':2})
            sage: W.index_set()
            ['a', 'b', 'c']
        """
        return [ self._index_set_inverse[i] for i in xrange(len(self._index_set)) ]

    def series(self):
        return [ self._type[i]['series'] for i in range(len(self._type)) ]

    @cached_method
    def hyperplane_index_set(self):
        r"""
        Returns the index set of the hyperplanes of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: W.hyperplane_index_set()
            [0, 1, 2, 3, 4, 5]
            sage: W = ComplexReflectionGroup((1,1,4),hyperplane_index_set=[1,3,'asdf',7,9,11])
            sage: W.hyperplane_index_set()
            [1, 3, 'asdf', 7, 9, 11]
            sage: W = ComplexReflectionGroup((1,1,4),hyperplane_index_set={'a':0,'b':1,'c':2,'d':3,'e':4,'f':5})
            sage: W.hyperplane_index_set()
            ['a', 'b', 'c', 'd', 'e', 'f']
        """
        return [ self._hyperplane_index_set_inverse[i] for i in xrange(len(self._hyperplane_index_set)) ]

    @cached_method
    def distinguished_reflections(self):
        r"""
        Returns a finite family containing the distinguished reflections of ``self``,
        indexed by ``self.hyperplane_index_set()``.
        These are the reflections in ``self`` acting on the complement
        of the fixed hyperplane `H` as `\operatorname{exp}(2 \pi i / n)`, where `n`
        is the order of the reflection subgroup fixing `H`.

       EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.distinguished_reflections()
            Finite family {0: (1,4)(2,3)(5,6), 1: (1,3)(2,5)(4,6), 2: (1,5)(2,4)(3,6)}

            sage: W = ComplexReflectionGroup((1,1,3),hyperplane_index_set=['a','b','c'])
            sage: W.distinguished_reflections()
            Finite family {'a': (1,4)(2,3)(5,6), 'c': (1,5)(2,4)(3,6), 'b': (1,3)(2,5)(4,6)}

            sage: W = ComplexReflectionGroup((3,1,1))
            sage: W.distinguished_reflections()
            Finite family {0: (1,2,3)}

            sage: W = ComplexReflectionGroup((1,1,3),(3,1,2))
            sage: W.distinguished_reflections()
            Finite family {0: (1,6)(2,5)(7,8), 1: (1,5)(2,7)(6,8), 2: (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30), 3: (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30), 4: (1,7)(2,6)(5,8), 5: (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26), 6: (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29), 7: (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28)}
        """
        from sage.sets.family import Family
        # imports all distinguished reflections from gap, the Set is used as every such appears multiple times.
        T = [ self(str(r)) for r in self._gap_group.Reflections() ]
        # makes sure that the simple reflections come first
        gens = self.gens()
        R = [ t for t in gens ]
        for t in T:
            if t not in R:
                R.append(t)
        return Family(self._hyperplane_index_set.keys(), lambda i: R[self._hyperplane_index_set[i]] )

    def distinguished_reflection(self, i):
        r"""
        Returns the ``i``-th distinguished reflection of ``self``.
        These are the reflections in ``self`` acting on the complement
        of the fixed hyperplane `H` as `\operatorname{exp}(2 \pi i / n)`, where `n`
        is the order of the reflection subgroup fixing `H`.

       EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.distinguished_reflection(0)
            (1,4)(2,3)(5,6)
            sage: W.distinguished_reflection(1)
            (1,3)(2,5)(4,6)
            sage: W.distinguished_reflection(2)
            (1,5)(2,4)(3,6)

            sage: W = ComplexReflectionGroup((3,1,1),hyperplane_index_set=['a'])
            sage: W.distinguished_reflection('a')
            (1,2,3)

            sage: W = ComplexReflectionGroup((1,1,3),(3,1,2))
            sage: for i in range(W.nr_reflecting_hyperplanes()): print W.distinguished_reflection(i)
            (1,6)(2,5)(7,8)
            (1,5)(2,7)(6,8)
            (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30)
            (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30)
            (1,7)(2,6)(5,8)
            (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26)
            (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29)
            (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28)
        """
        assert i in self.hyperplane_index_set()
        return self.distinguished_reflections()[i]

    @cached_method
    def reflection_index_set(self):
        r"""
        Returns the index set of the reflections of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: W.reflection_index_set()
            [0, 1, 2, 3, 4, 5]
            sage: W = ComplexReflectionGroup((1,1,4),reflection_index_set=[1,3,'asdf',7,9,11])
            sage: W.reflection_index_set()
            [1, 3, 'asdf', 7, 9, 11]
            sage: W = ComplexReflectionGroup((1,1,4),reflection_index_set={'a':0,'b':1,'c':2,'d':3,'e':4,'f':5})
            sage: W.reflection_index_set()
            ['a', 'b', 'c', 'd', 'e', 'f']
        """
        return [ self._reflection_index_set_inverse[i] for i in xrange(len(self._reflection_index_set)) ]

    @cached_method
    def reflections(self):
        r"""
        Returns a finite family containing the reflections of ``self``,
        indexed by ``self.reflection_index_set()``.

       EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.reflections()
            Finite family {0: (1,4)(2,3)(5,6), 1: (1,3)(2,5)(4,6), 2: (1,5)(2,4)(3,6)}

            sage: W = ComplexReflectionGroup((1,1,3),reflection_index_set=['a','b','c'])
            sage: W.reflections()
            Finite family {'a': (1,4)(2,3)(5,6), 'c': (1,5)(2,4)(3,6), 'b': (1,3)(2,5)(4,6)}

            sage: W = ComplexReflectionGroup((3,1,1))
            sage: W.reflections()
            Finite family {0: (1,2,3), 1: (1,3,2)}

            sage: W = ComplexReflectionGroup((1,1,3),(3,1,2))
            sage: W.reflections()
            Finite family {0: (1,6)(2,5)(7,8), 1: (1,5)(2,7)(6,8), 2: (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30), 3: (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30), 4: (1,7)(2,6)(5,8), 5: (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26), 6: (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29), 7: (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28), 8: (3,15,9)(4,16,10)(12,23,17)(14,24,18)(20,29,25)(21,26,22)(27,30,28), 9: (4,27,21)(10,28,22)(11,19,13)(12,20,14)(16,30,26)(17,25,18)(23,29,24)}
        """
        from sage.sets.family import Family
        T = self.distinguished_reflections().values()
        for i in xrange(self.nr_reflecting_hyperplanes()):
            T.extend( [ T[i]**j for j in range(2,T[i].order()) ] )
        return Family(self._reflection_index_set.keys(), lambda i: T[self._reflection_index_set[i]] )

    def reflection(self,i):
        r"""
        Returns the ``i``-th reflection of ``self``

       EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.reflection(0)
            (1,4)(2,3)(5,6)
            sage: W.reflection(1)
            (1,3)(2,5)(4,6)
            sage: W.reflection(2)
            (1,5)(2,4)(3,6)

            sage: W = ComplexReflectionGroup((3,1,1),reflection_index_set=['a','b'])
            sage: W.reflection('a')
            (1,2,3)
            sage: W.reflection('b')
            (1,3,2)
        """
        assert i in self.reflection_index_set()
        return self.reflections()[i]

    def reflection_character(self):
        r"""
        Returns the reflection characters of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.reflection_character()
            [2, 0, -1]
        """
        return self._gap_group.ReflectionCharacter().sage()

    def is_crystallographic(self):
        r"""
        Returns True if self is crystallographic, i.e., if the reflection representation of ``self``
        is defined over the rationals.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3)); W
            Irreducible finite complex reflection group of rank 2 and type A2
            sage: W.is_crystallographic()
            True

            sage: W = ComplexReflectionGroup((2,1,3)); W
            Irreducible finite complex reflection group of rank 3 and type B3
            sage: W.is_crystallographic()
            True

            sage: W = ComplexReflectionGroup(23); W
            Irreducible finite complex reflection group of rank 3 and type H3
            sage: W.is_crystallographic()
            False

            sage: W = ComplexReflectionGroup((3,1,3)); W
            Irreducible finite complex reflection group of rank 3 and type G(3,1,3)
            sage: W.is_crystallographic()
            False
        """
        from sage.rings.all import QQ
        return all( t.as_matrix().base_ring() is QQ for t in self.simple_reflections() )

    def _element_class(self):
        r"""
        A temporary workaround for compatibility with Sage's
        permutation groups

        TESTS::

            sage: W = ComplexReflectionGroup(23)                         # optional (require chevie)
            sage: W._element_class() is W.element_class                  # optional (require chevie)
            True
        """
        return self.element_class

    def nr_irreducible_components(self):
        r"""
        Returns the number of irreducible components of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.nr_irreducible_components()
            1

            sage: W = ComplexReflectionGroup((1,1,3),(2,1,3))
            sage: W.nr_irreducible_components()
            2
        """
        return len(self._type)

    def irreducible_components(self):
        r"""
        Returns a list containing the irreducible components of ``self`` as finite reflection groups.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.irreducible_components()
            [Irreducible finite complex reflection group of rank 2 and type A2]

            sage: W = ComplexReflectionGroup((1,1,3),(2,1,3))
            sage: W.irreducible_components()
            [Irreducible finite complex reflection group of rank 2 and type A2,
            Irreducible finite complex reflection group of rank 3 and type B3]
        """
        if self.nr_irreducible_components() == 1:
            irr_comps = [self]
        else:
            irr_comps = []
            for W_type in self._type:
                if W_type["series"] == "A":
                    W_str = (1,1,W_type["rank"]+1)
                elif W_type["series"] == "B":
                    W_str = (2,1,W_type["rank"])
                elif W_type["series"] == "D":
                    W_str = (2,2,W_type["rank"])
                elif W_type["series"] == "E":
                    if W_type["rank"] == 6:
                        W_str = 35
                    elif W_type["rank"] == 7:
                        W_str = 36
                    elif W_type["rank"] == 8:
                        W_str = 37
                elif W_type["series"] == "F":
                    W_str = 28
                elif W_type["series"] == "G":
                    W_str = (6,6,2)
                elif W_type["series"] == "H":
                    if W_type["rank"] == 3:
                        W_str = 23
                    elif W_type["rank"] == 4:
                        W_str = 30
                elif W_type["series"] == "I":
                    W_str = (W_type["bond"],W_type["bond"],2)
                elif "ST" in W_type:
                    W_str = W_type["ST"]

                else:
                    raise ValueError, "not yet implemented"
                irr_comps.append( ComplexReflectionGroup(W_str) )
        return irr_comps

    @cached_method
    def conjugacy_classes_representatives(self):
        r"""
        Returns the shortest representatives of the conjugacy classes of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: [ w.reduced_word() for w in W.conjugacy_classes_representatives() ]
            [word: , word: 0, word: 01]

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: sorted([ w.reduced_word() for w in W.conjugacy_classes_representatives() ])
            [word: , word: 0, word: 01, word: 012, word: 20]

            sage: W = ComplexReflectionGroup((3,1,2))
            sage: [ w.reduced_word() for w in W.conjugacy_classes_representatives() ]
            [word: , word: 0, word: 00, word: 1010, word: 10100, word: 100100, word: 1, word: 01, word: 001]


            sage: W = ComplexReflectionGroup(23)
            sage: [ w.reduced_word() for w in W.conjugacy_classes_representatives() ]
            [word: , word: 0, word: 01, word: 20, word: 12, word: 012, word: 0101, word: 01012, word: 012010121, word: 210120101201010]
        """
        if self._conjugacy_classes_representatives is None:
            S = str(gap3('List(ConjugacyClasses(%s),Representative)'%self._gap_group._name))
            exec('self._conjugacy_classes_representatives='+gap_return(S))
        return self._conjugacy_classes_representatives

    @cached_method
    def conjugacy_classes(self):
        r"""
        Returns the conjugacy classes of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: for C in W.conjugacy_classes(): print C
            frozenset([()])
            frozenset([(1,3)(2,5)(4,6), (1,5)(2,4)(3,6), (1,4)(2,3)(5,6)])
            frozenset([(1,2,6)(3,4,5), (1,6,2)(3,5,4)])

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: sum( len(C) for C in  W.conjugacy_classes() ) == W.cardinality()
            True

            sage: W = ComplexReflectionGroup((3,1,2))
            sage: sum( len(C) for C in  W.conjugacy_classes() ) == W.cardinality()
            True

            sage: W = ComplexReflectionGroup(23)
            sage: sum( len(C) for C in  W.conjugacy_classes() ) == W.cardinality()
            True
       """
        from sage.sets.family import Family
        return Family(self.conjugacy_classes_representatives(), lambda w: w.conjugacy_class() )

    def rank(self):
        r"""
        Returns the rank of ``self``. This is the dimension of the underlying vector space.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.rank()
            2
            sage: W = ComplexReflectionGroup((2,1,3))
            sage: W.rank()
            3
            sage: W = ComplexReflectionGroup((4,1,3))
            sage: W.rank()
            3
            sage: W = ComplexReflectionGroup((4,2,3))
            sage: W.rank()
            3
        """
        return self._rank

    @cached_method
    def degrees(self):
        r"""
        Returns the degrees of ``self`` ordered within each irreducible component of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: W.degrees()
            [2, 3, 4]

            sage: W = ComplexReflectionGroup((2,1,4))
            sage: W.degrees()
            [2, 4, 6, 8]

            sage: W = ComplexReflectionGroup((4,1,4))
            sage: W.degrees()
            [4, 8, 12, 16]

            sage: W = ComplexReflectionGroup((4,2,4))
            sage: W.degrees()
            [4, 8, 12, 8]

            sage: W = ComplexReflectionGroup((4,4,4))
            sage: W.degrees()
            [4, 8, 12, 4]

        Examples of reducible types::

            sage: W = ComplexReflectionGroup((1,1,4),(3,1,2)); W
            Reducible finite complex reflection group of rank 5 and type A3 x G(3,1,2)
            sage: W.degrees()
            [2, 3, 4, 3, 6]

            sage: W = ComplexReflectionGroup((1,1,4),(6,1,12),23) # fails in GAP3
            sage: W.degrees()
            [2, 3, 4, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 2, 6, 10]
        """
        if self.is_irreducible():
            try:
                return self._gap_group.degrees.sage()
            except:
                return self._gap_group.ReflectionDegrees().sage()
        else:
            return add( [comp.degrees() for comp in self.irreducible_components()], [] )

    cardinality = ComplexReflectionGroups.Finite.ParentMethods.cardinality.__func__

    @cached_method
    def codegrees(self):
        r"""
        Returns the codegrees of ``self`` ordered within each irreducible component of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: W.codegrees()
            [0, 1, 2]

            sage: W = ComplexReflectionGroup((2,1,4))
            sage: W.codegrees()
            [0, 2, 4, 6]

            sage: W = ComplexReflectionGroup((4,1,4))
            sage: W.codegrees()
            [0, 4, 8, 12]

            sage: W = ComplexReflectionGroup((4,2,4))
            sage: W.codegrees()
            [0, 4, 8, 12]

            sage: W = ComplexReflectionGroup((4,4,4))
            sage: W.codegrees()
            [8, 0, 4, 8]

            sage: W = ComplexReflectionGroup((1,1,4),(3,1,2))
            sage: W.codegrees()
            [0, 1, 2, 0, 3]
            
            sage: W = ComplexReflectionGroup((1,1,4),(6,1,12),23) # fails in GAP3
            sage: W.codegrees()
            [0, 1, 2, 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 0, 4, 8]
        """
        if self.is_irreducible():
            if self.is_well_generated():
                h = self.coxeter_number()
                return [ h-d for d in reversed(self.degrees()) ]
            else:
                return self._gap_group.ReflectionCoDegrees().sage()
        else:
            return add( [comp.codegrees() for comp in self.irreducible_components()], [] )

    @cached_method
    def reflection_eigenvalues_family(self):
        r"""
        Returns the reflection eigenvalues of ``self`` as a finite family indexed by the class representatives of ``self``.
        An eigenvalue `\zeta_n^k` is returned as the quotient `k/n` in the rationals.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.reflection_eigenvalues_family()
            Finite family {(): [0, 0], (1,4)(2,3)(5,6): [1/2, 0], (1,6,2)(3,5,4): [1/3, 2/3]}

            sage: W = ComplexReflectionGroup((3,1,2))
            sage: W.reflection_eigenvalues_family()
            Finite family {(1,3,9)(2,16,24)(4,20,21)(5,7,13)(6,12,23)(8,19,17)(10,15,22)(11,18,14): [1/3, 1/3], (1,13,9,7,3,5)(2,14,24,18,16,11)(4,6,21,23,20,12)(8,22,17,15,19,10): [1/3, 5/6], (1,7,3,13,9,5)(2,8,16,19,24,17)(4,14,20,11,21,18)(6,15,12,22,23,10): [1/6, 2/3], (1,9,3)(2,10,4)(6,17,11)(8,18,12)(14,23,19)(15,20,16)(21,24,22): [2/3, 0], (1,3,9)(2,4,10)(6,11,17)(8,12,18)(14,19,23)(15,16,20)(21,22,24): [1/3, 0], (1,9,3)(2,20,22)(4,15,24)(5,7,13)(6,18,19)(8,23,11)(10,16,21)(12,14,17): [1/3, 2/3], (1,9,3)(2,24,16)(4,21,20)(5,13,7)(6,23,12)(8,17,19)(10,22,15)(11,14,18): [2/3, 2/3], (): [0, 0], (1,5)(2,6)(3,7)(4,8)(9,13)(10,14)(11,15)(12,16)(17,21)(18,22)(19,20)(23,24): [1/2, 0]}

            sage: W = ComplexReflectionGroup(23)
            sage: W.reflection_eigenvalues_family()
            Finite family {(1,16)(2,5)(4,7)(6,9)(8,10)(11,13)(12,14)(17,20)(19,22)(21,24)(23,25)(26,28)(27,29): [1/2, 0, 0], (): [0, 0, 0], (1,24,17,16,9,2)(3,12,13,18,27,28)(4,21,29,19,6,14)(5,25,26,20,10,11)(7,23,30,22,8,15): [1/6, 1/2, 5/6], (1,8,4)(2,21,3)(5,10,11)(6,18,17)(7,9,12)(13,14,15)(16,23,19)(20,25,26)(22,24,27)(28,29,30): [1/3, 2/3, 0], (1,29,8,7,26,16,14,23,22,11)(2,9,25,3,4,17,24,10,18,19)(5,30,6,13,27,20,15,21,28,12): [3/10, 1/2, 7/10], (1,19,20,2,7)(3,6,11,13,9)(4,5,17,22,16)(8,12,15,14,10)(18,21,26,28,24)(23,27,30,29,25): [1/5, 4/5, 0], (1,16)(2,9)(3,18)(4,10)(5,6)(7,8)(11,14)(12,13)(17,24)(19,25)(20,21)(22,23)(26,29)(27,28): [1/2, 1/2, 0], (1,20,7,19,2)(3,11,9,6,13)(4,17,16,5,22)(8,15,10,12,14)(18,26,24,21,28)(23,30,25,27,29): [2/5, 3/5, 0], (1,23,26,29,22,16,8,11,14,7)(2,10,4,9,18,17,25,19,24,3)(5,21,27,30,28,20,6,12,15,13): [1/10, 1/2, 9/10], (1,16)(2,17)(3,18)(4,19)(5,20)(6,21)(7,22)(8,23)(9,24)(10,25)(11,26)(12,27)(13,28)(14,29)(15,30): [1/2, 1/2, 1/2]}
            """
        from sage.sets.family import Family
        class_representatives = self.conjugacy_classes_representatives()
        Ev_list = self._gap_group.ReflectionEigenvalues().sage()
        return Family( class_representatives, lambda w: Ev_list[class_representatives.index(w)] )

    @cached_method
    def reflection_eigenvalues(self,w,test_class_repr=True):
        r"""
        Returns the reflection eigenvalue of ``w`` in ``self``.

        .. seealso:: :meth:`reflection_eigenvalues_family`

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: for w in W: print w.reduced_word(), W.reflection_eigenvalues(w)
             [0, 0]
            0 [1/2, 0]
            1 [1/2, 0]
            01 [1/3, 2/3]
            10 [1/3, 2/3]
            010 [1/2, 0]
        """
        if test_class_repr:
            w_repr = w.conjugacy_class_representative()
        else:
            w_repr = w
        return self.reflection_eigenvalues_family()[ w_repr ]

    @cached_method
    def simple_roots(self):
        r"""
        Returns the *simple roots* of ``self``. These are the roots
        corresponding to the simple reflections.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.simple_roots()
            [(1, 0), (0, 1)]

            sage: W = ComplexReflectionGroup((1,1,4),(2,1,2))
            sage: W.simple_roots()
            [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 0, 1)]

            sage: W = ComplexReflectionGroup((3,1,2))
            sage: W.simple_roots()
            [(1, 0), (-1, 1)]

            sage: W = ComplexReflectionGroup((1,1,4),(3,1,2))
            sage: W.simple_roots()
            [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, -1, 1)]
        """
        return self.roots()[:len(self.gens())]

    @cached_method
    def simple_coroots(self):
        r"""
        Returns the *simple coroots* of ``self``. These are the coroots
        corresponding to the simple reflections.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.simple_coroots()
            [[2, -1], [-1, 2]]

            sage: W = ComplexReflectionGroup((1,1,4),(2,1,2))
            sage: W.simple_coroots()
            [[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, 0, 0], [0, 0, 0, 2, -2], [0, 0, 0, -1, 2]]

            sage: W = ComplexReflectionGroup((3,1,2))
            sage: W.simple_coroots()
            [[-2*E(3) - E(3)^2, 0], [-1, 1]]

            sage: W = ComplexReflectionGroup((1,1,4),(3,1,2))
            sage: W.simple_coroots()
            [[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, 0, 0], [0, 0, 0, -2*E(3) - E(3)^2, 0], [0, 0, 0, -1, 1]]
        """
        return self._gap_group.simpleCoroots.sage()

    @cached_method
    def independent_roots(self):
        r"""
        Returns a collection of simple roots generating the underlying vector space of ``self``.
        For well-generated groups, these are all simple roots. Otherwise, a linearly independent
        subset of the simple roots is chosen.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.independent_roots()
            [(1, 0), (0, 1)]

            sage: W = ComplexReflectionGroup((4,2,3))
            sage: W.simple_roots()
            [(1, 0, 0), (-E(4), 1, 0), (-1, 1, 0), (0, -1, 1)]
            sage: W.independent_roots()
            [(1, 0, 0), (-E(4), 1, 0), (0, -1, 1)]
        """
        Delta = self.simple_roots()
        if len(Delta) == self.rank():
            basis = Delta
        else:
            basis = []
            for alpha in Delta:
                if Matrix(basis+[alpha]).rank() == len(basis) + 1:
                    basis.append(alpha)
        return basis

    @cached_method
    def base_change_matrix(self):
        r"""
        Returns the base change from the standard basis of the vector space of ``self`` to the basis given by the independent roots of ``self``.

        FIXME:

        - for non-well-generated groups there is a conflict with construction of the matrix for an element

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.base_change_matrix()
            [1 0]
            [0 1]

            sage: W = ComplexReflectionGroup(23)
            sage: W.base_change_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: W = ComplexReflectionGroup((3,1,2))
            sage: W.base_change_matrix()
            [1 0]
            [1 1]

            sage: W = ComplexReflectionGroup((4,2,2))
            sage: W.base_change_matrix()
            [   1    0]
            [E(4)    1]
        """
        return Matrix( self.independent_roots() ).inverse()

    @cached_method
    def roots(self):
        r"""
        Returns all roots corresponding to all reflections of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.roots()
            [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]

            sage: W = ComplexReflectionGroup((3,1,2))
            sage: W.roots()
            [(1, 0), (-1, 1), (E(3), 0), (-E(3), 1), (0, 1), (1, -1), (0, E(3)), (1, -E(3)), (E(3)^2, 0), (-E(3)^2, 1), (E(3), -1), (E(3), -E(3)), (0, E(3)^2), (1, -E(3)^2), (-1, E(3)), (-E(3), E(3)), (E(3)^2, -1), (E(3)^2, -E(3)), (E(3), -E(3)^2), (-E(3)^2, E(3)), (-1, E(3)^2), (-E(3), E(3)^2), (E(3)^2, -E(3)^2), (-E(3)^2, E(3)^2)]

            sage: W = ComplexReflectionGroup((4,2,2))
            sage: W.roots()
            [(1, 0), (-E(4), 1), (-1, 1), (-1, 0), (E(4), 1), (1, 1), (0, -E(4)), (E(4), -1), (E(4), E(4)), (0, E(4)), (E(4), -E(4)), (0, 1), (1, -E(4)), (1, -1), (0, -1), (1, E(4)), (-E(4), 0), (-1, E(4)), (E(4), 0), (-E(4), E(4)), (-E(4), -1), (-E(4), -E(4)), (-1, -E(4)), (-1, -1)]

            sage: W = ComplexReflectionGroup((1,1,4),(3,1,2))
            sage: W.roots()
            [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, -1, 1), (1, 1, 0, 0, 0), (0, 1, 1, 0, 0), (1, 1, 1, 0, 0), (-1, 0, 0, 0, 0), (0, -1, 0, 0, 0), (0, 0, -1, 0, 0), (-1, -1, 0, 0, 0), (0, -1, -1, 0, 0), (-1, -1, -1, 0, 0), (0, 0, 0, E(3), 0), (0, 0, 0, -E(3), 1), (0, 0, 0, 0, 1), (0, 0, 0, 1, -1), (0, 0, 0, 0, E(3)), (0, 0, 0, 1, -E(3)), (0, 0, 0, E(3)^2, 0), (0, 0, 0, -E(3)^2, 1), (0, 0, 0, E(3), -1), (0, 0, 0, E(3), -E(3)), (0, 0, 0, 0, E(3)^2), (0, 0, 0, 1, -E(3)^2), (0, 0, 0, -1, E(3)), (0, 0, 0, -E(3), E(3)), (0, 0, 0, E(3)^2, -1), (0, 0, 0, E(3)^2, -E(3)), (0, 0, 0, E(3), -E(3)^2), (0, 0, 0, -E(3)^2, E(3)), (0, 0, 0, -1, E(3)^2), (0, 0, 0, -E(3), E(3)^2), (0, 0, 0, E(3)^2, -E(3)^2), (0, 0, 0, -E(3)^2, E(3)^2)]
        """
        roots = [ vector(eval(str(root).replace("^","**"))) for root in self._gap_group.roots ]
        for v in roots:
            v.set_immutable()
        return roots

    @cached_method
    def braid_relations(self):
        r"""
        Returns the braid relations of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.braid_relations()
            [[[2, 1, 2], [1, 2, 1]]]

            sage: W = ComplexReflectionGroup((2,1,3))
            sage: W.braid_relations()
            [[[2, 1, 2, 1], [1, 2, 1, 2]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]

            sage: W = ComplexReflectionGroup((2,2,3))
            sage: W.braid_relations()
            [[[2, 1, 2], [1, 2, 1]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]
        """
        return self._gap_group.BraidRelations().sage()

    def fundamental_invariants(self):
        r"""
        Returns the fundamental invariants of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: W.fundamental_invariants()
            [-2*x0^2 + 2*x0*x1 - 2*x1^2, 6*x0^2*x1 - 6*x0*x1^2]
        """
        from sage.rings.polynomial.all import PolynomialRing

        I = [ str(p) for p in gap3('List(Invariants(%s),x->ApplyFunc(x,List([0..%s],i->Mvp(SPrint("x",i)))))'%(self._gap_group._name,self.rank()-1)) ]
        P = PolynomialRing(QQ,['x%s'%i for i in range(0,self.rank())])
        x = P.gens()
        for i in range(len(I)):
            I[i] = I[i].replace('^','**')
            for j in range(len(x)):
                I[i] = I[i].replace('x%s'%j,'*x[%s]'%j)
        I = [ eval(p) for p in I ]
        return I

    def set_reflection_representation(self,refl_repr):
        self.one().as_matrix.clear_cache()
        if set( refl_repr.keys() ) != set( self.index_set() ):
            raise ValueError, "The reflection representation must be defined for the complete index set."
        self._reflection_representation = refl_repr

    class Element(PermutationGroupElement):

        _reduced_word=None

        def apply_simple_reflection_right(self,i):
            r"""
            Returns the product of ``self`` with the ``i``-th simple reflection.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.apply_simple_reflection_right(0)
                (1,4)(2,3)(5,6)
                ()
                (1,2,6)(3,4,5)
                (1,5)(2,4)(3,6)
                (1,3)(2,5)(4,6)
                (1,6,2)(3,5,4)
            """
            assert i in self.parent().index_set()
            gen = self.parent().gens()[self.parent()._index_set[i]]
            return self*gen

        def apply_simple_reflection_left(self,i):
            r"""
            Returns the product of ``self`` with the ``i``-th simple reflection.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.apply_simple_reflection_left(0)
                (1,4)(2,3)(5,6)
                ()
                (1,6,2)(3,5,4)
                (1,3)(2,5)(4,6)
                (1,5)(2,4)(3,6)
                (1,2,6)(3,4,5)
            """
            assert i in self.parent().index_set()
            gen = self.parent().gens()[self.parent()._index_set[i]]
            return gen*self

        @cached_in_parent_method
        def conjugacy_class_representative(self):
            r"""
            Returns a representative of the conjugacy class of ``self``.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.reduced_word(), w.conjugacy_class_representative().reduced_word()
                0 0
                1 0
                01 01
                10 01
                010 0
            """
            W = self.parent()
            for w in W._conjugacy_classes.keys():
                if self in W._conjugacy_classes[w]:
                    return w
            return W.conjugacy_classes_representatives()[ gap3("PositionClass(%s,%s)"%(W._gap_group._name,self)).sage()-1 ]

        def conjugacy_class(self):
            r"""
            Returns the conjugacy class of ``self``.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.conjugacy_class()
                frozenset([()])
                frozenset([(1,3)(2,5)(4,6), (1,5)(2,4)(3,6), (1,4)(2,3)(5,6)])
                frozenset([(1,3)(2,5)(4,6), (1,5)(2,4)(3,6), (1,4)(2,3)(5,6)])
                frozenset([(1,2,6)(3,4,5), (1,6,2)(3,5,4)])
                frozenset([(1,2,6)(3,4,5), (1,6,2)(3,5,4)])
                frozenset([(1,3)(2,5)(4,6), (1,5)(2,4)(3,6), (1,4)(2,3)(5,6)])
            """
            W = self.parent()
            if self not in W.conjugacy_classes_representatives():
                self = self.conjugacy_class_representative()
            if self in W._conjugacy_classes.keys():
                return W._conjugacy_classes[self]
            gens = W.simple_reflections()
            count = 0
            orbit = [self]
            orbit_set = set(orbit);
            while count < len(orbit):
                w = orbit[count]
                count += 1
                for s in gens:
                    w_new = s*w*s**-1
                    if w_new not in orbit_set:
                        orbit.append(w_new)
                        orbit_set.add(w_new)
            orbit_set = frozenset(orbit_set)
            W._conjugacy_classes[self] = orbit_set
            return orbit_set

        @cached_in_parent_method
        def reduced_word(self):
            r"""
            Returns a word in the simple reflections to obtain ``self``

            EXAMPLES::

                sage: W = ComplexReflectionGroup((5,1,1),index_set=['a'],hyperplane_index_set=['x'],reflection_index_set=['A','B','C','D'])

                sage: [ w.reduced_word() for w in W ]
                [word: , word: a, word: aa, word: aaa, word: aaaa]
            
            .. seealso:: :meth:`reduced_word_in_reflections`
            """
            W = self.parent()
            if self._reduced_word is None:
                inv_dict = dict( (W._index_set[i],i) for i in W._index_set.keys() )                    
                gens = [ W.simple_reflection(i) for i in W.index_set() ]
                word = gap_factorization(self,gens,inv_dict)
                self._reduced_word = Word(word)
            return self._reduced_word

        @cached_in_parent_method
        def reduced_word_in_reflections(self):
            r"""
            Returns a word in the reflections to obtain ``self``

            EXAMPLES::

                sage: W = ComplexReflectionGroup((5,1,1),index_set=['a'],reflection_index_set=['A','B','C','D'])

                sage: [ w.reduced_word_in_reflections() for w in W ]
                [word: , word: A, word: B, word: C, word: D]
            
            .. seealso:: :meth:`reduced_word`
            """
            W = self.parent()
            if self == W.one():
                return Word([])
            if W.is_real():
                R = W.reflections()
                r = self.reflection_length()
                for i in range(len(R)):
                    t = R[i]
                    w = t*self
                    if w.reflection_length() < r:
                        return Word([i]) + w.reduced_word_in_reflections()
            else: 
                inv_dict = dict( (W._reflection_index_set[i],i) for i in W._reflection_index_set.keys() )
                gens = [ W.reflection(i) for i in W.reflection_index_set() ]
                return Word(gap_factorization(self,gens,inv_dict))

        @cached_in_parent_method
        def length(self):
            r"""
            Returns the length of ``self`` in generating reflections. This is
            the minimal numbers of generating reflections needed to obtain ``self``.

            EXAMPLES::

                tba
            """
            return len( self.reduced_word() )

        @cached_in_parent_method
        def reflection_length(self, in_unitary_group=False):
            r"""
            Returns the reflection length of ``self``. This is
            the minimal numbers of reflections needed to obtain ``self``.

            INPUT:

            - in_unitary_group -- (default:False) if True, the reflection length is computed in the unitary group which is the dimension of the move space of ``self``.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: sorted([ t.reflection_length() for t in W ])
                [0, 1, 1, 1, 2, 2]

                sage: W = ComplexReflectionGroup((2,1,2))
                sage: sorted([ t.reflection_length() for t in W ])
                [0, 1, 1, 1, 1, 2, 2, 2]

                sage: W = ComplexReflectionGroup((2,2,2))
                sage: sorted([ t.reflection_length() for t in W ])
                [0, 1, 1, 2]

                sage: W = ComplexReflectionGroup((3,1,2))
                sage: sorted([ t.reflection_length() for t in W ])
                [0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            """
            W = self.parent()
            if self in W.conjugacy_classes_representatives():
                if in_unitary_group or W.is_real():
                    return W.rank()-self.reflection_eigenvalues(test_class_repr=False).count(0)
                else:
                    return len(self.reduced_word_in_reflections())
            else:
                w = self.conjugacy_class_representative()
                assert w in self.parent().conjugacy_classes_representatives()
                return w.reflection_length(in_unitary_group=in_unitary_group)

        @cached_in_parent_method
        def conjugacy_class_representative(self):
            conj = self.parent().conjugacy_classes()
            for w in conj.keys():
                if self in conj[w]:
                    return w

        @cached_in_parent_method
        def right_coset_representatives(self):
            r"""
            Returns the right coset representatives of ``self``.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.reduced_word(), [ v.reduced_word() for v in w.right_coset_representatives() ]
                 [word: , word: 1, word: 0, word: 10, word: 01, word: 010]
                0 [word: , word: 1, word: 10]
                1 [word: , word: 0, word: 01]
                01 [word: ]
                10 [word: ]
                010 [word: , word: 1, word: 0]
            """
            W = self.parent()
            T = W.reflections()
            T_fix = [ i+1 for i in T.keys() if self.fix_space().is_subspace(T[i].fix_space()) ]
            S = str(gap3('ReducedRightCosetRepresentatives(%s,ReflectionSubgroup(%s,%s))'%(W._gap_group._name,W._gap_group._name,T_fix)))
            exec('L = '+gap_return(S,coerce_obj='W'))
            return L

        def left_coset_representatives(self):
            r"""
            Returns the left coset representatives of ``self``.

            .. seealso:: :meth:`right_coset_representatives`

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.reduced_word(), [ v.reduced_word() for v in w.left_coset_representatives() ]
                 [word: , word: 1, word: 0, word: 01, word: 10, word: 010]
                0 [word: , word: 1, word: 01]
                1 [word: , word: 0, word: 10]
                01 [word: ]
                10 [word: ]
                010 [word: , word: 1, word: 0]
            """
            return [ w**-1 for w in self.right_coset_representatives() ]

        @cached_in_parent_method
        def as_matrix(self):
            r"""
            Returns ``self`` as a matrix acting on the underlying vector space.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W:
                ...       print w.reduced_word()
                ...       print w.as_matrix()
                [1 0]
                [0 1]
                0
                [-1  0]
                [ 1  1]
                1
                [ 1  1]
                [ 0 -1]
                01
                [-1 -1]
                [ 1  0]
                10
                [ 0  1]
                [-1 -1]
                010
                [ 0 -1]
                [-1  0]
            """
            W = self.parent()
            if W._reflection_representation is None:
                Delta = W.simple_roots()
                Phi = W.roots()
                M = Matrix([ Phi[self(Phi.index(alpha)+1)-1] for alpha in Delta ])
                return W.base_change_matrix() * M
            else:
                refl_repr = W._reflection_representation
                id_mat = identity_matrix(QQ,refl_repr[W.index_set()[0]].nrows())
                return prod( [refl_repr[i] for i in self.reduced_word()],id_mat  )

        @cached_in_parent_method
        def fix_space(self):
            r"""
            Returns the fix space of ``self``. This is the sub vector space of the
            underlying vector space on which ``self`` acts trivially.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W:
                ...       print w.reduced_word()
                ...       print w.fix_space()                
                Vector space of degree 2 and dimension 2 over Universal Cyclotomic Field endowed with the Zumbroich basis
                Basis matrix:
                [1 0]
                [0 1]
                0
                Vector space of degree 2 and dimension 1 over Universal Cyclotomic Field endowed with the Zumbroich basis
                Basis matrix:
                [0 1]
                1
                Vector space of degree 2 and dimension 1 over Universal Cyclotomic Field endowed with the Zumbroich basis
                Basis matrix:
                [1 0]
                01
                Vector space of degree 2 and dimension 0 over Universal Cyclotomic Field endowed with the Zumbroich basis
                Basis matrix:
                []
                10
                Vector space of degree 2 and dimension 0 over Universal Cyclotomic Field endowed with the Zumbroich basis
                Basis matrix:
                []
                010
                Vector space of degree 2 and dimension 1 over Universal Cyclotomic Field endowed with the Zumbroich basis
                Basis matrix:
                [ 1 -1]
            """
            return (self.as_matrix()-identity_matrix(UCF,self.parent().rank())).right_kernel()

        @cached_in_parent_method
        def reflection_eigenvalues(self,test_class_repr=True):
            r"""
            Returns the reflection eigenvalues of ``self``.
            """
            return self.parent().reflection_eigenvalues(self,test_class_repr=test_class_repr)

        @cached_in_parent_method
        def galois_conjugates(self):
            r"""
            Returns all Galois conjugates of ``self``.

            EXAMPLES::

                tba
            """
            rk = self.parent().rank()
            M = self.as_matrix()
            L = [ UCF(x) for x in M.list() ]
            m = lcm([ x.field_order() for x in L ])
            L_gals = [ x.galois_conjugates(m) for x in L ]
            conjugates = []
            for i in range(len(L_gals[0])):
                conjugates.append( Matrix(rk, [ X[i] for X in L_gals ] ) )
            return conjugates

        def __cmp__(self, other):
            r"""
            Without this comparison method, the initialization of this
            permutation group fails ...
            """
            return super(FiniteComplexReflectionGroup.Element, self).__cmp__(other)

class IrreducibleFiniteComplexReflectionGroup(FiniteComplexReflectionGroup):

    def _repr_(self):
        r"""
        Returns the string representation of ``self``.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3)); W
            Irreducible finite complex reflection group of rank 2 and type A2
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible finite complex reflection group of rank %s and type %s'%(self._rank,type_str)

    @cached_method
    def a_coxeter_element(self):
        r"""
        Returns a Coxeter element of a well-generated, irreducible reflection group. This is an element
        having a regular eigenvector (a vector not contained in any recflecting hyperplane of ``self``).

        REMARK:

        - ``self`` is assumed to be well-generated.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: W.a_coxeter_element().reduced_word()
            word: 012

            sage: W = ComplexReflectionGroup((2,1,4))
            sage: W.a_coxeter_element().reduced_word()
            word: 0123

            sage: W = ComplexReflectionGroup((4,1,4))
            sage: W.a_coxeter_element().reduced_word()
            word: 0123

            sage: W = ComplexReflectionGroup((4,4,4))
            sage: W.a_coxeter_element().reduced_word()
            word: 0123
        """
        assert self.is_irreducible()
        assert self.is_well_generated()
        inverse_index = dict([(self._index_set[i],i) for i in self._index_set.keys()])
        return self.from_word( inverse_index[i] for i in sorted(self._index_set.values()) )

    @cached_method
    def coxeter_elements(self):
        r"""
        Returns the (unique) conjugacy class in ``self`` containing all Coxeter elements.

        REMARK:

        - ``self`` is assumed to be well-generated.
        - This works even beyond real reflection groups, but the conjugacy class is not unique and we only obtain one such class.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))
            sage: [ c.reduced_word() for c in W.coxeter_elements() ]
            [word: 10, word: 01]

            sage: W = ComplexReflectionGroup((1,1,4))
            sage: [ c.reduced_word() for c in W.coxeter_elements() ]
            [word: 01201, word: 12010, word: 210, word: 012, word: 201, word: 120]
        """
        return self.a_coxeter_element().conjugacy_class()

    @cached_method
    def standard_coxeter_elements(self):
        r"""
        Returns all standard Coxeter elements in ``self``. This is the set of all
        elements in self obtained from any product of the simple reflections in ``self``.

        REMARK:

        - ``self`` is assumed to be well-generated.
        - This works even beyond real reflection groups, but the conjugacy class is not unique and we only obtain one such class.

        EXAMPLES::

            tba
        """
        assert self.is_irreducible()
        assert self.is_well_generated()
        from sage.combinat.permutation import Permutations
        return set( self.from_word(w) for w in Permutations(self.index_set()) )

    def elements_below_coxeter_element(self, c=None):
        r"""
        Returns all elements in ``self`` in the interval `[1,c]` in the absolute order of ``self``.
        This order is defines by `\omega \leq_R \tau \Leftrightarrow \ell_R(\omega) + \ell_R(\omega^{-1} \tau) = \ell_R(\tau)``.

        REMARK:

        - ``self`` is assumed to be well-generated.

        INPUT:

        - c -- (default:None) if an element ``c`` is given, it is used as the maximal element in the interval. If a list is given,
                the union of the various maximal elements is computed.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))

            sage: [ w.reduced_word() for w in W.elements_below_coxeter_element() ]
            [word: , word: 0, word: 1, word: 01, word: 010]

            sage: [ w.reduced_word() for w in W.elements_below_coxeter_element(W.from_word([1,0])) ]
            [word: , word: 0, word: 1, word: 10, word: 010]

            sage: [ w.reduced_word() for w in W.elements_below_coxeter_element(W.from_word([1])) ]
            [word: , word: 1]
        """
        if c in self:
            cs = [c]
        elif c is None:
            cs = [self.a_coxeter_element()]
        else:
            cs = list(c)
        l = cs[0].reflection_length(in_unitary_group=True)
        R = self.reflections()
        f = lambda pi: any( pi.reflection_length(in_unitary_group=True) + (c*pi**-1).reflection_length(in_unitary_group=True) == l for c in cs )
        # first computing the conjugacy classes only needed if the interaction with gap3 is slow due to a bug
        #self.conjugacy_classes()
        return filter(f,self)
        #return set( self.magma_closure_iter(I=R, predicate=lambda pi: pi.reflection_length(in_unitary_group=True) + (c*pi**-1).reflection_length(in_unitary_group=True) == l) )

    @cached_method
    def noncrossing_partition_lattice(self, c=None, L=None):
        r"""
        Returns the the interval `[1,c]` in the absolute order of ``self`` as a finite lattice.

        .. seealso:: :meth:`elements_below_coxeter_element`

        REMARK:

        - ``self`` is assumed to be well-generated.

        INPUT:

        - c -- (default:None) if an element ``c`` in ``self`` is given, it is used as the maximal element in the interval.

        - L -- (default:None) if a subset ``L`` (must be hashable!) of ``self`` is given, it is used as the underlying set. Only cover relations are checked though.

        EXAMPLES::

            sage: W = ComplexReflectionGroup((1,1,3))

            sage: [ w.reduced_word() for w in W.noncrossing_partition_lattice() ]
            [word: , word: 0, word: 1, word: 010, word: 01]

            sage: [ w.reduced_word() for w in W.noncrossing_partition_lattice(W.from_word([1,0])) ]
            [word: , word: 0, word: 1, word: 010, word: 10]

            sage: [ w.reduced_word() for w in W.noncrossing_partition_lattice(W.from_word([1])) ]
            [word: , word: 1]
        """
        from sage.combinat.posets.all import Poset, LatticePoset
        R = self.reflections()
        if L is None:
            L = self.elements_below_coxeter_element(c=c)
        rels = []
        for pi in L:
            for t in R:
                tau = pi*t
                if tau in L and pi.reflection_length(in_unitary_group=True) + 1 == tau.reflection_length(in_unitary_group=True):
                    rels.append([pi,tau])
        P = Poset(([],rels), cover_relations=True, facade=True)
        if P.is_lattice():
            return LatticePoset(P)
        else:
            return P

    def generalized_noncrossing_partitions(self, m, c=None, positive=False):
        from sage.combinat.combination import Combinations
        NC = self.noncrossing_partition_lattice(c=c)
        one = self.one()
        if c is None:
            c = self.a_coxeter_element()
        chains = NC.chains()
        NCm = set()
        iter = chains.breadth_first_search_iterator()
        chain = iter.next()
        chain = iter.next()
        while len(chain) <= m:
            chain.append( c )
            for i in range(len(chain)-1,0,-1):
                chain[i] = chain[i-1]**-1 * chain[i]
            k = m+1 - len(chain)
            for positions in Combinations(range(m+1),k):
                ncm = []
                for l in range(m+1):
                    if l in positions:
                        ncm.append( one )
                    else:
                        l_prime = l - len( [ i for i in positions if i <= l ] )
                        ncm.append( chain[l_prime] )
                if not positive or prod(ncm[:-1]).has_full_support():
                    NCm.add(tuple(ncm))
            try:
                chain = iter.next()
            except StopIteration:
                chain = range(m+1)
        return NCm

    @cached_method
    def absolute_poset(self):
        r"""
        Returns the poset induced by the absolute order of ``self`` as a finite lattice.

        .. seealso:: :meth:`noncrossing_partition_lattice`

        EXAMPLES::

            sage: P = ComplexReflectionGroup((1,1,3)).absolute_poset(); P
            Finite poset containing 6 elements

            sage: [ w.reduced_word() for w in P ]
            [word: , word: 0, word: 1, word: 010, word: 01, word: 10]
        """
        return self.noncrossing_partition_lattice(L=self)

    class Element(FiniteComplexReflectionGroup.Element):

        @cached_in_parent_method
        def is_coxeter_element(self,which_primitive=1,test_class_repr=True):
            r"""
            Returns True if ``self`` is a Coxeter element.

            .. seealso:: :meth:`a_coxeter_element`

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.reduced_word(), w.is_coxeter_element()
                 False
                0 False
                1 False
                01 True
                10 True
                010 False
            """
            assert self.parent().is_well_generated()
            h = self.parent().coxeter_number()
            # to check regularity for a Coxeter number h, we get that an eigenvector is regular for free
            return any( QQ(ev).denom() == h and QQ(ev).numer() == which_primitive for ev in self.reflection_eigenvalues(test_class_repr=test_class_repr) )

        @cached_in_parent_method
        def is_h_regular(self,test_class_repr=True):
            r"""
            Returns True if self is regular. I.e., self has an eigenvector
            with eigenvalue `h` and which does not lie in any reflecting hyperplane.
            Here, `h` denotes the *Coxeter number* of ``self.parent()``.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.reduced_word(), w.is_h_regular()
                 False
                0 False
                1 False
                01 True
                10 True
                010 False
            """
            assert self.parent().is_well_generated()
            h = self.parent().coxeter_number()
            # to check regularity for a Coxeter number h, we get that an eigenvector is regular for free
            return any( QQ(ev).denom() == h for ev in self.reflection_eigenvalues(test_class_repr=test_class_repr) )

        @cached_in_parent_method
        def is_regular(self,h,test_class_repr=True):
            r"""
            Returns True if self is regular. I.e., self has an eigenvector
            with eigenvalue `h` and which does not lie in any reflecting hyperplane.
            Here, `h` denotes the *Coxeter number* of ``self.parent()``.

            EXAMPLES::

                sage: W = ComplexReflectionGroup((1,1,3))
                sage: for w in W: print w.reduced_word(), w.is_regular(W.coxeter_number())
                 False
                0 False
                1 False
                01 True
                10 True
                010 False
            """
            assert self.parent().is_well_generated()
            evs = self.reflection_eigenvalues(test_class_repr=test_class_repr)
            for ev in evs:
                ev = QQ(ev)
                if h == ev.denom():
                    M = Matrix(UCF,(self.as_matrix()-E(ev.denom(),ev.numer())*identity_matrix(self.parent().rank())))
                    V = M.right_kernel()
                    if not any( V.is_subspace(H) for H in self.parent().reflecting_hyperplanes() ):
                        return True
            return False

def ComplexReflectionGroup(*args,**kwds):
    r"""
    Construct a finite (complex) reflection group as a Sage permutation group by
    fetching the permutation representation of the generators from chevie's database.

    INPUT:

    can be one of the following:

    - (a) integer between 4 and 37, which denotes an exeptional irreducible complex reflection group
    - (b) triple (r,p,n) with p divides r, which denotes the group G(r,p,n)
    - (c) list containing objects in (a) and (b)

    EXAMPLES:
    
    Finite reflection groups can be constructed from

    the complex infinite family `G(r,p,n)` with `p` divides `r`::

        sage: W = ComplexReflectionGroup((1,1,4)); W
        Irreducible finite complex reflection group of rank 3 and type A3

        sage: W = ComplexReflectionGroup((2,1,3)); W
        Irreducible finite complex reflection group of rank 3 and type B3

    Chevalley-Shepard-Todd exceptional classification types::

        sage: W = ComplexReflectionGroup(23); W
         Irreducible finite complex reflection group of rank 3 and type H3

    lists containing the above::

        sage: W = ComplexReflectionGroup((1,1,4),(2,1,3)); W
        Reducible finite complex reflection group of rank 6 and type A3 x B3
    """
    assert is_chevie_available()
    gap3.load_package("chevie")

    W_types = []
    for arg in args:
        if type(arg) is list:
            X = tuple(arg)
        else:
            X = arg
        assert is_Matrix(X) or isinstance(X,CartanMatrix) or isinstance(X,tuple) or ( X in ZZ and 4 <= X <= 37 ), "The input is not valid."
        if X == (2,2,2):
            W_types.extend([(1,1,2),(1,1,2)])
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
        cls = IrreducibleFiniteComplexReflectionGroup
    else:
        cls = FiniteComplexReflectionGroup
    return cls(tuple(W_types), index_set=kwds.get('index_set', None),
                               hyperplane_index_set=kwds.get('hyperplane_index_set', None),
                               reflection_index_set=kwds.get('reflection_index_set', None) )

def gap_factorization(w,gens,inv_dict):
    gap3.execute('W := GroupWithGenerators(%s)'%str(gens))
    gap3.execute(gap_factorization_code)
    fac = gap3('MinimalWord(W,%s)'%str(w)).sage()
    return [ inv_dict[i-1] for i in fac ]

gap_factorization_code = '# MinimalWord(G,w) \n \
# given a permutation group G find some expression of minimal length in the \n \
# generators of G and their inverses of the element w (an inverse is \n \
# representated by a negative index). \n \
# To speed up  later calls to  the same function  the fields G.base, G.words, \n \
# G.nbwordslength are kept. \n \
MinimalWord:=function(G,w) \n \
  local decode,i,p,g,h,n,bag,nbe,nbf,new,gens,inds; \n \
# to save space elements of G are represented as image of the base, and \n \
# words are represented as: index of previous elt, last generator applied; \n \
  if not IsBound(G.base) then \n \
    StabChain(G);g:=G; G.base:=[]; \n \
    while IsBound(g.orbit) do Add(G.base,g.orbit[1]); g:=g.stabilizer; od; \n \
  fi; \n \
  w:=OnTuples(G.base,w); \n \
  if not IsBound(G.words) then \n \
    G.words:=[G.base]; G.lastmult:=[[0,0]]; \n \
    G.nbwordslength:=[1]; \n \
  fi; \n \
  gens:=ShallowCopy(G.generators);inds:=[1..Length(gens)]; \n \
#  for g in G.generators do \n \
#    if g<>g^-1 then Add(gens,g^-1);Add(inds,-Position(gens,g));fi; \n \
#  od; \n \
  bag:=Set(G.words); \n \
  nbe:=0;nbf:=0; \n \
  decode:=function(i)local w;w:=[]; \n \
    while i<>1 do Add(w,G.lastmult[i][2]); i:=G.lastmult[i][1];od; \n \
    return Reversed(w); \n \
  end; \n \
  while true do \n \
    if w in bag then return decode(Position(G.words,w));fi; \n \
    new:=Length(G.words); \n \
    for g in [1..Length(gens)] do \n \
      for h in [1+Sum(G.nbwordslength{[1..Length(G.nbwordslength)-1]})..new] do \n \
         n:=OnTuples(G.words[h],gens[g]); \n \
         if n in bag then \n \
           nbe:=nbe+1;# if nbe mod 500=1 then Print(".\c");fi; \n \
         else \n \
           nbf:=nbf+1;# if nbf mod 500=1 then Print("*\c");fi; \n \
       Add(G.words,n);Add(G.lastmult,[h,inds[g]]);AddSet(bag,n); \n \
         fi; \n \
       od; \n \
    od; \n \
    Add(G.nbwordslength,Length(G.words)-new); \n \
    Print("\n",G.nbwordslength[Length(G.nbwordslength)]," elements of length ", \n \
      Length(G.nbwordslength)-1); \n \
  od; \n \
end;'

def gap_return(S,coerce_obj='self'):
    S = S.replace(' ','').replace('\n','')
    S = S.replace(',(','\',check=False),%s(\'('%coerce_obj).replace('[','[%s(\''%coerce_obj).replace(']','\',check=False)]')
    return S

@cached_function
def is_chevie_available():
    r"""
    Tests whether the GAP3 Chevie package is available

    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_group import is_chevie_available
        sage: is_chevie_available() # random
        False
        sage: is_chevie_available() in [True, False]
        True
    """
    try:
        from sage.interfaces.gap3 import gap3
        gap3.load_package("chevie")
        return True
    except:
        return False
