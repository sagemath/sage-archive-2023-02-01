r"""
Type A crystal on affine factorizations
"""
#*****************************************************************************
#  Copyright (C) 2014 Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.k_tableau import WeakTableaux, WeakTableau

class AffineFactorizationCrystal(UniqueRepresentation, Parent):
    r"""
    This is an implementation of the weak crystal on affine factorizations
    with a cut-point, as introduced by Morse and Schilling.

    EXAMPLES::

        sage: from sage.combinat.crystals.affine_factorization import AffineFactorizationCrystal
        sage: W = AffineFactorizationCrystal([[3,1,1],[1]],3,3)
        sage: W.list()
        [[1, s2, s3*s2*s1],
        [1, s3*s2, s3*s1],
        [1, s3*s2*s1, s3],
        [s3, s2, s3*s1],
        [s3, s2*s1, s3],
        [s3*s2, s1, s3],
        [s3*s2*s1, 1, s3],
        [s3*s2*s1, s3, 1],
        [s3*s2, 1, s3*s1],
        [s3*s2, s3, s1],
        [s3*s2, s3*s1, 1],
        [s2, 1, s3*s2*s1],
        [s2, s3, s2*s1],
        [s2, s3*s2, s1],
        [s2, s3*s2*s1, 1]]

    We can compute the highest weight elements::

        sage: hw = [w for w in W if w.is_highest_weight()]
        sage: hw
        [[1, s2, s3*s2*s1]]
        sage: hw[0].weight()
        (3, 1, 0)

    And show that this crystal is isomorphic to the tableau model of the same weight::

        sage: C = CrystalOfTableaux(['A',2],shape=[3,1])
        sage: GC = C.digraph()
        sage: GW = W.digraph()
        sage: GC.is_isomorphic(GW, edge_labels=True)
        True

    The crystal operators themselves move elements between adjacent factors::

        sage: b = hw[0];b
        [1, s2, s3*s2*s1]
        sage: b.f(1)
        [1, s3*s2, s3*s1]

    TESTS::

        sage: TestSuite(W).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """
    @staticmethod
    def __classcall_private__(cls, skew_shape, k, n, x = None):
        """
        Classcall to mend the input.

        TESTS::

            sage: from sage.combinat.crystals.affine_factorization import AffineFactorizationCrystal
            sage: AffineFactorizationCrystal([[3,1],[1]],3,4)
            Crystal on affine factorizations of type A4 of shape [[3, 1], [1]]
        """
        skew_shape = tuple(tuple(l) for l in skew_shape)
        return super(AffineFactorizationCrystal, cls).__classcall__(cls, skew_shape, k, n, x)

    def __init__(self, skew_shape, k, n, x=None):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.affine_factorization import AffineFactorizationCrystal
            sage: W = AffineFactorizationCrystal([[3,2],[2]],3,4,0)
            sage: W.x
            0
            sage: W = AffineFactorizationCrystal([[3,2],[2]],3,4)
            sage: W.x
            1
        """
        Parent.__init__(self, category = ClassicalCrystals())
        self.n = n
        self.k = k
        self.skew_shape = skew_shape
        cartan_type = CartanType(['A',n-1])
        self._cartan_type = cartan_type
        generators = [ t.to_core_tableau().to_factorized_permutation_tableau()
                      for mu in IntegerVectors(sum(skew_shape[0])-sum(skew_shape[1]), n)
                      for t in WeakTableaux(k, skew_shape, mu, representation='bounded')]
        self.module_generators = [self(t) for t in generators]
        if x is None:
            if generators != []:
                x = min( set(range(k+1)).difference(set(
                            sum([i.reduced_word() for i in generators[0]],[]))))
            else:
                x = 0
        self.x = x

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.affine_factorization import AffineFactorizationCrystal
            sage: AffineFactorizationCrystal([[3,1],[1]],3,4)
            Crystal on affine factorizations of type A4 of shape [[3, 1], [1]]
        """
        return "Crystal on affine factorizations of type A%s of shape %s"%(self.n, [[i for i in w] for w in self.skew_shape])

    # temporary workaround while an_element is overriden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    class Element(ElementWrapper):

        def e(self, i):
            r"""
            Returns the action of `e_i` on ``self``.

            EXAMPLES::

                sage: from sage.combinat.crystals.affine_factorization import AffineFactorizationCrystal
                sage: W = AffineFactorizationCrystal(tuple([tuple([3,1]),tuple([1])]),3,4)
                sage: t = W(W.module_generators[1]); t
                [1, 1, s3, s2*s1]
                sage: t.e(1)
                [1, 1, 1, s3*s2*s1]
            """
            assert i in self.index_set()
            x = self.parent().x
            k = self.parent().k
            n = self.parent().n
            b = self.bracketing(i)
            if b[0] == []:
                return None
            a = min(b[0])
            left = [j for j in (self.value[n-i-1]).reduced_word() if j != (a+x)%(k+1)]
            right = [(j-x)%(k+1) for j in (self.value[n-i]).reduced_word()]
            m = max([j for j in range(a) if (j+x)%(k+1) not in left])
            right += [m+1]
            right.sort(reverse=True)
            right = [(j+x)%(k+1) for j in right]
            t = [self.value[j].reduced_word() for j in range(n-i-1)] + [left] + [right] + [self.value[j].reduced_word() for j in range(n-i+1,n)]
            t = WeakTableau([self.value[j].reduced_word() for j in range(n-i-1)]
                            + [left] + [right]
                            + [self.value[j].reduced_word() for j in range(n-i+1,n)],
                            self.parent().k, inner_shape = self.value._inner_shape,
                            representation='factorized_permutation')
            return self.parent()(t)

        def f(self, i):
            r"""
            Returns the action of `f_i` on ``self``.

            EXAMPLES::

                sage: from sage.combinat.crystals.affine_factorization import AffineFactorizationCrystal
                sage: W = AffineFactorizationCrystal(tuple([tuple([3,1]),tuple([1])]),3,4)
                sage: t = W(W.module_generators[1]); t
                [1, 1, s3, s2*s1]
                sage: t.f(2)
                [1, s3, 1, s2*s1]
                sage: t.f(1)
                [1, 1, s3*s2, s1]
            """
            assert i in self.index_set()
            x = self.parent().x
            k = self.parent().k
            n = self.parent().n
            b = self.bracketing(i)
            if b[1] == []:
                return None
            a = max(b[1])
            right = [j for j in (self.value[n-i]).reduced_word() if j != (a+x)%(k+1)]
            left = [(j-x)%(k+1) for j in (self.value[n-i-1]).reduced_word()]
            m = min([j for j in range(a+1,k+2) if (j+x)%(k+1) not in right])
            left += [m-1]
            left.sort(reverse=True)
            left = [(j+x)%(k+1) for j in left]
            t = [self.value[j].reduced_word() for j in range(n-i-1)] + [left] + [right] + [self.value[j].reduced_word() for j in range(n-i+1,n)]
            t = WeakTableau([self.value[j].reduced_word() for j in range(n-i-1)]
                            + [left] + [right]
                            + [self.value[j].reduced_word() for j in range(n-i+1,n)],
                            self.parent().k, inner_shape = self.value._inner_shape,
                            representation='factorized_permutation')
            return self.parent()(t)

        def bracketing(self, i):
            r"""
            Removes all bracketed letters between `i`-th and `i+1`-th entry.

            EXAMPLES::

                sage: from sage.combinat.crystals.affine_factorization import AffineFactorizationCrystal
                sage: W = AffineFactorizationCrystal(tuple([tuple([3,1]),tuple([1])]),3,4)
                sage: t = W(W.module_generators[1])
                sage: t.bracketing(1)
                [[3], [2, 1]]
            """
            n = self.parent().n
            x = self.parent().x
            k = self.parent().k
            right = (self.value[n-i]).reduced_word()
            left = (self.value[n-i-1]).reduced_word()
            right_n = [(j-x)%(k+1) for j in right]
            left_n = [(j-x)%(k+1) for j in left]
            left_unbracketed = []
            while left_n != []:
                m = max(left_n)
                left_n.remove(m)
                l = [j for j in right_n if j>m]
                if l != []:
                    right_n.remove(min(l))
                else:
                    left_unbracketed += [m]
            return [[j for j in left_unbracketed],[j for j in right_n]]
