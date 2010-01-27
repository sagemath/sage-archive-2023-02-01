r"""
Affine crystals
"""

#*****************************************************************************
#       Copyright (C) 2008 Brant Jones <brant at math.ucdavis.edu>
#                          Anne Schilling <anne at math.ucdavis.edu>
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
# Acknowledgment: most of the design and implementation of this
# library is heavily inspired from MuPAD-Combinat.
#****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.combinat.crystals.crystals import Crystal, CrystalElement
from sage.combinat.root_system.cartan_type import CartanType

class AffineCrystal(Crystal):
    r"""
    The abstract class of affine crystals
    """

class AffineCrystalFromClassical(Crystal, UniqueRepresentation):
    r"""
    This abstract class can be used for affine crystals that are constructed from a classical crystal.
    The zero arrows can be implemented using different methods (for example using a Dynkin diagram
    automorphisms or virtual crystals).

    For general information about crystals see :mod:`sage.combinat.crystals`.

    INPUT:

    - ``cartan_type`` -  The Cartan type of the resulting affine crystal

    - ``classical_crystal`` - instance of a classical crystal.

    EXAMPLES::

        sage: n=2
        sage: C=CrystalOfTableaux(['A',n],shape=[1])
        sage: pr=lambda x : C(x.to_tableau().promotion(n))
        sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
        sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
        sage: A.list()
        [[[1]], [[2]], [[3]]]
        sage: A.cartan_type()
        ['A', 2, 1]
        sage: A.index_set()
        [0, 1, 2]
        sage: b=A(rows=[[1]])
        sage: b.weight()
        -Lambda[0] + Lambda[1]
        sage: b.classical_weight()
        (1, 0, 0)
        sage: [x.s(0) for x in A.list()]
        [[[3]], [[2]], [[1]]]
        sage: [x.s(1) for x in A.list()]
        [[[2]], [[1]], [[3]]]
    """

    @staticmethod
    def __classcall__(cls, cartan_type, *args, **options):
        ct = CartanType(cartan_type)
        return super(AffineCrystalFromClassical, cls).__classcall__(cls, ct, *args, **options)

    def __init__(self, cartan_type, classical_crystal, category = None):
        """
        Input is an affine Cartan type 'cartan_type', a classical crystal 'classical_crystal', and automorphism and its
        inverse 'automorphism' and 'inverse_automorphism', and the Dynkin node 'dynkin_node'

        EXAMPLES::

            sage: n=1
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: def pr(x): return C(x.to_tableau().promotion(n))
            sage: def pr_inverse(x): return C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1) # indirect doctest
            sage: A.list()
            [[[1]], [[2]]]
            sage: A.cartan_type()
            ['A', 1, 1]
            sage: A.index_set()
            [0, 1]

        TESTS::

            sage: import __main__; __main__.pr = pr; __main__.pr_inverse = pr_inverse # standard pickling hack
            sage: TestSuite(A).run()
        """
        if category is None:
            category = FiniteEnumeratedSets()
        self._cartan_type = cartan_type
        super(AffineCrystalFromClassical, self).__init__(category = category)
        self.classical_crystal = classical_crystal;
        self.module_generators = map( self.retract, self.classical_crystal.module_generators )
	self.element_class._latex_ = lambda x: x.lift()._latex_()

    def _repr_(self):
        """
        EXAMPLES::

            sage: n=1
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: def pr(x): return C(x.to_tableau().promotion(n))
            sage: def pr_inverse(x): return C(x.to_tableau().promotion_inverse(n))
            sage: AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1) # indirect doctest
            An affine crystal for type ['A', 1, 1]
        """
        return "An affine crystal for type %s"%self.cartan_type()


    def __iter__(self):
        r"""
        Construct the iterator from the underlying classical crystal.

        TESTS::

            sage: n=1
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1) # indirect doctest
            sage: A.list() # indirect doctest
            [[[1]], [[2]]]
        """
        for x in self.classical_crystal:
            yield self.retract(x)

    # should be removed once crystal defines __iter__ instead of list
    def list(self):
        """
        Returns the list of all crystal elements using the underlying classical crystal

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: A.list()
            [[[1]], [[2]], [[3]]]
        """
        return map( self.retract, self.classical_crystal.list() )

    def lift(self, affine_elt):
        """
        Lifts an affine crystal element to the corresponding classical crystal element

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A.list()[0]
            sage: A.lift(b)
            [[1]]
            sage: A.lift(b).parent()
            The crystal of tableaux of type ['A', 2] and shape(s) [[1]]
        """
        return affine_elt.lift()

    def retract(self, classical_elt):
        """
        Transforms a classical crystal element to the corresponding affine crystal element

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: t=C(rows=[[1]])
            sage: t.parent()
            The crystal of tableaux of type ['A', 2] and shape(s) [[1]]
            sage: A.retract(t)
            [[1]]
            sage: A.retract(t).parent() is A
            True
        """
        return self.element_class(classical_elt, parent = self)

    def __call__(self, *value, **options):
        r"""
        Coerces value into self.

        EXAMPLES:

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: def pr(x): return C(x.to_tableau().promotion(n))
            sage: def pr_inverse(x): return C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A(rows=[[1]])
            sage: b
            [[1]]
            sage: b.parent()
            An affine crystal for type ['A', 2, 1]
            sage: A(b) is b
            True
        """
        if len(value) == 1 and isinstance(value[0], self.element_class) and value[0].parent() == self:
            return value[0]
        else: # Should do sanity checks!  (Including check for inconsistent parent.)
            return self.retract(self.classical_crystal(*value, **options))

    def __contains__(self, x):
        r"""
        Checks whether x is an element of self.

        EXAMPLES:

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A(rows=[[1]])
            sage: A.__contains__(b)
            True
        """
        return x.parent() is self


class AffineCrystalFromClassicalElement(ElementWrapper, CrystalElement):
    r"""
    Elements of crystals that are constructed from a classical crystal.
    The elements inherit many of their methods from the classical crystal
    using lift and retract.

    This class is not instantiated directly but rather __call__ed from
    AffineCrystalFromClassical.  The syntax of this is governed by the
    (classical) CrystalOfTableaux.

    EXAMPLES::

        sage: n=2
        sage: C=CrystalOfTableaux(['A',n],shape=[1])
        sage: pr=lambda x : C(x.to_tableau().promotion(n))
        sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
        sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
        sage: b=A(rows=[[1]])
        sage: b.__repr__()
        '[[1]]'
    """

    def classical_weight(self):
        """
        Returns the classical weight corresponding to self.

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A(rows=[[1]])
            sage: b.classical_weight()
            (1, 0, 0)
        """
        return self.lift().weight()

    def lift(self):
        """
        Lifts an affine crystal element to the corresponding classical crystal element

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A.list()[0]
            sage: b.lift()
            [[1]]
            sage: b.lift().parent()
            The crystal of tableaux of type ['A', 2] and shape(s) [[1]]
        """
        return self.value

    def pp(self):
        """
        Method for pretty printing

        EXAMPLES::

            sage: K = KirillovReshetikhinCrystal(['D',3,2],1,1)
            sage: t=K(rows=[[1]])
            sage: t.pp()
            1
        """
        return self.lift().pp()

    @abstract_method
    def e0(self):
        r"""
        Assumes that `e_0` is implemented separately.
        """

    @abstract_method
    def f0(self):
        r"""
        Assumes that `f_0` is implemented separately.
        """

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A(rows=[[1]])
            sage: b.e(0)
            [[3]]
            sage: b.e(1)
        """
        if i == 0:
            return self.e0()
        else:
            x = self.lift().e(i)
            if (x == None):
                return None
            else:
                return self.parent().retract(x)

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A(rows=[[3]])
            sage: b.f(0)
            [[1]]
            sage: b.f(2)
        """
        if i == 0:
            return self.f0()
        else:
            x = self.lift().f(i)
            if (x == None):
                return None
            else:
                return self.parent().retract(x)

    def epsilon0(self):
        r"""
        Uses `epsilon_0` from the super class, but should be implemented if a faster implementation exists.
        """
        return super(AffineCrystalFromClassicalElement, self).epsilon(0)

    def epsilon(self, i):
        """
        Returns the maximal time the crystal operator `e_i` can be applied to self.

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: [x.epsilon(0) for x in A.list()]
            [1, 0, 0]
            sage: [x.epsilon(1) for x in A.list()]
            [0, 1, 0]
        """
        if i == 0:
            return self.epsilon0()
        else:
            return self.lift().epsilon(i)

    def phi0(self):
        r"""
        Uses `phi_0` from the super class, but should be implemented if a faster implementation exists.
        """
        return super(AffineCrystalFromClassicalElement, self).phi(0)

    def phi(self, i):
        r"""
        Returns the maximal time the crystal operator `f_i` can be applied to self.

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: [x.phi(0) for x in A.list()]
            [0, 0, 1]
            sage: [x.phi(1) for x in A.list()]
            [1, 0, 0]
        """
        if i == 0:
            return self.phi0()
        else:
            return self.lift().phi(i)

AffineCrystalFromClassical.Element = AffineCrystalFromClassicalElement


class AffineCrystalFromClassicalAndPromotion(AffineCrystalFromClassical):
    r"""
    Crystals that are constructed from a classical crystal and a
    Dynkin diagram automorphism `\sigma`.  In type `A_n`, the Dynkin
    diagram automorphism is `i \to i+1 \pmod n+1` and the
    corresponding map on the crystal is the promotion operation
    `\mathrm{pr}` on tableaux. The affine crystal operators are given
    by `f_0= \mathrm{pr}^{-1} f_{\sigma(0)} \mathrm{pr}`.

    INPUT:

    - ``cartan_type`` -  The Cartan type of the resulting affine crystal

    - ``classical_crystal`` - instance of a classical crystal.

    - ``automorphism, inverse_automorphism`` - A function on the
      elements of the classical_crystal

    - ``dynkin_node`` - Integer specifying the classical node in the
      image of the zero node under the automorphism sigma.

    EXAMPLES::

        sage: n=2
        sage: C=CrystalOfTableaux(['A',n],shape=[1])
        sage: pr=lambda x : C(x.to_tableau().promotion(n))
        sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
        sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
        sage: A.list()
        [[[1]], [[2]], [[3]]]
        sage: A.cartan_type()
        ['A', 2, 1]
        sage: A.index_set()
        [0, 1, 2]
        sage: b=A(rows=[[1]])
        sage: b.weight()
        -Lambda[0] + Lambda[1]
        sage: b.classical_weight()
        (1, 0, 0)
        sage: [x.s(0) for x in A.list()]
        [[[3]], [[2]], [[1]]]
        sage: [x.s(1) for x in A.list()]
        [[[2]], [[1]], [[3]]]
    """

    def __init__(self, cartan_type, classical_crystal, p_automorphism, p_inverse_automorphism, dynkin_node):
        """
        Input is an affine Cartan type 'cartan_type', a classical crystal 'classical_crystal', and automorphism and its
        inverse 'automorphism' and 'inverse_automorphism', and the Dynkin node 'dynkin_node'

        EXAMPLES::

            sage: n=1
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: def pr(x): return C(x.to_tableau().promotion(n))
            sage: def pr_inverse(x): return C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: A.list()
            [[[1]], [[2]]]
            sage: A.cartan_type()
            ['A', 1, 1]
            sage: A.index_set()
            [0, 1]

        TESTS::

            sage: import __main__; __main__.pr = pr; __main__.pr_inverse = pr_inverse # standard pickling hack
            sage: TestSuite(A).run()
        """
        AffineCrystalFromClassical.__init__(self, cartan_type, classical_crystal)
        self.p_automorphism = p_automorphism
        self.p_inverse_automorphism = p_inverse_automorphism
        self.dynkin_node = dynkin_node

    def automorphism(self, x):
        """
        Gives the analogue of the affine Dynkin diagram automorphism on the level of crystals

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A.list()[0]
            sage: A.automorphism(b)
            [[2]]
        """
        return self.retract( self.p_automorphism( x.lift() ) )

    def inverse_automorphism(self, x):
        """
        Gives the analogue of the inverse of the affine Dynkin diagram automorphism on the level of crystals

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A.list()[0]
            sage: A.inverse_automorphism(b)
            [[3]]
        """
        return self.retract( self.p_inverse_automorphism( x.lift() ) )

class AffineCrystalFromClassicalAndPromotionElement(AffineCrystalFromClassicalElement):
    r"""
    Elements of crystals that are constructed from a classical crystal
    and a Dynkin diagram automorphism.  In type A, the automorphism is
    the promotion operation on tableaux.

    This class is not instantiated directly but rather __call__ed from
    AffineCrystalFromClassicalAndPromotion.  The syntax of this is governed by the
    (classical) CrystalOfTableaux.

    Since this class inherits from AffineClassicalFromClassicalElement, the methods
    that need to be implemented are e0, f0 and possibly epsilon0 and phi0 if more efficient
    algorithms exist.

    EXAMPLES::

        sage: n=2
        sage: C=CrystalOfTableaux(['A',n],shape=[1])
        sage: pr=lambda x : C(x.to_tableau().promotion(n))
        sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
        sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
        sage: b=A(rows=[[1]])
        sage: b.__repr__()
        '[[1]]'
    """

    def e0(self):
        r"""
        Implements `e_0` using the automorphism as
        `e_0 = \pr^{-1} e_{dynkin_node} \pr`

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A(rows=[[1]])
            sage: b.e0()
            [[3]]
        """
        x = self.parent().automorphism(self).e(self.parent().dynkin_node)
        if (x == None):
            return None
        else:
            return self.parent().inverse_automorphism(x)

    def f0(self):
        r"""
        Implements `f_0` using the automorphism as
        `f_0 = \pr^{-1} f_{dynkin_node} \pr`

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: b=A(rows=[[3]])
            sage: b.f0()
            [[1]]
        """
        x = self.parent().automorphism(self).f(self.parent().dynkin_node)
        if (x == None):
            return None
        else:
            return self.parent().inverse_automorphism(x)

    def epsilon0(self):
        r"""
        Implements `epsilon_0` using the automorphism.

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: [x.epsilon0() for x in A.list()]
            [1, 0, 0]
        """
        x = self.parent().automorphism(self)
        return x.lift().epsilon(self.parent().dynkin_node)

    def phi0(self):
        r"""
        Implements `phi_0` using the automorphism.

        EXAMPLES::

            sage: n=2
            sage: C=CrystalOfTableaux(['A',n],shape=[1])
            sage: pr=lambda x : C(x.to_tableau().promotion(n))
            sage: pr_inverse = lambda x : C(x.to_tableau().promotion_inverse(n))
            sage: A=AffineCrystalFromClassicalAndPromotion(['A',n,1],C,pr,pr_inverse,1)
            sage: [x.phi0() for x in A.list()]
            [0, 0, 1]
        """
        x = self.parent().automorphism(self)
        return x.lift().phi(self.parent().dynkin_node)

AffineCrystalFromClassicalAndPromotion.Element = AffineCrystalFromClassicalAndPromotionElement
