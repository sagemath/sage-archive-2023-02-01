"""
Weyl Character Rings
"""
#*****************************************************************************
#  Copyright (C) 2011 Daniel Bump <bump at match.stanford.edu>
#                     Nicolas Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.all import Category, Algebras, AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.dynkin_diagram import DynkinDiagram
from sage.combinat.root_system.root_system import RootSystem
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.functional import is_even, is_odd
from sage.modules.free_module_element import vector
from sage.rings.all import ZZ, QQ

class WeylCharacterRing(CombinatorialFreeModule):
    """
    A class for rings of Weyl characters.

    Let `K` be a compact Lie group, which we assume is semisimple and
    simply-connected. Its complexified Lie algebra `L` is the Lie algebra of a
    complex analytic Lie group `G`. The following three categories are
    equivalent: finite-dimensional representations of `K`; finite-dimensional
    representations of `L`; and finite-dimensional analytic representations of
    `G`. In every case, there is a parametrization of the irreducible
    representations by their highest weight vectors. For this theory of Weyl,
    see (for example):

    * Adams, *Lectures on Lie groups*
    * Broecker and Tom Dieck, *Representations of Compact Lie groups*
    * Bump, *Lie Groups*
    * Fulton and Harris, *Representation Theory*
    * Goodman and Wallach, *Representations and Invariants of the Classical Groups*
    * Hall, *Lie Groups, Lie Algebras and Representations*
    * Humphreys, *Introduction to Lie Algebras and their representations*
    * Procesi, *Lie Groups*
    * Samelson, *Notes on Lie Algebras*
    * Varadarajan, *Lie groups, Lie algebras, and their representations*
    * Zhelobenko, *Compact Lie Groups and their Representations*.

    Computations that you can do with these include computing their
    weight multiplicities, products (thus decomposing the tensor
    product of a representation into irreducibles) and branching
    rules (restriction to a smaller group).

    There is associated with `K`, `L` or `G` as above a lattice, the weight
    lattice, whose elements (called weights) are characters of a Cartan
    subgroup or subalgebra. There is an action of the Weyl group `W` on
    the lattice, and elements of a fixed fundamental domain for `W`, the
    positive Weyl chamber, are called dominant. There is for each
    representation a unique highest dominant weight that occurs with
    nonzero multiplicity with respect to a certain partial order, and
    it is called the highest weight vector.

    EXAMPLES::

        sage: L = RootSystem("A2").ambient_space()
        sage: [fw1,fw2] = L.fundamental_weights()
        sage: R = WeylCharacterRing(['A',2], prefix="R")
        sage: [R(1),R(fw1),R(fw2)]
        [R(0,0,0), R(1,0,0), R(1,1,0)]

    Here ``R(1)``, ``R(fw1)``, and ``R(fw2)`` are irreducible representations
    with highest weight vectors `0`, `\Lambda_1`, and `\Lambda_2` respecitively
    (the first two fundamental weights).

    For type `A` (also `G_2`, `F_4`, `E_6` and `E_7`) we will take as the
    weight lattice not the weight lattice of the semisimple group, but for a
    larger one. For type `A`, this means we are concerned with the
    representation theory of `K = U(n)` or `G = GL(n, \CC)` rather than `SU(n)`
    or `SU(n, \CC)`. This is useful since the representation theory of `GL(n)`
    is ubiquitous, and also since we may then represent the fundamental
    weights (in :mod:`sage.combinat.root_system.root_system`) by vectors
    with integer entries. If you are only interested in `SL(3)`, say, use
    ``WeylCharacterRing(['A',2])`` as above but be aware that ``R([a,b,c])``
    and ``R([a+1,b+1,c+1])`` represent the same character of `SL(3)` since
    ``R([1,1,1])`` is the determinant.

    For more information, see the thematic tutorial *Lie Methods and
    Related Combinatorics in Sage*, available at:

    http://www.sagemath.org/doc/thematic_tutorials/lie.html
    """
    @staticmethod
    def __classcall__(cls, ct, base_ring=ZZ, prefix=None, style="lattice"):
        """
        TESTS::

            sage: R = WeylCharacterRing("G2", style="coroots")
            sage: R.cartan_type() is CartanType("G2")
            True
            sage: R.base_ring() is ZZ
            True
        """
        ct = CartanType(ct)
        if prefix == None:
            if ct.is_atomic():
                prefix = ct[0]+str(ct[1])
            else:
                prefix = ct.__repr__()
        return super(WeylCharacterRing, cls).__classcall__(cls, ct, base_ring=base_ring, prefix=prefix, style=style)

    def __init__(self, ct, base_ring=ZZ, prefix=None, style="lattice"):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: TestSuite(A2).run()
        """
        ct = CartanType(ct)
        self._cartan_type = ct
        self._rank = ct.rank()
        self._base_ring = base_ring
        self._space = RootSystem(self._cartan_type).ambient_space()
        self._origin = self._space.zero()
        if prefix == None:
            if ct.is_atomic():
                prefix = ct[0]+str(ct[1])
            else:
                prefix = ct.__repr__()
        self._prefix = prefix
        self._style = style
        if style == "coroots":
            self._word = self._space.weyl_group().long_element().reduced_word()
        # TODO: remove the Category.join once not needed anymore (bug in CombinatorialFreeModule)
        # TODO: use GradedAlgebrasWithBasis
        category = Category.join([AlgebrasWithBasis(base_ring), Algebras(base_ring).Subobjects()])
        CombinatorialFreeModule.__init__(self, base_ring, self._space, category = category)

        # Register the embedding of self into ambient as a coercion
        self.lift.register_as_coercion()
        # Register the partial inverse as a conversion
        self.register_conversion(self.retract)

    @cached_method
    def ambient(self):
        """
        Returns the weight ring of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A2").ambient()
            The Weight ring attached to The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients
        """
        return WeightRing(self)

    # Eventually, one could want to put the cache_method here rather
    # than on _irr_weights. Or just to merge this method and _irr_weights
    def lift_on_basis(self, irr):
        """
        Expand the basis element indexed by the weight ``irr`` into the
        weight ring of ``self``.

        INPUT:

        - ``irr`` -- a dominant weight

        This is used to implement :meth:`lift`.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([2,1,0]); v
            (2, 1, 0)
            sage: A2.lift_on_basis(v)
            2*a2(1,1,1) + a2(1,2,0) + a2(1,0,2) + a2(2,1,0) + a2(2,0,1) + a2(0,1,2) + a2(0,2,1)

        This is consistent with the analoguous calculation with symmetric
        Schur functions::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s[2,1].expand(3)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
        """
        return self.ambient()._from_dict(self._irr_weights(irr))

    def demazure_character(self, hwv, word, debug=False):
        r"""
        Compute the Demazure character.

        INPUT:

        - ``hwv`` -- a (usually dominant) weight
        - ``word`` -- a Weyl group word

        Produces the Demazure character with highest weight ``hwv`` and
        ``word`` as an element of the weight ring. Only available if
        ``style="coroots"``. The Demazure operators are also available as
        methods of :class:`WeightRing` elements, and as methods of crystals.
        Given a :class:`CrystalOfTableaux` with given highest weight vector,
        the Demazure method on the crystal will give the equivalent of this
        method, except that the Demazure character of the crystal is given
        as a sum of monomials instead of an element of the :class:`WeightRing`.

        See :meth:`WeightRing.Element.demazure` and
        :meth:`sage.categories.classical_crystals.ClassicalCrystals.ParentMethods.demazure_character`

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2",style="coroots")
            sage: h=sum(A2.fundamental_weights()); h
            (2, 1, 0)
            sage: A2.demazure_character(h,word=[1,2])
            a2(0,0) + a2(-2,1) + a2(2,-1) + a2(1,1) + a2(-1,2)
            sage: A2.demazure_character((1,1),word=[1,2])
            a2(0,0) + a2(-2,1) + a2(2,-1) + a2(1,1) + a2(-1,2)
        """
        if self._style != "coroots":
            raise ValueError('demazure method unavailable. Use style="coroots".')
        hwv = self._space.from_vector_notation(hwv, style = "coroots")
        return self.ambient()._from_dict(self._demazure_weights(hwv, word=word, debug=debug))

    @lazy_attribute
    def lift(self):
        """
        The embedding of ``self`` into its weight ring.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: A2.lift
            Generic morphism:
              From: The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients
              To:   The Weight ring attached to The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients

        ::

            sage: x = -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
            sage: A2.lift(x)
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)

        As a shortcut, you may also do::

            sage: x.lift()
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)

        Or even::

            sage: a2 = WeightRing(A2)
            sage: a2(x)
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)
        """
        return self.module_morphism(self.lift_on_basis,
                                    codomain = self.ambient(),
                                    category = AlgebrasWithBasis(self.base_ring()))

    def _retract(self, chi):
        """
        Construct a Weyl character from an invariant element of the weight ring

        INPUT:

        - ``chi`` -- a linear combination of weights which
          shall be invariant under the action of the Weyl group

        OUTPUT: the corresponding Weyl character

        Please use instead the morphism :meth:`retract` which is
        implemented using this method.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)

        ::

            sage: v = A2._space([3,1,0]); v
            (3, 1, 0)
            sage: chi = a2.sum_of_monomials(v.orbit()); chi
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)
            sage: A2._retract(chi)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
        """
        return self.char_from_weights(dict(chi))

    @lazy_attribute
    def retract(self):
        """
        The partial inverse map from the weight ring into ``self``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: A2.retract
            Generic morphism:
              From: The Weight ring attached to The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients
              To:   The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients

        ::

            sage: v = A2._space([3,1,0]); v
            (3, 1, 0)
            sage: chi = a2.sum_of_monomials(v.orbit()); chi
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)
            sage: A2.retract(chi)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)

        The input should be invariant::

            sage: A2.retract(a2.monomial(v))
            Traceback (most recent call last):
            ...
            ValueError: multiplicity dictionary may not be Weyl group invariant

        As a shortcut, you may use conversion::

            sage: A2(chi)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
            sage: A2(a2.monomial(v))
            Traceback (most recent call last):
            ...
            ValueError: multiplicity dictionary may not be Weyl group invariant
        """
        from sage.categories.homset import Hom
        from sage.categories.morphism import SetMorphism
        category = Algebras(self.base_ring())
        return SetMorphism(Hom(self.ambient(), self, category), self._retract)

    def _repr_(self):
        """
        EXAMPLES::

            sage: WeylCharacterRing("A3")
            The Weyl Character Ring of Type ['A', 3] with Integer Ring coefficients
        """
        return "The Weyl Character Ring of Type {} with {} coefficients".format(self._cartan_type, self._base_ring)

    def __call__(self, *args):
        """
        Construct an element of ``self``.

        The input can either be an object that can be coerced or
        converted into ``self`` (an element of ``self``, of the base
        ring, of the weight ring), or a dominant weight. In the later
        case, the basis element indexed by that weight is returned.

        To specify the weight, you may give it explicitly. Alternatively,
        you may give a tuple of integers. Normally these are the
        components of the vector in the standard realization of
        the weight lattice as a vector space. Alternatively, if
        the ring is constructed with ``style = "coroots"``, you may
        specify the weight by giving a set of integers, one for each
        fundamental weight; the weight is then the linear combination
        of the fundamental weights with these coefficients.

        As a syntactical shorthand, for tuples of length at least two,
        the parenthesis may be omitted.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: [A2(x) for x in [-2,-1,0,1,2]]
            [-2*A2(0,0,0), -A2(0,0,0), 0, A2(0,0,0), 2*A2(0,0,0)]
            sage: [A2(2,1,0), A2([2,1,0]), A2(2,1,0)== A2([2,1,0])]
            [A2(2,1,0), A2(2,1,0), True]
            sage: A2([2,1,0]) == A2(2,1,0)
            True
            sage: l = -2*A2(0,0,0) - A2(1,0,0) + A2(2,0,0) + 2*A2(3,0,0)
            sage: [l in A2, A2(l) == l]
            [True, True]
            sage: P.<q> = QQ[]
            sage: A2 = WeylCharacterRing(['A',2], base_ring = P)
            sage: [A2(x) for x in [-2,-1,0,1,2,-2*q,-q,q,2*q,(1-q)]]
            [-2*A2(0,0,0), -A2(0,0,0), 0, A2(0,0,0), 2*A2(0,0,0), -2*q*A2(0,0,0), -q*A2(0,0,0),
            q*A2(0,0,0), 2*q*A2(0,0,0), (-q+1)*A2(0,0,0)]
            sage: R.<q> = ZZ[]
            sage: A2 = WeylCharacterRing(['A',2], base_ring = R, style="coroots")
            sage: q*A2(1)
            q*A2(0,0)
            sage: [A2(x) for x in [-2,-1,0,1,2,-2*q,-q,q,2*q,(1-q)]]
            [-2*A2(0,0), -A2(0,0), 0, A2(0,0), 2*A2(0,0), -2*q*A2(0,0), -q*A2(0,0), q*A2(0,0), 2*q*A2(0,0), (-q+1)*A2(0,0)]

        """
        # The purpose of this __call__ method is only to handle the
        # syntactical shorthand; otherwise it just delegates the work
        # to the coercion model, which itself will call
        # _element_constructor_ if the input is made of exactly one
        # object which can't be coerced into self
        if len(args) > 1:
            args = (args,)
        return super(WeylCharacterRing, self).__call__(*args)

    def _element_constructor_(self, weight):
        """
        Construct a monomial from a dominant weight.

        INPUT:

        - ``weight`` -- an element of the weight space, or a tuple

        This method is responsible for constructing an appropriate
        dominant weight from ``weight``, and then return the monomial
        indexed by that weight. See :meth:`__call__` and
        :meth:`sage.combinat.root_system.ambient_space.AmbientSpace.from_vector`.

        TESTS::

            sage: A2 = WeylCharacterRing("A2")
            sage: A2._element_constructor_([2,1,0])
            A2(2,1,0)
        """
        weight = self._space.from_vector_notation(weight, style = self._style)
        if not weight.is_dominant_weight():
            raise ValueError("{} is not a dominant element of the weight lattice".format(weight))
        return self.monomial(weight)

    def product_on_basis(self, a, b):
        r"""
        Compute the tensor product of two irreducible representations ``a``
        and ``b``.

        EXAMPLES::

            sage: D4 = WeylCharacterRing(['D',4])
            sage: spin_plus = D4(1/2,1/2,1/2,1/2)
            sage: spin_minus = D4(1/2,1/2,1/2,-1/2)
            sage: spin_plus * spin_minus # indirect doctest
            D4(1,0,0,0) + D4(1,1,1,0)
            sage: spin_minus * spin_plus
            D4(1,0,0,0) + D4(1,1,1,0)

        Uses the Brauer-Klimyk method.
        """
        # The method is asymmetrical, and as a rule of thumb
        # it is fastest to switch the factors so that the
        # smaller character is the one that is decomposed
        # into weights.
        if sum(a.coefficients()) > sum(b.coefficients()):
            a,b = b,a
        return self._product_helper(self._irr_weights(a), b)

    def _product_helper(self, d1, b):
        """
        Helper function for :meth:`product_on_basis`.

        INPUT:

        - ``d1`` -- a dictionary of weight multiplicities
        - ``b`` -- a dominant weight

        If ``d1`` is the dictionary of weight multiplicities of a character,
        returns the product of that character by the irreducible character
        with highest weight ``b``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: r = A2(1,0,0)
            sage: [A2._product_helper(r.weight_multiplicities(),x) for x in A2.space().fundamental_weights()]
            [A2(1,1,0) + A2(2,0,0), A2(1,1,1) + A2(2,1,0)]
        """
        d = {}
        for k in d1:
            [epsilon,g] = self.dot_reduce(b+k)
            if epsilon == 1:
                d[g] = d.get(g,0) + d1[k]
            elif epsilon == -1:
                d[g] = d.get(g,0)- d1[k]
        return self._from_dict(d)

    def dot_reduce(self, a):
        r"""
        Auxiliary function for :meth:`product_on_basis`.

        Return a pair `[\epsilon, b]` where `b` is a dominant weight and
        `\epsilon` is 0, 1 or -1. To describe `b`, let `w` be an element of
        the Weyl group such that `w(a + \rho)` is dominant. If
        `w(a + \rho) - \rho` is dominant, then `\epsilon` is the sign of
        `w` and `b` is `w(a + \rho) - \rho`. Otherwise, `\epsilon` is zero.

        INPUT:

        - ``a`` -- a weight

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2")
            sage: weights=A2(2,1,0).weight_multiplicities().keys(); weights
            [(1, 2, 0), (2, 1, 0), (0, 2, 1), (2, 0, 1), (0, 1, 2), (1, 1, 1), (1, 0, 2)]
            sage: [A2.dot_reduce(x) for x in weights]
            [[0, (0, 0, 0)], [1, (2, 1, 0)], [-1, (1, 1, 1)], [0, (0, 0, 0)], [0, (0, 0, 0)], [1, (1, 1, 1)], [-1, (1, 1, 1)]]
        """
        alphacheck = self._space.simple_coroots()
        alpha = self._space.simple_roots()
        sr = self._space.weyl_group().simple_reflections()
        [epsilon, ret] = [1,a]
        done = False
        while not done:
            done = True
            for i in self._space.index_set():
                c = ret.inner_product(alphacheck[i])
                if c == -1:
                    return [0, self._space.zero()]
                elif c < -1:
                    epsilon = -epsilon
                    ret -= (1+c)*alpha[i]
                    done = False
                    break
        return [epsilon, ret]

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A3").some_elements()
            [A3(1,0,0,0), A3(1,1,0,0), A3(1,1,1,0)]
        """
        return [self.monomial(x) for x in self.fundamental_weights()]

    def one_basis(self):
        """
        Return the index of 1 in ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A3").one_basis()
            (0, 0, 0, 0)
            sage: WeylCharacterRing("A3").one()
            A3(0,0,0,0)
        """
        return self._space.zero()

    @cached_method
    def _irr_weights(self, hwv):
        """
        Compute the weights of an irreducible as a dictionary.

        Given a dominant weight ``hwv``, this produces a dictionary of
        weight multiplicities for the irreducible representation
        with highest weight vector ``hwv``. This method is cached
        for efficiency.

        INPUT:

        - ``hwv`` -- a dominant weight

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2")
            sage: v = A2.fundamental_weights()[1]; v
            (1, 0, 0)
            sage: A2._irr_weights(v)
            {(0, 1, 0): 1, (1, 0, 0): 1, (0, 0, 1): 1}
        """
        if self._style == "coroots":
            return self._demazure_weights(hwv)
        else:
            return irreducible_character_freudenthal(hwv)

    def _demazure_weights(self, hwv, word="long", debug=False):
        """
        Computes the weights of a Demazure character.

        This method duplicates the functionality of :meth:`_irr_weights`, under
        the assumption that ``style = "coroots"``, but allows an optional
        parameter ``word``. (This is not allowed in :meth:`_irr_weights` since
        it would interfere with the ``@cached_method``.) Produces the
        dictionary of weights for the irreducible character with highest
        weight ``hwv`` when ``word`` is omitted, or for the Demazure character
        if ``word`` is included.

        INPUT:

        - ``hwv`` -- a dominant weight

        EXAMPLES::

            sage: B2=WeylCharacterRing("B2", style="coroots")
            sage: [B2._demazure_weights(v, word=[1,2]) for v in B2.fundamental_weights()]
            [{(0, 1): 1, (1, 0): 1}, {(-1/2, 1/2): 1, (1/2, -1/2): 1, (1/2, 1/2): 1}]
        """
        alphacheck = self._space.simple_coroots()
        alpha = self._space.simple_roots()
        dd = {}
        h = tuple(int(hwv.inner_product(alphacheck[j])) for j in self._space.index_set())
        dd[h] = int(1)
        return self._demazure_helper(dd, word=word, debug=debug)

    def _demazure_helper(self, dd, word="long", debug=False):
        r"""
        Assumes ``style = "coroots"``. If the optional parameter ``word`` is
        specified, produces a Demazure character (defaults to the long Weyl
        group element.

        INPUT:

        - ``dd`` -- a dictionary of weights

        - ``word`` -- (optional) a Weyl group reduced word

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2",style="coroots")
            sage: dd = {}; dd[(1,1)]=int(1)
            sage: A2._demazure_helper(dd,word=[1,2])
            {(1, -1, 0): 1, (-1, 1, 0): 1, (1, 0, -1): 1, (0, 0, 0): 1, (0, 1, -1): 1}
        """
        if self._style != "coroots":
            raise ValueError('_demazure_helper method unavailable. Use style="coroots".')
        index_set = self._space.index_set()
        alphacheck = self._space.simple_coroots()
        alpha = self._space.simple_roots()
        r = self.rank()
        cm = {}
        for i in index_set:
            cm[i] = tuple(int(alpha[i].inner_product(alphacheck[j])) for j in index_set)
            if debug:
                print "cm[%s]=%s"%(i,cm[i])
        accum = dd
        if word == "long":
            word = self._word
        for i in reversed(word):
            if debug:
                print "i=%s"%i
            next = {}
            for v in accum:
                coroot = v[i-1]
                if debug:
                    print "   v=%s, coroot=%s"%(v, coroot)
                if coroot >= 0:
                    mu = v
                    for j in range(coroot+1):
                        next[mu] = next.get(mu,0)+accum[v]
                        if debug:
                            print "     mu=%s, next[mu]=%s"%(mu, next[mu])
                        mu = tuple(mu[k] - cm[i][k] for k in range(r))
                else:
                    mu = v
                    for j in range(-1-coroot):
                        mu = tuple(mu[k] + cm[i][k] for k in range(r))
                        next[mu] = next.get(mu,0)-accum[v]
                        if debug:
                            print "     mu=%s, next[mu]=%s"%(mu, next[mu])
            accum = {}
            for v in next:
                accum[v] = next[v]
        ret = {}
        for v in accum:
            if accum[v]:
                ret[self._space.from_vector_notation(v, style="coroots")] = accum[v]
        return ret

    @cached_method
    def _weight_multiplicities(self, x):
        """
        Produce weight multiplicities for the (possibly reducible)
        WeylCharacter ``x``.

        EXAMPLES::

            sage: B2=WeylCharacterRing("B2",style="coroots")
            sage: chi=2*B2(1,0)
            sage: B2._weight_multiplicities(chi)
            {(0, 1): 2, (1, 0): 2, (0, 0): 2, (-1, 0): 2, (0, -1): 2}
        """
        d = {}
        m = x._monomial_coefficients
        for k in m:
            c = m[k]
            d1 = self._irr_weights(k)
            for l in d1:
                if l in d:
                    d[l] += c*d1[l]
                else:
                    d[l] = c*d1[l]
        for k in d.keys():
            if d[k] == 0:
                del d[k]
            else:
                d[k] = self._base_ring(d[k])
        return d

    def base_ring(self):
        """
        Return the base ring of ``self``.

        EXAMPLES::

            sage: R = WeylCharacterRing(['A',3], base_ring = CC); R.base_ring()
            Complex Field with 53 bits of precision
        """
        return self._base_ring

    def irr_repr(self, hwv):
        """
        Return a string representing the irreducible character with highest
        weight vector ``hwv``.

        EXAMPLES::

            sage: B3 = WeylCharacterRing("B3")
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(1,1,0)', 'B3(1/2,1/2,1/2)']
            sage: B3 = WeylCharacterRing("B3", style="coroots")
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(0,1,0)', 'B3(0,0,1)']
        """
        return self._prefix+self._wt_repr(hwv)

    def _wt_repr(self, wt):
        """
        Produce a representation of a vector in either coweight or
        lattice notation (following the appendices in Bourbaki, Lie Groups and
        Lie Algebras, Chapters 4,5,6), depending on whether the parent
        :class:`WeylCharacterRing` is created with ``style="coweights"``
        or not.

        EXAMPLES::

            sage: [fw1,fw2]=RootSystem("G2").ambient_space().fundamental_weights(); fw1,fw2
            ((1, 0, -1), (2, -1, -1))
            sage: [WeylCharacterRing("G2")._wt_repr(v) for v in [fw1,fw2]]
            ['(1,0,-1)', '(2,-1,-1)']
            sage: [WeylCharacterRing("G2",style="coroots")._wt_repr(v) for v in [fw1,fw2]]
            ['(1,0)', '(0,1)']
        """
        if self._style == "lattice":
            vec = wt.to_vector()
        elif self._style == "coroots":
            vec = [wt.inner_product(x) for x in self.simple_coroots()]
        else:
            raise ValueError("unknown style")
        hstring = str(vec[0])
        for i in range(1,len(vec)):
            hstring=hstring+","+str(vec[i])
        return "("+hstring+")"

    def _repr_term(self, t):
        """
        Representation of the monomial corresponding to a weight ``t``.

        EXAMPLES::

            sage: G2 = WeylCharacterRing("G2") # indirect doctest
            sage: [G2._repr_term(x) for x in G2.fundamental_weights()]
            ['G2(1,0,-1)', 'G2(2,-1,-1)']
        """
        return self.irr_repr(t)

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A2").cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def fundamental_weights(self):
        """
        Return the fundamental weights.

        EXAMPLES::

            sage: WeylCharacterRing("G2").fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Return the simple roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def simple_coroots(self):
        """
        Return the simple coroots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_coroots()
            Finite family {1: (0, 1, -1), 2: (1/3, -2/3, 1/3)}
        """
        return self._space.simple_coroots()

    def highest_root(self):
        """
        Return the highest_root.

        EXAMPLES::

            sage: WeylCharacterRing("G2").highest_root()
             (2, -1, -1)
        """
        return self._space.highest_root()

    def positive_roots(self):
        """
        Return the positive roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("E7").dynkin_diagram()
                    O 2
                    |
                    |
            O---O---O---O---O---O
            1   3   4   5   6   7
            E7
        """
        return self.space().dynkin_diagram()

    def extended_dynkin_diagram(self):
        """
        Return the extended Dynkin diagram, which is the Dynkin diagram
        of the corresponding untwisted affine type.

        EXAMPLES::

            sage: WeylCharacterRing("E7").extended_dynkin_diagram()
                        O 2
                        |
                        |
            O---O---O---O---O---O---O
            0   1   3   4   5   6   7
            E7~
        """
        return DynkinDiagram([self.cartan_type()[0],self.cartan_type()[1],1])

    def rank(self):
        """
        Return the rank.

        EXAMPLES::

            sage: WeylCharacterRing("G2").rank()
            2
        """
        return self._rank

    def space(self):
        """
        Return the weight space associated to ``self``.

        EXAMPLES::

            sage: WeylCharacterRing(['E',8]).space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def char_from_weights(self, mdict):
        """
        Construct a Weyl character from an invariant linear combination
        of weights.

        INPUT:

        - ``mdict`` -- a dictionary mapping weights to coefficients,
          and representing a linear combination of weights which
          shall be invariant under the action of the Weyl group

        OUTPUT: the corresponding Weyl character

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([3,1,0]); v
            (3, 1, 0)
            sage: d = dict([(x,1) for x in v.orbit()]); d
            {(3, 0, 1): 1, (1, 0, 3): 1, (0, 1, 3): 1, (1, 3, 0): 1, (3, 1, 0): 1, (0, 3, 1): 1}
            sage: A2.char_from_weights(d)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
        """
        return self._from_dict(self._char_from_weights(mdict), coerce=True)

    def _char_from_weights(self, mdict):
        """
        Helper method for :meth:'char_from_weights'.

        INPUT:

        - ``mdict`` -- a dictionary of weight multiplicities

        The output of this method is a dictionary whose keys are dominant
        weights that is the same as the :meth:`monomial_coefficients` method
        of ``self.char_from_weights()``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([3,1,0])
            sage: d = dict([(x,1) for x in v.orbit()])
            sage: A2._char_from_weights(d)
            {(3, 1, 0): 1, (2, 1, 1): -1, (2, 2, 0): -1}
        """
        hdict = {}
        ddict = mdict.copy()
        while len(ddict) != 0:
            highest = max((x.inner_product(self._space.rho()),x) for x in ddict)[1]
            if not highest.is_dominant():
                raise ValueError("multiplicity dictionary may not be Weyl group invariant")
            sdict = self._irr_weights(highest)
            c = ddict[highest]
            if highest in hdict:
                hdict[highest] += c
            else:
                hdict[highest] = c
            for k in sdict:
                if k in ddict:
                    if ddict[k] == c*sdict[k]:
                        del ddict[k]
                    else:
                        ddict[k] = ddict[k]-c*sdict[k]
                else:
                    ddict[k] = -c*sdict[k]
        return hdict

    class Element(CombinatorialFreeModule.Element):
        """
        A class for Weyl characters.
        """
        def cartan_type(self):
            """
            Return the Cartan type of ``self``.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2")
                sage: A2([1,0,0]).cartan_type()
                ['A', 2]
            """
            return self.parent()._cartan_type

        def degree(self):
            """
            The degree of ``self``, that is, the dimension of module.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: [B3(x).degree() for x in B3.fundamental_weights()]
                [7, 21, 8]
            """
            L = self.parent()._space
            return sum(L.weyl_dimension(k)*c for k,c in self)

        def branch(self, S, rule="default"):
            """
            Return the restriction of the character to the subalgebra. If no
            rule is specified, we will try to specify one.

            INPUT:

            - ``S`` -- a Weyl character ring for a Lie subgroup or subalgebra

            -  ``rule`` -- a branching rule

            See :func:`branch_weyl_character` for more information
            about branching rules.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: A2 = WeylCharacterRing(['A',2])
                sage: [B3(w).branch(A2,rule="levi") for w in B3.fundamental_weights()]
                [A2(0,0,0) + A2(1,0,0) + A2(0,0,-1),
                A2(0,0,0) + A2(1,0,0) + A2(1,1,0) + A2(1,0,-1) + A2(0,-1,-1) + A2(0,0,-1),
                A2(-1/2,-1/2,-1/2) + A2(1/2,-1/2,-1/2) + A2(1/2,1/2,-1/2) + A2(1/2,1/2,1/2)]
            """
            return branch_weyl_character(self, self.parent(), S, rule=rule)

        def __pow__(self, n):
            """
            We override the method in :module:`sage.monoids.monoids` since
            using the Brauer-Klimyk algorithm, it is more efficient to
            compute ``a*(a*(a*a))`` than ``(a*a)*(a*a)``.

            EXAMPLES::

                sage: B4 = WeylCharacterRing("B4",style="coroots")
                sage: spin = B4(0,0,0,1)
                sage: [spin^k for k in [0,1,3]]
                [B4(0,0,0,0), B4(0,0,0,1), 5*B4(0,0,0,1) + 4*B4(1,0,0,1) + 3*B4(0,1,0,1) + 2*B4(0,0,1,1) + B4(0,0,0,3)]
            """
            if n == 0:
                return self.parent().one()
            elif n == 1:
                return self
            else:
                return self*self.__pow__(n-1)

        def is_irreducible(self):
            """
            Return whether ``self`` is an irreducible character.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: [B3(x).is_irreducible() for x in B3.fundamental_weights()]
                [True, True, True]
                sage: sum(B3(x) for x in B3.fundamental_weights()).is_irreducible()
                False
            """
            return self.coefficients() == [1]

        @cached_method
        def symmetric_power(self, k):
            r"""
            Return the `k`-th symmetric power of ``self``.

            INPUT:

            - `k` -- a nonnegative integer

            The algorithm is based on the
            identity `k h_k = \sum_{r=1}^k p_k h_{k-r}` relating the power-sum
            and complete symmetric polynomials. Applying this to the
            eigenvalues of an element of the parent Lie group in the
            representation ``self``, the `h_k` become symmetric powers and
            the `p_k` become Adams operations, giving an efficient recursive
            implementation.

            EXAMPLES::

                sage: B3=WeylCharacterRing("B3",style="coroots")
                sage: spin=B3(0,0,1)
                sage: spin.symmetric_power(6)
                B3(0,0,0) + B3(0,0,2) + B3(0,0,4) + B3(0,0,6)
            """
            par = self.parent()
            if k == 0:
                return par.one()
            if k == 1:
                return self
            ret = par.zero()
            for r in range(1,k+1):
                adam_r = self._adams_operation_helper(r)
                ret += par.linear_combination( (par._product_helper(adam_r, l), c) for (l, c) in self.symmetric_power(k-r))
            dd = {}
            m = ret.weight_multiplicities()
            for l in m:
                dd[l] = m[l]/k
            return self.parent().char_from_weights(dd)

        @cached_method
        def exterior_power(self, k):
            r"""
            Return the `k`-th exterior power of ``self``.

            INPUT:

            - ``k`` -- a nonnegative integer

            The algorithm is based on the
            identity `k e_k = \sum_{r=1}^k (-1)^{k-1} p_k e_{k-r}` relating the
            power-sum and elementary symmetric polynomials. Applying this to
            the eigenvalues of an element of the parent Lie group in the
            representation ``self``, the `e_k` become exterior powers and
            the `p_k` become Adams operations, giving an efficient recursive
            implementation.

            EXAMPLES::

                sage: B3=WeylCharacterRing("B3",style="coroots")
                sage: spin=B3(0,0,1)
                sage: spin.exterior_power(6)
                B3(1,0,0) + B3(0,1,0)
            """
            par = self.parent()
            if k == 0:
                return par.one()
            if k == 1:
                return self
            ret = par.zero()
            for r in range(1,k+1):
                adam_r = self._adams_operation_helper(r)
                if is_even(r):
                    ret -= par.linear_combination( (par._product_helper(adam_r, l), c) for (l, c) in self.exterior_power(k-r))
                else:
                    ret += par.linear_combination( (par._product_helper(adam_r, l), c) for (l, c) in self.exterior_power(k-r))
            dd = {}
            m = ret.weight_multiplicities()
            for l in m:
                dd[l] = m[l]/k
            return self.parent().char_from_weights(dd)

        def adams_operation(self, r):
            """
            Return the `r`-th Adams operation of ``self``.

            INPUT:

            - ``r`` -- a positive integer

            This is a virtual character,
            whose weights are the weights of ``self``, each multiplied by `r`.

            EXAMPLES::

                sage: A2=WeylCharacterRing("A2")
                sage: A2(1,1,0).adams_operation(3)
                A2(2,2,2) - A2(3,2,1) + A2(3,3,0)
            """
            return self.parent().char_from_weights(self._adams_operation_helper(r))

        def _adams_operation_helper(self, r):
            """
            Helper function for Adams operations.

            INPUT:

            - ``r`` -- a positive integer

            Return the dictionary of weight multiplicities for the Adams
            operation, needed for internal use by symmetric and exterior powers.

            EXAMPLES::

                sage: A2=WeylCharacterRing("A2")
                sage: A2(1,1,0)._adams_operation_helper(3)
                {(3, 3, 0): 1, (0, 3, 3): 1, (3, 0, 3): 1}
            """
            d = self.weight_multiplicities()
            dd = {}
            for k in d:
                dd[r*k] = d[k]
            return dd

        def symmetric_square(self):
            """
            Return the symmetric square of the character.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2",style="coroots")
                sage: A2(1,0).symmetric_square()
                A2(2,0)
            """
            # Conceptually, this converts self to the weight ring,
            # computes its square there, and converts the result back.
            #
            # This implementation uses that this is a squaring (and not
            # a generic product) in the weight ring to optimize by
            # running only through pairs of weights instead of couples.
            c = self.weight_multiplicities()
            ckeys = c.keys()
            d = {}
            for j in range(len(ckeys)):
                for i in range(j+1):
                    ci = ckeys[i]
                    cj = ckeys[j]
                    t = ci + cj
                    if i < j:
                        coef = c[ci]*c[cj]
                    else:
                        coef = c[ci]*(c[ci]+1)/2
                    if t in d:
                        d[t] += coef
                    else:
                        d[t] = coef
            for k in d.keys():
                if d[k] == 0:
                    del d[k]
            return self.parent().char_from_weights(d)

        def exterior_square(self):
            """
            Return the exterior square of the character.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2",style="coroots")
                sage: A2(1,0).exterior_square()
                A2(0,1)
            """
            c = self.weight_multiplicities()
            ckeys = c.keys()
            d = {}
            for j in range(len(ckeys)):
                for i in range(j+1):
                    ci = ckeys[i]
                    cj = ckeys[j]
                    t = ci + cj
                    if i < j:
                        coef = c[ci]*c[cj]
                    else:
                        coef = c[ci]*(c[ci]-1)/2
                    if t in d:
                        d[t] += coef
                    else:
                        d[t] = coef
            for k in d.keys():
                if d[k] == 0:
                    del d[k]
            return self.parent().char_from_weights(d)

        def frobenius_schur_indicator(self):
            """
            Return:

            - `1` if the representation is real (orthogonal)

            - `-1` if the representation is quaternionic (symplectic)

            - `0` if the representation is complex (not self dual)

            The Frobenius-Schur indicator of a character `\chi`
            of a compact group `G` is the Haar integral over the
            group of `\chi(g^2)`. Its value is 1, -1 or 0. This
            method computes it for irreducible characters of
            compact Lie groups by checking whether the symmetric
            and exterior square characters contain the trivial
            character.

            .. TODO::

                Try to compute this directly without actually calculating
                the full symmetric and exterior squares.

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style="coroots")
                sage: B2(1,0).frobenius_schur_indicator()
                 1
                sage: B2(0,1).frobenius_schur_indicator()
                 -1
            """
            if not self.is_irreducible():
                raise ValueError("Frobenius-Schur indicator is only valid for irreducible characters")
            z = self.parent()._space.zero()
            if self.symmetric_square().coefficient(z) != 0:
                return 1
            if self.exterior_square().coefficient(z) != 0:
                return -1
            return 0

        def weight_multiplicities(self):
            """
            Produce the dictionary of weight multiplicities for the Weyl
            character ``self``. The character does not have to be irreducible.

            EXAMPLES::

                sage: B2=WeylCharacterRing("B2",style="coroots")
                sage: B2(0,1).weight_multiplicities()
                {(-1/2, 1/2): 1, (-1/2, -1/2): 1, (1/2, -1/2): 1, (1/2, 1/2): 1}
            """
            return self.parent()._weight_multiplicities(self)

        def inner_product(self, other):
            """
            Compute the inner product with another character.

            The irreducible characters are an orthonormal basis with respect
            to the usual inner product of characters, interpreted as functions
            on a compact Lie group, by Schur orthogonality.

            INPUT:

            - ``other`` -- another character

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2")
                sage: [f1,f2]=A2.fundamental_weights()
                sage: r1 = A2(f1)*A2(f2); r1
                A2(1,1,1) + A2(2,1,0)
                sage: r2 = A2(f1)^3; r2
                A2(1,1,1) + 2*A2(2,1,0) + A2(3,0,0)
                sage: r1.inner_product(r2)
                3
            """
            return sum(self.coefficient(x)*other.coefficient(x) for x in self.monomial_coefficients())

        def invariant_degree(self):
            """
            Return the multiplicity of the trivial representation in ``self``.

            Multiplicities of other irreducibles may be obtained
            using :meth:`multiplicity`.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2",style="coroots")
                sage: rep = A2(1,0)^2*A2(0,1)^2; rep
                2*A2(0,0) + A2(0,3) + 4*A2(1,1) + A2(3,0) + A2(2,2)
                sage: rep.invariant_degree()
                2
            """
            return self.coefficient(self.parent().space()(0))

        def multiplicity(self, other):
            """
            Return the multiplicity of the irreducible ``other`` in ``self``.

            INPUT:

            - ``other`` -- an irreducible character

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style="coroots")
                sage: rep = B2(1,1)^2; rep
                B2(0,0) + B2(1,0) + 2*B2(0,2) + B2(2,0) + 2*B2(1,2) + B2(0,4) + B2(3,0) + B2(2,2)
                sage: rep.multiplicity(B2(0,2))
                2
            """
            if not other.is_irreducible():
                raise ValueError("{} is not irreducible".format(other))
            return self.coefficient(other.support()[0])

def irreducible_character_freudenthal(hwv, debug=False):
    """
    Return the dictionary of multiplicities for the irreducible
    character with highest weight `\lambda`.

    The weight multiplicities are computed by the Freudenthal multiplicity
    formula. The algorithm is based on recursion relation that is stated,
    for example, in Humphrey's book on Lie Algebras. The multiplicities are
    invariant under the Weyl group, so to compute them it would be sufficient
    to compute them for the weights in the positive Weyl chamber. However
    after some testing it was found to be faster to compute every
    weight using the recursion, since the use of the Weyl group is
    expensive in its current implementation.

    INPUT:

    - ``hwv`` -- a dominant weight in a weight lattice.

    - ``L`` -- the ambient space

    EXAMPLES::

        sage: WeylCharacterRing("A2")(2,1,0).weight_multiplicities() # indirect doctest
        {(1, 2, 0): 1, (2, 1, 0): 1, (0, 2, 1): 1, (2, 0, 1): 1, (0, 1, 2): 1, (1, 1, 1): 2, (1, 0, 2): 1}
    """
    L = hwv.parent()
    rho = L.rho()
    mdict = {}
    current_layer = {hwv:1}

    simple_roots = L.simple_roots()
    positive_roots = L.positive_roots()

    while len(current_layer) > 0:
        next_layer = {}
        for mu in current_layer:
            if current_layer[mu] != 0:
                mdict[mu] = current_layer[mu]
                for alpha in simple_roots:
                    next_layer[mu-alpha] = None
        if debug:
            print next_layer

        for mu in next_layer:
            if next_layer[mu] is None:
                accum = 0
                for alpha in positive_roots:
                    mu_plus_i_alpha = mu + alpha
                    while mu_plus_i_alpha in mdict:
                        accum += mdict[mu_plus_i_alpha]*(mu_plus_i_alpha).inner_product(alpha)
                        mu_plus_i_alpha += alpha
                if accum == 0:
                    next_layer[mu] = 0
                else:
                    hwv_plus_rho = hwv + rho
                    mu_plus_rho  = mu  + rho
                    next_layer[mu] = ZZ(2*accum)/ZZ((hwv_plus_rho).inner_product(hwv_plus_rho)-(mu_plus_rho).inner_product(mu_plus_rho))
        current_layer = next_layer
    return mdict

def branch_weyl_character(chi, R, S, rule="default"):
    r"""
    A branching rule describes the restriction of representations from
    a Lie group or algebra `G` to a smaller one `H`. See for example, R. C.
    King, Branching rules for classical Lie groups using tensor and
    spinor methods. J. Phys. A 8 (1975), 429-449, Howe, Tan and
    Willenbring, Stable branching rules for classical symmetric pairs,
    Trans. Amer. Math. Soc. 357 (2005), no. 4, 1601-1626, McKay and
    Patera, Tables of Dimensions, Indices and Branching Rules for
    Representations of Simple Lie Algebras (Marcel Dekker, 1981),
    and Fauser, Jarvis, King and Wybourne, New branching rules induced
    by plethysm. J. Phys. A 39 (2006), no. 11, 2611--2655.

    INPUT:

    - ``chi`` -- a character of `G`

    - ``R`` -- the Weyl Character Ring of `G`

    - ``S`` -- the Weyl Character Ring of `H`

    - ``rule`` -- a set of `r` dominant weights in `H` where `r` is the rank
      of `G` or one of the following:

      * ``"levi"``
      * ``"automorphic"``
      * ``"symmetric"``
      * ``"extended"``
      * ``"orthogonal_sum"``
      * ``"tensor"``
      * ``"triality"``
      * ``"miscellaneous"``

    The use of the various input to ``rule`` will be explained next.
    After the examples we will explain how to write your own branching
    rules for cases that we have omitted.

    To explain the predefined rules, we survey the most important
    branching rules. These may be classified into several cases, and
    once this is understood, the detailed classification can be read
    off from the Dynkin diagrams. Dynkin classified the maximal
    subgroups of Lie groups in Mat. Sbornik N.S. 30(72):349-462 (1952).

    We will list give predefined rules that cover most cases where the
    branching rule is to a maximal subgroup. For convenience, we
    also give some branching rules to subgroups that are not maximal.
    For example, a Levi subgroup may or may not be maximal.

    You may try omitting the rule if it is "obvious". Default
    rules are provided for the following cases:

    .. MATH::

        \begin{aligned}
        A_{2s} & \to B_s,
        \\ A_{2s-1} & \to C_s,
        \\ A_{2*s-1} & \to D_s.
        \end{aligned}

    The above default rules correspond to embedding the group
    `SO(2s+1)`, `Sp(2s)` or `SO(2s)` into the corresponding general
    or special linear group by the standard representation. Default
    rules are also specified for the following cases:

    .. MATH::

        \begin{aligned}
        B_{s+1} & \to D_s,
        \\ D_s & \to B_s.
        \end{aligned}

    These correspond to the embedding of `O(n)` into `O(n+1)` where
    `n = 2s` or `2s + 1`. Finally, the branching rule for the embedding of
    a Levi subgroup is also implemented as a default rule.

    EXAMPLES::

        sage: A1 = WeylCharacterRing("A1", style="coroots")
        sage: A2 = WeylCharacterRing("A2", style="coroots")
        sage: D4 = WeylCharacterRing("D4", style="coroots")
        sage: B3 = WeylCharacterRing("B3", style="coroots")
        sage: B4 = WeylCharacterRing("B4", style="coroots")
        sage: A6 = WeylCharacterRing("A6", style="coroots")
        sage: A7 = WeylCharacterRing("A7", style="coroots")
        sage: def try_default_rule(R,S): return [R(f).branch(S) for f in R.fundamental_weights()]
        sage: try_default_rule(A2,A1)
        [A1(0) + A1(1), A1(0) + A1(1)]
        sage: try_default_rule(D4,B3)
        [B3(0,0,0) + B3(1,0,0), B3(1,0,0) + B3(0,1,0), B3(0,0,1), B3(0,0,1)]
        sage: try_default_rule(B4,D4)
        [D4(0,0,0,0) + D4(1,0,0,0), D4(1,0,0,0) + D4(0,1,0,0),
        D4(0,1,0,0) + D4(0,0,1,1), D4(0,0,1,0) + D4(0,0,0,1)]
        sage: try_default_rule(A7,D4)
        [D4(1,0,0,0), D4(0,1,0,0), D4(0,0,1,1), D4(0,0,2,0) + D4(0,0,0,2),
        D4(0,0,1,1),
        D4(0,1,0,0),
        D4(1,0,0,0)]
        sage: try_default_rule(A6,B3)
        [B3(1,0,0), B3(0,1,0), B3(0,0,2), B3(0,0,2), B3(0,1,0), B3(1,0,0)]

    If a default rule is not known, you may cue Sage as to what the
    Lie group embedding is by supplying a rule from the list of
    predefined rules. We will treat these next.

    .. RUBRIC:: Levi Type

    These can be read off from the Dynkin diagram. If
    removing a node from the Dynkin diagram produces another Dynkin
    diagram, there is a branching rule. Currently we require that the
    smaller diagram be connected. For these rules use the option
    ``rule="levi"``:

    .. MATH::

        \begin{aligned}
        A_r & \to A_{r-1}
        \\ B_r & \to A_{r-1}
        \\ B_r & \to B_{r-1}
        \\ C_r & \to A_{r-1}
        \\ C_r & \to C_{r-1}
        \\ D_r & \to A_{r-1}
        \\ D_r & \to D_{r-1}
        \\ E_r & \to A_{r-1} \quad r = 7,8
        \\ E_r & \to D_{r-1} \quad r = 6,7,8
        \\ E_r & \to E_{r-1}
        \\ F_4 & \to B_3
        \\ F_4 & \to C_3
        \\ G_2 & \to A_1 \text{(short root)}
        \end{aligned}

    Not all Levi subgroups are maximal subgroups. If the Levi is not
    maximal there may or may not be a preprogrammed ``rule="levi"`` for
    it. If there is not, the branching rule may still be obtained by going
    through an intermediate subgroup that is maximal using rule="extended".
    Thus the other Levi branching rule from `G_2 \to A_1` corresponding to the
    long root is available by first branching `G_2 \to A_2` then `A_2 \to A_1`.
    Similarly the branching rules to the Levi subgroup:

    .. MATH::

        E_r \to A_{r-1} \quad r = 6,7,8

    may be obtained by first branching `E_6 \to A_5 \times A_1`, `E_7 \to A_7`
    or `E_8 \to A_8`.

    EXAMPLES::

        sage: A1 = WeylCharacterRing("A1")
        sage: A2 = WeylCharacterRing("A2")
        sage: A3 = WeylCharacterRing("A3")
        sage: A4 = WeylCharacterRing("A4")
        sage: A5 = WeylCharacterRing("A5")
        sage: B2 = WeylCharacterRing("B2")
        sage: B3 = WeylCharacterRing("B3")
        sage: B4 = WeylCharacterRing("B4")
        sage: C2 = WeylCharacterRing("C2")
        sage: C3 = WeylCharacterRing("C3")
        sage: D3 = WeylCharacterRing("D3")
        sage: D4 = WeylCharacterRing("D4")
        sage: D5 = WeylCharacterRing("D5")
        sage: G2 = WeylCharacterRing("G2")
        sage: F4 = WeylCharacterRing("F4",style="coroots")
        sage: E6=WeylCharacterRing("E6",style="coroots")
        sage: D5=WeylCharacterRing("D5",style="coroots")
        sage: [B3(w).branch(A2,rule="levi") for w in B3.fundamental_weights()]
        [A2(0,0,0) + A2(1,0,0) + A2(0,0,-1),
        A2(0,0,0) + A2(1,0,0) + A2(1,1,0) + A2(1,0,-1) + A2(0,-1,-1) + A2(0,0,-1),
        A2(-1/2,-1/2,-1/2) + A2(1/2,-1/2,-1/2) + A2(1/2,1/2,-1/2) + A2(1/2,1/2,1/2)]

    The last example must be understood as follows. The representation
    of `B_3` being branched is spin, which is not a representation of
    `SO(7)` but of its double cover `\mathrm{spin}(7)`. The group `A_2` is
    really GL(3) and the double cover of `SO(7)` induces a cover of `GL(3)`
    that is trivial over `SL(3)` but not over the center of `GL(3)`. The weight
    lattice for this `GL(3)` consists of triples `(a,b,c)` of half integers
    such that `a - b` and `b - c` are in `\ZZ`, and this is reflected in the
    last decomposition.

    ::

        sage: [C3(w).branch(A2,rule="levi") for w in C3.fundamental_weights()]
        [A2(1,0,0) + A2(0,0,-1),
        A2(1,1,0) + A2(1,0,-1) + A2(0,-1,-1),
        A2(-1,-1,-1) + A2(1,-1,-1) + A2(1,1,-1) + A2(1,1,1)]
        sage: [D4(w).branch(A3,rule="levi") for w in D4.fundamental_weights()]
        [A3(1,0,0,0) + A3(0,0,0,-1),
        A3(0,0,0,0) + A3(1,1,0,0) + A3(1,0,0,-1) + A3(0,0,-1,-1),
        A3(1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,-1/2),
        A3(-1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,1/2)]
        sage: [B3(w).branch(B2,rule="levi") for w in B3.fundamental_weights()]
        [2*B2(0,0) + B2(1,0), B2(0,0) + 2*B2(1,0) + B2(1,1), 2*B2(1/2,1/2)]
        sage: C3 = WeylCharacterRing(['C',3])
        sage: [C3(w).branch(C2,rule="levi") for w in C3.fundamental_weights()]
        [2*C2(0,0) + C2(1,0),
         C2(0,0) + 2*C2(1,0) + C2(1,1),
         C2(1,0) + 2*C2(1,1)]
        sage: [D5(w).branch(D4,rule="levi") for w in D5.fundamental_weights()]
        [2*D4(0,0,0,0) + D4(1,0,0,0),
         D4(0,0,0,0) + 2*D4(1,0,0,0) + D4(1,1,0,0),
         D4(1,0,0,0) + 2*D4(1,1,0,0) + D4(1,1,1,0),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2)]
        sage: G2(1,0,-1).branch(A1,rule="levi")
        A1(1,0) + A1(1,-1) + A1(0,-1)
        sage: E6=WeylCharacterRing("E6",style="coroots")
        sage: D5=WeylCharacterRing("D5",style="coroots")
        sage: fw = E6.fundamental_weights()
        sage: [E6(fw[i]).branch(D5,rule="levi") for i in [1,2,6]] # long time (3s)
        [D5(0,0,0,0,0) + D5(0,0,0,0,1) + D5(1,0,0,0,0),
         D5(0,0,0,0,0) + D5(0,0,0,1,0) + D5(0,0,0,0,1) + D5(0,1,0,0,0),
         D5(0,0,0,0,0) + D5(0,0,0,1,0) + D5(1,0,0,0,0)]
        sage: E7=WeylCharacterRing("E7",style="coroots")
        sage: D6=WeylCharacterRing("D6",style="coroots")
        sage: fw = E7.fundamental_weights()
        sage: [E7(fw[i]).branch(D6,rule="levi") for i in [1,2,7]] # long time (26s)
        [3*D6(0,0,0,0,0,0) + 2*D6(0,0,0,0,1,0) + D6(0,1,0,0,0,0),
         3*D6(0,0,0,0,0,1) + 2*D6(1,0,0,0,0,0) + 2*D6(0,0,1,0,0,0) + D6(1,0,0,0,1,0),
         D6(0,0,0,0,0,1) + 2*D6(1,0,0,0,0,0)]
        sage: D7=WeylCharacterRing("D7",style="coroots")
        sage: E8=WeylCharacterRing("E8",style="coroots")
        sage: D7=WeylCharacterRing("D7",style="coroots")
        sage: E8(1,0,0,0,0,0,0,0).branch(D7,rule="levi") # not tested (very long time) (121s)
         3*D7(0,0,0,0,0,0,0) + 2*D7(0,0,0,0,0,1,0) + 2*D7(0,0,0,0,0,0,1) + 2*D7(1,0,0,0,0,0,0)
         + D7(0,1,0,0,0,0,0) + 2*D7(0,0,1,0,0,0,0) + D7(0,0,0,1,0,0,0) + D7(1,0,0,0,0,1,0) + D7(1,0,0,0,0,0,1) + D7(2,0,0,0,0,0,0)
        sage: E8(0,0,0,0,0,0,0,1).branch(D7,rule="levi") # long time (3s)
         D7(0,0,0,0,0,0,0) + D7(0,0,0,0,0,1,0) + D7(0,0,0,0,0,0,1) + 2*D7(1,0,0,0,0,0,0) + D7(0,1,0,0,0,0,0)
        sage: [F4(fw).branch(B3,rule="levi") for fw in F4.fundamental_weights()] # long time (36s)
         [B3(0,0,0) + 2*B3(1/2,1/2,1/2) + 2*B3(1,0,0) + B3(1,1,0),
         B3(0,0,0) + 6*B3(1/2,1/2,1/2) + 5*B3(1,0,0) + 7*B3(1,1,0) + 3*B3(1,1,1)
         + 6*B3(3/2,1/2,1/2) + 2*B3(3/2,3/2,1/2) + B3(2,0,0) + 2*B3(2,1,0) + B3(2,1,1),
         3*B3(0,0,0) + 6*B3(1/2,1/2,1/2) + 4*B3(1,0,0) + 3*B3(1,1,0) + B3(1,1,1) + 2*B3(3/2,1/2,1/2),
         3*B3(0,0,0) + 2*B3(1/2,1/2,1/2) + B3(1,0,0)]
        sage: [F4(fw).branch(C3,rule="levi") for fw in F4.fundamental_weights()] # long time (6s)
         [3*C3(0,0,0) + 2*C3(1,1,1) + C3(2,0,0),
         3*C3(0,0,0) + 6*C3(1,1,1) + 4*C3(2,0,0) + 2*C3(2,1,0) + 3*C3(2,2,0) + C3(2,2,2) + C3(3,1,0) + 2*C3(3,1,1),
         2*C3(1,0,0) + 3*C3(1,1,0) + C3(2,0,0) + 2*C3(2,1,0) + C3(2,1,1),
         2*C3(1,0,0) + C3(1,1,0)]
        sage: A1xA1 = WeylCharacterRing("A1xA1")
        sage: [A3(hwv).branch(A1xA1,rule="levi") for hwv in A3.fundamental_weights()]
        [A1xA1(1,0,0,0) + A1xA1(0,0,1,0),
        A1xA1(1,1,0,0) + A1xA1(1,0,1,0) + A1xA1(0,0,1,1),
        A1xA1(1,1,1,0) + A1xA1(1,0,1,1)]
        sage: A1xB1=WeylCharacterRing("A1xB1",style="coroots")
        sage: [B3(x).branch(A1xB1,rule="levi") for x in B3.fundamental_weights()]
        [2*A1xB1(1,0) + A1xB1(0,2),
        3*A1xB1(0,0) + 2*A1xB1(1,2) + A1xB1(2,0) + A1xB1(0,2),
        A1xB1(1,1) + 2*A1xB1(0,1)]

    .. RUBRIC:: Automorphic Type

    If the Dynkin diagram has a symmetry, then there
    is an automorphism that is a special case of a branching rule.
    There is also an exotic "triality" automorphism of `D_4` having order 3.
    Use ``rule="automorphic"`` (or for `D_4` ``rule="triality"``):

    .. MATH::

        \begin{aligned}
        A_r & \to A_r
        \\ D_r & \to D_r
        \\ E_6 & \to E_6
        \end{aligned}

    EXAMPLES::

        sage: [A3(chi).branch(A3,rule="automorphic") for chi in A3.fundamental_weights()]
        [A3(0,0,0,-1), A3(0,0,-1,-1), A3(0,-1,-1,-1)]
        sage: [D4(chi).branch(D4,rule="automorphic") for chi in D4.fundamental_weights()]
        [D4(1,0,0,0), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1/2,1/2,1/2,-1/2)]

    Here is an example with `D_4` triality::

        sage: [D4(chi).branch(D4,rule="triality") for chi in D4.fundamental_weights()]
        [D4(1/2,1/2,1/2,-1/2), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1,0,0,0)]

    .. RUBRIC:: Symmetric Type

    Related to the automorphic type, when `G` admits
    an outer automorphism (usually of degree 2) we may consider
    the branching rule to the isotropy subgroup `H`. In many cases
    the Dynkin diagram of `H` can be obtained by folding the Dynkin
    diagram of `G`. For such isotropy subgroups use ``rule="symmetric"``.
    The last branching rule, `D_4 \to G_2` is not to a maximal subgroup
    since `D_4 \to B_3 \to G_2`, but it is included for convenience.

    .. MATH::

        A_{2r+1} & \to B_r
        \\ A_{2r} & \to C_r
        \\ A_{2r} & \to D_r
        \\ D_r & \to B_{r-1}
        \\ E_6 & \to F_4
        \\ D_4 & \to G_2

    EXAMPLES::

        sage: [w.branch(B2,rule="symmetric") for w in [A4(1,0,0,0,0),A4(1,1,0,0,0),A4(1,1,1,0,0),A4(2,0,0,0,0)]]
        [B2(1,0), B2(1,1), B2(1,1), B2(0,0) + B2(2,0)]
        sage: [A5(w).branch(C3,rule="symmetric") for w in A5.fundamental_weights()]
        [C3(1,0,0), C3(0,0,0) + C3(1,1,0), C3(1,0,0) + C3(1,1,1), C3(0,0,0) + C3(1,1,0), C3(1,0,0)]
        sage: [A5(w).branch(D3,rule="symmetric") for w in A5.fundamental_weights()]
        [D3(1,0,0), D3(1,1,0), D3(1,1,-1) + D3(1,1,1), D3(1,1,0), D3(1,0,0)]
        sage: [D4(x).branch(B3,rule="symmetric") for x in D4.fundamental_weights()]
        [B3(0,0,0) + B3(1,0,0), B3(1,0,0) + B3(1,1,0), B3(1/2,1/2,1/2), B3(1/2,1/2,1/2)]
        sage: [D4(x).branch(G2,rule="symmetric") for x in D4.fundamental_weights()]
        [G2(0,0,0) + G2(1,0,-1), 2*G2(1,0,-1) + G2(2,-1,-1), G2(0,0,0) + G2(1,0,-1), G2(0,0,0) + G2(1,0,-1)]
        sage: [E6(fw).branch(F4,rule="symmetric") for fw in E6.fundamental_weights()] # long time (36s)
        [F4(0,0,0,0) + F4(0,0,0,1),
         F4(0,0,0,1) + F4(1,0,0,0),
         F4(0,0,0,1) + F4(1,0,0,0) + F4(0,0,1,0),
         F4(1,0,0,0) + 2*F4(0,0,1,0) + F4(1,0,0,1) + F4(0,1,0,0),
         F4(0,0,0,1) + F4(1,0,0,0) + F4(0,0,1,0),
         F4(0,0,0,0) + F4(0,0,0,1)]

    .. RUBRIC:: Extended Type

    If removing a node from the extended Dynkin diagram
    results in a Dynkin diagram, then there is a branching rule. Use
    ``rule="extended"`` for these. We will also use this classification
    for some rules that are not of this type, mainly involving type `B`,
    such as `D_6 \to `B_3 \times B_3`.

    Here is the extended Dynkin diagram for `D_6`::

            0       6
            O       O
            |       |
            |       |
        O---O---O---O---O
        1   2   3   4   6

    Removing the node 3 results in an embedding `D_3 \times D_3 \to D_6`.
    This corresponds to the embedding `SO(6) \times SO(6) \to SO(12)`, and
    is of extended type. On the other hand the embedding `SO(5) \times SO(7)
    \to SO(12)` (e.g. `B_2 \times B_3 \to D_6`) cannot be explained this way
    but for uniformity is implemented under ``rule="extended"``.

    The following rules are implemented as special cases
    of ``rule="extended"``:

    .. MATH::

        \begin{aligned}
        E_6 & \to A_5 \times A_1, A_2 \times A_2 \times A_2
        \\ E_7 & \to A_7, D_6 \times A_1, A_3 \times A_3 \times A_1
        \\ E_8 & \to A_8, D_8, E_7 \times A_1, A_4 \times A_4,
        D_5 \times A_3, E_6 \times A_2
        \\ F_4 & \to B_4, C_3 \times A_1, A_2 \times A_2, A_3 \times A_1
        \\ G_2 => A_1 \times A_1
        \end{aligned}

    Note that `E_8` has only a limited number of representations of
    reasonably low degree.

    EXAMPLES::

        sage: [B3(x).branch(D3,rule="extended") for x in B3.fundamental_weights()]
        [D3(0,0,0) + D3(1,0,0),
         D3(1,0,0) + D3(1,1,0),
         D3(1/2,1/2,-1/2) + D3(1/2,1/2,1/2)]
        sage: [G2(w).branch(A2, rule="extended") for w in G2.fundamental_weights()]
        [A2(0,0,0) + A2(1/3,1/3,-2/3) + A2(2/3,-1/3,-1/3),
         A2(1/3,1/3,-2/3) + A2(2/3,-1/3,-1/3) + A2(1,0,-1)]
        sage: [F4(fw).branch(B4,rule="extended") for fw in F4.fundamental_weights()] # long time (9s)
        [B4(1/2,1/2,1/2,1/2) + B4(1,1,0,0),
         B4(1,1,0,0) + B4(1,1,1,0) + B4(3/2,1/2,1/2,1/2) + B4(3/2,3/2,1/2,1/2) + B4(2,1,1,0),
         B4(1/2,1/2,1/2,1/2) + B4(1,0,0,0) + B4(1,1,0,0) + B4(1,1,1,0) + B4(3/2,1/2,1/2,1/2),
         B4(0,0,0,0) + B4(1/2,1/2,1/2,1/2) + B4(1,0,0,0)]

        sage: E6 = WeylCharacterRing("E6", style="coroots")
        sage: A2xA2xA2=WeylCharacterRing("A2xA2xA2",style="coroots")
        sage: A5xA1=WeylCharacterRing("A5xA1",style="coroots")
        sage: G2 = WeylCharacterRing("G2", style="coroots")
        sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
        sage: F4 = WeylCharacterRing("F4",style="coroots")
        sage: A3xA1 = WeylCharacterRing("A3xA1", style="coroots")
        sage: A2xA2 = WeylCharacterRing("A2xA2", style="coroots")
        sage: A1xC3 = WeylCharacterRing("A1xC3",style="coroots")

        sage: E6(1,0,0,0,0,0).branch(A5xA1,rule="extended") # (0.7s)
         A5xA1(0,0,0,1,0,0) + A5xA1(1,0,0,0,0,1)
        sage: E6(1,0,0,0,0,0).branch(A2xA2xA2, rule="extended") # (0.7s)
        A2xA2xA2(0,1,1,0,0,0) + A2xA2xA2(1,0,0,0,0,1) + A2xA2xA2(0,0,0,1,1,0)
        sage: E7=WeylCharacterRing("E7",style="coroots")
        sage: A7=WeylCharacterRing("A7",style="coroots")
        sage: E7(1,0,0,0,0,0,0).branch(A7,rule="extended") # long time (5s)
         A7(0,0,0,1,0,0,0) + A7(1,0,0,0,0,0,1)
        sage: E8=WeylCharacterRing("E8",style="coroots")
        sage: D8=WeylCharacterRing("D8",style="coroots")
        sage: E8(0,0,0,0,0,0,0,1).branch(D8,rule="extended") # long time (19s)
         D8(0,0,0,0,0,0,1,0) + D8(0,1,0,0,0,0,0,0)
        sage: F4(1,0,0,0).branch(A1xC3,rule="extended") # (0.7s)
         A1xC3(1,0,0,1) + A1xC3(2,0,0,0) + A1xC3(0,2,0,0)
        sage: G2(0,1).branch(A1xA1, rule="extended")
         A1xA1(2,0) + A1xA1(3,1) + A1xA1(0,2)
        sage: F4(0,0,0,1).branch(A2xA2, rule="extended") # (0.4s)
         A2xA2(0,1,0,1) + A2xA2(1,0,1,0) + A2xA2(0,0,1,1)
        sage: F4(0,0,0,1).branch(A3xA1,rule="extended") # (0.34s)
         A3xA1(0,0,0,0) + A3xA1(0,0,1,1) + A3xA1(0,1,0,0) + A3xA1(1,0,0,1) + A3xA1(0,0,0,2)
        sage: D4=WeylCharacterRing("D4",style="coroots")
        sage: D2xD2=WeylCharacterRing("D2xD2",style="coroots") # We get D4 => A1xA1xA1xA1 by remembering that A1xA1 = D2.
        sage: [D4(fw).branch(D2xD2, rule="extended") for fw in D4.fundamental_weights()]
         [D2xD2(1,1,0,0) + D2xD2(0,0,1,1),
         D2xD2(2,0,0,0) + D2xD2(0,2,0,0) + D2xD2(1,1,1,1) + D2xD2(0,0,2,0) + D2xD2(0,0,0,2),
         D2xD2(1,0,0,1) + D2xD2(0,1,1,0),
         D2xD2(1,0,1,0) + D2xD2(0,1,0,1)]

    .. RUBRIC:: Orthogonal Sum

    Using ``rule="orthogonal_sum"``, for `n = a + b + c + \cdots`,
    you can get any branching rule

    .. MATH::

        \begin{aligned}
        SO(n) & \to SO(a) \times SO(b) \times SO(c) \times \cdots,
        \\ Sp(2n) & \to Sp(2a) \times Sp(2b) \times Sp(2c) x \times \cdots,
        \end{aligned}

    where `O(a)` is type `D_r` for `a = 2r` or `B_r` for `a = 2r+1`
    and `Sp(2r)` is type `C_r`. In some cases these are also of
    extended type, as in the case `D_3 \times D_3 \to D_6` discussed above.
    But in other cases, for example `B_3 \times B_3 \to D_7`, they are not
    of extended type.

    .. RUBRIC:: Tensor

    There are branching rules:

    .. MATH::

        \begin{aligned}
        A_{rs-1} & \to A_{r-1} \times A_{s-1},
        \\ B_{2rs+r+s} & \to B_r \times B_s,
        \\ D_{2rs+s} & \to B_r \times D_s,
        \\ D_{2rs} & \to D_r \times D_s,
        \\ D_{2rs} & \to C_r \times C_s,
        \\ C_{2rs+s} & \to B_r \times C_s,
        \\ C_{2rs} & \to C_r \times D_s.
        \end{aligned}

    corresponding to the tensor product homomorphism. For type
    `A`, the homomorphism is `GL(r) \times GL(s) \to GL(rs)`. For the
    classical types, the relevant fact is that if `V, W` are
    orthogonal or symplectic spaces, that is, spaces endowed
    with symmetric or skew-symmetric bilinear forms, then `V \otimes W`
    is also an orthogonal space (if `V` and `W` are both
    orthogonal or both symplectic) or symplectic (if one of
    `V` and `W` is orthogonal and the other symplectic).

    The corresponding branching rules are obtained using ``rule="tensor"``.

    EXAMPLES::

        sage: A5=WeylCharacterRing("A5", style="coroots")
        sage: A2xA1=WeylCharacterRing("A2xA1", style="coroots")
        sage: [A5(hwv).branch(A2xA1, rule="tensor") for hwv in A5.fundamental_weights()]
        [A2xA1(1,0,1),
        A2xA1(0,1,2) + A2xA1(2,0,0),
        A2xA1(1,1,1) + A2xA1(0,0,3),
        A2xA1(1,0,2) + A2xA1(0,2,0),
        A2xA1(0,1,1)]
        sage: B4=WeylCharacterRing("B4",style="coroots")
        sage: B1xB1=WeylCharacterRing("B1xB1",style="coroots")
        sage: [B4(f).branch(B1xB1,rule="tensor") for f in B4.fundamental_weights()]
        [B1xB1(2,2),
        B1xB1(2,0) + B1xB1(2,4) + B1xB1(4,2) + B1xB1(0,2),
        B1xB1(2,0) + B1xB1(2,2) + B1xB1(2,4) + B1xB1(4,2) + B1xB1(4,4) + B1xB1(6,0) + B1xB1(0,2) + B1xB1(0,6),
        B1xB1(1,3) + B1xB1(3,1)]
        sage: D4=WeylCharacterRing("D4",style="coroots")
        sage: C2xC1=WeylCharacterRing("C2xC1",style="coroots")
        sage: [D4(f).branch(C2xC1,rule="tensor") for f in D4.fundamental_weights()]
        [C2xC1(1,0,1),
        C2xC1(0,1,2) + C2xC1(2,0,0) + C2xC1(0,0,2),
        C2xC1(1,0,1),
        C2xC1(0,1,0) + C2xC1(0,0,2)]
        sage: C3=WeylCharacterRing("C3",style="coroots")
        sage: B1xC1=WeylCharacterRing("B1xC1",style="coroots")
        sage: [C3(f).branch(B1xC1,rule="tensor") for f in C3.fundamental_weights()]
        [B1xC1(2,1), B1xC1(2,2) + B1xC1(4,0), B1xC1(4,1) + B1xC1(0,3)]

    .. RUBRIC:: Symmetric Power

    The `k`-th symmetric and exterior power homomorphisms map

    .. MATH::

        GL(n) \to GL\left(\binom{n+k-1}{k}\right)
        \times GL\left(\binom{n}{k}\right).

    The corresponding branching rules are not implemented but a special
    case is. The `k`-th symmetric power homomorphism `SL(2) \to GL(k+1)`
    has its image inside of `SO(2r+1)` if `k = 2r` and inside of `Sp(2r)` if
    `k = 2r - 1`. Hence there are branching rules:

    .. MATH::

        \begin{aligned}
        B_r & \to A_1
        \\ C_r & \to A_1
        \end{aligned}

    and these may be obtained using the rule "symmetric_power".

    EXAMPLES::

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: B3=WeylCharacterRing("B3",style="coroots")
        sage: C3=WeylCharacterRing("C3",style="coroots")
        sage: [B3(fw).branch(A1,rule="symmetric_power") for fw in B3.fundamental_weights()]
        [A1(6), A1(2) + A1(6) + A1(10), A1(0) + A1(6)]
        sage: [C3(fw).branch(A1,rule="symmetric_power") for fw in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    .. RUBRIC:: Miscellaneous

    Use ``rule="miscellaneous"`` for the following rules:

    .. MATH::

        \begin{aligned}
        B_3 & \to G_2,
        \\ F_4 & \to G_2 \times A_1 \text{(not implemented yet)}.
        \end{aligned}

    EXAMPLES::

        sage: G2 = WeylCharacterRing("G2")
        sage: [fw1, fw2, fw3] = B3.fundamental_weights()
        sage: B3(fw1+fw3).branch(G2, rule="miscellaneous")
        G2(1,0,-1) + G2(2,-1,-1) + G2(2,0,-2)

    .. RUBRIC:: Branching Rules From Plethysms

    Nearly all branching rules `G \to H` where `G` is of type `A`, `B`, `C`
    or `D` are covered by the preceding rules. The function
    :func:`branching_rule_from_plethysm` covers the remaining cases.

    EXAMPLES:

    This is a general rule that includes any branching rule
    from types `A`, `B`, `C`, or `D` as a special case. Thus it could be
    used in place of the above rules and would give the same
    results. However it is most useful when branching from `G`
    to a maximal subgroup `H` such that
    `\mathrm{rank}(H) < \mathrm{rank}(G) - 1`.

    We consider a homomorphism `H \to G` where `G` is one of
    `SL(r+1)`, `SO(2r+1)`, `Sp(2r)` or `SO(2r)`. The function
    :func:`branching_rule_from_plethysm` produces the corresponding
    branching rule. The main ingredient is the character
    `\chi` of the representation of `H` that is the homomorphism
    to `GL(r+1)`, `GL(2r+1)` or `GL(2r)`.

    This rule is so powerful that it contains the other
    rules implemented above as special cases. First let
    us consider the symmetric fifth power representation
    of `SL(2)`.

    ::

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: chi=A1([5])
        sage: chi.degree()
         6
        sage: chi.frobenius_schur_indicator()
        -1

    This confirms that the character has degree 6 and
    is symplectic, so it corresponds to a homomorphism
    `SL(2) \to Sp(6)`, and there is a corresponding
    branching rule `C_3 \to A_1`.

    ::

        sage: C3 = WeylCharacterRing("C3",style="coroots")
        sage: sym5rule = branching_rule_from_plethysm(chi,"C3")
        sage: [C3(hwv).branch(A1,rule=sym5rule) for hwv in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    This is identical to the results we would obtain using
    ``rule="symmetric_power"``. The next example gives a branching
    not available by other standard rules.

    ::

        sage: G2 = WeylCharacterRing("G2",style="coroots")
        sage: D7 = WeylCharacterRing("D7",style="coroots")
        sage: ad=G2(0,1); ad.degree(); ad.frobenius_schur_indicator()
         14
         1
        sage: spin = D7(0,0,0,0,0,1,0); spin.degree()
         64
        sage: spin.branch(G2, rule=branching_rule_from_plethysm(ad, "D7"))
         G2(1,1)

    We have confirmed that the adjoint representation of `G_2`
    gives a homomorphism into `SO(14)`, and that the pullback
    of the one of the two 64 dimensional spin representations
    to `SO(14)` is an irreducible representation of `G_2`.

    .. RUBRIC:: Isomorphic Type

    Although not usually referred to as a branching
    rule, the effects of the accidental isomorphisms may be handled
    using ``rule="isomorphic"``:

    .. MATH::

        \begin{aligned}
        B_2 & \to C_2
        \\ C_2 & \to B_2
        \\ A_3 & \to D_3
        \\ D_3 & \to A_3
        \\ D_2 & \to A_1 \to A_1
        \\ B_1 & \to A_1
        \\ C_1 & \to A_1
        \end{aligned}

    EXAMPLES::

        sage: [B2(x).branch(C2, rule="isomorphic") for x in B2.fundamental_weights()]
        [C2(1,1), C2(1,0)]
        sage: [C2(x).branch(B2, rule="isomorphic") for x in C2.fundamental_weights()]
        [B2(1/2,1/2), B2(1,0)]
        sage: [A3(x).branch(D3,rule="isomorphic") for x in A3.fundamental_weights()]
        [D3(1/2,1/2,1/2), D3(1,0,0), D3(1/2,1/2,-1/2)]
        sage: [D3(x).branch(A3,rule="isomorphic") for x in D3.fundamental_weights()]
        [A3(1/2,1/2,-1/2,-1/2), A3(1/4,1/4,1/4,-3/4), A3(3/4,-1/4,-1/4,-1/4)]

    Here `A_3(x,y,z,w)` can be understood as a representation of `SL(4)`.
    The weights `x,y,z,w` and `x+t,y+t,z+t,w+t` represent the same
    representation of `SL(4)` - though not of `GL(4)` - since
    `A_3(x+t,y+t,z+t,w+t)` is the same as `A_3(x,y,z,w)` tensored with
    `\mathrm{det}^t`. So as a representation of `SL(4)`,
    ``A3(1/4,1/4,1/4,-3/4)`` is the same as ``A3(1,1,1,0)``. The exterior
    square representation `SL(4) \to GL(6)` admits an invariant symmetric
    bilinear form, so is a representation `SL(4) \to SO(6)` that lifts to
    an isomorphism `SL(4) \to \mathrm{Spin}(6)`. Conversely, there are two
    isomorphisms `SO(6) \to SL(4)`, of which we've selected one.

    In cases like this you might prefer ``style="coroots"``::

        sage: A3 = WeylCharacterRing("A3",style="coroots")
        sage: D3 = WeylCharacterRing("D3",style="coroots")
        sage: [D3(fw) for fw in D3.fundamental_weights()]
        [D3(1,0,0), D3(0,1,0), D3(0,0,1)]
        sage: [D3(fw).branch(A3,rule="isomorphic") for fw in D3.fundamental_weights()]
        [A3(0,1,0), A3(0,0,1), A3(1,0,0)]
        sage: D2 = WeylCharacterRing("D2", style="coroots")
        sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
        sage: [D2(fw).branch(A1xA1,rule="isomorphic") for fw in D2.fundamental_weights()]
        [A1xA1(1,0), A1xA1(0,1)]

    .. RUBRIC:: Branching From a Reducible Root System

    If you are branching from a reducible root system, the rule is
    a list of rules, one for each component type in the root system.
    The rules in the list are given in pairs ``[type, rule]``, where
    type is the root system to be branched to, and rule is the
    branching rule.

    EXAMPLES::

        sage: D4 = WeylCharacterRing("D4",style="coroots")
        sage: D2xD2 = WeylCharacterRing("D2xD2",style="coroots")
        sage: A1xA1xA1xA1 = WeylCharacterRing("A1xA1xA1xA1",style="coroots")
        sage: rr = [["A1xA1","isomorphic"],["A1xA1","isomorphic"]]
        sage: [D4(fw) for fw in D4.fundamental_weights()]
        [D4(1,0,0,0), D4(0,1,0,0), D4(0,0,1,0), D4(0,0,0,1)]
        sage: [D4(fw).branch(D2xD2,rule="extended").branch(A1xA1xA1xA1,rule=rr) for fw in D4.fundamental_weights()]
        [A1xA1xA1xA1(1,1,0,0) + A1xA1xA1xA1(0,0,1,1),
        A1xA1xA1xA1(1,1,1,1) + A1xA1xA1xA1(2,0,0,0) + A1xA1xA1xA1(0,2,0,0) + A1xA1xA1xA1(0,0,2,0) + A1xA1xA1xA1(0,0,0,2),
        A1xA1xA1xA1(1,0,0,1) + A1xA1xA1xA1(0,1,1,0),
        A1xA1xA1xA1(1,0,1,0) + A1xA1xA1xA1(0,1,0,1)]

    .. RUBRIC:: Writing Your Own (Branching) Rules

    Suppose you want to branch from a group `G` to a subgroup `H`.
    Arrange the embedding so that a Cartan subalgebra `U` of `H` is
    contained in a Cartan subalgebra `T` of `G`. There is thus
    a mapping from the weight spaces `\mathrm{Lie}(T)^* \to \mathrm{Lie}(U)^*`.
    Two embeddings will produce identical branching rules if they
    differ by an element of the Weyl group of `H`.

    The *rule* is this map `\mathrm{Lie}(T)^*`, which is ``G.space()``, to
    `\mathrm{Lie}(U)^*`, which is ``H.space()``,
    which you may implement as a function. As an example, let
    us consider how to implement the branching rule `A_3 \to C_2`.
    Here `H = C_2 = Sp(4)` embedded as a subgroup in `A_3 = GL(4)`. The
    Cartan subalgebra `U` consists of diagonal matrices with
    eigenvalues `u_1, u_2, -u_2, -u_1`. The ``C2.space()`` is the
    two dimensional vector spaces consisting of the linear
    functionals `u_1` and `u_2` on `U`. On the other hand `\mathrm{Lie}(T)` is
    `\RR^4`. A convenient way to see the restriction is to
    think of it as the adjoint of the map `(u_1, u_2) \to
    (u_1,u_2, -u_2, -u_1)`,
    that is, `(x_0, x_1, x_2, x_3) \to (x_0 - x_3, x_1 - x_2)`. Hence we may
    encode the rule as follows::

       def rule(x):
           return [x[0]-x[3],x[1]-x[2]]

    or simply::

        rule = lambda x : [x[0]-x[3],x[1]-x[2]]

    EXAMPLES::

        sage: A3 = WeylCharacterRing(['A',3])
        sage: C2 = WeylCharacterRing(['C',2])
        sage: rule = lambda x : [x[0]-x[3],x[1]-x[2]]
        sage: branch_weyl_character(A3([1,1,0,0]),A3,C2,rule)
        C2(0,0) + C2(1,1)
        sage: A3(1,1,0,0).branch(C2, rule) == C2(0,0) + C2(1,1)
        True
    """
    if type(rule) == str:
        rule = get_branching_rule(R._cartan_type, S._cartan_type, rule)
    elif R._cartan_type.is_compound():
        Rtypes = R._cartan_type.component_types()
        Stypes = [CartanType(l[0]) for l in rule]
        rules = [l[1] for l in rule]
        ntypes = len(Rtypes)
        rule_list = [get_branching_rule(Rtypes[i], Stypes[i], rules[i]) for i in range(ntypes)]
        shifts = R._cartan_type._shifts
        def rule(x):
            yl = []
            for i in range(ntypes):
                yl.append(rule_list[i](x[shifts[i]:shifts[i+1]]))
            return flatten(yl)
    mdict = {}
    for k in chi.weight_multiplicities():
        # TODO: Could this use the new from_vector of ambient_space ?
        if S._style == "coroots":
            if S._cartan_type.is_atomic() and S._cartan_type[0] == 'E':
                if S._cartan_type[1] == 6:
                    h = S._space(rule(list(k.to_vector()))).coerce_to_e6()
                elif S._cartan_type[1] == 7:
                    h = S._space(rule(list(k.to_vector()))).coerce_to_e7()
            else:
                h = (S._space(rule(list(k.to_vector())))).coerce_to_sl()
        else:
            h = S._space(rule(list(k.to_vector())))
        chi_mdict = chi.weight_multiplicities()
        if h in mdict:
            mdict[h] += chi_mdict[k]
        else:
            mdict[h] = chi_mdict[k]
    return S.char_from_weights(mdict)

def get_branching_rule(Rtype, Stype, rule):
    """
    Creates a branching rule.

    INPUT:

    - ``R`` -- the Weyl Character Ring of `G`

    - ``S`` -- the Weyl Character Ring of `H`

    - ``rule`` -- a string describing the branching rule as a map from
      the weight space of `S` to the weight space of `R`.

    If the rule parameter is omitted, in a very few cases, a default
    rule is supplied. See
    :func:`~sage.combinat.root_system.weyl_characters.branch_weyl_character`.

    EXAMPLES::

       sage: rule = get_branching_rule(CartanType("A3"),CartanType("C2"),"symmetric")
       sage: [rule(x) for x in WeylCharacterRing("A3").fundamental_weights()]
       [[1, 0], [1, 1], [1, 0]]
    """
    r = Rtype.rank()
    s = Stype.rank()
    rdim = Rtype.root_system().ambient_space().dimension()
    sdim = Stype.root_system().ambient_space().dimension()
    if Stype.is_compound():
        stypes = Stype.component_types()
    if rule == "default":
        if not Rtype.is_compound():
            if Stype.is_compound() and s == r-1:
                try:
                    return get_branching_rule(Rtype, Stype, rule="levi")
                except StandardError:
                    pass

            if Rtype[0] == "A":
                if Stype[0] == "B" and r == 2*s:
                    return get_branching_rule(Rtype, Stype, rule="symmetric")
                elif Stype[0] == "C" and r == 2*s-1:
                    return get_branching_rule(Rtype, Stype, rule="symmetric")
                elif Stype[0] == "D" and r == 2*s-1:
                    return get_branching_rule(Rtype, Stype, rule="symmetric")
            elif Rtype[0] == "B" and Stype[0] == "D" and r == s:
                return get_branching_rule(Rtype, Stype, rule="extended")
            elif Rtype[0] == "D" and Stype[0] == "B" and r == s+1:
                return get_branching_rule(Rtype, Stype, rule="symmetric")

            if s == r-1:
                try:
                    return get_branching_rule(Rtype, Stype, rule="levi")
                except StandardError:
                    pass

        raise ValueError("No default rule found (you must specify the rule)")
    elif rule == "levi":
        if not s == r-1:
            raise ValueError("Incompatible ranks")
        if Rtype[0] == 'A':
            if Stype.is_compound():
                if all(ct[0]=='A' for ct in stypes) and rdim == sdim:
                    return lambda x : x
                else:
                    raise ValueError("Rule not found")
            elif Stype[0] == 'A':
                return lambda x : list(x)[:r]
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] in ['B', 'C', 'D']:
            if Stype.is_atomic():
                if Stype[0] == 'A':
                    return  lambda x : x
                elif Stype[0] == Rtype[0]:
                    return lambda x : list(x)[1:]
            elif stypes[-1][0] == Rtype[0] and all(t[0] == 'A' for t in stypes[:-1]):
                return lambda x : x
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] == 'E':
            if Stype.is_atomic():
                if Stype[0] == 'D':
                    if r == 6:
                        return lambda x : [-x[4],-x[3],-x[2],-x[1],-x[0]]
                    if r == 7:
                        return lambda x : [-x[5],-x[4],-x[3],-x[2],-x[1],-x[0]]
                    if r == 8:
                        return lambda x : [-x[6],-x[5],-x[4],-x[3],-x[2],-x[1],-x[0]]
                elif r in [7,8] and Stype[0] == 'E':
                    return lambda x : x
                elif Stype[0] == 'A':
                    if r == 6:
                        raise NotImplementedError('A5 Levi is not maximal. Branch to A5xA1 (rule="extended").')
                    if r == 7:
                        raise NotImplementedError('A5 Levi is not maximal. Branch to A5xA1 (rule="extended").')
            raise NotImplementedError("Not implemented yet")
        elif Rtype[0] == 'F' and s == 3:
            if Stype[0] == 'B':
                return lambda x : list(x)[1:]
            elif Stype[0] == 'C':
                return lambda x : [x[1]-x[0],x[2]+x[3],x[2]-x[3]]
            else:
                raise NotImplementedError("Not implemented yet")
        elif Rtype[0] == 'G' and Stype[0] == 'A':
            return lambda x : list(x)[1:]
        else:
            raise ValueError("Rule not found")
    elif rule == "automorphic":
        if not Rtype == Stype:
            raise ValueError("Cartan types must agree for automorphic branching rule")
        elif Rtype[0] == 'A':
            def rule(x) : y = [-i for i in x]; y.reverse(); return y
            return rule
        elif Rtype[0] == 'D':
            def rule(x) : x[len(x)-1] = -x[len(x)-1]; return x
            return rule
        elif Rtype[0] == 'E' and r == 6:
            M = matrix(QQ,[(3, 3, 3, -3, 0, 0, 0, 0), \
                           (3, 3, -3, 3, 0, 0, 0, 0), \
                           (3, -3, 3, 3, 0, 0, 0, 0), \
                           (-3, 3, 3, 3, 0, 0, 0, 0), \
                           (0, 0, 0, 0, -3, -3, -3, 3), \
                           (0, 0, 0, 0, -3, 5, -1, 1), \
                           (0, 0, 0, 0, -3, -1, 5, 1), \
                           (0, 0, 0, 0, 3, 1, 1, 5)])/6
            return lambda x : tuple(M*vector(x))
        else:
            raise ValueError("No automorphism found")
    elif rule == "triality":
        if not Rtype == Stype:
            raise ValueError("Triality is an automorphic type (for D4 only)")
        elif not Rtype[0] == 'D' and r == 4:
            raise ValueError("Triality is for D4 only")
        else:
            return lambda x : [(x[0]+x[1]+x[2]+x[3])/2,(x[0]+x[1]-x[2]-x[3])/2,(x[0]-x[1]+x[2]-x[3])/2,(-x[0]+x[1]+x[2]-x[3])/2]
    elif rule == "symmetric":
        if Rtype[0] == 'A':
            if (Stype[0] == 'C' or Stype[0] == 'D' and r == 2*s-1) or (Stype[0] == 'B' and r == 2*s):
                return lambda x : [x[i]-x[r-i] for i in range(s)]
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] == 'D' and Stype[0] == 'B' and s == r-1:
            return lambda x : x[:s]
        elif Rtype[0] == 'D' and r == 4 and Stype[0] == 'G':
            return lambda x : [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]]
        elif Rtype[0] == 'E' and Stype[0] == 'F' and r == 6 and s == 4:
            return lambda x : [(x[4]-3*x[5])/2,(x[0]+x[1]+x[2]+x[3])/2,(-x[0]-x[1]+x[2]+x[3])/2,(-x[0]+x[1]-x[2]+x[3])/2]
        else:
            raise ValueError("Rule not found")
    elif rule == "extended" or rule == "orthogonal_sum":
        if rule == "extended" and not s == r:
            raise ValueError('Ranks should be equal for rule="extended"')
        if Stype.is_compound():
            if Rtype[0] in ['B','D'] and all(t[0] in ['B','D'] for t in stypes):
                if Rtype[0] == 'D':
                    rdeg = 2*r
                else:
                    rdeg = 2*r+1
                sdeg = 0
                for t in stypes:
                    if t[0] == 'D':
                        sdeg += 2*t[1]
                    else:
                        sdeg += 2*t[1]+1
                if rdeg == sdeg:
                    return lambda x : x[:s]
                else:
                    raise ValueError("Rule not found")
            elif Rtype[0] == 'C':
                if all(t[0] == Rtype[0] for t in stypes):
                    return lambda x : x
            if rule == "orthogonal_sum":
                raise ValueError("Rule not found")
            elif Rtype[0] == 'E':
                if r == 6:
                    if stypes[0][0] == 'A' and stypes[0][1] == 5: # need doctest
                        if stypes[1][0] == 'A' and stypes[1][1] == 1:
                            M = matrix(QQ,[(-3, -3, -3, -3, -3, -5, -5, 5), \
                                           (-9, 3, 3, 3, 3, 1, 1, -1), \
                                           (3, -9, 3, 3, 3, 1, 1, -1), \
                                           (3, 3, -9, 3, 3, 1, 1, -1), \
                                           (3, 3, 3, -9, 3, 1, 1, -1), \
                                           (3, 3, 3, 3, -9, 9, -3, 3), \
                                           (-3, -3, -3, -3, -3, -1, 11, 1), \
                                           (3, 3, 3, 3, 3, 1, 1, 11)])/12
                            return lambda x : tuple(M*vector(x))
                    if len(stypes) == 3 and all(x[0] == 'A' and x[1] == 2 for x in stypes): # need doctest
                        M = matrix(QQ,[(0, 0, -2, -2, -2, -2, -2, 2), \
                                       (-3, 3, 1, 1, 1, 1, 1, -1), \
                                       (3, -3, 1, 1, 1, 1, 1, -1), \
                                       (0, 0, -2, -2, 4, 0, 0, 0), \
                                       (0, 0, -2, 4, -2, 0, 0, 0), \
                                       (0, 0, 4, -2, -2, 0, 0, 0), \
                                       (0, 0, -2, -2, -2, 2, 2, -2), \
                                       (3, 3, 1, 1, 1, -1, -1, 1), \
                                       (-3, -3, 1, 1, 1, -1, -1, 1)])/6
                        return lambda x : tuple(M*vector(x))
                elif r == 7:
                    if stypes[0][0] == 'D' and stypes[0][1] == 6 and stypes[1][0] == 'A' and stypes[1][1] == 1:
                        return lambda x : [x[5],x[4],x[3],x[2],x[1],x[0],x[6],x[7]] # need doctest
                    elif stypes[0][0] == 'A' and stypes[1][0] == 'A':
                        if stypes[0][1] == 5 and stypes[1][1] == 2:
                            M = matrix(QQ,[(5, 1, 1, 1, 1, 1, 0, 0), \
                                           (-1, -5, 1, 1, 1, 1, 0, 0), \
                                           (-1, 1, -5, 1, 1, 1, 0, 0), \
                                           (-1, 1, 1, -5, 1, 1, 0, 0), \
                                           (-1, 1, 1, 1, -5, 1, 0, 0), \
                                           (-1, 1, 1, 1, 1, -5, 0, 0), \
                                           (1, -1, -1, -1, -1, -1, 0, -6), \
                                           (1, -1, -1, -1, -1, -1, -6, 0), \
                                           (-2, 2, 2, 2, 2, 2, -3, -3)])/6
                            return lambda x : tuple(M*vector(x))
                        if len(stypes) == 3 and stypes[2][0] == 'A' and [stypes[0][1],stypes[1][1],stypes[2][1]] == [3,3,1]: # need doctest
                            M = matrix(QQ, [(0, 0, -1, -1, -1, -1, 2, -2), \
                                            (0, 0, -1, -1, -1, -1, -2, 2), \
                                            (-2, 2, 1, 1, 1, 1, 0, 0), \
                                            (2, -2, 1, 1, 1, 1, 0, 0), \
                                            (0, 0, -1, -1, -1, 3, 0, 0), \
                                            (0, 0, -1, -1, 3, -1, 0, 0), \
                                            (0, 0, -1, 3, -1, -1, 0, 0), \
                                            (0, 0, 3, -1, -1, -1, 0, 0), \
                                            (2, 2, 0, 0, 0, 0, -2, -2), \
                                            (-2, -2, 0, 0, 0, 0, -2, -2)])/4
                            return lambda x : tuple(M*vector(x))
                elif r == 8:
                    if stypes[0][0] == 'A' and stypes[1][0] == 'A':
                        if stypes[0][1] == 7 and stypes[1][1] == 1:
                            raise NotImplementedError("Not maximal: first branch to E7xA1")
                        elif stypes[0][1] == 4 and stypes[1][1] == 4:
                            M = matrix(QQ,[(0, 0, 0, -4, -4, -4, -4, 4), \
                                           (-5, 5, 5, 1, 1, 1, 1, -1), \
                                           (5, -5, 5, 1, 1, 1, 1, -1), \
                                           (5, 5, -5, 1, 1, 1, 1, -1), \
                                           (-5, -5, -5, 1, 1, 1, 1, -1), \
                                           (0, 0, 0, -8, 2, 2, 2, -2), \
                                           (0, 0, 0, 2, -8, 2, 2, -2), \
                                           (0, 0, 0, 2, 2, -8, 2, -2), \
                                           (0, 0, 0, 2, 2, 2, -8, -2), \
                                           (0, 0, 0, 2, 2, 2, 2, 8)])/10
                            return lambda x : tuple(M*vector(x))
                        elif len(stypes)==3:
                            if stypes[0][1] == 5 and stypes[1][0] == 2 and stypes[2][0] == 1:
                                raise NotImplementedError("Not maximal: first branch to A7xA1")
                    elif stypes[0][0] == 'D' and stypes[1][0] == 'A':
                        if stypes[0][1] == 5 and stypes[1][1] == 3:
                            raise NotImplementedError("Not maximal: first branch to D8 then D5xD3=D5xA3")
                    elif stypes[0][0] == 'E' and stypes[1][0] == 'A':
                        if stypes[0][1] == 6 and stypes[1][1] == 2:
                            return lambda x : [x[0],x[1],x[2],x[3],x[4], \
                                               (x[5]+x[6]-x[7])/3,(x[5]+x[6]-x[7])/3,(-x[5]-x[6]+x[7])/3, \
                                               (-x[5]-x[6]-2*x[7])/3,(-x[5]+2*x[6]+x[7])/3,(2*x[5]-x[6]+x[7])/3]
                        elif stypes[0][1] == 7 and stypes[1][1] == 1:
                            return lambda x : [x[0],x[1],x[2],x[3],x[4],x[5],(x[6]-x[7])/2,(-x[6]+x[7])/2,(-x[6]-x[7])/2,(x[6]+x[7])/2]
                raise ValueError("Rule not found")
            elif Rtype[0] == 'F':
                if stypes[0][0] == 'C' and stypes[0][1] == 3:
                    if stypes[1][0] == 'A' and stypes[1][1] == 1:
                        return lambda x : [x[0]-x[1],x[2]+x[3],x[2]-x[3],(-x[0]-x[1])/2,(x[0]+x[1])/2]
                if stypes[0][0] == 'A' and stypes[0][1] == 1:
                    if stypes[1][0] == 'C' and stypes[1][1] == 3:
                        return lambda x : [(-x[0]-x[1])/2,(x[0]+x[1])/2,x[0]-x[1],x[2]+x[3],x[2]-x[3]]
                if stypes[0][0] == 'A' and stypes[1][0] == 'A':
                    if stypes[0][1] == 2 and stypes[1][1] == 2:
                        M = matrix(QQ,[(-2, -1, -1, 0), (1, 2, -1, 0), (1, -1, 2, 0), (1, -1, -1, 3), (1, -1, -1, -3), (-2, 2, 2, 0)])/3
                    elif stypes[0][1] == 3 and stypes[1][1] == 1:
                        M = matrix(QQ,[(-3, -1, -1, -1), (1, 3, -1, -1), (1, -1, 3, -1), (1, -1, -1, 3), (2, -2, -2, -2), (-2, 2, 2, 2)])/4
                    elif stypes[0][1] == 1 and stypes[1][1] == 3:
                        M = matrix(QQ,[(2, -2, -2, -2), (-2, 2, 2, 2), (-3, -1, -1, -1), (1, 3, -1, -1), (1, -1, 3, -1), (1, -1, -1, 3)])/4
                    return lambda x : tuple(M*vector(x))
                else:
                    raise ValueError("Rule not found")
            elif Rtype[0] == 'G':
                if all(t[0] == 'A' and t[1] == 1 for t in stypes):
                    return lambda x : [(x[1]-x[2])/2,-(x[1]-x[2])/2, x[0]/2, -x[0]/2]
            else:
                raise ValueError("Rule not found")
        else: # irreducible Stype
            if Rtype[0] == 'B' and Stype[0] == 'D':
                return lambda x : x
            elif Rtype[0] == 'E':
                if r == 7:
                    if Stype[0] == 'A':
                        M = matrix(QQ, [(-1, -1, -1, -1, -1, -1, 2, -2), \
                                        (-1, -1, -1, -1, -1, -1, -2, 2), \
                                        (-3, 1, 1, 1, 1, 1, 0, 0), \
                                        (1, -3, 1, 1, 1, 1, 0, 0), \
                                        (1, 1, -3, 1, 1, 1, 0, 0), \
                                        (1, 1, 1, -3, 1, 1, 0, 0), \
                                        (1, 1, 1, 1, -3, 1, 2, 2), \
                                        (1, 1, 1, 1, 1, -3, 2, 2)])/4
                        return lambda x : tuple(M*vector(x))
                elif r == 8:
                    if Stype[0] == 'D':
                        return lambda x : [-x[7],x[6],x[5],x[4],x[3],x[2],x[1],x[0]]
                    elif Stype[0] == 'A':
                        M = matrix([(-2, -2, -2, -2, -2, -2, -2, 2), \
                                    (-5, 1, 1, 1, 1, 1, 1, -1), \
                                    (1, -5, 1, 1, 1, 1, 1, -1), \
                                    (1, 1, -5, 1, 1, 1, 1, -1), \
                                    (1, 1, 1, -5, 1, 1, 1, -1), \
                                    (1, 1, 1, 1, -5, 1, 1, -1), \
                                    (1, 1, 1, 1, 1, -5, 1, -1), \
                                    (1, 1, 1, 1, 1, 1, -5, -1), \
                                    (1, 1, 1, 1, 1, 1, 1, 5)])/6 # doctest needed
                        return lambda x : tuple(M*vector(x))
            elif Rtype[0] == 'F' and Stype[0] == 'B' and s == r:
                return lambda x : [-x[0], x[1], x[2], x[3]]
            elif Rtype[0] == 'G' and Stype[0] == 'A' and s == r:
                return lambda x : [(x[0]-x[2])/3, (-x[1]+x[2])/3, (-x[0]+x[1])/3]
            else:
                raise ValueError("Rule not found")
    elif rule == "isomorphic":
        if r != s:
            raise ValueError("Incompatible ranks")
        if Rtype == Stype:
            return lambda x : x
        elif Rtype[0] == 'B' and r == 2 and Stype[0] == 'C':
            def rule(x) : [x1, x2] = x; return [x1+x2, x1-x2]
            return rule
        elif Rtype[0] == 'B' and r == 1 and Stype[0] == 'A':
            return lambda x : [x[0],-x[0]]
        elif Rtype[0] == 'C' and r == 2 and Stype[0] == 'B':
            def rule(x) : [x1, x2] = x; return [(x1+x2)/2, (x1-x2)/2]
            return rule
        elif Rtype[0] == 'C' and r == 1 and Stype[0] == 'A':
            return lambda x : [x[0]/2,-x[0]/2]
        elif Rtype[0] == 'A' and r == 3 and Stype[0] == 'D':
            def rule(x): [x1, x2, x3, x4] = x; return [(x1+x2-x3-x4)/2, (x1-x2+x3-x4)/2, (x1-x2-x3+x4)/2]
            return rule
        elif Rtype[0] == 'D' and r == 2 and Stype.is_compound() and \
                 all(t[0] == 'A' for t in stypes):
            def rule(x): [t1, t2] = x; return [(t1-t2)/2, -(t1-t2)/2, (t1+t2)/2, -(t1+t2)/2]
            return rule
        elif Rtype[0] == 'D' and r == 3 and Stype[0] == 'A':
            def rule(x): [t1, t2, t3] = x; return [(t1+t2+t3)/2, (t1-t2-t3)/2, (-t1+t2-t3)/2, (-t1-t2+t3)/2]
            return rule
        else:
            raise ValueError("Rule not found")
    elif rule == "tensor" or rule == "tensor-debug":
        if not Stype.is_compound():
            raise ValueError("Tensor product requires more than one factor")
        if len(stypes) is not 2:
            raise ValueError("Not implemented")
        if Rtype[0] is 'A':
            nr = Rtype[1]+1
        elif Rtype[0] is 'B':
            nr = 2*Rtype[1]+1
        elif Rtype[0] in ['C', 'D']:
            nr = 2*Rtype[1]
        else:
            raise ValueError("Rule not found")
        [s1, s2] = [stypes[i][1] for i in range(2)]
        ns = [s1, s2]
        for i in range(2):
            if stypes[i][0] is 'A':
                ns[i] = ns[i]+1
            if stypes[i][0] is 'B':
                ns[i] = 2*ns[i]+1
            if stypes[i][0] in ['C','D']:
                ns[i] = 2*ns[i]
        if nr != ns[0]*ns[1]:
            raise ValueError("Ranks don't agree with tensor product")
        if Rtype[0] == 'A':
            if all(t[0] == 'A' for t in stypes):
                def rule(x):
                    ret = [sum(x[i*ns[1]:(i+1)*ns[1]]) for i in range(ns[0])]
                    ret.extend([sum(x[ns[1]*j+i] for j in range(ns[0])) for i in range(ns[1])])
                    return ret
                return rule
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] == 'B':
            if not all(t[0] == 'B' for t in stypes):
                raise ValueError("Rule not found")
        elif Rtype[0] == 'C':
            if stypes[0][0] in ['B','D'] and stypes[1][0] is 'C':
                pass
            elif stypes[1][0] in ['B','D'] and stypes[0][0] is 'C':
                pass
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] == 'D':
            if stypes[0][0] in ['B','D'] and stypes[1][0] is 'D':
                pass
            elif stypes[1][0] is 'B' and stypes[0][0] is 'D':
                pass
            elif stypes[1][0] is 'C' and stypes[0][0] is 'C':
                pass
            else:
                raise ValueError("Rule not found")
        rows = []
        for i in range(s1):
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                nextrow[s1+j] = 1
                rows.append(nextrow)
        if stypes[1][0] == 'B':
            for i in range(s1):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                rows.append(nextrow)
        for i in range(s1):
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                nextrow[s1+j] = -1
                rows.append(nextrow)
        if stypes[0][0] == 'B':
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[s1+j] = 1
                rows.append(nextrow)
        mat = matrix(rows).transpose()
        if rule == "tensor-debug":
            print mat
        return lambda x : tuple(mat*vector(x))
    elif rule == "symmetric_power":
        if Stype[0] == 'A' and s == 1:
            if Rtype[0] == 'B':
                def rule(x):
                    a = sum((r-i)*x[i] for i in range(r))
                    return [a,-a]
                return rule
            elif Rtype[0] == 'C':
                def rule(x):
                    a = sum((2*r-2*i-1)*x[i] for i in range(r))
                    return [a/2,-a/2]
                return rule
            else:
                raise ValueError("Rule not found")
        else:
            raise ValueError("Rule not found")
    elif rule == "miscellaneous":
        if Rtype[0] == 'B' and Stype[0] == 'G' and r == 3:
            return lambda x : [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]]
        else:
            raise ValueError("Rule not found")

def branching_rule_from_plethysm(chi, cartan_type, return_matrix = False):
    r"""
    Create the branching rule of a plethysm.

    INPUT:

    - ``chi`` -- the character of an irreducible representation `\pi` of
      a group `G`
    - ``cartan_type`` -- a classical Cartan type (`A`,`B`,`C` or `D`).

    It is assumed that the image of the irreducible representation pi
    naturally has its image in the group `G`.

    Returns a branching rule for this plethysm.

    EXAMPLES:

    The adjoint representation `SL(3) \to GL(8)` factors
    through `SO(8)`. The branching rule in question will
    describe how representations of `SO(8)` composed with
    this homomorphism decompose into irreducible characters
    of `SL(3)`::

        sage: A2 = WeylCharacterRing("A2")
        sage: A2 = WeylCharacterRing("A2", style="coroots")
        sage: ad = A2(1,1)
        sage: ad.degree()
        8
        sage: ad.frobenius_schur_indicator()
        1

    This confirms that `ad` has degree 8 and is orthogonal,
    hence factors through `SO(8)` which is type `D_4`::

        sage: br = branching_rule_from_plethysm(ad,"D4")
        sage: D4 = WeylCharacterRing("D4")
        sage: [D4(f).branch(A2,rule = br) for f in D4.fundamental_weights()]
        [A2(1,1), A2(0,3) + A2(1,1) + A2(3,0), A2(1,1), A2(1,1)]
    """
    ct = CartanType(cartan_type)
    if ct[0] not in ["A","B","C","D"]:
        raise ValueError("not implemented for type {}".format(ct[0]))
    if ct[0] is "A":
        ret = []
        ml = chi.weight_multiplicities()
        for v in ml:
            n = ml[v]
            ret.extend(n*[v.to_vector()])
        M = matrix(ret).transpose()
        if len(M.columns()) != ct[1] + 1:
            raise ValueError("representation has wrong degree for type {}".format(ct))
        return lambda x : tuple(M*vector(x))
    if ct[0] in ["B","D"]:
        if chi.frobenius_schur_indicator() != 1:
            raise ValueError("character is not orthogonal")
    if ct[0] is "C":
        if chi.frobenius_schur_indicator() != -1:
            raise ValueError("character is not symplectic")
    if ct[0] is "B":
        if is_even(chi.degree()):
            raise ValueError("degree is not odd")
    if ct[0] is ["C","D"]:
        if is_odd(chi.degree()):
            raise ValueError("degree is not even")
    ret = []
    ml = chi.weight_multiplicities()
    for v in ml:
        n = ml[v]
        vec = v.to_vector()
        if all(x==0 for x in vec):
            if ct[0] is "B":
                n = (n-1)/2
            else:
                n = n/2
        elif [x for x in vec if x !=0][0] < 0:
            continue
        ret.extend(n*[vec])
    M = matrix(ret).transpose()
    if len(M.columns()) != ct.root_system().ambient_space().dimension():
        raise ValueError("representation has wrong degree for type {}".format(ct))
    if return_matrix:
        return M
    else:
        return lambda x : tuple(M*vector(x))

class WeightRing(CombinatorialFreeModule):
    """
    The weight ring, which is the group algebra over a weight lattice.

    A Weyl character may be regarded as an element of the weight ring.
    In fact, an element of the weight ring is an element of the
    :class:`weyl character ring <WeylCharacterRing>` if and only if it is
    invariant under the action of the Weyl group.

    The advantage of the weight ring over the Weyl character ring
    is that one may conduct calculations in the weight ring that
    involve sums of weights that are not Weyl group invariant.

    EXAMPLES::

        sage: A2 = WeylCharacterRing(['A',2])
        sage: a2 = WeightRing(A2)
        sage: wd = prod(a2(x/2)-a2(-x/2) for x in a2.space().positive_roots()); wd
        a2(-1,1,0) - a2(-1,0,1) - a2(1,-1,0) + a2(1,0,-1) + a2(0,-1,1) - a2(0,1,-1)
        sage: chi = A2([5,3,0]); chi
        A2(5,3,0)
        sage: a2(chi)
        a2(1,2,5) + 2*a2(1,3,4) + 2*a2(1,4,3) + a2(1,5,2) + a2(2,1,5)
        + 2*a2(2,2,4) + 3*a2(2,3,3) + 2*a2(2,4,2) + a2(2,5,1) + 2*a2(3,1,4)
        + 3*a2(3,2,3) + 3*a2(3,3,2) + 2*a2(3,4,1) + a2(3,5,0) + a2(3,0,5)
        + 2*a2(4,1,3) + 2*a2(4,2,2) + 2*a2(4,3,1) + a2(4,4,0) + a2(4,0,4)
        + a2(5,1,2) + a2(5,2,1) + a2(5,3,0) + a2(5,0,3) + a2(0,3,5)
        + a2(0,4,4) + a2(0,5,3)
        sage: a2(chi)*wd
        -a2(-1,3,6) + a2(-1,6,3) + a2(3,-1,6) - a2(3,6,-1) - a2(6,-1,3) + a2(6,3,-1)
        sage: sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.space().weyl_group())
        -a2(-1,3,6) + a2(-1,6,3) + a2(3,-1,6) - a2(3,6,-1) - a2(6,-1,3) + a2(6,3,-1)
        sage: a2(chi)*wd == sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.space().weyl_group())
        True
    """
    @staticmethod
    def __classcall__(cls, parent, prefix=None):
        """
        TESTS::

            sage: A3 = WeylCharacterRing("A3", style="coroots")
            sage: a3 = WeightRing(A3)
            sage: a3.cartan_type(), a3.base_ring(), a3.parent()
            (['A', 3], Integer Ring, The Weyl Character Ring of Type ['A', 3] with Integer Ring coefficients)
        """
        return super(WeightRing, cls).__classcall__(cls, parent, prefix=prefix)

    def __init__(self, parent, prefix):
        """
        EXAMPLES:

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: TestSuite(a2).run()

        TESTS::

            sage: A1xA1 = WeylCharacterRing("A1xA1")
            sage: a1xa1 = WeightRing(A1xA1)
            sage: TestSuite(a1xa1).run()
            sage: a1xa1.an_element()
            a1xa1(2,2,3,0)
        """
        self._parent = parent
        self._style = parent._style
        self._prefix = prefix
        self._space = parent._space
        self._cartan_type = parent._cartan_type
        self._rank = parent._rank
        self._origin = parent._origin
        self._base_ring = parent._base_ring
        if prefix is None:
            # TODO: refactor this fragile logic into CartanType's
            if self._parent._prefix.replace('x','_').isupper():
                # The 'x' workaround above is to support reducible Cartan types like 'A1xB2'
                prefix = self._parent._prefix.lower()
            elif self._parent._prefix.islower():
                prefix = self._parent._prefix.upper()
            else:
                # TODO: this only works for irreducible Cartan types!
                prefix = (self._cartan_type[0].lower()+str(self._rank))
        self._prefix = prefix
        category = AlgebrasWithBasis(self._base_ring)
        CombinatorialFreeModule.__init__(self, self._base_ring, self._space, category = category)


    def _repr_(self):
        """
        EXAMPLES::

            sage: P.<q>=QQ[]
            sage: G2 = WeylCharacterRing(['G',2], base_ring = P)
            sage: WeightRing(G2) # indirect doctest
            The Weight ring attached to The Weyl Character Ring of Type ['G', 2] with Univariate Polynomial Ring in q over Rational Field coefficients
        """
        return "The Weight ring attached to %s"%self._parent

    def __call__(self, *args):
        """
        Construct an element of ``self``.

        The input can either be an object that can be coerced or
        converted into ``self`` (an element of ``self``, of the base
        ring, of the weight ring), or a dominant weight. In the later
        case, the basis element indexed by that weight is returned.

        To specify the weight, you may give it explicitly. Alternatively,
        you may give a tuple of integers. Normally these are the
        components of the vector in the standard realization of
        the weight lattice as a vector space. Alternatively, if
        the ring is constructed with style="coroots", you may
        specify the weight by giving a set of integers, one for each
        fundamental weight; the weight is then the linear combination
        of the fundamental weights with these coefficients.

        As a syntactical shorthand, for tuples of length at least two,
        the parenthesis may be omitted.

        EXAMPLES::

            sage: a2 = WeightRing(WeylCharacterRing(['A',2]))
            sage: a2(-1)
            -a2(0,0,0)
        """
        # The purpose of this __call__ method is only to handle the
        # syntactical shorthand; otherwise it just delegates the work
        # to the coercion model, which itself will call
        # _element_constructor_ if the input is made of exactly one
        # object which can't be coerced into self
        if len(args) > 1:
            args = (args,)
        return super(WeightRing, self).__call__(*args)

    def _element_constructor_(self, weight):
        """
        Construct a monomial from a weight.

        INPUT:

        -  ``weight`` -- an element of the weight space, or a tuple

        This method is responsible for constructing an appropriate
        weight from the data in ``weight``, and then return the
        monomial indexed by that weight. See :meth:`__call__` and
        :meth:`sage.combinat.root_system.ambient_space.AmbientSpace.from_vector`.

        TESTS::

            sage: A2 = WeylCharacterRing("A2")
            sage: A2._element_constructor_([2,1,0])
            A2(2,1,0)
        """
        weight = self._space.from_vector_notation(weight, style = self._style)
        return self.monomial(weight)

    def product_on_basis(self, a, b):
        """
        Return the product of basis elements indexed by ``a`` and ``b``.

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2")
            sage: a2=WeightRing(A2)
            sage: a2(1,0,0) * a2(0,1,0) # indirect doctest
            a2(1,1,0)
        """
        return self(a+b)

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: A3=WeylCharacterRing("A3")
            sage: a3=WeightRing(A3)
            sage: a3.some_elements()
            [a3(1,0,0,0), a3(1,1,0,0), a3(1,1,1,0)]
        """
        return [self.monomial(x) for x in self.fundamental_weights()]

    def one_basis(self):
        """
        Return the index of `1`.

        EXAMPLES::

            sage: A3=WeylCharacterRing("A3")
            sage: WeightRing(A3).one_basis()
            (0, 0, 0, 0)
            sage: WeightRing(A3).one()
            a3(0,0,0,0)
        """
        return self._space.zero()

    def parent(self):
        """
        Return the parent Weyl character ring.

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2")
            sage: a2=WeightRing(A2)
            sage: a2.parent()
            The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients
            sage: a2.parent() == A2
            True

        """
        return self._parent

    def weyl_character_ring(self):
        """
        Return the parent Weyl Character Ring. A synonym for ``self.parent()``.

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2")
            sage: a2=WeightRing(A2)
            sage: a2.weyl_character_ring()
            The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients
        """
        return self._parent

    def cartan_type(self):
        """
        Return the Cartan type.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: WeightRing(A2).cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def space(self):
        """
        Return the weight space realization associated to ``self``.

        EXAMPLES::

            sage: E8 = WeylCharacterRing(['E',8])
            sage: e8 = WeightRing(E8)
            sage: e8.space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def fundamental_weights(self):
        """
        Return the fundamental weights.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Return the simple roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def positive_roots(self):
        """
        Return the positive roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def wt_repr(self, wt):
        r"""
        Return a string representing the irreducible character with
        highest weight vector ``wt``. Uses coroot notation if the associated
        Weyl character ring is defined with ``style="coroots"``.

        EXAMPLES::

            sage: G2 = WeylCharacterRing("G2")
            sage: [G2.ambient().wt_repr(x) for x in G2.fundamental_weights()]
            ['g2(1,0,-1)', 'g2(2,-1,-1)']
            sage: G2 = WeylCharacterRing("G2",style="coroots")
            sage: [G2.ambient().wt_repr(x) for x in G2.fundamental_weights()]
            ['g2(1,0)', 'g2(0,1)']
        """
        return self._prefix+self.parent()._wt_repr(wt)

    def _repr_term(self, t):
        """
        Representation of the monomial corresponding to a weight ``t``.

        EXAMPLES::

            sage: G2=WeylCharacterRing("G2")
            sage: g2=WeightRing(G2)
            sage: [g2(x) for x in g2.fundamental_weights()] # indirect doctest
            [g2(1,0,-1), g2(2,-1,-1)]
        """
        return self.wt_repr(t)

    class Element(CombinatorialFreeModule.Element):
        """
        A class for weight ring elements.
        """
        def cartan_type(self):
            """
            Return the Cartan type.

            EXAMPLES::

                sage: A2=WeylCharacterRing("A2")
                sage: a2 = WeightRing(A2)
                sage: a2([0,1,0]).cartan_type()
                ['A', 2]
            """
            return self.parent()._cartan_type

        def weyl_group_action(self, w):
            """
            Return the action of the Weyl group element ``w`` on ``self``.

            EXAMPLES::

                sage: G2 = WeylCharacterRing(['G',2])
                sage: g2 = WeightRing(G2)
                sage: L = g2.space()
                sage: [fw1, fw2] = L.fundamental_weights()
                sage: sum(g2(fw2).weyl_group_action(w) for w in L.weyl_group())
                2*g2(-2,1,1) + 2*g2(-1,-1,2) + 2*g2(-1,2,-1) + 2*g2(1,-2,1) + 2*g2(1,1,-2) + 2*g2(2,-1,-1)
            """
            return self.map_support(w.action)

        def character(self):
            """
            Assuming that ``self`` is invariant under the Weyl group, this will
            express it as a linear combination of characters. If ``self`` is
            not Weyl group invariant, this method will not terminate.

            EXAMPLES::

                sage: A2 = WeylCharacterRing(['A',2])
                sage: a2 = WeightRing(A2)
                sage: W = a2.space().weyl_group()
                sage: mu = a2(2,1,0)
                sage: nu = sum(mu.weyl_group_action(w) for w in W) ; nu
                a2(1,2,0) + a2(1,0,2) + a2(2,1,0) + a2(2,0,1) + a2(0,1,2) + a2(0,2,1)
                sage: nu.character()
                -2*A2(1,1,1) + A2(2,1,0)
            """
            return self.parent().parent().char_from_weights(self.monomial_coefficients())

        def scale(self, k):
            """
            Multiplies a weight by `k`. The operation is extended by linearity
            to the weight ring.

            INPUT:

            - ``k`` -- a nonzero integer

            EXAMPLES::

                sage: g2 = WeylCharacterRing("G2",style="coroots").ambient()
                sage: g2(2,3).scale(2)
                g2(4,6)
            """
            if k == 0:
                raise ValueError("parameter must be nonzero")
            d1 = self.monomial_coefficients()
            d2 = {}
            for mu in d1:
                d2[k*mu]=d1[mu]
            return self.parent()._from_dict(d2)

        def shift(self, mu):
            """
            Add `\mu` to any weight. Extended by linearity to the weight ring.

            INPUT:

            - ``mu`` -- a weight

            EXAMPLES::

                sage: g2 = WeylCharacterRing("G2",style="coroots").ambient()
                sage: [g2(1,2).shift(fw) for fw in g2.fundamental_weights()]
                [g2(2,2), g2(1,3)]
            """
            d1 = self.monomial_coefficients()
            d2 = {}
            for nu in d1:
                d2[mu+nu]=d1[nu]
            return self.parent()._from_dict(d2)

        def demazure(self, w, debug=False):
            r"""
            Return the result of applying the Demazure operator `\partial_w`
            to ``self``.

            INPUT:

            - ``w`` -- a Weyl group element, or its reduced word

            If `w = s_i` is a simple reflection, the operation `\partial_w`
            sends the weight `\lambda` to

            .. MATH::

                \frac{\lambda - s_i \cdot \lambda + \alpha_i}{1 + \alpha_i}

            where the numerator is divisible the denominator in the weight
            ring. This is extended by multiplicativity to all `w` in the
            Weyl group.

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style="coroots")
                sage: b2=WeightRing(B2)
                sage: b2(1,0).demazure([1])
                b2(1,0) + b2(-1,2)
                sage: b2(1,0).demazure([2])
                b2(1,0)
                sage: r=b2(1,0).demazure([1,2]); r
                b2(1,0) + b2(-1,2)
                sage: r.demazure([1])
                b2(1,0) + b2(-1,2)
                sage: r.demazure([2])
                b2(0,0) + b2(1,0) + b2(1,-2) + b2(-1,2)
            """
            if type(w) is list:
                word = w
            else:
                word = w.reduced_word()
            d1 = self.monomial_coefficients()
            d = {}
            alphacheck = self.parent()._space.simple_coroots()
            for v in d1:
                d[tuple(v.inner_product(alphacheck[j]) for j in self.parent().space().index_set())]=d1[v]
            return self.parent()._from_dict(self.parent().parent()._demazure_helper(d, word, debug=debug))

        def demazure_lusztig(self, i, v):
            r"""
            Return the result of applying the Demazure-Lusztig operator
            `T_i` to ``self``.

            INPUT:

            - ``i`` -- an element of the index set (or a reduced word or
              Weyl group element)
            - ``v`` -- an element of the base ring

            If `R` is the parent WeightRing, the Demazure-Lusztig operator
            `T_i` is the linear map `R \to R` that sends (for a weight
            `\lambda`) `R(\lambda)` to

            .. MATH::

                (R(\alpha_i)-1)^{-1} \bigl(R(\lambda) - R(s_i\lambda)
                - v(R(\lambda) - R(\alpha_i + s_i \lambda)) \bigr)

            where the numerator is divisible by the denominator in `R`.
            The Demazure-Lusztig operators give a representation of the
            Iwahori--Hecke algebra associated to the Weyl group. See

            * Lusztig, Equivariant `K`-theory and representations of Hecke
              algebras, Proc. Amer. Math. Soc. 94 (1985), no. 2, 337-342.
            * Cherednik, *Nonsymmetric Macdonald polynomials*. IMRN 10,
              483-515 (1995).

            In the examples, we confirm the braid and quadratic relations
            for type `B_2`.

            EXAMPLES::

                sage: P.<v> = PolynomialRing(QQ)
                sage: B2 = WeylCharacterRing("B2",style="coroots",base_ring=P); b2 = B2.ambient()
                sage: def T1(f) : return f.demazure_lusztig(1,v)
                sage: def T2(f) : return f.demazure_lusztig(2,v)
                sage: T1(T2(T1(T2(b2(1,-1)))))
                (v^2-v)*b2(0,-1) + v^2*b2(-1,1)
                sage: [T1(T1(f))==(v-1)*T1(f)+v*f for f in [b2(0,0), b2(1,0), b2(2,3)]]
                [True, True, True]
                sage: [T1(T2(T1(T2(b2(i,j))))) == T2(T1(T2(T1(b2(i,j))))) for i in [-2..2] for j in [-1,1]]
                [True, True, True, True, True, True, True, True, True, True]

            Instead of an index `i` one may use a reduced word or
            Weyl group element::

                sage: b2(1,0).demazure_lusztig([2,1],v)==T2(T1(b2(1,0)))
                True
                sage: W = B2.space().weyl_group(prefix="s")
                sage: [s1,s2]=W.simple_reflections()
                sage: b2(1,0).demazure_lusztig(s2*s1,v)==T2(T1(b2(1,0)))
                True
            """
            if i in self.parent().space().index_set():
                rho = self.parent().space().from_vector_notation(self.parent().space().rho(),style="coroots")
                inv = self.scale(-1)
                return (-inv.shift(-rho).demazure([i]).shift(rho)+v*inv.demazure([i])).scale(-1)
            elif type(i) is list:
                if len(i) == 0:
                    return self
                elif len(i) == 1:
                    return self.demazure_lusztig(i[0],v)
                else:
                    return self.demazure_lusztig(i[1:],v).demazure_lusztig(i[:1],v)
            else:
                try:
                    return self.demazure_lusztig(i.reduced_word(),v)
                except StandardError:
                    raise ValueError("unknown index {}".format(i))

