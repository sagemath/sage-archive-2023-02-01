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

    Let K be a compact Lie group, which we assume is semisimple and
    simply-connected. Its complexified Lie algebra L is the Lie algebra of a
    complex analytic Lie group G. The following three categories are
    equivalent: finite-dimensional representations of K; finite-dimensional
    representations of L; and finite-dimensional analytic representations of
    G. In every case, there is a parametrization of the irreducible
    representations by their highest weight vectors. For this theory of Weyl,
    see (for example):

    * J. F. Adams, Lectures on Lie groups
    * Broecker and Tom Dieck, Representations of Compact Lie groups
    * Bump, Lie Groups
    * Fulton and Harris, Representation Theory
    * Goodman and Wallach, Representations and Invariants of the Classical Groups
    * Hall, Lie Groups, Lie Algebras and Representations
    * Humphreys, Introduction to Lie Algebras and their representations
    * Procesi, Lie Groups
    * Samelson, Notes on Lie Algebras
    * Varadarajan, Lie groups, Lie algebras, and their representations
    * Zhelobenko, Compact Lie Groups and their Representations.

    Computations that you can do with these include computing their
    weight multiplicities, products (thus decomposing the tensor
    product of a representation into irreducibles) and branching
    rules (restriction to a smaller group).

    There is associated with K, L or G as above a lattice, the weight
    lattice, whose elements (called weights) are characters of a Cartan
    subgroup or subalgebra. There is an action of the Weyl group W on
    the lattice, and elements of a fixed fundamental domain for W, the
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

    Here R(1), R(fw1) and R(fw2) are irreducible representations with
    highest weight vectors 0, and the first two fundamental weights.

    For type A (also `G_2`, `F_4`, `E_6` and `E_7`) we will take as the weight
    lattice not the weight lattice of the semisimple group, but for a
    larger one. For type A, this means we are concerned with the
    representation theory of K=U(n) or G=GL(n,CC) rather than SU(n) or
    SU(n,CC). This is useful since the representation theory of GL(n)
    is ubiquitous, and also since we may then represent the fundamental
    weights (in root_system.py) by vectors with integer entries. If
    you are only interested in SL(3), say, use
    WeylCharacterRing(['A',2]) as above but be aware that R([a,b,c])
    and R([a+1,b+1,c+1]) represent the same character of SL(3) since
    R([1,1,1]) is the determinant.

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
        EXAMPLES:

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
        Expands the basis element indexed by the weight ``irr`` into
        the weight ring of ``self``.

        This is used to implement :meth:`lift`.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([2,1,0]); v
            (2, 1, 0)
            sage: A2.lift_on_basis(v)
            2*a2(1,1,1) + a2(1,2,0) + a2(1,0,2) + a2(2,1,0) + a2(2,0,1) + a2(0,1,2) + a2(0,2,1)

        This is consistent with the analoguous calculation with symmetric schur functions::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s[2,1].expand(3)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
        """
        return self.ambient()._from_dict(self._irr_weights(irr))

    @lazy_attribute
    def lift(self):
        """
        The embedding of ``self`` into its weight ring

        EXAMPLES ::

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

        EXAMPLES ::

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
        The partial inverse map from the weight ring into ``self``

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
        return "The Weyl Character Ring of Type %s with %s coefficients"%(self._cartan_type.__repr__(), self._base_ring.__repr__())

    def __call__(self, *args):
        """
        Construct an element of ``self``

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
        Construct a monomial from a weight

        INPUT:

        -  ``weight`` -- an element of the weight space, or a tuple

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
            raise ValueError, "%s is not a dominant element of the weight lattice"%weight
        return self.monomial(weight)

    def product_on_basis(self, a, b):
        """
        EXAMPLES::

            sage: D4 = WeylCharacterRing(['D',4])
            sage: spin_plus = D4(1/2,1/2,1/2,1/2)
            sage: spin_minus = D4(1/2,1/2,1/2,-1/2)
            sage: spin_plus*spin_minus
            D4(1,0,0,0) + D4(1,1,1,0)
            sage: spin_minus*spin_plus
            D4(1,0,0,0) + D4(1,1,1,0)
        """
        d = {}
        d1 = self._irr_weights(a)
        d2 = self._irr_weights(b)
        for k in d1:
            for l in d2:
                m = k+l
                if m in d:
                    d[m] += d1[k]*d2[l]
                else:
                    d[m] = d1[k]*d2[l]
        for k in d.keys():
            if d[k] == 0:
                del d[k]
        return self.char_from_weights(d)

    def some_elements(self):
        """
        EXAMPLE::

            sage: WeylCharacterRing("A3").some_elements()
            [A3(1,0,0,0), A3(1,1,0,0), A3(1,1,1,0)]

        """
        return [self.monomial(x) for x in self.fundamental_weights()]

    def one_basis(self):
        """
        EXAMPLE::

            sage: WeylCharacterRing("A3").one_basis()
            (0, 0, 0, 0)
            sage: WeylCharacterRing("A3").one()
            A3(0,0,0,0)
        """
        return self._space.zero()

    @cached_method
    def _irr_weights(self, hwv):
        """
        Given a dominant weight hwv, produce the dictionary of
        weight multiplicities for the irreducible representation
        with highest weight vector hwv. This method is cached
        for efficiency.

        EXAMPLE::

            sage: A2=WeylCharacterRing("A2")
            sage: v = A2.fundamental_weights()[1]; v
            (1, 0, 0)
            sage: A2._irr_weights(v)
            {(0, 1, 0): 1, (1, 0, 0): 1, (0, 0, 1): 1}
        """
        return irreducible_character_freudenthal(hwv, self._space)

    @cached_method
    def _weight_multiplicities(self, x):
        """
        Produces the dictionary of weight multiplicities for the
        WeylCharacter x. The character x does not have to be irreducible.

        EXAMPLE::

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
        Returns the base ring.

        EXAMPLES::

            sage: R = WeylCharacterRing(['A',3], base_ring = CC); R.base_ring()
            Complex Field with 53 bits of precision
        """
        return self._base_ring

    def irr_repr(self, hwv):
        """
        Return a string representing the irreducible character with highest
        weight vector hwv.

        EXAMPLES::

            sage: B3 = WeylCharacterRing("B3")
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(1,1,0)', 'B3(1/2,1/2,1/2)']
            sage: B3 = WeylCharacterRing("B3", style="coroots")
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(0,1,0)', 'B3(0,0,1)']
        """
        if self._style == "lattice":
            vec = hwv.to_vector()
        elif self._style == "coroots":
            vec = [hwv.inner_product(x) for x in self.simple_coroots()]
        else:
            raise ValueError, "unknown style"
        hstring = str(vec[0])
        for i in range(1,len(vec)):
            hstring=hstring+","+str(vec[i])
        return self._prefix+"("+hstring+")"

    def _repr_term(self, t):
        """
        Representation of the monomial corresponding to a weight t.

        EXAMPLES ::

            sage: G2=WeylCharacterRing("G2")
            sage: [G2._repr_term(x) for x in G2.fundamental_weights()]
            ['G2(1,0,-1)', 'G2(2,-1,-1)']
        """
        return self.irr_repr(t)

    def cartan_type(self):
        """
        Returns the Cartan type.

        EXAMPLES::

            sage: WeylCharacterRing("A2").cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def fundamental_weights(self):
        """
        Returns the fundamental weights.

        EXAMPLES::

            sage: WeylCharacterRing("G2").fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Returns the simple roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def simple_coroots(self):
        """
        Returns the simple coroots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_coroots()

    def highest_root(self):
        """
        Returns the highest_root.

        EXAMPLES::

            sage: WeylCharacterRing("G2").highest_root()
             (2, -1, -1)
        """
        return self._space.highest_root()

    def positive_roots(self):
        """
        Returns the positive roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def dynkin_diagram(self):
        """
        Returns the Dynkin diagram.

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
        Returns the extended Dynkin diagram, which is the Dynkin diagram
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
        Returns the rank.

        EXAMPLES::

            sage: WeylCharacterRing("G2").rank()
            2
        """
        return self._rank

    def space(self):
        """
        Returns the weight space associated to self.

        EXAMPLES::

            sage: WeylCharacterRing(['E',8]).space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def char_from_weights(self, mdict):
        """
        Construct a Weyl character from an invariant linear combination of weights

        INPUT:

         - ``mdict`` -- a dictionary mapping weights to coefficients,
           and representing a linear combination of weights which
           shall be invariant under the action of the Weyl group

        OUTPUT: the corresponding Weyl character

        EXAMPLES ::

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
        Helper method for char_from_weights. The output of this method is a
        dictionary whose keys are dominant weights that is the same as the
        monomial_coefficients method of self.char_from_weights().

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([3,1,0])
            sage: d = dict([(x,1) for x in v.orbit()])
            sage: A2._char_from_weights(d)
            {(3, 1, 0): 1, (2, 1, 1): -1, (2, 2, 0): -1}
        """
        hdict = {}
        ddict = mdict.copy()
        while not ddict == {}:
            highest = max((x.inner_product(self._space.rho()),x) for x in ddict)[1]
            if not highest.is_dominant():
                raise ValueError, "multiplicity dictionary may not be Weyl group invariant"
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
        A class for Weyl Characters.
        """
        def cartan_type(self):
            """
            Returns the Cartan type.

            EXAMPLES::

                sage: A2=WeylCharacterRing("A2")
                sage: A2([1,0,0]).cartan_type()
                ['A', 2]
            """
            return self.parent()._cartan_type

        def degree(self):
            """
            The degree of the character, that is, the dimension of module.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: [B3(x).degree() for x in B3.fundamental_weights()]
                [7, 21, 8]
            """
            L = self.parent()._space
            return sum(L.weyl_dimension(k)*c for k,c in self)

        def branch(self, S, rule="default"):
            """
            Returns the restriction of the character to the subalgebra. If no
            rule is specified, we will try to specify one.

            INPUT:

             - ``S`` - a Weyl character ring for a Lie subgroup or
               subalgebra

             -  ``rule`` - a branching rule.

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

        def is_irreducible(self):
            """
            Returns whether ``self`` is an irreducible character.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: [B3(x).is_irreducible() for x in B3.fundamental_weights()]
                [True, True, True]
                sage: sum(B3(x) for x in B3.fundamental_weights()).is_irreducible()
                False
            """
            return self.coefficients() == [1]

        def symmetric_square(self):
            """
            Returns the symmetric square of the character.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2",style="coroots")
                sage: A2(1,0).symmetric_square()
                A2(2,0)
            """
            # Conceptually, this converts self to the weight ring,
            # computes its square there, and convert the result back,
            # similar to what is done by product_on_basis. So in
            # principle, one could just return self^2.
            #
            # This implementation uses that this is a squaring (and not
            # a generic product) in the weight ring to optimize by
            # running only through pairs of weights instead of couples.
            # Worth it?
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
            Returns the exterior square of the character.

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
            Returns:

             - `1` if the representation is real (orthogonal)

             - `-1` if the representation is quaternionic (symplectic)

             - `0` if the representation is complex (not self dual)

            The Frobenius-Schur indicator of a character 'chi'
            of a compact group G is the Haar integral over the
            group of 'chi(g^2)'. Its value is 1,-1 or 0. This
            method computes it for irreducible characters of
            compact Lie groups by checking whether the symmetric
            and exterior square characters contain the trivial
            character.

            TODO: Try to compute this directly without actually calculating
            the full symmetric and exterior squares.

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style="coroots")
                sage: B2(1,0).frobenius_schur_indicator()
                 1
                sage: B2(0,1).frobenius_schur_indicator()
                 -1
            """
            if not self.is_irreducible():
                raise ValueError, "Frobenius-Schur indicator is only valid for irreducible characters"
            z = self.parent()._space.zero()
            if self.symmetric_square().coefficient(z) != 0:
                return 1
            if self.exterior_square().coefficient(z) != 0:
                return -1
            return 0

        def weight_multiplicities(self):
            """
            Produces the dictionary of weight multiplicities for the
            WeylCharacter self. The character does not have to be irreducible.

            EXAMPLE ::

                sage: B2=WeylCharacterRing("B2",style="coroots")
                sage: B2(0,1).weight_multiplicities()
                {(-1/2, 1/2): 1, (-1/2, -1/2): 1, (1/2, -1/2): 1, (1/2, 1/2): 1}
            """
            return self.parent()._weight_multiplicities(self)

        def inner_product(self, other):
            """
            Computes the inner product with another character. The
            irreducible characters are an orthonormal basis with
            respect to the usual inner product of characters,
            interpreted as functions on a compact Lie group,
            by Schur orthogonality.

            EXAMPLES ::

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
            Returns the multiplicity of the trivial representation
            in self. Multiplicities of other irreducibles may be
            obtained using the method ``multiplicity.``

            EXAMPLE ::

                sage: A2 = WeylCharacterRing("A2",style="coroots")
                sage: rep = A2(1,0)^2*A2(0,1)^2; rep
                2*A2(0,0) + A2(0,3) + 4*A2(1,1) + A2(3,0) + A2(2,2)
                sage: rep.invariant_degree()
                2
            """
            return self.coefficient(self.parent().space()(0))

        def multiplicity(self, other):
            """
            INPUT:

            - ``other`` - an irreducible character.

            Returns the multiplicity of the irreducible other in self.

            EXAMPLE::

                sage: B2 = WeylCharacterRing("B2",style="coroots")
                sage: rep = B2(1,1)^2; rep
                B2(0,0) + B2(1,0) + 2*B2(0,2) + B2(2,0) + 2*B2(1,2) + B2(0,4) + B2(3,0) + B2(2,2)
                sage: rep.multiplicity(B2(0,2))
                2
            """
            if not other.is_irreducible():
                raise ValueError, "%s is not irreducible"%other
            return self.coefficient(other.support()[0])

def irreducible_character_freudenthal(hwv, L, debug=False):
    """
    Returns the dictionary of multiplicities for the irreducible
    character with highest weight lamb. The weight multiplicities are
    computed by the Freudenthal multiplicity formula. The algorithm is
    based on recursion relation that is stated, for example, in
    Humphrey's book on Lie Algebras. The multiplicities are invariant
    under the Weyl group, so to compute them it would be sufficient to
    compute them for the weights in the positive Weyl chamber. However
    after some testing it was found to be faster to compute every
    weight using the recursion, since the use of the Weyl group is
    expensive in its current implementation.

    INPUT:

    - ``hwv`` - a dominant weight in a weight lattice.

    - ``L`` - the ambient space
    """
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
    A Branching rule describes the restriction of representations from
    a Lie group or algebra G to a smaller one. See for example, R. C.
    King, Branching rules for classical Lie groups using tensor and
    spinor methods. J. Phys. A 8 (1975), 429-449, Howe, Tan and
    Willenbring, Stable branching rules for classical symmetric pairs,
    Trans. Amer. Math. Soc. 357 (2005), no. 4, 1601-1626, McKay and
    Patera, Tables of Dimensions, Indices and Branching Rules for
    Representations of Simple Lie Algebras (Marcel Dekker, 1981),
    and Fauser, Jarvis, King and Wybourne, New branching rules induced
    by plethysm. J. Phys. A 39 (2006), no. 11, 2611--2655.

    INPUT:

    - ``chi`` - a character of G

    - ``R`` - the Weyl Character Ring of G

    - ``S`` - the Weyl Character Ring of H

    - ``rule`` - a set of r dominant weights in H where r is the rank
      of G.

    You may use a predefined rule by specifying rule = one of"levi",
    "automorphic", "symmetric", "extended", "orthogonal_sum", "tensor",
    "triality" or  "miscellaneous". The use of these rules will be
    explained next. After the examples we will explain how to write
    your own branching rules for cases that we have omitted.

    To explain the predefined rules we survey the most important
    branching rules. These may be classified into several cases, and
    once this is understood, the detailed classification can be read
    off from the Dynkin diagrams. Dynkin classified the maximal
    subgroups of Lie groups in Mat. Sbornik N.S. 30(72):349-462
    (1952).

    We will list give predefined rules that cover most cases where the
    branching rule is to a maximal subgroup. For convenience, we
    also give some branching rules to subgroups that are not maximal.
    For example, a Levi subgroup may or may not be maximal.

    LEVI TYPE. These can be read off from the Dynkin diagram. If
    removing a node from the Dynkin diagram produces another Dynkin
    diagram, there is a branching rule. Currently we require that the
    smaller diagram be connected. For these rules use the option
    rule="levi"::

       ['A',r] => ['A',r-1]
       ['B',r] => ['A',r-1]
       ['B',r] => ['B',r-1]
       ['C',r] => ['A',r-1]
       ['C',r] => ['C',r-1]
       ['D',r] => ['A',r-1]
       ['D',r] => ['D',r-1]
       ['E',r] => ['A',r-1] r = 7,8
       ['E',r] => ['D',r-1] r = 6,7,8
       ['E',r] => ['E',r-1]
       F4 => B3
       F4 => C3
       G2 => A1 (short root)

    Not all Levi subgroups are maximal subgroups. If the Levi is not
    maximal there may or may not be a preprogrammed rule="levi" for
    it. If there is not, the branching rule may still be obtained by going
    through an intermediate subgroup that is maximal using rule="extended".
    Thus the other Levi branching rule from G2 => A1 corresponding to the
    long root is available by first branching G2 => A_2 then A2 => A1.
    Similarly the branching rules to the Levi subgroup::

       ['E',r] => ['A',r-1] r = 6,7,8

    may be obtained by first branching E6=>A5xA1, E7=>A7 or E8=>A8.

    AUTOMORPHIC TYPE. If the Dynkin diagram has a symmetry, then there
    is an automorphism that is a special case of a branching rule.
    There is also an exotic "triality" automorphism of D4 having order
    3. Use rule="automorphic" or (for D4) rule="triality"::

        ['A',r] => ['A',r]
        ['D',r] => ['D',r]
        E6 => E6

    SYMMETRIC TYPE. Related to the automorphic type, when G admits
    an outer automorphism (usually of degree 2) we may consider
    the branching rule to the isotropy subgroup H. In many cases
    the Dynkin diagram of H can be obtained by folding the Dynkin
    diagram of G. For such isotropy subgroups use rule="symmetric".
    The last branching rule, D4=>G2 is not to a maximal subgroup
    since D4=>B3=>G2, but it is included for convenience. ::

        ['A',2r+1] => ['B',r]
        ['A',2r] => ['C',r]
        ['A',2r] => ['D',r]
        ['D',r] => ['B',r-1]
        E6 => F4
        D4 => G2

    EXTENDED TYPE. If removing a node from the extended Dynkin diagram
    results in a Dynkin diagram, then there is a branching rule. Use
    rule="extended" for these. We will also use this classification
    for some rules that are not of this type, mainly involving type B,
    such as D6 => B3xB3.

    Here is the extended Dynkin diagram for D6::

            0       6
            O       O
            |       |
            |       |
        O---O---O---O---O
        1   2   3   4   6

    Removing the node 3 results in an embedding D3xD3 -> D6. This
    corresponds to the embedding SO(6)xSO(6) -> SO(12), and is of
    extended type. On the other hand the embedding SO(5)xSO(7)-->SO(12)
    (e.g. B2xB3 -> D6) cannot be explained this way but for
    uniformity is implemented under rule="extended".

    The following rules are implemented as special cases of rule="extended". ::

        E6 => A5xA1, A2xA2xA2
        E7 => A7, D6xA1, A3xA3xA1
        E8 => A8, D8, E7xA1, A4xA4, D5xA3, E6xA2
        F4 => B4, C3xA1, A2xA2, A3xA1
        G2 => A1xA1

    Note that E8 has only a limited number of representations of reasonably low
    degree.

    ORTHOGONAL_SUM:

    Using rule="orthogonal_sum" you can get any branching rule
    SO(n) => SO(a) x SO(b) x SO(c) x ... where n = a+b+c+ ...
    Sp(2n) => Sp(2a) x Sp(2b) x Sp(2c) x ... where n = a+b+c+ ...
    where O(a) = ['D',r] (a=2r) or ['B',r] (a=2r+1)
    and Sp(2r)=['C',r]. In some cases these are also of
    extended type, as in the case ``D3xD3->D6`` discussed above.
    But in other cases, for example ``B3xB3->D7``, they are not
    of extended type.

    TENSOR: There are branching rules:
    ::

        ['A', rs-1] => ['A',r-1] x ['A',s-1]
        ['B',2rs+r+s] => ['B',r] x ['B',s]
        ['D',2rs+s] => ['B',r] x ['D',s]
        ['D',2rs] => ['D',r] x ['D',s]
        ['D',2rs] => ['C',r] x ['C',s]
        ['C',2rs+s] => ['B',r] x ['C',s]
        ['C',2rs] => ['C',r] x ['D',s].

    corresponding to the tensor product homomorphism. For type
    A, the homomorphism is GL(r) x GL(s) -> GL(rs). For the
    classical types, the relevant fact is that if V,W are
    orthogonal or symplectic spaces, that is, spaces endowed
    with symmetric or skew-symmetric bilinear forms, then V
    tensor W is also an orthogonal space (if V and W are both
    orthogonal or both symplectic) or symplectic (if one of
    V and W is orthogonal and the other symplectic).

    The corresponding branching rules are obtained using rule="tensor".

    SYMMETRIC POWER: The k-th symmetric and exterior power homomorphisms
    map GL(n) --> GL(binomial(n+k-1,k)) and GL(binomial(n,k)). The
    corresponding branching rules are not implemented but a special
    case is. The k-th symmetric power homomorphism SL(2) --> GL(k+1)
    has its image inside of SO(2r+1) if k=2r and inside of Sp(2r) if
    k=2r-1. Hence there are branching rules::

        ['B',r] => A1
        ['C',r] => A1

    and these may be obtained using the rule "symmetric_power".

    MISCELLANEOUS: Use rule="miscellaneous" for the following rules::

        B3 => G2
        F4 => G2xA1 (not implemented yet)

    BRANCHING RULES FROM PLETHYSMS

    Nearly all branching rules G => H where G is of type A,B,C or D
    are covered by the preceding rules. The function
    branching_rules_from_plethysm covers the remaining cases.

    ISOMORPHIC TYPE: Although not usually referred to as a branching
    rule, the effects of the accidental isomorphisms may be handled
    using rule="isomorphic"::

        B2 => C2
        C2 => B2
        A3 => D3
        D3 => A3
        D2 => A1xA1
        B1 => A1
        C1 => A1

    EXAMPLES: (Levi type)

    ::

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
    of B3 being branched is spin, which is not a representation of
    SO(7) but of its double cover spin(7). The group A2 is really GL(3)
    and the double cover of SO(7) induces a cover of GL(3) that is
    trivial over SL(3) but not over the center of GL(3). The weight
    lattice for this GL(3) consists of triples (a,b,c) of half integers
    such that a-b and b-c are in `\ZZ`, and this is reflected in the last
    decomposition.

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

    EXAMPLES: (Automorphic type, including D4 triality)

    ::

        sage: [A3(chi).branch(A3,rule="automorphic") for chi in A3.fundamental_weights()]
        [A3(0,0,0,-1), A3(0,0,-1,-1), A3(0,-1,-1,-1)]
        sage: [D4(chi).branch(D4,rule="automorphic") for chi in D4.fundamental_weights()]
        [D4(1,0,0,0), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1/2,1/2,1/2,-1/2)]
        sage: [D4(chi).branch(D4,rule="triality") for chi in D4.fundamental_weights()]
        [D4(1/2,1/2,1/2,-1/2), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1,0,0,0)]

    EXAMPLES: (Symmetric type)

    ::

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

    EXAMPLES: (Extended type)

    ::

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


    EXAMPLES: (Tensor type)

    ::

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

    EXAMPLES: (Symmetric Power)

    ::

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: B3=WeylCharacterRing("B3",style="coroots")
        sage: C3=WeylCharacterRing("C3",style="coroots")
        sage: [B3(fw).branch(A1,rule="symmetric_power") for fw in B3.fundamental_weights()]
        [A1(6), A1(2) + A1(6) + A1(10), A1(0) + A1(6)]
        sage: [C3(fw).branch(A1,rule="symmetric_power") for fw in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    EXAMPLES: (Miscellaneous type)

    ::

        sage: G2 = WeylCharacterRing("G2")
        sage: [fw1, fw2, fw3] = B3.fundamental_weights()
        sage: B3(fw1+fw3).branch(G2, rule="miscellaneous")
        G2(1,0,-1) + G2(2,-1,-1) + G2(2,0,-2)

    EXAMPLES: (Isomorphic type)

    ::

        sage: [B2(x).branch(C2, rule="isomorphic") for x in B2.fundamental_weights()]
        [C2(1,1), C2(1,0)]
        sage: [C2(x).branch(B2, rule="isomorphic") for x in C2.fundamental_weights()]
        [B2(1/2,1/2), B2(1,0)]
        sage: [A3(x).branch(D3,rule="isomorphic") for x in A3.fundamental_weights()]
        [D3(1/2,1/2,1/2), D3(1,0,0), D3(1/2,1/2,-1/2)]
        sage: [D3(x).branch(A3,rule="isomorphic") for x in D3.fundamental_weights()]
        [A3(1/2,1/2,-1/2,-1/2), A3(1/4,1/4,1/4,-3/4), A3(3/4,-1/4,-1/4,-1/4)]

    Here A3(x,y,z,w) can be understood as a representation of SL(4).
    The weights x,y,z,w and x+t,y+t,z+t,w+t represent the same
    representation of SL(4) - though not of GL(4) - since
    A3(x+t,y+t,z+t,w+t) is the same as A3(x,y,z,w) tensored with
    `det^t`. So as a representation of SL(4),
    A3(1/4,1/4,1/4,-3/4) is the same as A3(1,1,1,0). The exterior
    square representation SL(4) -> GL(6) admits an invariant symmetric
    bilinear form, so is a representation SL(4) -> SO(6) that lifts to
    an isomorphism SL(4) -> Spin(6). Conversely, there are two
    isomorphisms SO(6) -> SL(4), of which we've selected one.

    In cases like this you might prefer style="coroots".

    ::

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

    EXAMPLES: (Branching rules from plethysms)

    This is a general rule that includes any branching rule
    from types A,B,C or D as a special case. Thus it could be
    used in place of the above rules and would give the same
    results. However it is most useful when branching from G
    to a maximal subgroup H such that rank(H) < rank(G)-1.

    We consider a homomorphism H --> G where G is one of
    SL(r+1), SO(2r+1), Sp(2r) or SO(2r). The function
    branching_rule_from_plethysm produces the corresponding
    branching rule. The main ingredient is the character
    chi of the representation of H that is the homomorphism
    to GL(r+1), GL(2r+1) or GL(2r).

    This rule is so powerful that it contains the other
    rules implemented above as special cases. First let
    us consider the symmetric fifth power representation
    of SL(2).

    ::

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: chi=A1([5])
        sage: chi.degree()
         6
        sage: chi.frobenius_schur_indicator()
        -1

    This confirms that the character has degree 6 and
    is symplectic, so it corresponds to a homomorphism
    SL(2) --> Sp(6), and there is a corresponding
    branching rule C3 => A1.

    ::

        sage: C3 = WeylCharacterRing("C3",style="coroots")
        sage: sym5rule = branching_rule_from_plethysm(chi,"C3")
        sage: [C3(hwv).branch(A1,rule=sym5rule) for hwv in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    This is identical to the results we would obtain using
    rule="symmetric_power". The next example gives a branching
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

    We have confirmed that the adjoint representation of G2
    gives a homomorphism into SO(14), and that the pullback
    of the one of the two 64 dimensional spin representations
    to SO(14) is an irreducible representation of G2.

    BRANCHING FROM A REDUCIBLE ROOT SYSTEM

    If you are branching from a reducible root system, the rule is
    a list of rules, one for each component type in the root system.
    The rules in the list are given in pairs [type, rule], where
    type is the root system to be branched to, and rule is the
    branching rule.

    ::

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

    WRITING YOUR OWN RULES

    Suppose you want to branch from a group G to a subgroup H.
    Arrange the embedding so that a Cartan subalgebra U of H is
    contained in a Cartan subalgebra T of G. There is thus
    a mapping from the weight spaces Lie(T)* --> Lie(U)*.
    Two embeddings will produce identical branching rules if they
    differ by an element of the Weyl group of H.

    The RULE is this map Lie(T)* = G.space() to Lie(U)* = H.space(),
    which you may implement as a function. As an example, let
    us consider how to implement the branching rule A3 => C2.
    Here H = C2 = Sp(4) embedded as a subgroup in A3 = GL(4). The
    Cartan subalgebra U consists of diagonal matrices with
    eigenvalues u1, u2, -u2, -u1. The C2.space() is the
    two dimensional vector spaces consisting of the linear
    functionals u1 and u2 on U. On the other hand Lie(T) is
    RR^4. A convenient way to see the restriction is to
    think of it as the adjoint of the map [u1,u2] -> [u1,u2,-u2,-u1],
    that is, [x0,x1,x2,x3] -> [x0-x3,x1-x2]. Hence we may
    encode the rule:

    ::

       def rule(x):
           return [x[0]-x[3],x[1]-x[2]]

    or simply:

    ::

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
    INPUT:

    - ``R`` - the Weyl Character Ring of G

    - ``S`` - the Weyl Character Ring of H

    - ``rule`` - a string describing the branching rule as a map from
      the weight space of S to the weight space of R.
    """
    r = Rtype.rank()
    s = Stype.rank()
    rdim = Rtype.root_system().ambient_space().dimension()
    sdim = Stype.root_system().ambient_space().dimension()
    if Stype.is_compound():
        stypes = Stype.component_types()
    if rule == "default":
        return lambda x : x
    elif rule == "levi":
        if not s == r-1:
            raise ValueError, "Incompatible ranks"
        if Rtype[0] == 'A':
            if Stype.is_compound():
                if all(ct[0]=='A' for ct in stypes) and rdim == sdim:
                    return lambda x : x
                else:
                    raise ValueError, "Rule not found"
            elif Stype[0] == 'A':
                return lambda x : list(x)[:r]
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] in ['B', 'C', 'D']:
            if Stype.is_atomic():
                if Stype[0] == 'A':
                    return  lambda x : x
                elif Stype[0] == Rtype[0]:
                    return lambda x : list(x)[1:]
            elif stypes[-1][0] == Rtype[0] and all(t[0] == 'A' for t in stypes[:-1]):
                return lambda x : x
            else:
                raise ValueError, "Rule not found"
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
                        raise NotImplementedError, """A5 Levi is not maximal. Branch to A5xA1 (rule="extended")."""
                    if r == 7:
                        raise NotImplementedError, """A5 Levi is not maximal. Branch to A5xA1 (rule="extended")."""
            raise NotImplementedError, "Not implemented yet"
        elif Rtype[0] == 'F' and s == 3:
            if Stype[0] == 'B':
                return lambda x : list(x)[1:]
            elif Stype[0] == 'C':
                return lambda x : [x[1]-x[0],x[2]+x[3],x[2]-x[3]]
            else:
                raise NotImplementedError, "Not implemented yet"
        elif Rtype[0] == 'G' and Stype[0] == 'A':
            return lambda x : list(x)[1:]
        else:
            raise ValueError, "Rule not found"
    elif rule == "automorphic":
        if not Rtype == Stype:
            raise ValueError, "Cartan types must agree for automorphic branching rule"
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
            raise ValueError, "No automorphism found"
    elif rule == "triality":
        if not Rtype == Stype:
            raise ValueError, "Triality is an automorphic type (for D4 only)"
        elif not Rtype[0] == 'D' and r == 4:
            raise ValueError, "Triality is for D4 only"
        else:
            return lambda x : [(x[0]+x[1]+x[2]+x[3])/2,(x[0]+x[1]-x[2]-x[3])/2,(x[0]-x[1]+x[2]-x[3])/2,(-x[0]+x[1]+x[2]-x[3])/2]
    elif rule == "symmetric":
        if Rtype[0] == 'A':
            if (Stype[0] == 'C' or Stype[0] == 'D' and r == 2*s-1) or (Stype[0] == 'B' and r == 2*s):
                return lambda x : [x[i]-x[r-i] for i in range(s)]
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'D' and Stype[0] == 'B' and s == r-1:
            return lambda x : x[:s]
        elif Rtype[0] == 'D' and r == 4 and Stype[0] == 'G':
            return lambda x : [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]]
        elif Rtype[0] == 'E' and Stype[0] == 'F' and r == 6 and s == 4:
            return lambda x : [(x[4]-3*x[5])/2,(x[0]+x[1]+x[2]+x[3])/2,(-x[0]-x[1]+x[2]+x[3])/2,(-x[0]+x[1]-x[2]+x[3])/2]
        else:
            raise ValueError, "Rule not found"
    elif rule == "extended" or rule == "orthogonal_sum":
        if rule == "extended" and not s == r:
            raise ValueError, """Ranks should be equal for rule="extended" """
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
                    raise ValueError, "Rule not found"
            elif Rtype[0] == 'C':
                if all(t[0] == Rtype[0] for t in stypes):
                    return lambda x : x
            if rule == "orthogonal_sum":
                raise ValueError, "Rule not found"
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
                            raise NotImplementedError, "Not maximal: first branch to E7xA1"
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
                                raise NotImplementedError, "Not maximal: first branch to A7xA1"
                    elif stypes[0][0] == 'D' and stypes[1][0] == 'A':
                        if stypes[0][1] == 5 and stypes[1][1] == 3:
                            raise NotImplementedError, "Not maximal: first branch to D8 then D5xD3=D5xA3"
                    elif stypes[0][0] == 'E' and stypes[1][0] == 'A':
                        if stypes[0][1] == 6 and stypes[1][1] == 2:
                            return lambda x : [x[0],x[1],x[2],x[3],x[4], \
                                               (x[5]+x[6]-x[7])/3,(x[5]+x[6]-x[7])/3,(-x[5]-x[6]+x[7])/3, \
                                               (-x[5]-x[6]-2*x[7])/3,(-x[5]+2*x[6]+x[7])/3,(2*x[5]-x[6]+x[7])/3]
                        elif stypes[0][1] == 7 and stypes[1][1] == 1:
                            return lambda x : [x[0],x[1],x[2],x[3],x[4],x[5],(x[6]-x[7])/2,(-x[6]+x[7])/2,(-x[6]-x[7])/2,(x[6]+x[7])/2]
                raise ValueError, "Rule not found"
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
                    raise ValueError, "Rule not found"
            elif Rtype[0] == 'G':
                if all(t[0] == 'A' and t[1] == 1 for t in stypes):
                    return lambda x : [(x[1]-x[2])/2,-(x[1]-x[2])/2, x[0]/2, -x[0]/2]
            else:
                raise ValueError, "Rule not found"
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
                raise ValueError, "Rule not found"
    elif rule == "isomorphic":
        if r != s:
            raise ValueError, "Incompatible ranks"
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
            raise ValueError, "Rule not found"
    elif rule == "tensor" or rule == "tensor-debug":
        if not Stype.is_compound():
            raise ValueError, "Tensor product requires more than one factor"
        if len(stypes) is not 2:
            raise ValueError, "Not implemented"
        if Rtype[0] is 'A':
            nr = Rtype[1]+1
        elif Rtype[0] is 'B':
            nr = 2*Rtype[1]+1
        elif Rtype[0] in ['C', 'D']:
            nr = 2*Rtype[1]
        else:
            raise ValueError, "Rule not found"
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
            raise ValueError, "Ranks don't agree with tensor product"
        if Rtype[0] == 'A':
            if all(t[0] == 'A' for t in stypes):
                def rule(x):
                    ret = [sum(x[i*ns[1]:(i+1)*ns[1]]) for i in range(ns[0])]
                    ret.extend([sum(x[ns[1]*j+i] for j in range(ns[0])) for i in range(ns[1])])
                    return ret
                return rule
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'B':
            if not all(t[0] == 'B' for t in stypes):
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'C':
            if stypes[0][0] in ['B','D'] and stypes[1][0] is 'C':
                pass
            elif stypes[1][0] in ['B','D'] and stypes[0][0] is 'C':
                pass
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'D':
            if stypes[0][0] in ['B','D'] and stypes[1][0] is 'D':
                pass
            elif stypes[1][0] is 'B' and stypes[0][0] is 'D':
                pass
            elif stypes[1][0] is 'C' and stypes[0][0] is 'C':
                pass
            else:
                raise ValueError, "Rule not found"
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
                raise ValueError, "Rule not found"
        else:
            raise ValueError, "Rule not found"
    elif rule == "miscellaneous":
        if Rtype[0] == 'B' and Stype[0] == 'G' and r == 3:
            return lambda x : [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]]
        else:
            raise ValueError, "Rule not found"


def branching_rule_from_plethysm(chi, cartan_type, return_matrix = False):
    """
    INPUT:

    - ``chi`` - the character of an irreducible representation pi of a group G
    - ``cartan_type`` - a classical Cartan type (A,B,C or D).

    It is assumed that the image of the irreducible representation pi
    naturally has its image in the group G.

    Returns a branching rule for this plethysm.

    EXAMPLE:

    The adjoint representation SL(3) --> GL(8) factors
    through SO(8). The branching rule in question will
    describe how representations of SO(8) composed with
    this homomorphism decompose into irreducible characters
    of SL(3)::

        sage: A2 = WeylCharacterRing("A2")
        sage: A2 = WeylCharacterRing("A2", style="coroots")
        sage: ad = A2(1,1)
        sage: ad.degree()
        8
        sage: ad.frobenius_schur_indicator()
        1

    This confirms that `ad` has degree 8 and is orthogonal,
    hence factors through SO(8)=D4::

        sage: br = branching_rule_from_plethysm(ad,"D4")
        sage: D4 = WeylCharacterRing("D4")
        sage: [D4(f).branch(A2,rule = br) for f in D4.fundamental_weights()]
        [A2(1,1), A2(0,3) + A2(1,1) + A2(3,0), A2(1,1), A2(1,1)]

    """
    ct = CartanType(cartan_type)
    if ct[0] not in ["A","B","C","D"]:
        raise ValueError, "not implemented for type %s"%ct[0]
    if ct[0] is "A":
        ret = []
        ml = chi.weight_multiplicities()
        for v in ml:
            n = ml[v]
            ret.extend(n*[v.to_vector()])
        M = matrix(ret).transpose()
        if len(M.columns()) != ct[1] + 1:
            raise ValueError, "representation has wrong degree for type %s"%ct.__repr__()
        return lambda x : tuple(M*vector(x))
    if ct[0] in ["B","D"]:
        if chi.frobenius_schur_indicator() != 1:
            raise ValueError, "character is not orthogonal"
    if ct[0] is "C":
        if chi.frobenius_schur_indicator() != -1:
            raise ValueError, "character is not symplectic"
    if ct[0] is "B":
        if is_even(chi.degree()):
            raise ValueError, "degree is not odd"
    if ct[0] is ["C","D"]:
        if is_odd(chi.degree()):
            raise ValueError, "degree is not even"
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
        raise ValueError, "representation has wrong degree for type %s"%ct.__repr__()
    if return_matrix:
        return M
    else:
        return lambda x : tuple(M*vector(x))

class WeightRing(CombinatorialFreeModule):
    """
    The WeightRing is the group algebra over a weight lattice.

    A WeylCharacter may be regarded as an element of the WeightRing.
    In fact, an element of the WeightRing is an element of the
    WeylCharacterRing if and only if it is invariant under the
    action of the Weyl Group.

    The advantage of the WeightRing over the WeylCharacterRing
    is that one may conduct calculations in the WeightRing that
    involve sums of weights that are not Weyl Group invariant.

    EXAMPLE::

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
            sage: WeightRing(G2)
            The Weight ring attached to The Weyl Character Ring of Type ['G', 2] with Univariate Polynomial Ring in q over Rational Field coefficients
        """
        return "The Weight ring attached to %s"%self._parent

    def __call__(self, *args):
        """
        Construct an element of ``self``

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
        Construct a monomial from a weight

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
        EXAMPLE::

            sage: A2=WeylCharacterRing("A2")
            sage: a2=WeightRing(A2)
            sage: a2(1,0,0)*a2(0,1,0)
            a2(1,1,0)
        """
        return self(a+b)

    def some_elements(self):
        """
        EXAMPLE::

            sage: A3=WeylCharacterRing("A3")
            sage: a3=WeightRing(A3)
            sage: a3.some_elements()
            [a3(1,0,0,0), a3(1,1,0,0), a3(1,1,1,0)]
        """
        return [self.monomial(x) for x in self.fundamental_weights()]

    def one_basis(self):
        """
        EXAMPLE::

            sage: A3=WeylCharacterRing("A3")
            sage: WeightRing(A3).one_basis()
            (0, 0, 0, 0)
            sage: WeightRing(A3).one()
            a3(0,0,0,0)
        """
        return self._space.zero()

    def parent(self):
        """
        Returns the parent Weyl Character Ring.

        EXAMPLE::

            sage: A2=WeylCharacterRing("A2")
            sage: a2=WeightRing(A2)
            sage: a2.parent()
            The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients
            sage: a2.parent()==A2
            True

        """
        return self._parent

    def weyl_character_ring(self):
        """
        Returns the parent Weyl Character Ring. A synonym for self.parent().

        EXAMPLE::

            sage: A2=WeylCharacterRing("A2")
            sage: a2=WeightRing(A2)
            sage: a2.weyl_character_ring()
            The Weyl Character Ring of Type ['A', 2] with Integer Ring coefficients
        """
        return self._parent

    def cartan_type(self):
        """
        Returns the Cartan type.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: WeightRing(A2).cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def space(self):
        """
        Returns the weight space realization associated to self.

        EXAMPLES::

            sage: E8 = WeylCharacterRing(['E',8])
            sage: e8 = WeightRing(E8)
            sage: e8.space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def fundamental_weights(self):
        """
        Returns the fundamental weights.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Returns the simple roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def positive_roots(self):
        """
        Returns the positive roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def wt_repr(self, wt):
        """
        Returns a string representing the irreducible character with
        highest weight vector wt.

        EXAMPLES::

            sage: G2 = WeylCharacterRing(['G',2])
            sage: g2 = WeightRing(G2)
            sage: g2.wt_repr([1,0,0])
            'g2(1,0,0)'
        """
        hstring = str(wt[0])
        for i in range(1,self._space.n):
            hstring=hstring+","+str(wt[i])
        return self._prefix+"("+hstring+")"

    def _repr_term(self, t):
        """
        Representation of the monomial corresponding to a weight t.

        EXAMPLES ::

            sage: G2=WeylCharacterRing("G2")
            sage: g2=WeightRing(G2)
            sage: [g2(x) for x in g2.fundamental_weights()]
            [g2(1,0,-1), g2(2,-1,-1)]
        """
        return self.wt_repr(t)

    class Element(CombinatorialFreeModule.Element):
        """
        A class for Weight Ring Elements.
        """
        def cartan_type(self):
            """
            Returns the Cartan type.
            EXAMPLES::

                sage: A2=WeylCharacterRing("A2")
                sage: a2 = WeightRing(A2)
                sage: a2([0,1,0]).cartan_type()
                ['A', 2]
            """
            return self.parent()._cartan_type

        def weyl_group_action(self, w):
            """
            Returns the action of the Weyl group element w on self.

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
            Assuming that self is invariant under the Weyl group, this will
            express it as a linear combination of characters. If self is not
            Weyl group invariant, this method will not terminate.

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


