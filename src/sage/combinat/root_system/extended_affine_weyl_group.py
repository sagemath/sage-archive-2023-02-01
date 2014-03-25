"""
Extended Affine Weyl Groups

AUTHORS:

- Daniel Bump (2012): initial version
- Daniel Orr (2012): initial version
- Anne Schilling (2012): initial version
- Mark Shimozono (2012): initial version
- Nicolas Thiery (2012): initial version
- Mark Shimozono (2013): twisted affine root systems, multiple realizations

"""

#*****************************************************************************
#       Copyright (C) 2012 Daniel Bump <bump at match.stanford.edu>,
#                     2012 Daniel Orr <danorr at live.unc.edu>
#                     2012 Anne Schilling <anne at math.ucdavis.edu>
#                     2012 Mark Shimozono <mshimo at math.vt.edu>
#                     2012 Nicolas Thiery <nthiery at users.sf.net>
#
#                     2013 Mark Shimozono
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.categories.groups import Groups
from sage.categories.sets_cat import Sets
from sage.categories.all import WeylGroups, FiniteWeylGroups, AffineWeylGroups
from sage.misc.cachefunc import cached_method
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family
from sage.categories.realizations import Category_realization_of_parent
from sage.misc.bindable_class import BindableClass
from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
from sage.misc.abstract_method import abstract_method
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.groups.group_exp import GroupExp
from sage.groups.group_semidirect_product import GroupSemidirectProduct
from sage.combinat.root_system.root_system import RootSystem

def ExtendedAffineWeylGroup(cartan_type, style="PW0", **print_options):
    r"""
    The extended affine Weyl group.

    INPUT:

        - ``cartan_type`` -- An affine or finite Cartan type (a finite Cartan type is an abbreviation for its untwisted affinization)
        - ``style`` -- one of "PW0", "W0P", "PvW0", "W0Pv", "WF", "FW" (default: "PW0"); this describes the implementation.
        - ``print_options`` -- Special instructions for printing elements

    The mnemonics for ``style`` are:

    - "P" -- subgroup of translations
    - "Pv" -- subgroup of translations in a dual form
    - "W0" -- classical Weyl group
    - "W" -- affine Weyl group
    - "F" -- fundamental group of length zero elements

    Hence "PW0" means the semidirect product of P with action from W0 with P
    written on the left. "W0P" is similar but with translations on the right.
    "WF" is the semidirect product of W under the action of F with W on the left.
    "FW" is similar but with W on the right.

    Recognized arguments for ``print_options`` are:
    
    - print_tuple -- True or False (default: False) If True, elements are printed `(a,b)`, otherwise as `a * b`.
    - affine -- Prefix for simple reflections in the affine Weyl group
    - classical -- Prefix for simple reflections in the classical Weyl group
    - translation -- Prefix for the translation elements
    - fundamental -- Prefix for the elements of the fundamental group

    Unlike :class:`CombinatorialFreeModule`, these options are not mutable.

    OUTPUT:

        - A realization of :class:`ExtendedAffineWeylGroup_Class`.

    The *extended affine Weyl group* was introduced in Iwahori, Generalized Tits system
    (Bruhat decomposition) on p-adic semisimple groups. 1966 Algebraic Groups and Discontinuous
    Subgroups (AMS Proc. Symp. Pure Math.., 1965) pp. 71-83 Amer. Math. Soc.,
    Providence, R.I. See Kac, Bourbaki, *Lie Groups and Lie Algebras* IV.2 and
    Kac, *Infinite-dimensional Lie algebras*.

    ..RUBRIC: Notation

    - `R` -- An irreducible affine root system
    - `I` -- Set of nodes of the Dynkin diagram of `R`
    - `R_0` -- The classical subsystem of `R`
    - `I_0` -- Set of nodes of the Dynkin diagram of `R_0`
    - `E` -- Extended affine Weyl group of type `R`
    - `W` -- Affine Weyl group of type `R`
    - `W_0` -- finite (classical) Weyl group (of type `R_0`)
    - `M` -- translation lattice for `W`
    - `L` -- translation lattice for `E`
    - `F` -- Fundamental subgroup of `E` (the length zero elements)
    - `P` -- Finite weight lattice
    - `Q` -- Finite root lattice
    - `P^\vee` -- Finite coweight lattice
    - `Q^\vee` -- Finite coroot lattice

    ..RUBRIC: Translation lattices

    The styles "PW0" and "W0P" use the following lattices.

    - Untwisted affine -- `L = P^\vee`, `M = Q^\vee`
    - Dual of untwisted affine -- `L = P`, `M = Q`
    - `BC_n` (`A_{2n}^{(2)}`) -- `L = M = P`
    - Dual of `BC_n` (`A_{2n}^{(2)\dagger}`) --  `L = M = P^\vee`

    The styles "PvW0" and "W0Pv" use the following lattices for `E`.

    - Untwisted affine -- The weight lattice of the dual finite Cartan type.
    - Dual untwisted affine -- The same choice as above

    For mixed affine type (`A_{2n}^{(2)}`, aka  `\tilde{BC}_n`, and their affine duals)
    the styles "PvW0" and "W0Pv" are not implemented.

    ..RUBRIC: Finite and affine Weyl groups `W_0` and `W`

    The finite Weyl group `W_0` is generated by the simple reflections `s_i` for `i \in I_0` where
    `s_i` is the reflection across a suitable hyperplane `H_i` through the origin in the
    real span `V` of the lattice `M`.

    `R` specifies another (affine) hyperplane `H_0`. The affine Weyl group `W` is generated by `W_0`
    and the reflection `S_0` across `H_0`.

    ..RUBRIC: Extended affine Weyl group `E`

    The complement in `V` of the set `H` of hyperplanes obtained from the `H_i` by the action of
    `W`, has connected components called alcoves. `W` acts freely and transitively on the set
    of alcoves. After the choice of a certain alcove (the fundamental alcove),
    there is an induced bijection from `W` to the set of alcoves under which the identity
    in `W` maps to the fundamental alcove.

    Then `L` is the largest sublattice of `V`, whose translations stabilize the set of alcoves.

    There are isomorphisms

    `W \cong M \rtimes W_0 \cong W_0 \ltimes M`

    `E \cong L \rtimes W_0 \cong W_0 \ltimes L`

    ..RUBRIC: Fundamental group of affine Dynkin automorphisms

    Since `L` acts on the set of alcoves, the group `F = L/M` may be viewed as a
    subgroup of the symmetries of the fundamental alcove or equivalently the
    symmetries of the affine Dynkin diagram.
    `F` acts on the set of alcoves and hence on `W`. Conjugation by an element of `F`
    acts on `W` by permuting the indices of simple reflections.

    There are isomorphisms

    `E \cong F \ltimes W \cong W \rtimes F`

    An affine Dynkin node is *special* if it is conjugate to the zero node under some
    affine Dynkin automorphism.

    There is a bijection `i` `\mapsto` `\pi_i` from the set of special nodes
    to the group `F`, where `\pi_i` is the unique element of `F` that sends `0` to `i`.
    When `L=P` (resp. `L=P^\vee`) the element `\pi_i` is induced
    (under the isomorphism `F \cong L/M`) by addition of the coset of the
    `i`-th fundamental weight (resp. coweight).

    The length function of the Coxeter group `W` may be extended to `E` by
    `\ell(w \pi) = \ell(w)` where `w \in W` and `\pi\in F`.
    This is the number of hyperplanes in `H` separating the
    fundamental alcove from its image by `w \pi` (or equivalently `w`).

    It is known that if `G` is the compact Lie group of adjoint type with root
    system `R_0` then `F` is isomorphic to the fundamental group of `G`, or
    to the center of its simply-connected covering group. That is why we
    call `F` the *fundamental group*.

    In the future we may want to build an element of the group from an appropriate linear map f
    on some of the root lattice realizations for this cartan type: W.from_endomorphism(f).

    EXAMPLES::

        sage: PW0 = ExtendedAffineWeylGroup(["A",2,1]); PW0
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Multiplicative form of Coweight lattice of the Root system of type ['A', 2] acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)

    Style "PW0" is the default realization of
    :class:`~sage.combinat.root_system.extended_affine_weyl_group.ExtendedAffineWeylGroup_Class`
    which is stored in the variable ``E`` below.

    ::

        sage: E = PW0.realization_of(); E
        Extended affine Weyl group of type ['A', 2, 1]

        sage: type(E)
        <class 'sage.combinat.root_system.extended_affine_weyl_group.ExtendedAffineWeylGroup_Class_with_category'>

    The other realizations "W0P", "PvW0", "W0Pv", "WF", and "FW" can be accessed from the others via ``E``.

    ::

        sage: W0P = E.W0P(); W0P
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice) acting on Multiplicative form of Coweight lattice of the Root system of type ['A', 2]

        sage: PvW0 = E.PvW0(); PvW0
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['A', 2] acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)

        sage: W0Pv = E.W0Pv(); W0Pv
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice) acting on Multiplicative form of Weight lattice of the Root system of type ['A', 2]

        sage: WF = E.WF(); WF
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice) acted upon by Fundamental group of type ['A', 2, 1]

        sage: FW = E.FW(); FW
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Fundamental group of type ['A', 2, 1] acting on Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)        

    When the realizations are constructed from each other as above, there are built-in coercions between them.

    ::

        sage: x = WF.an_element(); x
        S0*S1*S2 * pi[2]
        sage: FW(x)
        pi[2] * S1*S2*S0
        sage: W0P(x)
        s1*s2*s1 * t[-2*Lambdacheck[1] - Lambdacheck[2]]
        sage: PW0(x)
        t[Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2*s1
        sage: PvW0(x)
        t[Lambda[1] + 2*Lambda[2]] * s1*s2*s1

    The translation lattice and its distinguished basis are obtained from ``E``.

    ::

        sage: L = E.lattice(); L
        Coweight lattice of the Root system of type ['A', 2]
        sage: b = E.lattice_basis(); b
        Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}

    Translation lattice elements can be coerced into any realization.

    ::

        sage: PW0(b[1]-b[2])
        t[Lambdacheck[1] - Lambdacheck[2]]
        sage: FW(b[1]-b[2])
        pi[2] * S0*S1

    The dual form of the translation lattice and its basis are similarly obtained.

    ::

        sage: Lv = E.dual_lattice(); Lv
        Weight lattice of the Root system of type ['A', 2]
        sage: bv = E.dual_lattice_basis(); bv
        Finite family {1: Lambda[1], 2: Lambda[2]}
        sage: FW(bv[1]-bv[2])
        pi[2] * S0*S1

    The abstract fundamental group is accessed from ``E``.

    ::

        sage: F = E.fundamental_group(); F
        Fundamental group of type ['A', 2, 1]

    Its elements are indexed by the set of special nodes of the affine Dynkin diagram.

    ::

        sage: E.special_nodes()
        [0, 1, 2]
        sage: F.family()
        Finite family {0: pi[0], 1: pi[1], 2: pi[2]}

    There is a coercion from the fundamental group into each realization:

    ::

        sage: F(2)
        pi[2]
        sage: WF(F(2))
        pi[2]
        sage: W0P(F(2))
        s2*s1 * t[-Lambdacheck[1]]
        sage: W0Pv(F(2))
        s2*s1 * t[-Lambda[1]]

    Using ``E`` one may access the classical and affine Weyl groups and their morphisms
    into each realization.

    ::

        sage: W0 = E.classical_weyl(); W0
        Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
        sage: v = W0.from_reduced_word([1,2,1]); v
        s1*s2*s1
        sage: PW0(v)
        s1*s2*s1
        sage: WF(v)
        S1*S2*S1
        sage: W = E.affine_weyl(); W
        Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
        sage: w = W.from_reduced_word([2,1,0]); w
        S2*S1*S0
        sage: WF(w)
        S2*S1*S0
        sage: PW0(w)
        t[Lambdacheck[1] - 2*Lambdacheck[2]] * s1

    Note that for untwisted affine type the dual form of the classical Weyl group
    is isomorphic to the usual one, but acts on a different lattice and is therefore different to sage::

        sage: W0v = E.dual_classical_weyl(); W0v
        Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
        sage: v = W0v.from_reduced_word([1,2])
        sage: x = PvW0(v); x
        s1*s2
        sage: y = PW0(v); y
        s1*s2
        sage: x == y
        False
        sage: x.parent() == y.parent()
        False

    An element can be created directly from a reduced word.

    ::

        sage: PW0.from_reduced_word([2,1,0])
        t[Lambdacheck[1] - 2*Lambdacheck[2]] * s1

    Here is a demonstration of the printing options.

    ::

        sage: PW0 = ExtendedAffineWeylGroup(["A",2,1], affine="sx", classical="Sx",translation="x",fundamental="pix")
        sage: E = PW0.realization_of()
        sage: y = PW0(E.lattice_basis()[1])
        sage: y
        x[Lambdacheck[1]]
        sage: FW = PW0.realization_of().FW()
        sage: FW(y)
        pix[1] * sx2*sx1
        sage: PW0.an_element()
        x[2*Lambdacheck[1] + 2*Lambdacheck[2]] * Sx1*Sx2

    .. TODO::

        - Implement a "slow" action of `E` on any affine root or weight lattice realization.
        - Implement the level `m` actions of `E` and `W` on the lattices of finite type.
        - Implement the relevant methods from the usual affine weyl group
        - Implementation by matrices: style "M".
        - Use case: implement the Hecke algebra on top of this

    The semidirect product construction in sage currently only
    admits multiplicative groups. Therefore for the styles involving "P" and "Pv", one must
    convert the additive group of translations `L` into a multiplicative group by
    applying the :class:`sage.groups.group_exp.GroupExp` functor.

    """

    cartan_type = CartanType(cartan_type)
    if cartan_type.is_reducible():
        raise ValueError, "Extended affine Weyl groups are only implemented for irreducible affine Cartan types"
    if cartan_type.is_finite(): # a finite Cartan type is an abbreviation for its untwisted affinization
        cartan_type = cartan_type.affine()
    elif not cartan_type.is_affine():
        raise ValueError, "Cartan type must be finite or affine"

    if style == "PW0":
        return ExtendedAffineWeylGroup_Class(cartan_type, **print_options).PW0()
    elif style == "W0P":
        return ExtendedAffineWeylGroup_Class(cartan_type, **print_options).W0P()
    elif style == "WF":
        return ExtendedAffineWeylGroup_Class(cartan_type, **print_options).WF()
    elif style == "FW":
        return ExtendedAffineWeylGroup_Class(cartan_type, **print_options).FW()
    elif style == "PvW0":
        return ExtendedAffineWeylGroup_Class(cartan_type, **print_options).PvW0()
    elif style == "W0Pv":
        return ExtendedAffineWeylGroup_Class(cartan_type, **print_options).W0Pv()
    else:
        raise ValueError, "Extended affine Weyl group of style %s is not implemented"%style
    return None

class ExtendedAffineWeylGroup_Class(UniqueRepresentation, Parent):
    r"""
    The parent-with-realization class of an extended affine Weyl group.
    """

    def __init__(self, cartan_type, **print_options):
        if not cartan_type.is_affine():
            raise ValueError, "%s is not affine"%cartan_type

        self._cartan_type = cartan_type

        self._prefixt = "t"
        self._prefixf = "pi"
        self._prefixcl = None
        self._prefixaf = None
        self._print_tuple = False

        for option in print_options.keys():
            if option == 'translation':
                self._prefixt = print_options['translation']
            elif option == 'fundamental':
                self._prefixf = print_options['fundamental']
            elif option == 'print_tuple':
                self._print_tuple = print_options['print_tuple']
            elif option == 'affine':
                self._prefixaf = print_options['affine']
            elif option == 'classical':
                self._prefixcl = print_options['classical']
            else:
                raise ValueError, "Print option %s is unrecognized"%option

        if self._prefixaf:
            if not self._prefixcl:
                if self._prefixaf.islower():
                    self._prefixcl = self._prefixaf.upper()
                else:
                    self._prefixcl = self._prefixaf.lower()
        elif self._prefixcl:
            if self._prefixcl.islower():
                self._prefixaf = self._prefixcl.upper()
            else:
                self._prefixaf = self._prefixcl.lower()
        else:
            self._prefixaf = "S"
            self._prefixcl = "s"

        self._ct0 = cartan_type.classical()
        self._R0  = self._ct0.root_system()
        self._I0 = self._ct0.index_set()
        self._ct0v = self._ct0.dual()
        self._R0v = self._ct0v.root_system()
        self._a0check = self._cartan_type.acheck()[0]

        # `BC` (`A_{2n}^{(2)\dagger}`) is considered untwisted and its dual is considered twisted
        self._untwisted = (self._cartan_type.is_untwisted_affine() or self._cartan_type.dual().type() == 'BC')

        # fundamental group
        self._fundamental_group = FundamentalGroupOfExtendedAffineWeylGroup(cartan_type, prefix=self._prefixf)

        # lattice data
        if self._untwisted:
            self._lattice = self._R0.coweight_lattice()
            self._basis = self._lattice.fundamental_weights()
            self._basis_name = 'Lambdacheck'
            self._simpleR0 = self._R0.root_lattice().simple_roots()
            if self._cartan_type.dual().type() == 'BC':
                # A_{2n}^{(2)} dual
                self._special_root = self._R0.coroot_lattice().highest_root()
                self._special_translation = self._lattice.fundamental_weight(1)
            else:
                # untwisted affine case
                self._special_root = self._R0.root_lattice().highest_root().associated_coroot()
                self._special_translation = self._special_root
            self._special_translation_covector = self._special_root.associated_coroot()
            # in the "Pv" realization, the weight lattice of dual type is used for translations
            self._dual_lattice = self._R0v.weight_lattice()
            self._dual_basis = self._dual_lattice.basis()
            self._dual_basis_name = 'Lambda'
        else:
            self._lattice = self._R0.weight_lattice()
            self._basis = self._lattice.fundamental_weights()
            self._basis_name = 'Lambda'
            self._simpleR0 = self._R0.coroot_lattice().simple_roots()
            if self._cartan_type.type() == 'BC':
                # A_{2n}^{(2)}
                self._special_root = self._R0.root_lattice().highest_root()
                self._special_translation = self._lattice.fundamental_weight(1)
                self._special_translation_covector = 2*self._special_root.associated_coroot()
            else:
                # dual untwisted case
                self._special_root = self._R0.coroot_lattice().highest_root().associated_coroot()
                self._special_translation = self._special_root
                self._special_translation_covector = self._special_root.associated_coroot()

            self._dual_lattice = self._lattice
            self._dual_basis = self._basis
            self._dual_basis_name = 'Lambda'

        # classical and affine Weyl groups
        self._W0 = WeylGroup(self._lattice, prefix=self._prefixcl)
        self._W = WeylGroup(self._cartan_type.root_system().root_lattice(), prefix=self._prefixaf)
        self._special_reflection = self._W0.from_reduced_word(self._special_root.associated_reflection())

        # "Pv" version of classical Weyl group; use same prefix as for W0
        self._W0v = WeylGroup(self._dual_lattice, prefix=self._prefixcl)

        # wrap the lattice into a multiplicative group for internal use in the semidirect product
        self._exp_lattice = GroupExp()(self._lattice)
        self._exp_dual_lattice = GroupExp()(self._dual_lattice)

        self._extended = True

        Parent.__init__(self, category = Groups().WithRealizations())

        # create the realizations (they are cached)
        PW0 = self.PW0()
        W0P = self.W0P()
        WF = self.WF()
        FW = self.FW()
        PvW0 = self.PvW0()
        W0Pv = self.W0Pv()

        # coercions between realizations

        W0P_to_PW0 = SetMorphism(Hom(W0P, PW0, Groups()),lambda x: PW0(x.to_opposite()))
        W0P_to_PW0.register_as_coercion()

        PW0_to_W0P = SetMorphism(Hom(PW0, W0P, Groups()),lambda x: W0P(x.to_opposite()))
        PW0_to_W0P.register_as_coercion()

        FW_to_WF = SetMorphism(Hom(FW, WF, Groups()),lambda x: WF(x.to_opposite()))
        FW_to_WF.register_as_coercion()

        WF_to_FW = SetMorphism(Hom(WF, FW, Groups()),lambda x: FW(x.to_opposite()))
        WF_to_FW.register_as_coercion()

        def PW0_to_WF_func(x):
            r"""
            Coercion from style "PW0" to "WF".
            """
            E = x.parent().realization_of()
            assert x.parent() == E.PW0()
            WF = E.WF()
            i = x.first_descent(side='left')
            if i is None:
                t = x.to_translation_left()
                # t must be zero or a special fundamental basis element
                if t == E.lattice().zero():
                    ispecial = 0
                else:
                    supp = t.support()
                    assert len(supp) == 1
                    ispecial = supp[0]
                return WF.fundamental_group_morphism(E.fundamental_group()(ispecial))
            sx = x.apply_simple_reflection(i, side='left')
            swf = PW0_to_WF_func(sx)
            return swf.apply_simple_reflection(i, side='left')

        PW0_to_WF = SetMorphism(Hom(PW0, WF, Groups()),PW0_to_WF_func)
        PW0_to_WF.register_as_coercion()

        def WF_to_PW0_func(x):
            r"""
            Coercion from style "WF" to "PW0".
            """
            E = x.parent().realization_of()
            assert x.parent() == E.WF(), "type should be WF but it is %s"%type(x.parent())  #!!! delete?
            PW0 = E.PW0()
            w = x.summand_projection(0)
            f = x.summand_projection(1)
            i = w.first_descent(side='left')
            if i is None:
                # the element is in the fundamental group; use its known formula
                ispecial = f.value()
                W=E.classical_weyl()
                if ispecial == 0:
                    weight = E.lattice().zero()
                    wo = W.one()
                else:
                    weight = E.lattice_basis()[ispecial]
                    wo = W.from_reduced_word(E.fundamental_group().finite_action()[ispecial])
                return PW0((weight,wo))
            sx = x.apply_simple_reflection(i, side='left')
            spw0 = WF_to_PW0_func(sx)
            return spw0.apply_simple_reflection(i, side='left')

        WF_to_PW0 = SetMorphism(Hom(WF, PW0, Groups()),WF_to_PW0_func)
        WF_to_PW0.register_as_coercion()

        PvW0_to_W0Pv = SetMorphism(Hom(PvW0, W0Pv, Groups()),lambda x: W0Pv(x.to_opposite()))
        PvW0_to_W0Pv.register_as_coercion()
        W0Pv_to_PvW0 = SetMorphism(Hom(W0Pv, PvW0, Groups()),lambda x: PvW0(x.to_opposite()))
        W0Pv_to_PvW0.register_as_coercion()

        if self._untwisted:
            PW0_to_PvW0 = SetMorphism(Hom(PW0, PvW0, Groups()), lambda x: PvW0((self.exp_dual_lattice()(x.summand_projection(0).value.to_dual_type_cospace()),self.dual_classical_weyl().from_reduced_word(x.summand_projection(1).reduced_word()))))
            PvW0_to_PW0 = SetMorphism(Hom(PvW0, PW0, Groups()), lambda x: PW0((self.exp_lattice()(x.summand_projection(0).value.to_dual_type_cospace()),self.classical_weyl().from_reduced_word(x.summand_projection(1).reduced_word()))))
            W0P_to_W0Pv = SetMorphism(Hom(W0P, W0Pv, Groups()), lambda x: W0Pv((self.dual_classical_weyl().from_reduced_word(x.summand_projection(0).reduced_word()),self.exp_dual_lattice()(x.summand_projection(1).value.to_dual_type_cospace()))))
            W0Pv_to_W0P = SetMorphism(Hom(W0Pv, W0P, Groups()), lambda x: W0P((self.classical_weyl().from_reduced_word(x.summand_projection(0).reduced_word()),self.exp_lattice()(x.summand_projection(1).value.to_dual_type_cospace()))))
        else:
            PW0_to_PvW0 = SetMorphism(Hom(PW0, PvW0, Groups()), lambda x: PvW0((x.summand_projection(0),self.dual_classical_weyl().from_reduced_word(x.summand_projection(1).reduced_word()))))
            PvW0_to_PW0 = SetMorphism(Hom(PvW0, PW0, Groups()), lambda x: PW0((x.summand_projection(0),self.classical_weyl().from_reduced_word(x.summand_projection(1).reduced_word()))))
            W0P_to_W0Pv = SetMorphism(Hom(W0P, W0Pv, Groups()), lambda x: W0Pv((self.dual_classical_weyl().from_reduced_word(x.summand_projection(0).reduced_word()),x.summand_projection(1))))
            W0Pv_to_W0P = SetMorphism(Hom(W0Pv, W0P, Groups()), lambda x: W0P((self.classical_weyl().from_reduced_word(x.summand_projection(0).reduced_word()),x.summand_projection(1))))

        PW0_to_PvW0.register_as_coercion()
        PvW0_to_PW0.register_as_coercion()
        W0P_to_W0Pv.register_as_coercion()
        W0Pv_to_W0P.register_as_coercion()

        # coercions of the translation lattice into the appropriate realizations
        P_to_PW0 = SetMorphism(Hom(self.lattice(), PW0, Sets()), lambda x: PW0.translation_group_morphism(x))
        P_to_PW0.register_as_coercion()
        P_to_W0P = SetMorphism(Hom(self.lattice(), W0P, Sets()), lambda x: W0P.translation_group_morphism(x))
        P_to_W0P.register_as_coercion()
        Pv_to_PvW0 = SetMorphism(Hom(self.dual_lattice(), PvW0, Sets()), lambda x: PvW0.dual_translation_group_morphism(x))
        Pv_to_PvW0.register_as_coercion()
        Pv_to_W0Pv = SetMorphism(Hom(self.dual_lattice(), W0Pv, Sets()), lambda x: W0Pv.dual_translation_group_morphism(x))
        Pv_to_W0Pv.register_as_coercion()

        # coercions of the classical Weyl group into the appropriate realizations

        W0_to_PW0 = SetMorphism(Hom(self.classical_weyl(), PW0, Groups()), lambda x: PW0.classical_weyl_morphism(x))
        W0_to_PW0.register_as_coercion()
        W0_to_W0P = SetMorphism(Hom(self.classical_weyl(), W0P, Groups()), lambda x: W0P.classical_weyl_morphism(x))
        W0_to_W0P.register_as_coercion()
        W0v_to_PvW0 = SetMorphism(Hom(self.dual_classical_weyl(), PvW0, Groups()), lambda x: PvW0.dual_classical_weyl_morphism(x))
        W0v_to_PvW0.register_as_coercion()
        W0v_to_W0Pv = SetMorphism(Hom(self.dual_classical_weyl(), W0Pv, Groups()), lambda x: W0Pv.dual_classical_weyl_morphism(x))
        W0v_to_W0Pv.register_as_coercion()

        # coercions of the fundamental group into the appropriate realizations

        F_to_WF = SetMorphism(Hom(self.fundamental_group(), WF, Groups()), lambda x: WF.fundamental_group_morphism(x))
        F_to_WF.register_as_coercion()
        F_to_FW = SetMorphism(Hom(self.fundamental_group(), FW, Groups()), lambda x: FW.fundamental_group_morphism(x))
        F_to_FW.register_as_coercion()

        # coercions of the affine Weyl group into the appropriate realizations

        W_to_WF = SetMorphism(Hom(self.affine_weyl(), WF, Groups()), lambda x: WF.affine_weyl_morphism(x))
        W_to_WF.register_as_coercion()

        W_to_FW = SetMorphism(Hom(self.affine_weyl(), FW, Groups()), lambda x: FW.affine_weyl_morphism(x))
        W_to_FW.register_as_coercion()

    @cached_method
    def PW0(self):
        r"""
        Realizes ``self`` in "PW0"-style.
        """
        return self.ExtendedAffineWeylGroupPW0()

    @cached_method
    def W0P(self):
        r"""
        Realizes ``self`` in "W0P"-style.
        """
        return self.ExtendedAffineWeylGroupW0P()

    @cached_method
    def WF(self):
        r"""
        Realizes ``self`` in "WF"-style.
        """
        return self.ExtendedAffineWeylGroupWF()

    @cached_method
    def FW(self):
        r"""
        Realizes ``self`` in "FW"-style.
        """
        return self.ExtendedAffineWeylGroupFW()

    @cached_method
    def PvW0(self):
        r"""
        Realizes ``self`` in "PvW0"-style.
        """
        return self.ExtendedAffineWeylGroupPvW0()

    @cached_method
    def W0Pv(self):
        r"""
        Realizes ``self`` in "W0Pv"-style.
        """
        return self.ExtendedAffineWeylGroupW0Pv()

    def cartan_type(self):
        return self._cartan_type

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
            sage: E = PW0.realization_of(); E
            Extended affine Weyl group of type ['A', 2, 1]

        """
        return "Extended affine Weyl group of type %s"%self.cartan_type()

    @cached_method
    def index_set(self):
        r"""
        Returns the index set of the affine Dynkin diagram.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup("B2").realization_of().index_set()
            (0, 1, 2)

        """
        return self.cartan_type().index_set()

    def fundamental_group(self):
        r"""
        Returns the abstract fundamental group.

        EXAMPLES::

            sage: F = ExtendedAffineWeylGroup(['D',5,1]).realization_of().fundamental_group()
            sage: f = F.family(); f
            Finite family {0: pi[0], 1: pi[1], 4: pi[4], 5: pi[5]}
            sage: f[0]
            pi[0]
            sage: f[1]^2
            pi[0]
            sage: F(1)^2
            pi[0]
            sage: f[4]*f[4]
            pi[1]

        """
        return self._fundamental_group

    @cached_method
    def special_nodes(self):
        r"""
        Returns the set of special nodes.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['D',4,2]).realization_of().special_nodes()
            [0, 3]

        """
        return self.fundamental_group().special_nodes()

    def is_finite(self):
        return False

    def lattice(self):
        r"""
        Returns the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().lattice()
            Coweight lattice of the Root system of type ['A', 2]
            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().lattice()
            Weight lattice of the Root system of type ['C', 3]
            sage: ExtendedAffineWeylGroup(['A',4,2]).realization_of().lattice()
            Weight lattice of the Root system of type ['C', 2]
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).realization_of().lattice()
            Coweight lattice of the Root system of type ['B', 2]

        """

        return self._lattice

    def exp_lattice(self):
        r"""
        Returns the multiplicative version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().exp_lattice()
            Multiplicative form of Coweight lattice of the Root system of type ['A', 2]

        """

        return self._exp_lattice

    def lattice_basis(self):
        r"""
        Returns the distinguished basis of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().lattice_basis()
            Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}
            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2], 3: Lambda[3]}
            sage: ExtendedAffineWeylGroup(['A',4,2]).realization_of().lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2]}
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).realization_of().lattice_basis()
            Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}

        """
        return self._basis

    def dual_lattice(self):
        r"""
        Returns the dual version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().dual_lattice()
            Weight lattice of the Root system of type ['A', 2]
            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().dual_lattice()
            Weight lattice of the Root system of type ['C', 3]

        """

        return self._dual_lattice

    def exp_dual_lattice(self):
        r"""
        Returns the multiplicative version of the dual version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().exp_dual_lattice()
            Multiplicative form of Weight lattice of the Root system of type ['A', 2]

        """

        return self._exp_dual_lattice

    def dual_lattice_basis(self):
        r"""
        Returns the distinguished basis of the dual version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().dual_lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2]}
            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().dual_lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2], 3: Lambda[3]}

        """
        return self._dual_basis

    def classical_weyl(self):
        r"""
        Returns the classical Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the weight lattice)
            sage: ExtendedAffineWeylGroup(['A',4,2]).realization_of().classical_weyl()
            Weyl Group of type ['C', 2] (as a matrix group acting on the weight lattice)
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).realization_of().classical_weyl()
            Weyl Group of type ['C', 2] (as a matrix group acting on the coweight lattice)

        """

        return self._W0

    def dual_classical_weyl(self):
        r"""
        Returns the dual version of the classical Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().dual_classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().dual_classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the weight lattice)

        """

        return self._W0v

    def affine_weyl(self):
        r"""
        Returns the affine Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).realization_of().affine_weyl()
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().affine_weyl()
            Weyl Group of type ['B', 3, 1]^* (as a matrix group acting on the root lattice)
            sage: ExtendedAffineWeylGroup(['A',4,2]).realization_of().affine_weyl()
            Weyl Group of type ['BC', 2, 2] (as a matrix group acting on the root lattice)
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).realization_of().affine_weyl()
            Weyl Group of type ['BC', 2, 2]^* (as a matrix group acting on the root lattice)

        """
        return self._W

    def classical_weyl_to_affine_morphism(self, w):
        r"""
        The image of `w` under the homomorphism from the classical Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: E = ExtendedAffineWeylGroup(['A',2,1]).realization_of()
            sage: W0 = E.classical_weyl()
            sage: w = W0.from_reduced_word([1,2]); w
            s1*s2
            sage: v = E.classical_weyl_to_affine_morphism(w); v
            S1*S2

        """
        return self.affine_weyl().from_reduced_word(w.reduced_word())

    def dual_classical_weyl_to_affine_morphism(self, w):
        r"""
        The image of `w` under the homomorphism from the dual version of the classical Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: E = ExtendedAffineWeylGroup(['A',2,1]).realization_of()
            sage: W0v = E.dual_classical_weyl()
            sage: w = W0v.from_reduced_word([1,2]); w
            s1*s2
            sage: v = E.dual_classical_weyl_to_affine_morphism(w); v
            S1*S2

        """
        return self.affine_weyl().from_reduced_word(w.reduced_word())

    def classical_root_system(self):
        r"""
        Returns the root system of the classical subrootsystem of ``self`` obtained by removing the affine `0` node.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',5,2]).realization_of().classical_root_system()
            Root system of type ['C', 3]

        """
        return self._R0

    def cardinality(self):
        """
        Returns infinity
        """
        return Infinity

    def a_realization(self):
        r"""
        Returns the default realization of an extended affine Weyl group.

        EXAMPLES::

            sage: PW0 = ExtendedAffineWeylGroup(['A',2,1]); PW0
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Multiplicative form of Coweight lattice of the Root system of type ['A', 2] acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
            sage: PW0 == PW0.realization_of().PW0()
            True

        """
        return self.PW0()

    class Realizations(Category_realization_of_parent):
        r"""
        The category of the realizations of an extended affine Weyl group
        """

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: R = ExtendedAffineWeylGroup(['A',2,1]).realization_of().Realizations(); R
                Category of realizations of Extended affine Weyl group of type ['A', 2, 1]
                sage: R.super_categories()
                [Category of associative inverse unital realizations of magmas]

            """
            return [Groups().Realizations()]

        class ParentMethods:

            @cached_method
            def one(self):
                r"""
                Returns the unit element.

                This default implementation takes the unit in the
                PW0 realization and coerces it into `self`.

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',2,1]).one()
                    1

                .. warning::

                    Must be implemented in style "PW0".

                """
                PW0 = self.realization_of().PW0()
                return self(PW0.one())

            @cached_method
            def fundamental_group_morphism(self, x):
                r"""
                Returns the image of `x` under the homomorphism from the fundamental group into
                ``self``.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: E = PW0.realization_of()
                    sage: Is = E.special_nodes()
                    sage: F = E.fundamental_group()
                    sage: [(i, PW0.fundamental_group_morphism(F(i))) for i in Is]
                    [(0, 1), (1, t[Lambdacheck[1]] * s1*s2*s3), (2, t[Lambdacheck[2]] * s2*s3*s1*s2), (3, t[Lambdacheck[3]] * s3*s2*s1)]
                    sage: [(i, E.W0P().fundamental_group_morphism((F(i)))) for i in Is]
                    [(0, 1), (1, s1*s2*s3 * t[-Lambdacheck[3]]), (2, s2*s3*s1*s2 * t[-Lambdacheck[2]]), (3, s3*s2*s1 * t[-Lambdacheck[1]])]
                    sage: [(i, E.WF().fundamental_group_morphism(F(i))) for i in Is]
                    [(0, 1), (1, pi[1]), (2, pi[2]), (3, pi[3])]

                .. warning::

                    This method must be implemented by the "WF" and "FW" realizations.

                """
                WF = self.realization_of().WF()
                return self(WF.fundamental_group_morphism(x))

            def translation_group_morphism(self, la):
                r"""
                Returns the image of `la` under the homomorphism of the translation lattice into
                ``self``.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
                    sage: E = PW0.realization_of()
                    sage: b = E.lattice_basis(); b
                    Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}
                    sage: x = PW0.translation_group_morphism(2*b[1]-b[2]); x
                    t[2*Lambdacheck[1] - Lambdacheck[2]]
                    sage: FW = E.FW()
                    sage: y = FW.translation_group_morphism(2*b[1]-b[2]); y
                    S0*S2*S0*S1
                    sage: FW(x) == y
                    True

                Since the implementation as a semidirect product requires
                wrapping the lattice group to make it multiplicative,
                we cannot declare that this map is a morphism for
                sage ``Groups()``.

                .. warning::

                    This method must be implemented by the "PW0" and "W0P" realizations.
                """
                PW0 = self.realization_of().PW0()
                return self(PW0.translation_group_morphism(la))

            def dual_translation_group_morphism(self, la):
                r"""
                Returns the image of `la` under the homomorphism of the dual version of the translation lattice into
                ``self``.

                EXAMPLES::

                    sage: PvW0 = ExtendedAffineWeylGroup(['A',2,1],style="PvW0")
                    sage: E = PvW0.realization_of()
                    sage: bv = E.dual_lattice_basis(); bv
                    Finite family {1: Lambda[1], 2: Lambda[2]}
                    sage: x = PvW0.dual_translation_group_morphism(2*bv[1]-bv[2]); x
                    t[2*Lambda[1] - Lambda[2]]
                    sage: FW = E.FW()
                    sage: y = FW.dual_translation_group_morphism(2*bv[1]-bv[2]); y
                    S0*S2*S0*S1
                    sage: FW(x) == y
                    True

                """

                return self(self.realization_of().PvW0().dual_translation_group_morphism(la))

            @abstract_method
            def simple_reflections(self):
                r"""
                Returns a family from the set of affine Dynkin nodes to the simple reflections
                in the realization of the extended affine Weyl group.

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1],style="W0P").simple_reflections()
                    Finite family {0: s1*s2*s3*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[3]], 1: s1, 2: s2, 3: s3}
                    sage: ExtendedAffineWeylGroup(['A',3,1],style="WF").simple_reflections()
                    Finite family {0: S0, 1: S1, 2: S2, 3: S3}
                    sage: ExtendedAffineWeylGroup(['A',3,1],style="FW",print_tuple=True).simple_reflections()
                    Finite family {0: (pi[0], S0), 1: (pi[0], S1), 2: (pi[0], S2), 3: (pi[0], S3)}
                    sage: ExtendedAffineWeylGroup(['A',3,1],style="FW",fundamental="f",print_tuple=True).simple_reflections()
                    Finite family {0: (f[0], S0), 1: (f[0], S1), 2: (f[0], S2), 3: (f[0], S3)}
                    sage: ExtendedAffineWeylGroup(['A',3,1],style="PvW0").simple_reflections()
                    Finite family {0: t[Lambda[1] + Lambda[3]] * s1*s2*s3*s2*s1, 1: s1, 2: s2, 3: s3}

                """

            def simple_reflection(self, i):
                r"""
                Returns the `i`-th simple reflection in the realization of the extended affine
                Weyl group.

                    sage: ExtendedAffineWeylGroup(['A',3,1],style="PW0").simple_reflection(0)
                    t[Lambdacheck[1] + Lambdacheck[3]] * s1*s2*s3*s2*s1
                    sage: ExtendedAffineWeylGroup(['A',3,1],style="WF").simple_reflection(0)
                    S0

                """
                return self.simple_reflections()[i]

            def classical_weyl_morphism(self, w):
                r"""
                Returns the image of `w` from the finite Weyl group into ``self``.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: E = PW0.realization_of()
                    sage: W0 = E.classical_weyl()
                    sage: w = W0.from_reduced_word([2,1,3])
                    sage: y = PW0.classical_weyl_morphism(w); y
                    s2*s3*s1
                    sage: y.parent() == PW0
                    True
                    sage: y.to_classical_weyl() == w
                    True
                    sage: W0P = E.W0P()
                    sage: z = W0P.classical_weyl_morphism(w); z
                    s2*s3*s1
                    sage: z.parent() == W0P
                    True
                    sage: W0P(y) == z
                    True
                    sage: FW = E.FW()
                    sage: x = FW.classical_weyl_morphism(w); x
                    S2*S3*S1
                    sage: x.parent() == FW
                    True
                    sage: FW(y) == x
                    True
                    sage: FW(z) == x
                    True

                .. warning::

                    Must be implemented in style "PW0" and "W0P".

                """
                PW0 = self.realization_of().PW0()
                return self(PW0.classical_weyl_morphism(w))

            def dual_classical_weyl_morphism(self, w):
                r"""
                Returns the image of `w` from the finite Weyl group of dual form into ``self``.

                EXAMPLES::

                    sage: PvW0 = ExtendedAffineWeylGroup(['A',3,1],style="PvW0")
                    sage: E = PvW0.realization_of()
                    sage: W0v = E.dual_classical_weyl()
                    sage: w = W0v.from_reduced_word([2,1,3])
                    sage: y = PvW0.dual_classical_weyl_morphism(w); y
                    s2*s3*s1
                    sage: y.parent() == PvW0
                    True
                    sage: y.to_dual_classical_weyl() == w
                    True
                    sage: x = E.FW().dual_classical_weyl_morphism(w); x
                    S2*S3*S1
                    sage: PvW0(x) == y
                    True

                .. warning::

                    Must be implemented in style "PvW0" and "W0Pv".

                """
                return self(self.realization_of().PvW0().dual_classical_weyl_morphism(w))

            def affine_weyl_morphism(self, w):
                r"""
                Returns the image of `w` under the homomorphism from the affine Weyl group
                into ``self``.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: E = PW0.realization_of()
                    sage: W = E.affine_weyl()
                    sage: w = W.from_reduced_word([2,1,3,0])
                    sage: x = PW0.affine_weyl_morphism(w); x
                    t[Lambdacheck[1] - 2*Lambdacheck[2] + Lambdacheck[3]] * s3*s1
                    sage: FW = E.FW()
                    sage: y = FW.affine_weyl_morphism(w); y
                    S2*S3*S1*S0
                    sage: FW(x) == y
                    True

                .. warning::

                    Must be implemented in style "WF" and "FW".

                """
                WF = self.realization_of().WF()
                return self(WF.affine_weyl_morphism(w))

            def from_reduced_word(self, word):
                r"""
                Converts an affine or finite reduced word into a group element.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
                    sage: PW0.from_reduced_word([1,0,1,2])
                    t[-Lambdacheck[1] + 2*Lambdacheck[2]]

                """

                return self.affine_weyl_morphism(self.realization_of().affine_weyl().from_reduced_word(word))

        class ElementMethods:

            @abstract_method
            def has_descent(self, i, side='right', positive=False):
                r"""
                Returns whether ``self`` * `s_i` < ``self`` where `s_i` is the `i`-th simple
                reflection in the realized group.

                INPUT:

                    - `i` -- An affine Dynkin node
                    - `side` -- 'right' or 'left' (default: 'right')
                    - `positive` -- True or False (default: False)

                If ``side``='left' then the reflection acts
                on the left. If ``positive`` = True then the inequality is reversed.

                EXAMPLES::

                    sage: WF = ExtendedAffineWeylGroup(['A',3,1],style="WF")
                    sage: x = WF.an_element(); x
                    S0*S1*S2*S3 * pi[2]
                    sage: I = WF.realization_of().index_set()
                    sage: [(i, x.has_descent(i)) for i in I]
                    [(0, False), (1, True), (2, False), (3, False)]
                    sage: [(i, x.has_descent(i,side='left')) for i in I]
                    [(0, True), (1, False), (2, False), (3, False)]
                    sage: [(i, x.has_descent(i,positive=True)) for i in I]
                    [(0, True), (1, False), (2, True), (3, True)]

                .. warning::

                    This method is abstract because it is used in the recursive coercions
                    between "PW0" and "WF" and other methods use this coercion.

                """

            def first_descent(self, side='right', positive=False, index_set=None):
                r"""
                Returns the first descent of ``self``.

                INPUT:

                    - ``side`` -- 'left' or 'right' (default: 'right')
                    - ``positive`` -- True or False (default: False)
                    - ``index_set`` -- an optional subset of Dynkin nodes

                If ``index_set`` is not None, then the descent must be in the ``index_set``.

                EXAMPLES::

                    sage: WF = ExtendedAffineWeylGroup(['A',3,1],style="WF")
                    sage: x = WF.an_element(); x
                    S0*S1*S2*S3 * pi[2]
                    sage: x.first_descent()
                    1
                    sage: x.first_descent(side='left')
                    0
                    sage: x.first_descent(positive=True)
                    0
                    sage: x.first_descent(side='left',positive=True)
                    1

                """
                if index_set is None:
                    index_set = self.parent().realization_of().index_set()
                for i in index_set:
                    if self.has_descent(i, side=side, positive=positive):
                        return i
                return None

            def apply_simple_reflection(self, i, side='right'):
                r"""
                Apply the `i`-th simple reflection to ``self``.

                EXAMPLES::

                    sage: WF = ExtendedAffineWeylGroup(['A',3,1],style="WF")
                    sage: x = WF.an_element(); x
                    S0*S1*S2*S3 * pi[2]
                    sage: x.apply_simple_reflection(1)
                    S0*S1*S2 * pi[2]
                    sage: x.apply_simple_reflection(1, side='left')
                    S0*S1*S2*S0*S3 * pi[2]

                """
                s = self.parent().simple_reflection(i)
                if side == 'right':
                    return self*s
                else:
                    return s*self

            def apply_simple_projection(self, i, side='right', length_increasing=True):
                r"""
                Returns ``self`` * `s_i` if it is greater than ``self`` and otherwise
                returns ``self``.

                If ``side`` = 'left' then `s_i` is multiplied on the left.
                If ``positive`` = True then replace the word "greater" by "less".

                EXAMPLES::

                    sage: WF = ExtendedAffineWeylGroup(['A',3,1],style="WF")
                    sage: x = WF.an_element(); x
                    S0*S1*S2*S3 * pi[2]
                    sage: x.apply_simple_projection(1)
                    S0*S1*S2*S3 * pi[2]
                    sage: x.apply_simple_projection(1, length_increasing=False)
                    S0*S1*S2 * pi[2]

                """
                if self.has_descent(i, side=side, positive=length_increasing):
                    return self.apply_simple_reflection(i, side=side)
                return self

            def to_fundamental_group(self):
                r"""
                Returns the image of ``self`` under the homomorphism to the fundamental group.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: b = PW0.realization_of().lattice_basis()
                    sage: [(x, PW0.translation_group_morphism(x).to_fundamental_group()) for x in b]
                    [(Lambdacheck[1], pi[1]), (Lambdacheck[2], pi[2]), (Lambdacheck[3], pi[3])]

                .. warning::

                    Must be implemented in style "WF".

                """
                WF = self.parent().realization_of().WF()
                return WF(self).to_fundamental_group()

            def to_classical_weyl(self):
                r"""
                Returns the image of ``self`` under the homomorphism to the classical Weyl group.

                    sage: WF = ExtendedAffineWeylGroup(['A',3,1],style="WF")
                    sage: WF.simple_reflection(0).to_classical_weyl()
                    s1*s2*s3*s2*s1

                .. warning::

                    Must be implemented in style "PW0".

                """
                PW0 = self.parent().realization_of().PW0()
                return PW0(self).to_classical_weyl()

            def to_dual_classical_weyl(self):
                r"""
                Returns the image of ``self`` under the homomorphism to the dual form of the classical Weyl group.

                    sage: WF = ExtendedAffineWeylGroup(['A',3,1],style="WF")
                    sage: x = WF.simple_reflection(0).to_dual_classical_weyl(); x
                    s1*s2*s3*s2*s1
                    sage: x.parent()
                    Weyl Group of type ['A', 3] (as a matrix group acting on the weight lattice)

                .. warning::

                    Must be implemented in style "PvW0".

                """
                PvW0 = self.parent().realization_of().PvW0()
                return PvW0(self).to_dual_classical_weyl()

            def to_affine_weyl_left(self):
                r"""
                Returns the projection of ``self`` to the affine Weyl group on the left,
                after factorizing using the style "WF".

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: b = PW0.realization_of().lattice_basis()
                    sage: [(x,PW0.translation_group_morphism(x).to_affine_weyl_left()) for x in b]
                    [(Lambdacheck[1], S0*S3*S2), (Lambdacheck[2], S0*S3*S1*S0), (Lambdacheck[3], S0*S1*S2)]

                .. warning::
                    Must be implemented in style "WF".

                """
                WF = self.parent().realization_of().WF()
                return WF(self).to_affine_weyl_left()

            def to_affine_weyl_right(self):
                r"""
                Returns the projection of ``self`` to the affine Weyl group on the right,
                after factorizing using the style "FW".

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: b = PW0.realization_of().lattice_basis()
                    sage: [(x,PW0.translation_group_morphism(x).to_affine_weyl_right()) for x in b]
                    [(Lambdacheck[1], S3*S2*S1), (Lambdacheck[2], S2*S3*S1*S2), (Lambdacheck[3], S1*S2*S3)]

                .. warning::

                    Must be implemented in style "FW".

                """
                FW = self.parent().realization_of().FW()
                return FW(self).to_affine_weyl_right()

            def to_translation_left(self):
                r"""
                Returns the projection of ``self`` to the translation lattice after factorizing
                it to the left using the style "PW0".

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: PW0.simple_reflection(0).to_translation_left()
                    Lambdacheck[1] + Lambdacheck[3]

                .. warning::

                    Must be implemented in style "PW0".

                """
                PW0 = self.parent().realization_of().PW0()
                return PW0(self).to_translation_left()

            def to_translation_right(self):
                r"""
                Returns the projection of ``self`` to the translation lattice after factorizing
                it to the right using the style "W0P".

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: PW0.simple_reflection(0).to_translation_right()
                    -Lambdacheck[1] - Lambdacheck[3]

                .. warning::

                    Must be implemented in style "W0P".

                """
                W0P = self.parent().realization_of().W0P()
                return W0P(self).to_translation_right()

            def to_dual_translation_left(self):
                r"""
                Returns the projection of ``self`` to the dual translation lattice after factorizing
                it to the left using the style "PvW0".

                EXAMPLES::

                    sage: PvW0 = ExtendedAffineWeylGroup(['A',3,1],style="PvW0")
                    sage: PvW0.simple_reflection(0).to_dual_translation_left()
                    Lambda[1] + Lambda[3]

                .. warning::

                    Must be implemented in style "PvW0".

                """
                PvW0 = self.parent().realization_of().PvW0()
                return PvW0(self).to_dual_translation_left()

            def to_dual_translation_right(self):
                r"""
                Returns the projection of ``self`` to the dual translation lattice after factorizing
                it to the right using the style "W0Pv".

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                    sage: PW0.simple_reflection(0).to_dual_translation_right()
                    -Lambda[1] - Lambda[3]

                .. warning::

                    Must be implemented in style "W0Pv".

                """
                W0Pv = self.parent().realization_of().W0Pv()
                return W0Pv(self).to_dual_translation_right()

            def length(self):
                r"""
                Returns the length of ``self`` in the Coxeter group sense.

                EXAMPLES::

                     sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                     sage: E = PW0.realization_of()
                     sage: I0 = E.classical_root_system().index_set()
                     sage: [PW0.translation_group_morphism(E.lattice_basis()[i]).length() for i in I0]
                     [3, 4, 3]

                 """
                return self.to_affine_weyl_left().length()

            def coset_representative(self, index_set, side='right'):
                r"""
                Returns the minimum length representative in the coset of ``self`` with respect to
                the subgroup generated by the reflections given by ``index_set``. If ``side`` is
                'left', the subgroup is on the left.

                EXAMPLES::

                     sage: WF = ExtendedAffineWeylGroup(['A',3,1],style="WF")
                     sage: E = WF.realization_of()
                     sage: b = E.lattice_basis()
                     sage: I0 = E.classical_root_system().index_set()
                     sage: [WF.translation_group_morphism(x).coset_representative(index_set=I0) for x in b]
                     [pi[1], pi[2], pi[3]]

                """

                while True:
                    i = self.first_descent(index_set=index_set, side=side)
                    if i is None:
                        return self
                    self = self.apply_simple_reflection(i,side=side)

            def is_grassmannian(self, index_set, side='right'):
                r"""
                Returns whether ``self`` is of minimum length in its coset with respect to the
                subgroup generated by the reflections of ``index_set``.

                EXAMPLES::

                     sage: PW0 = ExtendedAffineWeylGroup(['A',3,1])
                     sage: E = PW0.realization_of()
                     sage: x = PW0.translation_group_morphism(E.lattice_basis()[1]); x
                     t[Lambdacheck[1]]
                     sage: I = E.cartan_type().index_set()
                     sage: [(i, x.is_grassmannian(index_set=[i])) for i in I]
                     [(0, True), (1, False), (2, True), (3, True)]
                     sage: [(i, x.is_grassmannian(index_set=[i], side='left')) for i in I]
                     [(0, False), (1, True), (2, True), (3, True)]

                """

                return self == self.coset_representative(index_set=index_set,side=side)

            def to_affine_grassmannian(self):
                r"""
                Returns the unique affine Grassmannian element in the same coset of ``self``
                with respect to the finite Weyl group acting on the right.

                EXAMPLES::

                     sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
                     sage: elts = PW0.some_elements()
                     sage: [(x, x.to_affine_grassmannian()) for x in elts]
                     [(t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2, t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2*s1)]

                """

                return self.coset_representative(index_set=self.parent().realization_of().classical_root_system().index_set())

            def is_affine_grassmannian(self):
                r"""
                Returns whether ``self`` is affine Grassmannian.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
                    sage: E = PW0.realization_of()
                    sage: F = E.fundamental_group().family()
                    sage: [(x,PW0.fundamental_group_morphism(x).is_affine_grassmannian()) for x in F]
                    [(pi[0], True), (pi[1], True), (pi[2], True)]
                    sage: b = E.lattice_basis()
                    sage: [(-x,PW0.translation_group_morphism(-x).is_affine_grassmannian()) for x in b]
                    [(-Lambdacheck[1], True), (-Lambdacheck[2], True)]

                """

                return self == self.to_affine_grassmannian()

            def bruhat_le(self, x):
                r"""
                Returns whether ``self`` <= `x` in Bruhat order.

                EXAMPLES::

                    sage: WF = ExtendedAffineWeylGroup(['A',2,1],print_tuple=True,style="WF")
                    sage: E = WF.realization_of()
                    sage: W = E.affine_weyl()
                    sage: v = W.from_reduced_word([2,1,0])
                    sage: w = W.from_reduced_word([2,0,1,0])
                    sage: v.bruhat_le(w)
                    True
                    sage: vx = WF.affine_weyl_morphism(v); vx
                    (S2*S1*S0, pi[0])
                    sage: wx = WF.affine_weyl_morphism(w); wx
                    (S2*S0*S1*S0, pi[0])
                    sage: vx.bruhat_le(wx)
                    True
                    sage: F = E.fundamental_group()
                    sage: f = WF.fundamental_group_morphism(F(2))
                    sage: vx.bruhat_le(wx*f)
                    False
                    sage: (vx*f).bruhat_le(wx*f)
                    True

                .. warning::

                    Must be implemented by "WF".

                """
                WF = self.parent().realization_of().WF()
                return WF(self).bruhat_le(WF(x))

            def is_translation(self):
                r"""
                Returns whether ``self`` is a translation element or not.

                EXAMPLES::

                    sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW")
                    sage: E = FW.realization_of()
                    sage: F = E.fundamental_group()
                    sage: FW.affine_weyl_morphism(E.affine_weyl().from_reduced_word([1,2,1,0])).is_translation()
                    True
                    sage: FW.translation_group_morphism(E.lattice_basis()[1]).is_translation()
                    True
                    sage: FW.simple_reflection(0).is_translation()
                    False

                """
                w = self.to_classical_weyl()
                return w == w.parent().one()

            def action(self, la):
                r"""
                Action of ``self`` on a lattice element ``la``.

                EXAMPLES::

                    sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW",affine="s")
                    sage: E = FW.realization_of()
                    sage: x = FW.an_element(); x
                    pi[2] * s0*s1*s2
                    sage: la = E.lattice().an_element(); la
                    2*Lambdacheck[1] + 2*Lambdacheck[2]
                    sage: x.action(la)
                    5*Lambdacheck[1] - 3*Lambdacheck[2]

                .. warning::

                    Must be implemented by style "PW0".

                """
                PW0 = self.parent().realization_of().PW0()
                return PW0(self).action(la)

            def dual_action(self, la):
                r"""
                Action of ``self`` on a dual lattice element ``la``.

                EXAMPLES::

                    sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW",affine="s")
                    sage: E = FW.realization_of()
                    sage: x = FW.an_element(); x
                    pi[2] * s0*s1*s2
                    sage: la = E.dual_lattice().an_element(); la
                    2*Lambda[1] + 2*Lambda[2]
                    sage: x.dual_action(la)
                    5*Lambda[1] - 3*Lambda[2]

                .. warning::

                    Must be implemented by style "PvW0".

                """
                PvW0 = self.parent().realization_of().PvW0()
                return PvW0(self).dual_action(la)

            def action_on_affine_roots(self, beta):
                r"""
                Act by ``self`` on the affine root lattice element ``beta``.

                .. warning::

                    Must be implemented by style "FW".
                """
                WR = self.parent().realization_of()
                assert beta in RootSystem(WR.cartan_type()).root_lattice()
                return WR.FW()(self).action_on_affine_roots(beta)

            def face_data(self, i):
                r"""
                Returns data for the alcove of ``self`` and the face labeled by the Dynkin node ``i``.

                Let the extended affine Weyl group element ``self`` act on the affine simple root `\alpha_i`:

                ..MATH::

                    ``self`` `\alpha_i` = m `\delta` + `\beta`

                Returns the 2-tuple `(m,\beta)` where `m` is the height of the bounding hyperplane and `\beta` is
                the classical root giving the normal vector.

                EXAMPLES::

                    sage: We=ExtendedAffineWeylGroup(['A',2,1])
                    sage: x = We.an_element(); x
                    t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2
                    sage: x.face_data(0)
                    (-1, alpha[1])

                """
                Qaf = RootSystem(self.parent().realization_of().cartan_type()).root_lattice()
                gamma = self.action_on_affine_roots(Qaf.simple_root(i))
                return gamma[0], Qaf.classical()(gamma)

            def alcove_walk_signs(self):
                r"""
                Returns a signed alcove walk for ``self``.

                The returned data is a 3-tuple `g`, ``rw``, ``signs``.
                Let ``self`` = `g` * `w` where `g` has length zero and `w` has reduced word ``rw``.
                Starting with `g` and applying simple reflections from `rw` on the right, one obtains
                a sequence of extended affine Weyl group elements (that is, alcoves) and simple roots.
                The signs give the sequence of sides on which the alcoves lie, relative to the face
                indicated by the simple roots.

                EXAMPLES::

                    sage: We=ExtendedAffineWeylGroup(['A',3,1],style="FW")
                    sage: w = We.from_reduced_word([0,2,1,3,0])*We.fundamental_group_morphism(1); w
                    pi[1] * S3*S1*S2*S0*S3
                    sage: w.alcove_walk_signs()
                    (pi[1], [3, 1, 2, 0, 3], [-1, 1, -1, -1, 1])

                """
                We = self.parent()
                gw = We.realization_of().FW()(self)
                g = gw.summand_projection(0)
                w = gw.summand_projection(1)
                rw = w.reduced_word()
                u_curr = We.fundamental_group_morphism(g.value())
                signs=[]
                for i in rw:
                    m, beta = u_curr.face_data(i)
                    if beta.is_positive_root():
                        signs = signs + [1]
                    else:
                        signs = signs + [-1]
                    u_curr = u_curr * We.simple_reflection(i)
                return g, rw, signs

    class ExtendedAffineWeylGroupPW0Element(GroupSemidirectProduct.Element):
        r"""
        The element class for the "PW0" realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            Whether `self` has `i` as a descent.

            INPUT:

                - `i` - an index.

            OPTIONAL:

                - side -- 'left' or 'right' (default: 'right')
                - positive -- True or False (default: False)

            EXAMPLES::

                sage: We = ExtendedAffineWeylGroup(['A',4,2])
                sage: w = We.from_reduced_word([0,1]); w
                t[Lambda[1]] * s1*s2
                sage: w.has_descent(0, side='left')
                True

            """

            E = self.parent().realization_of()
            if side == 'right':
                self = ~self
            if positive:
                return not self.has_descent(i, side='left')
            la = self.summand_projection(0).value
            w = self.summand_projection(1)
            if i == 0:
                ip = la.scalar(E._special_translation_covector) * E._a0check
                if ip > 1:
                    return True
                if ip < 1:
                    return False
                return E._special_root.weyl_action(w, inverse=True).is_positive_root()
            ip = la.scalar(E._simpleR0[i]) # test height versus simple (co)root
            if ip < 0:
                return True
            if ip > 0:
                return False
            return w.has_descent(i, side='left')

        def action(self, la):
            r"""
            Returns the action of ``self`` on an element ``la`` of the translation lattice.

            EXAMPLES::

                sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
                sage: E = PW0.realization_of()
                sage: x = PW0.an_element(); x
                t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2
                sage: la = E.lattice().an_element(); la
                2*Lambdacheck[1] + 2*Lambdacheck[2]
                sage: x.action(la)
                -2*Lambdacheck[1] + 4*Lambdacheck[2]

            """
            w = self.summand_projection(1)
            assert la in w.parent().domain()
            return self.summand_projection(0).value + w.action(la)

        def to_translation_left(self):
            r"""
            The image of ``self`` under the map that projects to the translation lattice
            factor after factoring it to the left as in style "PW0".

            EXAMPLES::

                sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
                sage: s = PW0.S0(); s
                t[Lambdacheck[1] + Lambdacheck[2]] * s1*s2*s1
                sage: s.to_translation_left()
                Lambdacheck[1] + Lambdacheck[2]

            """
            return self.summand_projection(0).value # undo the GroupExp

        def to_classical_weyl(self):
            r"""
            Returns the image of ``self`` under the homomorphism that projects to the classical
            Weyl group factor after rewriting it in either style "PW0" or "W0P".

            EXAMPLES::

                sage: PW0 = ExtendedAffineWeylGroup(['A',2,1])
                sage: s = PW0.S0(); s
                t[Lambdacheck[1] + Lambdacheck[2]] * s1*s2*s1
                sage: s.to_classical_weyl()
                s1*s2*s1

            """
            return self.summand_projection(1)

    class ExtendedAffineWeylGroupPW0(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the translation lattice
        by the finite Weyl group.

        INPUT:

        - `E` -- A parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1],style="PW0")
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Multiplicative form of Coweight lattice of the Root system of type ['A', 2] acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)

        """

        def __init__(self, E):
            # note that we have to use the multiplicative version of the translation lattice
            # and change the twist to deal with this
            def twist(w,l):
                return E.exp_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.exp_lattice(), E.classical_weyl(), twist = twist, act_to_right=False, prefix0=E._prefixt, print_tuple = E._print_tuple, category=E.Realizations())
            self._style = "PW0"

        def _repr_(self):
            return self.realization_of()._repr_() + " realized by " + super(ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPW0, self)._repr_()

        def translation_group_morphism(self, la):
            r"""
            Map the translation lattice element ``la`` into ``self``.

            EXAMPLES::

                sage: PW0 = ExtendedAffineWeylGroup(['A',2,1], style="PW0", translation="tau", print_tuple = True)
                sage: la = PW0.realization_of().lattice().an_element(); la
                2*Lambdacheck[1] + 2*Lambdacheck[2]
                sage: PW0.translation_group_morphism(la)
                (tau[2*Lambdacheck[1] + 2*Lambdacheck[2]], 1)

            """
            E = self.realization_of()
            return self((E.exp_lattice()(la),self.summands()[1].one()))

        @cached_method
        def S0(self):
            """
            Returns the affine simple reflection.

            EXAMPLE::

                sage: ExtendedAffineWeylGroup("B2").S0()
                t[Lambdacheck[2]] * s2*s1*s2

            """
            E = self.realization_of()
            return self((E.exp_lattice()(E.lattice()(E._special_translation)), E._special_reflection))

        @cached_method
        def simple_reflection(self, i):
            r"""
            Returns the `i`-th simple reflection in ``self``.

            EXAMPLES::
               sage: PW0 = ExtendedAffineWeylGroup("G2")
               sage: E = PW0.realization_of()
               sage: [(i, PW0.simple_reflection(i)) for i in E.index_set()]
               [(0, t[Lambdacheck[2]] * s2*s1*s2*s1*s2), (1, s1), (2, s2)]

            """
            if i == 0:
                return self.S0()
            else:
                E = self.realization_of()
                return self.classical_weyl_morphism(E.classical_weyl().simple_reflection(i))

        @cached_method
        def simple_reflections(self):
            r"""
            Returns a family for the simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup("A3").simple_reflections()
                Finite family {0: t[Lambdacheck[1] + Lambdacheck[3]] * s1*s2*s3*s2*s1, 1: s1, 2: s2, 3: s3}

            """
            return Family(self.realization_of().index_set(), lambda i: self.simple_reflection(i))

        def classical_weyl_morphism(self, w):
            r"""
            Returns the image of `w` under the homomorphism of the classical Weyl group into ``self``.

            EXAMPLES::

                sage: PW0 = ExtendedAffineWeylGroup("A3",print_tuple=True)
                sage: E = PW0.realization_of()
                sage: PW0.classical_weyl_morphism(E.classical_weyl().from_reduced_word([1,2]))
                (t[0], s1*s2)

            """
            return self((self.summands()[0].one(),w))

    class ExtendedAffineWeylGroupW0PElement(GroupSemidirectProduct.Element):
        r"""
        The element class for the W0P realization.
        """
        def has_descent(self, i, side='right', positive=False):
            r"""
            INPUT:

            - `i` - an index.

            OPTIONAL:

            - side - 'left' or 'right' (default: 'right')
            - positive - True or False (default: False)

            EXAMPLES::

                sage: We = ExtendedAffineWeylGroup(['A',4,2],style="W0P")
                sage: w = We.from_reduced_word([0,1]); w
                s1*s2 * t[Lambda[1] - Lambda[2]]
                sage: w.has_descent(0, side='left')
                True

            """
            E = self.parent().realization_of()
            if side == 'left':
                self = ~self
            if positive:
                return not self.has_descent(i, side='right')
            w = self.summand_projection(0)
            la = self.summand_projection(1).value
            if i == 0:
                ip = la.scalar(E._special_translation_covector) * E._a0check
                if ip < -1:
                    return True
                if ip > -1:
                    return False
                return E._special_root.weyl_action(w).is_positive_root()
            ip = la.scalar(E._simpleR0[i]) # test height versus simple (co)root
            if ip > 0:
                return True
            if ip < 0:
                return False
            return w.has_descent(i, side='right')

        def to_classical_weyl(self):
            r"""
            Project ``self`` into the classical Weyl group.

            EXAMPLES::

                sage: W0P = ExtendedAffineWeylGroup(['A',2,1],style="W0P")
                sage: x = W0P.simple_reflection(0); x
                s1*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[2]]
                sage: x.to_classical_weyl()
                s1*s2*s1

            """
            return self.summand_projection(0)

        def to_translation_right(self):
            r"""
            Project onto the right (translation) factor in the "W0P" style.

            EXAMPLES::

                sage: W0P = ExtendedAffineWeylGroup(['A',2,1],style="W0P")
                sage: x = W0P.simple_reflection(0); x
                s1*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[2]]
                sage: x.to_translation_right()
                -Lambdacheck[1] - Lambdacheck[2]

            """
            return self.summand_projection(1).value

    class ExtendedAffineWeylGroupW0P(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the finite Weyl group
        by the translation lattice.

        INPUT:

        - `E` -- A parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1], style="W0P")
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice) acting on Multiplicative form of Coweight lattice of the Root system of type ['A', 2]

        """

        def __init__(self, E):
            def twist(w,l):
                return E.exp_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.classical_weyl(), E.exp_lattice(), twist=twist, act_to_right=True, prefix1=E._prefixt, print_tuple=E._print_tuple, category=E.Realizations())
            self._style = "W0P"

        def _repr_(self):
            return self.realization_of()._repr_() + " realized by " + super(ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0P, self)._repr_()

        def S0(self):
            r"""
            Returns the zero-th simple reflection in style "W0P".

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(["A",3,1],style="W0P").S0()
                s1*s2*s3*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[3]]

            """
            E = self.realization_of()
            return self((E._special_reflection,E.exp_lattice()(E.lattice()(-E._special_translation))))

        def simple_reflection(self, i):
            r"""
            Returns the `i`-th simple reflection in ``self``.
            """
            if i == 0:
                return self.S0()
            E = self.realization_of()
            return self.classical_weyl_morphism(E.classical_weyl().simple_reflection(i))

        @cached_method
        def simple_reflections(self):
            r"""
            Returns the family of simple reflections.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(["A",3,1],style="W0P").simple_reflections()
                Finite family {0: s1*s2*s3*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[3]], 1: s1, 2: s2, 3: s3}

            """

            return Family(self.realization_of().index_set(), lambda i: self.simple_reflection(i))

        def classical_weyl_morphism(self, w):
            r"""
            Returns the image of the classical Weyl group element `w` in ``self``.

            EXAMPLES::

                sage: W0P = ExtendedAffineWeylGroup(['A',2,1],style="W0P",print_tuple=True)
                sage: W0P.classical_weyl_morphism(W0P.realization_of().classical_weyl().from_reduced_word([2,1]))
                (s2*s1, t[0])

            """
            return self((w,self.summands()[1].one()))

        def translation_group_morphism(self, la):
            r"""
            Returns the image of the lattice element ``la`` in ``self``.

            EXAMPLES::

                sage: W0P = ExtendedAffineWeylGroup(['A',2,1],style="W0P",print_tuple=True)
                sage: W0P.translation_group_morphism(W0P.realization_of().lattice().an_element())
                (1, t[2*Lambdacheck[1] + 2*Lambdacheck[2]])

            """
            return self((self.summands()[0].one(),self.realization_of().exp_lattice()(la)))

    class ExtendedAffineWeylGroupWFElement(GroupSemidirectProduct.Element):
        r"""
        Element class for the "WF" realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            INPUT:

            - `i` -- an index.

            OPTIONAL:

            - ``side`` -- 'left' or 'right' (default: 'right')
            - ``positive`` -- True or False (default: False)

            """
            E = self.parent().realization_of()
            if side == 'right':
                self = ~self
            if positive:
                return not self.has_descent(i, side='left')
            return self.summand_projection(0).has_descent(i, side='left')

        def to_fundamental_group(self):
            r"""
            Project ``self`` to the right (fundamental group) factor in the "WF" style.

            EXAMPLES::

                sage: WF = ExtendedAffineWeylGroup(['A',2,1],style="WF")
                sage: E = WF.realization_of()
                sage: x = WF.translation_group_morphism(E.lattice_basis()[1]); x
                S0*S2 * pi[1]
                sage: x.to_fundamental_group()
                pi[1]

            """
            return self.summand_projection(1)

        def to_affine_weyl_left(self):
            r"""
            Project ``self`` to the left (affine Weyl group) factor in the "WF" style.

            EXAMPLES::

                sage: WF = ExtendedAffineWeylGroup(['A',2,1],style="WF")
                sage: E = WF.realization_of()
                sage: x = WF.translation_group_morphism(E.lattice_basis()[1]); x
                S0*S2 * pi[1]
                sage: x.to_affine_weyl_left()
                S0*S2

            """
            return self.summand_projection(0)

        def bruhat_le(self, x):
            r"""
            Return whether ``self`` is less than or equal to `x` in the Bruhat order.

            EXAMPLES::

                sage: WF = ExtendedAffineWeylGroup(['A',2,1],style="WF",affine="s", print_tuple=True)
                sage: E = WF.realization_of()
                sage: r = E.affine_weyl().from_reduced_word
                sage: v = r([1,0])
                sage: w = r([1,2,0])
                sage: v.bruhat_le(w)
                True
                sage: vv = WF.affine_weyl_morphism(v); vv
                (s1*s0, pi[0])
                sage: ww = WF.affine_weyl_morphism(w); ww
                (s1*s2*s0, pi[0])
                sage: vv.bruhat_le(ww)
                True
                sage: f = E.fundamental_group()(2); f
                pi[2]
                sage: ff = WF.fundamental_group_morphism(f); ff
                (1, pi[2])
                sage: vv.bruhat_le(ww*ff)
                False
                sage: (vv*ff).bruhat_le(ww*ff)
                True

            """
            if self.summand_projection(1) != x.summand_projection(1):
                return False
            return self.summand_projection(0).bruhat_le(x.summand_projection(0))

    class ExtendedAffineWeylGroupWF(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the affine Weyl group
        by the fundamental group.

        INPUT:

        - `E` -- A parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1], style="WF")
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice) acted upon by Fundamental group of type ['A', 2, 1]            

        """

        def __init__(self, E):
            def twist(g,w):
                return g.act_on_affine_weyl(w)

            GroupSemidirectProduct.__init__(self, E.affine_weyl(), E.fundamental_group(), twist = twist, act_to_right=False, print_tuple = E._print_tuple, category=E.Realizations())
            self._style = "WF"

        def _repr_(self):
            return self.realization_of()._repr_() + " realized by " + super(ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupWF, self)._repr_()

        def affine_weyl_morphism(self, w):
            r"""
            Returns the image of the affine Weyl group element `w` in ``self``.

            EXAMPLES::

                sage: WF = ExtendedAffineWeylGroup(['C',2,1],style="WF",print_tuple=True)
                sage: WF.affine_weyl_morphism(WF.realization_of().affine_weyl().from_reduced_word([1,2,1,0]))
                (S1*S2*S1*S0, pi[0])

            """
            return self((w,self.summands()[1].one()))

        @cached_method
        def simple_reflections(self):
            r"""
            Returns the family of simple reflections.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(["A",3,1],style="WF",affine="r").simple_reflections()
                Finite family {0: r0, 1: r1, 2: r2, 3: r3}

            """
            E = self.realization_of()
            W = E.affine_weyl()
            return Family(E.index_set(), lambda i: self.affine_weyl_morphism(W.simple_reflection(i)))

        @cached_method
        def fundamental_group_morphism(self, f):
            r"""
            Returns the image of `f` under the homomorphism from the fundamental group into
            the right (fundamental group) factor in "WF" style.

            EXAMPLES::

                sage: WF = ExtendedAffineWeylGroup(['E',6,1],style="WF",print_tuple=True)
                sage: F = WF.realization_of().fundamental_group().family()
                sage: [(x,WF.fundamental_group_morphism(x)) for x in F]
                [(pi[0], (1, pi[0])), (pi[1], (1, pi[1])), (pi[6], (1, pi[6]))]

            """
            return self((self.summands()[0].one(),f))

    class ExtendedAffineWeylGroupFWElement(GroupSemidirectProduct.Element):
        r"""
        The element class for the "FW" realization.
        """
        def has_descent(self, i, side='right', positive=False):
            r"""
            INPUT:

            - `i` -- an index.

            OPTIONAL:

            - ``side`` -- 'left' or 'right' (default: 'right')
            - ``positive`` -- True or False (default: False)

            """
            E = self.parent().realization_of()
            if side == 'left':
                self = ~self
            if positive:
                return not self.has_descent(i, side='right')
            return self.summand_projection(1).has_descent(i, side='right')

        def to_fundamental_group(self):
            r"""
            Returns the projection of ``self`` to the fundamental group in the "FW" style.

            EXAMPLES::

                sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW")
                sage: E = FW.realization_of()
                sage: x = FW.translation_group_morphism(E.lattice_basis()[2]); x
                pi[2] * S1*S2
                sage: x.to_fundamental_group()
                pi[2]

            """
            return self.summand_projection(0)

        def to_affine_weyl_right(self):
            r"""
            Project ``self`` to the right (affine Weyl group) factor in the "FW" style.

            EXAMPLES::

                sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW")
                sage: E = FW.realization_of()
                sage: x = FW.translation_group_morphism(E.lattice_basis()[1]); x
                pi[1] * S2*S1
                sage: x.to_affine_weyl_right()
                S2*S1

            """
            return self.summand_projection(1)

        def action_on_affine_roots(self, beta):
            r"""
            Act by ``self`` on the affine root lattice element ``beta``.

            EXAMPLES::

                sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW",affine="s")
                sage: E = FW.realization_of()
                sage: w = FW.an_element(); w
                pi[2] * s0*s1*s2
                sage: Qaf = RootSystem(['A',2,1]).root_lattice()
                sage: v = Qaf.an_element(); v
                2*alpha[0] + 2*alpha[1] + 3*alpha[2]
                sage: w.action_on_affine_roots(v)
                alpha[0] + alpha[1]

            """
            g = self.summand_projection(0)
            w = self.summand_projection(1)
            return g.act_on_affine_lattice(w.action(beta))

    class ExtendedAffineWeylGroupFW(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the affine Weyl group
        by the fundamental group.

        INPUT:

        - `E` -- A parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1],style="FW")
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Fundamental group of type ['A', 2, 1] acting on Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)            

        """

        def __init__(self, E):
            def twist(g,w):
                return g.act_on_affine_weyl(w)

            GroupSemidirectProduct.__init__(self, E.fundamental_group(), E.affine_weyl(), twist = twist, act_to_right=True, print_tuple = E._print_tuple, category=E.Realizations())
            self._style = "FW"

        def _repr_(self):
            return self.realization_of()._repr_() + " realized by " + super(ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupFW, self)._repr_()

        @cached_method
        def simple_reflections(self):
            r"""
            Returns the family of simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',2,1],style="FW",print_tuple=True).simple_reflections()
                Finite family {0: (pi[0], S0), 1: (pi[0], S1), 2: (pi[0], S2)}

            """
            E = self.realization_of()
            W = E.affine_weyl()
            return Family(E.index_set(), lambda i: self.affine_weyl_morphism(W.simple_reflection(i)))

        def affine_weyl_morphism(self, w):
            r"""
            Returns the image of `w` under the map of the affine Weyl group into the right
            (affine Weyl group) factor in the "FW" style.

            EXAMPLES::

                sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW",print_tuple=True)
                sage: E = FW.realization_of()
                sage: FW.affine_weyl_morphism(E.affine_weyl().from_reduced_word([0,2,1]))
                (pi[0], S0*S2*S1)

            """
            return self((self.summands()[0].one(),w))

        @cached_method
        def fundamental_group_morphism(self, f):
            r"""
            Returns the image of the fundamental group element `f` into ``self``.

            EXAMPLES::

                sage: FW = ExtendedAffineWeylGroup(['A',2,1],style="FW",print_tuple=True)
                sage: E = FW.realization_of()
                sage: FW.fundamental_group_morphism(E.fundamental_group()(2))
                (pi[2], 1)

            """
            return self((f,self.summands()[1].one()))

    class ExtendedAffineWeylGroupPvW0Element(GroupSemidirectProduct.Element):
        r"""
        The element class for the "PvW0" realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            Whether `self` has `i` as a descent.

            INPUT:

                - `i` - an index.

            OPTIONAL:

                - side -- 'left' or 'right' (default: 'right')
                - positive -- True or False (default: False)

            EXAMPLES::

                sage: PvW0 = ExtendedAffineWeylGroup(['A',4,2],style="PvW0")
                sage: w = PvW0.from_reduced_word([0,1]); w
                t[Lambda[1]] * s1*s2
                sage: w.has_descent(0, side='left')
                True

            """

            return self.parent().realization_of().PW0()(self).has_descent(i, side=side, positive=positive)

        def dual_action(self, la):
            r"""
            Returns the action of ``self`` on an element ``la`` of the dual version of the translation lattice.

            EXAMPLES::

                sage: PvW0 = ExtendedAffineWeylGroup(['A',2,1],style="PvW0")
                sage: E = PvW0.realization_of()
                sage: x = PvW0.an_element(); x
                t[2*Lambda[1] + 2*Lambda[2]] * s1*s2
                sage: la = E.dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: x.dual_action(la)
                -2*Lambda[1] + 4*Lambda[2]

            """
            w = self.summand_projection(1)
            assert la in w.parent().domain()
            return self.summand_projection(0).value + w.action(la)

        def to_dual_translation_left(self):
            r"""
            The image of ``self`` under the map that projects to the dual translation lattice
            factor after factoring it to the left as in style "PvW0".

            EXAMPLES::

                sage: PvW0 = ExtendedAffineWeylGroup(['A',2,1],style="PvW0")
                sage: s = PvW0.simple_reflection(0); s
                t[Lambda[1] + Lambda[2]] * s1*s2*s1
                sage: s.to_dual_translation_left()
                Lambda[1] + Lambda[2]

            """
            return self.summand_projection(0).value # undo the GroupExp

        def to_dual_classical_weyl(self):
            r"""
            Returns the image of ``self`` under the homomorphism that projects to the dual classical
            Weyl group factor after rewriting it in either style "PvW0" or "W0Pv".

            EXAMPLES::

                sage: PvW0 = ExtendedAffineWeylGroup(['A',2,1],style="PvW0")
                sage: s = PvW0.simple_reflection(0); s
                t[Lambda[1] + Lambda[2]] * s1*s2*s1
                sage: s.to_dual_classical_weyl()
                s1*s2*s1

            """
            return self.summand_projection(1)

            def is_translation(self):
                r"""
                Returns whether ``self`` is a translation element or not.

                EXAMPLES::

                    sage: PvW0 = ExtendedAffineWeylGroup(['A',2,1],style="PvW0")
                    sage: E = PvW0.realization_of()
                    sage: t = PvW0.from_reduced_word([1,2,1,0])
                    sage: t.is_translation()
                    True
                    sage: PvW0.simple_reflection(0).is_translation()
                    False

                """
                w = self.to_dual_classical_weyl()
                return w == w.parent().one()

    class ExtendedAffineWeylGroupPvW0(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the dual form of the translation lattice
        by the finite Weyl group.

        INPUT:

        - `E` -- A parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1],style="PvW0")
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['A', 2] acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)

        """

        def __init__(self, E):
            # note that we have to use the multiplicative version of the translation lattice
            # and change the twist to deal with this
            def twist(w,l):
                return E.exp_dual_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.exp_dual_lattice(), E.dual_classical_weyl(), twist = twist, act_to_right=False, prefix0=E._prefixt, print_tuple = E._print_tuple, category=E.Realizations())
            self._style = "PvW0"

        def _repr_(self):
            return self.realization_of()._repr_() + " realized by " + super(ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPvW0, self)._repr_()

        def dual_translation_group_morphism(self, la):
            r"""
            Map the dual translation lattice element ``la`` into ``self``.

            EXAMPLES::

                sage: PvW0 = ExtendedAffineWeylGroup(['A',2,1], style="PvW0", translation="tau", print_tuple = True)
                sage: la = PvW0.realization_of().dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: PvW0.dual_translation_group_morphism(la)
                (tau[2*Lambda[1] + 2*Lambda[2]], 1)

            """
            E = self.realization_of()
            return self((E.exp_dual_lattice()(la),self.summands()[1].one()))

        @cached_method
        def simple_reflections(self):
            r"""
            Returns a family for the simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',3,1],style="PvW0").simple_reflections()
                Finite family {0: t[Lambda[1] + Lambda[3]] * s1*s2*s3*s2*s1, 1: s1, 2: s2, 3: s3}

            """
            WR = self.realization_of()
            return Family(WR.index_set(), lambda i: self(WR.PW0().simple_reflection(i)))

        def dual_classical_weyl_morphism(self, w):
            r"""
            Returns the image of `w` under the homomorphism of the dual form of the classical Weyl group into ``self``.

            EXAMPLES::

                sage: PvW0 = ExtendedAffineWeylGroup(['A',3,1],style="PvW0",print_tuple=True)
                sage: E = PvW0.realization_of()
                sage: PvW0.dual_classical_weyl_morphism(E.dual_classical_weyl().from_reduced_word([1,2]))
                (t[0], s1*s2)

            """
            return self((self.summands()[0].one(),w))

    class ExtendedAffineWeylGroupW0PvElement(GroupSemidirectProduct.Element):
        r"""
        The element class for the "W0Pv" realization.
        """

        def dual_action(self, la):
            r"""
            Returns the action of ``self`` on an element ``la`` of the dual version of the translation lattice.

            EXAMPLES::

                sage: W0Pv = ExtendedAffineWeylGroup(['A',2,1],style="W0Pv")
                sage: E = W0Pv.realization_of()
                sage: x = W0Pv.an_element(); x
                s1*s2 * t[2*Lambda[1] + 2*Lambda[2]]
                sage: la = E.dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: x.dual_action(la)
                -8*Lambda[1] + 4*Lambda[2]

            """
            w = self.summand_projection(0)
            assert la in w.parent().domain()
            return w.action(self.summand_projection(1).value + la)

        def has_descent(self, i, side='right', positive=False):
            r"""
            INPUT:

            - `i` - an index.

            OPTIONAL:

            - side - 'left' or 'right' (default: 'right')
            - positive - True or False (default: False)

            EXAMPLES::

                sage: W0Pv = ExtendedAffineWeylGroup(['A',4,2],style="W0Pv")
                sage: w = W0Pv.from_reduced_word([0,1]); w
                s1*s2 * t[Lambda[1] - Lambda[2]]
                sage: w.has_descent(0, side='left')
                True

            """
            return self.parent().realization_of().W0P()(self).has_descent(i, side=side, positive=positive)

        def to_dual_translation_right(self):
            r"""
            The image of ``self`` under the map that projects to the dual translation lattice
            factor after factoring it to the right as in style "W0Pv".

            EXAMPLES::

                sage: W0Pv = ExtendedAffineWeylGroup(['A',2,1],style="W0Pv")
                sage: s = W0Pv.simple_reflection(0); s
                s1*s2*s1 * t[-Lambda[1] - Lambda[2]]
                sage: s.to_dual_translation_right()
                -Lambda[1] - Lambda[2]

            """
            return self.summand_projection(1).value # undo the GroupExp

        def to_dual_classical_weyl(self):
            r"""
            Returns the image of ``self`` under the homomorphism that projects to the dual classical
            Weyl group factor after rewriting it in either style "PvW0" or "W0Pv".

            EXAMPLES::

                sage: W0Pv = ExtendedAffineWeylGroup(['A',2,1],style="W0Pv")
                sage: s = W0Pv.simple_reflection(0); s
                s1*s2*s1 * t[-Lambda[1] - Lambda[2]]
                sage: s.to_dual_classical_weyl()
                s1*s2*s1

            """
            return self.summand_projection(0)

            def is_translation(self):
                r"""
                Returns whether ``self`` is a translation element or not.

                EXAMPLES::

                    sage: W0Pv = ExtendedAffineWeylGroup(['A',2,1],style="W0Pv")
                    sage: E = W0Pv.realization_of()
                    sage: t = W0Pv.from_reduced_word([1,2,1,0])
                    sage: t.is_translation()
                    True

                """
                w = self.to_dual_classical_weyl()
                return w == w.parent().one()

    class ExtendedAffineWeylGroupW0Pv(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the finite Weyl group, acting on the
        dual form of the translation lattice.

        INPUT:

        - `E` -- A parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1],style="W0Pv")
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice) acting on Multiplicative form of Weight lattice of the Root system of type ['A', 2]

        """

        def __init__(self, E):
            # note that we have to use the multiplicative version of the translation lattice
            # and change the twist to deal with this
            def twist(w,l):
                return E.exp_dual_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.dual_classical_weyl(), E.exp_dual_lattice(), twist = twist, act_to_right=True, prefix1=E._prefixt, print_tuple = E._print_tuple, category=E.Realizations())
            self._style = "W0Pv"

        def _repr_(self):
            return self.realization_of()._repr_() + " realized by " + super(ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0Pv, self)._repr_()

        def dual_translation_group_morphism(self, la):
            r"""
            Map the dual translation lattice element ``la`` into ``self``.

            EXAMPLES::

                sage: W0Pv = ExtendedAffineWeylGroup(['A',2,1], style="W0Pv", translation="tau", print_tuple = True)
                sage: la = W0Pv.realization_of().dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: W0Pv.dual_translation_group_morphism(la)
                (1, tau[2*Lambda[1] + 2*Lambda[2]])

            """
            E = self.realization_of()
            return self((self.summands()[0].one(),E.exp_dual_lattice()(la)))

        @cached_method
        def simple_reflections(self):
            r"""
            Returns a family for the simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',3,1],style="W0Pv").simple_reflections()
                Finite family {0: s1*s2*s3*s2*s1 * t[-Lambda[1] - Lambda[3]], 1: s1, 2: s2, 3: s3}

            """
            WR = self.realization_of()
            return Family(WR.index_set(), lambda i: self(WR.PW0().simple_reflection(i)))

        def dual_classical_weyl_morphism(self, w):
            r"""
            Returns the image of `w` under the homomorphism of the dual form of the classical Weyl group into ``self``.

            EXAMPLES::

                sage: W0Pv = ExtendedAffineWeylGroup(['A',3,1],style="W0Pv",print_tuple=True)
                sage: E = W0Pv.realization_of()
                sage: W0Pv.dual_classical_weyl_morphism(E.dual_classical_weyl().from_reduced_word([1,2]))
                (s1*s2, t[0])

            """
            return self((w,self.summands()[1].one()))

ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPW0.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPW0Element
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0P.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0PElement
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupWF.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupWFElement
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupFW.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupFWElement
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPvW0.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPvW0Element
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0Pv.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0PvElement
