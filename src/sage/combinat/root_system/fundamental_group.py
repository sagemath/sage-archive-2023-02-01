"""
Fundamental Group of an Extended Affine Weyl Group
"""

#*****************************************************************************
#       Copyright (C) 2013 Mark Shimozono <mshimo at math.vt.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.root_system.cartan_type import CartanType
from sage.categories.groups import Groups
from sage.misc.cachefunc import cached_method
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family
from sage.combinat.root_system.root_system import RootSystem
from sage.categories.finite_sets import FiniteSets

def FundamentalGroupOfExtendedAffineWeylGroup(cartan_type, prefix='pi', general_linear=None):
    r"""
    Factory for the fundamental group of an extended affine Weyl group.

    INPUTS:

    - `cartan_type` -- a Cartan type that is either affine or finite, with the latter being a
      shorthand for the untwisted affinization
    - `prefix` (default: 'pi') -- string that labels the elements of the group
    - `general_linear` -- (default: None, meaning False) In untwisted type A, if True, use the
      universal central extension

    ..RUBRIC::

    Associated to each affine Cartan type `\tilde{X}` is an extended affine Weyl group `E`.
    Its subgroup of length-zero elements is called the fundamental group `F`.
    The group `F` can be identified with a subgroup of the group of automorphisms of the
    affine Dynkin diagram. As such, every element of `F` can be viewed as a permutation of the
    set `I` of affine Dynkin nodes.

    Let `0 \in I` be the distinguished affine node; it is the one whose removal produces the associated
    finite Cartan type (call it `X`). A node `i \in I` is called *special*
    if some automorphism of the affine Dynkin diagram, sends `0` to `i`.
    The node `0` is always special due to the identity automorphism.
    There is a bijection of the set of special nodes with the fundamental group. We denote the image of `i` by
    `\pi_i`. The structure of `F` is determined as follows.

    - `\tilde{X}` is untwisted -- `F` is isomorphic to `P^\vee/Q^\vee` where `P^\vee` and `Q^\vee` are the
    coweight and coroot lattices of type `X`. The group `P^\vee/Q^\vee` consists of the cosets `\omega_i^\vee + Q^\vee`
    for special nodes `i`, where `\omega_0^\vee = 0` by convention. In this case the special nodes `i`
    are the *cominuscule* nodes, the ones such that ``omega_i^\vee(\alpha_j)`` is `0` or `1` for all `j\in I_0 = I \setminus \{0\}`.
    For `i` special, addition by `\omega_i^\vee+Q^\vee` permutes `P^\vee/Q^\vee` and therefore permutes the set of special nodes.
    This permutation extends uniquely to an automorphism of the affine Dynkin diagram.

    - `\tilde{X}` is dual untwisted -- (that is, the dual of `\tilde{X}` is untwisted) `F` is isomorphic to `P/Q`
    where `P` and `Q` are the weight and root lattices of type `X`. The group `P/Q` consists of the cosets
    `\omega_i + Q` for special nodes `i`, where `\omega_0 = 0` by convention. In this case the special nodes `i`
    are the *minuscule* nodes, the ones such that ``\alpha_j^\vee(omega_i)`` is `0` or `1` for all `j \in I_0`.
    For `i` special, addition by `\omega_i+Q` permutes `P/Q` and therefore permutes the set of special nodes.
    This permutation extends uniquely to an automorphism of the affine Dynkin diagram.

    - `\tilde{X}` is mixed -- (that is, not of the above two types) `F` is the trivial group.

    EXAMPLES::

        sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
        sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1]); F
        Fundamental group of type ['A', 3, 1]
        sage: F.cartan_type().dynkin_diagram()
        0
        O-------+
        |       |
        |       |
        O---O---O
        1   2   3
        A3~
        sage: F.special_nodes()
        (0, 1, 2, 3)
        sage: F(1)^2
        pi[2]
        sage: F(1)*F(2)
        pi[3]
        sage: F(3)^(-1)
        pi[1]

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("B3"); F
        Fundamental group of type ['B', 3, 1]
        sage: F.cartan_type().dynkin_diagram()
            O 0
            |
            |
        O---O=>=O
        1   2   3
        B3~
        sage: F.special_nodes()
        (0, 1)

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("C2"); F
        Fundamental group of type ['C', 2, 1]
        sage: F.cartan_type().dynkin_diagram()
        O=>=O=<=O
        0   1   2
        C2~
        sage: F.special_nodes()
        (0, 2)

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("D4"); F
        Fundamental group of type ['D', 4, 1]
        sage: F.cartan_type().dynkin_diagram()
            O 4
            |
            |
        O---O---O
        1   |2  3
            |
            O 0
        D4~
        sage: F.special_nodes()
        (0, 1, 3, 4)
        sage: (F(4), F(4)^2)
        (pi[4], pi[0])

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("D5"); F
        Fundamental group of type ['D', 5, 1]
        sage: F.cartan_type().dynkin_diagram()
          0 O   O 5
            |   |
            |   |
        O---O---O---O
        1   2   3   4
        D5~
        sage: F.special_nodes()
        (0, 1, 4, 5)
        sage: (F(5), F(5)^2, F(5)^3, F(5)^4)
        (pi[5], pi[1], pi[4], pi[0])
        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("E6"); F
        Fundamental group of type ['E', 6, 1]
        sage: F.cartan_type().dynkin_diagram()
                O 0
                |
                |
                O 2
                |
                |
        O---O---O---O---O
        1   3   4   5   6
        E6~
        sage: F.special_nodes()
        (0, 1, 6)
        sage: F(1)^2
        pi[6]

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['D',4,2]); F
        Fundamental group of type ['C', 3, 1]^*
        sage: F.cartan_type().dynkin_diagram()
        O=<=O---O=>=O
        0   1   2   3
        C3~*
        sage: F.special_nodes()
        (0, 3)

    We also implement a fundamental group for `GL_n`. It is defined to be the group of integers, which is the
    covering group of the fundamental group Z/nZ for affine `SL_n`::

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True); F
        Fundamental group of GL(3)
        sage: x = F.an_element(); x
        pi[5]
        sage: x*x
        pi[10]
        sage: x.inverse()
        pi[-5]
        sage: wt = F.cartan_type().classical().root_system().ambient_space().an_element(); wt
        (2, 2, 3)
        sage: x.act_on_classical_ambient(wt)
        (2, 3, 2)
        sage: w = WeylGroup(F.cartan_type(),prefix="s").an_element(); w
        s0*s1*s2
        sage: x.act_on_affine_weyl(w)
        s2*s0*s1

    """
    cartan_type = CartanType(cartan_type)
    if cartan_type.is_finite():
        cartan_type = cartan_type.affine()
    if not cartan_type.is_affine():
        raise NotImplementedError
    if general_linear is True:
        if cartan_type.is_untwisted_affine() and cartan_type.type() == "A":
            return FundamentalGroupGL(cartan_type, prefix)
        else:
            raise ValueError, "General Linear Fundamental group is untwisted type A"
    return FundamentalGroupOfExtendedAffineWeylGroup_Class(cartan_type,prefix)

class FundamentalGroupElement(MultiplicativeGroupElement):
    def __init__(self, parent, x):
        r"""
        This should not be called directly
        """
        if x not in parent.special_nodes():
            raise ValueError, "%s is not a special node"%x
        self._value = x
        MultiplicativeGroupElement.__init__(self, parent)

    def value(self):
        r"""
        Returns the special node which indexes the special automorphism `self`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1], prefix="f")
            sage: x = F.an_element(); x
            f[2]
            sage: x.value()
            2

        """
        return self._value

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1], prefix="f")
            sage: F(2)^3 # indirect doctest
            f[1]

        """
        return self.parent()._prefix + "[" + repr(self.value()) + "]"

    def inverse(self):
        r"""
        Returns the inverse element of `self`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: F(1).inverse()
            pi[3]
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['E',6,1], prefix="f")
            sage: F(1).inverse()
            f[6]

        """
        par = self.parent()
        return self.__class__(par, par.dual_node(self.value()))

    __invert__ = inverse

    def __cmp__(self, x):
        r"""
        Compare `self` with `x`.
        """
        if self.__class__  != x.__class__:
            return cmp(self.__class__,x.__class__)
        return cmp(self.value(), x.value())

    def act_on_affine_weyl(self, w):
        r"""
        Act by `self` on the element `w` of the affine Weyl group.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: W = WeylGroup(F.cartan_type(),prefix="s")
            sage: w = W.from_reduced_word([2,3,0])
            sage: F(1).act_on_affine_weyl(w).reduced_word()
            [3, 0, 1]

        """
        par = self.parent()
        if self == par.one():
            return w
        action = par.action(self.value())
        return w.parent().from_reduced_word([action(j) for j in w.reduced_word()])

    def act_on_affine_lattice(self, wt):
        r"""
        Act by `self` on the element `wt` of an affine root/weight lattice realization.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: la = RootSystem(F.cartan_type()).weight_lattice().an_element(); la
            2*Lambda[0] + 2*Lambda[1] + 3*Lambda[2]
            sage: F(3).act_on_affine_lattice(la)
            2*Lambda[0] + 3*Lambda[1] + 2*Lambda[3]

        warning::

            Doesn't work on ambient spaces.

        """
        return wt.map_support(self.parent().action(self.value()))

class FundamentalGroupOfExtendedAffineWeylGroup_Class(UniqueRepresentation, Parent):
    r"""
    The group of length zero elements in the extended affine Weyl group.
    """
    Element = FundamentalGroupElement

    def __init__(self, cartan_type, prefix):
        def get_the_index(beta):
            r"""
            Given a dictionary with one key, return this key
            """
            supp = beta.support()
            assert len(supp) == 1
            return supp[0]

        self._cartan_type = cartan_type
        self._prefix = prefix

        if cartan_type.dual().is_untwisted_affine():
            cartan_type = cartan_type.dual()
        if cartan_type.is_untwisted_affine():
            cartan_type_classical = cartan_type.classical()
            I = [i for i in cartan_type_classical.index_set()]
            Q = RootSystem(cartan_type_classical).root_lattice()
            alpha = Q.simple_roots()
            theta = Q.highest_root()
            self._special_nodes = tuple([0] + [i for i in I if theta[i] == 1])
            om = RootSystem(cartan_type_classical).weight_lattice().fundamental_weights()
            W = Q.weyl_group(prefix="s")
            w0 = W.long_element()
            fg_dict = {}
            for j in cartan_type.index_set():
                fg_dict[0,j] = j
            inverse_dict = {}
            inverse_dict[0] = 0
            finite_action_dict = {}
            finite_action_dict[0] = tuple([])
            for i in self._special_nodes:
                if i == 0:
                    continue
                antiwt, red = om[i].to_dominant_chamber(reduced_word=True, positive=False)
                finite_action_dict[i] = tuple(red)
                w0i = W.from_reduced_word(red)
                idual = get_the_index(-antiwt)
                inverse_dict[i] = idual
                fg_dict[i,0] = i
                for j in I:
                    if j != idual:
                        fg_dict[i,j] = get_the_index(w0i.action(alpha[j]))
                    else:
                        fg_dict[i,j] = 0
            self._action = Family(self._special_nodes, lambda i: Family(cartan_type.index_set(), lambda j: fg_dict[i,j]))
            self._dual_node = Family(self._special_nodes, lambda i: inverse_dict[i])
            self._finite_action = Family(self._special_nodes, lambda i: finite_action_dict[i])

        if cartan_type.type() == "BC" or cartan_type.dual().type() == "BC":
            self._special_nodes = tuple([0])
            self._action = Family({0:Family(cartan_type.index_set(), lambda i: i)})
            self._dual_node = Family({0:0})
            self._finite_action = Family({0:tuple([])})

        Parent.__init__(self, category = Groups().Commutative())

    def _element_constructor_(self, x):
        if isinstance(x, self.element_class) and x.parent() is self:
            return x
        return self.element_class(self, x)

    @cached_method
    def one(self):
        r"""
        Returns the identity element of the fundamental group.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: F.one()
            pi[0]

        """
        return self._element_constructor_(0)

    def product(self, x, y):
        r"""
        Returns the product of `x` and `y`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: F.special_nodes()
            (0, 1, 2, 3)
            sage: F(2)*F(3)
            pi[1]
            sage: F(1)*F(3)^(-1)
            pi[2]

        """
        return self(self.action(x.value())(y.value()))

    def cartan_type(self):
        return self._cartan_type

    def _repr_(self):
        return "Fundamental group of type %s"%self.cartan_type()

    def special_nodes(self):
        return self._special_nodes

    def family(self):
        r"""
        Returns a family indexed by special nodes whose values are the corresponding fundamental
        group elements.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['E',6,1],prefix="f")
            sage: fam = F.family(); fam
            Finite family {0: f[0], 1: f[1], 6: f[6]}
            sage: fam[0] == F(0)
            True
            sage: fam[6]^2
            f[1]

        """
        return Family(self.special_nodes(), lambda i: self(i))

    @cached_method
    def index_set(self):
        return self.cartan_type().index_set()

    def action(self, i):
        r"""
        Returns a function which permutes the affine Dynkin node set by the `i`-th special automorphism.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1])
            sage: [[(i, j, F.action(i)(j)) for j in F.index_set()] for i in F.special_nodes()]
            [[(0, 0, 0), (0, 1, 1), (0, 2, 2)], [(1, 0, 1), (1, 1, 2), (1, 2, 0)], [(2, 0, 2), (2, 1, 0), (2, 2, 1)]]
            sage: G = FundamentalGroupOfExtendedAffineWeylGroup(['D',4,1])
            sage: [[(i, j, G.action(i)(j)) for j in G.index_set()] for i in G.special_nodes()]
            [[(0, 0, 0), (0, 1, 1), (0, 2, 2), (0, 3, 3), (0, 4, 4)], [(1, 0, 1), (1, 1, 0), (1, 2, 2), (1, 3, 4), (1, 4, 3)], [(3, 0, 3), (3, 1, 4), (3, 2, 2), (3, 3, 0), (3, 4, 1)], [(4, 0, 4), (4, 1, 3), (4, 2, 2), (4, 3, 1), (4, 4, 0)]]
        """
        return lambda j: self._action[i][j]

    def dual_node(self, i):
        r"""
        Returns the node that indexes the inverse of the `i`-th element.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1])
            sage: [(i, F.dual_node(i)) for i in F.special_nodes()]
            [(0, 0), (1, 4), (2, 3), (3, 2), (4, 1)]
            sage: G = FundamentalGroupOfExtendedAffineWeylGroup(['E',6,1])
            sage: [(i, G.dual_node(i)) for i in G.special_nodes()]
            [(0, 0), (1, 6), (6, 1)]
            sage: H = FundamentalGroupOfExtendedAffineWeylGroup(['D',5,1])
            sage: [(i, H.dual_node(i)) for i in H.special_nodes()]
            [(0, 0), (1, 1), (4, 5), (5, 4)]

        """
        return self._dual_node[i]

    def finite_action(self, i):
        r"""
        Returns a reduced word for the finite Weyl group element associated with the `i`-th special automorphism.

        More precisely, for each special node `i`, `self.finite_action(i)` is a reduced word for
        the element `v` in the finite Weyl group such that in the extended affine Weyl group,
        the `i`-th special automorphism is equal to `t v` where `t` is a translation element.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: [(i, F.finite_action(i)) for i in F.special_nodes()]
            [(0, ()), (1, (1, 2, 3)), (2, (2, 1, 3, 2)), (3, (3, 2, 1))]

        """
        return self._finite_action[i]

from sage.rings.finite_rings.integer_mod import Mod

class FundamentalGroupGLElement(FundamentalGroupElement):
    def act_on_classical_ambient(self, wt):
        return wt.map_support(self.parent().action(self.value()))

class FundamentalGroupGL(FundamentalGroupOfExtendedAffineWeylGroup_Class):
    r"""
    Fundamental group of `GL_n`. It is just the integers with extra privileges.
    """
    Element = FundamentalGroupGLElement

    def __init__(self, cartan_type, prefix='pi'):
        FundamentalGroupOfExtendedAffineWeylGroup_Class.__init__(self, cartan_type, prefix)
        from sage.rings.integer_ring import ZZ
        self._special_nodes = ZZ
        self._n = cartan_type.n + 1

    def _element_constructor_(self, x):
        if isinstance(x, self.element_class) and x.parent() is self:
            return x
        return self.element_class(self, x)

    @cached_method
    def one(self):
        r"""
        Returns the identity element of the fundamental group.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).one()
            pi[0]

        """
        from sage.rings.integer_ring import ZZ
        return self(ZZ(0))

    def product(self, x, y):
        r"""
        Returns the product of `x` and `y`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.special_nodes()
            Integer Ring
            sage: F(2)*F(3)
            pi[5]
            sage: F(1)*F(3)^(-1)
            pi[-2]

        """
        return self(x.value()+y.value())

    def _repr_(self):
        return "Fundamental group of GL(%s)"%self._n

    def family(self):
        return NotImplemented

    @cached_method
    def an_element(self):
        r"""
        An element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).an_element()
            pi[5]

        """
        from sage.rings.integer_ring import ZZ
        return self(ZZ(5))

    def action(self, i):
        r"""
        The action of the `i`-th automorphism on the affine Dynkin node set.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.action(4)(2)
            0
            sage: F.action(-4)(2)
            1

        """
        return lambda j: Mod(i + j, self._n)

    def dual_node(self, i):
        r"""
        The node whose special automorphism is inverse to that of `i`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.dual_node(2)
            -2
        """
        return -i

    @cached_method
    def finite_action(self, i):
        r"""
        A reduced word for the finite permutation part of the
        special automorphism indexed by `i`.

        More precisely, return a reduced word for the finite Weyl group element `u`
        where `i`-th automorphism (expressed in the extended affine Weyl group)
        has the form `t u` where `t` is a translation element.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.finite_action(10)
            (1, 2)

        """
        i = Mod(i, self._n)
        if i == 0:
            return tuple([])
        om = self.cartan_type().classical().root_system().weight_lattice().fundamental_weight(i)
        return tuple((-om).reduced_word())

