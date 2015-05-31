from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.groups import Groups
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.morphism import SetMorphism
from sage.categories.functor import Functor
from sage.categories.homset import Hom
from sage.structure.element_wrapper import ElementWrapper


class GroupExp(Functor):
    r"""
    A functor that wraps a commutative additive group to become a
    multiplicative group.

    EXAMPLES::

        sage: E = GroupExp(); E
        Functor from Category of commutative additive groups to Category of groups
        sage: EZ = E(ZZ); EZ
        Multiplicative form of Integer Ring
        sage: x = EZ(-3); x
        -3
        sage: x.parent()
        Multiplicative form of Integer Ring
        sage: EZ(-1)*EZ(6) == EZ(5)
        True
        sage: EZ(3)^(-1)
        -3
        sage: EZ.one()
        0
        sage: L = RootSystem(['A',2]).ambient_space()
        sage: EL = E(L)
        sage: W = L.weyl_group(prefix="s")
        sage: s2 = W.simple_reflection(2)
        sage: def my_action(mu):
        ....:     return s2.action(mu)
        sage: from sage.categories.morphism import SetMorphism
        sage: from sage.categories.homset import Hom
        sage: f = SetMorphism(Hom(L,L,CommutativeAdditiveGroups()), my_action)
        sage: F = E(f); F
        Generic endomorphism of Multiplicative form of Ambient space of the Root system of type ['A', 2]
        sage: v = L.an_element(); v
        (2, 2, 3)
        sage: y = F(EL(v)); y
        (2, 3, 2)
        sage: y.parent()
        Multiplicative form of Ambient space of the Root system of type ['A', 2]

    """
    def __init__(self):
        Functor.__init__(self, CommutativeAdditiveGroups(), Groups())

    def _apply_functor(self, x):
        r"""
        Given a commutative additive group `x`, returns the isomorphic
        multiplicative group.

        EXAMPLES::

            sage: EQ2 = GroupExp()(QQ^2)
            sage: x = EQ2(vector(QQ,(-2,1))); x
            (-2, 1)
            sage: x^(-1)
            (2, -1)
            sage: x*x
            (-4, 2)
            sage: EQ2(vector(QQ,(-1,1)))*EQ2(vector(QQ,(3,4))) == EQ2(vector(QQ,(2,5)))
            True
            sage: EQ2.one()
            (0, 0)
        """
        return GroupExp_Class(x)

    def _apply_functor_to_morphism(self, f):
        r"""
        Given a morphism of commutative additive groups, returns the corresponding morphism
        of multiplicative groups.

        EXAMPLES::

            sage: def double(x):
            ....:     return x + x
            sage: from sage.categories.morphism import SetMorphism
            sage: from sage.categories.homset import Hom
            sage: f = SetMorphism(Hom(ZZ,ZZ,CommutativeAdditiveGroups()),double)
            sage: E = GroupExp()
            sage: EZ = E._apply_functor(ZZ)
            sage: F = E._apply_functor_to_morphism(f)
            sage: F.domain() == EZ
            True
            sage: F.codomain() == EZ
            True
            sage: F(EZ(3)) == EZ(3)*EZ(3)
            True
        """
        new_domain = self._apply_functor(f.domain())
        new_codomain = self._apply_functor(f.codomain())
        new_f = lambda a: new_codomain(f(a.value))
        return SetMorphism(Hom(new_domain, new_codomain, Groups()), new_f)


class GroupExpElement(ElementWrapper, MultiplicativeGroupElement):
    r"""
    The Element class for a GroupExp_Class object.

    INPUT:

    - ``self`` -- the instance being created
    - ``parent`` -- the parent of ``self``
    - ``x`` -- the additive group element being wrapped

    EXAMPLES::

        sage: G = QQ^2
        sage: EG = GroupExp()(G)
        sage: x = GroupExpElement(EG, vector(QQ, (1,-3))); x
        (1, -3)
        sage: x.parent()
        Multiplicative form of Vector space of dimension 2 over Rational Field
    """
    def __init__(self, parent, x):
        if x not in parent._G:
            return ValueError("%s is not an element of %s" % (x, parent._G))
        ElementWrapper.__init__(self, parent, x)

    def inverse(self):
        r"""
        Returns the inverse of the element `self`.

        EXAMPLES::

            sage: EZ = GroupExp()(ZZ)
            sage: EZ(-3)^(-1)
            3
        """
        return GroupExpElement(self.parent(), -self.value)

    __invert__ = inverse

    def __mul__(self, x):
        r"""
        Returns the product of `self` and `x`.

        EXAMPLES::

            sage: G = GroupExp()(ZZ)
            sage: x = G(2)
            sage: x.__mul__(G(3))
            5
            sage: G.product(G(2),G(3))
            5
        """
        return GroupExpElement(self.parent(), self.value + x.value)


class GroupExp_Class(UniqueRepresentation, Parent):
    r"""
    The multiplicative form of a commutative additive group.

    INPUT:

    - `G`: a commutative additive group
    """
    def __init__(self, G):
        if not G in CommutativeAdditiveGroups():
            raise TypeError("%s must be a commutative additive group" % G)
        self._G = G
        Parent.__init__(self, category=Groups())

    def _repr_(self):
        r"""
        Returns a string describing ``self``.

        EXAMPLES::

            sage: GroupExp()(ZZ) # indirect doctest
            Multiplicative form of Integer Ring
        """
        return "Multiplicative form of %s" % self._G

    def _element_constructor_(self, x):
        r"""
        Constructs the element of ``self`` that wraps the additive
        group element `x`.

        EXAMPLES::

            sage: G = GroupExp()(ZZ)
            sage: G(4) # indirect doctest
            4
        """
        return GroupExpElement(self, x)

    def one(self):
        r"""
        Returns the identity element of ``self``.

        EXAMPLES::

            sage: G = GroupExp()(ZZ^2)
            sage: G.one()
            (0, 0)
        """
        return GroupExpElement(self, self._G.zero())

    def an_element(self):
        r"""
        Returns an element of ``self``.

        EXAMPLES::

            sage: L = RootSystem(['A',2]).weight_lattice()
            sage: EL = GroupExp()(L)
            sage: x = EL.an_element(); x
            2*Lambda[1] + 2*Lambda[2]
            sage: x.parent()
            Multiplicative form of Weight lattice of the Root system of type ['A', 2]
        """
        return GroupExpElement(self, self._G.an_element())

    def product(self, x, y):
        r"""
        Returns the product of `x` and `y` in ``self``.

        EXAMPLES::

            sage: G = GroupExp()(ZZ)
            sage: G.product(G(2),G(7))
            9
            sage: x = G(2)
            sage: x.__mul__(G(7))
            9
        """
        return GroupExpElement(self, x.value + y.value)


GroupExp_Class.Element = GroupExpElement
