from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP
from sage.groups.group import AbelianGroup as AbelianGroupBase
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.groups import Groups

def AbelianGroupGap(generator_orders):
    r"""
    Create the multiplicative abelian group with given orders of generators.
    
    INPUT:
    
    - ``generator_orders`` -- a list of nonnegative integers where `0` will give a factor isomorphic to `\ZZ`.
    
    OUTPUT:
    
    - an abelian group 
    
    EXAMPLES::
    
            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: AbelianGroupGap([3,6])    
            Abelian group with gap, generator orders (3, 6)
            sage: AbelianGroupGap([3,6,5])
            Abelian group with gap, generator orders (3, 6, 5)
            sage: AbelianGroupGap([3,6,0])      # optional gap_packages
            Abelian group with gap, generator orders (3, 6, 0)
    """
    generator_orders = tuple(ZZ(e) for e in generator_orders)
    if not all([e >= 0 for e in generator_orders]):
        return ValueError("Generator orders must be nonnegative")
    category = Groups().Commutative()
    if 0 in generator_orders:
        category = category.Finite().Enumerated()
    else:
        category = category.Infinite()
    polycyclic_package = libgap.LoadPackage("Polycyclic")
    return AbelianGroupAmbient_gap(generator_orders, polycyclic_package=polycyclic_package, category=None)

class AbelianGroupElement_gap(ElementLibGAP):
    r"""
    An element of an abelian group via libgap
    """
    def __init__(self, parent, x, check=True):
        """
        The Python constructor.

        See :class:`AbelianGroupElement_gap` for details.
        
        INPUT:
        
        - ``parent`` -- an instance of :class:`AbelianGroup_gap`
        - ``x`` -- an instance of :class:`sage.libs.gap.element.GapElement`
        - ``check`` -- boolean (default: ``True``) check 
          if ``x`` is an element  of the group

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,6])
            sage: g = G.an_element()
            sage: TestSuite(g).run()
        """
        if check:
            if not x in parent.gap():
                raise ValueError("%s is not in the group %s" %(x, parent))
        ElementLibGAP.__init__(self, parent, x)

    def __hash__(self):
        r"""
        Return the hash of this element.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,4])
            sage: g = G.an_element()
            sage: g.__hash__()    # random
            1693277541873681615
        """
        return hash(self.parent()) ^ hash(self.exponents())

    def __reduce__(self):
        r"""
        Implement pickling.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,4])
            sage: g = G.an_element()
            sage: g == loads(dumps(g))
            True
            sage: g.__reduce__()
            (Abelian group with gap, generator orders (3, 2, 4), ((1, 1, 1),))
        """
        return self.parent(), (self.exponents(),)

    def _repr_(self):
        """
        The string representation of this element.
        
        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,4])
            sage: g = G.an_element()
            sage: g._repr_()    
            'g1*g2*g3'
        """
        rep = self.gap()._repr_()
        return rep.replace('f','g')

    def exponents(self):
        r"""
        Return the tuple of exponents.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4,7,9])
            sage: gens = G.gens()
            sage: g = gens[0]^2 * gens[1]^4 * gens[2]^8
            sage: g.exponents()
            (2, 4, 8)
            sage: G = AbelianGroupGap([4,7,0])         # optional - gap_packages
            sage: gens = G.gens()
            sage: g = gens[0]^2 * gens[1]^4 * gens[2]^8
            sage: g.exponents()
            (2, 4, 8)
        """
        if self.parent()._with_pc:
            exp = self.gap().Exponents().sage()
        else:
            # works only for small groups
            # as gap has problems to solve the word problem
            P = self.parent()
            x = libgap.Factorization(P.gap(), self.gap())
            L = x.ExtRepOfObj().sage()
            Lgens = L[::2]
            Lexpo = L[1::2]
            exp = []
            orders = P.gens_orders()
            i = 0
            for k in range(len(P.gens())):
                if not k+1 in Lgens:
                    exp.append(0)
                else:
                    i = Lgens.index(k+1)
                    exp.append(Lexpo[i] % orders[k])
        return tuple(exp)

    def order(self):
        r"""
        Return the order of this element.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4])
            sage: g = G.gens()[0]
            sage: g.order()
            4
            sage: G = AbelianGroupGap([0])          # optional - gap_packages
            sage: g = G.gens()[0]
            sage: g.order()                         # optional - gap_packages
            +Infinity
        """
        return self.gap().Order().sage()

#    def sage(self):
#        r"""
#        Convert this element to the corresponding abelian group in sage.
#
#        EXAMPLES::
#
#            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
#            sage: G = AbelianGroupGap([2,3,7,4])
#            sage: all([a == a.gap().sage() for a in A])
#            True
#            sage: all([g == g.sage().gap() for g in G])
#            True
#        """
#        P = self.parent()
#        gens_sage = P._A.gens()
#        e = P._A.identity()
#        exp = self.exponents()
#        for i in range(len(exp)):
#            e *= gens_sage[i]**exp[i]
#        return e

class AbelianGroup_gap(UniqueRepresentation,GroupMixinLibGAP, ParentLibGAP, AbelianGroupBase):
    r"""
    Python wrapper for finitely generated abelian groups in gap.

    Needs the gap package "Polycyclic" in case the group is infinite.

    INPUT:

    - ``G`` -- (default:``None``) a gap group
    - ``ambient`` -- (default:``None``) an :class:`AbelianGroup_gap`
    - ``polycyclic_package`` -- (default: ``False``) boolean
    - ``category` -- a category`

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: G = AbelianGroupGap([3,2,5])
        sage: G    
        Abelian group with gap, generator orders (3, 2, 5)

    """
    def __init__(self, G, ambient=None, polycyclic_package=False, category=None):
        r"""
        Create an instance of this class.

        See :class:`AbelianGroup_gap` for details

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,5])
            sage: TestSuite(G).run()
        """
        self._with_pc = polycyclic_package
        AbelianGroupBase.__init__(self, category=category)
        ParentLibGAP.__init__(self, G, ambient=ambient)

    Element = AbelianGroupElement_gap

    def __reduce__(self):
        r"""
        """
        if self is self.ambient():
            return AbelianGroupGap, (self.gens_orders(),)
        else:
            return self.ambient().subgroup, (self.gens(),)

    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` canonically coerces to ``self``.

        INPUT:

        - ``S`` -- anything.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: G._coerce_map_from_(S)
            True
            sage: S._coerce_map_from_(G)
        """
        try:
            if S.ambient() is self:
                return True
        except AttributeError:
            pass

    def _element_constructor_(self,x,check=True):
        r"""
        Defines coercions and conversions.
        
        INPUT:
        
        - ``x` -- an element of this group, a gap element `
        """
        if isinstance(x, AbelianGroupElement_gap):
            x = x.gap()
        elif x == 1:
            x = self.gap().Identity()
        elif not isinstance(x, GapElement):
            try: 
                exp = x.exponents()
            except AttributeError:
                exp = tuple(x) 
            # turn the exponents into a gap element
            gens_gap = self.gens()
            if len(exp) != len(gens_gap):
                raise ValueError("Input does not match the number of generators.")
            x = gens_gap[0]**0
            for i in range(len(exp)):
                x *= gens_gap[i]**exp[i]
            x = x.gap()
        return self.element_class(self, x, check=check)

    def is_trivial(self):
        r"""
        Return if this group is the trivial group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([])
            sage: G
            Abelian group with gap, generator orders ()
        """
        return 1 == self.order()

    def identity(self):
        r"""
        Return the identity element of this group.

        EXAMPLES:

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4,10])
            sage: G.identity()
            id
        """
        return self(self.gap().Identity())

    @cached_method
    def elementary_divisors(self):
        r"""
        Return the elementary divisors of the group.

        See :meth:`sage.groups.abelian_gps.abelian_group_gap.elementary_divisors`

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: G.elementary_divisors()
            (2, 60)
        """
        ediv = self.gap().AbelianInvariants().sage()
        from sage.matrix.constructor import diagonal_matrix
        ed = diagonal_matrix(ZZ, ediv).elementary_divisors()
        return tuple(d for d in ed if d!=1)

    @cached_method
    def exponent(self):
        r"""
        Return the exponent of this abelian group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,7]) 
            sage: G
            Abelian group with gap, generator orders (2, 3, 7)
            sage: G = AbelianGroupGap([2,4,6])
            sage: G    
            Abelian group with gap, generator orders (2, 4, 6)
            sage: G.exponent()
            12
        """
        return self.gap().Exponent().sage()

    @cached_method
    def gens_orders(self):
        r"""
        Return the orders of the generators.

        Use :meth:`elementary_divisors` if you are looking for an
        invariant of the group.

        OUTPUT:

        - a tuple of integers

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: Z2xZ3 = AbelianGroupGap([2,3])
            sage: Z2xZ3.gens_orders()
            (2, 3)
            sage: Z2xZ3.elementary_divisors()
            (6,)
            sage: Z6 = AbelianGroupGap([6])
            sage: Z6.gens_orders()
            (6,)
            sage: Z6.elementary_divisors()
            (6,)
            sage: Z2xZ3.is_isomorphic(Z6)
            True
            sage: Z2xZ3 is Z6
            False
        """
        from sage.rings.infinity import Infinity
        orders = []
        for g in self.gens():
            order = g.order()
            if order == Infinity:
                order = 0
            orders.append(order)
        return tuple(orders)

    def subgroup(self, gens):
        r"""
        Return the subgroup of this group generated by ``gens``.

        INPUT:

        - ``gens`` -- a list of elements coercible into this group

        OUTPUT:

        - a subgroup

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: S    
            Subgroup of Abelian group with gap, generator orders (2, 3, 4, 5) generated by (g1, g2)
            sage: g = G.an_element()
            sage: s = S.an_element()
            sage: g*s
            g2^2*g3*g4
            sage: G = AbelianGroupGap([3,4,0,2])     # optional - gap_packages
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: g = G.an_element()
            sage: s = S.an_element()
            sage: g*s                       # optional - gap_packages
            g1^2*g2^2*g3*g4

        TESTS::

            sage: h = G.gens()[3]
            sage: h in S
            False
        """
        gens = tuple(self(g) for g in gens)
        return AbelianGroupSubgroup_gap(self, gens)
    
class AbelianGroupAmbient_gap(AbelianGroup_gap):
    r"""
    Ambient abelian groups with gap.
    
    Do not use this class directly. Instead use :function:`AbelianGroupGap`.
    Needs the gap package "Polycyclic" in case the group is infinite.

    INPUT:

    - ``generator_orders`` - a tuple of nonnegative integers
    - ``polycyclic_package`` -- (default: ``False``) boolean
    - ``category` -- a category`
    
    EXAMPLES::
    
        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupAmbient_gap
        sage: AbelianGroupAmbient_gap((2,3,4))
        Abelian group with gap, generator orders (2, 3, 4)
    """
    def __init__(self, generator_orders, polycyclic_package=False, category=None):
        if polycyclic_package:
            G = libgap.eval("AbelianPcpGroup(%s)"%list(generator_orders))
        else:
            G = libgap.AbelianGroup(generator_orders)
        AbelianGroup_gap.__init__(self, G, ambient=None, polycyclic_package=polycyclic_package, category=category)
        
    def _latex_(self):
        """
        Return the latex representation of this group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,6])
            sage: G._latex_()    
            'Abelian group with gap, generator orders $(2, 6)$'
        """
        return "Abelian group with gap, generator orders $" + str(self.gens_orders()) + "$"

    def _repr_(self):
        r"""
        Return the string representation of this group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,6])
            sage: G._repr_()    
            'Abelian group with gap, generator orders (2, 6)'
        """
        return "Abelian group with gap, generator orders " + str(self.gens_orders())

class AbelianGroupSubgroup_gap(AbelianGroup_gap):
    r"""
    Subgroups of abelian groups with gap.
    
    Do not use this class directly. Instead use :meth:`subgroup`.
    
    INPUT:
    
    - ``ambient`` -- the ambient group
    - ``gens`` -- generators of the subgroup
    
    EXAMPLES::
        
            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)                # indirect doctest
    """
    def __init__(self, ambient, gens):
        r"""
        Initialize this module
        
        TESTS::
        
            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap,AbelianGroupSubgroup_gap
            sage: G = AbelianGroupGap([])
            sage: gen = G.gens()
            sage: AbelianGroupSubgroup_gap(G, gen)
            Subgroup of Abelian group with gap, generator orders () generated by ()
        """
        polycyclic_package = ambient._with_pc
        category = ambient.category()
        gens_gap = tuple(g.gap() for g in gens)
        G = ambient.gap().Subgroup(gens_gap)
        AbelianGroup_gap.__init__(self, G, ambient=ambient, polycyclic_package=polycyclic_package, category=category)
        
    def __repr__(self):
        r"""
        Return the string representation of this subgroup.
        
        EXAMPLES::
        
            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)                # indirect doctest
            sage: S.__repr__()      
            'Subgroup of Abelian group with gap, generator orders (2, 3, 4, 5) generated by (g1, g2)'
        """
        s = "Subgroup of %s generated by %s"%(self.ambient(),self.gens())
        return s
