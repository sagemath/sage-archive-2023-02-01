from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP
from sage.groups.group import AbelianGroup as AbelianGroupBase
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.arith.all import GCD, LCM
from sage.structure.unique_representation import UniqueRepresentation



class AbelianGroupElement_gap(ElementLibGAP):
    r"""
    An element of an abelian group via libgap
    """
    def __init__(self, parent, x, check=True):
        """
        The Python constructor.

        See :class:`AbelianGroupElement_gap` for details.

        TESTS::

            sage: A = AbelianGroup([3,6])
            sage: G = A.gap()
            sage: g = G.an_element()
            sage: TestSuite(g).run()
        """
        if isinstance(x, AbelianGroupElement_gap):
            x = x.gap()
        if not isinstance(x, GapElement):
            A = parent._A
            x = A(x)
            # turn this into a gap element
            gens_gap = parent.gens()
            exp = x.exponents()
            x = gens_gap[0]**0
            for i in range(len(exp)):
                x *= gens_gap[i]**exp[i]
            x = x.gap()
        if check:
            if not x in parent.gap():
                raise ValueError("%s is not in the group %s" %(x, parent))
        ElementLibGAP.__init__(self, parent, x)

    def __hash__(self):
        r"""
        Return the hash of this element

        EXAMPLES::

            sage: A = AbelianGroup([3,2,4])
            sage: G = A.gap()
            sage: g = G.an_element()
            sage: g.__hash__()    # random
            1693277541873681615
        """
        return hash(self.parent()) ^ hash(self.exponents())

    def __reduce__(self):
        r"""
        Implement pickling

            sage: A = AbelianGroup([3,2,4])
            sage: G = A.gap()
            sage: g = G.an_element()
            sage: g == loads(dumps(g))
            True
            sage: g.__reduce__()
            (Multiplicative Abelian group isomorphic to C3 x C2 x C4 with gap, (f0*f1*f2,))
        """
        return self.parent(), (self.sage(),)

    def exponents(self):
        r"""
        Return the tuple of exponents.

        EXAMPLES::

            sage: A = AbelianGroup([4,7,9])
            sage: G = A.gap()
            sage: gens = G.gens()
            sage: g = gens[0]^2 * gens[1]^4 * gens[2]^8
            sage: g.exponents()
            (2, 4, 8)
            sage: A = AbelianGroup([4,7,0])         # optional - gap_packages
            sage: G = A.gap()
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

            sage: G = AbelianGroup([4]).gap()
            sage: g = G.gens()[0]
            sage: g.order()
            4
            sage: G = AbelianGroup([0]).gap()       # optional - gap_packages
            sage: g = G.gens()[0]
            sage: g.order()
            +Infinity
        """
        return self.gap().Order().sage()

    def sage(self):
        r"""
        Convert this element to the corresponding abelian group in sage.

        EXAMPLES::

            sage: A = AbelianGroup([2,3,7,4])
            sage: G = A.gap()
            sage: all([a == a.gap().sage() for a in A])
            True
            sage: all([g == g.sage().gap() for g in G])
            True
        """
        P = self.parent()
        gens_sage = P._A.gens()
        e = P._A.identity()
        exp = self.exponents()
        for i in range(len(exp)):
            e *= gens_sage[i]**exp[i]
        return e

class AbelianGroup_gap(UniqueRepresentation, GroupMixinLibGAP, ParentLibGAP, AbelianGroupBase):
    r"""
    Class for finitely generated abelian groups implemented with gap.

    Needs the gap package "Polycyclic" in case the group is infinite.

    INPUT:

        - ``A`` -- an `AbelianGroup`
        - ``G`` -- (default:``None``) a gap group
        - ``ambient`` -- (default:``None``) an AbelianGroup_gap

    EXAMPLES::

        sage: A = AbelianGroup([3,2,5])
        sage: G = A.gap()
        sage: TestSuite(G).run()
    """
    def __init__(self, A, G=None, ambient=None):
        r"""
        """
        AbelianGroupBase.__init__(self, category=A.category())
        self._with_pc = libgap.LoadPackage("Polycyclic")
        if G is None:
            if self._with_pc:
                G = libgap.eval("AbelianPcpGroup(%s)"%list(A.gens_orders()))
            else:
                G = libgap(A)
        ParentLibGAP.__init__(self, G, ambient=ambient)
        self._A = A

    Element = AbelianGroupElement_gap

    def _latex_(self):
        """
        Return the latex representation of this group.

        EXAMPLES::

            sage: A = AbelianGroup([2,6])
            sage: G = A.gap()
            sage: G._latex_()
            '$\\mathrm{AbelianGroup}( 2, (2, 6) )$ with gap'
        """
        return self._A._latex_() + " with gap"

    def _repr_(self):
        r"""
        Return the string representation of this group.

        EXAMPLES::

            sage: A = AbelianGroup([2,6])
            sage: G = A.gap()
            sage: G._repr_()
            'Multiplicative Abelian group isomorphic to C2 x C6 with gap'
        """
        s = self._A._repr_()
        s += " with gap"
        return s

    def __hash__(self):
        r"""
        A hash function.

        EXAMPLES::

            sage: A = AbelianGroup([2,6])
            sage: G = A.gap()
            sage: G.__hash__()      # random
            -9223363266866470866
        """
        return hash(self._A) ^ hash(type(self))

    def _coerce_map_from_(self, S):
        r"""
        """
        try:
            if S.ambient() is self:
                return True
        except AttributeError:
            pass
        if self._A is S:
            return True

    def is_trivial(self):
        r"""
        Return if this group is the trivial group.

        EXAMPLES::

            sage: A = AbelianGroup([]).gap()
            sage: A.is_trivial()
            True
        """
        return 1 == self.order()

    def identity(self):
        r"""
        Return the identity element of this group.

        EXAMPLES:

            sage: G = AbelianGroup([4,10]).gap()
            sage: G.identity()
            1
        """
        return self(self.gap().Identity())

    @cached_method
    def elementary_divisors(self):
        r"""
        Return the elementary divisors of the group.

        See :meth:`sage.groups.abelian_gps.abelian_group_gap.elementary_divisors`
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

            sage: G = AbelianGroup([2,3,7]).gap() 
            sage: G
            Multiplicative Abelian group isomorphic to C2 x C3 x C7 with gap
            sage: G.exponent()
            42
            sage: G = AbelianGroup([2,4,6]).gap()
            sage: G
            Multiplicative Abelian group isomorphic to C2 x C4 x C6 with gap
            sage: G.exponent()
            12
        """
        return self.gap().Exponent().sage()

    @cached_method
    def gens_orders(self):
        r"""
        Return the orders of the cyclic factors that this group has
        been defined with.

        Use :meth:`elementary_divisors` if you are looking for an
        invariant of the group.

        OUTPUT:

            - A tuple of integers.

        EXAMPLES::

            sage: Z2xZ3 = AbelianGroup([2,3]).gap()
            sage: Z2xZ3.gens_orders()
            (2, 3)
            sage: Z2xZ3.elementary_divisors()
            (6,)
            sage: Z6 = AbelianGroup([6]).gap()
            sage: Z6.gens_orders()
            (6,)
            sage: Z6.elementary_divisors()
            (6,)
            sage: Z2xZ3.is_isomorphic(Z6)
            True
            sage: Z2xZ3 is Z6
            False
        """
        return tuple(g.order() for g in self.gens())

    def sage(self):
        r"""
        Return the sage pendant of this abelian group.

        EXAMPLES::

            sage: A = AbelianGroup([])
            sage: G = A.gap()
            sage: A is G.sage()
            True
        """
        return self._A

    def subgroup(self, gens):
        r"""
        Return the subgroup of this group generated by ``gens``.

        INPUT:

            - ``gens`` -- a list of elements coercible into this group

        OUTPUT:

            - a subgroup which remembers that it is a subgroup

        EXAMPLES::

            sage: A = AbelianGroup([2,3,4,5])
            sage: G = A.gap()
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: S
            Multiplicative Abelian subgroup isomorphic to C2 x C3 generated by {f0, f1} with gap
            sage: g = G.an_element()
            sage: s = S.an_element()
            sage: a = A.an_element()
            sage: g*s
            g2^2*g3*g4
            sage: g*a
            g2^2*g3^2*g4^2
            sage: A = AbelianGroup([3,4,0,2])
            sage: G = A.gap()
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: g = G.an_element()
            sage: s = S.an_element()
            sage: a = A.an_element()
            sage: g*s
            g1^2*g2^2*g3*g4
            sage: g*a
            g1^2*g2^2*g3^2

        TESTS::

            sage: h = G.gens()[3]
            sage: h in S
            False
        """
        gens_gap = [self(g).gap() for g in gens]
        gens_sage = [g.sage() for g in gens]
        G = self.gap().Subgroup(gens_gap)
        A = self._A.subgroup(gens_sage)
        return AbelianGroup_gap(A, G=G, ambient=self)
