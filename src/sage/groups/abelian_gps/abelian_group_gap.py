from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP
from sage.groups.group import AbelianGroup as AbelianGroupBase
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.arith.all import GCD, LCM



class AbelianGroupElement_gap(ElementLibGAP):
    r"""
    An element of an abelian group via libgap
    """
    def __init__(self, parent, x, check=False):
        """
        The Python constructor.

        See :class:`AbelianGroupElement_gap` for details.

        TESTS::
            
            sage: A = AbelianGroup([3,6])
            sage: G = A.gap()

        """
        if not isinstance(x, GapElement):
            A = parent._A
            x = A(x)
            # turn this into a gap element
            gens_gap = parent.gens()
            exp = x.exponents()
            x = gens_gap[0]**0
            for i in range(len(exp)):
                x *= gens_gap[i]**exp[i]
        if check:
            x in parent.gap()
        ElementLibGAP.__init__(self, parent, x)
        
    def __hash__(self):
        r"""
        """
        return hash(self.parent()) ^ hash(self._exponents)
    
    def exponents(self):
        r"""
        """
        P = self.parent()
        x = libgap.Factorization(P.gap(), self.gap())
        L = x.ExtRepOfObj()
        Lgens = L[::2]
        Lexpo = L[1::2]
        exp = []
        for k in range(len(P.gens())):
            if k in Lgens:
                exp.append(0)
            else:
                exp.append(Lexpo[k])
        return exp
    
    def order(self):
        r"""
        """
        if self.parent().is_finite():
            return self.gap().Order()
        else:
            # is there a way to do this in gap?
            P = self.parent()
            order = P.gens_orders()
            L = self.exponents()
            N = LCM([order[i]/GCD(order[i],L[i]) for i in range(len(order)) if L[i]!=0])
        return N
    
    def sage(self):
        r"""
        """
        P = self.parent()
        x = libgap.Factorization(P, self.gap())
        L = x.ExtRepOfObj()
        gens_sage = P._A.gens()
        e = P._A.identity()
        n = len(L)
        for k in range(0,n,2):
            e *= gens_sage[L[k]-1]**L[k+1]
        return e

class AbelianGroup_gap(GroupMixinLibGAP, ParentLibGAP, AbelianGroupBase):
    r"""
    """
    def __init__(self, A):
        r"""
        """
        AbelianGroupBase.__init__(self, category=A.category())
        ParentLibGAP.__init__(self, libgap(A))
        
        self._A = A
    
    Element = AbelianGroupElement_gap
           
    def _latex_(self):
        return self._A.latex()
    
    def _gap_init(self):
        return self._A._gap_init()
    
    def _repr_(self):
        r"""
        Return a string representation
        
        EXAMPLES::
        
            sage: A = AbelianGroup([2,6])
            sage: G = A.gap()
            sage: G._repr_()
        """
        s = self._A._repr_()
        s += " implemented with gap"
        return s

    def dual_group(self):
        return self._A.dual_group()
        
    def is_trivial(self):
        return 1 == self.order()
        
    def identity(self):
        return self(self.gap().Identity())
    
    @cached_method
    def elementary_divisors(self):
        ediv = self.gap().AbelianInvariants()
        from sage.matrix.constructor import diagonal_matrix
        ed = diagonal_matrix(ZZ, ediv).elementary_divisors()
        return tuple(d for d in ed if d!=1)
        
    @cached_method
    def exponent(self):
        r"""
        """
        return self.gap().Exponent()

    
    def gens_orders(self):
        r"""
        """
        return (g.order() for g in self.gens())
    
    def sage(self):
        r"""
        """
        return self._A
