r"""
"Named" Permutation groups (such as the symmetric group, S_n)

You can construct the following permutation groups:

-- SymmetricGroup, $S_n$ of order $n!$ (n can also be a list $X$ of distinct
                   positive integers, in which case it returns $S_X$)

-- AlternatingGroup, $A_n$ or order $n!/2$ (n can also be a list $X$
                   of distinct positive integers, in which case it returns
                   $A_X$)

-- DihedralGroup, $D_n$ of order $2n$

-- CyclicPermutationGroup, $C_n$ of order $n$

-- TransitiveGroup, $i^{th}$ transitive group of degree $n$
                      from the GAP tables of transitive groups (requires
                      the "optional" package database_gap)

-- MathieuGroup(degree), Mathieu group of degree 9, 10, 11, 12, 21, 22, 23, or 24.

-- KleinFourGroup, subgroup of $S_4$ of order $4$ which is not $C_2 \times C_2$

-- PGL(n,q), projective general linear group of $n\times n$ matrices over
             the finite field GF(q)

-- PSL(n,q), projective special linear group of $n\times n$ matrices over
             the finite field GF(q)

-- PSp(2n,q), projective symplectic linear group of $2n\times 2n$ matrices
              over the finite field GF(q)

-- PSU(n,q), projective special unitary group of $n\times n$ matrices having
             coefficients in the finite field $GF(q^2)$ that respect a
             fixed nondegenerate sesquilinear form, of determinant 1.

-- PGU(n,q), projective general unitary group of $n\times n$ matrices having
             coefficients in the finite field $GF(q^2)$ that respect a
             fixed nondegenerate sesquilinear form, modulo the centre.

-- SuzukiGroup(q), Suzuki group over GF(q), $^2 B_2(2^{2k+1}) = Sz(2^{2k+1})$.


AUTHOR:
    - David Joyner (2007-06): split from permgp.py (suggested by Nick Alexander)

REFERENCES:
    Cameron, P., Permutation Groups. New York: Cambridge University Press, 1999.
    Wielandt, H., Finite Permutation Groups. New York: Academic Press, 1964.
    Dixon, J. and Mortimer, B., Permutation Groups, Springer-Verlag, Berlin/New York, 1996.

NOTE:
    Though Suzuki groups are okay, Ree groups should *not* be wrapped as
    permutation groups - the onstruction is too slow - unless (for
    small values or the parameter) they are made using explicit generators.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                          David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import random
from types import ListType

import sage.structure.element as element
import sage.groups.group as group

from sage.rings.all      import RationalField, Integer
from sage.interfaces.all import gap, is_GapElement, is_ExpectElement
import sage.structure.coerce as coerce
from sage.rings.finite_field import FiniteField as GF
from sage.rings.arith import factor
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.misc.functional import is_even, log
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic

def permutation_list_to_degree(v):
    try:
        v = [Integer(z) for z in v]
    except TypeError:
        raise ValueError, "each entry of list must be an integer"
    if min(v) < 1:
        raise ValueError, "each element of list must be positive"
    return max(v), v


class SymmetricGroup(PermutationGroup_generic):
    """
    The full symmetric group of order $n!$, as a permutation group.
    (If n is a list of positive integers then it returns the
    symmetric group of the associated set.)

    INPUT:
       n -- a positive integer

    EXAMPLE:
        sage: G = SymmetricGroup(8)
        sage: G.order()
        40320
        sage: G
        Symmetric group of order 8! as a permutation group
        sage: G.degree()
        8
        sage: S8 = SymmetricGroup(8)
        sage: loads(dumps(S8)) == S8
        True
        sage: G = SymmetricGroup([1,2,4,5])
        sage: G
        Symmetric group of order 4! as a permutation group
        sage: G.set()
        [1, 2, 4, 5]
        sage: G = SymmetricGroup(4)
        sage: G
        Symmetric group of order 4! as a permutation group
        sage: G.set()
        [1, 2, 3, 4]
    """
    def __init__(self, n):
        if isinstance(n, (list, tuple)):
            self._deg, n = permutation_list_to_degree(n)
            PermutationGroup_generic.__init__(self, 'SymmetricGroup(%s)'%n, from_group = True)
            self._set = n
        else:
            try:
                self._deg = n = Integer(n)
                self._set = n
                if n < 1:
                    raise ValueError, "n (=%s) must be >= 1"%n
                PermutationGroup_generic.__init__(self, 'SymmetricGroup(%s)'%n, from_group = True)
            except TypeError, msg:
                raise ValueError, "%s\nn (=%s) must be an integer >= 1 or a list (but n has type %s)"%(msg, n,type(n))

    def _num_symbols(self):
        try:
            return self.__num_symbols
        except AttributeError:
            self.__num_symbols = len(self._set) if isinstance(self._set,ListType) else self.degree()
        return self.__num_symbols

    def _repr_(self):
        return "Symmetric group of order %s! as a permutation group"%self._num_symbols()

    def __str__(self):
        """
        EXAMPLES:
            sage: S = SymmetricGroup([2,3,7]); S
            Symmetric group of order 3! as a permutation group
            sage: str(S)
            'SymmetricGroup([2, 3, 7])'
            sage: S = SymmetricGroup(5); S
            Symmetric group of order 5! as a permutation group
            sage: str(S)
            'SymmetricGroup(5)'
        """
        if isinstance(self._set, ListType):
            x = self._set
        else:
            x = self.degree()
        return "SymmetricGroup(%s)"%x


    def set(self):
        if isinstance(self._set, list):
            return self._set
        else:
            return range(1, self._deg + 1)

class AlternatingGroup(PermutationGroup_generic):
    """
    The alternating group of order $n!/2$, as a permutation group.

    INPUT:
        n -- integer $n \geq 1$

    EXAMPLE:
        sage: G = AlternatingGroup(8)
        sage: G.order()
        20160
        sage: G
        Alternating group of order 8!/2 as a permutation group
        sage: loads(G.dumps()) == G
        True
        sage: G = AlternatingGroup([1,2,4,5])
        sage: G
        Alternating group of order 4!/2 as a permutation group
        sage: G.set()
        [1, 2, 4, 5]
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n):
        if isinstance(n, (list, tuple)):
            deg, n = permutation_list_to_degree(n)
            PermutationGroup_generic.__init__(self, 'AlternatingGroup(%s)'%n, from_group = True)
            self._set = n
        else:
            try:
                self._deg = n = Integer(n)
                self._set = n
                if n < 1:
                    raise ValueError, "n (=%s) must be >= 1"%n
                PermutationGroup_generic.__init__(self, 'AlternatingGroup(%s)'%n, from_group = True)
            except TypeError, msg:
                raise ValueError, "n (=%s) must be an integer >= 1 or a list"%n

    def _num_symbols(self):
        try:
            return self.__num_symbols
        except AttributeError:
            self.__num_symbols = len(self._set) if isinstance(self._set,ListType) else self.degree()
        return self.__num_symbols

    def _repr_(self):
        """
        EXAMPLES:
            sage: A = AlternatingGroup([2,3,7]); A
            Alternating group of order 3!/2 as a permutation group
        """
        z = self._num_symbols()
        return "Alternating group of order %s!/2 as a permutation group" % z

    def __str__(self):
        """
        EXAMPLES:
            sage: A = AlternatingGroup([2,3,7]); A
            Alternating group of order 3!/2 as a permutation group
            sage: str(A)
            'AlternatingGroup([2, 3, 7])'
        """
        if isinstance(self._set, ListType):
            x = self._set
        else:
            x = self.degree()
        return "AlternatingGroup(%s)" % x

    def set(self):
        if isinstance(self._set, list):
            return self._set
        else:
            return range(1, self._deg + 1)

class CyclicPermutationGroup(PermutationGroup_generic):
    """
    A cyclic group of order n, as a permutation group.

    INPUT:
        n -- a positive integer

    EXAMPLE:
        sage: G = CyclicPermutationGroup(8)
        sage: G.order()
        8
        sage: G
        Cyclic group of order 8 as a permutation group
        sage: loads(G.dumps()) == G
        True
        sage: C = CyclicPermutationGroup(10)
        sage: C.is_abelian()
        True
        sage: C = CyclicPermutationGroup(10)
        sage: C.as_AbelianGroup()
        Multiplicative Abelian Group isomorphic to C2 x C5

    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n):
        """
        """
        n = Integer(n)
        if n < 1:
            raise ValueError, "n (=%s) must be >= 1"%n
        gens = tuple(range(1, n+1))
        PermutationGroup_generic.__init__(self, [gens], n)

    def _repr_(self):
        return "Cyclic group of order %s as a permutation group"%self.order()

    def is_commutative(self):
        return True

    def is_abelian(self):
        return True

    def as_AbelianGroup(self):
	"""
	Returns the corresponding Abelian Group instance.

	EXAMPLES:

	"""
	n = self.order()
        a = list(factor(n))
        invs = [x[0]**x[1] for x in a]
        G = AbelianGroup(len(a),invs)
	return G

class KleinFourGroup(PermutationGroup_generic):
    r"""
    The Klein 4 Group, which has order $4$ and exponent $2$, viewed
    as a subgroup of $S_4$.

    OUTPUT:
        -- the Klein 4 group of order 4, as a permutation group of degree 4.

    EXAMPLES:
        sage: G = KleinFourGroup(); G
        The Klein 4 group of order 4, as a permutation group
        sage: list(G)
        [(), (3,4), (1,2), (1,2)(3,4)]

    AUTHOR:
        -- Bobby Moretti (2006-10)
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self):
        gens = ((1,2),(3,4))
        PermutationGroup_generic.__init__(self, gens, from_group=True)

    def _repr_(self):
        return 'The Klein 4 group of order 4, as a permutation group'


class DihedralGroup(PermutationGroup_generic):
    """
    The Dihedral group of order $2n$ for any integer $n\geq 1$.

    INPUT:
        n -- a positive integer

    OUTPUT:
        -- the dihedral group of order 2*n, as a permutation group

    EXAMPLE:
        sage: DihedralGroup(1)
        Dihedral group of order 2 as a permutation group

        sage: DihedralGroup(2)
        Dihedral group of order 4 as a permutation group
        sage: DihedralGroup(2).gens()
        ((1,2), (3,4))

        sage: DihedralGroup(5).gens()
        ((1,2,3,4,5), (1,5)(2,4))
        sage: list(DihedralGroup(5))
        [(), (2,5)(3,4), (1,2)(3,5), (1,2,3,4,5), (1,3)(4,5), (1,3,5,2,4), (1,4)(2,3), (1,4,2,5,3), (1,5,4,3,2), (1,5)(2,4)]

        sage: G = DihedralGroup(6)
        sage: G.order()
        12
        sage: G = DihedralGroup(5)
        sage: G.order()
        10
        sage: G
        Dihedral group of order 10 as a permutation group
        sage: loads(G.dumps()) == G
        True
        sage: G.gens()
        ((1,2,3,4,5), (1,5)(2,4))

        sage: DihedralGroup(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n):
        n = Integer(n)
        if n <= 0:
            raise ValueError, "n must be positive"

        # the first generator generates the cyclic subgroup of D_n, <(1...n)> in
        # cycle notation
        gen0 = range(1,n+1)

        if n < 1:
            raise ValueError, "n (=%s) must be >= 1"%n

        # D_1 is a subgroup of S_2, we need the cyclic group of order 2
        if n == 1:
            gens = CyclicPermutationGroup(2).gens()
        elif n == 2:
            gens = ((1,2),(3,4))
        else:
            #gens = tuple(gen0.append([(i, n+1-i) for i in range(1, n//2 + 1)]))
            gen1 = [(i, n-i+1) for i in range(1, n//2 +1)]
            gens = tuple([tuple(gen0),tuple(gen1)])
        # send this off to the parent's class __init__()
        PermutationGroup_generic.__init__(self, gens, from_group = True)

    def _repr_(self):
        return "Dihedral group of order %s as a permutation group"%self.order()

class MathieuGroup(PermutationGroup_generic):
    """
    The Mathieu group of degree $n$.

    INPUT:
        n -- a positive integer in  {9, 10, 11, 12, 21, 22, 23, 24}.

    OUTPUT:
        -- the Mathieu group of degree n, as a permutation group

    EXAMPLE:
        sage: G = MathieuGroup(12)
        sage: G
        Mathieu group of degree 12 and order 95040 as a permutation group
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n):
        n = Integer(n)
        self._n = n
        if not(n in [9, 10, 11, 12, 21, 22, 23, 24]):
            raise ValueError,"argument must belong to {9, 10, 11, 12, 21, 22, 23, 24}."
        id = 'MathieuGroup(%s)'%n
        PermutationGroup_generic.__init__(self, id, from_group=True, check=False)

    def _repr_(self):
        return "Mathieu group of degree %s and order %s as a permutation group"%(self._n,self.order())

class TransitiveGroup(PermutationGroup_generic):
    """
    The transitive group from the GAP tables of transitive groups.
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, d, n):
        """
        INPUT:
            d -- positive integer; the degree
            n -- positive integer; the number

        OUTPUT:
            the n-th transitive group of degree d

        EXAMPLE:
            sage: G = TransitiveGroup(1,1); G
            Transitive group number 1 of degree 1
            sage: G = TransitiveGroup(5, 2); G         # requires optional database_gap
            Transitive group number 2 of degree 5
            sage: G.gens()                             # requires optional database_gap
            ((1,2,3,4,5), (1,4)(2,3))

            sage: loads(G.dumps()) == G                # requires optional database_gap
            True
        """
        if d == 1:
            id = 'Group([()])'
        else:
            id = 'TransitiveGroup(%s,%s)'%(d,n)
        try:
            PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        except RuntimeError:
            from sage.misc.misc import verbose
            verbose("Warning: Computing with TransitiveGroups requires the optional database_gap package. Please install it.", level=0)

        self._d = d
        self._n = n

    def _repr_(self):
        return "Transitive group number %s of degree %s"%(self._n, self._d)

class PGL(PermutationGroup_generic):
    """
    The projective general linear groups over GF(q).
    """
    def __init__(self, n, q, name='a'):
        """
        INPUT:
            n -- positive integer; the degree
            q -- prime power; the size of the ground field
            name -- (default: 'a') variable name of indeterminate of finite field GF(q)

        OUTPUT:
            PGL(n,q)

        EXAMPLE:
            sage: G = PGL(2,3); G
            Permutation Group with generators [(3,4), (1,2,4)]
            sage: print G
            The projective general linear group of degree 2 over Finite Field of size 3
            sage: G.base_ring()
            Finite Field of size 3
            sage: G.order()
            24


            sage: G = PGL(2, 9, 'b'); G
            Permutation Group with generators [(3,10,9,8,4,7,6,5), (1,2,4)(5,6,8)(7,9,10)]
            sage: G.base_ring()
            Finite Field in b of size 3^2
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
        if n == 1:
            id = 'Group([()])'
        else:
            id = 'PGL(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._n = n

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The projective general linear group of degree %s over %s"%(self._n, self.base_ring())

class PSL(PermutationGroup_generic):
    """
    The projective special linear groups over GF(q).
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n, q, name='a'):
        """
        INPUT:
            n -- positive integer; the degree
            q -- prime power; the size of the ground field
            name -- (default: 'a') variable name of indeterminate of finite field GF(q)

        OUTPUT:
            PSL(n,q)

        EXAMPLE:
            sage: G = PSL(2,3); G
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
            sage: G.order()
            12
            sage: G.base_ring()
            Finite Field of size 3
            sage: print G
            The projective special linear group of degree 2 over Finite Field of size 3

        We create two groups over nontrivial finite fields:
            sage: G = PSL(2, 4, 'b'); G
            Permutation Group with generators [(3,4,5), (1,2,3)]
            sage: G.base_ring()
            Finite Field in b of size 2^2
            sage: G = PSL(2, 8); G
            Permutation Group with generators [(3,8,6,4,9,7,5), (1,2,3)(4,7,5)(6,9,8)]
            sage: G.base_ring()
            Finite Field in a of size 2^3

        """
        if n == 1:
            id = 'Group([()])'
        else:
            id = 'PSL(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._n = n

    def matrix_degree(self):
        return self._n

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The projective special linear group of degree %s over %s"%(self._n, self.base_ring())

    def ramification_module_decomposition_hurwitz_curve(self):
        """
        Helps compute the decomposition of the ramification module
        for the Hurwitz curves X (over CC say) with automorphism group
        G = PSL(2,q), q a "Hurwitz prime" (ie, p is $\pm 1 \pmod 7$).
        Using this computation and Borne's formula helps determine the
        G-module structure of the RR spaces of equivariant
        divisors can be determined explicitly.

        The output is a list of integer multiplicities: [m1,...,mn],
        where n is the number of conj classes of G=PSL(2,p) and mi is the
        multiplicity of pi_i in the ramification module of a
        Hurwitz curve with automorphism group G.
        Here IrrRepns(G) = [pi_1,...,pi_n] (in the order listed in the
        output of self.character_table()).

        REFERENCE: David Joyner, Amy Ksir, Roger Vogeler,
                   "Group representations on Riemann-Roch spaces of some
                   Hurwitz curves," preprint, 2006.

        EXAMPLES:
            sage: G = PSL(2,13)
            sage: G.ramification_module_decomposition_hurwitz_curve() #random
            [0, 7, 7, 12, 12, 12, 13, 15, 14]

        This means, for example, that the trivial representation does not
        occur in the ramification module of a Hurwitz curve with automorphism
        group PSL(2,13), since the trivial representation is listed first
        and that entry has multiplicity 0. The "randomness" is due to the
        fact that GAP randomly orders the conjugacy classes of the same order
        in the list of all conjugacy classes. Similarly, there is some
        randomness to the ordering of the characters.

        If you try to use this function on a group PSL(2,q) where q is
        not a (smallish) "Hurwitz prime", an error mesage will be printed.
        """
        if self.matrix_degree()!=2:
            return ValueError, "Degree must be 2."
        F = self.base_ring()
        q = F.order()
        from sage.misc.misc import SAGE_EXTCODE
        gapcode = SAGE_EXTCODE + '/gap/joyner/hurwitz_crv_rr_sp.gap'
        gap.eval('Read("'+gapcode+'")')
        mults = gap.eval("ram_module_hurwitz("+str(q)+")")
        return eval(mults)

    def ramification_module_decomposition_modular_curve(self):
        """
        Helps compute the decomposition of the ramification module
        for the modular curve X(p) (over CC say) with automorphism group G = PSL(2,q),
        q a prime > 5. Using this computation and Borne's formula helps determine the
        G-module structure of the RR spaces of equivariant
        divisors can be determined explicitly.

        The output is a list of integer multiplicities: [m1,...,mn],
        where n is the number of conj classes of G=PSL(2,p) and mi is the
        multiplicity of pi_i in the ramification module of a
        modular curve with automorphism group G.
        Here IrrRepns(G) = [pi_1,...,pi_n] (in the order listed in the
        output of self.character_table()).

        REFERENCE: D. Joyner and A. Ksir, 'Modular representations
                   on some Riemann-Roch spaces of modular curves
                   $X(N)$', Computational Aspects of Algebraic Curves,
                   (Editor: T. Shaska) Lecture Notes in Computing, WorldScientific,
                   2005.)

        EXAMPLES:
            sage: G = PSL(2,7)
            sage: G.ramification_module_decomposition_modular_curve() ## random
            [0, 4, 3, 6, 7, 8]

        This means, for example, that the trivial representation does not
        occur in the ramification module of X(7), since the trivial representation
        is listed first and that entry has multiplicity 0. The "randomness" is due to the
        fact that GAP randomly orders the conjugacy classes of the same order
        in the list of all conjugacy classes. Similarly, there is some
        randomness to the ordering of the characters.
        """
        if self.matrix_degree()!=2:
            return ValueError, "Degree must be 2."
        F = self.base_ring()
        q = F.order()
        from sage.misc.misc import SAGE_EXTCODE
        gapcode = SAGE_EXTCODE + '/gap/joyner/modular_crv_rr_sp.gap'
        gap.eval('Read("'+gapcode+'")')
        mults = gap.eval("ram_module_X("+str(q)+")")
        return eval(mults)

class PSp(PermutationGroup_generic):
    """
    The projective symplectic linear groups over GF(q).
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n, q, name='a'):
        """
        INPUT:
            n -- positive integer; the degree
            q -- prime power; the size of the ground field
            name -- (default: 'a') variable name of indeterminate of finite field GF(q)

        OUTPUT:
            PSp(n,q)

        EXAMPLE:
            sage: G = PSp(2,3); G
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
            sage: G.order()
            12
            sage: G = PSp(4,3); G
            Permutation Group with generators [(3,4)(6,7)(9,10)(12,13)(17,20)(18,21)(19,22)(23,32)(24,33)(25,34)(26,38)(27,39)(28,40)(29,35)(30,36)(31,37), (1,5,14,17,27,22,19,36,3)(2,6,32)(4,7,23,20,37,13,16,26,40)(8,24,29,30,39,10,33,11,34)(9,15,35)(12,25,38)(21,28,31)]
            sage: G.order()
            25920
            sage: print G
            The projective symplectic linear group of degree 4 over Finite Field of size 3
            sage: G.base_ring()
            Finite Field of size 3

            sage: G = PSp(2, 8, name='alpha'); G
            Permutation Group with generators [(3,8,6,4,9,7,5), (1,2,3)(4,7,5)(6,9,8)]
            sage: G.base_ring()
            Finite Field in alpha of size 2^3
        """
        if n%2 == 1:
            raise TypeError, "The degree n must be even"
        else:
            id = 'PSp(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._n = n

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The projective symplectic linear group of degree %s over %s"%(self._n, self.base_ring())

PSP = PSp

class PSU(PermutationGroup_generic):
    """
    The projective special unitary groups over GF(q).

    INPUT:
        n -- positive integer; the degree
        q -- prime power; the size of the ground field
        name -- (default: 'a') variable name of indeterminate of finite field GF(q)
    OUTPUT:
        PSU(n,q)

    EXAMPLE:
        sage: PSU(2,3)
        The projective special unitary group of degree 2 over Finite Field of size 3

        sage: G = PSU(2, 8, name='alpha'); G
        The projective special unitary group of degree 2 over Finite Field in alpha of size 2^3
        sage: G.base_ring()
        Finite Field in alpha of size 2^3
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n, q, name='a'):
        id = 'PSU(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._field_of_definition = GF(q**2, name)
        self._n = n

    def field_of_definition(self):
        return self._field_of_definition

    def base_ring(self):
        return self._base_ring

    def _repr_(self):
        return "The projective special unitary group of degree %s over %s"%(self._n, self.base_ring())

class PGU(PermutationGroup_generic):
    """
    The projective general unitary groups over GF(q).

    INPUT:
        n -- positive integer; the degree
        q -- prime power; the size of the ground field
        name -- (default: 'a') variable name of indeterminate of finite field GF(q)

    OUTPUT:
        PGU(n,q)

    EXAMPLE:
        sage: PGU(2,3)
        The projective general unitary group of degree 2 over Finite Field of size 3

        sage: G = PGU(2, 8, name='alpha'); G
        The projective general unitary group of degree 2 over Finite Field in alpha of size 2^3
        sage: G.base_ring()
        Finite Field in alpha of size 2^3
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, n, q, name='a'):
        id = 'PGU(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q, name=name)
        self._field_of_definition = GF(q**2, name)
        self._n = n

    def field_of_definition(self):
        return self._field_of_definition

    def base_ring(self):
        return self._base_ring

    def _repr_(self):
        return "The projective general unitary group of degree %s over %s"%(self._n, self.base_ring())


class SuzukiGroup(PermutationGroup_generic):
    r"""
    The Suzuki group over GF(q), $^2 B_2(2^{2k+1}) = Sz(2^{2k+1})$. A wrapper for the GAP function SuzukiGroup.

    INPUT:
        q -- 2^n, an odd power of 2; the size of the ground
             field. (Strictly speaking, n should be greater than 1, or
             else this group os not simple.)
        name -- (default: 'a') variable name of indeterminate of
                finite field GF(q)

    OUTPUT:
        A Suzuki group.

    EXAMPLE:
        sage: SuzukiGroup(8)
        Permutation Group with generators [(1,28,10,44)(3,50,11,42)(4,43,53,64)(5,9,39,52)(6,36,63,13)(7,51,60,57)(8,33,37,16)(12,24,55,29)(14,30,48,47)(15,19,61,54)(17,59,22,62)(18,23,34,31)(20,38,49,25)(21,26,45,58)(27,32,41,65)(35,46,40,56), (1,2)(3,10)(4,42)(5,18)(6,50)(7,26)(8,58)(9,34)(12,28)(13,45)(14,44)(15,23)(16,31)(17,21)(19,39)(20,38)(22,25)(24,61)(27,60)(29,65)(30,55)(32,33)(35,52)(36,49)(37,59)(40,54)(41,62)(43,53)(46,48)(47,56)(51,63)(57,64)]
        sage: print SuzukiGroup(8)
        The Suzuki group over Finite Field in a of size 2^3

        sage: G = SuzukiGroup(32, name='alpha')
        sage: G.order()
        32537600
        sage: G.order().factor()
        2^10 * 5^2 * 31 * 41
        sage: G.base_ring()
        Finite Field in alpha of size 2^5

    REFERENCES:
        http://en.wikipedia.org/wiki/Group_of_Lie_type\#Suzuki-Ree_groups
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
    def __init__(self, q, name='a'):
        q = Integer(q)
        from sage.rings.arith import valuation
        t = valuation(q, 2)
        if 2**t != q or is_even(t):
	    raise ValueError,"The ground field size %s must be an odd power of 2."%q
        id = 'SuzukiGroup(IsPermGroup,%s)'%q
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q, name=name)

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The Suzuki group over %s"%self.base_ring()







