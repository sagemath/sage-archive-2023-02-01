r"""
The set `\mathbb{P}^1(\mathbb{Q})` of cusps

EXAMPLES::

    sage: Cusps
    Set P^1(QQ) of all cusps

::

    sage: Cusp(oo)
    Infinity
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#*****************************************************************************

from sage.rings.all import is_Infinite, Rational, Integer, ZZ, QQ
from sage.rings.integer_ring import IntegerRing
from sage.structure.parent_base import ParentWithBase
from sage.structure.element import Element, is_InfinityElement
from sage.modular.modsym.p1list import lift_to_sl2z_llong

class Cusps_class(ParentWithBase):
    """
    The set of cusps.

    EXAMPLES::

        sage: C = Cusps; C
        Set P^1(QQ) of all cusps
        sage: loads(C.dumps()) == C
        True
    """
    def __init__(self):
        r"""
        The set of cusps, i.e. `\mathbb{P}^1(\mathbb{Q})`.

        EXAMPLES::

            sage: C = sage.modular.cusps.Cusps_class() ; C
            Set P^1(QQ) of all cusps
            sage: Cusps == C
            True
        """
        ParentWithBase.__init__(self, self)

    def __cmp__(self, right):
        """
        Return equality only if right is the set of cusps.

        EXAMPLES::

            sage: Cusps == Cusps
            True
            sage: Cusps == QQ
            False
        """
        return cmp(type(self), type(right))

    def _repr_(self):
        """
        String representation of the set of cusps.

        EXAMPLES::

            sage: Cusps
            Set P^1(QQ) of all cusps
            sage: Cusps._repr_()
            'Set P^1(QQ) of all cusps'
            sage: Cusps.rename('CUSPS'); Cusps
            CUSPS
            sage: Cusps.rename(); Cusps
            Set P^1(QQ) of all cusps
            sage: Cusps
            Set P^1(QQ) of all cusps
        """
        return "Set P^1(QQ) of all cusps"

    def _latex_(self):
        """
        Return latex representation of self.

        EXAMPLES::

            sage: latex(Cusps)
            \mathbf{P}^1(\mathbf{Q})
            sage: latex(Cusps) == Cusps._latex_()
            True
        """
        return "\\mathbf{P}^1(\\mathbf{Q})"

    def __call__(self, x):
        """
        Coerce x into the set of cusps.

        EXAMPLES::

            sage: a = Cusps(-4/5); a
            -4/5
            sage: Cusps(a) is a
            False
            sage: Cusps(1.5)
            3/2
            sage: Cusps(oo)
            Infinity
            sage: Cusps(I)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert I to a Cusp
        """
        return Cusp(x, parent=self)

    def _coerce_impl(self, x):
        """
        Canonical coercion of x into the set of cusps.

        EXAMPLES::

            sage: Cusps._coerce_(7/13)
            7/13
            sage: Cusps._coerce_(GF(7)(3))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self
            sage: Cusps(GF(7)(3))
            3
            sage: Cusps._coerce_impl(GF(7)(3))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self
        """
        if is_Infinite(x):
            return Cusp(x, parent=self)
        else:
            return self._coerce_try(x, QQ)

Cusps = Cusps_class()


class Cusp(Element):
    """
    A cusp.

    A cusp is either a rational number or infinity, i.e., an element of
    the projective line over Q. A Cusp is stored as a pair (a,b), where
    gcd(a,b)=1 and a,b are of type Integer.

    EXAMPLES::

        sage: a = Cusp(2/3); b = Cusp(oo)
        sage: a.parent()
        Set P^1(QQ) of all cusps
        sage: a.parent() is b.parent()
        True
    """

    def __init__(self, a, b=None, parent=None, check=True):
        r"""
        Create the cusp a/b in `\mathbb{P}^1(\mathbb{Q})`, where if b=0
        this is the cusp at infinity.

        When present, b must either be Infinity or coercible to an
        Integer.

        EXAMPLES::

            sage: Cusp(2,3)
            2/3
            sage: Cusp(3,6)
            1/2
            sage: Cusp(1,0)
            Infinity
            sage: Cusp(infinity)
            Infinity
            sage: Cusp(5)
            5
            sage: Cusp(1/2)
            1/2
            sage: Cusp(1.5)
            3/2
            sage: Cusp(int(7))
            7
            sage: Cusp(1, 2, check=False)
            1/2
            sage: Cusp('sage', 2.5, check=False)          # don't do this!
            sage/2.50000000000000

        ::

            sage: I**2
            -1
            sage: Cusp(I)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert I to a Cusp

        ::

            sage: a = Cusp(2,3)
            sage: loads(a.dumps()) == a
            True

        ::

            sage: Cusp(1/3,0)
            Infinity
            sage: Cusp((1,0))
            Infinity

        TESTS::

            sage: Cusp("1/3", 5)
            1/15
            sage: Cusp(Cusp(3/5), 7)
            3/35
            sage: Cusp(5/3, 0)
            Infinity
            sage: Cusp(3,oo)
            0
            sage: Cusp((7,3), 5)
            7/15
            sage: Cusp(int(5), 7)
            5/7

        ::

            sage: Cusp(0,0)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert (0, 0) to a Cusp

        ::

            sage: Cusp(oo,oo)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert (+Infinity, +Infinity) to a Cusp

        ::

            sage: Cusp(Cusp(oo),oo)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert (Infinity, +Infinity) to a Cusp
        """
        if parent is None:
            parent = Cusps
        Element.__init__(self, parent)

        if not check:
            self.__a = a; self.__b = b
            return

        if b is None:
            if isinstance(a, Integer):
                self.__a = a
                self.__b = ZZ(1)
            elif isinstance(a, Rational):
                self.__a = a.numer()
                self.__b = a.denom()
            elif is_InfinityElement(a):
                self.__a = ZZ(1)
                self.__b = ZZ(0)
            elif isinstance(a, Cusp):
                self.__a = a.__a
                self.__b = a.__b
            elif isinstance(a, (int, long)):
                self.__a = ZZ(a)
                self.__b = ZZ(1)
            elif isinstance(a, (tuple, list)):
                if len(a) != 2:
                    raise TypeError, "Unable to convert %s to a Cusp"%a
                if ZZ(a[1]) == 0:
                    self.__a = ZZ(1)
                    self.__b = ZZ(0)
                    return
                try:
                    r = QQ((a[0], a[1]))
                    self.__a = r.numer()
                    self.__b = r.denom()
                except (ValueError, TypeError):
                    raise TypeError, "Unable to convert %s to a Cusp"%a
            else:
                try:
                    r = QQ(a)
                    self.__a = r.numer()
                    self.__b = r.denom()
                except (ValueError, TypeError):
                    raise TypeError, "Unable to convert %s to a Cusp"%a
            return

        if is_InfinityElement(b):
            if is_InfinityElement(a) or (isinstance(a, Cusp) and a.is_infinity()):
                raise TypeError, "Unable to convert (%s, %s) to a Cusp"%(a, b)
            self.__a = ZZ(0)
            self.__b = ZZ(1)
            return
        elif not b:
            if not a:
                raise TypeError, "Unable to convert (%s, %s) to a Cusp"%(a, b)
            self.__a = ZZ(1)
            self.__b = ZZ(0)
            return

        if isinstance(a, Integer) or isinstance(a, Rational):
            r = a / ZZ(b)
        elif is_InfinityElement(a):
            self.__a = ZZ(1)
            self.__b = ZZ(0)
            return
        elif isinstance(a, Cusp):
            if a.__b:
                r = a.__a / (a.__b * b)
            else:
                self.__a = ZZ(1)
                self.__b = ZZ(0)
                return
        elif isinstance(a, (int, long)):
            r = ZZ(a) / b
        elif isinstance(a, (tuple, list)):
            if len(a) != 2:
                raise TypeError, "Unable to convert (%s, %s) to a Cusp"%(a, b)
            r = ZZ(a[0]) / (ZZ(a[1]) * b)
        else:
            try:
                r = QQ(a) / b
            except (ValueError, TypeError):
                raise TypeError, "Unable to convert (%s, %s) to a Cusp"%(a, b)

        self.__a = r.numer()
        self.__b = r.denom()


    def __hash__(self):
        """
        EXAMPLES:
            sage: hash(Cusp(1/3))
            1298787075             # 32-bit
            3713081631933328131    # 64-bit
            sage: hash(Cusp(oo))
            1302034650             # 32-bit
            3713081631936575706    # 64-bit
        """
        return hash((self.__a, self.__b))

    def __cmp__(self, right):
        """
        Compare the cusps self and right. Comparison is as for rational
        numbers, except with the cusp oo greater than everything but
        itself.

        The ordering in comparison is only really meaningful for infinity
        or elements that coerce to the rationals.

        EXAMPLES::

            sage: Cusp(2/3) == Cusp(oo)
            False

        ::

            sage: Cusp(2/3) < Cusp(oo)
            True

        ::

            sage: Cusp(2/3)> Cusp(oo)
            False

        ::

            sage: Cusp(2/3) > Cusp(5/2)
            False

        ::

            sage: Cusp(2/3) < Cusp(5/2)
            True

        ::

            sage: Cusp(2/3) == Cusp(5/2)
            False

        ::

            sage: Cusp(oo) == Cusp(oo)
            True

        ::

            sage: 19/3 < Cusp(oo)
            True

        ::

            sage: Cusp(oo) < 19/3
            False

        ::

            sage: Cusp(2/3) < Cusp(11/7)
            True

        ::

            sage: Cusp(11/7) < Cusp(2/3)
            False

        ::

            sage: 2 < Cusp(3)
            True
        """
        if not self.__b:
            # self is oo, which is bigger than everything but oo.
            if not right.__b:
                return 0
            else:
                return 1
        elif not right.__b:
            if not self.__b:
                return 0
            else:
                return -1
        return cmp(self._rational_(), right._rational_())

    def is_infinity(self):
        """
        Returns True if this is the cusp infinity.

        EXAMPLES::

            sage: Cusp(3/5).is_infinity()
            False
            sage: Cusp(1,0).is_infinity()
            True
            sage: Cusp(0,1).is_infinity()
            False
        """
        return not self.__b

    def numerator(self):
        """
        Return the numerator of the cusp a/b.

        EXAMPLES::

            sage: x=Cusp(6,9); x
            2/3
            sage: x.numerator()
            2
            sage: Cusp(oo).numerator()
            1
            sage: Cusp(-5/10).numerator()
            -1
        """
        return self.__a

    def denominator(self):
        """
        Return the denominator of the cusp a/b.

        EXAMPLES::

            sage: x=Cusp(6,9); x
            2/3
            sage: x.denominator()
            3
            sage: Cusp(oo).denominator()
            0
            sage: Cusp(-5/10).denominator()
            2
        """
        return self.__b

    def _rational_(self):
        """
        Coerce to a rational number.

        EXAMPLES::

            sage: QQ(Cusp(oo))
            Traceback (most recent call last):
            ...
            TypeError: cusp Infinity is not a rational number
            sage: QQ(Cusp(-3,7))
            -3/7
            sage: Cusp(11,2)._rational_()
            11/2
        """
        try:
            return self.__rational
        except AttributeError:
            pass

        if not self.__b:
            raise TypeError, "cusp %s is not a rational number"%self
        self.__rational = self.__a / self.__b
        return self.__rational

    def _integer_(self, ZZ=None):
        """
        Coerce to an integer.

        EXAMPLES::

            sage: ZZ(Cusp(-19))
            -19
            sage: Cusp(4,2)._integer_()
            2

        ::

            sage: ZZ(Cusp(oo))
            Traceback (most recent call last):
            ...
            TypeError: cusp Infinity is not an integer
            sage: ZZ(Cusp(-3,7))
            Traceback (most recent call last):
            ...
            TypeError: cusp -3/7 is not an integer
        """
        if self.__b != 1:
            raise TypeError, "cusp %s is not an integer"%self
        return self.__a

    def _repr_(self):
        """
        String representation of this cusp.

        EXAMPLES::

            sage: a = Cusp(2/3); a
            2/3
            sage: a._repr_()
            '2/3'
            sage: a.rename('2/3(cusp)'); a
            2/3(cusp)
        """
        if self.__b.is_zero():
            return "Infinity"
        if self.__b != 1:
            return "%s/%s"%(self.__a,self.__b)
        else:
            return str(self.__a)

    def _latex_(self):
        r"""
        Latex representation of this cusp.

        EXAMPLES::

            sage: latex(Cusp(-2/7))
            \frac{-2}{7}
            sage: latex(Cusp(oo))
            \infty
            sage: latex(Cusp(oo)) == Cusp(oo)._latex_()
            True
        """
        if self.__b.is_zero():
            return "\\infty"
        if self.__b != 1:
            return "\\frac{%s}{%s}"%(self.__a,self.__b)
        else:
            return str(self.__a)

    def __neg__(self):
        """
        The negative of this cusp.

        EXAMPLES::

            sage: -Cusp(2/7)
            -2/7
            sage: -Cusp(oo)
            Infinity
        """
        return Cusp(-self.__a, self.__b)

    def is_gamma0_equiv(self, other, N, transformation = False):
        r"""
        Return whether self and other are equivalent modulo the action of
        `\Gamma_0(N)` via linear fractional transformations.

        INPUT:


        -  ``other`` - Cusp

        -  ``N`` - an integer (specifies the group
           Gamma_0(N))

        -  ``transformation`` - bool (default: False), if True,
           also return upper left entry of a matrix in Gamma_0(N) that sends
           self to other.


        OUTPUT:


        -  ``bool`` - True if self and other are equivalent

        -  ``integer`` - returned only if transformation is
           True


        EXAMPLES::

            sage: x = Cusp(2,3)
            sage: y = Cusp(4,5)
            sage: x.is_gamma0_equiv(y, 2)
            True
            sage: x.is_gamma0_equiv(y, 2, True)
            (True, 1)
            sage: x.is_gamma0_equiv(y, 3)
            False
            sage: x.is_gamma0_equiv(y, 3, True)
            (False, None)
            sage: Cusp(1,0)
            Infinity
            sage: z = Cusp(1,0)
            sage: x.is_gamma0_equiv(z, 3, True)
            (True, 2)

        ALGORITHM: See Proposition 2.2.3 of Cremona's book "Algorithms for
        Modular Elliptic Curves", or Prop 2.27 of Stein's Ph.D. thesis.
        """
        if not isinstance(other, Cusp):
            other = Cusp(other)
        N = ZZ(N)
        u1 = self.__a
        v1 = self.__b
        u2 = other.__a
        v2 = other.__b
        if (u1,v1) != (ZZ(0),ZZ(1)):
            if v1 in [ZZ(0),ZZ(1)]:
                s1 = ZZ(1)
            else:
                s1 = u1.inverse_mod(v1)
        else:
            s1 = 0
        if (u2,v2) != (ZZ(0),ZZ(1)):
            if v2 in [ZZ(0),ZZ(1)]:
                s2 = 1
            else:
                s2 = u2.inverse_mod(v2)
        else:
            s2 = 0
        g = (v1*v2).gcd(N)
        a = s1*v2 - s2*v1
        if a%g != ZZ(0):
            if transformation:
                return False, None
            else:
                return False

        if not transformation:
            return True

        # Now we know the cusps are equivalent.  Use the proof of Prop 2.2.3
        # of Cremona to find a matrix in Gamma_0(N) relating them.
        dum, s2, r2 = u2.xgcd(-v2)
        assert dum.is_one()
        dum, s1, r1 = u1.xgcd(-v1)
        assert dum.is_one()
        a = s1*v2 - s2*v1
        assert (a%g).is_zero()
        # solve x*v1*v2 + a = 0 (mod N).
        d,x0,y0 = (v1*v2).xgcd(N)          # x0*v1*v2 + y0*N = d = g.
        # so x0*v1*v2 - g = 0 (mod N)
        x = -x0 * (a/g)
        # now  x*v1*v2 + a = 0 (mod N)
        s1p = s1+x*v1
        return (True, (u2*s1p-r2*v1)%N)

    def is_gamma1_equiv(self, other, N):
        """
        Return whether self and other are equivalent modulo the action of
        Gamma_1(N) via linear fractional transformations.

        INPUT:


        -  ``other`` - Cusp

        -  ``N`` - an integer (specifies the group
           Gamma_1(N))


        OUTPUT:


        -  ``bool`` - True if self and other are equivalent

        -  ``int`` - 0, 1 or -1, gives further information
           about the equivalence: If the two cusps are u1/v1 and u2/v2, then
           they are equivalent if and only if v1 = v2 (mod N) and u1 = u2 (mod
           gcd(v1,N)) or v1 = -v2 (mod N) and u1 = -u2 (mod gcd(v1,N)) The
           sign is +1 for the first and -1 for the second. If the two cusps
           are not equivalent then 0 is returned.


        EXAMPLES::

            sage: x = Cusp(2,3)
            sage: y = Cusp(4,5)
            sage: x.is_gamma1_equiv(y,2)
            (True, 1)
            sage: x.is_gamma1_equiv(y,3)
            (False, 0)
            sage: z = Cusp(QQ(x) + 10)
            sage: x.is_gamma1_equiv(z,10)
            (True, 1)
            sage: z = Cusp(1,0)
            sage: x.is_gamma1_equiv(z, 3)
            (True, -1)
            sage: Cusp(0).is_gamma1_equiv(oo, 1)
            (True, 1)
            sage: Cusp(0).is_gamma1_equiv(oo, 3)
            (False, 0)
        """
        if not isinstance(other, Cusp):
            other = Cusp(other)
        N = ZZ(N)
        u1 = self.__a
        v1 = self.__b
        u2 = other.__a
        v2 = other.__b
        g = v1.gcd(N)
        if ((v2 - v1) % N == 0 and (u2 - u1)%g== 0):
            return True, 1
        elif ((v2 + v1) % N == 0 and (u2 + u1)%g== 0):
            return True, -1
        return False, 0

    def is_gamma_h_equiv(self, other, G):
        """
        Return a pair (b, t), where b is True or False as self and other
        are equivalent under the action of G, and t is 1 or -1, as
        described below.

        Two cusps `u1/v1` and `u2/v2` are equivalent modulo
        Gamma_H(N) if and only if `v1 =  h*v2 (\mathrm{mod} N)` and
        `u1 =  h^{(-1)}*u2 (\mathrm{mod} gcd(v1,N))` or
        `v1 = -h*v2 (mod N)` and
        `u1 = -h^{(-1)}*u2 (\mathrm{mod} gcd(v1,N))` for some
        `h \in H`. Then t is 1 or -1 as c and c' fall into the
        first or second case, respectively.

        INPUT:


        -  ``other`` - Cusp

        -  ``G`` - a congruence subgroup Gamma_H(N)


        OUTPUT:


        -  ``bool`` - True if self and other are equivalent

        -  ``int`` - -1, 0, 1; extra info


        EXAMPLES::

            sage: x = Cusp(2,3)
            sage: y = Cusp(4,5)
            sage: x.is_gamma_h_equiv(y,GammaH(13,[2]))
            (True, 1)
            sage: x.is_gamma_h_equiv(y,GammaH(13,[5]))
            (False, 0)
            sage: x.is_gamma_h_equiv(y,GammaH(5,[]))
            (False, 0)
            sage: x.is_gamma_h_equiv(y,GammaH(23,[4]))
            (True, -1)

        Enumerating the cusps for a space of modular symbols uses this
        function.

        ::

            sage: G = GammaH(25,[6]) ; M = G.modular_symbols() ; M
            Modular Symbols space of dimension 11 for Congruence Subgroup Gamma_H(25) with H generated by [6] of weight 2 with sign 0 and over Rational Field
            sage: M.cusps()
            [37/75, 1/2, 31/125, 1/4, -2/5, 2/5, -1/5, 1/10, -3/10, 1/15, 7/15, 9/20]
            sage: len(M.cusps())
            12

        This is always one more than the associated space of weight 2 Eisenstein
        series.

        ::

            sage: G.dimension_eis(2)
            11
            sage: M.cuspidal_subspace()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 11 for Congruence Subgroup Gamma_H(25) with H generated by [6] of weight 2 with sign 0 and over Rational Field
            sage: G.dimension_cusp_forms(2)
            0
        """
        from sage.modular.arithgroup.all import is_GammaH
        if not isinstance(other, Cusp):
            other = Cusp(other)
        if not is_GammaH(G):
            raise TypeError, "G must be a group GammaH(N)."

        H = G._list_of_elements_in_H()
        N = ZZ(G.level())
        u1 = self.__a
        v1 = self.__b
        u2 = other.__a
        v2 = other.__b
        g = v1.gcd(N)

        for h in H:
            v_tmp = (h*v1)%N
            u_tmp = (h*u2)%N
            if (v_tmp - v2)%N == 0 and (u_tmp - u1)%g == 0:
                return True, 1
            if (v_tmp + v2)%N == 0 and (u_tmp + u1)%g == 0:
                return True, -1
        return False, 0

    def apply(self, g):
        """
        Return g(self), where g=[a,b,c,d] is a list of length 4, which we
        view as a linear fractional transformation.

        EXAMPLES: Apply the identity matrix::

            sage: Cusp(0).apply([1,0,0,1])
            0
            sage: Cusp(0).apply([0,-1,1,0])
            Infinity
            sage: Cusp(0).apply([1,-3,0,1])
            -3
        """
        return Cusp(g[0]*self.__a + g[1]*self.__b, g[2]*self.__a + g[3]*self.__b)

    def galois_action(self, t, N):
        r"""
        Suppose this cusp is `\alpha`, `G` is a congruence subgroup of
        level `N`, and `\sigma` is the automorphism in the Galois
        group of `\QQ(\zeta_N)/\QQ` that sends `\zeta_N` to
        `\zeta_N^t`.  Then this function computes a cusp `\beta` such
        that `\sigma([\alpha]) = [\beta]`, where `[\alpha]` is the
        equivalence class of `\alpha` modulo `G`.

        INPUT:

           - `t` -- integer that is coprime to N

           - `N` -- positive integer (level)

        OUTPUT:

           - a cusp

        EXAMPLES::

            sage: Cusp(1/10).galois_action(3, 50)
            1/170
            sage: Cusp(oo).galois_action(3, 50)
            Infinity
            sage: Cusp(0).galois_action(3, 50)
            0

        Here we compute explicitly the permutations of the action for
        t=3 on cusps for Gamma0(50)::

            sage: N = 50; t=3; G = Gamma0(N); C = G.cusps()
            sage: cl = lambda z: exists(C, lambda y:y.is_gamma0_equiv(z, N))[1]
            sage: for i in range(5): print i, t^i, [cl(alpha.galois_action(t^i,N)) for alpha in C]
            0 1 [0, 1/25, 1/10, 1/5, 3/10, 2/5, 1/2, 3/5, 7/10, 4/5, 9/10, Infinity]
            1 3 [0, 1/25, 7/10, 2/5, 1/10, 4/5, 1/2, 1/5, 9/10, 3/5, 3/10, Infinity]
            2 9 [0, 1/25, 9/10, 4/5, 7/10, 3/5, 1/2, 2/5, 3/10, 1/5, 1/10, Infinity]
            3 27 [0, 1/25, 3/10, 3/5, 9/10, 1/5, 1/2, 4/5, 1/10, 2/5, 7/10, Infinity]
            4 81 [0, 1/25, 1/10, 1/5, 3/10, 2/5, 1/2, 3/5, 7/10, 4/5, 9/10, Infinity]

        REFERENCES:

            - Section 1.3 of Glenn Stevens, "Arithmetic on Modular Curves"

            - There is a long comment about our algorithm in the source code for this function.

        WARNING: In some cases `N` must fit in a long long, i.e., there
                 are cases where this algorithm isn't fully implemented.

        AUTHORS:

            - William Stein, 2009-04-18

        """
        if self.is_infinity() or not self.__a: return self
        if not isinstance(t, Integer): t = Integer(t)

        # Our algorithm for computing the Galois action works as
        # follows (see Section 1.3 of Glenn Stevens "Arithmetic on
        # Modular Curves" for a proof that the action given below is
        # correct).  We alternatively view the set of cusps as the
        # Gamma-equivalence classes of column vectors [a;b] with
        # gcd(a,b,N)=1, and the left action of Gamma by matrix
        # multiplication.  The action of t is induced by [a;b] |-->
        # [a;t'*b], where t' is an inverse mod N of t.  For [a;t'*b]
        # with gcd(a,t'*b)==1, the cusp corresponding to [a;t'*b] is
        # just the rational number a/(t'*b).  Thus in this case, to
        # compute the action of t we just do a/b <--> [a;b] |--->
        # [a;t'*b] <--> a/(t'*b).  IN the other case when we get
        # [a;t'*b] with gcd(a,t'*b) != 1, which can and does happen,
        # we have to work a bit harder.  We need to find [c;d] such
        # that [c;d] is congruent to [a;t'*b] modulo N, and
        # gcd(c,d)=1.  There is a standard lifting algorithm that is
        # implemented for working with P^1(Z/NZ) [it is needed for
        # modular symbols algorithms], so we just apply it to lift
        # [a,t'*b] to a matrix [A,B;c,d] in SL_2(Z) with lower two
        # entries congruent to [a,t'*b] modulo N.  This exactly solves
        # our problem, since gcd(c,d)=1.

        a = self.__a
        b = self.__b * t.inverse_mod(N)
        if b.gcd(a) != 1:
            _,_,a,b = lift_to_sl2z_llong(a,b,N)
            a = Integer(a); b = Integer(b)

        # Now that we've computed the Galois action, we efficiently
        # construct the corresponding cusp as a Cusp object.
        return Cusp(a,b,check=False)
