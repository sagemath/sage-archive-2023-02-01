r"""
$p$-adic Numbers

We represent a $p$-adic number as a product $p^r \cdot u$, where $r$ is
an integer or $\infty$ and $u$ is a $p$-adic unit known to some
precision.  If $r = \infty$, then $p^r\cdot u$ is the $0$ element.

Binary operations on two elements of $\Q_p$ reduce the precision of
the unit part of the argument with larger precision to that of the one
with lesser precision.  This applies to all operations, including
equality testing, so, e.g., the element $O(p)$ is equal to every
$p$-adic integer that is divisible by $p$, since comparison will truncate
the other $p$-adic integer to precision $O(p)$.

AUTHORS:
    - William Stein (2004): first version
    - David Joyner (2005-12-24): examples
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

from sage.libs.all import pari, pari_gen, PariError
from sage.structure.element import Element
from sage.misc.latex import latex
from infinity import infinity

import arith
import integer
import integer_mod
import padic_field
import sage.rings.arith
import rational
import field_element
from rational_field import frac, QQ

class pAdic(field_element.FieldElement):
    r"""
    A $p$-adic numbers of either finite or infinite precision.

    EXAMPLES:
      Coercion from rational field using the default precision of the
      p-adic field:
        sage: K = pAdicField(5, 4)
        sage: K(3)         # has an exact representation
          3 + ... + O(5^Infinity)
        sage: K(3/5)
          3*5^-1 + ... + O(5^Infinity)
        sage: K(1/3)       # needs to be approximated; use default precision
          2 + 3*5 + 5^2 + 3*5^3 + O(5^4)
        sage: K(5/3)
          2*5 + 3*5^2 + 5^3 + O(5^4)
        sage: K(1/15)
          2*5^-1 + 3 + 5 + 3*5^2 + 5^3 + O(5^4)

      Coercion from rational field with explicitly selected precision:
        sage: K = pAdicField(5, 4)
        sage: K(3, infinity)
          3 + ... + O(5^Infinity)
        sage: K(3, 8)
          3 + O(5^8)
        sage: K(3/5, 8)
          3*5^-1 + O(5^8)
        sage: K(1/3, 8)
          2 + 3*5 + 5^2 + 3*5^3 + 5^4 + 3*5^5 + 5^6 + 3*5^7 + O(5^8)
        sage: K(5/3, 8)
          2*5 + 3*5^2 + 5^3 + 3*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^8)
        sage: K(1/15, 8)
          2*5^-1 + 3 + 5 + 3*5^2 + 5^3 + 3*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^8)

    Saving and loading p-adics:
        sage: a = 1 + 17 + 3*17^2 + O(17^3)
        sage: a
        1 + 17 + 3*17^2 + O(17^3)
        sage: b = loads(a.dumps())
        sage: b == a
        True
        sage: b
        1 + 17 + 3*17^2 + O(17^3)
    """
    def __init__(self, parent, x, big_oh=infinity, ordp=None, construct=False):
        r"""
        INPUT:
            parent -- a p-adic field.
            x -- anything that can be coerced to a p-adic number.
            big_oh -- is such that, e.g. 3^(-1) + 1 + 2 * 3 + (3^2)
                       has big_oh equal to 2.
        """
        field_element.FieldElement.__init__(self, parent)
        if construct:
            self.__parent = parent
            (self.__p, self.__unit, self.__ordp, self.__prec) = x
            return

        if isinstance(x, pAdic):
            self.__parent = x.__parent
            self.__p = x.__p
            self.__unit = x.__unit
            self.__ordp = x.__ordp
            if big_oh == infinity:
                self.__prec = x.__prec
            else:
                if x.__ordp != infinity:
                    self.__prec = min(x.__prec, big_oh-x.__ordp)
                else:
                    self.__prec = min(x.__prec, big_oh)
            return

        self.__p = int(parent.residue_characteristic())
        self.__parent = parent

        if isinstance(x, pari_gen) and x.type() == "t_PADIC":
            t = x.lift()
            big_oh = x.padicprec(self.__p)
            if t.type() == 't_INT':
                x = int(t)
            else:
                x = QQ(t)

            # we then use code below (so don't make the next line elif)

        if isinstance(x, integer.Integer):
            if x >= 0:
                # Can represent the integer exactly
                x = int(x)
            else:
                # Need to approximate.
                # todo: I'm being lazy here, just use the rational
                # code for the moment.... all this will get rewritten
                # anyway....
                x = rational.Rational(x)

        if isinstance(x, rational.Rational):
            if x == 0:
                x=0
            else:
                ordp = x.valuation(self.__p)
                ppow = rational.Rational(self.__p)**ordp
                x = x/ppow
                if x.denominator() == 1 and x.numerator() >= 0:
                    # Can represent the rational exactly
                    x = int(x)
                else:
                    # Need to approximate, using either default precision,
                    # or requested precision
                    if big_oh is infinity:
                        pr = parent.prec() - ordp
                    else:
                        pr = big_oh - ordp

                    modulus = self.__p**pr
                    x = x.numerator() * \
                          arith.inverse_mod(x.denominator(), modulus)
                    # x = x % modulus
                    big_oh = pr + ordp
        if x==0:
            unit = 1
            ordp = infinity
        elif ordp != None:
            unit = int(x)
        else:  # ord_p is None, so we compute it from x, which need not be a unit.
            p = self.__p
            if isinstance(x, (int, long)):
                ordp = sage.rings.arith.valuation(x,p)
                unit = x/(p**ordp)
            elif isinstance(x, rational.Rational):
                ordp = x.valuation(p)
                unit = int(x/(p**ordp))
            else:
                raise TypeError, "unable to compute ordp of element of type %s"%type(x)
        self.__unit = unit
        self.__ordp = ordp
        if big_oh is infinity and ordp is infinity:
            self.__prec = infinity
        else:
            if ordp != infinity:
                self.__prec = big_oh - ordp  # since internal prec is of the unit part
            else:
                self.__prec = big_oh
        self.__order = None

    def __reduce__(self):
        return reduce_version0, (self.__ordp, self.__parent, self.__p, self.__unit, self.__prec)

    def prec(self):
        """
        Return the precision to which this p-adic number is known.
        """
        return self.__prec

    def _pari_init_(self):
        """
        PARI representation of a p-adic is the same as in SAGE.
        """
        if self.__prec == infinity:
            if self.__ordp == infinity:
                return '0'
            else:
                return str(self.__unit * (self.__p**self.__ordp))
        else:
            return str(self)

    def sqrt(self):
        """
        Return a square root of self.  See the documentation of
        self.square_root for examples.
        """
        return self.square_root()

    def square_root(self):
        """
        Return a square root of self in the $p$-adics.

        EXAMPLES:
            sage: n = 1 +17 + O(17^2); n
            1 + 17 + O(17^2)
            sage: m = n.square_root(); m
            1 + 9*17 + O(17^2)
            sage: m^2
            1 + 17 + O(17^2)

            sage: K = pAdicField(3,prec=10)
            sage: x = K(22/7); x
            1 + 2*3 + 2*3^3 + 3^4 + 2*3^5 + 3^7 + 2*3^9 + O(3^10)
            sage: y = x.sqrt(); y
            1 + 3 + 3^2 + 3^3 + 3^5 + 3^6 + 2*3^7 + 2*3^8 + 3^9 + O(3^10)
            sage: y^2
            1 + 2*3 + 2*3^3 + 3^4 + 2*3^5 + 3^7 + 2*3^9 + O(3^10)

            sage: x = x/9; x
            3^-2 + 2*3^-1 + 2*3 + 3^2 + 2*3^3 + 3^5 + 2*3^7 + O(3^8)
            sage: y = x.sqrt(); y
            3^-1 + 1 + 3 + 3^2 + 3^4 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^9)
            sage: y^2
            3^-2 + 2*3^-1 + 2*3 + 3^2 + 2*3^3 + 3^5 + 2*3^7 + O(3^8)
        """
        try:
            return self.parent()(self._pari_().sqrt())
        except PariError:
            raise ValueError, "square root not a padic number"


    def denominator(self):
        """
        Return the denominator of this p-adic number, which is an integer
        that is a power of $p$.

        EXAMPLES:
            sage: K = Qp(11, 10)
            sage: a = K(211/17); a
            4 + 4*11 + 11^2 + 7*11^3 + 9*11^5 + 5*11^6 + 4*11^7 + 8*11^8 + 7*11^9 + O(11^10)
            sage: a.denominator()
            1
            sage: b = K(3211/11^2); b
            10*11^-2 + 5*11^-1 + 4 + 2*11 + ... + O(11^Infinity)
            sage: b.denominator()
            121
        """
        return self.__p**(-min(self.__ordp,0))

    def ordp(self):
        """
        Return the $p$-adic valuation at this $p$-adic number, normalized
        so that the valuation of $p$ is $1$.

        EXAMPLES:
            sage: K = Qp(11)
            sage: K.prec()
            20
            sage: a = K(211/17); a
            4 + 4*11 + 11^2 + 7*11^3 + 9*11^5 + 5*11^6 + 4*11^7 + 8*11^8 + 7*11^9 + 9*11^10 + 3*11^11 + 10*11^12 + 11^13 + 5*11^14 + 6*11^15 + 2*11^16 + 3*11^17 + 11^18 + 7*11^19 + O(11^20)
            sage: a.ordp()
            0
            sage: b = K(3211/11^2); b
            10*11^-2 + 5*11^-1 + 4 + 2*11 + ... + O(11^Infinity)
            sage: b.ordp()
            -2

            sage: K = Qp(11, prec=5)
            sage: a = K(211/17); a
            4 + 4*11 + 11^2 + 7*11^3 + O(11^5)
        """
        return self.__ordp

    def valuation(self):
        """
        Same as ordp.
        """
        return self.__ordp

    def additive_order(self):
        """
        The additive order of self as an element of the additive group.

        EXAMPLES:
            sage: K = Qp(11)
            sage: a = K(1); a
            1
            sage: a.additive_order()
            Infinity
            sage: b = zero(K); b
            0
            sage: b.additive_order()
            1
        """
        if self.is_zero():
            return integer.Integer(1)
        else:
            return infinity

    def multiplicative_order(self):
        """
        The multiplicative order of this as an element, if defined.

        If the element is known to infinite precision, then it is
        truncated to the parent default precision before the order is
        computed.

            sage: K = Qp(13)
            sage: a = K(-1)
            sage: a.multiplicative_order()
            2

        We immediately know that 13 has infinite order, since it is 0
        modulo 13.
            sage: b = K(13)
            sage: b.multiplicative_order()
            Infinity

        The following element has finite multiplicative order modulo $5^2$:
            sage: c = 3 + 3*5 + 2*5**2 + O(5**3)
            sage: c.multiplicative_order()
            4
        """
        if self.__ordp != 0:
            return infinity
        if self == 1:
            return integer.Integer(1)
        if self == -1:
            return integer.Integer(2)
        if self.__prec == infinity:
            prec = min(self.__prec, self.parent().prec())
        else:
            prec = self.__prec
        return integer_mod.Mod(self.__unit, self.__p**prec).multiplicative_order()

    def copy(self):
        return pAdic(self.__parent, self.__unit, self.__prec, self.__ordp)

    def big_oh(self):
        if self.__ordp == infinity:
            return self.__prec
        return self.__prec + self.__ordp

    def _repr_(self, do_latex=False):
        if not self.parent().series_print():
            if do_latex:
                return "%s^{%s} * (%s + O(%s^{%s}))"%(self.__p, self.__ordp, \
                                                  self.__unit, self.__p, latex(self.__prec))
            else:
                return "%s^%s * (%s + O(%s^%s))"%(self.__p, self.__ordp, \
                                                  self.__unit, self.__p, self.__prec)
        # series printing
        if self.__ordp == infinity:
            o = self.big_oh()
            if o == infinity:
                return "0"
            else:
                return "O(%s^%s)"%(self.__p, o)
        if self.__ordp == 0 and self.__prec == infinity and self.__unit == 1:
            return "1"
        s     = ""
        u     = self.__unit
        exp   = self.__ordp
        p     = self.__p
        if self.__ordp == infinity:
            prec = self.__prec
        else:
            prec  = min(self.__prec, self.parent().print_prec() - self.__ordp)
        if prec == infinity:
            prec = self.parent().prec()
        u   %= self.__p ** prec
        while u != 0:
            coeff = u % p
            if coeff != 0:
                if exp == 0:
                    s += "%s + "%coeff
                else:
                    var = "%s"%p
                    if exp != 1:
                        if do_latex:
                            var += "^{%s}"%exp
                        else:
                            var += "^%s"%exp
                    if coeff != 1:
                        if do_latex:
                            s += "%s\\cdot %s + "%(coeff,var)
                        else:
                            s += "%s*%s + "%(coeff,var)
                    else:
                        s += "%s + "%var
            exp += 1
            u = (u-coeff)/p
        if self.big_oh() == infinity and exp != infinity:
            s += '... + '
        s += "O(%s"%(p)
        if self.big_oh() == 1:
            s += ")"
        else:
            if do_latex:
                s += "^{%s})"%latex(self.big_oh())
            else:
                s += "^%s)"%self.big_oh()
        return s

    def _latex_(self):
        return self._repr_(do_latex=True)

    def _add_(self, right):
        """
        EXAMPLES:
            sage: K = Qp(11, 10); K.print_prec(5)
            sage: a = K(-1); a
            10 + 10*11 + 10*11^2 + 10*11^3 + 10*11^4 + O(11^10)
            sage: b = K(1); b
            1
            sage: a+b
            O(11^10)
        """
        if self.__ordp <= right.__ordp:
            x = self; y = right
        else:
            x = right; y = self
        big_oh = min(x.big_oh(), y.big_oh())
        if y == 0:
            return pAdic(self.__parent, x, big_oh)
        p = x.__p
        a = x.__unit + p**(y.__ordp - x.__ordp)*y.__unit
        n = arith.valuation(a,p)
        if n == infinity:    # 0
            return pAdic(self.__parent, 0, big_oh)
        a /= p**n
        prec = big_oh - n - x.__ordp
        if prec != infinity:
            a %= p**prec
        if a == 0:
            return pAdic(self.__parent, 0, big_oh)
        return pAdic(self.__parent, a, big_oh, x.__ordp + n)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: K = Qp(19, 5)
            sage: zero(K) - one(K)
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
        """
        return self + (-right)

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: K = Qp(19, 5)
            sage: (-1)*one(K)
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a = K(2/19); a
            2*19^-1 + ... + O(19^Infinity)
            sage: b = K(3/19); b
            3*19^-1 + ... + O(19^Infinity)
            sage: a*b
            6*19^-2 + ... + O(19^Infinity)
        """
        prec = min(self.__prec, right.__prec)
        ordp = self.__ordp + right.__ordp
        unit = self.__unit * right.__unit
        if prec != infinity:
            unit %= self.__p**prec
        return pAdic(self.__parent, unit, prec+ordp, ordp)

    def _div_(self, right):
        """
        EXAMPLES:
            sage: K = Qp(19, 5)
            sage: a = K(2/19); a
            2*19^-1 + ... + O(19^Infinity)
            sage: b = K(3/19); b
            3*19^-1 + ... + O(19^Infinity)
            sage: a/b
            7 + 6*19 + 6*19^2 + 6*19^3 + 6*19^4 + O(19^5)
            sage: (5 + O(5)) / 1
            O(5^1)
        """
        return self * right.__invert__(self.__prec + max(self.__ordp,0))

    def __mod__(self, right):
        """
        EXAMPLES:
            sage: a = 2 + 13*17 + 15*17^2 + O(17^3)
            sage: a % 17^2
            2 + 13*17 + O(17^2)

            sage: a % 17
            2 + O(17)

        We reduce a 16-th
            sage: zeta = pAdicField(17,8).zeta(16); zeta
            14 + 3*17 + 14*17^2 + 13*17^3 + 16*17^4 + 5*17^5 + 12*17^6 + 16*17^7 + O(17^8)
            sage: zeta^16
            1 + O(17^8)
            sage: zeta % 17
            14 + O(17)
            sage: zeta % 17^3
            14 + 3*17 + 14*17^2 + O(17^3)
        """
        n = int(right)
        if n != right:
            raise TypeError, "modulus must be an integer"
        p = self.__p
        v = arith.valuation(n, p)
        if p**v != n:
            raise ValueError, "modulus must be a power of p (=%s)"%p
        return self + self.__parent(0, v)

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: K = Qp(19, 5)
            sage: a = K(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a^2
            1 + O(19^5)
            sage: a^3
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: K(5)^30
            11 + 14*19 + 19^2 + 7*19^3 + ... + O(19^Infinity)
        """
        right = integer.Integer(right)
        if self == 0:
            if right == 0:
                raise ValueError, "0^0 not defined"
            return 0
        if right < 0:
            inv = 1/self
            return inv**(-right)
        if right == 0:
            return pAdic(self.__parent, 1)
        ordp = right * self.__ordp
        if self.__prec == infinity:
            z = pAdic(self.__parent, self.__unit**right)
            z.__ordp = ordp
            return z
        else:
            prec = self.__prec
        unit = arith.power_mod(self.__unit, right, self.__p**prec)
        return pAdic(self.__parent, unit, prec + ordp, ordp)

    def __neg__(self):
        """
        EXAMPLES:
           sage: K = Qp(5, 3)
           sage: K(1)
           1
           sage: -K(1)
           4 + 4*5 + 4*5^2 + O(5^3)
        """
        u = - self.__unit
        if self.__prec is infinity:
            prec = self.__parent.prec()
        else:
            prec = self.__prec
        return pAdic(self.__parent, u, prec + self.__ordp, self.__ordp)

    def __pos__(self):
        return self

    def __invert__(self, prec=infinity):
        r"""
        EXAMPLES:
            sage: K = Qp(19, 5)
            sage: a = K(20); a
            1 + 19 + ... + O(19^Infinity)
            sage: b = ~a    # calls __invert__
            sage: b
            1 + 18*19 + 18*19^3 + O(19^5)
            sage: a*b
            1 + O(19^5)

        One can pass an optional argument to __invert__ to
        affect the precision, which is especially useful when
        the input has infinite precision:
            sage: K = pAdicField(5,10)
            sage: 1/K(2)
            3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + O(5^10)
            sage: K(2).__invert__()
            3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + O(5^10)
            sage: K(2).__invert__(15)
            3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + 2*5^10 + 2*5^11 + 2*5^12 + 2*5^13 + 2*5^14 + O(5^15)

        Check we don't add extra precision:
            sage: 1/(1+O(5^2))
             1 + O(5^2)

        AUTHORS:
            -- David Harvey (2006-09-27): fixed precision bug (trac \#89)
        """
        if self.__prec is infinity:
            # If we have infinite precision, approximate the answer to the
            # requested precision (or use the parent's default).
            if prec is infinity:
                prec = self.parent().prec()
        else:
            # Otherwise, use our current precision, or honour a request
            # for lower precision.
            prec = min(self.__prec, prec)

        if prec <= 0:
            raise ZeroDivisionError, "Can not invert self"

        unit = arith.inverse_mod(self.__unit, self.__p**prec)
        big_oh = prec - self.__ordp
        return pAdic(self.__parent, unit, big_oh,  -self.__ordp)

    def lift(self):
        """
        Return rational number with denominator only divisible by
        $p$ that lifts this $p$-adic number, to the given precision.

        EXAMPLES:
            sage: a = 4596/7^2 + O(7^4); a
            4*7^-2 + 5*7^-1 + 2 + 6*7 + 7^2 + O(7^4)
            sage: a.lift()
            4596/49

            sage: K = Qp(19, 5)

            sage: a = K(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a.lift()
            2476098

            sage: a = 4596/18 + O(7^4); a
            1 + 6*7 + 2*7^2 + 5*7^3 + O(7^4)
            sage: a.lift()
            1856

            sage: z = 7 + 4*17 + 2*17^2 + 16*17^3 + 17^4 + 15*17^5 + O(17^6)
            sage: lift(z)
            21460637
            sage: lift(z%17)
            7
            sage: w = z/17^2; w
            7*17^-2 + 4*17^-1 + 2 + 16*17 + 17^2 + 15*17^3 + O(17^4)
            sage: lift(w%17)
            653/289
            sage: 7*17^-2 + 4*17^-1 + 2
            653/289
        """
        if self.is_zero():
            return frac(0,1)
        p = self.__p
        if self.__prec == infinity:
            return frac(p,1)**self.__ordp * frac(self.__unit,1)
        return frac(p,1)**self.__ordp * (self.__unit % (p**self.__prec))

    def _integer_(self):
        return self.lift()

    def __cmp__(self, other):
        """
        First compare valuations, then compare normalized
        residue of unit part.

        EXAMPLES:
            sage: K = Qp(19, 5)
            sage: a = K(2); a
            2 + ... + O(19^Infinity)
            sage: b = K(3); b
            3 + ... + O(19^Infinity)
            sage: a < b
            True
        """
        m = min(self.big_oh(), other.big_oh())
        x_ordp = self.__ordp
        if x_ordp >= m :
            x_ordp = infinity
        y_ordp = other.__ordp
        if y_ordp >= m :
            y_ordp = infinity
        if x_ordp < y_ordp:
            return -1
        elif x_ordp > y_ordp:
            return 1
        else:  # equal ordp
            if x_ordp == infinity:
                return 0 # since both are zero
        p = self.__p
        a = self.__unit
        b = other.__unit
        prec = min(self.__prec, other.__prec)
        if prec != infinity:
            ppow = p**prec
            a %= ppow
            b %= ppow
        if a < b:
            return -1
        elif a == b:
            return 0
        else:
            return 1

    def unit_part(self):
        r"""
        Return the unit part of $p^r\cdot u$, which is simply u.

        EXAMPLES:
            sage: x = 9*(2+3+O(3**7))
            sage: x.unit_part()
            2 + 3 + O(3^7)
            sage: K = Qp(19, 5)
            sage: a = K(2)/19; a
            2*19^-1 + O(19^4)
            sage: a.unit_part()
            2 + O(19^5)
        """
        return pAdic(self.__parent, self.__unit, self.__prec, 0)

    def is_zero(self):
        """
            EXAMPLES:
            sage: K = Qp(11)
            sage: a = K(2); a
            2 + ... + O(11^Infinity)
            sage: a.is_zero()
            False

            sage: K = Qp(11)
            sage: K(0).is_zero()
            True
        """
        return self.__ordp >= self.big_oh()

    def is_unit(self):
        """
        EXAMPLES:
            sage: K = Qp(11)
            sage: a = K(2); a
            2 + ... + O(11^Infinity)
            sage: a.is_unit()
            True
            sage: K(121).is_unit()
            False
            sage: K(0).is_unit()
            False
        """
        return self.__ordp == 0

    def rational_reconstruction(self):
        r"""
        Try to lift the p-adic number to the rationals using rational
        reconstruction, as follows: Suppose the p-adic number is
        $p^{r}\cdot (u+O(p^n))$, where u is a unit.  Using rational
        reconstruction, try to find the unique rational number a/b
        such that a/b is congruent to u modulo $p^n$, and abs(a),
        abs(b) are both at most $\sqrt{p/2}$.  If such $a/b$ exists,
        return $p^r \cdot (a/b)$, otherwise raises ValueError.

        EXAMPLES:
            sage: K = Qp(11, 10); K.print_prec(10)
            sage: a = K(-1); a
            10 + 10*11 + 10*11^2 + 10*11^3 + 10*11^4 + 10*11^5 + 10*11^6 + 10*11^7 + 10*11^8 + 10*11^9 + O(11^10)
            sage: a.rational_reconstruction()
            -1
            sage: a = K(2); a
            2 + ... + O(11^Infinity)
            sage: a.rational_reconstruction()
            2
        """
        if self.is_zero():
            return frac(0,1)
        p = self.__p
        if self.__prec == infinity:
            return frac(p,1)**self.__ordp * frac(self.__unit,1)
        alpha = self.__unit
        m = integer.Integer(p**self.__prec)
        r = arith.rational_reconstruction(alpha, m)
        return (frac(p,1)**self.__ordp)*r

    def log(self):
        r"""
        Compute the p-adic logarithm of any nonzero element of $\Q_p$.
        (See below for normalization.)

        The usual power series for log with values in the additive
        group of $\Q_p$ only converges for 1-units (units congruent to
        1 modulo p).  However, there is a unique extension of log to a
        homomorphism defined on all the units.  If u = a*v is a unit
        with v = 1 (mod p), then we define log(u) = log(v).  This is
        the correct extension because the units U of Z_p splits as a
        product U = V x <w>, where V is the subgroup of 1-units and w
        is a (p-1)st root of unity.  The <w> factor is torsion, so
        must go to 0 under any homomorphism to the torsion free group
        $(\Q_p, +)$.

        Notes -- What some other systems do:
           PARI:  Seems to define log the same way as we do.
           MAGMA: Gives an error when unit is not a 1-unit.

        Algorithm:
           Input: Some p-adic unit u.
           1. Check that the input p-adic number is really a unit
              (i.e., valuation 0)
           2. Let $1-x = u^{p-1}$, which is a 1-unit.
           3. Use the series expansion
              $$
                \log(1-x) = F(x) = -x - 1/2*x^2 - 1/3*x^3 - 1/4*x^4 - 1/5*x^5 - ...
              $$
              to compute the logarithm log(u**(p-1)).  Use enough
              terms so that terms added on are zero (to the default
              precision, if the input has infinite precision).
           4. Then $$\log(u) = log(u^{p-1})/(p-1) = F(1-u^{p-1})/(p-1).$$

        EXAMPLES:
            sage: Q13 = Qp(13, 10)
            sage: a = Q13(14); a
            1 + 13 + ... + O(13^Infinity)
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)

        The next few examples illustrate precision when computing $p$-adic logs.
        First we create a field with {\em default} precision 10.
            sage: K = pAdicField(5,10)
            sage: e = K(389); e
            4 + 2*5 + 3*5^3 + ... + O(5^Infinity)

        If the input has infinite precision, the output is computed to the
        default precision of the parent field.
            sage: e.log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)

        In contrast, if the input has larger precision than the default, then
        the output is also computed to larger precision.
            sage: (e+O(5^30)).log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + 2*5^10 + 2*5^11 + 5^12 + 5^14 + 3*5^15 + 5^16 + 4*5^18 + 4*5^19 + 5^20 + 3*5^21 + 2*5^22 + 2*5^23 + 5^24 + 5^25 + 3*5^27 + 2*5^28 + 3*5^29 + O(5^30)

        Check that results are consistent over a range of precision:
            sage: max_prec = 40
            sage: p = 3
            sage: K = pAdicField(p, max_prec)
            sage: full_log = (1 + p + O(p**max_prec)).log()
            sage: for prec in range(2, max_prec):
            ...       assert (1 + p + O(p**prec)).log() == full_log


        AUTHORS:
            -- David Harvey (2006-09-13): corrected subtle precision bug
               (need to take denominators into account! -- see trac \#53)

        TODO:
            -- Surely there must exist an $O(N^{1+\epsilon})$ time algorithm
            for this? The current one is pretty much quadratic in the
            precision. (Yes: Kiran pointed out the paper by Bernstein:
            http://cr.yp.to/lineartime/multapps-20041007.pdf -- David)

        """

        # Step 1 -- a unit?
        if self.is_unit() and self.__unit % self.__p == 1:
            # It's already a 1-unit, so just use the series
            # (base case of "induction")

            if self.__prec == infinity:
                prec = self.parent().prec()
            else:
                prec = self.__prec

            # Need extra precision to take into account powers of p
            # in the denominators of the series. (Indeed, it's a
            # not-entirely-trivial fact that if x is given mod p^n, that
            # log(x) is well-defined mod p^n !) Specifically:
            # we are only guaranteed that $x^j/j$ is zero mod $p^n$ if
            # j >= floor(log_p(j)) + n.
            extra_prec = integer.Integer(0)
            while extra_prec < (prec + extra_prec).exact_log(self.__p):
                extra_prec += 1

            x = 1 - self
            p = self.__p
            R = sage.rings.all.Integers(p**(prec + extra_prec))
            x = R(x.lift())

            xpow = x
            ans = R(0)
            for n in range(1, prec + extra_prec):
                ans = ans - R(xpow.lift()/n)
                xpow = xpow * x

            return pAdic(self.__parent, ans.lift(), prec)
        else:
            z = self.unit_part()
            p = self.__p
            return (z**(p-1)).log()/(p-1)

    def algdep(self, n):
        """
        Returns a polynomial of degree at most $n$ which is approximately
        satisfied by this number.  Note that the returned polynomial
        need not be irreducible, and indeed usually won't be if this number
        is a good approximation to an algebraic number of degree less than $n$.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLE:
        sage: K = pAdicField(3)
        sage: a = K(7/19); a
        1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
        sage: a.algdep(1)
        19*x - 7
        """
        return sage.rings.arith.algdep(self, n)

pAdic.algebraic_dependency = pAdic.algdep


def reduce_version0(ordp, parent, p, unit, prec):
    return pAdic(parent, (p, unit, ordp, prec), construct=True)
