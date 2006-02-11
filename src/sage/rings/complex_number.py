"""
Complex Numbers

AUTHOR:
    -- William Stein (2006-01-26): complete rewrite
"""

#################################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import math
import operator

from coerce import bin_op
import complex_field
import sage.rings.ring_element as ring_element
import real_field
import sage.misc.misc
import sage.libs.pari.all as pari
import integer
import infinity
#from sage.databases.odlyzko import zeta_zeroes

def is_ComplexNumber(x):
    return isinstance(x, ComplexNumber)

class ComplexNumber(ring_element.RingElement):
    """
    A complex number.

    EXAMPLES:
        sage: C = ComplexField()
        sage: I = C.0
        sage: b = 1.5 + 2.5*I
        sage: loads(b.dumps()) == b
        True
    """
    def __init__(self, parent, real, imag=None):
        ring_element.RingElement.__init__(self, parent)
        if imag is None:
            if isinstance(real, ComplexNumber):
                real, imag = real.__re, real.__im
            elif isinstance(real, list) or isinstance(real,tuple):
                real, imag = real
            elif isinstance(real, complex):
                real, imag = real.real, real.imag
            elif isinstance(real, pari.gen):
                x = real
                real = str(x.real()).replace(' E','e')
                imag = str(x.imag()).replace(' E','e')
            else:
                imag = 0
        try:
            R = parent._real_field()
            self.__re = R(real)
            self.__im = R(imag)
        except TypeError:
            raise TypeError, "unable to coerce (%s,%s) to a ComplexNumber"%(real,imag)
        self.__repr = None

    def _repr_(self):
        s = ""
        if self.__re != 0:
            s = str(self.__re)
        if self.__im != 0:
            y  =  self.__im
            if s!="":
                if y < 0:
                    s = s+" - "
                    y = -y
                else:
                    s = s+" + "
            s = s+"%s*I"%y
        if len(s) == 0:
            s = "0"
        return s

    def _pari_(self):
        try:
            return self.__pari
        except AttributeError:
            self.__pari = pari.pari.new_with_bits_prec(str(self), self.prec())
            return self.__pari

    def _add_(self, right):
        return ComplexNumber(self.parent(), self.__re + right.__re, self.__im + right.__im)

    def _sub_(self, right):
        return ComplexNumber(self.parent(), self.__re - right.__re, self.__im - right.__im)

    def _mul_(self, right):
        return ComplexNumber(self.parent(), self.__re * right.__re - self.__im * right.__im,
                             self.__re * right.__im + self.__im * right.__re)

    def norm(self):
        return self.__re*self.__re + self.__im*self.__im

    def _div_(self, right):
        right_nm = right.norm()
        a = right.__re/right_nm
        b = right.__im/right_nm
        re =  a * self.__re + b * self.__im
        im = -b * self.__re + a * self.__im
        return ComplexNumber(self.parent(), re, im)

    def __rdiv__(self, left):
        return ComplexNumber(self.parent(), left)/self

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: C, i = ComplexField(20).objgen()
            sage: a = i^2; a
            -1.0000000
            sage: a.parent()
            Complex Field with 20 bits of precision
            sage: a = (1+i)^i; a
            0.42882872 + 0.15487170*I
            sage: (1+i)^(1+i)
            0.27395725 + 0.58370018*I
            sage: a.parent()
            Complex Field with 20 bits of precision
            sage: i^i
            0.20787954
            sage: (2+i)^(0.5)
            1.4553452 + 0.34356070*I
        """
        if isinstance(right, (int, long, integer.Integer)):
            return ring_element.RingElement.__pow__(self, right)
        z = self._pari_()
        P = self.parent()
        w = P(right)._pari_()
        m = z**w
        return P(m)

    def prec(self):
        """
        Return square root, which is a complex number.

        EXAMPLES:
            sage: i = ComplexField(2000).0
            sage: i.prec()
            2000
        """
        return self.__re.prec()

    def real(self):
        """
        Return real part of self.

        EXAMPLES:
            sage: i = ComplexField(100).0
            sage: z = 2 + 3*i
            sage: x = z.real(); x
            2.0000000000000000000000000000000
            sage: x.parent()
            Real Field with 100 bits of precision
        """
        return self.__re

    def imag(self):
        """
        Return imaginary part of self.

        EXAMPLES:
            sage: i = ComplexField(100).0
            sage: z = 2 + 3*i
            sage: x = z.imag(); x
            3.0000000000000000000000000000000
            sage: x.parent()
            Real Field with 100 bits of precision
        """
        return self.__im

    def __neg__(self):
        return ComplexNumber(self.parent(), -self.__re, -self.__im)

    def __pos__(self):
        return self

    def __abs__(self):
        R = self.__re.parent()
        x = R(self.__re**2 + self.__im**2)
        return x.sqrt()

    def __invert__(self):
        a = abs(self)
        return ComplexNumber(self.parent(), self.__re/a, -self.__im/a)

    def __int__(self):
        if self.__im == 0:
            return int(self.__re)
        raise TypeError

    def __long__(self):
        if self.__im == 0:
            return long(self.__re)
        raise TypeError

    def __float__(self):
        if self.__im == 0:
            return float(self.__re)
        raise TypeError

    def __complex__(self):
        return complex(float(self.__re), float(self.__im))

    def __cmp__(self, other):
        if not isinstance(other, ComplexNumber):
            other = self.parent()(other)
        return sage.misc.misc.generic_cmp((self.__re, self.__im) , (other.__re, other.__im))

    def multiplicative_order(self):
        """
        Return the multiplicative order of this complex number, if known, or raise
        a NotImplementedError.

        EXAMPLES:
            sage: C, i = ComplexField().objgen()
            sage: i.multiplicative_order()
            4
            sage: C(1).multiplicative_order()
            1
            sage: C(-1).multiplicative_order()
            2
            sage: C(i^2).multiplicative_order()
            2
            sage: C(-i).multiplicative_order()
            4
            sage: C(2).multiplicative_order()
            Infinity
            sage: w = (1+sqrt(-3))/2; w
            0.50000000000000000 + 0.86602540378443860*I
            sage: abs(w)
            0.99999999999999989
            sage: w.multiplicative_order()
            Traceback (most recent call last):
            ...
            NotImplementedError: order of 0.50000000000000000 + 0.86602540378443860*I not known
        """
        if self == 1:
            return integer.Integer(1)
        elif self == -1:
            return integer.Integer(2)
        elif self == self.parent().gen():
            return integer.Integer(4)
        elif self == -self.parent().gen():
            return integer.Integer(4)
        elif hasattr(self, "_multiplicative_order"):
            return self._multiplicative_order
        elif abs(abs(self) - 1) > 0.1:  # clearly not a root of unity
            return infinity.infinity
        raise NotImplementedError, "order of %s not known"%self


    ########################################################################
    # Transcendental (and other) functions
    ########################################################################


    # Trig functions
    def acos(self):
        """
        EXAMPLES:
            sage: (1+I).acos()
            0.90455689430238129 - 1.0612750619050355*I
        """
        return self.parent()(self._pari_().acos())

    def acosh(self):
        """
        EXAMPLES:
            sage: (1+I).acosh()
            1.0612750619050355 + 0.90455689430238129*I
        """
        return self.parent()(self._pari_().acosh())

    def asin(self):
        """
        EXAMPLES:
            sage: (1+I).asin()
            0.66623943249251516 + 1.0612750619050355*I
        """
        return self.parent()(self._pari_().asin())

    def asinh(self):
        """
        EXAMPLES:
            sage: (1+I).asinh()
            1.0612750619050355 + 0.66623943249251516*I
        """
        return self.parent()(self._pari_().asinh())

    def atan(self):
        """
        EXAMPLES:
            sage: (1+I).atan()
            1.0172219678978514 + 0.40235947810852507*I
        """
        return self.parent()(self._pari_().atan())

    def atanh(self):
        """
        EXAMPLES:
            sage: (1+I).atanh()
            0.40235947810852507 + 1.0172219678978514*I
        """
        return self.parent()(self._pari_().atanh())

    def cotan(self):
        """
        EXAMPLES:
            sage: (1+I).cotan()
            0.21762156185440268 - 0.86801414289592493*I
            sage: i = ComplexField(200).0
            sage: (1+i).cotan()
            0.21762156185440268136513424360523807352075436916785404091068113 - 0.86801414289592494863584920891627388827343874994609327121115055*I   # 32-bit
            0.21762156185440268136513424360523807352075436916785404090999986 - 0.86801414289592494863584920891627388827343874994609327120999992*I   # 64-bit
            sage: i = ComplexField(220).0
            sage: (1+i).cotan()
             0.21762156185440268136513424360523807352075436916785404091068124239236 - 0.86801414289592494863584920891627388827343874994609327121115071646176*I     # 32-bit
             0.21762156185440268136513424360523807352075436916785404091068124239992 - 0.86801414289592494863584920891627388827343874994609327121115071649974*I   # 64-bit
        """
        return self.parent()(self._pari_().cotan())

    def cos(self):
        """
        EXAMPLES:
            sage: (1+I).cos()
            0.83373002513114902 - 0.98889770576286506*I
        """
        return self.parent()(self._pari_().cos())

    def cosh(self):
        """
        EXAMPLES:
            sage: (1+I).cosh()
            0.83373002513114902 + 0.98889770576286506*I
        """
        return self.parent()(self._pari_().cosh())

    def sin(self):
        """
        EXAMPLES:
            sage: (1+I).sin()
            1.2984575814159771 + 0.63496391478473602*I
        """
        return self.parent()(self._pari_().sin())

    def sinh(self):
        """
        EXAMPLES:
            sage: (1+I).sinh()
            0.63496391478473602 + 1.2984575814159771*I
        """
        return self.parent()(self._pari_().sinh())

    def tan(self):
        """
        EXAMPLES:
            sage: (1+I).tan()
            0.27175258531951169 + 1.0839233273386943*I
        """
        return self.parent()(self._pari_().tan())

    def tanh(self):
        """
        EXAMPLES:
            sage: (1+I).tanh()
            1.0839233273386943 + 0.27175258531951169*I
        """
        return self.parent()(self._pari_().tanh())

    # Other special functions
    def agm(self, right):
        """
        EXAMPLES:
            sage: (1+I).agm(2-I)
            1.6278054848727064 + 0.13682754839736855*I
        """
        t = self.parent()(right)._pari_()
        return self.parent()(self._pari_().agm(t))


    def argument(self):
        r"""
        The argument (angle) of the complex number, normalized
        so that $-\pi < \theta \leq \pi$.

        EXAMPLES:
            sage: i = ComplexField().gen()
            sage: (i^2).argument()
            3.1415926535897931
            sage: (1+i).argument()
            0.78539816339744828
            sage: i.argument()
            1.5707963267948966
            sage: (-i).argument()
            -1.5707963267948966
            sage: (-0.001 - i).argument()
            -1.5717963264615633
        """
        return self.parent()(self._pari_().arg())

    def arg(self):
        """
        Same as argument.

        EXAMPLES:
            sage: i = ComplexField().gen()
            sage: (i^2).arg()
            3.1415926535897931
        """
        return self.argument()

    def conjugate(self):
        """
        Return the complex conjugate of this complex number.

        EXAMPLES:
            sage: i = ComplexField().gen()
            sage: (1+i).conjugate()
            1.0000000000000000 - 1.0000000000000000*I
        """
        return ComplexNumber(self.parent(), self.__re, -self.__im)

    def dilog(self):
        return self.parent()(self._pari_().dilog())

    def exp(self):
        """
        Compute exp(z).

        EXAMPLES:
            sage: i = ComplexField(300).gen()
            sage: z = 1 + i
            sage: z.exp()
            1.4686939399158851571389675973266042613269567366290087227976756763109369658595121387227244973 + 2.2873552871788423912081719067005018089555862566683556809386581141036471601893454092673448502*I   # 32-bit
            1.4686939399158851571389675973266042613269567366290087227976756763109369658595121387199999995 + 2.2873552871788423912081719067005018089555862566683556809386581141036471601893454092699999982*I   # 64-bit
        """
        return self.parent()(self._pari_().exp())

    def gamma(self):
        """
        Return the Riemann zeta function evaluated at this complex number.

        EXAMPLES:
            sage: i = ComplexField(30).gen()
            sage: (1+i).gamma()
            0.49801566778 - 0.15494982828*I
        """
        return self.parent()(self._pari_().gamma())

    def gamma_inc(self, t):
        """
        Return the Riemann zeta function evaluated at this complex number.

        EXAMPLES:
            sage: C, i = ComplexField(30).objgen()
            sage: (1+i).gamma_inc(2+3*i)
            0.0020969148609 - 0.059981913655*I
            sage: (1+i).gamma_inc(5)
            -0.0013781309353 + 0.0065198200173*I
            sage: C(2).gamma_inc(1+i)
            0.70709209610 - 0.42035364080*I
            sage: gamma_inc(2, 1+i)
            0.70709209610 - 0.42035364080*I
            sage: gamma_inc(2, 5)
            0.040427681994512799
        """
        t = self.parent()(t)._pari_()
        return self.parent()(self._pari_().incgam(t))

    def log(self):
        """
        Complex logarithm of z with branch chosen as follows:
        Write z = rho*exp(i*theta) with -pi <= theta < pi.  Then
               log(z) = log(rho) + i*theta.

        WARNING: Currently the real log is computed using floats, so there is
        potential precision loss.
        """
        theta = self.argument()
        rho = abs(self)
        return ComplexNumber(self.parent(), rho.log(), theta)

    def order(self):
        if self.__re == 1 and self.__im == 0:
            return 1
        if self.__re == -1 and self.__im == 0:
            return 2
        raise NotImplementedError, "general multiplicative order over complex field not implemented."

    def sqrt(self):
        """
        EXAMPLES:
            sage: C, i = ComplexField(30).objgen()
            sage: i.sqrt()
            0.70710678026 + 0.70710678026*I
            sage: (1+i).sqrt()
            1.0986841116 + 0.45508986013*I
            sage: (C(-1)).sqrt()
            1.0000000000*I
            sage: i = ComplexField(200).0
            sage: i.sqrt()
            0.70710678118654752440084436210484903928483593768847403658833981 + 0.70710678118654752440084436210484903928483593768847403658833981*I   # 32-bit
            0.70710678118654752440084436210484903928483593768847403658999949 + 0.70710678118654752440084436210484903928483593768847403658999949*I   # 64-bit

        """
        return self.parent()(self._pari_().sqrt())

    def square_root(self):
        """
        Return square root, which is a complex number.

        EXAMPLES:
            sage: i = ComplexField(100).0
            sage: (-i).sqrt()
            0.70710678118654752440084436210459 - 0.70710678118654752440084436210459*I
        """
        return self.sqrt()

    def zeta(self):
        """
        Return the Riemann zeta function evaluated at this complex number.

        EXAMPLES:
            sage: i = ComplexField(30).gen()
            sage: z = 1 + i
            sage: z.zeta()
            0.58215805888 - 0.92684856430*I
            sage: zeta(z)
            0.58215805888 - 0.92684856430*I
        """
        return self.parent()(self._pari_().zeta())

