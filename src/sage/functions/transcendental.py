"""
Transcendental Functions
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

import  sage.libs.pari.all
pari = sage.libs.pari.all.pari
import sage.rings.complex_field as complex_field
import sage.rings.real_mpfr as real_field
import sage.rings.complex_number

from sage.rings.all import is_RealNumber, RealField, is_ComplexNumber, ComplexField

def __prep_num(x):
    if isinstance(x, sage.rings.complex_number.ComplexNumber):
        x = str(x).replace("i","I")
    else:
        x = str(x)
    return x

CC = complex_field.ComplexField()
I = CC.gen(0)
def __eval(x):
    return eval(x)

def exponential_integral_1(x, n=0):
    r"""
    Returns the exponential integral $E_1(x)$. If the optional argument
    $n$ is given, computes list of the first $n$ values of the exponential
    integral $E_1(x m)$.

    The exponential integral $E_1(x)$ is
    $$
             E_1(x) = \int_{x}^{\infty} e^{-t}/t dt
    $$

    INPUT:
        x -- a positive real number

        n -- (default: 0) a nonnegative integer; if nonzero,
             then return a list of values E_1(x*m) for
             m = 1,2,3,...,n.   This is useful, e.g., when
             computing derivatives of L-functions.

    OUTPUT:
        float -- if n is 0 (the default)
      or
        list -- list of floats if n > 0

    EXAMPLES:
        sage: exponential_integral_1(2)
        0.048900510708061118
        sage: w = exponential_integral_1(2,4); w
        [0.048900510708061118, 0.0037793524098489067, 0.00036008245216265867, 3.76656228439249e-05]


    IMPLEMENTATION: We use the PARI C-library functions eint1 and
    veceint1.

    REFERENCE: See page 262, Prop 5.6.12, of Cohen's book "A Course
    in Computational Algebraic Number Theory".

    REMARKS: When called with the optional argument n, the PARI
    C-library is fast for values of n up to some bound, then very
    very slow.  For example, if x=5, then the computation takes less
    than a second for n=800000, and takes "forever" for n=900000.

    """
    if n <= 0:
        return float(pari(x).eint1())
    else:
        return [float(z) for z in pari(x).eint1(n)]

def gamma(s):
    """
    Gamma function at s.

    EXAMPLES:
        sage: gamma(CDF(0.5,14))
        -4.05370307804e-10 - 5.77329961615e-10*I
        sage: gamma(I)
        -0.154949828301811 - 0.498015668118356*I
        sage: gamma(6)
        120.000000000000
    """
    try:
        return s.gamma()
    except AttributeError:
        return CC(s).gamma()

def gamma_inc(s, t):
    """
    Incomplete Gamma function Gamma(s,t).

    EXAMPLES:
        sage: gamma_inc(CDF(0,1), 3)
        0.00320857499337 + 0.0124061862007*I
        sage: gamma_inc(3, 3)
        0.846380162253687 + 2.52435489670724e-29*I
        sage: gamma_inc(RDF(1), 3)
        0.0497870683678639
    """
    try:
        return s.gamma_inc(t)
    except AttributeError:
        if not (is_ComplexNumber(s)):
            if is_ComplexNumber(t):
                C = t.parent()
            else:
                C = ComplexField()
            s = C(s)
        return s.gamma_inc(t)


# synonym.
incomplete_gamma = gamma_inc

def zeta(s):
    """
    Riemann zeta function at s with s a real or complex number.

    INPUT:
        s -- real or complex number

    If s is a real number the computation is done using the MPFR
    library.  When the input is not real, the computation is done
    using the PARI C library.

    EXAMPLES:
        sage: zeta(2)
        1.64493406684823
        sage: RR = RealField(200)
        sage: zeta(RR(2))
        1.6449340668482264364724151666460251892189499012067984377356
    """
    try:
        return s.zeta()
    except AttributeError:
        return RealField()(s).zeta()

##     prec = s.prec()
##     s = pari.new_with_prec(s, prec)
##     z = s.zeta()._sage_()
##     if z.prec() < prec:
##         raise RuntimeError, "Error computing zeta(%s) -- precision loss."%s
##     return z

def zeta_symmetric(s):
    r"""
    Completed function $\xi(s)$ that satisfies $\xi(s) = \xi(1-s)$ and
    has zeros at the same points as the Riemann zeta function.

    INPUT:
        s -- real or complex number

    If s is a real number the computation is done using the MPFR
    library.  When the input is not real, the computation is done
    using the PARI C library.

    More precisely,
    $$
       xi(s) = \gamma(s/2 + 1) * (s-1) * \pi^{-s/2} * \zeta(s).
    $$

    EXAMPLES:
        sage: zeta_symmetric(0.7)
        0.497580414651127
        sage: zeta_symmetric(1-0.7)
        0.497580414651127
        sage: RR = RealField(200)
        sage: zeta_symmetric(RR('0.7'))
        0.49758041465112690357779107525638385212657443284080589766062
        sage: I = CC.0
        sage: zeta_symmetric(RR('0.5') + I*RR('14.0'))
        0.000201294444235258 + 4.74338450462408e-20*I
        sage: zeta_symmetric(RR('0.5') + I*RR('14.1'))
        0.0000489893483255687 + 1.18584612615602e-20*I
        sage: zeta_symmetric(RR('0.5') + I*RR('14.2'))
        -0.0000868931282620101 - 2.03287907341032e-20*I

    REFERENCE:
      I copied the definition of xi from
        \url{http://www.math.ubc.ca/~pugh/RiemannZeta/RiemannZetaLong.html}
    """
    if not (is_ComplexNumber(s) or is_RealNumber(s)):
        s = RealField()(s)

    if s == 1:  # deal with poles, hopefully
        return s.parent()(1/2)

    R = s.parent()
    return (s/2 + 1).gamma()   *    (s-1)   * (R.pi()**(-s/2))  *  s.zeta()


##     # Use PARI on complex nubmer
##     prec = s.prec()
##     s = pari.new_with_bits_prec(s, prec)
##     pi = pari.pi()
##     w = (s/2 + 1).gamma() * (s-1) * pi **(-s/2) * s.zeta()
##     z = w._sage_()
##     if z.prec() < prec:
##         raise RuntimeError, "Error computing zeta_symmetric(%s) -- precision loss."%s
##     return z


#def pi_approx(prec=53):
#    """
#    Return pi computed to prec bits of precision.
#    """
#   return real_field.RealField(prec).pi()
