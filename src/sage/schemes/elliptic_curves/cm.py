"""
Complex multiplication for elliptic curves

This module implements the functions

- ``hilbert_class_polynomial``
- ``cm_j_invariants``
- ``cm_j_invariants_and_orders``

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

from sage.interfaces.all import magma
from sage.rings.all import (Integer,
                            RationalField,
                            IntegerRing,
                            is_fundamental_discriminant,
                            PolynomialRing)

def hilbert_class_polynomial(D, algorithm=None):
    r"""
    Returns the Hilbert class polynomial for discriminant `D`.

    INPUT:

    - ``D`` (int) -- a negative integer congruent to 0 or 1 modulo 4.

    - ``algorithm`` (string, default None) -- if "sage" then use the Sage implementation; if "magma" then call Magma (if available).

    OUTPUT:

    (integer polynomial) The Hilbert class polynomial for the
    discriminant `D`.

    ALGORITHM:

    - If ``algorithm`` = "sage": Use complex approximations to the roots.

    - If ``algorithm`` = "magma": Call the appropriate Magma function.

    AUTHORS:

    - Sage implementation originally by Eduardo Ocampo Alvarez and
      AndreyTimofeev

    - Sage implementation corrected by John Cremona (using corrected precision bounds from Andreas Enge)

    - Magma implementation by David Kohel

    EXAMPLES::

        sage: hilbert_class_polynomial(-4)
        x - 1728
        sage: hilbert_class_polynomial(-7)
        x + 3375
        sage: hilbert_class_polynomial(-23)
        x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
        sage: hilbert_class_polynomial(-37*4)
        x^2 - 39660183801072000*x - 7898242515936467904000000
        sage: hilbert_class_polynomial(-37*4, algorithm="magma") # optional - magma
        x^2 - 39660183801072000*x - 7898242515936467904000000
        sage: hilbert_class_polynomial(-163)
        x + 262537412640768000
        sage: hilbert_class_polynomial(-163, algorithm="magma") # optional - magma
        x + 262537412640768000

    """
    if algorithm is None:
        algorithm = "sage"

    D = Integer(D)
    if D >= 0:
        raise ValueError, "D (=%s) must be negative"%D
    if not (D%4 in [0,1]):
         raise ValueError, "D (=%s) must be a discriminant"%D

    if algorithm == "magma":
        magma.eval("R<x> := PolynomialRing(IntegerRing())")
        f = str(magma.eval("HilbertClassPolynomial(%s)"%D))
        return IntegerRing()['x'](f)

    if algorithm != "sage":
        raise ValueError, "%s is not a valid algorithm"%algorithm

    from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
    from sage.rings.all import RR, ZZ, ComplexField
    from sage.functions.all import elliptic_j

    # get all primitive reduced quadratic forms, (necessary to exclude
    # imprimitive forms when D is not a fundamental discriminant):

    rqf = BinaryQF_reduced_representatives(D, primitive_only=True)

    # compute needed precision
    #
    # NB: [http://arxiv.org/abs/0802.0979v1], quoting Enge (2006), is
    # incorrect.  Enge writes (2009-04-20 email to John Cremona) "The
    # source is my paper on class polynomials
    # [http://hal.inria.fr/inria-00001040] It was pointed out to me by
    # the referee after ANTS that the constant given there was
    # wrong. The final version contains a corrected constant on p.7
    # which is consistent with your example. It says:

    # "The logarithm of the absolute value of the coefficient in front
    # of X^j is bounded above by
    #
    # log (2*k_2) * h + pi * sqrt(|D|) * sum (1/A_i)
    #
    # independently of j", where k_2 \approx 10.163.

    h = len(rqf) # class number
    c1 = 3.05682737291380 # log(2*10.63)
    c2 = sum([1/RR(qf.a) for qf in rqf], RR(0))
    prec =  c2*RR(3.142)*RR(D).abs().sqrt() + h*c1  # bound on log
    prec = prec * 1.45   # bound on log_2 (1/log(2) = 1.44..)
    prec = 10 + prec.ceil()  # allow for rounding error
#    print "prec = ",prec

    # set appropriate precision for further computing

    Dsqrt = D.sqrt(prec=prec)
    R = ComplexField(prec)['t']
    t = R.gen()
    pol = R(1)
    for qf in rqf:
        a, b, c = list(qf)
        tau = (b+Dsqrt)/(a<<1)
        pol *=  (t - elliptic_j(tau))

    coeffs = [cof.real().round() for cof in pol.coeffs()]
    return IntegerRing()['x'](coeffs)


def cm_j_invariants(K):
    r"""
    Return a list of all CM `j`-invariants in the field `K`.

    INPUT:

    - ``K`` -- a number field (currently only implemented for `K=\QQ`)

    OUTPUT:

    (list) -- A list of CM `j`-invariants in the field `K`.

    .. note::

       This is currently only implemented for the rationals.  David
       Kohel has large tables for other fields, but they are not in
       Sage yet.

    EXAMPLE::

        sage: cm_j_invariants(QQ)
        [0, 54000, -12288000, 1728, 287496, -3375, 16581375, 8000, -32768, -884736, -884736000, -147197952000, -262537412640768000]

    ::

        sage: cm_j_invariants(QuadraticField(-1, 'i'))
        Traceback (most recent call last):
        ...
        NotImplementedError: Enumeration of CM j-invariants over Number Field in i with defining polynomial x^2 + 1 not yet implemented

    """
    if K == RationalField():
        return [Integer(x) for x in [0, 54000, -12288000, 1728, \
                               287496, -3375, 16581375, 8000, \
                               -32768,  -884736, -884736000,\
                               -147197952000, -262537412640768000]]
    else:
        raise NotImplementedError, "Enumeration of CM j-invariants over %s not yet implemented"%K


def cm_j_invariants_and_orders(K):
    r"""
    Return a list of all CM `j`-invariants in the field `K`, together with the associated orders.

    INPUT:

    - ``K`` -- a number field (currently only implemented for `K=\QQ`)

    OUTPUT:

    (list) A list of 3-tuples `(D,f,j)` where `j` is a CM
    `j`-invariant with quadratic fundamental discriminant `D` and
    conductor `f`

    EXAMPLE::

        sage: cm_j_invariants_and_orders(QQ)
        [(-3, 3, -12288000), (-3, 2, 54000), (-3, 1, 0), (-4, 2, 287496), (-4, 1, 1728), (-7, 2, 16581375), (-7, 1, -3375), (-8, 1, 8000), (-11, 1, -32768), (-19, 1, -884736), (-43, 1, -884736000), (-67, 1, -147197952000), (-163, 1, -262537412640768000)]

    ::

        sage: cm_j_invariants_and_orders(QuadraticField(-1, 'i'))
        Traceback (most recent call last):
        ...
        NotImplementedError: Enumeration of CM j-invariants over Number Field in i with defining polynomial x^2 + 1 not yet implemented

    """
    if K == RationalField():
        T = [ (0,-3, 1), (54000,-3,2), (-12288000, -3,3), (1728,-1, 1), \
               (287496,-1, 2), (-3375,-7,1), (16581375, -7, 2), (8000,-2,1), \
               (-32768, -11, 1),  (-884736, -19,1), (-884736000,-43,1),\
               (-147197952000, -67,1), (-262537412640768000,-163,1)
               ]
        dis = lambda D:  Integer(D) if D%4 in [0,1] else Integer(4*D)
        T = [(dis(D),Integer(f),Integer(j)) for (j,D,f) in T]
        T.sort()
        T.reverse()
        return T
    else:
        raise NotImplementedError, "Enumeration of CM j-invariants over %s not yet implemented"%K

