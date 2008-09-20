r"""
Finite Fields

\SAGE supports arithmetic in finite prime and extension fields.
Several implementation for prime fields are implemented natively in
\SAGE for several sizes of primes $p$. These implementations are
\begin{itemize}
\item \code{sage.rings.integer_mod.IntegerMod_int},
\item \code{sage.rings.integer_mod.IntegerMod_int64}, and
\item \code{sage.rings.integer_mod.IntegerMod_gmp}.
\end{itemize}
Small extension fields
of cardinality $< 2^{16}$ are implemented using tables of Zech logs
via the Givaro C++ library
(\code{sage.rings.finite_field_givaro.FiniteField_givaro}). While this
representation is very fast it is limited to finite fields of small
cardinality. Larger finite extension fields of order $q >= 2^{16}$ are
internally represented as polynomials over a smaller finite prime
fields. If the characteristic of such a field is 2 then NTL is used
internally to represent the field
(\code{sage.rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e}). In all
other case the PARI C library is used
(\code{sage.rings.finite_field_ext_pari.FiniteField_ext_pari}).

However, this distinction is internal only and the user usually does
not have to worry about it because consistency across all
implementations is aimed for. In all extension field implementations
the user may either specify a minimal polynomial or leave the choice
to \SAGE.

For small finite fields the default choice are Conway polynomials.

The Conway polynomial $C_n$ is the lexicographically first monic
irreducible, primitive polynomial of degree $n$ over $GF(p)$ with the
property that for a root $\alpha$ of $C_n$ we have that $\beta=
\alpha^{(p^n - 1)/(p^m - 1)}$ is a root of $C_m$ for all $m$ dividing
$n$. \SAGE contains a database of conway polynomials which also can be
queried independendtly of finite field construction.

While \SAGE supports basic arithmetic in finite fields some more
advanced features for computing with finite fields are still not
implemented. For instance, \SAGE does not calculate embeddings of
finite fields yet.

EXAMPLES:
    sage: k = GF(5); type(k)
    <class 'sage.rings.finite_field_prime_modn.FiniteField_prime_modn'>

    sage: k = GF(5^2,'c'); type(k)
    <type 'sage.rings.finite_field_givaro.FiniteField_givaro'>

    sage: k = GF(2^16,'c'); type(k)
    <type 'sage.rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e'>

    sage: k = GF(3^16,'c'); type(k)
    <class 'sage.rings.finite_field_ext_pari.FiniteField_ext_pari'>

Finite Fields support iteration, starting with 0.

    sage: k = GF(9, 'a')
    sage: for i,x in enumerate(k):  print i,x
    0 0
    1 2*a
    2 a + 1
    3 a + 2
    4 2
    5 a
    6 2*a + 2
    7 2*a + 1
    8 1
    sage: for a in GF(5):
    ...    print a
    0
    1
    2
    3
    4

We output the base rings of several finite fields.

    sage: k = GF(3); type(k)
    <class 'sage.rings.finite_field_prime_modn.FiniteField_prime_modn'>
    sage: k.base_ring()
    Finite Field of size 3

    sage: k = GF(9,'alpha'); type(k)
    <type 'sage.rings.finite_field_givaro.FiniteField_givaro'>
    sage: k.base_ring()
    Finite Field of size 3

    sage: k = GF(3^40,'b'); type(k)
    <class 'sage.rings.finite_field_ext_pari.FiniteField_ext_pari'>
    sage: k.base_ring()
    Finite Field of size 3

Further examples:
    sage: GF(2).is_field()
    True
    sage: GF(next_prime(10^20)).is_field()
    True
    sage: GF(19^20,'a').is_field()
    True
    sage: GF(8,'a').is_field()
    True

AUTHORS:
     -- William Stein: initial version
     -- Robert Bradshaw: prime field implementation
     -- Martin Albrecht: Givaro and ntl.GF2E implementations
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
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

import weakref

from ring import is_FiniteField
from sage.structure.parent_gens import normalize_names

import arith
import integer

import polynomial.polynomial_element as polynomial_element
import polynomial.multi_polynomial_element as multi_polynomial_element

# We don't late import this because this means trouble with the Givaro library
# TODO: figure out why
from finite_field_givaro import FiniteField_givaro

import sage.interfaces.gap
import sage.databases.conway

cache = {}

def FiniteField(order, name=None, modulus=None, names=None,
                elem_cache=False, check_irreducible=True, *args, **kwds):
    """
    Return the globally unique finite field of given order with generator
    labeled by the given name and possibly with given modulus.

    INPUT:
        order --   int
        name --    string; must be specified if not a prime field
        modulus -- (optional) defining polynomial for field, i.e.,
                   generator of the field will be a root of this
                   polynomial; if not specified the choice of
                   definining polynomials can be arbitrary.
        elem_cache -- cache all elements to avoid creation time  (default: order<500)
        check_irreducible -- verify that the polynomial modulus is irreducible
        args -- additional parameters passed to finite field implementations
        kwds -- additional keyword parameters passed to finite field implementations

    ALIAS:
        You can also use GF instead of FiniteField -- they are identical.

    EXAMPLES:
        sage: k.<a> = FiniteField(9); k
        Finite Field in a of size 3^2
        sage: parent(a)
        Finite Field in a of size 3^2
        sage: charpoly(a, 'y')
        y^2 + 2*y + 2

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x +1 )
        sage: f = K.modulus(); f
        x^5 + 4*x + 1
        sage: type(f)
        <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_mod_p'>

    The modulus must be irreducible:
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x )
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not

    You can't accidently fool the constructor into thinking the
    modulus is irreducible when it isn't mod p, since it actually
    tests irreducibility modulo p.

        sage: F.<x> = QQ[]
        sage: factor(x^5+2)
        x^5 + 2
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 + 2 )
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not

    If you wish to live dangerously, you can tell the constructor not
    to test irreducibility using check_irreducible=False, but this
    can easily lead to crashes and hangs -- so do not do it unless
    you know that the modulus really is irreducible!

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**2, name='a', modulus=x^2 + 2, check_irreducible=False)

    For example, you may print finite field elements as integers. This
    currently only works if the order of field is $<2^{16}$, though.

        sage: k.<a> = GF(2^8,repr='int')
        sage: a
        2

    The order of a finite field must be a prime power:

        sage: GF(100)
        Traceback (most recent call last):
        ...
        ValueError: order of finite field must be a prime power

    Finite fields with random modulus are not cached:
        sage: k.<a> = GF(2^17,modulus='random')
        sage: n.<a> = GF(2^17,modulus='random')
        sage: n is k
        False
    """
    if not names is None: name = names
    order = int(order)
    name = normalize_names(1,name)

    if elem_cache is None:
        elem_cache = order < 500

    key = (order, name, modulus, str([args, kwds]))
    if modulus != 'random' and cache.has_key(key):
        K = cache[key]()
        if not K is None:
            return K
    if arith.is_prime(order):
        from finite_field_prime_modn import FiniteField_prime_modn
        K = FiniteField_prime_modn(order,*args,**kwds)
    else:
        if not arith.is_prime_power(order):
            raise ValueError, "order of finite field must be a prime power"
        if check_irreducible and polynomial_element.is_Polynomial(modulus):
            if modulus.parent().base_ring().characteristic() == 0:
                p = arith.factor(order)[0][0]
                modulus = modulus.change_ring(FiniteField(p))
            if not modulus.is_irreducible():
                raise ValueError, "finite field modulus must be irreducible but it is not"
        if name is None:
            raise TypeError, "you must specify the generator name"
        if order < zech_log_bound:
            # DO *NOT* use for prime subfield, since that would lead to
            # a circular reference in the call to ParentWithGens in the
            # __init__ method.
            K = FiniteField_givaro(order, name, modulus, cache=elem_cache, *args,**kwds)
        else:
            if integer.Integer(order).factor()[0][0] == 2:
                from finite_field_ntl_gf2e import FiniteField_ntl_gf2e
                K = FiniteField_ntl_gf2e(order, name, modulus, *args, **kwds)
            else:
                from finite_field_ext_pari import FiniteField_ext_pari
                K = FiniteField_ext_pari(order, name, modulus, *args, **kwds)

    if modulus != 'random':
        cache[key] = weakref.ref(K)
    return K


def is_PrimeFiniteField(x):
    """
    Returns True if x is a prime finite field.

    EXAMPLES:
        sage: is_PrimeFiniteField(QQ)
        False
        sage: is_PrimeFiniteField(GF(7))
        True
        sage: is_PrimeFiniteField(GF(7^2,'a'))
        False
        sage: is_PrimeFiniteField(GF(next_prime(10^90,proof=False)))
        True
    """
    from finite_field_prime_modn import FiniteField_prime_modn
    from ring import FiniteField as FiniteField_generic

    return isinstance(x, FiniteField_prime_modn) or \
           (isinstance(x, FiniteField_generic) and x.degree() == 1)

##################################################################

def conway_polynomial(p, n):
    r"""
    Return the Conway polynomial of degree n over GF(p), which is
    loaded from a table.

    If the requested polynomial is not known, this function raises a
    RuntimeError exception.

    INPUT:
        p -- int
        n -- int

    OUTPUT:
        Polynomial -- a polynomial over the prime finite field GF(p).

    NOTE: The first time this function is called a table is read from
    disk, which takes a fraction of a second.  Subsequent calls do not
    require reloading the table.

    See also the \code{ConwayPolynomials()} object, which is a table of
    Conway polynomials.   For example, if \code{c=ConwayPolynomials}, then
    \code{c.primes()} is a list of all primes for which the polynomials are
    known, and for a given prime $p$,  \code{c.degree(p)} is a list of all
    degrees for which the Conway polynomials are known.

    EXAMPLES:
        sage: conway_polynomial(2,5)
        x^5 + x^2 + 1
        sage: conway_polynomial(101,5)
        x^5 + 2*x + 99
        sage: conway_polynomial(97,101)
        Traceback (most recent call last):
        ...
        RuntimeError: requested conway polynomial not in database.
    """
    (p,n)=(int(p),int(n))
    R = FiniteField(p)['x']
    try:
        return R(sage.databases.conway.ConwayPolynomials()[p][n])
    except KeyError:
        raise RuntimeError, "requested conway polynomial not in database."

def exists_conway_polynomial(p, n):
    r"""
    Return True if the Conway polynomial over $F_p$ of degree $n$ is in the
    database and False otherwise.

    If the Conway polynomial is in the database, to obtain it use the
    command \code{conway_polynomial(p,n)}.

    EXAMPLES:
        sage: exists_conway_polynomial(2,3)
        True
        sage: exists_conway_polynomial(2,-1)
        False
        sage: exists_conway_polynomial(97,200)
        False
        sage: exists_conway_polynomial(6,6)
        False
    """
    return sage.databases.conway.ConwayPolynomials().has_polynomial(p,n)


zech_log_bound = 2**16
