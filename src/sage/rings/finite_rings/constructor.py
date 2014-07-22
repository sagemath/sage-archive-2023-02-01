r"""
Finite Fields

Sage supports arithmetic in finite prime and extension fields.
Several implementation for prime fields are implemented natively in
Sage for several sizes of primes `p`. These implementations
are


-  ``sage.rings.finite_rings.integer_mod.IntegerMod_int``,

-  ``sage.rings.finite_rings.integer_mod.IntegerMod_int64``, and

-  ``sage.rings.finite_rings.integer_mod.IntegerMod_gmp``.


Small extension fields of cardinality `< 2^{16}` are
implemented using tables of Zech logs via the Givaro C++ library
(``sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro``).
While this representation is very fast it is limited to finite
fields of small cardinality. Larger finite extension fields of
order `q >= 2^{16}` are internally represented as
polynomials over smaller finite prime fields. If the
characteristic of such a field is 2 then NTL is used internally to
represent the field
(``sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e``).
In all other case the PARI C library is used
(``sage.rings.finite_rings.finite_field_ext_pari.FiniteField_ext_pari``).

However, this distinction is internal only and the user usually
does not have to worry about it because consistency across all
implementations is aimed for. In all extension field
implementations the user may either specify a minimal polynomial or
leave the choice to Sage.

For small finite fields the default choice are Conway polynomials.

The Conway polynomial `C_n` is the lexicographically first
monic irreducible, primitive polynomial of degree `n` over
`GF(p)` with the property that for a root `\alpha`
of `C_n` we have that
`\beta=
\alpha^{(p^n - 1)/(p^m - 1)}` is a root of
`C_m` for all `m` dividing `n`. Sage
contains a database of Conway polynomials which also can be queried
independently of finite field construction.

While Sage supports basic arithmetic in finite fields some more
advanced features for computing with finite fields are still not
implemented. For instance, Sage does not calculate embeddings of
finite fields yet.

EXAMPLES::

    sage: k = GF(5); type(k)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>

::

    sage: k = GF(5^2,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>

::

    sage: k = GF(2^16,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>

::

    sage: k = GF(3^16,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>

Finite Fields support iteration, starting with 0.

::

    sage: k = GF(9, 'a')
    sage: for i,x in enumerate(k):  print i,x
    0 0
    1 a
    2 a + 1
    3 2*a + 1
    4 2
    5 2*a
    6 2*a + 2
    7 a + 2
    8 1
    sage: for a in GF(5):
    ...    print a
    0
    1
    2
    3
    4

We output the base rings of several finite fields.

::

    sage: k = GF(3); type(k)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

::

    sage: k = GF(9,'alpha'); type(k)
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

::

    sage: k = GF(3^40,'b'); type(k)
    <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

Further examples::

    sage: GF(2).is_field()
    True
    sage: GF(next_prime(10^20)).is_field()
    True
    sage: GF(19^20,'a').is_field()
    True
    sage: GF(8,'a').is_field()
    True

AUTHORS:

- William Stein: initial version

- Robert Bradshaw: prime field implementation

- Martin Albrecht: Givaro and ntl.GF2E implementations
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

import random

from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.structure.parent_gens import normalize_names

import sage.rings.arith as arith
import sage.rings.integer as integer

import sage.rings.polynomial.polynomial_element as polynomial_element
import sage.rings.polynomial.multi_polynomial_element as multi_polynomial_element

# We don't late import this because this means trouble with the Givaro library
# On a Macbook Pro OSX 10.5.8, this manifests as a Bus Error on exiting Sage.
# TODO: figure out why
from finite_field_givaro import FiniteField_givaro

import sage.interfaces.gap

from sage.structure.factory import UniqueFactory

class FiniteFieldFactory(UniqueFactory):
    """
    Return the globally unique finite field of given order with
    generator labeled by the given name and possibly with given
    modulus.

    INPUT:

    - ``order`` -- a prime power

    - ``name`` -- string; must be specified unless ``order`` is prime.

    - ``modulus`` -- (optional) either a defining polynomial for the
      field, or a string specifying an algorithm to use to generate
      such a polynomial.  If ``modulus`` is a string, it is passed to
      :meth:`~sage.rings.polynomial.irreducible_element()` as the
      parameter ``algorithm``; see there for the permissible values of
      this parameter.

    - ``elem_cache`` -- cache all elements to avoid creation time
      (default: order < 500)

    - ``check_irreducible`` -- verify that the polynomial modulus is
      irreducible

    - ``proof`` -- bool (default: ``True``): if ``True``, use provable
      primality test; otherwise only use pseudoprimality test.

    - ``args`` -- additional parameters passed to finite field
      implementations

    - ``kwds`` -- additional keyword parameters passed to finite field
      implementations

    ALIAS: You can also use ``GF`` instead of ``FiniteField`` -- they
    are identical.

    EXAMPLES::

        sage: k.<a> = FiniteField(9); k
        Finite Field in a of size 3^2
        sage: parent(a)
        Finite Field in a of size 3^2
        sage: charpoly(a, 'y')
        y^2 + 2*y + 2

    We illustrate the proof flag.  The following example would hang
    for a very long time if we didn't use ``proof=False``.

    .. NOTE::

        Magma only supports ``proof=False`` for making finite fields,
        so falsely appears to be faster than Sage -- see :trac:10975.

    ::

        sage: k = FiniteField(10^1000 + 453, proof=False)
        sage: k = FiniteField((10^1000 + 453)^2, 'a', proof=False)      # long time -- about 5 seconds

    ::

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x +1 )
        sage: f = K.modulus(); f
        x^5 + 4*x + 1
        sage: type(f)
         <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>

    The modulus must be irreducible::

        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x)
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not.

    You can't accidentally fool the constructor into thinking the
    modulus is irreducible when it is not, since it actually tests
    irreducibility modulo `p`.  Also, the modulus has to be of the
    right degree::

        sage: F.<x> = QQ[]
        sage: factor(x^5 + 2)
        x^5 + 2
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 + 2)
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not.
        sage: K.<a> = GF(5**5, name='a', modulus=x^3 + 3*x + 3)
        Traceback (most recent call last):
        ...
        ValueError: The degree of the modulus does not correspond to the
        cardinality of the field.

    If you wish to live dangerously, you can tell the constructor not
    to test irreducibility using ``check_irreducible=False``, but this
    can easily lead to crashes and hangs -- so do not do it unless you
    know that the modulus really is irreducible and has the correct
    degree!

    ::

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**2, name='a', modulus=x^2 + 2, check_irreducible=False)

        sage: L = GF(3**2, name='a', modulus=QQ[x](x - 1), check_irreducible=False)
        sage: L.list()  # random
        [0, a, 1, 2, 1, 2, 1, 2, 1]

    The order of a finite field must be a prime power::

        sage: GF(1)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be > 1.
        sage: GF(100)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be a prime power.

    Finite fields with explicit random modulus are not cached::

        sage: k.<a> = GF(5**10, modulus='random')
        sage: n.<a> = GF(5**10, modulus='random')
        sage: n is k
        False
        sage: GF(5**10, 'a') is GF(5**10, 'a')
        True

    We check that various ways of creating the same finite field yield
    the same object, which is cached::

        sage: K = GF(7, 'a')
        sage: L = GF(7, 'b')
        sage: K is L
        True
        sage: K = GF(4,'a'); K.modulus()
        x^2 + x + 1
        sage: L = GF(4,'a', K.modulus())
        sage: K is L
        True
        sage: M = GF(4,'a', K.modulus().change_variable_name('y'))
        sage: K is M
        True

    You may print finite field elements as integers. This currently
    only works if the order of field is `<2^{16}`, though::

        sage: k.<a> = GF(2^8, repr='int')
        sage: a
        2

    The following demonstrate coercions for finite fields using Conway
    polynomials::

        sage: k = GF(5^2, conway=True, prefix='z'); a = k.gen()
        sage: l = GF(5^5, conway=True, prefix='z'); b = l.gen()
        sage: a + b
        3*z10^5 + z10^4 + z10^2 + 3*z10 + 1

    Note that embeddings are compatible in lattices of such finite
    fields::

        sage: m = GF(5^3, conway=True, prefix='z'); c = m.gen()
        sage: (a+b)+c == a+(b+c)
        True
        sage: (a*b)*c == a*(b*c)
        True
        sage: from sage.categories.pushout import pushout
        sage: n = pushout(k, l)
        sage: o = pushout(l, m)
        sage: q = pushout(n, o)
        sage: q(o(b)) == q(n(b))
        True

    Another check that embeddings are defined properly::

        sage: k = GF(3**10, conway=True, prefix='z')
        sage: l = GF(3**20, conway=True, prefix='z')
        sage: l(k.gen()**10) == l(k.gen())**10
        True
    """
    def create_key_and_extra_args(self, order, name=None, modulus=None, names=None,
                                  impl=None, proof=None, **kwds):
        """
        EXAMPLES::

            sage: GF.create_key_and_extra_args(9, 'a')
            ((9, ('a',), x^2 + 2*x + 2, None, '{}', 3, 2, True), {})
            sage: GF.create_key_and_extra_args(9, 'a', foo='value')
            ((9, ('a',), x^2 + 2*x + 2, None, "{'foo': 'value'}", 3, 2, True), {'foo': 'value'})
        """
        from sage.structure.proof.all import WithProof, arithmetic
        if proof is None: proof = arithmetic()
        with WithProof('arithmetic', proof):
            order = int(order)
            if order <= 1:
                raise ValueError("the order of a finite field must be > 1.")

            if arith.is_prime(order):
                name = None
                modulus = None
                p = integer.Integer(order)
                n = integer.Integer(1)
            elif arith.is_prime_power(order):
                if not names is None: name = names
                name = normalize_names(1,name)

                p, n = arith.factor(order)[0]

                # The following is a temporary solution that allows us
                # to construct compatible systems of finite fields
                # until algebraic closures of finite fields are
                # implemented in Sage.  It requires the user to
                # specify two parameters:
                #
                # - `conway` -- boolean; if True, this field is
                #   constructed to fit in a compatible system using
                #   a Conway polynomial.
                # - `prefix` -- a string used to generate names for
                #   automatically constructed finite fields
                #
                # See the docstring of FiniteFieldFactory for examples.
                #
                # Once algebraic closures of finite fields are
                # implemented, this syntax should be superseded by
                # something like the following:
                #
                #     sage: Fpbar = GF(5).algebraic_closure('z')
                #     sage: F, e = Fpbar.subfield(3)  # e is the embedding into Fpbar
                #     sage: F
                #     Finite field in z3 of size 5^3
                #
                # This temporary solution only uses actual Conway
                # polynomials (no pseudo-Conway polynomials), since
                # pseudo-Conway polynomials are not unique, and until
                # we have algebraic closures of finite fields, there
                # is no good place to store a specific choice of
                # pseudo-Conway polynomials.
                if name is None:
                    if not ('conway' in kwds and kwds['conway']):
                        raise ValueError("parameter 'conway' is required if no name given")
                    if 'prefix' not in kwds:
                        raise ValueError("parameter 'prefix' is required if no name given")
                    name = kwds['prefix'] + str(n)

                if 'conway' in kwds and kwds['conway']:
                    from conway_polynomials import conway_polynomial
                    if 'prefix' not in kwds:
                        raise ValueError("a prefix must be specified if conway=True")
                    if modulus is not None:
                        raise ValueError("no modulus may be specified if conway=True")
                    # The following raises a RuntimeError if no polynomial is found.
                    modulus = conway_polynomial(p, n)

                if modulus is None or isinstance(modulus, str):
                    # A string specifies an algorithm to find a suitable modulus.
                    if modulus == "default":    # for backward compatibility
                        modulus = None
                    modulus = GF(p)['x'].irreducible_element(n, algorithm=modulus)
                elif isinstance(modulus, (list, tuple)):
                    modulus = GF(p)['x'](modulus)
                elif sage.rings.polynomial.polynomial_element.is_Polynomial(modulus):
                    modulus = modulus.change_variable_name('x')
                else:
                    raise TypeError("wrong type for modulus parameter")
            else:
                raise ValueError("the order of a finite field must be a prime power.")

            return (order, name, modulus, impl, str(kwds), p, n, proof), kwds

    def create_object(self, version, key, check_irreducible=True, elem_cache=None,
                      names=None, **kwds):
        """
        EXAMPLES::

            sage: K = GF(19) # indirect doctest
            sage: TestSuite(K).run()
        """
        # IMPORTANT!  If you add a new class to the list of classes
        # that get cached by this factor object, then you *must* add
        # the following method to that class in order to fully support
        # pickling:
        #
        #     def __reduce__(self):   # and include good doctests, please!
        #         return self._factory_data[0].reduce_data(self)
        #
        # This is not in the base class for finite fields, since some finite
        # fields need not be created using this factory object, e.g., residue
        # class fields.

        if len(key) == 5:
            # for backward compatibility of pickles (see trac 10975).
            order, name, modulus, impl, _ = key
            p, n = arith.factor(order)[0]
            proof = True
        else:
            order, name, modulus, impl, _, p, n, proof = key

        if elem_cache is None:
            elem_cache = order < 500

        if n == 1 and (impl is None or impl == 'modn'):
            from finite_field_prime_modn import FiniteField_prime_modn
            # Using a check option here is probably a worthwhile
            # compromise since this constructor is simple and used a
            # huge amount.
            K = FiniteField_prime_modn(order, check=False)
        else:
            # We have to do this with block so that the finite field
            # constructors below will use the proof flag that was
            # passed in when checking for primality, factoring, etc.
            # Otherwise, we would have to complicate all of their
            # constructors with check options (like above).
            from sage.structure.proof.all import WithProof
            with WithProof('arithmetic', proof):
                if check_irreducible and polynomial_element.is_Polynomial(modulus):
                    if modulus.parent().base_ring().characteristic() == 0:
                        modulus = modulus.change_ring(FiniteField(p))
                    if not modulus.is_irreducible():
                        raise ValueError("finite field modulus must be irreducible but it is not.")
                    if modulus.degree() != n:
                        raise ValueError("The degree of the modulus does not correspond to the cardinality of the field.")
                if name is None:
                    raise TypeError("you must specify the generator name.")
                if impl is None:
                    if order < zech_log_bound:
                        # DO *NOT* use for prime subfield, since that would lead to
                        # a circular reference in the call to ParentWithGens in the
                        # __init__ method.
                        impl = 'givaro'
                    elif order % 2 == 0:
                        impl = 'ntl'
                    else:
                        impl = 'pari_ffelt'
                if impl == 'givaro':
                    if 'repr' in kwds:
                        repr = kwds['repr']
                    else:
                        repr = 'poly'
                    K = FiniteField_givaro(order, name, modulus, repr=repr, cache=elem_cache)
                elif impl == 'ntl':
                    from finite_field_ntl_gf2e import FiniteField_ntl_gf2e
                    K = FiniteField_ntl_gf2e(order, name, modulus)
                elif impl == 'pari_ffelt':
                    from finite_field_pari_ffelt import FiniteField_pari_ffelt
                    K = FiniteField_pari_ffelt(p, modulus, name)
                elif (impl == 'pari_mod'
                      or impl == 'pari'):    # for unpickling old pickles
                    from finite_field_ext_pari import FiniteField_ext_pari
                    K = FiniteField_ext_pari(order, name, modulus)
                else:
                    raise ValueError("no such finite field implementation: %s" % impl)

            # Temporary; see create_key_and_extra_args() above.
            if 'prefix' in kwds:
                K._prefix = kwds['prefix']

        return K

    def other_keys(self, key, K):
        """
        EXAMPLES::

            sage: key, extra = GF.create_key_and_extra_args(9, 'a'); key
            (9, ('a',), x^2 + 2*x + 2, None, '{}', 3, 2, True)
            sage: K = GF.create_object(0, key); K
            Finite Field in a of size 3^2
            sage: GF.other_keys(key, K)
            [(9, ('a',), x^2 + 2*x + 2, None, '{}', 3, 2, True),
             (9, ('a',), x^2 + 2*x + 2, 'givaro', '{}', 3, 2, True)]
        """
        if len(key) == 5: # backward compat
            order, name, modulus, impl, _ = key
            p, n = arith.factor(order)[0]
            proof = True
        else:
            order, name, modulus, impl, _, p, n, proof = key

        from sage.structure.proof.all import WithProof
        with WithProof('arithmetic', proof):
            if K.degree() > 1:
                modulus = K.modulus().change_variable_name('x')
            new_keys = [(order, name, modulus, impl, _, p, n, proof)]
            from finite_field_prime_modn import FiniteField_prime_modn
            if isinstance(K, FiniteField_prime_modn):
                impl = 'modn'
            elif isinstance(K, FiniteField_givaro):
                impl = 'givaro'
            else:
                from finite_field_ntl_gf2e import FiniteField_ntl_gf2e
                from finite_field_ext_pari import FiniteField_ext_pari
                from finite_field_pari_ffelt import FiniteField_pari_ffelt
                if isinstance(K, FiniteField_ntl_gf2e):
                    impl = 'ntl'
                elif isinstance(K, FiniteField_ext_pari):
                    impl = 'pari_mod'
                elif isinstance(K, FiniteField_pari_ffelt):
                    impl = 'pari_ffelt'
            new_keys.append( (order, name, modulus, impl, _, p, n, proof) )
            return new_keys


GF = FiniteField = FiniteFieldFactory("FiniteField")


def is_PrimeFiniteField(x):
    """
    Returns True if x is a prime finite field.

    EXAMPLES::

        sage: from sage.rings.finite_rings.constructor import is_PrimeFiniteField
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
    from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic

    return isinstance(x, FiniteField_prime_modn) or \
           (isinstance(x, FiniteField_generic) and x.degree() == 1)

zech_log_bound = 2**16
