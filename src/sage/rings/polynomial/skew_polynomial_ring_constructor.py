r"""
Constructor for skew polynomial rings

This module provides the function :func:`SkewPolynomialRing`, which constructs
rings of univariate skew polynomials, and implements caching to prevent the
same ring being created in memory multiple times (which is wasteful and breaks
the general assumption in Sage that parents are unique).

AUTHOR:

- Xavier Caruso (2012-06-29)
"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

import cysignals
from sage.structure.category_object import normalize_names
import sage.rings.ring as ring
from sage.categories.morphism import Morphism, IdentityMorphism

def SkewPolynomialRing(base_ring, base_ring_automorphism=None, names=None, sparse=False):
    r"""
    Return the globally unique skew polynomial ring with given
    properties and variable name.

    Given a ring `R` and a ring automorphism `\base_ring_automorphism` of `R`, the ring of
    skew polynomials `R[X,\base_ring_automorphism]` is the usual abelian group polynomial
    `R[X]` equipped with the modification multiplication deduced from the
    rule `X*a = \base_ring_automorphism(a)*X`.

    .. SEEALSO::

        Class ``SkewPolynomialRing_general`` in sage.rings.polynomial.skew_polynomial_ring.py
        Class ``SkewPolynomial`` in sage.rings.polynomial.skew_polynomial_element

    INPUT:

    - ``base_ring`` -- a commutative ring

    - ``base_ring_automorphism`` -- an automorphism of the base ring (also called
      twisting map)

    - ``names`` -- a string or a list of one string

    - ``sparse`` -- a boolean (default: ``False``)

    .. NOTE::

        The current implementation of skew polynomial rings does not support derivations
        and the ring thus created is actually a special case of skew polynomials where the
        derivation is taken to be zero. Such skew polynomials are called Linearized Polynomials.

        Sparse skew polynomials and multivariate skew polynomials are
        currently not implemented.

    OUTPUT:

    ``SkewPolynomialRing(base_ring, base_ring_automorphism, names, sparse=False)``
    returns a univariate skew polynomial ring over `base_ring` twisted by
    `\base_ring_automorphism`.

    UNIQUENESS and IMMUTABILITY:

    In Sage there is exactly one skew polynomial ring for each
    triple (base ring, twisting map, name of the variable).

        EXAMPLES of VARIABLE NAME CONTEXT::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R, sigma); S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

        The names of the variables defined above cannot be arbitrarily modified because
        each skew polynomial ring is unique in Sage and other objects in Sage could have
        pointers to that skew polynomial ring.

        However, the variable can be changed within the scope of a ``with`` block using
        the localvars context::

            sage: with localvars(S, ['y']):
            ....:     print S
            Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    SQUARE BRACKETS NOTATION:

    You can alternatively create a skew polynomial ring over `R`
    twisted by `\base_ring_automorphism` by writing ``R['varname', base_ring_automorphism]``.

    EXAMPLES:

    We first define the base ring::

        sage: R.<t> = ZZ[]; R
        Univariate Polynomial Ring in t over Integer Ring

    and the twisting map::

        sage: base_ring_automorphism = R.hom([t+1]); base_ring_automorphism
        Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
          Defn: t |--> t + 1

    Now, we are ready to define the skew polynomial ring::

        sage: S = SkewPolynomialRing(R, base_ring_automorphism, names='x'); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    Use the diamond brackets notation to make the variable ready
    for use after you define the ring::

        sage: S.<x> = SkewPolynomialRing(R, base_ring_automorphism)
        sage: (x + t)^2
        x^2 + (2*t + 1)*x + t^2

    Here is an example with the square bracket notations::

        sage: S.<x> = R['x', base_ring_automorphism]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    Rings with different variables names are different::

        sage: R['x', base_ring_automorphism] == R['y', base_ring_automorphism]
        False

    TESTS:

    You must specify a variable name::

        sage: SkewPolynomialRing(R, base_ring_automorphism)
        Traceback (most recent call last):
        ...
        TypeError: you must specify the name of the variable.

    With this syntax, it is not possible to omit the name of the
    variable neither in LHS nor in RHS. If we omit it in LHS, the
    variable is not created::

        sage: Sy = R['y', base_ring_automorphism]; Sy
        Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        sage: y.parent()
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined

    If we omit it in RHS, sage tries to create a polynomial ring and fails::

        sage: Sz.<z> = R[base_ring_automorphism]
        Traceback (most recent call last):
        ...
        ValueError: variable name 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring\n  Defn: t |--> t + 1' is not alphanumeric

    Multivariate skew polynomial rings are not supported::

        sage: S = SkewPolynomialRing(R, base_ring_automorphism,names=['x','y'])
        Traceback (most recent call last):
        ...
        NotImplementedError: multivariate skew polynomials rings not supported.

    Sparse skew polynomial rings are not implemented::

        sage: S = SkewPolynomialRing(R, base_ring_automorphism, names='x', sparse=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: sparse skew polynomial rings are not implemented.

    TODO::

    - Sparse Skew Polynomial Ring
    - Multivariate Skew Polynomial Ring
    - Add derivations.
    """

    if not isinstance(base_ring, ring.CommutativeRing):
        raise TypeError('base_ring must be a commutative ring')
    if base_ring_automorphism is None:
        base_ring_automorphism = IdentityMorphism(base_ring)
    else:
        if not isinstance(base_ring_automorphism,Morphism) or \
                base_ring_automorphism.domain() != base_ring or \
                base_ring_automorphism.codomain() != base_ring:
            raise TypeError("base_ring_automorphism must be a ring automorphism of base_ring (=%s)" % base_ring)
    if sparse:
        raise NotImplementedError("sparse skew polynomial rings are not implemented.")
    if names is None:
        raise TypeError("you must specify the name of the variable.")
    try:
        names = normalize_names(1, names)[0]
    except IndexError:
        raise NotImplementedError("multivariate skew polynomials rings not supported.")

    import sage.rings.polynomial.skew_polynomial_ring
    R = sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general(base_ring, base_ring_automorphism, names, sparse)

    return R
