r"""
Constructors for sjew polynomial rings

This module provides the function :func:`SkewPolynomialRing`, which constructs
rings of univariate skew polynomials, and implements caching to prevent the
same ring being created in memory multiple times (which iswasteful and breaks
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

from sage.structure.parent_gens import normalize_names
from sage.structure.element import is_Element
import sage.rings.ring as ring

from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.categories.morphism import Morphism, IdentityMorphism
from sage.rings.morphism import RingHomomorphism

def SkewPolynomialRing(base_ring,sigma=None,name=None,names=None,sparse=False):
    r"""
    Return the globally unique skew polynomial ring with given
    properties and variable name.

    INPUT:

    - ``base_ring`` -- a commutative ring

    - ``sigma`` -- an endomorphism of the base ring (so-called
      twisting map)

    - ``name`` -- a string

    - ``names`` -- a string or a list of one string

    - ``sparse`` -- a boolean (default: False)

    .. NOTE::

        The input ``names`` is redundant with ``name`` but useful
        to allow the syntax ``S.<x> = SkewPolynomialRing(R,sigma)``.

        Sparse skew polynomials are currently not implemented.

    OUTPUT:

    The univariate skew polynomial ring over ``base_ring``
    twisted by `\sigma` (i.e. the ring `R` of usual polynomial
    ring over ``base_ring`` with the modified multiplication
    given by the rule `X * a = \sigma(a)` for all scalar `a`
    in the base ring).

    UNIQUENESS and IMMUTABILITY:

    In Sage there is exactly one skew polynomial ring for each
    triple (base ring, twisting map, name of the variable).

    SQUARE BRACKETS NOTATION:

    You can alternatively create a skew polynomial ring over `R`
    twisted by `\sigma` by writing ``R['varname',sigma]``.

    EXAMPLES::

    We first define the base ring::

        sage: R.<t> = ZZ[]; R
        Univariate Polynomial Ring in t over Integer Ring

    and the twisting map::

        sage: sigma = R.hom([t+1]); sigma
        Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
          Defn: t |--> t + 1

    Now, we are ready to define the skew polynomial ring::

        sage: S = SkewPolynomialRing(R,sigma,name='x'); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    Use the diamond brackets notation to make the variable ready
    for use after you define the ring::

        sage: S.<x> = SkewPolynomialRing(R,sigma)
        sage: (x + t)^2
        x^2 + (2*t + 1)*x + t^2

    You must specify a variable name::

        sage: SkewPolynomialRing(R,sigma)
        Traceback (most recent call last):
        ...
        TypeError: You must specify the name of the variable.

    Here is an example with the square bracket notations::

        sage: S.<x> = R['x',sigma]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    With this syntax, it is not possible to omit the name of the
    variable neither in LHS nor in RHS. If we omit it in LHS, the
    variable is not created::

        sage: Sy = R['y',sigma]; Sy
        Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        sage: y.parent()
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined

    If we omit it in RHS, sage tries to create a polynomial ring and fails::

        sage: Sz.<z> = R[sigma]
        Traceback (most recent call last):
        ...
        ValueError: variable names must be alphanumeric, but one is 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
          Defn: t |--> t + 1' which is not.


    Rings with different variables names are different::

        sage: R['x',sigma] == R['y',sigma]
        False
    """
    import sage.rings.polynomial.skew_polynomial_ring as m

    if not isinstance(base_ring, ring.CommutativeRing):
        raise TypeError('base_ring must be a commutative ring')
    if sigma is None:
        sigma = IdentityMorphism(base_ring)
        # should one return a polynomial ring in that case?
    else:
        if not isinstance(sigma,Morphism) or sigma.domain() != base_ring or sigma.codomain() != base_ring:
            raise TypeError("sigma must be a ring endomorphism of base_ring (=%s)" % base_ring)

    if name is None:
        name = names
    if name is None:
        raise TypeError("You must specify the name of the variable.")

    R = None
    try:
        name = normalize_names(1, name)[0]
    except IndexError:
        raise NotImplementedError("Multivariate skew polynomials rings not supported.")

    if is_FiniteField(base_ring):
        R = m.SkewPolynomialRing_finite_field(base_ring,sigma,name,sparse)
    else:
        R = m.SkewPolynomialRing_general(base_ring,sigma,name,sparse)

    return R
