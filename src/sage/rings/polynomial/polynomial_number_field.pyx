r"""
Univariate polynomials over number fields

AUTHOR:

- Luis Felipe Tabera Alonso (2014-02): initial version.

EXAMPLES:

Define a polynomial over an absolute number field and compute basic operations with them::

    sage: N.<a> = NumberField(x^2-2)
    sage: K.<x> = N[]
    sage: f = x - a
    sage: g = x^3 - 2*a + 1
    sage: f*(x + a)
    x^2 - 2
    sage: f + g
    x^3 + x - 3*a + 1
    sage: g // f
    x^2 + a*x + 2
    sage: g % f
    1
    sage: factor(x^3 - 2*a*x^2 - 2*x + 4*a)
    (x - 2*a) * (x - a) * (x + a)
    sage: gcd(f, x - a)
    x - a

Polynomials are aware of embeddings of the underlying field::

    sage: x = var('x')
    sage: Q7 = Qp(7)
    sage: r1 = Q7(3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8 +\
             6*7^9 + 6*7^10 + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16 + 7^17 +\
             4*7^18 + 6*7^19)
    sage: N.<b> = NumberField(x^2-2, embedding = r1)
    sage: K.<t> = N[]
    sage: f = t^3-2*t+1
    sage: f(r1)
    1 + O(7^20)

We can also construct polynomials over relative number fields::

    sage: N.<i, s2> = QQ[I, sqrt(2)]
    sage: K.<x> = N[]
    sage: f = x - s2
    sage: g = x^3 - 2*i*x^2 + s2*x
    sage: f*(x + s2)
    x^2 - 2
    sage: f + g
    x^3 - 2*I*x^2 + (sqrt2 + 1)*x - sqrt2
    sage: g // f
    x^2 + (-2*I + sqrt2)*x - 2*sqrt2*I + sqrt2 + 2
    sage: g % f
    -4*I + 2*sqrt2 + 2
    sage: factor(i*x^4 - 2*i*x^2 + 9*i)
    (I) * (x - I + sqrt2) * (x + I - sqrt2) * (x - I - sqrt2) * (x + I + sqrt2)
    sage: gcd(f, x-i)
    1
"""

#*****************************************************************************
#       Copyright (C) 2013 Luis Felipe Tabera Alonso <taberalf@unican.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from polynomial_element_generic import Polynomial_generic_dense_field

class Polynomial_absolute_number_field_dense(Polynomial_generic_dense_field):
    """
    Class of dense univariate polynomials over an absolute number field.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new polynomial of the polynomial ring ``parent``.

        INPUT:

             - ``parent`` -- the underlying Polynomial ring.

             - ``x`` -- (default: None) An object representing the polynomial.
               e.g. a list of coefficients. See
               :meth:`sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field.__init__` for more details.

             - ``check`` -- boolean (default: True) if True make sure that the
               coefficients of the polynomial are the ring of coefficients.

             - ``is_gen`` -- boolean (defaul: False) A boolean, True if `x` is
               the disinguished generator of the polynomial ring.

             - ``construct`` -- (default: False) A boolean, unused.

        EXAMPLES::

            sage: f = QQ[I][x].random_element()
            sage: type(f)
            <class 'sage.rings.polynomial.polynomial_number_field.Polynomial_absolute_number_field_dense'>
            sage: a = QQ[I][x](x)
            sage: a.is_gen()
            True
        """
        Polynomial_generic_dense_field.__init__(self, parent, x, check, is_gen, construct)


class Polynomial_relative_number_field_dense(Polynomial_generic_dense_field):
    """
    Class of dense univariate polynomials over a relative number field.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new polynomial of the polynomial ring ``parent``.

        INPUT:

             - ``parent`` -- the underlying Polynomial ring.

             - ``x`` -- (default: None) An object representing the polynomial.
               e.g. a list of coefficients. See
               :meth:`sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field.__init__` for more details.

             - ``check`` -- boolean (default: True) if True make sure that the
               coefficients of the polynomial are the ring of coefficients.

             - ``is_gen`` -- boolean (defaul: False) A boolean, True if ``x`` is
               the disinguished generator of the polynomial ring.

             - ``construct`` -- (default: False) A boolean, unused.

        EXAMPLES::

            sage: f = NumberField([x^2-2, x^2-3], 'a')[x].random_element()
            sage: type(f)
            <class 'sage.rings.polynomial.polynomial_number_field.Polynomial_relative_number_field_dense'>
        """
        Polynomial_generic_dense_field.__init__(self, parent, x, check, is_gen, construct)
