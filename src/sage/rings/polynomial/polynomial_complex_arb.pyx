# -*- coding: utf-8
r"""
Univariate polynomials over `\CC` with interval coefficients using Arb.

This is a binding to the `Arb library <http://fredrikj.net/arb/>`_; it
may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

.. SEEALSO::

    - :mod:`Complex balls using Arb <sage.rings.complex_arb>`

TESTS:

    sage: type(polygen(ComplexBallField(140)))
    <type 'sage.rings.polynomial.polynomial_complex_arb.Polynomial_complex_arb'>
    sage: Pol.<x> = CBF[]
    sage: (x+1/2)^3
    x^3 + 1.500000000000000*x^2 + 0.7500000000000000*x + 0.1250000000000000
"""

# temporary base, will switch to Polynomial as soon as we have a different data
# structure
cdef class Polynomial_complex_arb(Polynomial_generic_dense):
    pass
