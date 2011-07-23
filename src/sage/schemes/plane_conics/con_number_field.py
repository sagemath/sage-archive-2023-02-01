r"""
Projective plane conics over a number field

AUTHORS:

- Marco Streng (2010-07-20)

"""
#*****************************************************************************
#       Copyright (C) 2009/2010 Marco Streng <marco.streng@gmail.com>
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

from sage.rings.all import (is_RationalField,
                            is_RingHomomorphism, is_RealIntervalField,
                            is_ComplexField, is_ComplexIntervalField,
                            RDF, CDF, AA, QQbar, PolynomialRing)

from con_field import ProjectiveConic_field

class ProjectiveConic_number_field(ProjectiveConic_field):
    r"""
    Create a projective plane conic curve over a number field.
    See ``Conic`` for full documentation.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3 - 2, 'a')
        sage: P.<X, Y, Z> = K[]
        sage: Conic(X^2 + Y^2 - a*Z^2)
        Projective Conic Curve over Number Field in a with defining polynomial x^3 - 2 defined by X^2 + Y^2 + (-a)*Z^2

    TESTS::

        sage: K.<a> = NumberField(x^3 - 3, 'a')
        sage: Conic([a, 1, -1])._test_pickling()
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        EXAMPLES ::

            sage: Conic([1, 1, 1])
            Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2
        """
        ProjectiveConic_field.__init__(self, A, f)

        # a single prime such that self has no point over the completion
        self._local_obstruction = None
        # all finite primes such that self has no point over the completion
        self._finite_obstructions = None
        # all infinite primes such that self has no point over the completion
        self._infinite_obstructions = None



