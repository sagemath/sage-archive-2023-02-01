r"""
Projective plane conics over a field

AUTHORS:

- Lennart Ackermans (2015-08-06)

"""
#*****************************************************************************
#       Copyright (C) 2015 Lennart Ackermans <lennart@ackermans.info>
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

from sage.rings.all import PolynomialRing

from con_field import ProjectiveConic_field

class ProjectiveConic_rational_function_field(ProjectiveConic_field):
    r"""
    Create a projective plane conic curve over a rational function field.
    See ``Conic`` for full documentation.

    EXAMPLES::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: P.<X, Y, Z> = K[]
        sage: Conic(X^2 + Y^2 - Z^2)
        Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Rational Field defined by X^2 + Y^2 - Z^2

    TESTS::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: Conic([K(1), 1, -1])._test_pickling()
    """
    
    def __init__(self, A, f):
        ProjectiveConic_field.__init__(self, A, f)
