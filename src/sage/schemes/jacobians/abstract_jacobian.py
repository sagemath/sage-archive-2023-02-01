"""
Base class for Jacobians of curves
"""

#*******************************************************************************
#  Copyright (C) 2005 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.rings.all import is_Field
from sage.schemes.generic.scheme import Scheme, is_Scheme

def is_Jacobian(J):
    """
    Return True if `J` is of type Jacobian_generic.

    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian, is_Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: J = Jacobian(C)
        sage: is_Jacobian(J)
        True

    ::

        sage: E = EllipticCurve('37a1')
        sage: is_Jacobian(E)
        False
    """
    return isinstance(J, Jacobian_generic)

def Jacobian(C):
    """
    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: Jacobian(C)
        Jacobian of Projective Curve over Rational Field defined by x^3 + y^3 + z^3
    """
    try:
        return C.jacobian()
    except AttributeError:
        return Jacobian_generic(C)

class Jacobian_generic(Scheme):
    """
    Base class for Jacobians of projective curves.

    The input must be a projective curve over a field.

    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: J = Jacobian(C); J
        Jacobian of Projective Curve over Rational Field defined by x^3 + y^3 + z^3
    """
    def __init__(self, C):
        """
        TESTS::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve(x^3 + y^3 + z^3)
            sage: J = Jacobian_generic(C); J
            Jacobian of Projective Curve over Rational Field defined by x^3 + y^3 + z^3
            sage: type(J)
            <class 'sage.schemes.jacobians.abstract_jacobian.Jacobian_generic'>
            sage: J == loads(dumps(J))
            True

        ::

            sage: Jacobian_generic(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Argument (=Integer Ring) must be a scheme.
            sage: Jacobian_generic(P2)
            Traceback (most recent call last):
            ...
            ValueError: C (=Projective Space of dimension 2 over Rational Field) must have dimension 1.
            sage: P2.<x, y, z> = ProjectiveSpace(Zmod(6), 2)
            sage: C = Curve(x + y + z)
            sage: Jacobian_generic(C)
            Traceback (most recent call last):
            ...
            TypeError: C (=Projective Curve over Ring of integers modulo 6 defined by x + y + z) must be defined over a field.
        """
        if not is_Scheme(C):
            raise TypeError, "Argument (=%s) must be a scheme."%C
        if not is_Field(C.base_ring()):
            raise TypeError, "C (=%s) must be defined over a field."%C
        if C.dimension() != 1:
            raise ValueError, "C (=%s) must have dimension 1."%C
        self.__curve = C
        Scheme.__init__(self, C.base_scheme())

    def __cmp__(self, J):
        """
        Compare the Jacobian self to `J`.  If `J` is a Jacobian, then
        self and `J` are equal if and only if their curves are equal.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J1 = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J1 == J1
            True
            sage: J1 == P2
            False
            sage: J1 < P2
            False
            sage: J1 > P2
            True
            sage: J2 = Jacobian(Curve(x + y + z))
            sage: J1 == J2
            False
            sage: J1 > J2
            True
            sage: J1 < J2
            False
        """
        if not is_Jacobian(J):
            return cmp(type(self), type(J))
        return cmp(self.curve(), J.curve())

    def _repr_(self):
        """
        Return a string representation of this Jacobian.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3)); J
            Jacobian of Projective Curve over Rational Field defined by x^3 + y^3 + z^3
            sage: J._repr_()
            'Jacobian of Projective Curve over Rational Field defined by x^3 + y^3 + z^3'
        """
        return "Jacobian of %s"%self.__curve

    def _point_class(self):
        raise NotImplementedError

    def curve(self):
        """
        Return the curve of which self is the Jacobian.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J.curve()
            Projective Curve over Rational Field defined by x^3 + y^3 + z^3
        """
        return self.__curve
