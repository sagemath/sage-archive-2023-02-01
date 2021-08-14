"""
Elliptic-curve morphisms

This class serves as a common parent for various specializations
of morphisms between elliptic curves. In the future, some of the
code currently contained in the EllipticCurveIsogeny class should
be moved here, and new code for cases not currently covered by the
EllipticCurveIsogeny class should inherit from this class in order
to eventually provide a uniform interface for all elliptic-curve
maps — regardless of differences in internal representations.
"""

from sage.categories.morphism import Morphism

class EllipticCurveHom(Morphism):

    def _repr_type(self):           # used by Morphism._repr_
        return 'Elliptic-curve'

