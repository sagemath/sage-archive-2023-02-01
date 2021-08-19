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

    @staticmethod
    def _composition_impl(left, right):
        """
        Called by :meth:`_composition_`.
        """
        raise NotImplementedError('subclasses should implement _composition_impl')

    def _composition_(self, other, homset):
        """
        Return the composition of this elliptic-curve morphism
        with another elliptic-curve morphism.
        """
        if not isinstance(self, EllipticCurveHom) or not isinstance(other, EllipticCurveHom):
            raise TypeError(f'cannot compose {type(self)} with {type(other)}')

        try:
            return self._composition_impl(self, other)
        except NotImplementedError:
            pass

        try:
            return other._composition_impl(self, other)
        except NotImplementedError:
            pass

        # fall back to generic formal composite map
        return Morphism._composition_(self, other, homset)

