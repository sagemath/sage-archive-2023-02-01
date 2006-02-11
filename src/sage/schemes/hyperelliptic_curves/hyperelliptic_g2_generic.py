#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import hyperelliptic_generic
import jacobian_g2

from sage.schemes.generic.projective_space import ProjectiveSpace

class HyperellipticCurve_g2_generic(hyperelliptic_generic.HyperellipticCurve_generic):
    def igusa_invariants(self):
        raise NotImplementedError

    def is_odd_degree(self):
        """
        Returns True if the curve is an odd degree model.
        """
        f, h = self.hyperelliptic_polynomials()
        df = f.degree()
        if h.degree() < 3:
            return df%2 == 1
        elif df < 6:
            return False
        else:
            a0 = f.leading_coefficient()
            c0 = h.leading_coefficient()
            return (c0^2 + 4*a0) == 0

    def jacobian(self):
        return jacobian_g2.HyperellipticJacobian_g2(self)

    def kummer_morphism(self):
        """
        Returns the morphism of an odd degree hyperelliptic curve to the Kummer
        surface of its Jacobian.  (This could be extended to an even degree model
        if a prescribed embedding in its Jacobian is fixed.)
        """
        try:
            return self._kummer_morphism
        except AttributeError:
            pass
        if not self.is_odd_degree():
            raise TypeError, \
                  "Kummer embedding not determined for even degree model curves."
        J = self.jacobian()
        K = J.kummer_surface()
        return self._kummer_morphism
