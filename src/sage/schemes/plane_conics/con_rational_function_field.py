r"""
Projective plane conics over a rational function field

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

from sage.schemes.plane_conics.con_field import ProjectiveConic_field

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
        
    def has_rational_point(self, point = False, algorithm = 'default', read_cache = True):
        from sage.libs.pari.pari_instance import PariInstance
        from constructor import Conic
        
        pari = PariInstance()
        if read_cache:
            if self._rational_point is not None:
                if point:
                    return True, self._rational_point
                else:
                    return True
        else:
            self._rational_point = None
        
        if algorithm == 'default':
            if self.base_ring().characteristic() == 2:
                raise NotImplementedError("has_rational_point not implemented for field of characteristic 2.")
            new_conic = self.diagonalization()[0]
            coeff = new_conic.coefficients()
            if coeff[0] == 0:
                self.point([1,0,0])
            elif coeff[3] == 0:
                self.point([0,1,0])
            elif coeff[5] == 0:
                self.point([0,0,1])
            if self._rational_point is not None:
                if point:
                    return True, self._rational_point
                else:
                    return True
            
            (coeff, multipliers) = new_conic._reduce_conic()
            if coeff[0].degree() % 2 == coeff[1].degree() % 2 and coeff[1].degree() % 2 == coeff[2].degree() % 2:
                case = 0
            else:
                case = 1
            
            t, = self.base_ring().base().gens()
            supp = []
            roots = [[], [], []]
            for i in (0,1,2):
                supp.append(coeff[i].factor())
                for p in supp[i]:
                    if p[1] != 1 or p[0].leading_coefficient() != 1:
                        raise ValueError("Expected monic factor of degree 1.")
                    u = pari('u')
                    f = pari(coeff[2] * u**2 + coeff[0]).Mod(p[0])
                    factor = f.factor()
                    roots[i].append(PolynomialRing(self.base_ring().base().quotient(p[0]), 'u')(factor[0][0]))
            
            if case == 0:
                leading_conic = Conic(self.base_ring().base_ring(), [coeff[0].leading_coefficient(), coeff[1].leading_coefficient(), coeff[2].leading_coefficient()])
                pt = leading_conic.has_rational_point(True)
                if pt:
                    self._find_point(coeff, roots, pt[1])
                else:
                    return False
            
            if point:
                return (True, self._find_point())
            else:
                return True
        
        
    
    def _reduce_conic(self):
        from sage.rings.arith import lcm, gcd
        from sage.modules.free_module_element import vector
        from sage.rings.fraction_field import is_FractionField
        
        # start with removing fractions
        coeff = [self.coefficients()[0], self.coefficients()[3], self.coefficients()[5]]
        coeff = lcm(lcm(coeff[0].denominator(), coeff[1].denominator()), coeff[2].denominator()) * vector(coeff)
        # go to base ring of fraction field
        coeff = [x.numerator() for x in coeff]
        
        # remove common divisors
        labda = mu = nu = 1
        g1 = g2 = g3 = 0
        ca, cb, cc = coeff
        while g1 != 1 or g2 != 1 or g3 != 1:
            g1 = gcd(ca,cb); ca = ca/g1; cb = cb/g1; cc = cc*g1; nu = g1*nu
            g2 = gcd(ca,cc); ca = ca/g2; cc = cc/g2; cb = cb*g2; mu = g2*mu
            g3 = gcd(cb,cc); cb = cb/g3; cc = cc/g3; ca = ca*g3; labda = g3*labda
        coeff = [ca, cb, cc]
        multipliers = [labda, mu, nu]
        
        # remove squares
        for i, x in enumerate(coeff):
            if (is_FractionField(x.parent())):
                # go to base ring of fraction field
                x = x.numerator()
            
            decom = x.squarefree_decomposition()
            x1 = 1; x2 = 1
            for factor in decom:
                if factor[1] % 2 == 0:
                    x2 = x2 * factor[0] ** (factor[1] / 2)
                else:
                    x1 = x1 * factor[0] ** factor[1]
            if (len(decom) != 0):
                x = x1
                for j, y in enumerate(multipliers):
                    if j != i:
                        multipliers[j] = y * x2
            coeff[i] = self.base_ring().base().coerce(x);
        
        return (coeff, multipliers)
        
    def _find_point(self, coefficients, certificate, solution):
        #
        return
