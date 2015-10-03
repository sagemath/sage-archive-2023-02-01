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

from sage.rings.all import PolynomialRing, NumberField

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
        from constructor import Conic
        
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
                supp.append(list(coeff[i].factor()))
                for p in supp[i]:
                    if p[1] != 1 or p[0].leading_coefficient() != 1:
                        raise ValueError("Expected monic factor of degree 1.")
                    N = NumberField(p[0], 'tbar')
                    R = PolynomialRing(N, 'u')
                    u, = R.gens()
                    if i == 0:
                        f = N(coeff[1]) * u**2 + N(coeff[2])
                    elif i == 1:
                        f = N(coeff[2]) * u**2 + N(coeff[0])
                    else:
                        f = N(coeff[0]) * u**2 + N(coeff[1])
                    if f.is_irreducible():
                        return False
                    roots[i].append(f.roots()[0][0])
            import pdb;pdb.set_trace()
            if case == 0:
                leading_conic = Conic(self.base_ring().base_ring(), [coeff[0].leading_coefficient(), coeff[1].leading_coefficient(), coeff[2].leading_coefficient()])
                has_point = leading_conic.has_rational_point(True)
                if has_point[0]:
                    if point:
                        pt = self._find_point(coeff, roots, supp, has_point[1])
                        pt = self.point([pt[0] * multipliers[0], pt[1] * multipliers[1], pt[2] * multipliers[2]])
                        return (True, pt)
                    else:
                        return True
                else:
                    return False
            
            if point:
                pt = self._find_point(coeff, roots, supp)
                pt = self.point([pt[0] * multipliers[0], pt[1] * multipliers[1], pt[2] * multipliers[2]])
                return (True, pt)
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
            x = decom.unit(); x2 = 1
            for factor in decom:
                if factor[1] > 1:
                    x2 = x2 * factor[0] ** (factor[1] - 1)
                x = x * factor[0]
            for j, y in enumerate(multipliers):
                if j != i:
                    multipliers[j] = y * x2
            coeff[i] = self.base_ring().base().coerce(x);
        
        return (coeff, multipliers)
        
    def _find_point(self, coefficients, roots, supports, solution = 0):
        from sage.matrix.constructor import matrix
        Ft = self.base_ring().base()
        deg = [coefficients[0].degree(), coefficients[1].degree(), coefficients[2].degree()]
        if deg[0] % 2 == deg[1] % 2 and deg[1] % 2 == deg[2] % 2:
            case = 0
            for (n, support) in enumerate(supports):
                for (i, factor) in enumerate(support):
                    if factor[1] == 1:
                        case = 1
                        support.pop(i)
                        roots[n].pop(i)
                        break
                if case == 1:
                    break
        else:
            case = 1
        A = ((deg[1] + deg[2]) / 2).round('up') - case
        B = ((deg[2] + deg[0]) / 2).round('up') - case
        C = ((deg[0] + deg[1]) / 2).round('up') - case
        var_names = [Ft.gens()[0]] + ['X%d' % i for i in range(A+1)] + ['Y%d' % i for i in range(B+1)] + ['Z%d' % i for i in range(C+1)] + ['W']
        R = PolynomialRing(self.base_ring().base_ring(), A+B+C+5, var_names)
        var_names = R.gens()
        XX = var_names[1:A+2]
        YY = var_names[A+2:A+B+3]
        ZZ = var_names[A+B+3:A+B+C+4]
        t = var_names[0]
        X = sum([XX[n]*t**n for n in range(A+1)])
        Y = sum([YY[n]*t**n for n in range(B+1)])
        Z = sum([ZZ[n]*t**n for n in range(C+1)])
        E = [] # list that will contain linear polynomials (set E from the article)
        
        if case == 0:
            W = var_names[A+B+C+4]
            (x,y,z) = solution
            E += [XX[A] - x*W, YY[B] - y*W, ZZ[C] - z*W]
        for (i, p) in enumerate(supports[0]):
            alpha = roots[0][i].lift().parent().hom([t])(roots[0][i].lift())
            r = (Y - alpha*Z).quo_rem(p[0])[1]
            for i in range(r.degree(t) + 1):
                E.append(r.coefficient({t:i}))
        for (i, p) in enumerate(supports[1]):
            alpha = roots[1][i].lift().parent().hom([t])(roots[1][i].lift())
            r = (Z - alpha*X).quo_rem(p[0])[1]
            for i in range(r.degree(t) + 1):
                E.append(r.coefficient({t:i}))
        for (i, p) in enumerate(supports[2]):
            alpha = roots[2][i].lift().parent().hom([t])(roots[2][i].lift())
            r = (X - alpha*Y).quo_rem(p[0])[1]
            for i in range(r.degree(t) + 1):
                E.append(r.coefficient({t:i}))
        E2 = []
        for f in E:
            column = []
            for i in range(1,A+B+C+5):
                column.append(f.coefficient(var_names[i]))
            E2.append(column)
        M = matrix(Ft.base(), E2)
        solution_space = M.right_kernel()
        for v in solution_space.basis():
            if v[:A+B+C+3] != 0:
                X = X.subs({XX[i]:v[i] for i in range(A+1)})
                Y = Y.subs({YY[i]:v[A+1+i] for i in range(B+1)})
                Z = Z.subs({ZZ[i]:v[A+B+2+i] for i in range(C+1)})
                return [X,Y,Z]
        
        raise RuntimeError("No solution has been found though a solution exists.")
