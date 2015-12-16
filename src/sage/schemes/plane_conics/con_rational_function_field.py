# -*- coding: utf-8 -*-
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
    Create a projective plane conic curve over a rational function field
    `F(t)`, where `F` is any field.
    See :func:`Conic` for full documentation.

    The algorithms used in this class come mostly from [HC2006]_.

    EXAMPLES::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: P.<X, Y, Z> = K[]
        sage: Conic(X^2 + Y^2 - Z^2)
        Projective Conic Curve over Fraction Field of Univariate
        Polynomial Ring in t over Rational Field defined by
        X^2 + Y^2 - Z^2

    TESTS::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: Conic([K(1), 1, -1])._test_pickling()

    REFERENCES:

    .. [HC2006] Mark van Hoeij and John Cremona, Solving Conics over
        function fields. J. Théor. Nombres Bordeaux, 2006.
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        EXAMPLES::

            sage: c = Conic([1, 1, 1]); c
            Projective Conic Curve over Rational Field defined by
            x^2 + y^2 + z^2
        """
        ProjectiveConic_field.__init__(self, A, f)
    
    def has_rational_point(self, point = False, algorithm = 'default',
        read_cache = True):
        r"""
        Returns True if and only if the conic ``self``
        has a point over its base field `F(t)`, which is a field of rational
        functions.

        If ``point`` is True, then returns a second output, which is
        a rational point if one exists.

        Points are cached whenever they are found. Cached information
        is used if and only if ``read_cache`` is True.
        
        The default algorithm does not (yet) work for all base fields `F`.
        In particular, sage is required to have:
        
        * an algorithm for finding the square root of elements in finite
          extensions of `F`;
        
        * a factorization algorithm for `F[t]`, as well as a square-free
          decomposition and gcd function;
        
        * an algorithm for solving conics over `F`.
        
        ALGORITHM:
        
        The parameter ``algorithm`` specifies the algorithm
        to be used:

        * ``'default'`` -- use a native Sage implementation, based on the
           algorithm Conic in [HC2006]_.

        * ``'magma'`` (requires Magma to be installed) --
          delegates the task to the Magma computer algebra
          system.
        
        EXAMPLES:
        
        We can find points for function fields over `\QQ` and finite
        fields::
        
            sage: K.<t> = FractionField(PolynomialRing(QQ, 't'))
            sage: C = Conic(K, [t^2-2, 2*t^3, -2*t^3-13*t^2-2*t+18])
            sage: C.has_rational_point(point=True)
            (True, (3 : (1/3*t + 1/3)/(1/3*t) : 1))
            sage: R.<t> = FiniteField(23)[]
            sage: C = Conic([2, t^2+1, t^2+5])
            sage: C.has_rational_point()
            True
            sage: C.has_rational_point(point=True)
            (True, (5*t : 8 : 1))
        
        Due to limitations in Sage of algorithms we depend on, it is not
        yet possible to find points on conics over multivariate function fields
        (see the requirements above)::
        
            sage: F.<t1> = FractionField(QQ['t1'])
            sage: K.<t2> = FractionField(F['t2'])
            sage: a = K(1)
            sage: b = 2*t2^2+2*t1*t2-t1^2
            sage: c = -3*t2^4-4*t1*t2^3+8*t1^2*t2^2+16*t1^3-t2-48*t1^4
            sage: C = Conic([a,b,c])
            sage: C.has_rational_point()
            ...
            Traceback (most recent call last):
            ...
            NotImplementedError: is_square() not implemented for elements of
            Univariate Quotient Polynomial Ring in tbar over Fraction Field
            of Univariate Polynomial Ring in t1 over Rational Field with
            modulus tbar^2 + t1*tbar - 1/2*t1^2
        
        In some cases, the algorithm requires us to be
        able to solve conics over `F`. In particular, the following does not
        work, and will still not work after square-free decomposition has
        been implemented for this field::

            sage: P.<u> = QQ[]
            sage: E = P.fraction_field()
            sage: Q.<Y> = E[]
            sage: F.<v> = E.extension(Y^2 - u^3 - 1)
            sage: R.<t> = F[]
            sage: K = R.fraction_field()
            sage: C = Conic(K, [u, v, 1])
            sage: C.has_rational_point()
            ...
            Traceback (most recent call last):
            ...
            NotImplementedError: square-free decomposition not implemented
            for this polynomial
            
        """
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
                raise NotImplementedError("has_rational_point not implemented \
for field of characteristic 2.")
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
            
            # We save the coefficients of the reduced form in coeff
            # A zero of the reduced conic can be multiplied by multipliers
            # to get a zero of the old conic
            (coeff, multipliers) = new_conic._reduce_conic()
            if coeff[0].degree() % 2 == coeff[1].degree() % 2 and \
                    coeff[1].degree() % 2 == coeff[2].degree() % 2:
                case = 0
            else:
                case = 1
            
            t, = self.base_ring().base().gens() # t in F[t]
            supp = []
            roots = [[], [], []]
            # loop through the coefficients and find a root of f_i (as in
            # [HC2006]) modulo each element in the coefficients support
            for i in (0,1,2):
                supp.append(list(coeff[i].factor()))
                for p in supp[i]:
                    if p[1] != 1:
                        raise ValueError("Expected factor of degree 1.")
                    if p[0].leading_coefficient() != 1:
                        x = p[0]/list(p[0])[-1]
                    else:
                        x = p[0]
                    N = p[0].base_ring().extension(x, 'tbar')
                    R = PolynomialRing(N, 'u')
                    u, = R.gens()
                    # If p[0] has degree 1, sage might forget the "defining
                    # polynomial" of N, so we define our own modulo operation
                    if p[0].degree() == 1:
                        mod = t.parent().hom([-x+t])
                    else:
                        mod = N
                    if i == 0:
                        x = -mod(coeff[2])/mod(coeff[1])
                    elif i == 1:
                        x = -mod(coeff[0])/mod(coeff[2])
                    else:
                        x = -mod(coeff[1])/mod(coeff[0])
                    if x.is_square():
                        root = N(x.sqrt())
                    else:
                        if point:
                            return False, False
                        else:
                            return False
                    roots[i].append(root)
            
            if case == 0:
            # Find a solution of (5) in [HC2006]
                leading_conic = Conic(self.base_ring().base_ring(),
                            [coeff[0].leading_coefficient(),
                            coeff[1].leading_coefficient(),
                            coeff[2].leading_coefficient()])
                has_point = leading_conic.has_rational_point(True)
                if has_point[0]:
                    if point:
                        pt = _find_point(coeff, roots, supp, has_point[1])
                        pt = self.point([pt[0] * multipliers[0],
                        pt[1] * multipliers[1], pt[2] * multipliers[2]])
                        return (True, pt)
                    else:
                        return True
                else:
                    if point:
                        return False, False
                    else:
                        return False
            
            if point:
                pt = _find_point(coeff, roots, supp)
                pt = self.point([pt[0] * multipliers[0],
                pt[1] * multipliers[1], pt[2] * multipliers[2]])
                return (True, pt)
            else:
                return True
        else:
            return ProjectiveConic_field.has_rational_point(self, point,
                algorithm, read_cache)
        
    
    def _reduce_conic(self):
        r"""
        Return the reduced form of the conic, i.e. a conic with base field
        `K=F(t)` and coefficients `a,b,c` such that `a,b,c \in F[t]`,
        `\gcd(a,b)=\gcd(b,c)=\gcd(c,a)=1` and `abc` is square-free.
        
        OUTPUT:
        
        A tuple (coefficients, multipliers), the coefficients of the conic
        in reduced form and multipliers `\lambda, \mu, \nu \in F(t)^*` such
        that `(x,y,z) \in F(t)` is a solution of the reduced conic if and only
        if `(\lambda x, \mu y, \nu z)` is a solution of the original conic.
        
        ALGORITMH:
        
        The algorithm used is the algorithm ReduceConic in [HC2006]_.
        
        EXAMPLES::
        
            sage: K.<t> = FractionField(PolynomialRing(QQ, 't'))
            sage: C = Conic(K, [t^2-2, 2*t^3, -2*t^3-13*t^2-2*t+18])
            sage: C._reduce_conic()
            ([t^2 - 2, 2*t, -2*t^3 - 13*t^2 - 2*t + 18], [t, 1, t])
        """
        from sage.rings.arith import lcm, gcd
        from sage.modules.free_module_element import vector
        from sage.rings.fraction_field import is_FractionField
        
        # start with removing fractions
        coeff = [self.coefficients()[0], self.coefficients()[3],
                self.coefficients()[5]]
        coeff = lcm(lcm(coeff[0].denominator(), coeff[1].denominator()),
                coeff[2].denominator()) * vector(coeff)
        # go to base ring of fraction field
        coeff = [x.numerator() for x in coeff]
        
        # remove common divisors
        labda = mu = nu = 1
        g1 = g2 = g3 = 0
        ca, cb, cc = coeff
        while g1 != 1 or g2 != 1 or g3 != 1:
            g1 = gcd(ca,cb); ca = ca/g1; cb = cb/g1; cc = cc*g1; nu = g1*nu
            g2 = gcd(ca,cc); ca = ca/g2; cc = cc/g2; cb = cb*g2; mu = g2*mu
            g3 = gcd(cb,cc); cb = cb/g3; cc = cc/g3; ca = ca*g3;
            labda = g3*labda
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
                    if factor[1] % 2 == 0:
                        x2 = x2 * factor[0] ** (factor[1] / 2)
                    else:
                        x = x * factor[0]
                        x2 = x2 * factor[0] ** ((factor[1]-1) / 2)
                else:
                    x = x * factor[0]
            for j, y in enumerate(multipliers):
                if j != i:
                    multipliers[j] = y * x2
            coeff[i] = self.base_ring().base().coerce(x);
        
        return (coeff, multipliers)
    
def _find_point(coefficients, roots, supports, solution = 0):
    r"""
    Given a solubility certificate like in [HC2006]_, find a point on
    the conic given by ``coefficients``.
    
    INPUT:
    
    - ``coefficients`` -- list of three defining coefficients of a
    conic in reduced form.
    - ``roots`` -- roots of a solubility certificate like in [HC2006]_
    - ``supports`` -- 3-tuple where supports[i] is a list of all monic
    irreducible `p \in F[t]` that divide the `i`'th coefficient.
    - ``solution`` -- (default: 0) a solution of (5) in [HC2006]_, if
    one is required, 0 otherwise.
    
    OUTPUT:
    
    A point `(x,y,z) \in F(t)` of the conic given by ``coefficients``.
    
    ALGORITMH:
    
    The algorithm used is the algorithm FindPoint in [HC2006]_.
    
    EXAMPLES::
        
        sage: K.<t> = FractionField(QQ['t'])
        sage: C = Conic(K, [t^2-2, 2*t^3, -2*t^3-13*t^2-2*t+18])
        sage: C.has_rational_point(point=True)
        (True, (3 : (1/3*t + 1/3)/(1/3*t) : 1))
    """
    from sage.matrix.constructor import matrix
    Ft = coefficients[0].parent()
    deg = [coefficients[0].degree(), coefficients[1].degree(),
            coefficients[2].degree()]
    if deg[0] % 2 == deg[1] % 2 and deg[1] % 2 == deg[2] % 2:
        case = 0
        for (n, support) in enumerate(supports):
            for (i, factor) in enumerate(support):
                if factor[0].degree() == 1:
                    case = 1
                    support.pop(i)
                    roots[n].pop(i)
                    break
                if case == 1:
                    break
    else:
        case = 1
    # definitions as in [HC2006]
    A = max(0, ((deg[1] + deg[2]) / 2).round('up') - case)
    B = max(0, ((deg[2] + deg[0]) / 2).round('up') - case)
    C = max(0, ((deg[0] + deg[1]) / 2).round('up') - case)
    # generate the names of all variables we require
    var_names = list(Ft.variable_names()) + ['X%d' % i for i in range(A+1)] + \
                ['Y%d' % i for i in range(B+1)] + \
                ['Z%d' % i for i in range(C+1)] + ['W']
    R = PolynomialRing(coefficients[0].parent().base_ring(), var_names)
    var_names = R.gens()
    XX = var_names[1:A+2]
    YY = var_names[A+2:A+B+3]
    ZZ = var_names[A+B+3:A+B+C+4]
    t = var_names[0]
    X = sum([XX[n]*t**n for n in range(A+1)])
    Y = sum([YY[n]*t**n for n in range(B+1)])
    Z = sum([ZZ[n]*t**n for n in range(C+1)])
    E = [] # list that will contain linear polynomials (set E in [HC2006])
    
    if case == 0:
        # We need an extra variable and linear equation
        W = var_names[A+B+C+4]
        (x,y,z) = solution
        E += [XX[A] - x*W, YY[B] - y*W, ZZ[C] - z*W]
    # For all roots as calculated by has_rational_point(), do the
    # calculations from step (6) of FindPoint in [HC2006]
    for (i, p) in enumerate(supports[0]):
        # lift to F[t] and map to R, with R as defined above
        if roots[0][i].parent().is_finite():
            root = roots[0][i].polynomial()
        else:
            root = roots[0][i].lift()
        alpha = root.parent().hom([t])(root)
        r = (Y - alpha*Z).quo_rem(p[0])[1]
        for i in range(r.degree(t) + 1):
            E.append(r.coefficient({t:i}))
    for (i, p) in enumerate(supports[1]):
        if roots[1][i].parent().is_finite():
            root = roots[1][i].polynomial()
        else:
            root = roots[1][i].lift()
        alpha = root.parent().hom([t])(root)
        r = (Z - alpha*X).quo_rem(p[0])[1]
        for i in range(r.degree(t) + 1):
            E.append(r.coefficient({t:i}))
    for (i, p) in enumerate(supports[2]):
        if roots[2][i].parent().is_finite():
            root = roots[2][i].polynomial()
        else:
            root = roots[2][i].lift()
        alpha = root.parent().hom([t])(root)
        r = (X - alpha*Y).quo_rem(p[0])[1]
        for i in range(r.degree(t) + 1):
            E.append(r.coefficient({t:i}))
    # Create a matrix to compute a solution of the system obtained by
    # equating all elements of E to 0
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
