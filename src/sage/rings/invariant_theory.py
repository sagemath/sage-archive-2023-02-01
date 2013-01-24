r"""
Classical Invariant Theory

This module lists classical invariants and covariants of homogeneous
polynomials (also called algebraic forms) under the action of the
special linear group. That is, we are dealing with polynomials of
degree `d` in `n` variables. The special linear group `SL(n,\CC)` acts
on the variables `(x_1,\dots, x_n)` linearly,

.. math::

    (x_1,\dots, x_n)^t \to A (x_1,\dots, x_n)^t
    ,\qquad
    A \in SL(n,\CC)

The linear action on the variables transforms a polynomial `p`
generally into a different polynomial `gp`. We can think of it as an
action on the space of coefficients in `p`. An invariant is a
polynomial in the coefficients that is invariant under this action. A
covariant is a polynomial in the coefficients and the variables
`(x_1,\dots, x_n)` that is invariant under the combined action.

For example, the binary quadric `p(x,y) = a x^2 + b x y + c y^2` has
as its invariant the discriminant `\mathop{disc}(p) = b^2 - 4 a
c`. This means that for any `SL(2,\CC)` coordinate change

.. math::

    \begin{pmatrix} x' \\ y' \end{pmatrix}
    =
    \begin{pmatrix} \alpha & \beta \\ \gamma & \delta \end{pmatrix}
    \begin{pmatrix} x \\ y \end{pmatrix}
    \qquad
    \alpha\delta-\beta\gamma=1

the discriminant is invariant, `\mathop{disc}\big(p(x',y')\big) =
\mathop{disc}\big(p(x,y)\big)`.

To use this module, you should use the factory object
:class:`invariant_theory <InvariantTheoryFactory>`. For example, take
the quartic::

    sage: R.<x,y> = QQ[]
    sage: q = x^4 + y^4
    sage: quartic = invariant_theory.binary_quartic(q);  quartic
    binary quartic with coefficients (1, 0, 0, 0, 1)


One invariant of a quartic is known as the Eisenstein
D-invariant. Since it is an invariant, it is a polynomial in the
coefficients (which are integers in this example)::

    sage: quartic.EisensteinD()
    1

One example of a covariant of a quartic is the so-called g-covariant
(actually, the Hessian). As with all covariants, it is a polynomial in
`x`, `y` and the coefficients::

    sage: quartic.g_covariant()
    -x^2*y^2

As usual, use tab completion and the online help to discover the
implemented invariants and covariants.

In general, the variables of the defining polynomial cannot be
guessed. For example, the zero polynomial can be thought of as a
homogeneous polynomial of any degree. Also, since we also want to
allow polynomial coefficients we cannot just take all variables of the
polynomial ring as the variables of the form. This is why you will
have to specify the variables explicitly if there is any potential
ambiguity. For example::

    sage: invariant_theory.binary_quartic(R.zero(), [x,y])
    binary quartic with coefficients (0, 0, 0, 0, 0)

    sage: invariant_theory.binary_quartic(x^4, [x,y])
    binary quartic with coefficients (0, 0, 0, 0, 1)

    sage: R.<x,y,t> = QQ[]
    sage: invariant_theory.binary_quartic(x^4 + y^4 + t*x^2*y^2, [x,y])
    binary quartic with coefficients (1, 0, t, 0, 1)

Finally, it is often convenient to use inhomogeneous polynomials where
it is understood that one wants to homogenize them. This is also
supported, just define the form with an inhomogeneous polynomial and
specify one less variable::

    sage: R.<x,t> = QQ[]
    sage: invariant_theory.binary_quartic(x^4 + 1 + t*x^2, [x])
    binary quartic with coefficients (1, 0, t, 0, 1)

REFERENCES:

..  [WpInvariantTheory]
    http://en.wikipedia.org/wiki/Glossary_of_invariant_theory
"""

########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.rings.all import QQ
from sage.misc.functional import is_odd
from sage.matrix.constructor import matrix
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method



######################################################################
def _guess_variables(polynomial, *args):
    """
    Return the polynomial variables

    INPUT:

    - ``polynomial`` -- a polynomial.

    - ``*args`` -- the variables. If none are specified, all variables
      in ``polynomial`` are returned. If a list or tuple is passed,
      the content is returned. If multiple arguments are passed, they
      are returned.

    OUTPUT:

    A tuple of variables in the parent ring of the polynoimal.

    EXAMPLES::

       sage: from sage.rings.invariant_theory import _guess_variables
       sage: R.<x,y> = QQ[]
       sage: _guess_variables(x^2+y^2)
       (x, y)
       sage: _guess_variables(x^2+y^2, x)
       (x,)
       sage: _guess_variables(x^2+y^2, x,y)
       (x, y)
       sage: _guess_variables(x^2+y^2, [x,y])
       (x, y)
    """
    if len(args)==0 or (len(args)==1 and args[0] is None):
        return polynomial.variables()
    elif len(args) == 1 and isinstance(args[0], (tuple, list)):
        return tuple(args[0])
    else:
        return tuple(args)


######################################################################

class AlgebraicForm(SageObject):
    """
    The base class of algebraic forms (i.e. homogeneous polynomials)

    You should only instantiate the derived classes of this base
    class.

    Derived classes must implement ``coeffs()`` and
    ``scaled_coeffs()``

    INPUT:

    - ``n`` -- The number of variables.

    - ``d`` -- The degree of the polynomial.

    - ``polynomial`` -- The polynomial.

    - ``*args`` -- The variables, as a single list/tuple, multiple
      arguments, or ``None`` to use all variables of the polynomial.

    EXAMPLES::

        sage: from sage.rings.invariant_theory import AlgebraicForm
        sage: R.<x,y> = QQ[]
        sage: p = x^2 + y^2
        sage: AlgebraicForm(2, 2, p).variables()
        (x, y)
        sage: AlgebraicForm(2, 2, p, None).variables()
        (x, y)
        sage: AlgebraicForm(3, 2, p).variables()
        (x, y, None)
        sage: AlgebraicForm(3, 2, p, None).variables()
        (x, y, None)

        sage: from sage.rings.invariant_theory import AlgebraicForm
        sage: R.<x,y,s,t> = QQ[]
        sage: p = s*x^2 + t*y^2
        sage: AlgebraicForm(2, 2, p, [x,y]).variables()
        (x, y)
        sage: AlgebraicForm(2, 2, p, x,y).variables()
        (x, y)

        sage: AlgebraicForm(3, 2, p, [x,y,None]).variables()
        (x, y, None)
        sage: AlgebraicForm(3, 2, p, x,y,None).variables()
        (x, y, None)

        sage: AlgebraicForm(2, 1, p, [x,y]).variables()
        Traceback (most recent call last):
        ...
        ValueError: polynomial is of the wrong degree

        sage: AlgebraicForm(2, 2, x^2+y, [x,y]).variables()
        Traceback (most recent call last):
        ...
        ValueError: polynomial is not homogeneous
    """

    def __init__(self, n, d, polynomial, *args, **kwds):
        """
        The Python constructor

        INPUT:

        See the class documentation.

        TESTS::

            sage: from sage.rings.invariant_theory import AlgebraicForm
            sage: R.<x,y> = QQ[]
            sage: form = AlgebraicForm(2, 2, x^2 + y^2)
        """
        self._n = n
        self._d = d
        self._polynomial = polynomial
        variables = _guess_variables(polynomial, *args)
        if len(variables) == n:
            pass
        elif len(variables) == n-1:
            variables = variables + (None,)
        else:
            raise ValueError('need '+str(n)+' or '+
                             str(n-1)+' variables, got '+str(variables))
        self._variables = variables
        self._ring = polynomial.parent()
        self._homogeneous = self._variables[-1] is not None
        self._check()


    def _check(self):
        """
        Check that the input is of the correct degree and number of
        variables

        EXAMPLES::

            sage: from sage.rings.invariant_theory import AlgebraicForm
            sage: R.<x,y,t> = QQ[]
            sage: p = x^2 + y^2
            sage: inv = AlgebraicForm(3, 2, p, [x,y,None])
            sage: inv._check()
        """
        degrees = set()
        R = self._ring
        for e in self._polynomial.exponents():
            deg = sum([ e[R.gens().index(x)]
                        for x in self._variables if x is not None ])
            degrees.add(deg)
        if self._homogeneous and len(degrees)>1:
            raise ValueError('polynomial is not homogeneous')
        if degrees == set() or \
                (self._homogeneous and degrees == set([self._d])) or \
                (not self._homogeneous and max(degrees) <= self._d):
            return
        else:
            raise ValueError('polynomial is of the wrong degree')


    def _check_covariant(self, method_name, g=None, invariant=False):
        """
        Test whether ``method_name`` actually returns a covariant

        INPUT:

        - ``method_name`` -- string. The name of the method that
          returns the invariant / covariant to test.

        - ``g`` -- a `SL(n,\CC)` matrix or ``None`` (default). The
          test will be to check that the covariant transforms
          corrently under this special linear group element acting on
          the homogeneous variables. If ``None``, a random matrix will
          be picked.

        - ``invariant`` -- boolean. Whether to additionaly test that
          it is an invariant.

        EXAMPLES::

            sage: R.<a0, a1, a2, a3, a4, x0, x1> = QQ[]
            sage: p = a0*x1^4 + a1*x1^3*x0 + a2*x1^2*x0^2 + a3*x1*x0^3 + a4*x0^4
            sage: quartic = invariant_theory.binary_quartic(p, x0, x1)

            sage: quartic._check_covariant('EisensteinE', invariant=True)
            sage: quartic._check_covariant('h_covariant')

            sage: quartic._check_covariant('h_covariant', invariant=True)
            Traceback (most recent call last):
            ...
            AssertionError: Not invariant
        """
        F = self._ring.base_ring()
        from sage.matrix.constructor import vector, random_matrix
        g = random_matrix(F, self._n, algorithm='unimodular')
        assert self._homogeneous
        v = vector(self.variables())
        g_v = g*v
        transform = dict( (v[i], g_v[i]) for i in range(self._n) )
        # The covariant of the transformed polynomial
        g_self = self.__class__(self.form().subs(transform), self.variables())
        cov_g = getattr(g_self, method_name)()
        # The transform of the covariant
        g_cov = getattr(self, method_name)().subs(transform)
        # they must be the same
        assert (g_cov - cov_g).is_zero(),  'Not covariant'
        if invariant:
            cov = getattr(self, method_name)()
            assert (cov - cov_g).is_zero(), 'Not invariant'


    def __cmp__(self, other):
        """
        Compare ``self`` with ``other``

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: quartic = invariant_theory.binary_quartic(x^4+y^4)
            sage: cmp(quartic, 'foo') == 0
            False
            sage: cmp(quartic, quartic)
            0
            sage: quartic.__cmp__(quartic)
            0
        """
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        return cmp(self.coeffs(), other.coeffs())


    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: quartic = invariant_theory.binary_quartic(x^4+y^4)
            sage: quartic._repr_()
            'binary quartic with coefficients (1, 0, 0, 0, 1)'
        """
        s = ''
        ary = ['unary', 'binary', 'ternary', 'quaternary', 'quinary',
               'senary', 'septenary', 'octonary', 'nonary', 'denary']
        try:
            s += ary[self._n-1]
        except IndexError:
            s += 'algebraic'
        ic = ['monic', 'quadratic', 'cubic', 'quartic', 'quintic',
              'sextic', 'septimic', 'octavic', 'nonic', 'decimic',
              'undecimic', 'duodecimic']
        s += ' '
        try:
            s += ic[self._d-1]
        except IndexError:
            s += 'form'
        s += ' with coefficients ' + str(self.coeffs())
        return s


    def form(self):
        """
        Return the defining polynomial

        OUTPUT:

        The polynomial used to define the algebraic form.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: quartic = invariant_theory.binary_quartic(x^4+y^4)
            sage: quartic.form()
            x^4 + y^4
        """
        return self._polynomial


    def variables(self):
        """
        Return the variables of the form

        OUTPUT:

        A tuple of variables. If inhomogeneous notation is use for the
        defining polynomial then the last entry will be ``None``.

        EXAMPLES::

            sage: R.<x,y,t> = QQ[]
            sage: quartic = invariant_theory.binary_quartic(x^4+y^4+t*x^2*y^2, [x,y])
            sage: quartic.variables()
            (x, y)

            sage: R.<x,y,t> = QQ[]
            sage: quartic = invariant_theory.binary_quartic(x^4+1+t*x^2, [x])
            sage: quartic.variables()
            (x, None)
        """
        return self._variables


    def _extract_coefficients(self, monomials):
        """
        Return the coefficients of ``monomials``.

        INPUT:

        - ``polynomial`` -- the input polynomial

        - ``monomials`` -- a list of monomials in the polynomial ring

        - ``variables`` -- a list of variables in the polynomial ring

        OUTPUT:

        A tuple containing the coefficients of the monomials in the given
        polynomial.

        EXAMPLES::

            sage: from sage.rings.invariant_theory import AlgebraicForm
            sage: R.<x,y,z,a30,a21,a12,a03,a20,a11,a02,a10,a01,a00> = QQ[]
            sage: p = ( a30*x^3 + a21*x^2*y + a12*x*y^2 + a03*y^3 + a20*x^2*z +
            ...         a11*x*y*z + a02*y^2*z + a10*x*z^2 + a01*y*z^2 + a00*z^3 )
            sage: base = AlgebraicForm(3, 3, p, [x,y,z])
            sage: m = [x^3, y^3, z^3, x^2*y, x^2*z, x*y^2, y^2*z, x*z^2, y*z^2, x*y*z]
            sage: base._extract_coefficients(m)
            (a30, a03, a00, a21, a20, a12, a02, a10, a01, a11)

            sage: base = AlgebraicForm(3, 3, p.subs(z=1), [x,y])
            sage: m = [x^3, y^3, 1, x^2*y, x^2, x*y^2, y^2, x, y, x*y]
            sage: base._extract_coefficients(m)
            (a30, a03, a00, a21, a20, a12, a02, a10, a01, a11)
        """
        R = self._ring
        if self._homogeneous:
            variables = self._variables
        else:
            variables = self._variables[0:-1]
        indices = [ R.gens().index(x) for x in variables ]
        def index(monomial):
            if monomial in R.base_ring():
                return tuple(0 for i in indices)
            e = monomial.exponents()[0]
            return tuple(e[i] for i in indices)
        coeffs = dict()
        for c,m in self._polynomial:
            i = index(m)
            coeffs[i] = c*m + coeffs.pop(i, R.zero())
        result = tuple(coeffs.pop(index(m), R.zero()) // m for m in monomials)
        if len(coeffs) != 0:
            raise ValueError('The polynomial contains more monomials than given: '+str(coeffs))
        return result


    def coefficients(self):
        """
        Alias for ``coeffs()``

        See the documentation for ``coeffs()`` for details.

        EXAMPLES::

            sage: R.<a,b,c,d,e,f,g, x,y,z> = QQ[]
            sage: p = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*x*z + f*y*z
            sage: q = invariant_theory.quadratic_form(p, x,y,z)
            sage: q.coefficients()
            (a, b, c, d, e, f)
            sage: q.coeffs()
            (a, b, c, d, e, f)
        """
        return self.coeffs()


######################################################################

class QuadraticForm(AlgebraicForm):
    """
    Invariant theory of a multivariate quadratic form

    You should use the :class:`invariant_theory
    <InvariantTheoryFactory>` factory object to contstruct instances
    of this class. See :meth:`~InvariantTheoryFactory.quadratic_form`
    for details.

    TESTS::

        sage: R.<a,b,c,d,e,f,g, x,y,z> = QQ[]
        sage: p = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*x*z + f*y*z
        sage: invariant_theory.quadratic_form(p, x,y,z)
        ternary quadratic with coefficients (a, b, c, d, e, f)
        sage: type(_)
        <class 'sage.rings.invariant_theory.TernaryQuadratic'>

        sage: R.<a,b,c,d,e,f,g, x,y,z> = QQ[]
        sage: p = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*x*z + f*y*z
        sage: invariant_theory.quadratic_form(p, x,y,z)
        ternary quadratic with coefficients (a, b, c, d, e, f)
        sage: type(_)
        <class 'sage.rings.invariant_theory.TernaryQuadratic'>

    Since we can not always decide whether the form is homogeneous or
    not based on the number of variables, you need to explicitly
    specify it if you want the variables to be treated as
    inhomogeneous::

        sage: invariant_theory.inhomogeneous_quadratic_form(p.subs(z=1), x,y)
        ternary quadratic with coefficients (a, b, c, d, e, f)
   """

    def __init__(self, n, polynomial, *args):
        """
        The Python constructor

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: from sage.rings.invariant_theory import QuadraticForm
            sage: QuadraticForm(2, x^2+2*y^2+3*x*y)
            binary quadratic with coefficients (1, 2, 3)
            sage: QuadraticForm(3, x^2+y^2)
            ternary quadratic with coefficients (1, 1, 0, 0, 0, 0)
        """
        super(QuadraticForm, self).__init__(n, 2, polynomial, *args)


    @cached_method
    def monomials(self):
        """
        List the basis monomials in the form.

        OUTPUT:

        A tuple of monomials. They are in the same order as
        :meth:`coeffs`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: quadratic = invariant_theory.quadratic_form(x^2+y^2)
            sage: quadratic.monomials()
            (x^2, y^2, x*y)

            sage: quadratic = invariant_theory.inhomogeneous_quadratic_form(x^2+y^2)
            sage: quadratic.monomials()
            (x^2, y^2, 1, x*y, x, y)
        """
        var = self._variables
        def prod(a,b):
            if a is None and b is None:
                return self._ring.one()
            elif a is None:
                return b
            elif b is None:
                return a
            else:
                return a*b
        squares = tuple( prod(x,x) for x in var )
        mixed = []
        for i in range(self._n):
            for j in range(i+1, self._n):
                mixed.append(prod(var[i], var[j]))
        mixed = tuple(mixed)
        return squares + mixed


    @cached_method
    def coeffs(self):
        r"""
        The coefficients of a quadratic form.

        Given

        .. math::

            f(x) = \sum_{0\leq i<n} a_i x_i^2 + \sum_{0\leq j <k<n}
            a_{jk} x_j x_k

        this function returns `a = (a_0, \dots, a_n, a_{00}, a_{01}, \dots, a_{n-1,n})`

        EXAMPLES::

            sage: R.<a,b,c,d,e,f,g, x,y,z> = QQ[]
            sage: p = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*x*z + f*y*z
            sage: inv = invariant_theory.quadratic_form(p, x,y,z); inv
            ternary quadratic with coefficients (a, b, c, d, e, f)
            sage: inv.coeffs()
            (a, b, c, d, e, f)
            sage: inv.scaled_coeffs()
            (a, b, c, 1/2*d, 1/2*e, 1/2*f)
        """
        return self._extract_coefficients(self.monomials())


    def scaled_coeffs(self):
        """
        The scaled coefficients of a quadratic form.

        Given

        .. math::

            f(x) = \sum_{0\leq i<n} a_i x_i^2 + \sum_{0\leq j <k<n}
            2 a_{jk} x_j x_k

        this function returns `a = (a_0, \cdots, a_n, a_{00}, a_{01}, \dots, a_{n-1,n})`

        EXAMPLES::

            sage: R.<a,b,c,d,e,f,g, x,y,z> = QQ[]
            sage: p = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*x*z + f*y*z
            sage: inv = invariant_theory.quadratic_form(p, x,y,z); inv
            ternary quadratic with coefficients (a, b, c, d, e, f)
            sage: inv.coeffs()
            (a, b, c, d, e, f)
            sage: inv.scaled_coeffs()
            (a, b, c, 1/2*d, 1/2*e, 1/2*f)
        """
        coeff = self.coeffs()
        squares = coeff[0:self._n]
        mixed = tuple( c/2 for c in coeff[self._n:] )
        return squares + mixed


    @cached_method
    def matrix(self):
        """
        Return the quadratic form as a symmetric matrix

        OUTPUT:

        This method returns a symmetric matrix `A` such that the
        quadratic `Q` equals

        .. math::

            Q(x,y,z,\dots) = (x,y,\dots) A (x,y,\dots)^t

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: quadratic = invariant_theory.ternary_quadratic(x^2+y^2+z^2+x*y)
            sage: matrix(quadratic)
            [  1 1/2   0]
            [1/2   1   0]
            [  0   0   1]
            sage: quadratic._matrix_() == matrix(quadratic)
            True
        """
        coeff = self.scaled_coeffs()
        A = matrix(self._ring, self._n)
        for i in range(self._n):
            A[i,i] = coeff[i]
        ij = self._n
        for i in range(self._n):
            for j in range(i+1, self._n):
                A[i,j] = coeff[ij]
                A[j,i] = coeff[ij]
                ij += 1
        return A

    _matrix_ = matrix


    def discriminant(self):
        """
        Return the discriminant of the quadratic form.

        Up to an overall constant factor, this is just the determinant
        of the defining matrix, see :meth:`matrix`. For a quadratic
        form in `n` variables, the overall constant is `2^{n-1}` if
        `n` is odd and `(-1)^{n/2} 2^n` if `n` is even.

        EXAMPLES::

            sage: R.<a,b,c, x,y> = QQ[]
            sage: p = a*x^2+b*x*y+c*y^2
            sage: quadratic = invariant_theory.quadratic_form(p, x,y)
            sage: quadratic.discriminant()
            b^2 - 4*a*c

            sage: R.<a,b,c,d,e,f,g, x,y,z> = QQ[]
            sage: p = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*x*z + f*y*z
            sage: quadratic = invariant_theory.quadratic_form(p, x,y,z)
            sage: quadratic.discriminant()
            4*a*b*c - c*d^2 - b*e^2 + d*e*f - a*f^2
        """
        A = 2*self._matrix_()
        if is_odd(self._n):
            return A.det() / 2
        else:
            return (-1)**(self._n/2) * A.det()



    def as_QuadraticForm(self):
        """
        Convert into a :class:`~sage.quadratic_forms.quadratic_form.QuadraticForm`.

        OUTPUT:

        Sage has a special quadratic forms subsystem. This method
        converts ``self`` into this
        :class:`~sage.quadratic_forms.quadratic_form.QuadraticForm`
        representation.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: p = x^2+y^2+z^2+2*x*y+3*x*z
            sage: quadratic = invariant_theory.ternary_quadratic(p)
            sage: matrix(quadratic)
            [  1   1 3/2]
            [  1   1   0]
            [3/2   0   1]
            sage: quadratic.as_QuadraticForm()
            Quadratic form in 3 variables over Multivariate Polynomial
            Ring in x, y, z over Rational Field with coefficients:
            [ 1/2 1 3/2 ]
            [ * 1/2 0 ]
            [ * * 1/2 ]
            sage: _.polynomial('X,Y,Z')
            X^2 + 2*X*Y + Y^2 + 3*X*Z + Z^2
       """
        R = self._ring
        B = self._matrix_()
        import sage.quadratic_forms.quadratic_form
        return sage.quadratic_forms.quadratic_form.QuadraticForm(R, B)

######################################################################

class BinaryQuartic(AlgebraicForm):
    """
    Invariant theory of a binary quartic

    You should use the :class:`invariant_theory
    <InvariantTheoryFactory>` factory object to contstruct instances
    of this class. See :meth:`~InvariantTheoryFactory.binary_quartic`
    for details.

    TESTS::

        sage: R.<a0, a1, a2, a3, a4, x0, x1> = QQ[]
        sage: p = a0*x1^4 + a1*x1^3*x0 + a2*x1^2*x0^2 + a3*x1*x0^3 + a4*x0^4
        sage: quartic = invariant_theory.binary_quartic(p, x0, x1)
        sage: quartic._check_covariant('form')
        sage: quartic._check_covariant('EisensteinD', invariant=True)
        sage: quartic._check_covariant('EisensteinE', invariant=True)
        sage: quartic._check_covariant('g_covariant')
        sage: quartic._check_covariant('h_covariant')
        sage: TestSuite(quartic).run()
    """

    def __init__(self, polynomial, *args):
        """
        The Python constructor

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: from sage.rings.invariant_theory import BinaryQuartic
            sage: BinaryQuartic(x^4+y^4)
            binary quartic with coefficients (1, 0, 0, 0, 1)
        """
        super(BinaryQuartic, self).__init__(2, 4, polynomial, *args)
        self._x = self._variables[0]
        self._y = self._variables[1]


    @cached_method
    def monomials(self):
        """
        List the basis monomials in the form.

        OUTPUT:

        A tuple of monomials. They are in the same order as
        :meth:`coeffs`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: quartic = invariant_theory.binary_quartic(x^4+y^4)
            sage: quartic.monomials()
            (y^4, x*y^3, x^2*y^2, x^3*y, x^4)
        """
        quartic = self._polynomial
        x0 = self._x
        x1 = self._y
        if self._homogeneous:
            return (x1**4, x1**3*x0, x1**2*x0**2, x1*x0**3, x0**4)
        else:
            return (self._ring.one(), x0, x0**2, x0**3, x0**4)


    @cached_method
    def coeffs(self):
        """
        The coefficients of a binary quartic.

        Given

        .. math::

            f(x) = a_0 x_1^4 + a_1 x_0 x_1^3 + a_2 x_0^2 x_1^2 +
                   a_3 x_0^3 x_1 + a_4 x_0^4

        this function returns `a = (a_0, a_1, a_2, a_3, a_4)`

        EXAMPLES::

            sage: R.<a0, a1, a2, a3, a4, x0, x1> = QQ[]
            sage: p = a0*x1^4 + a1*x1^3*x0 + a2*x1^2*x0^2 + a3*x1*x0^3 + a4*x0^4
            sage: quartic = invariant_theory.binary_quartic(p, x0, x1)
            sage: quartic.coeffs()
            (a0, a1, a2, a3, a4)

            sage: R.<a0, a1, a2, a3, a4, x> = QQ[]
            sage: p = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4
            sage: quartic = invariant_theory.binary_quartic(p, x)
            sage: quartic.coeffs()
            (a0, a1, a2, a3, a4)
        """
        return self._extract_coefficients(self.monomials())


    def scaled_coeffs(self):
        """
        The coefficients of a binary quartic.

        Given

        .. math::

            f(x) = a_0 x_1^4 + 4 a_1 x_0 x_1^3 + 6 a_2 x_0^2 x_1^2 +
                   4 a_3 x_0^3 x_1 + a_4 x_0^4

        this function returns `a = (a_0, a_1, a_2, a_3, a_4)`

        EXAMPLES::

            sage: R.<a0, a1, a2, a3, a4, x0, x1> = QQ[]
            sage: quartic = a0*x1^4 + 4*a1*x1^3*x0 + 6*a2*x1^2*x0^2 + 4*a3*x1*x0^3 + a4*x0^4
            sage: inv = invariant_theory.binary_quartic(quartic, x0, x1)
            sage: inv.scaled_coeffs()
            (a0, a1, a2, a3, a4)

            sage: R.<a0, a1, a2, a3, a4, x> = QQ[]
            sage: quartic = a0 + 4*a1*x + 6*a2*x^2 + 4*a3*x^3 + a4*x^4
            sage: inv = invariant_theory.binary_quartic(quartic, x)
            sage: inv.scaled_coeffs()
            (a0, a1, a2, a3, a4)
        """
        coeff = self.coeffs()
        return (coeff[0], coeff[1]/4, coeff[2]/6, coeff[3]/4, coeff[4])


    @cached_method
    def EisensteinD(self):
        r"""
        One of the Eisenstein invariants of a binary quartic.

        OUTPUT:

        The Eisenstein D-invariant of the quadric.

        .. math::

              f(x) = a_0 x_1^4 + 4 a_1 x_0 x_1^3 + 6 a_2 x_0^2 x_1^2 +
                      4 a_3 x_0^3 x_1 + a_4 x_0^4
              \\
              \Rightarrow
              D(f) = a_0 a_4+3 a_2^2-4 a_1 a_3

        EXAMPLES::

            sage: R.<a0, a1, a2, a3, a4, x0, x1> = QQ[]
            sage: f = a0*x1^4+4*a1*x0*x1^3+6*a2*x0^2*x1^2+4*a3*x0^3*x1+a4*x0^4
            sage: inv = invariant_theory.binary_quartic(f, x0, x1)
            sage: inv.EisensteinD()
            3*a2^2 - 4*a1*a3 + a0*a4
        """
        a = self.scaled_coeffs()
        assert len(a) == 5
        return a[0]*a[4]+3*a[2]**2-4*a[1]*a[3]


    @cached_method
    def EisensteinE(self):
        r"""
        One of the Eisenstein invariants of a binary quartic.

        OUTPUT:

        The Eisenstein E-invariant of the quadric.

        .. math::

            f(x) = a_0 x_1^4 + 4 a_1 x_0 x_1^3 + 6 a_2 x_0^2 x_1^2 +
                   4 a_3 x_0^3 x_1 + a_4 x_0^4
            \\ \Rightarrow
            E(f) = a_0 a_3^2 +a_1^2 a_4 -a_0 a_2 a_4
            -2 a_1 a_2 a_3 + a_2^3

        EXAMPLES::

            sage: R.<a0, a1, a2, a3, a4, x0, x1> = QQ[]
            sage: f = a0*x1^4+4*a1*x0*x1^3+6*a2*x0^2*x1^2+4*a3*x0^3*x1+a4*x0^4
            sage: inv = invariant_theory.binary_quartic(f, x0, x1)
            sage: inv.EisensteinE()
            a2^3 - 2*a1*a2*a3 + a0*a3^2 + a1^2*a4 - a0*a2*a4
        """
        a = self.scaled_coeffs()
        assert len(a) == 5
        return a[0]*a[3]**2 +a[1]**2*a[4] -a[0]*a[2]*a[4] -2*a[1]*a[2]*a[3] +a[2]**3


    @cached_method
    def g_covariant(self):
        r"""
        The g-covariant of a binary quartic.

        OUTPUT:

        The g-covariant of the quadric.

        .. math::

              f(x) = a_0 x_1^4 + 4 a_1 x_0 x_1^3 + 6 a_2 x_0^2 x_1^2 +
                      4 a_3 x_0^3 x_1 + a_4 x_0^4
              \\
              \Rightarrow
              D(f) = \frac{1}{144}
              \begin{pmatrix}
                 \frac{\partial^2 f}{\partial x \partial x}
              \end{pmatrix}

        EXAMPLES::

            sage: R.<a0, a1, a2, a3, a4, x, y> = QQ[]
            sage: p = a0*x^4+4*a1*x^3*y+6*a2*x^2*y^2+4*a3*x*y^3+a4*y^4
            sage: inv = invariant_theory.binary_quartic(p, x, y)
            sage: g = inv.g_covariant();  g
            a1^2*x^4 - a0*a2*x^4 + 2*a1*a2*x^3*y - 2*a0*a3*x^3*y + 3*a2^2*x^2*y^2
            - 2*a1*a3*x^2*y^2 - a0*a4*x^2*y^2 + 2*a2*a3*x*y^3
            - 2*a1*a4*x*y^3 + a3^2*y^4 - a2*a4*y^4

            sage: inv_inhomogeneous = invariant_theory.binary_quartic(p.subs(y=1), x)
            sage: inv_inhomogeneous.g_covariant()
            a1^2*x^4 - a0*a2*x^4 + 2*a1*a2*x^3 - 2*a0*a3*x^3 + 3*a2^2*x^2
            - 2*a1*a3*x^2 - a0*a4*x^2 + 2*a2*a3*x - 2*a1*a4*x + a3^2 - a2*a4

            sage: g == 1/144 * (p.derivative(x,y)^2 - p.derivative(x,x)*p.derivative(y,y))
            True
        """
        a4, a3, a2, a1, a0 = self.scaled_coeffs()
        x0 = self._x
        x1 = self._y
        if self._homogeneous:
            xpow = [x0**4, x0**3 * x1, x0**2 * x1**2, x0 * x1**3, x1**4]
        else:
            xpow = [x0**4, x0**3, x0**2, x0, self._ring.one()]
        return (a1**2 - a0*a2)*xpow[0] + \
            (2*a1*a2 - 2*a0*a3)*xpow[1] + \
            (3*a2**2 - 2*a1*a3 - a0*a4)*xpow[2] + \
            (2*a2*a3 - 2*a1*a4)*xpow[3] + \
            (a3**2 - a2*a4)*xpow[4]


    @cached_method
    def h_covariant(self):
        r"""
        The h-covariant of a binary quartic.

        OUTPUT:

        The h-covariant of the quadric.

        .. math::

              f(x) = a_0 x_1^4 + 4 a_1 x_0 x_1^3 + 6 a_2 x_0^2 x_1^2 +
                      4 a_3 x_0^3 x_1 + a_4 x_0^4
              \\
              \Rightarrow
              D(f) = \frac{1}{144}
              \begin{pmatrix}
                 \frac{\partial^2 f}{\partial x \partial x}
              \end{pmatrix}

        EXAMPLES::

            sage: R.<a0, a1, a2, a3, a4, x, y> = QQ[]
            sage: p = a0*x^4+4*a1*x^3*y+6*a2*x^2*y^2+4*a3*x*y^3+a4*y^4
            sage: inv = invariant_theory.binary_quartic(p, x, y)
            sage: h = inv.h_covariant();   h
            -2*a1^3*x^6 + 3*a0*a1*a2*x^6 - a0^2*a3*x^6 - 6*a1^2*a2*x^5*y + 9*a0*a2^2*x^5*y
            - 2*a0*a1*a3*x^5*y - a0^2*a4*x^5*y - 10*a1^2*a3*x^4*y^2 + 15*a0*a2*a3*x^4*y^2
            - 5*a0*a1*a4*x^4*y^2 + 10*a0*a3^2*x^3*y^3 - 10*a1^2*a4*x^3*y^3
            + 10*a1*a3^2*x^2*y^4 - 15*a1*a2*a4*x^2*y^4 + 5*a0*a3*a4*x^2*y^4
            + 6*a2*a3^2*x*y^5 - 9*a2^2*a4*x*y^5 + 2*a1*a3*a4*x*y^5 + a0*a4^2*x*y^5
            + 2*a3^3*y^6 - 3*a2*a3*a4*y^6 + a1*a4^2*y^6

            sage: inv_inhomogeneous = invariant_theory.binary_quartic(p.subs(y=1), x)
            sage: inv_inhomogeneous.h_covariant()
            -2*a1^3*x^6 + 3*a0*a1*a2*x^6 - a0^2*a3*x^6 - 6*a1^2*a2*x^5 + 9*a0*a2^2*x^5
            - 2*a0*a1*a3*x^5 - a0^2*a4*x^5 - 10*a1^2*a3*x^4 + 15*a0*a2*a3*x^4
            - 5*a0*a1*a4*x^4 + 10*a0*a3^2*x^3 - 10*a1^2*a4*x^3 + 10*a1*a3^2*x^2
            - 15*a1*a2*a4*x^2 + 5*a0*a3*a4*x^2 + 6*a2*a3^2*x - 9*a2^2*a4*x
            + 2*a1*a3*a4*x + a0*a4^2*x + 2*a3^3 - 3*a2*a3*a4 + a1*a4^2

            sage: g = inv.g_covariant()
            sage: h == 1/8 * (p.derivative(x)*g.derivative(y)-p.derivative(y)*g.derivative(x))
            True
        """
        a0, a1, a2, a3, a4 = self.scaled_coeffs()
        x0 = self._x
        x1 = self._y
        if self._homogeneous:
            xpow = [x0**6, x0**5 * x1, x0**4 * x1**2, x0**3 * x1**3,
                    x0**2 * x1**4, x0 * x1**5, x1**6]
        else:
            xpow = [x0**6, x0**5, x0**4, x0**3, x0**2, x0, x0.parent().one()]
        return (-2*a3**3 + 3*a2*a3*a4 - a1*a4**2) * xpow[0] + \
            (-6*a2*a3**2 + 9*a2**2*a4 - 2*a1*a3*a4 - a0*a4**2) * xpow[1] + \
            5 * (-2*a1*a3**2 + 3*a1*a2*a4 - a0*a3*a4) * xpow[2] + \
            10 * (-a0*a3**2 + a1**2*a4) * xpow[3] + \
            5 * (2*a1**2*a3 - 3*a0*a2*a3 + a0*a1*a4) * xpow[4] + \
            (6*a1**2*a2 - 9*a0*a2**2 + 2*a0*a1*a3 + a0**2*a4) * xpow[5] + \
            (2*a1**3 - 3*a0*a1*a2 + a0**2*a3) * xpow[6]


######################################################################
def _covariant_conic(A_scaled_coeffs, B_scaled_coeffs, monomials):
    """
    Helper funciton for :meth:`TernaryQuadratic.covariant_conic`

    INPUT:

    - ``A_scaled_coeffs``, ``B_scaled_coeffs`` -- The scaled
      coefficients of the two ternary quadratics.

    - ``monomials`` -- The monomials :meth:`~TernaryQuadratic.monomials`.

    OUTPUT:

    The so-called covariant conic, a ternary quadratic. It is
    symmetric under exchange of ``A`` and ``B``.

    EXAMPLES::

        sage: ring.<x,y,z> = QQ[]
        sage: A = invariant_theory.ternary_quadratic(x^2+y^2+z^2)
        sage: B = invariant_theory.ternary_quadratic(x*y+x*z+y*z)
        sage: from sage.rings.invariant_theory import _covariant_conic
        sage: _covariant_conic(A.scaled_coeffs(), B.scaled_coeffs(), A.monomials())
        -x*y - x*z - y*z
    """
    a0, b0, c0, h0, g0, f0 = A_scaled_coeffs
    a1, b1, c1, h1, g1, f1 = B_scaled_coeffs
    return (
        (b0*c1+c0*b1-2*f0*f1) * monomials[0] +
        (a0*c1+c0*a1-2*g0*g1) * monomials[1] +
        (a0*b1+b0*a1-2*h0*h1) * monomials[2] +
        2*(f0*g1+g0*f1 -c0*h1-h0*c1) * monomials[3] +
        2*(h0*f1+f0*h1 -b0*g1-g0*b1) * monomials[4] +
        2*(g0*h1+h0*g1 -a0*f1-f0*a1) * monomials[5]  )


######################################################################
class TernaryQuadratic(QuadraticForm):
    """
    Invariant theory of a ternary quadratic

    You should use the :class:`invariant_theory
    <InvariantTheoryFactory>` factory object to contstruct instances
    of this class. See
    :meth:`~InvariantTheoryFactory.ternary_quadratic` for details.

    TESTS::

        sage: R.<x,y,z> = QQ[]
        sage: quadratic = invariant_theory.ternary_quadratic(x^2+y^2+z^2)
        sage: quadratic
        ternary quadratic with coefficients (1, 1, 1, 0, 0, 0)
        sage: TestSuite(quadratic).run()
    """

    def __init__(self, polynomial, *args):
        """
        The Python constructor.

        INPUT:

        See :meth:`~InvariantTheoryFactory.ternary_quadratic`.

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: from sage.rings.invariant_theory import TernaryQuadratic
            sage: TernaryQuadratic(x^2+y^2+z^2)
            ternary quadratic with coefficients (1, 1, 1, 0, 0, 0)
        """
        super(QuadraticForm, self).__init__(3, 2, polynomial, *args)
        self._x = self._variables[0]
        self._y = self._variables[1]
        self._z = self._variables[2]


    @cached_method
    def monomials(self):
        """
        List the basis monomials in the form.

        OUTPUT:

        A tuple of monomials. They are in the same order as
        :meth:`coeffs`.

        EXAMPLES::


            sage: R.<x,y,z> = QQ[]
            sage: quadratic = invariant_theory.ternary_quadratic(x^2+y*z)
            sage: quadratic.monomials()
            (x^2, y^2, z^2, x*y, x*z, y*z)
        """
        R = self._ring
        x,y,z = self._x, self._y, self._z
        if self._homogeneous:
            return (x**2, y**2, z**2, x*y, x*z, y*z)
        else:
            return (x**2, y**2, R.one(), x*y, x, y)


    @cached_method
    def coeffs(self):
        """
        Return the coefficients of a quadratic

        Given

        .. math::

            p(x,y) =&\;
            a_{20} x^{2} + a_{11} x y + a_{02} y^{2} +
            a_{10} x + a_{01} y + a_{00}

        this function returns
        `a = (a_{20}, a_{02}, a_{00}, a_{11}, a_{10}, a_{01} )`

        EXAMPLES::

            sage: R.<x,y,z,a20,a11,a02,a10,a01,a00> = QQ[]
            sage: p = ( a20*x^2 + a11*x*y + a02*y^2 +
            ...         a10*x*z + a01*y*z + a00*z^2 )
            sage: invariant_theory.ternary_quadratic(p, x,y,z).coeffs()
            (a20, a02, a00, a11, a10, a01)
            sage: invariant_theory.ternary_quadratic(p.subs(z=1), x, y).coeffs()
            (a20, a02, a00, a11, a10, a01)
        """
        return self._extract_coefficients(self.monomials())


    def scaled_coeffs(self):
        """
        Return the scaled coefficients of a quadratic

        Given

        .. math::

            p(x,y) =&\;
            a_{20} x^{2} + a_{11} x y + a_{02} y^{2} +
            a_{10} x + a_{01} y + a_{00}

        this function returns
        `a = (a_{20}, a_{02}, a_{00}, a_{11}/2, a_{10}/2, a_{01}/2, )`

        EXAMPLES::

            sage: R.<x,y,z,a20,a11,a02,a10,a01,a00> = QQ[]
            sage: p = ( a20*x^2 + a11*x*y + a02*y^2 +
            ...         a10*x*z + a01*y*z + a00*z^2 )
            sage: invariant_theory.ternary_quadratic(p, x,y,z).scaled_coeffs()
            (a20, a02, a00, 1/2*a11, 1/2*a10, 1/2*a01)
            sage: invariant_theory.ternary_quadratic(p.subs(z=1), x, y).scaled_coeffs()
            (a20, a02, a00, 1/2*a11, 1/2*a10, 1/2*a01)
        """
        F = self._ring.base_ring()
        a200, a020, a002, a110, a101, a011 = self.coeffs()
        return (a200, a020, a002, a110/F(2), a101/F(2), a011/F(2))


    @cached_method
    def dual(self):
        """
        Return the dual ternary quadratic

        EXAMPLES::


            sage: R.<a,b,c,x,y,z> = QQ[]
            sage: cubic = x^2+y^2+z^2
            sage: quadratic = invariant_theory.ternary_quadratic(a*x^2+b*y^2+c*z^2, [x,y,z])
            sage: quadratic.form()
            a*x^2 + b*y^2 + c*z^2
            sage: quadratic.dual().form()
            b*c*x^2 + a*c*y^2 + a*b*z^2

            sage: R.<x,y,z, t> = QQ[]
            sage: cubic = x^2+y^2+z^2
            sage: quadratic = invariant_theory.ternary_quadratic(x^2+y^2+z^2 + t*x*y, [x,y,z])
            sage: quadratic.dual()
            ternary quadratic with coefficients (1, 1, -1/4*t^2 + 1, -t, 0, 0)

            sage: R.<x,y, t> = QQ[]
            sage: quadratic = invariant_theory.ternary_quadratic(x^2+y^2+1 + t*x*y, [x,y])
            sage: quadratic.dual()
            ternary quadratic with coefficients (1, 1, -1/4*t^2 + 1, -t, 0, 0)

        TESTS::

            sage: R = PolynomialRing(QQ, 'a20,a11,a02,a10,a01,a00,x,y,z', order='lex')
            sage: R.inject_variables()
            Defining a20, a11, a02, a10, a01, a00, x, y, z
            sage: p = ( a20*x^2 + a11*x*y + a02*y^2 +
            ...         a10*x*z + a01*y*z + a00*z^2 )
            sage: quadratic = invariant_theory.ternary_quadratic(p, x,y,z)
            sage: transformation = {x:x + 3*y, y:4*x + 13*y - z, z:-5*y + 6*z}
            sage: inverse = {x:73*x - 18*y - 3*z, y:-24*x + 6*y + z, z:-20*x + 5*y + z}

            quadratic.dual().form().subs(transformation)
            quadratic.dual().form().subs(inverse)
        """
        A = self.matrix()
        Aadj = A.adjoint()
        if self._homogeneous:
            var = [self._x, self._y, self._z]
        else:
            var = [self._x, self._y, 1]
        R = self._ring
        p = sum([ sum([ Aadj[i,j]*var[i]*var[j] for i in range(3) ]) for j in range(3)])
        return TernaryQuadratic(p, self.variables())


    def covariant_conic(self, other):
        """
        Returns the ternary quadratic covariant to ``self`` and ``other``

        INPUT:

        - ``other`` -- Another ternary quadratic.

        OUTPUT:

        The so-called covariant conic, a ternary quadratic. It is
        symmetric under exchange of ``self`` and ``other``.

        EXAMPLES::

            sage: ring.<x,y,z> = QQ[]
            sage: Q = invariant_theory.ternary_quadratic(x^2+y^2+z^2)
            sage: R = invariant_theory.ternary_quadratic(x*y+x*z+y*z)
            sage: Q.covariant_conic(R)
            -x*y - x*z - y*z
            sage: R.covariant_conic(Q)
            -x*y - x*z - y*z

        TESTS::

            sage: R.<a,a_,b,b_,c,c_,f,f_,g,g_,h,h_,x,y,z> = QQ[]
            sage: p = ( a*x^2 + 2*h*x*y + b*y^2 +
            ...         2*g*x*z + 2*f*y*z + c*z^2 )
            sage: Q = invariant_theory.ternary_quadratic(p, [x,y,z])
            sage: Q.matrix()
            [a h g]
            [h b f]
            [g f c]
            sage: p = ( a_*x^2 + 2*h_*x*y + b_*y^2 +
            ...         2*g_*x*z + 2*f_*y*z + c_*z^2 )
            sage: Q_ = invariant_theory.ternary_quadratic(p, [x,y,z])
            sage: Q_.matrix()
            [a_ h_ g_]
            [h_ b_ f_]
            [g_ f_ c_]
            sage: QQ_ = Q.covariant_conic(Q_)
            sage: invariant_theory.ternary_quadratic(QQ_, [x,y,z]).matrix()
            [      b_*c + b*c_ - 2*f*f_  f_*g + f*g_ - c_*h - c*h_ -b_*g - b*g_ + f_*h + f*h_]
            [ f_*g + f*g_ - c_*h - c*h_       a_*c + a*c_ - 2*g*g_ -a_*f - a*f_ + g_*h + g*h_]
            [-b_*g - b*g_ + f_*h + f*h_ -a_*f - a*f_ + g_*h + g*h_       a_*b + a*b_ - 2*h*h_]
        """
        return _covariant_conic(self.scaled_coeffs(), other.scaled_coeffs(),
                                self.monomials())


######################################################################

class TernaryCubic(AlgebraicForm):
    """
    Invariant theory of a ternary cubic

    You should use the :class:`invariant_theory
    <InvariantTheoryFactory>` factory object to contstruct instances
    of this class. See :meth:`~InvariantTheoryFactory.ternary_cubic`
    for details.

    TESTS::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+z^3)
        sage: cubic
        ternary cubic with coefficients (1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
        sage: TestSuite(cubic).run()
    """

    def __init__(self, polynomial, *args):
        """
        The Python constructor

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: p = 2837*x^3 + 1363*x^2*y + 6709*x^2*z + \
            ...     5147*x*y^2 + 2769*x*y*z + 912*x*z^2 + 4976*y^3 + \
            ...     2017*y^2*z + 4589*y*z^2 + 9681*z^3
            sage: cubic = invariant_theory.ternary_cubic(p)
            sage: cubic._check_covariant('S_invariant', invariant=True)
            sage: cubic._check_covariant('T_invariant', invariant=True)
            sage: cubic._check_covariant('form')
            sage: cubic._check_covariant('Hessian')
            sage: cubic._check_covariant('Theta_covariant')
            sage: cubic._check_covariant('J_covariant')
        """
        super(TernaryCubic, self).__init__(3, 3, polynomial, *args)
        self._x = self._variables[0]
        self._y = self._variables[1]
        self._z = self._variables[2]


    @cached_method
    def monomials(self):
        """
        List the basis monomials in the form.

        OUTPUT:

        A tuple of monomials. They are in the same order as
        :meth:`coeffs`.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y*z^2)
            sage: cubic.monomials()
            (x^3, y^3, z^3, x^2*y, x^2*z, x*y^2, y^2*z, x*z^2, y*z^2, x*y*z)
        """
        R = self._ring
        x,y,z = self._x, self._y, self._z
        if self._homogeneous:
            return (x**3, y**3, z**3, x**2*y, x**2*z, x*y**2,
                    y**2*z, x*z**2, y*z**2, x*y*z)
        else:
            return (x**3, y**3, R.one(), x**2*y, x**2, x*y**2,
                    y**2, x, y, x*y)


    @cached_method
    def coeffs(self):
        r"""
        Return the coefficients of a cubic

        Given

        .. math::

            \begin{split}
              p(x,y) =&\;
              a_{30} x^{3} + a_{21} x^{2} y + a_{12} x y^{2} +
              a_{03} y^{3} + a_{20} x^{2} +
              \\ &\;
              a_{11} x y +
              a_{02} y^{2} + a_{10} x + a_{01} y + a_{00}
            \end{split}

        this function returns
        `a = (a_{30}, a_{03}, a_{00}, a_{21}, a_{20}, a_{12}, a_{02}, a_{10}, a_{01}, a_{11})`

        EXAMPLES::

            sage: R.<x,y,z,a30,a21,a12,a03,a20,a11,a02,a10,a01,a00> = QQ[]
            sage: p = ( a30*x^3 + a21*x^2*y + a12*x*y^2 + a03*y^3 + a20*x^2*z +
            ...         a11*x*y*z + a02*y^2*z + a10*x*z^2 + a01*y*z^2 + a00*z^3 )
            sage: invariant_theory.ternary_cubic(p, x,y,z).coeffs()
            (a30, a03, a00, a21, a20, a12, a02, a10, a01, a11)
            sage: invariant_theory.ternary_cubic(p.subs(z=1), x, y).coeffs()
            (a30, a03, a00, a21, a20, a12, a02, a10, a01, a11)
        """
        return self._extract_coefficients(self.monomials())


    def scaled_coeffs(self):
        r"""
        Return the coefficients of a cubic

        Compared to :meth:`coeffs`, this method returns rescaled
        coefficiens that are often used in invariant theory.

        Given

        .. math::

            \begin{split}
              p(x,y) =&\;
              a_{30} x^{3} + a_{21} x^{2} y + a_{12} x y^{2} +
              a_{03} y^{3} + a_{20} x^{2} +
              \\ &\;
              a_{11} x y +
              a_{02} y^{2} + a_{10} x + a_{01} y + a_{00}
            \end{split}

        this function returns
        `a = (a_{30}, a_{03}, a_{00}, a_{21}/3, a_{20}/3, a_{12}/3, a_{02}/3, a_{10}/3, a_{01}/3, a_{11}/6)`

        EXAMPLES::

            sage: R.<x,y,z,a30,a21,a12,a03,a20,a11,a02,a10,a01,a00> = QQ[]
            sage: p = ( a30*x^3 + a21*x^2*y + a12*x*y^2 + a03*y^3 + a20*x^2*z +
            ...         a11*x*y*z + a02*y^2*z + a10*x*z^2 + a01*y*z^2 + a00*z^3 )
            sage: invariant_theory.ternary_cubic(p, x,y,z).scaled_coeffs()
            (a30, a03, a00, 1/3*a21, 1/3*a20, 1/3*a12, 1/3*a02, 1/3*a10, 1/3*a01, 1/6*a11)
        """
        a = self.coeffs()
        F = self._ring.base_ring()
        return (a[0], a[1], a[2],
                1/F(3)*a[3], 1/F(3)*a[4], 1/F(3)*a[5],
                1/F(3)*a[6], 1/F(3)*a[7], 1/F(3)*a[8],
                1/F(6)*a[9])


    def S_invariant(self):
        """
        Return the S-invariant

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^2*y+y^3+z^3+x*y*z)
            sage: cubic.S_invariant()
            -1/1296
        """
        a,b,c,a2,a3,b1,b3,c1,c2,m = self.scaled_coeffs()
        S = ( a*b*c*m-(b*c*a2*a3+c*a*b1*b3+a*b*c1*c2)
              -m*(a*b3*c2+b*c1*a3+c*a2*b1)
              +(a*b1*c2**2+a*c1*b3**2+b*a2*c1**2+b*c2*a3**2+c*b3*a2**2+c*a3*b1**2)
              -m**4+2*m**2*(b1*c1+c2*a2+a3*b3)
              -3*m*(a2*b3*c1+a3*b1*c2)
              -(b1**2*c1**2+c2**2*a2**2+a3**2*b3**2)
              +(c2*a2*a3*b3+a3*b3*b1*c1+b1*c1*c2*a2) )
        return S


    def T_invariant(self):
        """
        Return the T-invariant

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+z^3)
            sage: cubic.T_invariant()
            1

            sage: R.<x,y,z,t> = GF(7)[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+z^3+t*x*y*z, [x,y,z])
            sage: cubic.T_invariant()
            -t^6 - t^3 + 1
        """
        a,b,c,a2,a3,b1,b3,c1,c2,m = self.scaled_coeffs()
        T = ( a**2*b**2*c**2-6*a*b*c*(a*b3*c2+b*c1*a3+c*a2*b1)
              -20*a*b*c*m**3+12*a*b*c*m*(b1*c1+c2*a2+a3*b3)
              +6*a*b*c*(a2*b3*c1+a3*b1*c2)+
              4*(a**2*b*c2**3+a**2*c*b3**3+b**2*c*a3**3+
                 b**2*a*c1**3+c**2*a*b1**3+c**2*b*a2**3)
              +36*m**2*(b*c*a2*a3+c*a*b1*b3+a*b*c1*c2)
              -24*m*(b*c*b1*a3**2+b*c*c1*a2**2+c*a*c2*b1**2+c*a*a2*b3**2+a*b*a3*c2**2+
                     a*b*b3*c1**2)
              -3*(a**2*b3**2*c2**2+b**2*c1**2*a3**2+c**2*a2**2*b1**2)+
              18*(b*c*b1*c1*a2*a3+c*a*c2*a2*b3*b1+a*b*a3*b3*c1*c2)
              -12*(b*c*c2*a3*a2**2+b*c*b3*a2*a3**2+c*a*c1*b3*b1**2+
                   c*a*a3*b1*b3**2+a*b*a2*c1*c2**2+a*b*b1*c2*c1**2)
              -12*m**3*(a*b3*c2+b*c1*a3+c*a2*b1)
              +12*m**2*(a*b1*c2**2+a*c1*b3**2+b*a2*c1**2+
                        b*c2*a3**2+c*b3*a2**2+c*a3*b1**2)
              -60*m*(a*b1*b3*c1*c2+b*c1*c2*a2*a3+c*a2*a3*b1*b3)
              +12*m*(a*a2*b3*c2**2+a*a3*c2*b3**2+b*b3*c1*a3**2+
                     b*b1*a3*c1**2+c*c1*a2*b1**2+c*c2*b1*a2**2)
              +6*(a*b3*c2+b*c1*a3+c*a2*b1)*(a2*b3*c1+a3*b1*c2)
              +24*(a*b1*b3**2*c1**2+a*c1*c2**2*b1**2+b*c2*c1**2*a2**2
                   +b*a2*a3**2*c2**2+c*a3*a2**2*b3**2+c*b3*b1**2*a3**2)
              -12*(a*a2*b1*c2**3+a*a3*c1*b3**3+b*b3*c2*a3**3+b*b1*a2*c1**3
                   +c*c1*a3*b1**3+c*c2*b3*a2**3)
              -8*m**6+24*m**4*(b1*c1+c2*a2+a3*b3)-36*m**3*(a2*b3*c1+a3*b1*c2)
              -12*m**2*(b1*c1*c2*a2+c2*a2*a3*b3+a3*b3*b1*c1)
              -24*m**2*(b1**2*c1**2+c2**2*a2**2+a3**2*b3**2)
              +36*m*(a2*b3*c1+a3*b1*c2)*(b1*c1+c2*a2+a3*b3)
              +8*(b1**3*c1**3+c2**3*a2**3+a3**3*b3**3)
              -27*(a2**2*b3**2*c1**2+a3**2*b1**2*c2**2)-6*b1*c1*c2*a2*a3*b3
              -12*(b1**2*c1**2*c2*a2+b1**2*c1**2*a3*b3+c2**2*a2**2*a3*b3+
                   c2**2*a2**2*b1*c1+a3**2*b3**2*b1*c1+a3**2*b3**2*c2*a2) )
        return T


    @cached_method
    def polar_conic(self):
        """
        Return the polar conic of the cubic

        OUTPUT:

        Given the ternary cubic `f(X,Y,Z)`, this method returns the
        symmetric matrix `A(x,y,z)` defined by

        .. math::

            x f_X + y f_Y + z f_Z = (X,Y,Z) \cdot A(x,y,z) \cdot (X,Y,Z)^t

        EXAMPLES::

            sage: R.<x,y,z,X,Y,Z,a30,a21,a12,a03,a20,a11,a02,a10,a01,a00> = QQ[]
            sage: p = ( a30*x^3 + a21*x^2*y + a12*x*y^2 + a03*y^3 + a20*x^2*z +
            ...         a11*x*y*z + a02*y^2*z + a10*x*z^2 + a01*y*z^2 + a00*z^3 )
            sage: cubic = invariant_theory.ternary_cubic(p, x,y,z)
            sage: cubic.polar_conic()
            [  3*x*a30 + y*a21 + z*a20 x*a21 + y*a12 + 1/2*z*a11 x*a20 + 1/2*y*a11 + z*a10]
            [x*a21 + y*a12 + 1/2*z*a11   x*a12 + 3*y*a03 + z*a02 1/2*x*a11 + y*a02 + z*a01]
            [x*a20 + 1/2*y*a11 + z*a10 1/2*x*a11 + y*a02 + z*a01   x*a10 + y*a01 + 3*z*a00]

            sage: polar_eqn = X*p.derivative(x) + Y*p.derivative(y) + Z*p.derivative(z)
            sage: polar = invariant_theory.ternary_quadratic(polar_eqn, [x,y,z])
            sage: polar.matrix().subs(X=x,Y=y,Z=z) == cubic.polar_conic()
            True
        """
        a30, a03, a00, a21, a20, a12, a02, a10, a01, a11 = self.coeffs()
        if self._homogeneous:
            x,y,z = self.variables()
        else:
            x,y,z = (self._x, self._y, 1)
        F = self._ring.base_ring()
        A00 = 3*x*a30 + y*a21 + z*a20
        A11 = x*a12 + 3*y*a03 + z*a02
        A22 = x*a10 + y*a01 + 3*z*a00
        A01 = x*a21 + y*a12 + 1/F(2)*z*a11
        A02 = x*a20 + 1/F(2)*y*a11 + z*a10
        A12 = 1/F(2)*x*a11 + y*a02 + z*a01
        polar = matrix(self._ring, [[A00, A01, A02],[A01, A11, A12],[A02, A12, A22]])
        return polar


    @cached_method
    def Hessian(self):
        """
        Return the Hessian covariant

        OUTPUT:

        The Hessian matrix multiplied with the conventional
        normalization factor `1/216`.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+z^3)
            sage: cubic.Hessian()
            x*y*z

            sage: R.<x,y> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+1)
            sage: cubic.Hessian()
            x*y
        """
        a30, a03, a00, a21, a20, a12, a02, a10, a01, a11 = self.coeffs()
        if self._homogeneous:
            x, y, z = self.variables()
        else:
            x, y, z = self._x, self._y, 1
        Uxx = 6*x*a30 + 2*y*a21 + 2*z*a20
        Uxy = 2*x*a21 + 2*y*a12 + z*a11
        Uxz = 2*x*a20 + y*a11 + 2*z*a10
        Uyy = 2*x*a12 + 6*y*a03 + 2*z*a02
        Uyz = x*a11 + 2*y*a02 + 2*z*a01
        Uzz = 2*x*a10 + 2*y*a01 + 6*z*a00
        H = matrix(self._ring, [[Uxx, Uxy, Uxz],[Uxy, Uyy, Uyz],[Uxz, Uyz, Uzz]])
        F = self._ring.base_ring()
        return 1/F(216) * H.det()


    def Theta_covariant(self):
        """
        Return the `\Theta` covariant

        EXAMPLES::


            sage: R.<x,y,z> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+z^3)
            sage: cubic.Theta_covariant()
            -x^3*y^3 - x^3*z^3 - y^3*z^3

            sage: R.<x,y> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+1)
            sage: cubic.Theta_covariant()
            -x^3*y^3 - x^3 - y^3

            sage: R.<x,y,z,a30,a21,a12,a03,a20,a11,a02,a10,a01,a00> = QQ[]
            sage: p = ( a30*x^3 + a21*x^2*y + a12*x*y^2 + a03*y^3 + a20*x^2*z +
            ...         a11*x*y*z + a02*y^2*z + a10*x*z^2 + a01*y*z^2 + a00*z^3 )
            sage: cubic = invariant_theory.ternary_cubic(p, x,y,z)
            sage: len(list(cubic.Theta_covariant()))
            6952
        """
        U_conic = self.polar_conic().adjoint()
        U_coeffs = ( U_conic[0,0], U_conic[1,1], U_conic[2,2],
                     U_conic[0,1], U_conic[0,2], U_conic[1,2] )
        H_conic = TernaryCubic(self.Hessian(), self.variables()).polar_conic().adjoint()
        H_coeffs = ( H_conic[0,0], H_conic[1,1], H_conic[2,2],
                     H_conic[0,1], H_conic[0,2], H_conic[1,2] )
        quadratic = TernaryQuadratic(self._ring.zero(), self.variables())
        F = self._ring.base_ring()
        return 1/F(9) * _covariant_conic(U_coeffs, H_coeffs, quadratic.monomials())


    def J_covariant(self):
        """
        Return the J-covariant of the ternary cubic

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+z^3)
            sage: cubic.J_covariant()
            x^6*y^3 - x^3*y^6 - x^6*z^3 + y^6*z^3 + x^3*z^6 - y^3*z^6

            sage: R.<x,y> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+1)
            sage: cubic.J_covariant()
            x^6*y^3 - x^3*y^6 - x^6 + y^6 + x^3 - y^3
        """
        def diff(p,d):
            px = p.derivative(self._x)
            py = p.derivative(self._y)
            if self._homogeneous:
                pz = p.derivative(self._z)
            else:
                pz = d*p -self._x*px -self._y*py
            return (px, py, pz)

        grad_U = diff(self.form() ,3)
        grad_H = diff(self.Hessian(), 3)
        grad_T = diff(self.Theta_covariant(), 6)
        F = self._ring.base_ring()
        return 1/F(9) * matrix(self._ring, [grad_U, grad_H, grad_T]).det()




######################################################################

class InvariantTheoryFactory(SageObject):
    """
    Factory object for invariants of multilinear forms

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: invariant_theory.ternary_cubic(x^3+y^3+z^3)
        ternary cubic with coefficients (1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
    """

    def quadratic_form(self, polynomial, *args):
        """
        Invariants of a homogeneous quadratic form

        INPUT:

        - ``polynomial`` -- a homogeneous or inhomogenous quadratic form.

        - ``*args`` -- the variables as multiple arguments, or as a
          single list/tuple. If the last argument is ``None``, the
          cubic is assumed to be inhomogeneous.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: quadratic = x^2+y^2+z^2
            sage: inv = invariant_theory.quadratic_form(quadratic)
            sage: type(inv)
            <class 'sage.rings.invariant_theory.TernaryQuadratic'>

        If some of the ring variables are to be treated as coefficients
        you need to specify the polynomial variables::

            sage: R.<x,y,z, a,b> = QQ[]
            sage: quadratic = a*x^2+b*y^2+z^2+2*y*z
            sage: invariant_theory.quadratic_form(quadratic, x,y,z)
            ternary quadratic with coefficients (a, b, 1, 0, 0, 2)
            sage: invariant_theory.quadratic_form(quadratic, [x,y,z])   # alternate syntax
            ternary quadratic with coefficients (a, b, 1, 0, 0, 2)

        Inhomogeneous quadratic forms (see also
        :meth:`inhomogeneous_quadratic_form`) can be specified by
        passing ``None`` as the last variable::

            sage: inhom = quadratic.subs(z=1)
            sage: invariant_theory.quadratic_form(inhom, x,y,None)
            ternary quadratic with coefficients (a, b, 1, 0, 0, 2)
        """
        variables = _guess_variables(polynomial, *args)
        n = len(variables)
        if n == 3:
            return TernaryQuadratic(polynomial, *args)
        else:
            return QuadraticForm(n, polynomial, *args)

    def inhomogeneous_quadratic_form(self, polynomial, *args):
        """
        Invariants of an inhomogeneous quadratic form

        INPUT:

        - ``polynomial`` -- an inhomogenous quadratic form.

        - ``*args`` -- the variables as multiple arguments, or as a
          single list/tuple.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: quadratic = x^2+2*y^2+3*x*y+4*x+5*y+6
            sage: inv = invariant_theory.inhomogeneous_quadratic_form(quadratic)
            sage: type(inv)
            <class 'sage.rings.invariant_theory.TernaryQuadratic'>
        """
        variables = _guess_variables(polynomial, *args)
        n = len(variables) + 1
        if n == 3:
            return TernaryQuadratic(polynomial, *args)
        else:
            return QuadraticForm(n, polynomial, *args, **kwds)

    def binary_quartic(self, quartic, *args, **kwds):
        """
        Invariant theory of a quartic in two variables

        The algebra of invariants of a quartic form is generated by
        invariants `i`, `j` of degrees 2, 3.  This ring is naturally
        isomorphic to the ring of modular forms of level 1, with the
        two generators corresponding to the Eisenstein series `E_4`
        (see
        :meth:`~sage.rings.invariant_theory.BinaryQuartic.EisensteinD`)
        and `E_6` (see
        :meth:`~sage.rings.invariant_theory.BinaryQuartic.EisensteinE`). The
        algebra of covariants is generated by these two invariants
        together with the form `f` of degree 1 and order 4, the
        Hessian `g` (see :meth:`~BinaryQuartic.g_covariant`) of degree
        2 and order 4, and a covariant `h` (see
        :meth:`~BinaryQuartic.h_covariant`) of degree 3 and order
        6. They are related by a syzygy

        .. math::

            j f^3 - g f^2 i + 4 g^3 + h^2 = 0

        of degree 6 and order 12.

        INPUT:

        - ``quartic`` -- a quartic.

        - ``x``, ``y`` -- the homogeneous variables. If ``y`` is
          ``None``, the quartic is assumed to be inhomogeneous.

        REFERENCES:

        ..  [WpBinaryForm]
            http://en.wikipedia.org/wiki/Invariant_of_a_binary_form

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: quartic = invariant_theory.binary_quartic(x^4+y^4)
            sage: quartic
            binary quartic with coefficients (1, 0, 0, 0, 1)
            sage: type(quartic)
            <class 'sage.rings.invariant_theory.BinaryQuartic'>
        """
        return BinaryQuartic(quartic, *args, **kwds)

    def ternary_quadratic(self, quadratic, *args, **kwds):
        """
        Invariants of a quadratic in three variables

        INPUT:

        - ``quadratic`` -- a homogeneous quadratic in 3 homogeneous
          variables, or an inhomogeneous quadratic in 2 variables.

        - ``x``, ``y``, ``z`` -- the variables. If ``z`` is ``None``,
          the quadratic is assumed to be inhomogeneous.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: cubic = x^2+y^2+z^2
            sage: inv = invariant_theory.ternary_quadratic(cubic)
            sage: type(inv)
            <class 'sage.rings.invariant_theory.TernaryQuadratic'>
        """
        return TernaryQuadratic(quadratic, *args, **kwds)

    def ternary_cubic(self, cubic, *args, **kwds):
        r"""
        Invariants of a cubic in three variables

        The algebra of invariants of a ternary cubic under `SL_3(\CC)`
        is a polynomial algebra generated by two invariants `S` (see
        :meth:`~sage.rings.invariant_theory.TernaryCubic.S_invariant`)
        and T (see
        :meth:`~sage.rings.invariant_theory.TernaryCubic.T_invariant`)
        of degrees 4 and 6, called Aronhold invariants.

        The ring of covariants is given as follows. The identity
        covariant U of a ternary cubic has degree 1 and order 3.  The
        Hessian `H` (see
        :meth:`~sage.rings.invariant_theory.TernaryCubic.Hessian`)
        is a covariant of ternary cubics of degree 3 and order 3.
        There is a covariant `\Theta` (see
        :meth:`~sage.rings.invariant_theory.TernaryCubic.Theta_covariant`)
        of ternary cubics of degree 8 and order 6 that vanishes on
        points `x` lying on the Salmon conic of the polar of `x` with
        respect to the curve and its Hessian curve. The Brioschi
        covariant `J` (see
        :meth:`~sage.rings.invariant_theory.TernaryCubic.J_covariant`)
        is the Jacobian of `U`, `\Theta`, and `H` of degree 12, order
        9.  The algebra of covariants of a ternary cubic is generated
        over the ring of invariants by `U`, `\Theta`, `H`, and `J`,
        with a relation

        .. math::

            \begin{split}
              J^2 =& 4 \Theta^3 + T U^2 \Theta^2 +
              \Theta (-4 S^3 U^4 + 2 S T U^3 H
              - 72 S^2 U^2 H^2
              \\ &
              - 18 T U H^3 +  108 S H^4)
              -16 S^4 U^5 H - 11 S^2 T U^4 H^2
              \\ &
              -4 T^2 U^3 H^3
              +54 S T U^2 H^4 -432 S^2 U H^5 -27 T H^6
            \end{split}


        REFERENCES:

        ..  [WpTernaryCubic]
            http://en.wikipedia.org/wiki/Ternary_cubic

        INPUT:

        - ``cubic`` -- a homogeneous cubic in 3 homogeneous variables,
          or an inhomogeneous cubic in 2 variables.

        - ``x``, ``y``, ``z`` -- the variables. If ``z`` is ``None``, the
          cubic is assumed to be inhomogeneous.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: cubic = invariant_theory.ternary_cubic(x^3+y^3+z^3)
            sage: type(cubic)
            <class 'sage.rings.invariant_theory.TernaryCubic'>

        There is one syzygy of the covariants in degree 18, we use this to
        check that the expressions for the covariants are correct::

            sage: R.<x,y,z> = QQ[]
            sage: monomials = (x^3, y^3, z^3, x^2*y, x^2*z, x*y^2,
            ...                y^2*z, x*z^2, y*z^2, x*y*z)
            sage: random_poly = sum([ randint(0,10000) * m for m in monomials ])
            sage: cubic = invariant_theory.ternary_cubic(random_poly)
            sage: U = cubic.form()
            sage: S = cubic.S_invariant()
            sage: T = cubic.T_invariant()
            sage: H = cubic.Hessian()
            sage: Theta = cubic.Theta_covariant()
            sage: J = cubic.J_covariant()

            sage: ( -J^2 + 4*Theta^3 + T*U^2*Theta^2 +
            ...     Theta*(-4*S^3*U^4 + 2*S*T*U^3*H - 72*S^2*U^2*H^2
            ...            - 18*T*U*H^3 +  108*S*H^4)
            ...     -16*S^4*U^5*H - 11*S^2*T*U^4*H^2 -4*T^2*U^3*H^3
            ...     +54*S*T*U^2*H^4 -432*S^2*U*H^5 -27*T*H^6 )
            0
        """
        return TernaryCubic(cubic, *args, **kwds)


invariant_theory = InvariantTheoryFactory()

