"""
Quaternion Algebra Elements.


"""

#*****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#       Copyright (C) 2009 Jonathon Bober <jwbober@gmail.com>
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

from sage.structure.element cimport AlgebraElement, RingElement, ModuleElement, Element
from sage.algebras.quaternion_algebra_element cimport QuaternionAlgebraElement_abstract
from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.arith import lcm
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
from sage.rings.number_field.number_field_element cimport NumberFieldElement

include "../ext/gmp.pxi"
include "../ext/stdsage.pxi"

include "../libs/flint/fmpz.pxi"
include "../libs/flint/fmpz_poly.pxi"
include "../libs/flint/ntl_interface.pxd"

# variables for holding temporary values computed in
# QuaternionAlgebraElement_rational_field._mul_()
cdef mpz_t T1, T2, t3, t4, t5, t6, t7, t8, s1, s2, U1, U2
cdef fmpz_poly_t fT1, fT2, ft3, ft4, ft5, ft6, ft7, ft8, fs1, fs2, fU1, fU2

def _init_globals():
    """
    Clear all global variables allocated for optimization of
    quaternion algebra arithmetic.

    Do *not* call this unless you want to leak memory.

    EXAMPLES::

        sage: sage.algebras.quaternion_algebra_element._clear_globals()
        sage: sage.algebras.quaternion_algebra_element._init_globals()
    """
    # over QQ
    mpz_init(T1)
    mpz_init(T2)
    mpz_init(t3)
    mpz_init(t4)
    mpz_init(t5)
    mpz_init(t6)
    mpz_init(t7)
    mpz_init(t8)

    mpz_init(s1)
    mpz_init(s2)

    mpz_init(U1)
    mpz_init(U2)

    # Number fields
    fmpz_poly_init(fT1)
    fmpz_poly_init(fT2)
    fmpz_poly_init(ft3)
    fmpz_poly_init(ft4)
    fmpz_poly_init(ft5)
    fmpz_poly_init(ft6)
    fmpz_poly_init(ft7)
    fmpz_poly_init(ft8)

    fmpz_poly_init(fs1)
    fmpz_poly_init(fs2)

    fmpz_poly_init(fU1)
    fmpz_poly_init(fU2)


# Initialize module-scope global C variables.
_init_globals()

def _clear_globals():
    """
    Clear all global variables allocated for optimization of
    quaternion algebra arithmetic.

    Do *not* call this except on exit of Sage, unless you want to see segfaults.

    This is called in the function quit_sage(), which is defined in all.py.

    EXAMPLES::

        sage: sage.algebras.quaternion_algebra_element._clear_globals()
        sage: sage.algebras.quaternion_algebra_element._init_globals()
    """
    mpz_clear(T1)
    mpz_clear(T2)
    mpz_clear(t3)
    mpz_clear(t4)
    mpz_clear(t5)
    mpz_clear(t6)
    mpz_clear(t7)
    mpz_clear(t8)

    mpz_clear(s1)
    mpz_clear(s2)

    mpz_clear(U1)
    mpz_clear(U2)

    fmpz_poly_clear(fT1)
    fmpz_poly_clear(fT2)
    fmpz_poly_clear(ft3)
    fmpz_poly_clear(ft4)
    fmpz_poly_clear(ft5)
    fmpz_poly_clear(ft6)
    fmpz_poly_clear(ft7)
    fmpz_poly_clear(ft8)

    fmpz_poly_clear(fs1)
    fmpz_poly_clear(fs2)

    fmpz_poly_clear(fU1)
    fmpz_poly_clear(fU2)

cdef to_quaternion(R, x):
    """
    Internal function used implicitly by quaternion algebra creation.

    INPUT::

        - R -- callable
        - x -- element or 4-tuple

    Given a callable R and an x that defines a quaternion, which can be a
    4-tuple, list of length 4, or something that coerces to R, return
    4-tuple of elements of R.

    EXAMPLES::
        sage: Q.<i,j,kkkk> = QuaternionAlgebra(QQ,-7, 13)
        sage: kkkk._repr_()   # implicit doctest
        'kkkk'
    """
    if isinstance(x, (list, tuple)):
        return R(x[0]), R(x[1]), R(x[2]), R(x[3])
    else:
        return R(x), R(0), R(0), R(0)

cdef inline print_coeff(y, i, bint atomic):
    """
    Internal function used implicitly by all quaternion algebra printing.

    INPUT::

        - y -- coefficient
        - i -- string (name of a generator)
        - atomic -- boolean int; whether or not elements of base ring
          print atomically

    EXAMPLES::
        sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-7, 13)
        sage: i._repr_()   # implicit doctest
        'i'
    """
    if not y:
        return ''
    if y == 1:
        return i
    elif y == -1:
        return "-%s"%i
    y = str(y)
    if not atomic and ('+' in y or '-' in y):
        return '(%s)*%s'%(y, i)
    else:
        return '%s*%s'%(y, i)

cdef class QuaternionAlgebraElement_abstract(AlgebraElement):
    cdef _do_print(self, x,y,z,w):
        """
        Used internally by the print function.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(-17,-19)
            sage: str(i+j+k-3/4)                            # indirect doctest
            '-3/4 + i + j + k'
        """
        cdef bint atomic = self._parent._base.is_atomic_repr()
        v = []
        i,j,k = self._parent.variable_names()
        if x:
            v.append(str(x))
        c = print_coeff(y,i,atomic)
        if c: v.append(c)
        c = print_coeff(z,j,atomic)
        if c: v.append(c)
        c = print_coeff(w,k,atomic)
        if c: v.append(c)
        if len(v) == 0: return '0'
        return ' + '.join(v).replace('+ -','- ')

    def _repr_(self):
        """
        Return string representation of this quaternion:

        EXAMPLES::
            sage: R.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(R,-5*x,-2)
            sage: a = x + i*x^3 + j*x^2 + k*x
            sage: a._repr_()
            'x + x^3*i + x^2*j + x*k'
            sage: a = x+2/3 + i*x^3 + j*(x^2-5/2) + k*x
            sage: a._repr_()
            'x + 2/3 + x^3*i + (x^2 - 5/2)*j + x*k'
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: Q(0)._repr_()
            '0'
        """
        return self._do_print(self[0], self[1], self[2], self[3])


    cdef int _cmp_c_impl(self, sage.structure.element.Element right) except -2:
        """
        Comparing elements.

        TESTS::

            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(QQ['x']),-5,-2)
            sage: a = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: a == a
            True
            sage: a == a + 1
            False
            sage: a < a + 1
            True

            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: (1/a)*a == 1
            True
            sage: (1/a)*a
            1024/1024
        """
        cdef int i
        for i in range(4):
            c = cmp(self[i], right[i])
            if c: return c
        return 0

    cpdef conjugate(self):
        """
        Return the conjugate of this quaternion.

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: a = 3*i - j + 2
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: a.conjugate()
            2 - 3*i + j

        A generic test::

            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a.conjugate()
            x + (-2*x)*i + (-3)*j + (-x + 2)*k
        """
        return self.__class__(self._parent, (self[0], -self[1], -self[2], -self[3]))

    cpdef trace(self):
        """
        Return the *reduced* trace of this quaternion.

        This is 2*x.

        EXAMPLES::
            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a.trace()
            2*x
        """
        return 2*self[0]

    def reduced_trace(self):
        """
        Return the reduced trace of self.

        This is an alias for self.trace().

        EXAMPLES::
            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a.reduced_trace()
            2*x
        """
        # do not overload this in a derived class!
        return self.trace()

    cpdef norm(self):
        """
        Return the *reduced* norm of this quaternion.

        This is w^2*a*b - y^2*a - z^2*b + x^2.

        EXAMPLES::
            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a.norm()
            2*x^4 - 12*x^3 + 9*x^2 - 18*x
        """
        a,b,x,y,z,w = self._parent._a,self._parent._b,self[0],self[1],self[2],self[3]
        return w*w*a*b - y*y*a - z*z*b + x*x

    def reduced_norm(self):
        """
        Return the reduced norm of self.

        This is an alias for self.norm().

        EXAMPLES::
            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a.reduced_norm()
            2*x^4 - 12*x^3 + 9*x^2 - 18*x
        """
        # do not overload this in a derived class!
        return self.norm()

    def __invert__(self):
        """
        Return inverse of self.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-7,-13)
            sage: a = 1/3 - 2/3*i + 4/19*j - 17/3*k
            sage: (1/a) * a
            1
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>

        Note that the quaternion algebra need not be a division
        algebra, in which case we can get a ZeroDivisionException::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,4,9)
            sage: a = 2-i
            sage: a.norm()
            0
            sage: 1/a
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        A test with a generic type:
            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(QQ['x,y']),4,9)
            sage: a = 1/3 - 2/3*i + 4/19*j - 17/3*k
            sage: (1/a) * a
            1
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: 1/(2-i)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert 0
        """
        return self.norm().__invert__() * self.conjugate()

    cpdef RingElement _div_(self, RingElement right):
        """
        Return quotient of self by right.

        EXAMPLES::
            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a._div_(a)   # funny reduced form of Frac(QQ['x'])...
            64/64
            sage: a._div_(a) == 1
            True
        """
        return self * right.__invert__()

cdef class QuaternionAlgebraElement_generic(QuaternionAlgebraElement_abstract):
    """
    TESTS::

    We test pickling::

        sage: R.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(R,-5*x,-2)
        sage: a = x + i*x^3 + j*x^2 + k*x
        sage: a == loads(dumps(a))
        True
    """
    def __init__(self, parent, v):
        """
        Create a quaternion over some general base field.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic(Q, (x,1,-7,2/3*x^3))
            x + i + (-7)*j + 2/3*x^3*k
        """
        self._parent = parent
        self.x, self.y, self.z, self.w = to_quaternion(parent._base, v)

    def __getitem__(self, int i):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(QQ['x']),-5,-2)
            sage: a = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: list(a)
            [1/2, 2/3, -3/4, 5/7]
        """
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        elif i == 3:
            return self.w
        else:
            raise IndexError, "quaternion element index out of range"

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::
            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: a = 1/x + x*i - (x+1)*j + 2/(3*x^3+5)*k
            sage: loads(dumps(a)) == a
            True
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
        """
        return (unpickle_QuaternionAlgebraElement_generic_v0,
                (self._parent, (self.x, self.y, self.z, self.w)))

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        Return the sum of self and _right.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: (x+i+j+x^3*k)._add_(x-i-j+ (2/3*x^3+x)*k)
            2*x + (5/3*x^3 + x)*k
            sage: type(i)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
        """
        cdef QuaternionAlgebraElement_generic right = _right
        # TODO -- make this, etc. use PY_NEW
        return QuaternionAlgebraElement_generic(self._parent, (self.x + right.x, self.y + right.y, self.z + right.z, self.w + right.w))

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        Return the difference of self and _right.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: type(i)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: (x+i+j+x^3*k)._sub_(x-i-j+ (2/3*x^3+x)*k)
            2*i + 2*j + (1/3*x^3 - x)*k
        """
        cdef QuaternionAlgebraElement_generic right = _right
        return QuaternionAlgebraElement_generic(self._parent, (self.x - right.x, self.y - right.y, self.z - right.z, self.w - right.w))

    cpdef RingElement _mul_(self, RingElement _right):
        """
        Return the product of self and _right.

        EXAMPLES::
            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: type(i)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: (x+i+j+x^3*k)._mul_(x-i-j+ (2/3*x^3+x)*k)
            -20/3*x^6 - 10*x^4 + x^2 + 7 + (10/3*x^3 + 2*x)*i + (-25/3*x^3 - 5*x)*j + (5/3*x^4 + x^2)*k
        """
        cdef QuaternionAlgebraElement_generic right = _right

        a = self._parent._a
        b = self._parent._b

        x1, y1, z1, w1 = self.x, self.y, self.z, self.w
        x2, y2, z2, w2 = right.x, right.y, right.z, right.w

        #x = x1*x2 + y1*y2*a + z1*z2*b - w1*w2*a*b
        #y = x1*y2 + y1*x2 - z1*w2*b + w1*z2*b
        #z = x1*z2 + y1*w2 + z1*x2 - w1*y2*a
        #w = x1*w2 + y1*z2 - z1*y2 + w1*x2

        t1 = x1 * x2
        t2 = y1 * y2
        t3 = z1 * z2
        t4 = w1 * w2
        t5 = x2 * z1
        t6 = y2 * w1
        t7 = x1 * z2
        t8 = y1 * w2

        x = t1 + a * t2 + b * (t3 - a*t4)
        y = (x1 + y1)*(x2 + y2) - t1 - t2 + b*( (z1 + w1)*(z2 - w2) - t3 + t4)
        z = t5 - a*t6 + t7 + a*t8
        w = (x2 - y2)*(z1 + w1) - t5 + t6 + (x1 + y1)*(z2 + w2) - t7 - t8

        return QuaternionAlgebraElement_generic(self._parent, (x, y, z, w))

    def _repr_(self):
        """
        Print representation of self.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: a = 1/x + x*i - (x+1)*j + 2/(3*x^3+5)*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a._repr_()
            '1/x + x*i + (-x - 1)*j + (2/(3*x^3 + 5))*k'
        """
        return self._do_print(self.x, self.y, self.z, self.w)


cdef class QuaternionAlgebraElement_rational_field(QuaternionAlgebraElement_abstract):
    """

    TESTS::

    We test pickling::

        sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
        sage: i + j + k == loads(dumps(i+j+k))
        True
    """

    # Implementation Notes:
    #
    # A Quaternion algebra element (call it a) over Q are implemented as a 4-tuple of
    # integers x, y, z, w and a demoninator d, all of type mpz_t, such that
    #
    #           a = (1/d) * (x + y * i + z * j + w * k)
    #
    # (although different letters may be specified instead of i, j, and k, if desired).
    #
    # Inside the element we also store mpz_t integers a and b, where
    #
    #       i^2 = a   and   j^2 = b

    def __cinit__(self):
        """
        Initialize C variables.

        EXAMPLES::
            sage: QuaternionAlgebra(QQ,-5,-2)([1/2,-1/3,2/3,4/5])  # implicit doctest
            1/2 - 1/3*i + 2/3*j + 4/5*k
        """
        mpz_init(self.x)
        mpz_init(self.y)
        mpz_init(self.z)
        mpz_init(self.w)
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.d)

    def  __dealloc__(self):
        mpz_clear(self.x)
        mpz_clear(self.y)
        mpz_clear(self.z)
        mpz_clear(self.w)
        mpz_clear(self.a)
        mpz_clear(self.b)
        mpz_clear(self.d)

    cdef int _cmp_c_impl(self, sage.structure.element.Element _right) except -2:
        """
        Compare two quaternions.  The comparison is fairly arbitrary
        -- first the denominators are compared and if equal then each
        of the other coefficients are compared.

        TESTS::
            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: i < j
            False
            sage: -i < j
            True
            sage: i == i
            True
        """
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef int i
        i = mpz_cmp(self.d, right.d)
        if i < 0: return -1
        elif i > 0: return 1
        i = mpz_cmp(self.x, right.x)
        if i < 0: return -1
        elif i > 0: return 1
        i = mpz_cmp(self.y, right.y)
        if i < 0: return -1
        elif i > 0: return 1
        i = mpz_cmp(self.z, right.z)
        if i < 0: return -1
        elif i > 0: return 1
        i = mpz_cmp(self.w, right.w)
        if i < 0: return -1
        elif i > 0: return 1
        return 0

    def __init__(self, parent, v):
        """
        Setup element data from parent and coordinates.

        EXAMPLES::

            sage: QuaternionAlgebra(QQ,-5,-2)([-1/2,-10/3,-2/3,-4/5])     # implicit doctest
            -1/2 - 10/3*i - 2/3*j - 4/5*k
        """
        self._parent = parent
        _x, _y, _z, _w = to_quaternion(Rational, v)

        # cache _a and _b
        _a = Integer(parent._a)
        _b = Integer(parent._b)

        common_denominator = lcm( [A.denominator() for A in (_x, _y, _z, _w)])
        xx, yy, zz, ww = [Integer(A * common_denominator) for A in (_x, _y, _z, _w)]

        mpz_set(self.x, (<Integer>xx).value)
        mpz_set(self.y, (<Integer>yy).value)
        mpz_set(self.z, (<Integer>zz).value)
        mpz_set(self.w, (<Integer>ww).value)
        mpz_set(self.a, (<Integer>_a).value)
        mpz_set(self.b, (<Integer>_b).value)
        mpz_set(self.d, (<Integer>common_denominator).value)

    def __getitem__(self, int i):
        """
        TESTS::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: a = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: list(a)
            [1/2, 2/3, -3/4, 5/7]
        """
        cdef Rational r = Rational()
        if i == 0:
            mpq_set_num(r.value, self.x)
        elif i == 1:
            mpq_set_num(r.value, self.y)
        elif i == 2:
            mpq_set_num(r.value, self.z)
        elif i == 3:
            mpq_set_num(r.value, self.w)
        else:
            raise IndexError, "quaternion element index out of range"
        mpq_set_den(r.value, self.d)
        mpq_canonicalize(r.value)
        return r

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: K.<x> = QQ[]
            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(K),-5,-19)
            sage: a = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: loads(dumps(a)) == a
            True

        """
        return (unpickle_QuaternionAlgebraElement_rational_field_v0,
                (self._parent, (self[0], self[1], self[2], self[3])))

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(15)
            sage: type(i)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: (2/3 + 3/4*i + 5/6*j + 7/8*k)._add_(-2/3 - 3/4*i + 5/6*j + 7/8*k)
            5/3*j + 7/4*k
        """

        #   Given two quaternion algebra elements
        #       a = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #       b = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #   we compute their sum as
        #
        #   (a + b) = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #   with    d3 = d1 * d2
        #           x3 = d1 * x2 + d2 * x1
        #           y3 = d1 * y2 + d2 * y1
        #           z3 = d1 * z2 + d2 * z1
        #           w3 = d1 * w2 + d2 * w1
        #
        #   and then we reduce the sum by dividing everything
        #   by the gcd of d3, x3, y3, z3, and w3

        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_mul(U1, self.x, right.d)        # U1 = x1 * d2
        mpz_mul(U2, right.x, self.d)        # U2 = x2 * d1
        mpz_add(result.x, U1, U2)           # x3 = x1 * d2 + x2 * d1

        mpz_mul(U1, self.y, right.d)        # U1 = y1 * d2
        mpz_mul(U2, right.y, self.d)        # U2 = y2 * d1
        mpz_add(result.y, U1, U2)           # x3 = y1 * d2 + y2 * d1

        mpz_mul(U1, self.z, right.d)        # U1 = z1 * d2
        mpz_mul(U2, right.z, self.d)        # U2 = z2 * d1
        mpz_add(result.z, U1, U2)           # z3 = z1 * d2 + z2 * d1

        mpz_mul(U1, self.w, right.d)        # U1 = w1 * d2
        mpz_mul(U2, right.w, self.d)        # U2 = w2 * d1
        mpz_add(result.w, U1, U2)           # w3 = w1 * d2 + w2 * d1

        mpz_mul(result.d, self.d, right.d) # d3 = d1 * d2

        result.canonicalize()

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        return result

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(15)
            sage: type(i)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: (2/3 + 3/4*i + 5/6*j + 7/8*k)._sub_(-2/3 - 3/4*i + 5/6*j + 7/8*k)
            4/3 + 3/2*i
        """
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        # Implementation Note: To obtain _sub_, we simply replace every occurrence of
        # "add" in _add_ with "sub"; that is, we s/add/sub to get _sub_

        mpz_mul(U1, self.x, right.d)
        mpz_mul(U2, right.x, self.d)
        mpz_sub(result.x, U1, U2)

        mpz_mul(U1, self.y, right.d)
        mpz_mul(U2, right.y, self.d)
        mpz_sub(result.y, U1, U2)

        mpz_mul(U1, self.z, right.d)
        mpz_mul(U2, right.z, self.d)
        mpz_sub(result.z, U1, U2)

        mpz_mul(U1, self.w, right.d)
        mpz_mul(U2, right.w, self.d)
        mpz_sub(result.w, U1, U2)

        mpz_mul(result.d, self.d, right.d)

        result.canonicalize()

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        return result




    cpdef RingElement _mul_(self, RingElement _right):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(15)
            sage: type(i)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: (2/3 + 3/4*i + 5/6*j + 7/8*k)._mul_(-2/3 - 3/4*i + 5/6*j + 7/8*k)
            9331/576 - i - 63/16*j + 5/4*k
        """

        # We use the following formula for multiplication:
        #
        #    Given two quaternion algebra elements
        #
        #        a = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #        b = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #    we compute their product as
        #
        #    ab = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #    where
        #       x3 = t1 + a * t2 + b * (t3 - a*t4)
        #       y3 = s1*(x2 + y2) - t1 - t2 + b*( s2*(z2 - w2) - t3 + t4)
        #       z3 = t5 - a*t6 + t7 + a*t8
        #       w3 = (x2 - y2)*s2 - t5 + t6 + s1*(z2 + w2) - t7 - t8
        #
        #       and where
        #           t1 = x1 * x2
        #           t2 = y1 * y2
        #           t3 = z1 * z2
        #           t4 = w1 * w2
        #           t5 = x2 * z1
        #           t6 = y2 * w1
        #           t7 = x1 * z2
        #           t8 = y1 * w2
        #
        #           s1 = x1 + y1
        #           s2 = z1 + w1
        #
        # This takes more integer addition operations but fewer integer multiplication
        # operations than the "straightforward" multiplication method.
        #
        # There might be a way to optimize this formula further.


        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)

        mpz_mul(T1, self.x, right.x)    # t1 = x1 * x2
        mpz_mul(T2, self.y, right.y)    # t2 = y1 * y2
        mpz_mul(t3, self.z, right.z)    # t3 = z1 * z2
        mpz_mul(t4, self.w, right.w)    # t4 = w1 * w2
        mpz_mul(t5, right.x, self.z)    # t5 = x2 * z1
        mpz_mul(t6, right.y, self.w)    # t6 = y2 * w1
        mpz_mul(t7, self.x, right.z)    # t7 = x1 * z2
        mpz_mul(t8, self.y, right.w)    # t8 = y1 * w2

        mpz_add(s1, self.x, self.y)     # s1 = x1 + y1
        mpz_add(s2, self.z, self.w)     # s2 = z1 + w1

        #------------------

        mpz_mul(U1, self.a, t4)
        mpz_sub(U1, t3, U1)
        mpz_mul(U1, U1, self.b)
        mpz_mul(U2, self.a, T2)
        mpz_add(result.x, T1, U2)
        mpz_add(result.x, result.x, U1)

        #------------------

        mpz_sub(U1, right.z, right.w)
        mpz_mul(U1, U1, s2)
        mpz_sub(U1, U1, t3)
        mpz_add(U1, U1, t4)
        mpz_mul(U1, U1, self.b)
        mpz_sub(U1, U1, T2)
        mpz_sub(U1, U1, T1)
        mpz_add(U2, right.x, right.y)
        mpz_mul(U2, s1, U2)
        mpz_add(result.y, U1, U2)

        #------------------

        mpz_mul(U1, self.a, t8)
        mpz_add(U1, U1, t7)
        mpz_mul(U2, self.a, t6)
        mpz_sub(U1, U1, U2)
        mpz_add(result.z, U1, t5)

        #------------------

        mpz_add(U1, right.z, right.w)
        mpz_mul(U1, U1, s1)
        mpz_sub(U1, U1, t7)
        mpz_sub(U1, U1, t8)
        mpz_add(U1, U1, t6)
        mpz_sub(U1, U1, t5)
        mpz_sub(U2, right.x, right.y)
        mpz_mul(U2, U2, s2)
        mpz_add(result.w, U1, U2)



        mpz_mul(result.d, self.d, right.d)

        result.canonicalize()

        return result

    cpdef norm(self):
        """
        Return the *reduced* norm of this quaternion.

        This is w^2*a*b - y^2*a - z^2*b + x^2.

        EXAMPLES::
            sage: K.<i,j,k> = QuaternionAlgebra(QQ, -5, -2)
            sage: i.norm()
            5
            sage: j.norm()
            2
            sage: a = 1/3 + 1/5*i + 1/7*j + k
            sage: a.norm()
            22826/2205
        """
        mpz_mul(U1, self.x, self.x)         # U1 = x*x
        mpz_mul(U2, self.b, self.z)         # U2 = b*x
        mpz_mul(U2, U2, self.z)             # U2 = b*z*z
        mpz_sub(U2, U1, U2)                 # U2 = -b*z*z + x*x

        mpz_mul(U1, self.y, self.a)         # U1 = a*y
        mpz_mul(U1, U1, self.y)             # U1 = a*y*y
        mpz_sub(U2, U2, U1)                 # U2 = -a*y*y - b*z*z + x*x

        mpz_mul(U1, self.w, self.w)         # U1 = w*w
        mpz_mul(U1, U1, self.a)             # U1 = w*w*a
        mpz_mul(U1, U1, self.b)             # U1 = w*w*a*b

        mpz_add(U1, U1, U2)                 # U1 = w*w*a*b - a*y*y - b*z*z + x*x

        mpz_mul(U2, self.d, self.d)

        cdef Rational result = PY_NEW(Rational)
        mpq_set_num(result.value, U1)
        mpq_set_den(result.value, U2)
        mpq_canonicalize(result.value)

        return result

    cpdef conjugate(self):
        """
        Return the conjugate of this quaternion.

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: a = 3*i - j + 2
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: a.conjugate()
            2 - 3*i + j
            sage: b = 1 + 1/3*i + 1/5*j - 1/7*k
            sage: b.conjugate()
            1 - 1/3*i - 1/5*j + 1/7*k
        """

#        return self.__class__(self._parent, (self[0], -self[1], -self[2], -self[3]))
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        mpz_set(result.d, self.d)

        mpz_set(result.x, self.x)
        mpz_mul_si(result.y, self.y, -1)
        mpz_mul_si(result.z, self.z, -1)
        mpz_mul_si(result.w, self.w, -1)

        return result

    cpdef trace(self):
        """
        Return the *reduced* trace of this quaternion.

        This is 2 * x.

        EXAMPLES::
            sage: K.<i,j,k> = QuaternionAlgebra(QQ, -5, -2)
            sage: i.trace()
            0
            sage: j.trace()
            0
            sage: a = 1/3 + 1/5*i + 1/7*j + k
            sage: a.trace()
            2/3
        """
        #return 2*self[0]

        mpz_mul_si(U1, self.x, 2)
        cdef Rational result = PY_NEW(Rational)
        mpq_set_num(result.value, U1)
        mpq_set_den(result.value, self.d)
        mpq_canonicalize(result.value)
        return result

    cdef inline canonicalize(self):
        """
        Put the representation of this quaternion element into
        smallest form. For a = (1/d)*(x + yi + zj + wk) we
        divide a, x, y, z, and w by the gcd of all of them.

        TESTS::
            sage: K.<i,j,k> = QuaternionAlgebra(QQ, -10, -7)
            sage: (1/4 + 1/2 * i + 1/7 * j + 1/28 * k)*14*i     # implicit doctest
            -70 + 7/2*i + 5*j - 2*k
        """

        # Note: this function changes the module-level global variable
        # U1, so it isn't always safe to use this in the middle of
        # another function. Normally this function is called
        # at the end of an arithmetic routine, so this is fine.

        # Implemenationwise, we compute the gcd's one at a time,
        # and quit if it ever becomes one


        mpz_gcd(U1, self.d, self.x)
        if mpz_cmp_ui(U1, 1) != 0:
            mpz_gcd(U1, U1, self.y)
            if mpz_cmp_ui(U1, 1) != 0:
                mpz_gcd(U1, U1, self.z)
                if mpz_cmp_ui(U1, 1) != 0:
                    mpz_gcd(U1, U1, self.w)
                    if mpz_cmp_ui(U1, 1) != 0:
                        # at this point U1 actually contains the gcd of all the terms, and we divide
                        mpz_divexact(self.d, self.d, U1)
                        mpz_divexact(self.x, self.x, U1)
                        mpz_divexact(self.y, self.y, U1)
                        mpz_divexact(self.z, self.z, U1)
                        mpz_divexact(self.w, self.w, U1)



cdef class QuaternionAlgebraElement_number_field(QuaternionAlgebraElement_abstract):
    def __new__(self):
        """
        Allocate memory for this quaternion over a number field.
        """
        fmpz_poly_init(self.x)
        fmpz_poly_init(self.y)
        fmpz_poly_init(self.z)
        fmpz_poly_init(self.w)
        fmpz_poly_init(self.a)
        fmpz_poly_init(self.b)
        fmpz_poly_init(self.modulus)
        mpz_init(self.d)

    def __dealloc__(self):
        """
        Free memory used by this quaternion over a number field.
        """
        fmpz_poly_clear(self.x)
        fmpz_poly_clear(self.y)
        fmpz_poly_clear(self.z)
        fmpz_poly_clear(self.w)
        fmpz_poly_clear(self.a)
        fmpz_poly_clear(self.b)
        fmpz_poly_clear(self.modulus)
        mpz_clear(self.d)

    def __init__(self, parent, v):
        """
        EXAMPLES::

            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K,-a,a+1)
            sage: Q([a,-2/3,a^2-1/2,a*2])           # implicit doctest
            a + (-2/3)*i + (a^2 - 1/2)*j + 2*a*k

        """
        self._parent = parent
        x, y, z, w = to_quaternion(parent._base, v)
        cdef NumberFieldElement a = <NumberFieldElement>(parent._base(parent._a))
        cdef NumberFieldElement b = <NumberFieldElement>(parent._base(parent._b))
        ZZX_to_fmpz_poly(self.x, (<NumberFieldElement>x).__numerator)
        ZZX_to_fmpz_poly(self.y, (<NumberFieldElement>y).__numerator)
        ZZX_to_fmpz_poly(self.z, (<NumberFieldElement>z).__numerator)
        ZZX_to_fmpz_poly(self.w, (<NumberFieldElement>w).__numerator)

        ZZ_to_mpz(&T1, &(<NumberFieldElement>x).__denominator)
        ZZ_to_mpz(&T2, &(<NumberFieldElement>y).__denominator)
        ZZ_to_mpz(&t3, &(<NumberFieldElement>z).__denominator)
        ZZ_to_mpz(&t4, &(<NumberFieldElement>w).__denominator)

        mpz_lcm(self.d, T1, T2)
        mpz_lcm(self.d, self.d, t3)
        mpz_lcm(self.d, self.d, t4)

        mpz_divexact(T1, self.d, T1)
        mpz_divexact(T2, self.d, T2)
        mpz_divexact(t3, self.d, t3)
        mpz_divexact(t4, self.d, t4)

        fmpz_poly_scalar_mul_mpz(self.x, self.x, T1)
        fmpz_poly_scalar_mul_mpz(self.y, self.y, T2)
        fmpz_poly_scalar_mul_mpz(self.z, self.z, t3)
        fmpz_poly_scalar_mul_mpz(self.w, self.w, t4)

        ZZX_to_fmpz_poly(self.a, a.__numerator)     # we will assume that the denominafor of a and b are 1
        ZZX_to_fmpz_poly(self.b, b.__numerator)

        ZZX_to_fmpz_poly(self.modulus, (<NumberFieldElement>x).__fld_numerator.x) # and same for the modulus

    def __getitem__(self, int i):
        """
        EXAMPLES::

            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K,-a,a+1)
            sage: Q([a,-2/3,a^2-1/2,a*2])
            a + (-2/3)*i + (a^2 - 1/2)*j + 2*a*k
            sage: x = Q([a,-2/3,a^2-1/2,a*2])
            sage: type(x)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: x[0]
            a
            sage: x[1]
            -2/3
            sage: x[2]
            a^2 - 1/2
            sage: x[3]
            2*a
            sage: list(x)
            [a, -2/3, a^2 - 1/2, 2*a]
        """
        # general number -- this code assumes that the number field is not quadratic!!

        cdef NumberFieldElement el = <NumberFieldElement>(self._parent.base_ring().an_element())
        cdef NumberFieldElement item = el._new()

        if i == 0:
            fmpz_poly_to_ZZX(item.__numerator, self.x)
        elif i == 1:
            fmpz_poly_to_ZZX(item.__numerator, self.y)
        elif i == 2:
            fmpz_poly_to_ZZX(item.__numerator, self.z)
        elif i == 3:
            fmpz_poly_to_ZZX(item.__numerator, self.w)
        else:
            raise IndexError, "quaternion element index out of range"

        mpz_to_ZZ(&item.__denominator, &self.d)

        return item


    def __reduce__(self):
        """
        EXAMPLES::

            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = (i+j+k+a)^2; z
            a^2 + 4*a - 3 + 2*a*i + 2*a*j + 2*a*k
            sage: f, t = z.__reduce__()
            sage: f(*t)
            a^2 + 4*a - 3 + 2*a*i + 2*a*j + 2*a*k
            sage: loads(dumps(z)) == z
            True
        """
        return (unpickle_QuaternionAlgebraElement_number_field_v0,
                (self._parent, (self[0], self[1], self[2], self[3])))

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        Add self and _right:

        EXAMPLES::
            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = a + i + (2/3)*a^3*j + (1+a)*k; w = a - i - (2/3)*a^3*j + (1/3+a)*k
            sage: type(z)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: z._add_(w)
            2*a + (2*a + 4/3)*k
        """

        #   Given two quaternion algebra elements
        #       a = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #       b = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #   we compute their sum as
        #
        #   (a + b) = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #   with    d3 = d1 * d2
        #           x3 = d1 * x2 + d2 * x1
        #           y3 = d1 * y2 + d2 * y1
        #           z3 = d1 * z2 + d2 * z1
        #           w3 = d1 * w2 + d2 * w1
        #
        #   and then we reduce the sum by calling canonicalize().

        # Note: We are assuming in this routine that the modulus is monic. This shouldn't
        # currently be an issue because it is impossible to create a number field with
        # a modulus that is not monic.

        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> PY_NEW(QuaternionAlgebraElement_number_field)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)
        result._parent = self._parent

        fmpz_poly_scalar_mul_mpz(fU1, self.x, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.x, self.d)
        fmpz_poly_add(result.x, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.y, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.y, self.d)
        fmpz_poly_add(result.y, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.w, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.w, self.d)
        fmpz_poly_add(result.w, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.z, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.z, self.d)
        fmpz_poly_add(result.z, fU1, fU2)

        mpz_mul(result.d, self.d, right.d)

        self.canonicalize()

        return result





    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        Subtract _right from self.

        EXAMPLES::

            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = a + i + (2/3)*a^3*j + (1+a)*k; w = a - i - (2/3)*a^3*j + (1/3+a)*k
            sage: type(z)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: z._sub_(w)
            2*i + 8/3*j + 2/3*k

        """
        # Implementation Note: To obtain _sub_, we simply replace every occurrence of
        # "add" in _add_ with "sub"; that is, we s/add/sub to get _sub_

        # Note: We are assuming in this routine that the modulus is monic. This shouldn't
        # currently be an issue because it is impossible to create a number field with
        # a modulus that is not monic.
        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> PY_NEW(QuaternionAlgebraElement_number_field)


        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)
        result._parent = self._parent

        fmpz_poly_scalar_mul_mpz(fU1, self.x, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.x, self.d)
        fmpz_poly_sub(result.x, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.y, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.y, self.d)
        fmpz_poly_sub(result.y, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.w, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.w, self.d)
        fmpz_poly_sub(result.w, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.z, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.z, self.d)
        fmpz_poly_sub(result.z, fU1, fU2)

        mpz_mul(result.d, self.d, right.d)

        self.canonicalize()

        return result

    cpdef RingElement _mul_(self, RingElement _right):
        """
        Multiply self and _right.

        EXAMPLES::

            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = a + i + (2/3)*a^3*j + (1+a)*k; w = a - i - (2/3)*a^3*j + (1/3+a)*k
            sage: type(z)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: z._mul_(w)
            5*a^2 - 7/9*a + 9 + (-8/3*a^2 - 16/9*a)*i + (-6*a - 4)*j + (2*a^2 + 4/3*a)*k
        """

        # We use the following formula for multiplication:
        #
        #    Given two quaternion algebra elements
        #
        #        a = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #        b = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #    we compute their product as
        #
        #    ab = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #    where
        #       x3 = t1 + a * t2 + b * (t3 - a*t4)
        #       y3 = s1*(x2 + y2) - t1 - t2 + b*( s2*(z2 - w2) - t3 + t4)
        #       z3 = t5 - a*t6 + t7 + a*t8
        #       w3 = (x2 - y2)*s2 - t5 + t6 + s1*(z2 + w2) - t7 - t8
        #
        #       and where
        #           t1 = x1 * x2
        #           t2 = y1 * y2
        #           t3 = z1 * z2
        #           t4 = w1 * w2
        #           t5 = x2 * z1
        #           t6 = y2 * w1
        #           t7 = x1 * z2
        #           t8 = y1 * w2
        #
        #           s1 = x1 + y1
        #           s2 = z1 + w1
        #
        # This takes more polynomial addition operations but fewer polynomial multiplication
        # operations than the "straightforward" multiplication method.
        #
        # There might be a way to optimize this formula further.

        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> PY_NEW(QuaternionAlgebraElement_number_field)

        mpz_set_si(result.d, 1)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)
        result._parent = self._parent

        fmpz_poly_mul(fT1, self.x, right.x)    # t1 = x1 * x2
        fmpz_poly_mul(fT2, self.y, right.y)    # t2 = y1 * y2
        fmpz_poly_mul(ft3, self.z, right.z)    # t3 = x1 * x2
        fmpz_poly_mul(ft4, self.w, right.w)    # t4 = w1 * w2
        fmpz_poly_mul(ft5, right.x, self.z)    # t5 = x2 * z1
        fmpz_poly_mul(ft6, right.y, self.w)    # t6 = y2 * w1
        fmpz_poly_mul(ft7, self.x, right.z)    # t7 = x1 * z2
        fmpz_poly_mul(ft8, self.y, right.w)    # t8 = y1 * w2

        fmpz_poly_add(fs1, self.x, self.y)     # s1 = x1 + y1
        fmpz_poly_add(fs2, self.z, self.w)     # s2 = z1 + w

        #------------------

        fmpz_poly_mul(fU1, self.a, ft4)
        fmpz_poly_sub(fU1, ft3, fU1)
        fmpz_poly_mul(fU1, fU1, self.b)
        fmpz_poly_mul(fU2, self.a, fT2)
        fmpz_poly_add(result.x, fT1, fU2)
        fmpz_poly_add(result.x, result.x, fU1)

        #------------------

        fmpz_poly_sub(fU1, right.z, right.w)
        fmpz_poly_mul(fU1, fU1, fs2)
        fmpz_poly_sub(fU1, fU1, ft3)
        fmpz_poly_add(fU1, fU1, ft4)
        fmpz_poly_mul(fU1, fU1, self.b)
        fmpz_poly_sub(fU1, fU1, fT2)
        fmpz_poly_sub(fU1, fU1, fT1)
        fmpz_poly_add(fU2, right.x, right.y)
        fmpz_poly_mul(fU2, fs1, fU2)
        fmpz_poly_add(result.y, fU1, fU2)

        #------------------

        fmpz_poly_mul(fU1, self.a, ft8)
        fmpz_poly_add(fU1, fU1, ft7)
        fmpz_poly_mul(fU2, self.a, ft6)
        fmpz_poly_sub(fU1, fU1, fU2)
        fmpz_poly_add(result.z, fU1, ft5)

        #------------------

        fmpz_poly_add(fU1, right.z, right.w)
        fmpz_poly_mul(fU1, fU1, fs1)
        fmpz_poly_sub(fU1, fU1, ft7)
        fmpz_poly_sub(fU1, fU1, ft8)
        fmpz_poly_add(fU1, fU1, ft6)
        fmpz_poly_sub(fU1, fU1, ft5)
        fmpz_poly_sub(fU2, right.x, right.y)
        fmpz_poly_mul(fU2, fU2, fs2)
        fmpz_poly_add(result.w, fU1, fU2)

        # At this point we have essentially computed the product, but we still
        # need to reduce modulo the modulus, which is what the following 12 lines do.
        #
        # When this was written, the version of flint in Sage had problems with
        # fpmz_poly_divrem(). This should be fixed in the newest version of
        # flint, which also should have some new functions which should do
        # this faster (Bill Hart sent an email to Bober and William about this).
        #
        # This should be fixed in the near future. (I don't know how much
        # faster it will be when it is updated, but the following code is
        # currently quite a bottleneck.

        fmpz_poly_div(fT1, result.x, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.x, result.x, fT1)

        fmpz_poly_div(fT1, result.y, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.y, result.y, fT1)

        fmpz_poly_div(fT1, result.z, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.z, result.z, fT1)

        fmpz_poly_div(fT1, result.w, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.w, result.w, fT1)

        mpz_mul(result.d, self.d, right.d)

        self.canonicalize()

        return result

    cdef inline canonicalize(self):
        """
        Put the representation of this quaternion element into
        smallest form. For a = (1/d)*(x + yi + zj + wk) we
        divide a, x, y, z, and w by the gcd of all of them.

        TESTS::
            sage: F = QQ[3^(1/3)]
            sage: a = F.gen()
            sage: K.<i,j,k> = QuaternionAlgebra(F, -10 + a, -7 - a)
            sage: ((1/4 + 1/2 * i + a^3/7 * j + a/28 * k)*14*i)^3   # implicit doctest
            34503/2*a^2 + 132195/2*a + 791399/4 + (203/8*a^2 - 10591*a + 169225/4)*i + (-84695/4*a^2 + 483413/8*a + 18591/4)*j + (-87/2*a^2 + 18156*a - 72525)*k
        """

        # Note: this function changes the module-level global variables
        # U1 and U2, so it isn't always safe to use this in the middle of
        # another function. Normally this function is called
        # at the end of an arithmetic routine, so this is fine.

        # Implemenationwise, we compute the gcd's one at a time,
        # and quit if it ever becomes one

        cdef fmpz_t content = fmpz_init(fmpz_poly_max_limbs(self.x)) # TODO: think about how big this should be (probably the size of d)
                                                                     # Note that we have to allocate this here, and not
                                                                     # as a global variable, because fmpz_t's do not
                                                                     # self allocate memory
        fmpz_poly_content(content, self.x)
        fmpz_to_mpz(U1, content)
        mpz_gcd(U1, self.d, U1)
        if mpz_cmp_ui(U1, 1) != 0:
            fmpz_poly_content(content, self.y)
            fmpz_to_mpz(U2, content)
            mpz_gcd(U1, U1, U2)
            if mpz_cmp_ui(U1, 1) != 0:
                fmpz_poly_content(content, self.z)
                fmpz_to_mpz(U2, content)
                mpz_gcd(U1, U1, U2)
                if mpz_cmp_ui(U1, 1) != 0:
                    fmpz_poly_content(content, self.w)
                    fmpz_to_mpz(U2, content)
                    mpz_gcd(U1, U1, U2)
                    if mpz_cmp_ui(U1, 1) != 0:
                        fmpz_poly_scalar_div_mpz(self.x, self.x, U1)
                        fmpz_poly_scalar_div_mpz(self.y, self.y, U1)
                        fmpz_poly_scalar_div_mpz(self.z, self.z, U1)
                        fmpz_poly_scalar_div_mpz(self.w, self.w, U1)
                        mpz_divexact(self.d, self.d, U1)

        fmpz_clear(content)




#######################################################################
# Versioned unpickle functions
#######################################################################

def unpickle_QuaternionAlgebraElement_generic_v0(*args):
    """
    EXAMPLES::
        sage: K.<X> = QQ[]
        sage: Q.<i,j,k> = QuaternionAlgebra(Frac(K), -5,-19); z = 2/3 + i*X - X^2*j + X^3*k
        sage: f, t = z.__reduce__()
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_generic_v0(*t)
        2/3 + X*i + (-X^2)*j + X^3*k
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_generic_v0(*t) == z
        True
    """
    return QuaternionAlgebraElement_generic(*args)

def unpickle_QuaternionAlgebraElement_rational_field_v0(*args):
    """
    EXAMPLES::

        sage: Q.<i,j,k> = QuaternionAlgebra(-5,-19); a = 2/3 + i*5/7 - j*2/5 +19/2
        sage: f, t = a.__reduce__()
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_rational_field_v0(*t)
        61/6 + 5/7*i - 2/5*j
    """
    return QuaternionAlgebraElement_rational_field(*args)

def unpickle_QuaternionAlgebraElement_number_field_v0(*args):
    """
    EXAMPLES::

        sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a); z = i + j
        sage: f, t = z.__reduce__()
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_number_field_v0(*t)
        i + j
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_number_field_v0(*t) == z
        True
    """
    return QuaternionAlgebraElement_number_field(*args)
