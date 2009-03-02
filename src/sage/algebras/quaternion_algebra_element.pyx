from sage.structure.element cimport AlgebraElement, RingElement, ModuleElement, Element
from sage.algebras.quaternion_algebra_element cimport QuaternionAlgebraElement_abstract
from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.arith import lcm
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

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
        return ' + '.join(v).replace('+ -','- ')

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

    def conjugate(self):
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

        This is w^2*a*b - y^2*a - z^2*b + x^2.

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
            sage: Q.<i,j,k> = QuaternionAlgebra(QQ['x,y'],4,9)
            sage: a = 1/3 - 2/3*i + 4/19*j - 17/3*k
            sage: (1/a) * a
            1
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: 1/(2-i)
            Traceback (most recent call last):
            ...
            ZeroDivisionError
        """
        return self.norm().__invert__() * self.conjugate()

    cpdef RingElement _div_(self, RingElement right):
        """
        Return quotient of self by right.

        EXAMPLES::

        """
        return self * right.__invert__()

cdef class QuaternionAlgebraElement_generic(QuaternionAlgebraElement_abstract):
    def __init__(self, parent, v):
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
        # TODO -- make this versioned!!!
        return QuaternionAlgebraElement_generic, (self._parent, (self.x, self.y, self.z, self.w))

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        Return the sum of self and _right.

        EXAMPLES:
            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-7,-11)
            sage: i + k
            i + k
            sage: 1/2 + i + j
            1/2 + i + j
        """
        cdef QuaternionAlgebraElement_generic right = _right
        # TODO -- make this, etc. use PY_NEW
        return QuaternionAlgebraElement_generic(self._parent, (self.x + right.x, self.y + right.y, self.z + right.z, self.w + right.w))

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef QuaternionAlgebraElement_generic right = _right
        return QuaternionAlgebraElement_generic(self._parent, (self.x - right.x, self.y - right.y, self.z - right.z, self.w - right.w))

    cpdef RingElement _mul_(self, RingElement _right):
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
    def __cinit__(self):
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

    def __cmp__(self, _right):
        """
        TESTS::

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
        return 0

    def __init__(self, parent, v):
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
            sage: Q.<i,j,k> = QuaternionAlgebra(K,-5,-19)
            sage: a = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: loads(dumps(a)) == a
            True

        """
        # TODO -- make this versioned!!!
        return QuaternionAlgebraElement_rational_field, (self._parent, (self[0], self[1], self[2], self[3]))

    cpdef ModuleElement _add_(self, ModuleElement _right):
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_mul(U1, self.x, right.d)
        mpz_mul(U2, right.x, self.d)
        mpz_add(result.x, U1, U2)

        mpz_mul(U1, self.y, right.d)
        mpz_mul(U2, right.y, self.d)
        mpz_add(result.y, U1, U2)

        mpz_mul(U1, self.z, right.d)
        mpz_mul(U2, right.z, self.d)
        mpz_add(result.z, U1, U2)

        mpz_mul(U1, self.w, right.d)
        mpz_mul(U2, right.w, self.d)
        mpz_add(result.w, U1, U2)

        mpz_mul(result.d, self.d, right.d)

        mpz_gcd(U1, result.d, result.x)
        if mpz_cmp_ui(U1, 1) != 0:
            mpz_gcd(U1, U1, result.y)
            if mpz_cmp_ui(U1, 1) != 0:
                mpz_gcd(U1, U1, result.z)
                if mpz_cmp_ui(U1, 1) != 0:
                    mpz_gcd(U1, U1, result.w)
                    if mpz_cmp_ui(U1, 1) != 0:
                        mpz_divexact(result.d, result.d, U1)
                        mpz_divexact(result.x, result.x, U1)
                        mpz_divexact(result.y, result.y, U1)
                        mpz_divexact(result.z, result.z, U1)
                        mpz_divexact(result.w, result.w, U1)

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        return result

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

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

        mpz_gcd(U1, result.d, result.x)
        if mpz_cmp_ui(U1, 1) != 0:
            mpz_gcd(U1, U1, result.y)
            if mpz_cmp_ui(U1, 1) != 0:
                mpz_gcd(U1, U1, result.z)
                if mpz_cmp_ui(U1, 1) != 0:
                    mpz_gcd(U1, U1, result.w)
                    if mpz_cmp_ui(U1, 1) != 0:
                        mpz_divexact(result.d, result.d, U1)
                        mpz_divexact(result.x, result.x, U1)
                        mpz_divexact(result.y, result.y, U1)
                        mpz_divexact(result.z, result.z, U1)
                        mpz_divexact(result.w, result.w, U1)




        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        return result




    def _repr_(_self):
        cdef QuaternionAlgebraElement_rational_field self = _self
        cdef Integer _x = Integer()
        cdef Integer _y = Integer()
        cdef Integer _z = Integer()
        cdef Integer _w = Integer()
        cdef Integer d = Integer()

        mpz_set(_x.value, self.x)
        mpz_set(_y.value, self.y)
        mpz_set(_z.value, self.z)
        mpz_set(_w.value, self.w)
        mpz_set(d.value, self.d)

        x = _x/d
        y = _y/d
        z = _z/d
        w = _w/d

        i,j,k = 'i', 'j', 'k'
        v = []
        return self._do_print(x,y,z,w)


    cpdef RingElement _mul_(self, RingElement _right):
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

        mpz_gcd(U1, result.d, result.x)
        if mpz_cmp_ui(U1, 1) != 0:
            mpz_gcd(U1, U1, result.y)
            if mpz_cmp_ui(U1, 1) != 0:
                mpz_gcd(U1, U1, result.z)
                if mpz_cmp_ui(U1, 1) != 0:
                    mpz_gcd(U1, U1, result.w)
                    if mpz_cmp_ui(U1, 1) != 0:
                        mpz_divexact(result.d, result.d, U1)
                        mpz_divexact(result.x, result.x, U1)
                        mpz_divexact(result.y, result.y, U1)
                        mpz_divexact(result.z, result.z, U1)
                        mpz_divexact(result.w, result.w, U1)

        return result




cdef class QuaternionAlgebraElement_number_field(QuaternionAlgebraElement_abstract):
    def __new__(self):
        fmpz_poly_init(self.x)
        fmpz_poly_init(self.y)
        fmpz_poly_init(self.z)
        fmpz_poly_init(self.w)
        fmpz_poly_init(self.a)
        fmpz_poly_init(self.b)
        fmpz_poly_init(self.modulus)
        mpz_init(self.d)

    def __dealloc__(self):
        fmpz_poly_clear(self.x)
        fmpz_poly_clear(self.y)
        fmpz_poly_clear(self.z)
        fmpz_poly_clear(self.w)
        fmpz_poly_clear(self.a)
        fmpz_poly_clear(self.b)
        fmpz_poly_clear(self.modulus)
        mpz_clear(self.d)

    def __init__(self, parent, v):
        self._parent = parent
        x, y, z, w = to_quaternion(parent._base, v)
##         fmpz_poly_set_from_number_field_element(self.x, d_x, x)
##         fmpz_poly_set_from_number_field_element(self.y, d_y, y)
##         fmpz_poly_set_from_number_field_element(self.z, d_z, z)
##         fmpz_poly_set_from_number_field_element(self.w, d_w, w)
##         fmpz_poly_set_from_number_field_element(self.a, parent._a)
##         fmpz_poly_set_from_number_field_element(self.b, parent._b)

##         #TODO: #fmpz_poly_set_from_number_field_element(self.b, parent.modulus())

##         mpz_set(self.d, (<Integer>_d).value)


    cpdef ModuleElement _add_(self, ModuleElement _right):
        # Note: We are assuming in the routine that the modulus is monic. If it isn't there could be some problems.
        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> PY_NEW(QuaternionAlgebraElement_number_field)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)

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

        cdef fmpz_t content = fmpz_init(fmpz_poly_max_limbs(result.x)) # TODO: think about how big this should be (probably the size of d)

        fmpz_poly_content(content, result.x)
        fmpz_to_mpz(U1, content)
        mpz_gcd(U1, result.d, U1)
        if mpz_cmp_ui(U1, 1) != 0:
            fmpz_poly_content(content, result.y)
            fmpz_to_mpz(U2, content)
            mpz_gcd(U1, U1, U2)
            if mpz_cmp_ui(U1, 1) != 0:
                fmpz_poly_content(content, result.z)
                fmpz_to_mpz(U2, content)
                mpz_gcd(U1, U1, U2)
                if mpz_cmp_ui(U1, 1) != 0:
                    fmpz_poly_content(content, result.w)
                    fmpz_to_mpz(U2, content)
                    mpz_gcd(U1, U1, U2)
                    if mpz_cmp_ui(U1, 1) != 0:
                        fmpz_poly_scalar_div_mpz(result.x, result.x, U1)
                        fmpz_poly_scalar_div_mpz(result.y, result.y, U1)
                        fmpz_poly_scalar_div_mpz(result.z, result.z, U1)
                        fmpz_poly_scalar_div_mpz(result.w, result.w, U1)
                        mpz_divexact(result.d, result.d, U1)

        fmpz_clear(content)

        return result





    cpdef ModuleElement _sub_(self, ModuleElement _right):
        # Note: We are assuming in the routine that the modulus is monic. If it isn't there could be some problems.
        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> PY_NEW(QuaternionAlgebraElement_number_field)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)

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

        cdef fmpz_t content = fmpz_init(fmpz_poly_max_limbs(result.x)) # TODO: think about how big this should be (probably the size of d)

        fmpz_poly_content(content, result.x)
        fmpz_to_mpz(U1, content)
        mpz_gcd(U1, result.d, U1)
        if mpz_cmp_ui(U1, 1) != 0:
            fmpz_poly_content(content, result.y)
            fmpz_to_mpz(U2, content)
            mpz_gcd(U1, U1, U2)
            if mpz_cmp_ui(U1, 1) != 0:
                fmpz_poly_content(content, result.z)
                fmpz_to_mpz(U2, content)
                mpz_gcd(U1, U1, U2)
                if mpz_cmp_ui(U1, 1) != 0:
                    fmpz_poly_content(content, result.w)
                    fmpz_to_mpz(U2, content)
                    mpz_gcd(U1, U1, U2)
                    if mpz_cmp_ui(U1, 1) != 0:
                        fmpz_poly_scalar_div_mpz(result.x, result.x, U1)
                        fmpz_poly_scalar_div_mpz(result.y, result.y, U1)
                        fmpz_poly_scalar_div_mpz(result.z, result.z, U1)
                        fmpz_poly_scalar_div_mpz(result.w, result.w, U1)
                        mpz_divexact(result.d, result.d, U1)

        fmpz_clear(content)

        return result

    def __getitem__(self, int i):
        # quadratic case

        # general number
        print "stub for getitem"
        return self._base(1)

    cpdef RingElement _mul_(self, RingElement _right):
        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> PY_NEW(QuaternionAlgebraElement_number_field)

        mpz_set_si(result.d, 1)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)

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

        cdef fmpz_t content = fmpz_init(fmpz_poly_max_limbs(result.x)) # TODO: think about how big this should be (probably the size of d)

        fmpz_poly_content(content, result.x)
        fmpz_to_mpz(U1, content)
        mpz_gcd(U1, result.d, U1)
        if mpz_cmp_ui(U1, 1) != 0:
            fmpz_poly_content(content, result.y)
            fmpz_to_mpz(U2, content)
            mpz_gcd(U1, U1, U2)
            if mpz_cmp_ui(U1, 1) != 0:
                fmpz_poly_content(content, result.z)
                fmpz_to_mpz(U2, content)
                mpz_gcd(U1, U1, U2)
                if mpz_cmp_ui(U1, 1) != 0:
                    fmpz_poly_content(content, result.w)
                    fmpz_to_mpz(U2, content)
                    mpz_gcd(U1, U1, U2)
                    if mpz_cmp_ui(U1, 1) != 0:
                        fmpz_poly_scalar_div_mpz(result.x, result.x, U1)
                        fmpz_poly_scalar_div_mpz(result.y, result.y, U1)
                        fmpz_poly_scalar_div_mpz(result.z, result.z, U1)
                        fmpz_poly_scalar_div_mpz(result.w, result.w, U1)
                        mpz_divexact(result.d, result.d, U1)

        fmpz_clear(content)

        return result
