from sage.structure.element cimport AlgebraElement, RingElement, ModuleElement, Element
from sage.algebras.quaternion_algebra_element cimport QuaternionAlgebraElement_abstract
from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.arith import lcm

include "../ext/gmp.pxi"
include "../ext/stdsage.pxi"


# variables for holding temporary values computed in
# QuaternionAlgebraElement_rational_field._mul_()
#

cdef mpz_t T1
cdef mpz_t T2
cdef mpz_t t3
cdef mpz_t t4
cdef mpz_t t5
cdef mpz_t t6
cdef mpz_t t7
cdef mpz_t t8
cdef mpz_t s1
cdef mpz_t s2
cdef mpz_t U1
cdef mpz_t U2
#
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



cdef class QuaternionAlgebraElement_generic(QuaternionAlgebraElement_abstract):
    def __init__(self, parent, _x, _y, _z, _w):
        self._parent = parent
        self.x = _x
        self.y = _y
        self.z = _z
        self.w = _w

    def __reduce__(self):
        return QuaternionAlgebraElement_generic, (self._parent, self.x, self.y, self.z, self.w)

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
        return QuaternionAlgebraElement_generic(self._parent, self.x + right.x, self.y + right.y, self.z + right.z, self.w + right.w)

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef QuaternionAlgebraElement_generic right = _right
        return QuaternionAlgebraElement_generic(self._parent, self.x - right.x, self.y - right.y, self.z - right.z, self.w - right.w)

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
        y = (x1 + y1)*(x2 + y2) - t1 - y2 + b*( (z1 + w1)*(z2 - w2) - t3 + t4)
        z = t5 - a*t6 + t7 + a*t8
        w = (x2 - y2)*(z1 + w1) - t5 + t6 + (x1 + y1)*(z2 + w2) - t7 - t8

        return QuaternionAlgebraElement_generic(self._parent, x, y, z, w)

    def _repr_(self):
        i,j,k = self._parent.variable_names()
        v = []
        if self.x:
            v.append(str(self.x))
        if self.y:
            if self.y == 1:
                v.append(i)
            elif self.y == -1:
                v.append("-%s"%i)
            else:
                v.append('%s*%s'%(self.y, i))
        if self.z:
            if self.z == 1:
                v.append(j)
            elif self.z == -1:
                v.append("-%s"%j)
            else:
                v.append('%s*%s'%(self.z, j))
        if self.w:
            if self.w == 1:
                v.append(k)
            elif self.w == -1:
                v.append("-%s"%k)
            else:
                v.append('%s*%s'%(self.w, k))

        return ' + '.join(v).replace('+ -','- ')

cdef class QuaternionAlgebraElement_rational_field(QuaternionAlgebraElement_abstract):
    def __new__(self):
        mpz_init(self.x)
        mpz_init(self.y)
        mpz_init(self.z)
        mpz_init(self.w)
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.d)

    def __init__(self, _x, _y, _z, _w, _a, _b):
        _x = Rational(_x)
        _y = Rational(_y)
        _z = Rational(_z)
        _w = Rational(_w)
        _a = Integer(_a)
        _b = Integer(_b)

        common_denominator = lcm( [A.denominator() for A in (_x, _y, _z, _w)])
        xx, yy, zz, ww = [Integer(A * common_denominator) for A in (_x, _y, _z, _w)]

        mpz_set(self.x, (<Integer>xx).value)
        mpz_set(self.y, (<Integer>yy).value)
        mpz_set(self.z, (<Integer>zz).value)
        mpz_set(self.w, (<Integer>ww).value)
        mpz_set(self.a, (<Integer>_a).value)
        mpz_set(self.b, (<Integer>_b).value)
        mpz_set(self.d, (<Integer>common_denominator).value)

    cpdef ModuleElement _add_(self, ModuleElement _right):
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)

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


        return result

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)

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
        if x:
            v.append(str(x))
        if y:
            if y == 1:
                v.append(i)
            elif y == -1:
                v.append("-%s"%i)
            else:
                v.append('%s*%s'%(y, i))
        if z:
            if z == 1:
                v.append(j)
            elif z == -1:
                v.append("-%s"%j)
            else:
                v.append('%s*%s'%(z, j))
        if w:
            if w == 1:
                v.append(k)
            elif w == -1:
                v.append("-%s"%k)
            else:
                v.append('%s*%s'%(w, k))

        return ' + '.join(v).replace('+ -','- ')



    cpdef RingElement _mul_(self, RingElement _right):
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> PY_NEW(QuaternionAlgebraElement_rational_field)

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)

        mpz_mul(T1, self.x, right.x)    # t1 = x1 * x2
        mpz_mul(T2, self.y, right.y)    # t2 = y1 * y2
        mpz_mul(t3, self.z, right.z)    # t3 = x1 * x2
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
