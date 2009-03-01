from sage.structure.element cimport AlgebraElement, RingElement, ModuleElement, Element
from sage.algebras.quaternion_algebra_element cimport QuaternionAlgebraElement_abstract


cdef class QuaternionAlgebraElement_generic(QuaternionAlgebraElement_abstract):
    def __init__(self, parent, _x, _y, _z, _w):
        self._parent = parent
        self.x = _x
        self.y = _y
        self.z = _z
        self.w = _w

    cpdef ModuleElement _add_(self, ModuleElement _right):
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

        x = x1*x2 + y1*y2*a + z1*z2*b - w1*w2*a*b
        y = x1*y2 + y1*x2 - z1*w2*b + w1*z2*b
        z = x1*z2 + y1*w2 + z1*x2 - w1*y2*a
        w = x1*w2 + y1*z2 - z1*y2 + w1*x2

        return QuaternionAlgebraElement_generic(self._parent, x, y, z, w)

    def _repr_(self):
        return "%s + %s*i +%s*j + %s*k"%(self.x,self.y,self.z,self.w)
