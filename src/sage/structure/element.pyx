"""
Elements

PYREX: sage.structure.element

AUTHORS:
   -- David Harvey (2006-10-16): changed CommutativeAlgebraElement to derive
   from CommutativeRingElement instead of AlgebraElement
"""

##################################################################
# Generic element, so all this functionality must be defined
# by any element.  Derived class must call __init__
##################################################################

include "../ext/cdefs.pxi"
include '../ext/stdsage.pxi'

import operator

import coerce

# This classes uses element.pxd.  To add data members, you
# must change that file.

def make_element(_class, _dict, parent):
    """
    Used for unpickling Element objects (and subclasses).

    This should work for any Python class deriving from Element, as long
    as it doesn't implement some screwy __new__() method.

    See also Element.__reduce__().
    """
    new_object = _class.__new__(_class)
    new_object._set_parent(parent)
    new_object.__dict__ = _dict
    return new_object


cdef class Element(sage_object.SageObject):
    """
    Generic element of a structure. All other types of elements
    (RingElement, ModuleElement, etc) derive from this type.

    Subtypes must either call __init__() to set _parent, or may
    set _parent themselves if that would be more efficient.
    """

    def __init__(self, parent):
        r"""
        INPUT:
            parent -- a SageObject
        """
        self._parent = parent

    def _set_parent(self, parent):
        r"""
        INPUT:
            parent -- a SageObject
        """
        self._parent = parent


    def _repr_(self):
        return "Generic element of a structure"

    def __reduce__(self):
        return (make_element, (self.__class__, self.__dict__, self._parent))

    def __hash__(self):
        return hash(str(self))

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of self in codomain under the map that sends
        the images of the generators of the parent of self to the
        tuple of elements of im_gens.
        """
        raise NotImplementedError


    def base_ring(self):
        """
        Returns the base ring of this element's parent (if that makes sense).
        """
        return self._parent.base_ring()

    def category(self):
        from sage.categories.category import Elements
        return Elements(self._parent)

    def parent(self, x=None):
        """
        Returns parent of this element; or, if the optional argument x is
        supplied, the result of coercing x into the parent of this element.
        """
        if x is None:
            return self._parent
        else:
            return self._parent(x)

    def __xor__(self, right):
        raise RuntimeError, "Use ** for exponentiation, not '^', which means xor\n"+\
              "in Python, and has the wrong precedence."

    def _coeff_repr(self, no_space=True):
        if self._is_atomic():
            s = str(self)
        else:
            s = "(%s)"%self
        if no_space:
            return s.replace(' ','')
        return s

    def _latex_coeff_repr(self):
        try:
            s = self._latex_()
        except AttributeError:
            s = str(self)
        if self._is_atomic():
            return s
        else:
            return "\\left(%s\\right)"%s

    def _is_atomic(self):
        if self._parent.is_atomic_repr():
            return True
        s = str(self)
        return bool(s.find("+") == -1 and s.find("-") == -1 and s.find(" ") == -1)

    def _rich_to_bool(self, int op, int n):
        if op == 0:
            return bool(n < 0)
        elif op == 1:
            return bool(n <= 0)
        elif op == 2:
            return bool(n == 0)
        elif op == 3:
            return bool(n != 0)
        elif op == 4:
            return bool(n > 0)
        elif op == 5:
            return bool(n >= 0)
        assert False, "bug in _rich_to_bool"

    #def __richcmp__(self, right, int op):
        # Warning -- this is *not* inherited by derived Pyrex types,
        # so you must copy it to them.  It works fine for Python
        # derived classes (even of derived Pyrex types).  Weird.
        # For some important details about this weirdness, see
        # the email by Lenard Lindstrom in the notes subdirectory.
    #    print "hi", self, right, op
    #    cdef int n
    #    if not isinstance(right, Element) or right.parent() != self.parent():
    #        try:
    #            n = coerce.cmp(self, right)
    #        except TypeError:
    #            n = -1
    #    else:
    #        n = self.__cmp__(right)
    #    return self._rich_to_bool(op, n)

    def __cmp__(self, other):
        raise NotImplementedError
    #    return self._cmp(other)

    #def act(self, right):
    #    try:
    #        return right.r_act(self)
    #    except TypeError:
    #        raise TypeError, "right action of %s on %s not defined"%(right, self)


    def is_zero(self):
        return bool(self == self._parent(0))

cdef class ModuleElement(Element):
    """
    Generic element of a module.
    """
    ##################################################
    def is_zero(self):
        return bool(self == self._parent(0))

##     def is_nonzero(self):
##         return not self.is_zero()

    def __add__(self, right):
        if not isinstance(self, Element) or \
               not isinstance(right, Element) or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.add)
        return self._add_(right)

    def _add_(self, right):
        raise NotImplementedError

    def __sub__(self, right):
        if not isinstance(self, Element) or \
               not isinstance(right, Element) or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.sub)
        return self._sub_(right)

    def _sub_(self, right):
        return self + (-right)

    def __neg__(self):
        return coerce.bin_op(-1, self, operator.mul)

    def __pos__(self):
        return self

    # addition is commutative in a module:
    def __rsub__(self, left):
        return self.parent()(left) - self

    def __radd__(self, left):
        return self.parent()(left) + self

    def __rmul__(self, n):
        if not isinstance(n, (int, long)):
            raise TypeError, "multiplication on by %s not defined"%n
        return self*n

    ##################################################
    def order(self):
        """
        Return the additive order of self.
        """
        return self.additive_order()

    def additive_order(self):
        """
        Return the additive order of self.
        """
        raise NotImplementedError


def is_ModuleElement(x):
    return isinstance(x, ModuleElement)


cdef class MonoidElement(Element):
    """
    Generic element of a monoid.
    """
    def order(self):
        """
        Return the multiplicative order of self.
        """
        return self.multiplicative_order()

    def multiplicative_order(self):
        """
        Return the multiplicative order of self.
        """
        raise NotImplementedError

    def __mul__(self, right):
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.mul)
        return self._mul_(right)

    def _mul_(self, right):
        return NotImplementedError

    def __pow__(self, n, dummy):
        cdef int i

        if isinstance(n, float):
            raise TypeError, "raising %s to the power of the float %s not defined"%(self, n)

        n = int(n)

        a = self
        power = self.parent()(1)
        if n < 0:
            n = -n
            a = ~self
        elif n == 0:
            return power

        power = self.parent()(1)
        apow = a
        while True:
            if n&1 > 0: power = power*apow
            n = n >> 1
            if n != 0:
                apow = apow*apow
            else:
                break

        return power


cdef class AdditiveGroupElement(ModuleElement):
    """
    Generic element of an additive group.
    """
    def order(self):
        """
        Return additive order of element
        """
        return self.additive_order()

    #def _mul_(self, x):
    #    raise ArithmeticError, "multiplication not defined in an additive group"

    def __invert__(self):
        raise NotImplementedError, "multiplicative inverse not defined for additive group elements"

    def _mul_(self, m):
        m = int(m)
        if m<0:
            return (-self)*(-m)
        if m==1:
            return self
        P = self.scheme()(0)
        if m==0:
            return P
        power = P
        i = 0
        apow2 = self
        while ((m>>i) > 0):
            if((m>>i) & 1):
                power = power + apow2
            apow2 = apow2 + apow2
            i = i + 1
        return power


cdef class MultiplicativeGroupElement(MonoidElement):
    """
    Generic element of a multiplicative group.
    """
    def order(self):
        """
        Return the multiplicative order of self.
        """
        return self.multiplicative_order()

    def _add_(self, x):
        raise ArithmeticError, "addition not defined in a multiplicative group"

    def __truediv__(self, right):
        if not isinstance(self, Element):
            return coerce.bin_op(self, right, operator.div)
        return self.__div__(right) # in sage all divs are true

    def __div__(self, right):
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.div)
        return self._div_(right)

    def _div_(self, right):
        raise NotImplementedError

    def __invert__(self):
        return 1/self


cdef class RingElement(Element):
    ##################################################
    def is_zero(self):
        return bool(self == self.parent()(0))

    def is_one(self):
        return bool(self == self.parent()(1))

##     def is_nonzero(self):
##         return not self.is_zero()

    def __add__(self, right):
        """
        Addition for RingElements.

        This method tries to call _add_sibling_cdef or _add_sibling. If the
        parents disagree, it tries canonical coercion to make the parents
        agree.
        """
        # We know at least one of the arguments is a RingElement. So if their
        # types are *equal* (fast to check) then they are both RingElements.
        # Otherwise use the slower test via PY_TYPE_CHECK.

        if (PY_TYPE(self) is PY_TYPE(right)) or \
                   (PY_TYPE_CHECK(right, RingElement) and \
                    PY_TYPE_CHECK(self, RingElement)):
        #if (PY_TYPE(self) is PY_TYPE(right)):
            # If parents agree, then we can use _add_sibling_cdef.
            if (<RingElement>right)._parent is (<RingElement>self)._parent:
                return (<RingElement>self)._add_sibling_cdef(
                                                     <RingElement>right)

        # Fast pathway didn't work. Fall back on the old coercion code.
        # todo:
        # (1) optimise this too
        # (2) change "parent != parent" to "parent is not parent" and see
        #     what breaks
        if not PY_TYPE_CHECK(self, RingElement) or \
                  not PY_TYPE_CHECK(right, RingElement) \
                  or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.add)
        return self._add_(right)


    cdef RingElement _add_sibling_cdef(RingElement self, RingElement right):
        """
        Adds two RingElements. They are assumed to belong to the SAME parent.
        The return value MUST belong to the SAME parent.
        """
        # This default implementation tries to dispatch to a Python version
        # (in case cdef version not available)
        if self.__class__._add_sibling is RingElement._add_sibling:
            # This is to prevent an infinite loop between _add_sibling_cdef
            # and _add_sibling:
            # raise NotImplementedError

            # Actually, we can't do that right now, because most objects
            # don't have _add_sibling yet. So we need to look up _add_
            # as well:
            if self.__class__._add_ is RingElement._add_:
                raise NotImplementedError
            else:
                return self._add_(right)

        else:
            return self._add_sibling(right)


    def _add_sibling(RingElement self, RingElement right):
        """
        Adds two RingElements. They are assumed to belong to the SAME parent.
        The return value MUST belong to the SAME parent.
        """
        # This is the default Python implementation, which dispatches to any
        # available cdef version.
        return self._add_sibling_cdef(right)


    def _add_(RingElement self, RingElement right):
        """
        See docstring for _add_sibling.

        NOTE: The plan is to replace all '_add_' methods in SAGE by
        '_add_sibling'. Eventually this method will be removed.
        """
        self._add_sibling_cdef(right)


    def __sub__(self, right):
        if not isinstance(self, Element) or \
               not isinstance(right, Element) or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.sub)
        return self._sub_(right)

    def _sub_(self, right):
        return self + (-right)

    def __mul__(self, right):
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.mul)
        return self._mul_(right)

    def _mul_(self, right):
        raise NotImplementedError

    def __truediv__(self, right):
        if not isinstance(self, Element):
            return coerce.bin_op(self, right, operator.div)
        return self.__div__(right) # in sage all divs are true

    def __div__(self, right):
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or right.parent() != self.parent():
            return coerce.bin_op(self, right, operator.div)
        return self._div_(right)

    def _div_(self, right):
        return self.parent().fraction_field()(self, right)

    def __neg__(self):
        return coerce.bin_op(-1, self, operator.mul)

    def __pos__(self):
        return self

    def __invert__(self):
        return 1/self

    ##################################################

    def order(self):
        """
        Return the additive order of self.
        """
        return self.additive_order()

    def additive_order(self):
        """
        Return the additive order of self.
        """
        raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of self, if self is a unit, or raise
        \code{ArithmeticError} otherwise.
        """
        if not self.is_unit():
            raise ArithmeticError, "self (=%s) must be a unit to have a multiplicative order."
        raise NotImplementedError

    def is_unit(self):
        if self == 1 or self == -1:
            return True
        raise NotImplementedError

    def __pow__(self, n, dummy):
        cdef int i
        if isinstance(n, float):
            raise TypeError, "raising %s to the power of the float %s not defined"%(self, n)

        n = int(n)
        try:
            return self._pow(n)
        except AttributeError:
            pass

        a = self
        power = self.parent()(1)
        if n < 0:
            n = -n
            a = ~self
        elif n == 0:
            return power
        i = 0
        apow2 = a
        while (n>>i) > 0:
            if (n>>i) & 1:
                power = power * apow2
            if n == 0: break   # to not waste time doing an extra multiplication/increment
            apow2 = apow2 * apow2
            i = i+1
        return power




cdef class CommutativeRingElement(RingElement):
    def _im_gens_(self, codomain, im_gens):
        if len(im_gens) == 1 and self.parent().gen(0) == 1:
            return codomain(self)
        raise NotImplementedError

    def inverse_mod(self, I):
        r"""
        Return an inverse of self modulo the ideal $I$, if defined,
        i.e., if $I$ and self together generate the unit ideal.
        """
        raise NotImplementedError

    def mod(self, I):
        r"""
        Return a representative for self modulo the ideal I (or the ideal
        generated by the elements of I if I is not an ideal.)

        EXAMPLE:  Integers
        Reduction of 5 modulo an ideal:
            sage: n = 5
            sage: n.mod(3*ZZ)
            2

        Reduction of 5 modulo the ideal generated by 3.
            sage: n.mod(3)
            2

        Reduction of 5 modulo the ideal generated by 15 and 6, which is $(3)$.
            sage: n.mod([15,6])
            2


        EXAMPLE: Univiate polynomials
            sage: x = PolynomialRing(QQ).0
            sage: f = x^3 + x + 1
            sage: f.mod(x + 1)
            -1

        When little is implemented about a given ring, then mod may
        return simply return $f$.  For example, reduction is not
        implemented for $\Z[x]$ yet. (TODO!)

            sage: x = PolynomialRing(ZZ).0
            sage: f = x^3 + x + 1
            sage: f.mod(x + 1)
            x^3 + x + 1



        EXAMPLE: Multivariate polynomials
        We reduce a polynomial in two variables modulo a polynomial
        and an ideal:
            sage: x,y,z = PolynomialRing(QQ, 3, 'xyz').gens()
            sage: (x^2 + y^2 + z^2).mod(x+y+z)
            2*z^2 + 2*y*z + 2*y^2

        Notice above that $x$ is eliminated.  In the next example,
        both $y$ and $z$ are eliminated.

            sage: (x^2 + y^2 + z^2).mod( (x - y, y - z) )
            3*z^2
            sage: f = (x^2 + y^2 + z^2)^2; f
            z^4 + 2*y^2*z^2 + y^4 + 2*x^2*z^2 + 2*x^2*y^2 + x^4
            sage: f.mod( (x - y, y - z) )
            9*z^4

        In this example $y$ is eliminated.
            sage: (x^2 + y^2 + z^2).mod( (x^3, y - z) )
            2*z^2 + x^2
        """
        from sage.rings.all import is_Ideal
        if not is_Ideal(I) or not I.ring() is self.parent():
            I = self.parent().ideal(I)
            #raise TypeError, "I = %s must be an ideal in %s"%(I, self.parent())
        return I.reduce(self)

cdef class IntegralDomainElement(CommutativeRingElement):
    pass

cdef class DedekindDomainElement(IntegralDomainElement):

    pass

cdef class PrincipalIdealDomainElement(DedekindDomainElement):
    def lcm(self, right):
        """
        Returns the least common multiple of self and right.
        """
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or right.parent() != self.parent():
            return coerce.bin_op(self, right, coerce.lcm)
        return self._lcm(right)

    def gcd(self, right):
        """
        Returns the gcd of self and right, or 0 if both are 0.
        """
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or right.parent() != self.parent():
            return coerce.bin_op(self, right, coerce.gcd)
        return self._gcd(right)

    def xgcd(self, right):
        r"""
        Return the extended gcd of self and other, i.e., elements $r, s, t$ such that
        $$
           r = s \cdot self + t \cdot other.
        $$
        """
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or right.parent() != self.parent():
            return coerce.bin_op(self, right, coerce.xgcd)
        return self._xgcd(right)


# This is pretty nasty low level stuff. The idea is to speed up construction
# of EuclideanDomainElements (in particular Integers) by skipping some tp_new
# calls up the inheritance tree.
PY_SET_TP_NEW(EuclideanDomainElement, Element)

cdef class EuclideanDomainElement(PrincipalIdealDomainElement):

    def degree(self):
        raise NotImplementedError

    def _gcd(self, other):
        """
        Return the greatest common divisor of self and other.

        Algorithm 3.2.1 in Cohen, GTM 138.
        """
        A = self
        B = other
        while not B.is_zero():
            Q, R = A.quo_rem(B)
            A = B
            B = R
        return A

    def leading_coefficient(self):
        raise NotImplementedError

    def quo_rem(self, other):
        raise NotImplementedError

    def __floordiv__(self,right):
        """
        Quotient of division of self by other.  This is denoted //.
        """
        Q, _ = self.quo_rem(right)
        return Q

    def __mod__(self, other):
        """
        Remainder of division of self by other.
        EXAMPLES:
            sage: x = PolynomialRing(IntegerRing()).gen()
            sage: x % (x+1)
            -1
            sage: (x**3 + x - 1) % (x**2 - 1)
            2*x - 1
        """
        _, R = self.quo_rem(other)
        return R



cdef class FieldElement(CommutativeRingElement):

    def is_unit(self):
        return not self.is_zero()

    def _gcd(self, FieldElement other):
        """
        Return the greatest common divisor of self and other.
        """
        if self.is_zero() and other.is_zero():
            return self
        else:
            return self.parent()(1)

    def _lcm(self, FieldElement other):
        """
        Return the least common multiple of self and other.
        """
        if self.is_zero() and other.is_zero():
            return self
        else:
            return self.parent()(1)

    def _xgcd(self, FieldElement other):
        R = self.parent()
        if not self.is_zero():
            return R(1), ~self, R(0)
        elif not other.is_zero():
            return R(1), R(0), ~self
        else: # both are 0
            return self, self, self


    def quo_rem(self, right):
        if not isinstance(right, FieldElement) or right.parent() != self:
            right = self.parent()(right)
        return self/right, 0


cdef class AlgebraElement(RingElement):
    pass

cdef class CommutativeAlgebraElement(CommutativeRingElement):
    pass

cdef class InfinityElement(RingElement):
    pass


# A Python class
class Element_cmp_:
    """
    Class for defining comparisons between elements.
    """
    def __cmp__(self, right):
        try:
            if right.parent() != self.parent():
                return coerce.cmp(self, right)
        except AttributeError:
            return coerce.cmp(self, right)
        return self._cmp_(right)

    def _cmp_(self, right):
        raise NotImplementedError
