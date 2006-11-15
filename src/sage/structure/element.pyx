r"""
Elements

AUTHORS:
   -- David Harvey (2006-10-16): changed CommutativeAlgebraElement to derive
   from CommutativeRingElement instead of AlgebraElement
   -- David Harvey (2006-10-29): implementation and documentation of new
   arithmetic architecture
   -- William Stein (2006-11): arithmetic architecture -- pushing it through to completion.


\subsection{The Abstract Element Class Heierarchy}
This is the abstract class heierchary, i.e., these are all
abstract base classes.
\begin{verbatim}
SageObject
    Element
        ModuleElement
            AdditiveGroupElement

        MonoidElement
            MultiplicativeGroupElement

        RingElement
            CommutativeRingElement
                IntegralDomainElement
                    DedekindDomainElement
                        PrincipalIdealDomainElement
                            EuclideanDomainElement
            FieldElement
                FiniteFieldElement
            AlgebraElement   (note -- can't derive from module, since no multiple inheritence)
                CommutativeAlgebraElement
            InfinityElement
\end{verbatim}

\subsection{How to Define a New Element Class}
SAGE has a special system in place for handling arithmetic operations
for all Element subclasses. There are various rules that must be
followed by both arithmetic implementors and callers.

A quick summary for the impatient:
\begin{itemize}
 \item DO NOT OVERRIDE _add_c. EVER. THANKS.
 \item DO NOT CALL _add_c_impl DIRECTLY.
 \item To implement addition for a python class, override def _add_().
 \item To implement addition for a pyrex class, override cdef _add_c_impl().
 \item If you want to add x and y, whose parents you know are IDENTICAL,
   you may call _add_(x, y) (from python or pyrex) or _add_c(x, y) (from
   pyrex -- this will be faster). This will be the fastest way to guarantee
   that the correct implementation gets called. Of course you can still
   always use "x + y".
\end{itemize}

Now in more detail. The aims of this system are to provide (1) an efficient
calling protocol from both python and pyrex, (2) uniform coercion semantics
across SAGE, (3) ease of use, (4) readability of code.

We will take addition of RingElements as an example; all other operators
and classes are similar. There are four relevant functions.

{\bf def RingElement.__add__}

   This function is called by python or pyrex when the binary "+" operator
   is encountered. It ASSUMES that at least one of its arguments is a
   RingElement; only a really twisted programmer would violate this
   condition. It has a fast pathway to deal with the most common case
   where the arguments have the same parent. Otherwise, it uses the coercion
   module to work out how to make them have the same parent. After any
   necessary coercions have been performed, it calls _add_c to dispatch to
   the correct underlying addition implementation.

   Note that although this function is declared as def, it doesn't have the
   usual overheads associated with python functions (either for the caller
   or for __add__ itself). This is because python has optimised calling
   protocols for such special functions.

{\bf cdef RingElement._add_c}

   DO ***NOT*** OVERRIDE THIS FUNCTION.

   The two arguments to this function MUST have the SAME PARENT.
   Its return value MUST have the SAME PARENT as its arguments.

   If you want to add two objects from pyrex, and you know that their
   parents are the same object, you are encouraged to call this function
   directly, instead of using "x + y".

   This function dispatches to either _add_ or _add_c_impl as appropriate.
   It takes necessary steps to decide whether a pyrex implementation of
   _add_c_impl has been overridden by some python implementation of _add_.
   The code is optimised in favour of reaching _add_c_impl as soon as
   possible.

{\bf def RingElement._add_}

   This is the function you should override to implement addition in a
   python subclass of RingElement.

   WARNING: if you override this in a *pyrex* class, it won't get called.
   You should override _add_c_impl instead. It is especially important to
   keep this in mind whenever you move a class down from python to pyrex.

   The two arguments to this function are guaranteed to have the
   SAME PARENT. Its return value MUST have the SAME PARENT as its
   arguments.

   If you want to add two objects from python, and you know that their
   parents are the same object, you are encouraged to call this function
   directly, instead of using "x + y".

   The default implementation of this function is to call _add_c_impl,
   so if no-one has defined a python implementation, the correct pyrex
   implementation will get called.

{\bf cdef RingElement._add_c_impl}

   This is the function you should override to implement addition in a
   pyrex subclass of RingElement.

   The two arguments to this function are guaranteed to have the
   SAME PARENT. Its return value MUST have the SAME PARENT as its
   arguments.

   DO ***NOT*** CALL THIS FUNCTION DIRECTLY.

   (Exception: you know EXACTLY what you are doing, and you know exactly
   which implementation you are trying to call; i.e. you're not trying to
   write generic code. In particular, if you call this directly, your code
   will not work correctly if you run it on a python class derived from
   a pyrex class where someone has redefined _add_ in python.)

   The default implementation of this function is to raise a
   NotImplementedError, which will happen if no-one has supplied
   implementations of either _add_ or _add_c_impl.

"""


##################################################################
# Generic element, so all this functionality must be defined
# by any element.  Derived class must call __init__
##################################################################

include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"
include "../ext/python.pxi"

import operator

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



def is_Element(x):
    """
    Return True if x is of type Element.

    EXAMPLES:
        sage: is_Element(2/3)
        True
        sage: is_Element(QQ^3)
        False
    """
    return IS_INSTANCE(x, Element)

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
        self._parent = <sage.structure.parent.Parent> parent

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

    def __cmp__(left, right):
        if not have_same_parent(left, right):
            # TODO: can make faster using the cdef interface to coerce
            return cmp_c(left, right)

        if HAS_DICTIONARY(left):
            return left._cmp_(right)
        else:
            return left._cmp_c_impl(right)


    def is_zero(self):
        return bool(self == self._parent(0))


    cdef _richcmp(left, right, int op):
        """
        Compare left and right.
        """
        cdef int r

        if not have_same_parent(left, right):
            try:
                left, right = canonical_coercion_c(left, right)
                r = cmp(left, right)
            except TypeError:
                r = cmp(type(left), type(right))
                if r == 0:
                    r = -1
        else:
            if HAS_DICTIONARY(left):   # fast check
                r = left._cmp_(right)
            else:
                r = left._cmp_c_impl(right)

        if op == 0:  #<
            return bool(r  < 0)
        elif op == 2: #==
            return bool(r == 0)
        elif op == 4: #>
            return bool(r  > 0)
        elif op == 1: #<=
            return bool(r <= 0)
        elif op == 3: #!=
            return bool(r != 0)
        elif op == 5: #>=
            return bool(r >= 0)

    ####################################################################
    # For a derived Python class, you **must** put the following in
    # your subclasses, in order for it to take advantage of the
    # above generic comparison code.  You must also define
    # _cmp_c_impl.
    ####################################################################
    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        raise NotImplementedError, "sort algorithm for elements of type %s not implemented"%(type(left))


    ####################################################################
    # For a derived Python class, put the following:
    ####################################################################
    # import sage.structure.element as element
    #def __cmp__(left, right):
    #    if not isinstance(right, Element) or not (left.parent() is right.parent()):
    #        return element.cmp(self, other)
    #    else:
    #        ...
    # TODO: can we do something more systematic (and faster)?!



def is_ModuleElement(x):
    """
    Return True if x is of type ModuleElement.

    This is even faster than using isinstance inline.

    EXAMPLES:
        sage: is_ModuleElement(2/3)
        False
        sage: is_ModuleElement((QQ^3).0)
        True
    """
    return IS_INSTANCE(x, ModuleElement)

cdef class ModuleElement(Element):
    """
    Generic element of a module.
    """
    ##################################################
    def is_zero(self):
        return bool(self == self._parent(0))

    ##################################################
    # Addition
    ##################################################
    def __add__(self, right):
        """
        Top-level addition operator for ModuleElements.

        See extensive documentation at the top of element.pyx.
        """

        # Try fast pathway if they are both ModuleElements and the parents
        # match.

        # (We know at least one of the arguments is a ModuleElement. So if
        # their types are *equal* (fast to check) then they are both
        # ModuleElements. Otherwise use the slower test via PY_TYPE_CHECK.)
        if have_same_parent(self, right):
            return (<ModuleElement>self)._add_c(<ModuleElement>right)
        else:
            return bin_op_c(self, right, operator.add)
        return self._add_(right)


    cdef ModuleElement _add_c(self, ModuleElement right):
        """
        Addition dispatcher for ModuleElements.

        DO NOT OVERRIDE THIS FUNCTION.

        See extensive documentation at the top of element.pyx.
        """
        if HAS_DICTIONARY(self):   # fast check
            # TODO: this bit will be unnecessarily slow if someone derives
            # from the pyrex class *without* overriding _add_, since then
            # we'll be making an unnecessary python call to _add_, which will
            # end up in _add_c_impl anyway. There must be a simple way to
            # distinguish this situation. It's complicated because someone
            # can even override it at the instance level (without overriding
            # it in the class.)
            return self._add_(right)
        else:
            # Must be a pure Pyrex class.
            return self._add_c_impl(right)


    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Pyrex classes should override this function to implement addition.

        DO NOT CALL THIS FUNCTION DIRECTLY.

        See extensive documentation at the top of element.pyx.
        """
        raise NotImplementedError


    def _add_(ModuleElement self, ModuleElement right):
        """
        Python classes should override this function to implement addition.

        See extensive documentation at the top of element.pyx.
        """
        return self._add_c_impl(right)

    ##################################################
    # Subtraction
    ##################################################

    def __sub__(self, right):
        """
        Top-level subtraction operator for ModuleElements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent(self, right):
            return (<ModuleElement>self)._sub_c(<ModuleElement>right)
        else:
            return bin_op_c(self, right, operator.sub)
        return self._sub_(right)

    cdef ModuleElement _sub_c(self, ModuleElement right):
        """
        Subtraction dispatcher for ModuleElements.

        DO NOT OVERRIDE THIS FUNCTION.

        See extensive documentation at the top of element.pyx.
        """

        if HAS_DICTIONARY(self):   # fast check
            return self._sub_(right)
        else:
            # Must be a pure Pyrex class.
            return self._sub_c_impl(right)


    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Pyrex classes should override this function to implement subtraction.

        DO NOT CALL THIS FUNCTION DIRECTLY.

        See extensive documentation at the top of element.pyx.
        """
        # default implementation is to use the negation and addition
        # dispatchers:
        return self._add_c(right._neg_c())


    def _sub_(ModuleElement self, ModuleElement right):
        """
        Python classes should override this function to implement subtraction.

        See extensive documentation at the top of element.pyx.
        """
        return self._sub_c_impl(right)

    ##################################################
    # Negation
    ##################################################

    def __neg__(self):
        """
        Top-level negation operator for ModuleElements.
        See extensive documentation at the top of element.pyx.
        """
        # We ASSUME that self is a ModuleElement. No type checks.
        return (<ModuleElement>self)._neg_c()


    cdef ModuleElement _neg_c(self):
        """
        Negation dispatcher for ModuleElements.
        DO NOT OVERRIDE THIS FUNCTION.
        See extensive documentation at the top of element.pyx.
        """

        if HAS_DICTIONARY(self):   # fast check
            return self._neg_()
        else:
            # Must be a pure Pyrex class.
            return self._neg_c_impl()


    cdef ModuleElement _neg_c_impl(self):
        """
        Pyrex classes should override this function to implement negation.
        DO NOT CALL THIS FUNCTION DIRECTLY.
        See extensive documentation at the top of element.pyx.
        """
        # default implementation is to try multiplying by -1.
        return bin_op_c(-1, self, operator.mul)


    def _neg_(ModuleElement self):
        """
        Python classes should override this function to implement negation.

        See extensive documentation at the top of element.pyx.
        """
        return self._neg_c_impl()

    ##################################################
    # Scalar multiplication
    ##################################################
    cdef ModuleElement _lmul_c(self, RingElement left):
        """
        DO NOT OVERRIDE THIS FUNCTION.  OK to call.
        """
        if HAS_DICTIONARY(self):
            return self._lmul_c_impl(left)
        else:
            return self._lmul_(left)

    cdef ModuleElement _lmul_c_impl(self, RingElement left):
        """
        Default left multiplication, which is to try to canonically
        coerce the scalar to the integers and do that multiplication,
        which is always defined.
        """
        from sage.rings.all import ZZ
        n = ZZ._coerce_(left)
        a = self
        if n < 0:
            a = -a
            n = -n
        prod = self._parent(0)
        aprod = a
        while True:
            if n&1 > 0: prod = prod + aprod
            n = n >> 1
            if n != 0:
                aprod = aprod + aprod
            else:
                break
        return prod

    def _lmul_(self, left):
        return self._lmul_c_impl(left)

    cdef ModuleElement _rmul_c(self, RingElement right):
        """
        DO NOT OVERRIDE THIS FUNCTION.
        """
        if HAS_DICTIONARY(self):
            return self._rmul_c_impl(right)
        else:
            return self._rmul_(right)

    cdef ModuleElement _rmul_c_impl(self, RingElement right):
        """
        Default right multiplication, which is to try to canonically
        coerce the scalar to the integers and do that multiplication,
        which is always defined.
        """
        return self._lmul_c(right)

    def _rmul_(self, right):
        return self._rmul_c_impl(right)


    ##################################################
    # Other properties
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



def is_MonoidElement(x):
    """
    Return True if x is of type MonoidElement.
    """
    return IS_INSTANCE(x, MonoidElement)

cdef class MonoidElement(Element):
    """
    Generic element of a monoid.
    """

    #############################################################
    # Multiplication
    #############################################################
    def __mul__(self, right):
        """
        Top-level multiplication operator for ring elements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent(self, right):
            return (<MonoidElement>self)._mul_c(<MonoidElement>right)
        return bin_op_c(self, right, operator.mul)

    cdef MonoidElement _mul_c(self, MonoidElement right):
        """
        Multiplication dispatcher for RingElements.
        DO NOT OVERRIDE THIS FUNCTION.
        See extensive documentation at the top of element.pyx.
        """
        if HAS_DICTIONARY(self):   # fast check
            return self._mul_(right)
        else:
            return self._mul_c_impl(right)


    cdef MonoidElement _mul_c_impl(self, MonoidElement right):
        """
        Pyrex classes should override this function to implement multiplication.
        DO NOT CALL THIS FUNCTION DIRECTLY.
        See extensive documentation at the top of element.pyx.
        """
        raise NotImplementedError

    def _mul_(self, right):
        return self._mul_c_impl(right)

    #############################################################
    # Other generic functions that should be available to
    # any monoid element.
    #############################################################
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

        power = self._parent(1)
        apow = a
        while True:
            if n&1 > 0: power = power*apow
            n = n >> 1
            if n != 0:
                apow = apow*apow
            else:
                break

        return power


def is_AdditiveGroupElement(x):
    """
    Return True if x is of type AdditiveGroupElement.
    """
    return IS_INSTANCE(x, AdditiveGroupElement)

cdef class AdditiveGroupElement(ModuleElement):
    """
    Generic element of an additive group.
    """
    def order(self):
        """
        Return additive order of element
        """
        return self.additive_order()

    def __invert__(self):
        raise NotImplementedError, "multiplicative inverse not defined for additive group elements"

    cdef ModuleElement _lmul_c_impl(self, RingElement left):
        return self._rmul_c_impl(left)

    cdef ModuleElement _rmul_c_impl(self, RingElement right):
        cdef int m
        m = int(right)  # a little worrisome.
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


def is_MultiplicativeGroupElement(x):
    """
    Return True if x is of type MultiplicativeGroupElement.
    """
    return IS_INSTANCE(x, MultiplicativeGroupElement)

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
            return bin_op_c(self, right, operator.div)
        return self.__div__(right) # in sage all divs are true

    def __div__(self, right):
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or not (right.parent() != self.parent()):
            return bin_op_c(self, right, operator.div)
        return self._div_(right)

    def _div_(self, right):
        raise NotImplementedError

    def __invert__(self):
        return 1/self


def is_RingElement(x):
    """
    Return True if x is of type RingElement.
    """
    return IS_INSTANCE(x, RingElement)

cdef class RingElement(ModuleElement):
    ##################################################
    def is_zero(self):
        return bool(self == self.parent()(0))

    def is_one(self):
        return bool(self == self.parent()(1))

    ##################################
    # Multiplication
    ##################################

    def __mul__(self, right):
        """
        Top-level multiplication operator for ring elements.
        See extensive documentation at the top of element.pyx.
        """
        # Try fast pathway if they are both RingElements and the parents match.
        # (We know at least one of the arguments is a RingElement. So if their
        # types are *equal* (fast to check) then they are both RingElements.
        # Otherwise use the slower test via PY_TYPE_CHECK.)
        if have_same_parent(self, right):
            return (<RingElement>self)._mul_c(<RingElement>right)
        return bin_op_c(self, right, operator.mul)

    cdef RingElement _mul_c(self, RingElement right):
        """
        Multiplication dispatcher for RingElements.
        DO NOT OVERRIDE THIS FUNCTION.
        See extensive documentation at the top of element.pyx.
        """
        if HAS_DICTIONARY(self):   # fast check
            return self._mul_(right)
        else:
            return self._mul_c_impl(right)

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Pyrex classes should override this function to implement multiplication.
        DO NOT CALL THIS FUNCTION DIRECTLY.
        See extensive documentation at the top of element.pyx.
        """
        raise NotImplementedError

    def _mul_(RingElement self, RingElement right):
        """
        Python classes should override this function to implement multiplication.
        See extensive documentation at the top of element.pyx.
        """
        return self._mul_c_impl(right)


    ##################################
    # Division
    ##################################

    def __truediv__(self, right):
        if not isinstance(self, Element):
            return bin_op_c(self, right, operator.div)
        return self.__div__(right) # in sage all divs are true

    def __div__(self, right):
        """
        Top-level multiplication operator for ring elements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent(self, right):
            return (<RingElement>self)._div_c(<RingElement>right)
        else:
            return bin_op_c(self, right, operator.div)

    cdef RingElement _div_c(self, RingElement right):
        """
        Multiplication dispatcher for RingElements.
        DO NOT OVERRIDE THIS FUNCTION.
        See extensive documentation at the top of element.pyx.
        """
        if HAS_DICTIONARY(self):   # fast check
            return self._div_(right)
        else:
            return self._div_c_impl(right)

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Pyrex classes should override this function to implement division.
        DO NOT CALL THIS FUNCTION DIRECTLY.
        See extensive documentation at the top of element.pyx.
        """
        return self._parent.fraction_field()(self, right)

    def _div_(RingElement self, RingElement right):
        """
        Python classes should override this function to implement division.
        """
        return self._div_c_impl(right)

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




def is_CommutativeRingElement(x):
    """
    Return True if x is of type CommutativeRingElement.
    """
    return IS_INSTANCE(x, CommutativeRingElement)

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
            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^3 + x + 1
            sage: f.mod(x + 1)
            -1

        When little is implemented about a given ring, then mod may
        return simply return $f$.  For example, reduction is not
        implemented for $\Z[x]$ yet. (TODO!)

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = x^3 + x + 1
            sage: f.mod(x + 1)
            x^3 + x + 1



        EXAMPLE: Multivariate polynomials
        We reduce a polynomial in two variables modulo a polynomial
        and an ideal:
            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
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

def is_IntegralDomainElement(x):
    """
    Return True if x is of type IntegralDomainElement.
    """
    return IS_INSTANCE(x, IntegralDomainElement)

cdef class IntegralDomainElement(CommutativeRingElement):
    pass


def is_DedekindDomainElement(x):
    """
    Return True if x is of type DedekindDomainElement.
    """
    return IS_INSTANCE(x, DedekindDomainElement)

cdef class DedekindDomainElement(IntegralDomainElement):
    pass

def is_PrincipalIdealDomainElement(x):
    """
    Return True if x is of type PrincipalIdealDomainElement.
    """
    return IS_INSTANCE(x, PrincipalIdealDomainElement)

cdef class PrincipalIdealDomainElement(DedekindDomainElement):
    def lcm(self, right):
        """
        Returns the least common multiple of self and right.
        """
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or not (right.parent() is self.parent()):
            return bin_op_c(self, right, lcm)
        return self._lcm(right)

    def gcd(self, right):
        """
        Returns the gcd of self and right, or 0 if both are 0.
        """
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or not (right.parent() is self.parent()):
            return bin_op_c(self, right, gcd)
        return self._gcd(right)

    def xgcd(self, right):
        r"""
        Return the extended gcd of self and other, i.e., elements $r, s, t$ such that
        $$
           r = s \cdot self + t \cdot other.
        $$
        """
        if not isinstance(self, Element) or not isinstance(right, Element) \
               or not (right.parent() is self.parent()):
            return bin_op_c(self, right, xgcd)
        return self._xgcd(right)


# This is pretty nasty low level stuff. The idea is to speed up construction
# of EuclideanDomainElements (in particular Integers) by skipping some tp_new
# calls up the inheritance tree.
PY_SET_TP_NEW(EuclideanDomainElement, Element)

def is_EuclideanDomainElement(x):
    """
    Return True if x is of type EuclideanDomainElement.
    """
    return IS_INSTANCE(x, EuclideanDomainElement)

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

def is_FieldElement(x):
    """
    Return True if x is of type FieldElement.
    """
    return IS_INSTANCE(x, FieldElement)

cdef class FieldElement(CommutativeRingElement):

    def is_unit(self):
        """
        Return True if self is a unit in its parent ring.

        EXAMPLES:
            sage: a = 2/3; a.is_unit()
            True

        On the other hand, 2 is not a unit, since its parent is ZZ.
            sage: a = 2; a.is_unit()
            False
            sage: parent(a)
            Integer Ring

        However, a is a unit when viewed as an element of QQ:
            sage: a = QQ(2); a.is_unit()
            True
        """
        return bool(not self.is_zero())

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
        if not isinstance(right, FieldElement) or not (right.parent() is self.parent()):
            right = self.parent()(right)
        return self/right, 0

def is_FiniteFieldElement(x):
    """
    Return True if x is of type FiniteFieldElement.
    """
    return IS_INSTANCE(x, FiniteFieldElement)

cdef class FiniteFieldElement(FieldElement):
    pass

def is_AlgebraElement(x):
    """
    Return True if x is of type AlgebraElement.
    """
    return IS_INSTANCE(x, AlgebraElement)

cdef class AlgebraElement(RingElement):
    pass

def is_CommutativeAlgebraElement(x):
    """
    Return True if x is of type CommutativeAlgebraElement.
    """
    return IS_INSTANCE(x, CommutativeAlgebraElement)

cdef class CommutativeAlgebraElement(CommutativeRingElement):
    pass

def is_InfinityElement(x):
    """
    Return True if x is of type InfinityElement.
    """
    return IS_INSTANCE(x, InfinityElement)

cdef class InfinityElement(RingElement):
    pass


cdef int have_same_parent(left, right):
    """
    Return nonzero true value if and only if left and right are
    elements and have the same parent.
    """
    # (We know at least one of the arguments is an Element. So if
    # their types are *equal* (fast to check) then they are both
    # Elements.  Otherwise use the slower test via PY_TYPE_CHECK.)
    if PY_TYPE(left) is PY_TYPE(right):
        return (<Element>left)._parent is (<Element>right)._parent

    if PY_TYPE_CHECK(right, Element) and PY_TYPE_CHECK(left, Element):
        return (<Element>left)._parent is (<Element>right)._parent

    return 0




#################################################################################
#
#  Coercion of elements
#
#################################################################################
import __builtin__
import operator

cimport sage.modules.module
import  sage.modules.module

#################################################################################
# parent
#################################################################################
cdef parent_c(x):
    if PY_TYPE_CHECK(x,Element):
        return (<Element>x)._parent
    return <object>PY_TYPE(x)

def parent(x):
    return parent_c(x)

#################################################################################
# coerce
#################################################################################
def coerce(p, x):
    try:
        return p._coerce_(x)
    except AttributeError:
        return p(x)

#################################################################################
# canonical coercion of two ring elements into one of their parents.
#################################################################################
def canonical_coercion(x, y):
    return canonical_coercion_c(x,y)

cdef canonical_coercion_c(x, y):
    cdef int i
    xp = parent_c(x)
    yp = parent_c(y)
    if xp is yp:
        return x, y
    if PY_IS_NUMERIC(x):
        return yp(x), y
    elif PY_IS_NUMERIC(y):
        return x, xp(y)
    try:
        if xp.has_coerce_map_from(yp):
            return x, xp(y)
    except AttributeError:
        pass
    try:
        if yp.has_coerce_map_from(xp):
            return yp(x), y
    except AttributeError:
        pass
    raise TypeError, "unable to find a common canonical parent for x and y"

def canonical_base_coercion(x, y):
    try:
        xb = x.base_ring()
    except AttributeError:
        #raise TypeError, "unable to find base ring for %s (parent: %s)"%(x,x.parent())
        raise TypeError, "unable to find base ring"
    try:
        yb = y.base_ring()
    except AttributeError:
        raise TypeError, "unable to find base ring"
        #raise TypeError, "unable to find base ring for %s (parent: %s)"%(y,y.parent())
    try:
        b = canonical_coercion_c(xb(0),yb(0))[0].parent()
    except TypeError:
        raise TypeError, "unable to find base ring"
        #raise TypeError, "unable to find a common base ring for %s (base ring: %s) and %s (base ring %s)"%(x,xb,y,yb)
    return x.change_ring(b), y.change_ring(b)


def bin_op(x, y, op):
    return bin_op_c(x,y,op)

cdef bin_op_c(x, y, op):
    """
    Compute x op y, where coercion of x and y works according to
    SAGE's coercion rules.
    """

    # 1. Try canonical ring element coercion.
    try:
        x, y = canonical_coercion_c(x, y)
        return op(x,y)
    except TypeError, msg:
        pass

    # If the op is multiplication, then some other algebra multiplications
    # may be defined
    if not op is operator.mul:
        raise TypeError, "%s\nno canonical coercion."%msg

    # 2. Try scalar multiplication.
    # No way to multiply x and y using the ``coerce into a canonical
    # parent'' rule.
    # The next rule to try is scalar multiplication by coercing
    # into the base ring.
    cdef int x_is_modelt, y_is_modelt
    y_is_modelt = PY_TYPE_CHECK(y, ModuleElement)
    if y_is_modelt:
        # First try to coerce x into the base ring of y if y is an element.
        try:
            x = (<ModuleElement> y)._parent._base(x)
            return (<ModuleElement> y)._lmul_c(x)     # the product x * y
        except TypeError, msg:
            pass

    x_is_modelt = PY_TYPE_CHECK(x, ModuleElement)
    if x_is_modelt:
        # That did not work.  Try to coerce y into the base ring of x.
        try:
            y = (<ModuleElement> x)._parent._base(y)
            return (<ModuleElement> x)._rmul_c(y)    # the product x * y
        except TypeError:
            pass

    if y_is_modelt and x_is_modelt:
        # 3. Both coercion failed, but both are module elements.
        # Try base extending the right object by the parent of the left

        ## TODO -- WORRY -- only unambiguous if one succeeds!
        if  PY_TYPE_CHECK(x, RingElement):
            try:
                return x * y.base_extend((<RingElement>x)._parent)
            except (TypeError, AttributeError), msg:
                pass
        # Also try to base extending the left object by the parent of the right
        if  PY_TYPE_CHECK(y, RingElement):
            try:
                return y * x.base_extend((<Element>y)._parent)
            except (TypeError, AttributeError), msg:
                pass

    # 4. Try _l_action or _r_action.
    # Test to see if an _r_action or _l_action is
    # defined on either side.
    try:
        return y._r_action(x)
    except (AttributeError, TypeError):
        pass
    try:
        return x._l_action(y)
    except (AttributeError, TypeError):
        pass

    raise TypeError, "unable to multiply %s times %s"%(x,y)


def coerce_cmp(x,y):  # external interface to cmp_cdef
    try:
        x, y = canonical_coercion_c(x, y)
        return cmp(left, right)
    except TypeError:
        return cmp(type(x), type(y))



###############################################################################

def lcm(x,y):
    from sage.rings.arith import lcm
    return lcm(x,y)

def gcd(x,y):
    from sage.rings.arith import gcd
    return gcd(x,y)

def xgcd(x,y):
    from sage.rings.arith import xgcd
    return xgcd(x,y)

