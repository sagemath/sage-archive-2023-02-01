r"""
Elements

AUTHORS:
   -- David Harvey (2006-10-16): changed CommutativeAlgebraElement to derive
   from CommutativeRingElement instead of AlgebraElement
   -- David Harvey (2006-10-29): implementation and documentation of new
   arithmetic architecture
   -- William Stein (2006-11): arithmetic architecture -- pushing it through to completion.

\subsection{Some technical implementation notes about arithmetic architecture}
(This documentation should get moved into a more useful place.)

SAGE has a special system in place for handling arithmetic operations
for all Element subclasses. There are various rules that must be followed
by both arithmetic implementors and callers.

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

\subsection{The Element Class Heierarchy}
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

            # TODO: can make faster using the cdef interface to coerce
            r = cmp_c(left, right)

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



cdef class ModuleElement(Element):
    """
    Generic element of a module.
    """
    ##################################################
    def is_zero(self):
        return bool(self == self._parent(0))

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


    def __sub__(self, right):
        """
        Top-level subtraction operator for ModuleElements.

        See extensive documentation at the top of element.pyx.
        """

        # Try fast pathway if they are both ModuleElements and the parents
        # match.

        # (We know at least one of the arguments is a ModuleElement. So if
        # their types are *equal* (fast to check) then they are both
        # ModuleElements. Otherwise use the slower test via PY_TYPE_CHECK.)

        if (PY_TYPE(self) is PY_TYPE(right)) or \
                   (PY_TYPE_CHECK(right, ModuleElement) and \
                    PY_TYPE_CHECK(self, ModuleElement)):

            if (<ModuleElement>right)._parent is (<ModuleElement>self)._parent:

                return (<ModuleElement>self)._sub_c(<ModuleElement>right)

        # Fast pathway didn't work.

        # todo:
        # For now we are falling back on the old coercion code.
        # This needs to be optimised and re-thought-out.
        # In particular, the coercion code doesn't yet know about _add_c
        # and all that.
        if not PY_TYPE_CHECK(self, ModuleElement) or \
                  not PY_TYPE_CHECK(right, ModuleElement) \
                  or not (right.parent() is self.parent()):
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

    # addition is commutative in a module:
    def __rsub__(self, left):
        return self.parent()(left) - self

    def __radd__(self, left):
        return self.parent()(left) + self

    def __rmul__(self, n):
        if not isinstance(n, (int, long)):
            raise TypeError, "multiplication on by %s not defined"%n
        return self*n

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

    cdef ModuleElement _mul_left_scalar_c(self, RingElement m):
        """
        DO NOT OVERRIDE THIS FUNCTION.
        """
        if HAS_DICTIONARY(self):
            return self._mul_left_scalar_c_impl(right)
        else:
            return self._mul_left_scalar(right)

    cdef ModuleElement _mul_left_scalar_c_impl(self, RingElement m):
        raise NotImplementedError
    def _mul_left_scalar(self, m):
        return self._mul_left_scalar_c_impl(m)


    cdef ModuleElement _mul_right_scalar_c(self, RingElement m):
        """
        DO NOT OVERRIDE THIS FUNCTION.
        """
        if HAS_DICTIONARY(self):
            return self._mul_right_scalar_c_impl(right)
        else:
            return self._mul_right_scalar(right)

    cdef ModuleElement _mul_right_scalar_c_impl(self, RingElement m):
        raise NotImplementedError
    def _mul_right_scalar(self, m):
        return self._mul_right_scalar_c_impl(m)


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
               or not (right.parent() is self.parent()):
            return bin_op_c(self, right, operator.mul)
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

    def __invert__(self):
        raise NotImplementedError, "multiplicative inverse not defined for additive group elements"

    cdef ModuleElement _mul_left_scalar_c_impl(self, RingElement m):
        return self._mul_right_scalar_c_impl(m)

    cdef ModuleElement _mul_right_scalar_c_impl(self, RingElement m):
        m = int(m)  # a little worrisome.
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


cdef class RingElement(Element):
    ##################################################
    def is_zero(self):
        return bool(self == self.parent()(0))

    def is_one(self):
        return bool(self == self.parent()(1))

    ##################################
    # Addition
    ##################################
    def __add__(self, right):
        """
        Top-level addition operator for RingElements.
        See extensive documentation at the top of element.pyx.
        """
        # Try fast pathway if they are both RingElements and the parents match.
        # (We know at least one of the arguments is a RingElement. So if their
        # types are *equal* (fast to check) then they are both RingElements.
        # Otherwise use the slower test via PY_TYPE_CHECK.)
        if have_same_parent(self, right):
            return (<RingElement>self)._add_c(<RingElement>right)
        else:
            return bin_op_c(self, right, operator.add)

    cdef RingElement _add_c(self, RingElement right):
        """
        Addition dispatcher for RingElements.
        DO NOT OVERRIDE THIS FUNCTION.
        See extensive documentation at the top of element.pyx.
        """

        if HAS_DICTIONARY(self):   # fast check
            # TODO: this bit will be unnecessarily slow if someone derives
            # from the pyrex class *without* overriding _add_, since then
            # we'll be making an unnecessary call to _add_, which will
            # end up in _add_c_impl anyway. There must be a simple way to
            # distinguish this situation. It's complicated because someone
            # can even override it at the instance level (without overriding
            # it in the class.)
            # TODO: if you fix this up here, fix it up for the other operators
            # (sub, neg, etc.) as well. Thanks.
            return self._add_(right)
        else:
            # Must be a pure Pyrex class.
            return self._add_c_impl(right)


    cdef RingElement _add_c_impl(self, RingElement right):
        """
        Pyrex classes should override this function to implement addition.
        DO NOT CALL THIS FUNCTION DIRECTLY.
        See extensive documentation at the top of element.pyx.
        """
        raise NotImplementedError


    def _add_(RingElement self, RingElement right):
        """
        Python classes should override this function to implement addition.
        See extensive documentation at the top of element.pyx.
        """
        return self._add_c_impl(right)


    ##################################
    # Subtraction
    ##################################

    def __sub__(self, right):
        """
        Top-level subtraction operator for RingElements.

        See extensive documentation at the top of element.pyx.
        """

        # Try fast pathway if they are both RingElements and the parents match.

        # (We know at least one of the arguments is a RingElement. So if their
        # types are *equal* (fast to check) then they are both RingElements.
        # Otherwise use the slower test via PY_TYPE_CHECK.)

        if (PY_TYPE(self) is PY_TYPE(right)) or \
                   (PY_TYPE_CHECK(right, RingElement) and \
                    PY_TYPE_CHECK(self, RingElement)):

            if (<RingElement>right)._parent is (<RingElement>self)._parent:

                return (<RingElement>self)._sub_c(<RingElement>right)

        # Fast pathway didn't work.

        # todo:
        # For now we are falling back on the old coercion code.
        # This needs to be optimised and re-thought-out.
        # In particular, the coercion code doesn't yet know about _add_c
        # and all that.
        if not PY_TYPE_CHECK(self, RingElement) or \
                  not PY_TYPE_CHECK(right, RingElement) \
                  or not (right.parent() is self.parent()):
            return bin_op_c(self, right, operator.sub)
        return self._sub_(right)


    cdef RingElement _sub_c(self, RingElement right):
        """
        Subtraction dispatcher for RingElements.

        DO NOT OVERRIDE THIS FUNCTION.

        See extensive documentation at the top of element.pyx.
        """

        if HAS_DICTIONARY(self):   # fast check
            return self._sub_(right)
        else:
            # Must be a pure Pyrex class.
            return self._sub_c_impl(right)


    cdef RingElement _sub_c_impl(self, RingElement right):
        """
        Pyrex classes should override this function to implement subtraction.

        DO NOT CALL THIS FUNCTION DIRECTLY.

        See extensive documentation at the top of element.pyx.
        """
        # default implementation is to use the negation and addition
        # dispatchers:
        return self._add_c(right._neg_c())


    def _sub_(RingElement self, RingElement right):
        """
        Python classes should override this function to implement subtraction.

        See extensive documentation at the top of element.pyx.
        """
        return self._sub_c_impl(right)

    ##################################
    # Negation
    ##################################

    def __neg__(self):
        """
        Top-level negation operator for RingElements.

        See extensive documentation at the top of element.pyx.
        """
        # We ASSUME that self is a RingElement. No type checks.
        return (<RingElement>self)._neg_c()


    cdef RingElement _neg_c(self):
        """
        Negation dispatcher for RingElements.

        DO NOT OVERRIDE THIS FUNCTION.

        See extensive documentation at the top of element.pyx.
        """

        if HAS_DICTIONARY(self):   # fast check
            return self._neg_()
        else:
            # Must be a pure Pyrex class.
            return self._neg_c_impl()


    cdef RingElement _neg_c_impl(self):
        """
        Pyrex classes should override this function to implement negation.

        DO NOT CALL THIS FUNCTION DIRECTLY.

        See extensive documentation at the top of element.pyx.
        """
        # default implementation is to try multiplying by -1.
        return bin_op_c(-1, self, operator.mul)


    def _neg_(RingElement self):
        """
        Python classes should override this function to implement negation.

        See extensive documentation at the top of element.pyx.
        """
        return self._neg_c_impl()

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

cdef class FiniteFieldElement(FieldElement):
    pass

cdef class AlgebraElement(RingElement):
    cdef AlgebraElement _mul_left_scalar_c(self, RingElement m):
        """
        DO NOT OVERRIDE THIS FUNCTION.
        """
        if HAS_DICTIONARY(self):
            return self._mul_left_scalar_c_impl(right)
        else:
            return self._mul_left_scalar(right)

    cdef AlgebraElement _mul_left_scalar_c_impl(self, RingElement m):
        raise NotImplementedError

    def _mul_left_scalar(self, m):
        return self._mul_left_scalar_c_impl(m)


    cdef AlgebraElement _mul_right_scalar_c(self, RingElement m):
        """
        DO NOT OVERRIDE THIS FUNCTION.
        """
        if HAS_DICTIONARY(self):
            return self._mul_right_scalar_c_impl(right)
        else:
            return self._mul_right_scalar(right)

    cdef AlgebraElement _mul_right_scalar_c_impl(self, RingElement m):
        raise NotImplementedError

    def _mul_right_scalar(self, m):
        return self._mul_right_scalar_c_impl(m)



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
            if not (right.parent() is self.parent()):
                return cmp(self, right)
        except AttributeError:
            return cmp(self, right)
        return self._cmp_(right)

    def _cmp_(self, right):
        raise NotImplementedError


cdef int have_same_parent(left, right):
    # (We know at least one of the arguments is an Element. So if
    # their types are *equal* (fast to check) then they are both
    # Elements.  Otherwise use the slower test via PY_TYPE_CHECK.)
    return (PY_TYPE(left) is PY_TYPE(right)) or \
           (PY_TYPE_CHECK(right, Element) and PY_TYPE_CHECK(left, Element) and \
            (<Element>left)._parent is (<Element>right)._parent)



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
    if xp.has_coerce_map_from(yp):
        return x, xp(y)
    elif yp.has_coerce_map_from(xp):
        return yp(x), y
    raise TypeError, "unable to find a common canonical parent for x and y"

def x2_canonical_coercion(x, y):
    return canonical_coercion_c(x,y)

cdef x2_canonical_coercion_c(x, y):
    cdef int i
    xp = parent_c(x)
    yp = parent_c(y)
    if xp is yp:
        return x, y

    try:
        if PY_IS_NUMERIC(x):
            try:
                x = yp(x)
            except TypeError:
                y = xp(y)
                return x, y
        elif PY_IS_NUMERIC(y):
            try:
                y = xp(y)
            except TypeError:
                x = yp(x)
                return x, y
        else:
            i = 0
            try:
                x0 = x
                x = coerce(yp, x)
            except TypeError, msg:
                i = i + 1
            try:
                y = coerce(xp, y)
            except TypeError, msg:
                i = i + 1
            if i == 0:
                # Both succeed.  But we must be careful to take x before
                # it was coerced, or we end up *switching* to the parents,
                # which is no good.
                return x0, y
            elif i == 2:
                import  sage.rings.ring
                if isinstance(x, sage.rings.ring.Ring) or isinstance(y, sage.rings.ring.Ring):
                    raise TypeError, "you cannot +,*,/ a ring with a number."
                raise TypeError, "unable to find a common parent for %s (parent: %s) and %s (parent: %s)"%(x,xp, y, yp)
        return x, y
    except AttributeError:
        raise TypeError, "unable to find a common canonical parent"

########################################

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
    return functions.bin_op_c(x,y,op)
cdef bin_op_c(x, y, op):
    """
    Compute x op y, where coercion of x and y works according to
    SAGE's coercion rules.
    """
    #print "bin_op(%s,%s,%s)"%(x,y,op)   # debug
    if isinstance(y, InfinityElement):
        return op(y,x)
    if op == operator.mul and \
           isinstance(y, (\
                      ModuleElement,
                      AlgebraElement,
                      sage.modules.module.Module)) and \
           isinstance(x, (RingElement, int, long, float)):
        return op(y,x)
    try:
        #print 1, x, y, x.parent(), y.parent()
        x, y = canonical_coercion_c(x, y)
        #print 2, x, y, x.parent(), y.parent()
    except TypeError, mesg:
        try:
            return y._r_action(x)
        except AttributeError:
            raise TypeError, '%s: x parent: %s, y parent: %s'%(mesg, parent(x), parent(y))
        except TypeError:
            raise TypeError, "No right action defined"
        try:
            return x._l_action(y)
        except AttributeError:
            raise TypeError, mesg
        except TypeError:
            raise TypeError, "No left action defined"
            #raise TypeError, "No left action of %s on %s defined"%(x,y)

    return op(x,y)


def coerce_cmp(x,y):  # external interface to cmp_cdef
    return functions.cmp_c(x,y)

N = type(None)  # todo -- use Python/C API

cdef cmp_c(x, y):
    tx = type(x); ty = type(y)
    if (tx == N and ty != N) or (tx != N and ty == N):
        return -1
    elif isinstance(y, InfinityElement):
        return -y.__cmp__(x)

    xp = parent(x)
    yp = parent(y)
    if xp is yp:
        return __builtin__.cmp(x,y)

    cdef int fails
    fails = 0
    if isinstance(x, (int, long)):

        return __builtin__.cmp(yp(x), y)

    elif isinstance(y, (int, long)):

        return __builtin__.cmp(x, xp(y))

    else:

        fails = 0
        try:
            x0 = x
            x = coerce(yp, x)
        except (TypeError, ValueError):
        #    print "coercing %s to %s failed"%(x0,yp)
            fails = fails + 1
        #else:
        #    print "coercion %s to %s suceeded with %s (parent=%s)"%(x0,yp,x,parent(x))

        try:
            y0 = y
            y = coerce(xp, y)
        except (TypeError, ValueError):
        #    print "coercing %s to %s failed"%(y0,xp)
            fails = fails + 1
        #else:
        #    print "coercion %s to %s suceeded with %s (parent=%s)"%(y0,xp,y,parent(y))


        if fails == 0:
            assert (parent(x0) is parent(y))  # debug
            return __builtin__.cmp(x0,y)

        elif fails == 2:

            return -1

        else:
            if not (parent(x) is parent(y)):
                raise RuntimeError, "There is a bug in coercion: x=%s (parent=%s), y=%s (parent=%s)"%(x, parent(x), y, parent(y))
            return __builtin__.cmp(x,y)


def x_canonical_coercion(x, y):
    return x_canonical_coercion_c(x,y)
cdef x_canonical_coercion_c(x, y):
    cdef int i
    xp = parent(x)
    yp = parent(y)
    if xp is yp:
        return x, y

    try:
        if x.__class__ in [int, long, float, complex]:
            try:
                x = yp(x)
            except TypeError:
                y = x.__class__(y)
                return x, y
        elif y.__class__ in [int, long, float, complex]:
            try:
                y = xp(y)
            except TypeError:
                x = y.__class__(x)
                return x, y
        else:
            i = 0
            try:
                x0 = x
                x = coerce(yp, x)
            except TypeError, msg:
                i = i + 1
            try:
                y = coerce(xp, y)
            except TypeError, msg:
                i = i + 1
            if i == 0:
                # Both succeed.  But we must be careful to take x before
                # it was coerced, or we end up *switching* to the parents,
                # which is no good.
                return x0, y
            elif i == 2:
                import  sage.rings.ring
                if isinstance(x, sage.rings.ring.Ring) or isinstance(y, sage.rings.ring.Ring):
                    raise TypeError, "you cannot +,*,/ a ring with a number."
                raise TypeError, "unable to find a common parent for %s (parent: %s) and %s (parent: %s)"%(x,xp, y, yp)
        return x, y
    except AttributeError:
        raise TypeError, "unable to find a common canonical parent"

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

