r"""
Elements

AUTHORS:
   -- David Harvey (2006-10-16): changed CommutativeAlgebraElement to derive
   from CommutativeRingElement instead of AlgebraElement
   -- David Harvey (2006-10-29): implementation and documentation of new
   arithmetic architecture
   -- William Stein (2006-11): arithmetic architecture -- pushing it through to completion.
   -- Gonzalo Tornaria (2007-06): recursive base extend for coercion -- lots of tests


\subsection{The Abstract Element Class Hierarchy}
This is the abstract class hierarchy, i.e., these are all
abstract base classes.
\begin{verbatim}
SageObject
    Element
        ModuleElement
            RingElement
                CommutativeRingElement
                    IntegralDomainElement
                        DedekindDomainElement
                            PrincipalIdealDomainElement
                                EuclideanDomainElement
                    FieldElement
                        FiniteFieldElement
                    CommutativeAlgebraElement
                AlgebraElement   (note -- can't derive from module, since no multiple inheritence)
                    CommutativeAlgebra ??? (should be removed from element.pxd)
                    Matrix
                InfinityElement
                    PlusInfinityElement
                    MinusInfinityElement
            AdditiveGroupElement
            Vector

        MonoidElement
            MultiplicativeGroupElement

\end{verbatim}

\subsection{How to Define a New Element Class}
Elements typically define a method \code{_new_c}, e.g.,
\begin{verbatim}
    cdef _new_c(self, defining data):
        cdef FreeModuleElement_generic_dense x
        x = PY_NEW(FreeModuleElement_generic_dense)
        x._parent = self._parent
        x._entries = v
\end{verbatim}
that creates a new sibling very quickly from defining data
with assumed properties.

SAGE has a special system in place for handling arithmetic operations
for all Element subclasses. There are various rules that must be
followed by both arithmetic implementors and callers.

A quick summary for the impatient:
\begin{itemize}
 \item To implement addition for any Element class, override def _add_().
 \item If you want to add x and y, whose parents you know are IDENTICAL,
   you may call _add_(x, y). This will be the fastest way to guarantee
   that the correct implementation gets called. Of course you can still
   always use "x + y".
\end{itemize}

Now in more detail. The aims of this system are to provide (1) an efficient
calling protocol from both python and cython, (2) uniform coercion semantics
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
   necessary coercions have been performed, it calls _add_ to dispatch to
   the correct underlying addition implementation.

   Note that although this function is declared as def, it doesn't have the
   usual overheads associated with python functions (either for the caller
   or for __add__ itself). This is because python has optimised calling
   protocols for such special functions.

{\bf def RingElement._add_}

   This is the function you should override to implement addition in a
   python subclass of RingElement.

   WARNING: if you override this in a *pyrex* class, it won't get called.
   You should override _add_ instead. It is especially important to
   keep this in mind whenever you move a class down from python to pyrex.

   The two arguments to this function are guaranteed to have the
   SAME PARENT. Its return value MUST have the SAME PARENT as its
   arguments.

   If you want to add two objects from python, and you know that their
   parents are the same object, you are encouraged to call this function
   directly, instead of using "x + y".

   The default implementation of this function is to call _add_,
   so if no-one has defined a python implementation, the correct pyrex
   implementation will get called.

{\bf cpdef RingElement._add_}

   This is the function you should override to implement addition in a
   pyrex subclass of RingElement.

   The two arguments to this function are guaranteed to have the
   SAME PARENT. Its return value MUST have the SAME PARENT as its
   arguments.

   The default implementation of this function is to raise a
   NotImplementedError, which will happen if no-one has supplied
   implementations of either _add_.


For speed, there are also {\bf inplace} version of the arithmetic commands.
DD NOT call them directly, they may mutate the object and will be called
when and only when it has been determined that the old object will no longer
be accessible from the calling function after this operation.

{\bf def RingElement._iadd_}

   This is the function you should override to inplace implement addition
   in a python subclass of RingElement.

   The two arguments to this function are guaranteed to have the
   SAME PARENT. Its return value MUST have the SAME PARENT as its
   arguments.

   The default implementation of this function is to call _add_,
   so if no-one has defined a python implementation, the correct pyrex
   implementation will get called.
"""


cdef extern from *:
    ctypedef struct RefPyObject "PyObject":
        int ob_refcnt




##################################################################
# Generic element, so all this functionality must be defined
# by any element.  Derived class must call __init__
##################################################################

include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"
include "../ext/python.pxi"

import operator
import sys
import traceback

cdef MethodType
from types import MethodType

from sage.structure.parent      cimport Parent

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

def py_scalar_to_element(py):
    from sage.rings.integer_ring import ZZ
    from sage.rings.real_double import RDF
    from sage.rings.complex_double import CDF
    if PyInt_Check(py) or PyLong_Check(py):
        return ZZ(py)
    elif PyFloat_Check(py):
        return RDF(py)
    elif PyComplex_Check(py):
        return CDF(py)
    else:
        raise TypeError, "Not a scalar"

def is_Element(x):
    """
    Return True if x is of type Element.

    EXAMPLES:
        sage: from sage.structure.element import is_Element
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
        #if parent is None:
        #    raise RuntimeError, "bug -- can't set parent to None"
        self._parent = parent

    def _set_parent(self, parent):
        r"""
        INPUT:
            parent -- a SageObject
        """
        self._parent = parent

    cdef _set_parent_c(self, Parent parent):
        self._parent = parent

    def _make_new_with_parent_c(self, Parent parent):
        self._parent = parent
        return self

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

    cdef base_extend_c(self, Parent R):
        if HAS_DICTIONARY(self):
            return self.base_extend(R)
        else:
            return self.base_extend_c_impl(R)

    cdef base_extend_c_impl(self, Parent R):
        cdef Parent V
        V = self._parent.base_extend(R)
        return (<Parent>V)._coerce_c(self)

    def base_extend(self, R):
        return self.base_extend_c_impl(R)

    def base_ring(self):
        """
        Returns the base ring of this element's parent (if that makes sense).
        """
        return self._parent.base_ring()

    def category(self):
        from sage.categories.category_types import Elements
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


    def subs(self, in_dict=None, **kwds):
        """
        Substitutes given generators with given values while not touching
        other generators. This is a generic wrapper around __call__.
        The syntax is meant to be compatible with the corresponding method
        for symbolic expressions.

        INPUT:
            in_dict -- (optional) dictionary of inputs
            **kwds  -- named parameters

        OUTPUT:
            new object if substitution is possible, otherwise self.

        EXAMPLES:
            sage: x, y = PolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.subs({x:5})
            25*y^2 + y + 30
            sage: f.subs(x=5)
            25*y^2 + y + 30
            sage: (1/f).subs(x=5)
            1/(25*y^2 + y + 30)
            sage: Integer(5).subs(x=4)
            5
        """
        if not hasattr(self,'__call__'):
            return self
        parent=self._parent
        from sage.structure.parent_gens import is_ParentWithGens
        if not is_ParentWithGens(parent):
            return self
        variables=[]
        # use "gen" instead of "gens" as a ParentWithGens is not
        # required to have the latter
        for i in xrange(0,parent.ngens()):
            gen=parent.gen(i)
            if kwds.has_key(str(gen)):
                variables.append(kwds[str(gen)])
            elif in_dict and in_dict.has_key(gen):
                variables.append(in_dict[gen])
            else:
                variables.append(gen)
        return self(*variables)

    def n(self, prec=None, digits=None):
        """
        Return a numerical approximation of x with at least prec bits of
        precision.

        EXAMPLES:
            sage: (2/3).n()
            0.666666666666667
            sage: a = 2/3
            sage: pi.n(digits=10)
            3.141592654
            sage: pi.n(prec=20)   # 20 bits
            3.1416
        """
        import sage.misc.functional
        return sage.misc.functional.numerical_approx(self, prec=prec, digits=digits)

    def substitute(self,in_dict=None,**kwds):
        """
        This is an alias for self.subs().

        INPUT:
            in_dict -- (optional) dictionary of inputs
            **kwds  -- named parameters

        OUTPUT:
            new object if substitution is possible, otherwise self.

        EXAMPLES:
            sage: x, y = PolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.substitute({x:5})
            25*y^2 + y + 30
            sage: f.substitute(x=5)
            25*y^2 + y + 30
            sage: (1/f).substitute(x=5)
            1/(25*y^2 + y + 30)
            sage: Integer(5).substitute(x=4)
            5
         """
        return self.subs(in_dict,**kwds)




    def __xor__(self, right):
        raise RuntimeError, "Use ** for exponentiation, not '^', which means xor\n"+\
              "in Python, and has the wrong precedence."

    def __pos__(self):
        return self

    def _coeff_repr(self, no_space=True):
        if self._is_atomic():
            s = repr(self)
        else:
            s = "(%s)"%repr(self)
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
        """
        Return True if and only if parenthesis are not required when
        *printing* out any of $x - s$, $x + s$, $x^s$ and $x/s$.

        EXAMPLES:
            sage: n = 5; n._is_atomic()
            True
            sage: n = x+1; n._is_atomic()
            False
        """
        if self._parent.is_atomic_repr():
            return True
        s = str(self)
        return s.find("+") == -1 and s.find("-") == -1 and s.find(" ") == -1

    def __nonzero__(self):
        """
        Return True if self does not equal self.parent()(0).
        """
        return self != self._parent(0)

    def is_zero(self):
        """
        Return True if self equals self.parent()(0). The default
        implementation is to fall back to 'not self.__nonzero__'.

        NOTE: Do not re-implement this method in your subclass but
        implement __nonzero__ instead.
        """
        return not self

    def _cmp_(left, right):
        return left._cmp_c_impl(right)

    cdef _cmp(left, right):
        """
        Compare left and right.
        """
        global coercion_model
        cdef int r
        if not have_same_parent(left, right):
            if left is None or left is Ellipsis:
                return -1
            elif right is None or right is Ellipsis:
                return 1
            try:
                _left, _right = coercion_model.canonical_coercion(left, right)
                if PY_IS_NUMERIC(_left):
                    return cmp(_left, _right)
                else:
                    return _left._cmp_(_right)
            except TypeError:
                r = cmp(type(left), type(right))
                if r == 0:
                    r = -1
                return r
        else:
            if HAS_DICTIONARY(left):
                left_cmp = left._cmp_
                if PY_TYPE_CHECK(left_cmp, MethodType):
                    # it must have been overriden
                    return left_cmp(right)

            return left._cmp_c_impl(right)

    def _richcmp_(left, right, op):
        return left._richcmp(right, op)

    cdef _richcmp(left, right, int op):
        """
        Compare left and right, according to the comparison operator op.
        """
        global coercion_model
        cdef int r
        if not have_same_parent(left, right):
            if left is None or left is Ellipsis:
                return _rich_to_bool(op, -1)
            elif right is None or right is Ellipsis:
                return _rich_to_bool(op, 1)
            try:
                _left, _right = coercion_model.canonical_coercion(left, right)
                if PY_IS_NUMERIC(_left):
                    return _rich_to_bool(op, cmp(_left, _right))
                else:
                    return _left._richcmp_(_right, op)
            except (TypeError, NotImplementedError):
                r = cmp(type(left), type(right))
                if r == 0:
                    r = -1
                # Often things are compared against 0 (or 1), even when there
                # is not a cannonical coercion ZZ -> other
                # Things should implement and/or use __nonzero__ and is_one()
                # but we can't do that here as that calls this.
                # The old coercion model would declare a coercion if 0 went in.
                # (Though would fail with a TypeError for other values, thus
                # contaminating the _has_coerce_map_from cache.)
                from sage.rings.integer import Integer
                try:
                    if PY_TYPE_CHECK(left, Element) and isinstance(right, (int, float, Integer)) and not right:
                        right = (<Element>left)._parent(right)
                    elif PY_TYPE_CHECK(right, Element) and isinstance(left, (int, float, Integer)) and not left:
                        left = (<Element>right)._parent(left)
                    else:
                        return _rich_to_bool(op, r)
                except TypeError:
                    return _rich_to_bool(op, r)

        if HAS_DICTIONARY(left):   # fast check
            left_cmp = left.__cmp__
            if PY_TYPE_CHECK(left_cmp, MethodType):
                # it must have been overriden
                return _rich_to_bool(op, left_cmp(right))

        return left._richcmp_c_impl(right, op)

    cdef _rich_to_bool(self, int op, int r):
        return _rich_to_bool(op, r)

    ####################################################################
    # For a derived Cython class, you **must** put the following in
    # your subclasses, in order for it to take advantage of the
    # above generic comparison code.  You must also define
    # either _cmp_c_impl (if your subclass is totally ordered),
    # _richcmp_c_impl (if your subclass is partially ordered), or both
    # (if your class has both a total order and a partial order;
    # then the total order will be available with cmp(), and the partial
    # order will be available with the relation operators; in this case
    # you must also define __cmp__ in your subclass).
    # This is simply how Python works.
    #
    # For a *Python* class just define __cmp__ as always.
    # But note that when this gets called you can assume that
    # both inputs have identical parents.
    #
    # If your __cmp__ methods are not getting called, verify that the
    # canonical_coercion(x,y) is not throwing errors.
    #
    ####################################################################
    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    ####################################################################
    # If your subclass has both a partial order (available with the
    # relation operators) and a total order (available with cmp()),
    # you **must** put the following in your subclass.
    #
    # Note that in this case the total order defined by cmp() will
    # not properly respect coercions.
    ####################################################################
    def __cmp__(left, right):
        return (<Element>left)._cmp(right)

    cdef _richcmp_c_impl(left, Element right, int op):
        if (<Element>left)._richcmp_c_impl == Element._richcmp_c_impl and \
               (<Element>left)._cmp_c_impl == Element._cmp_c_impl:
            # Not implemented, try some basic defaults
            if op == Py_EQ:
                 return left is right
            elif op == Py_NE:
                return left is not right
        return left._rich_to_bool(op, left._cmp_c_impl(right))

    cdef int _cmp_c_impl(left, Element right) except -2:
        ### For derived SAGEX code, you *MUST* ALSO COPY the __richcmp__ above
        ### into your class!!!  For Python code just use __cmp__.
        raise NotImplementedError, "BUG: sort algorithm for elements of '%s' not implemented"%right.parent()

cdef inline bint _rich_to_bool(int op, int r):
    if op == Py_LT:  #<
        return (r  < 0)
    elif op == Py_EQ: #==
        return (r == 0)
    elif op == Py_GT: #>
        return (r  > 0)
    elif op == Py_LE: #<=
        return (r <= 0)
    elif op == Py_NE: #!=
        return (r != 0)
    elif op == Py_GE: #>=
        return (r >= 0)


def is_ModuleElement(x):
    """
    Return True if x is of type ModuleElement.

    This is even faster than using isinstance inline.

    EXAMPLES:
        sage: from sage.structure.element import is_ModuleElement
        sage: is_ModuleElement(2/3)
        True
        sage: is_ModuleElement((QQ^3).0)
        True
        sage: is_ModuleElement('a')
        False
    """
    return IS_INSTANCE(x, ModuleElement)


cdef class ModuleElement(Element):
    """
    Generic element of a module.
    """

    ##################################################
    # Addition
    ##################################################
    def __add__(left, right):
        """
        Top-level addition operator for ModuleElements.

        See extensive documentation at the top of element.pyx.
        """
        # Try fast pathway if they are both ModuleElements and the parents
        # match.

        # (We know at least one of the arguments is a ModuleElement. So if
        # their types are *equal* (fast to check) then they are both
        # ModuleElements. Otherwise use the slower test via PY_TYPE_CHECK.)
        if have_same_parent(left, right):
            # If we hold the only references to this object, we can
            # safely mutate it. NOTE the threshold is different by one
            # for __add__ and __iadd__.
            if  (<RefPyObject *>left).ob_refcnt < inplace_threshold:
                return (<ModuleElement>left)._iadd_(<ModuleElement>right)
            else:
                return (<ModuleElement>left)._add_(<ModuleElement>right)

        global coercion_model
        return coercion_model.bin_op(left, right, add)

    cpdef ModuleElement _add_(left, ModuleElement right):
        raise TypeError, arith_error_message(left, right, add)

    def __iadd__(ModuleElement self, right):
        if have_same_parent(self, right):
            if  (<RefPyObject *>self).ob_refcnt <= inplace_threshold:
                return self._iadd_(<ModuleElement>right)
            else:
                return self._add_(<ModuleElement>right)
        else:
            global coercion_model
            return coercion_model.bin_op(self, right, iadd)

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        return self._add_(right)

    ##################################################
    # Subtraction
    ##################################################

    def __sub__(left, right):
        """
        Top-level subtraction operator for ModuleElements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent(left, right):
            if  (<RefPyObject *>left).ob_refcnt < inplace_threshold:
                return (<ModuleElement>left)._isub_(<ModuleElement>right)
            else:
                return (<ModuleElement>left)._sub_(<ModuleElement>right)
        global coercion_model
        return coercion_model.bin_op(left, right, sub)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        # default implementation is to use the negation and addition
        # dispatchers:
        return left._add_(-right)

    def __isub__(ModuleElement self, right):
        if have_same_parent(self, right):
            if  (<RefPyObject *>self).ob_refcnt <= inplace_threshold:
                return self._isub_(<ModuleElement>right)
            else:
                return self._sub_(<ModuleElement>right)
        else:
            global coercion_model
            return coercion_model.bin_op(self, right, isub)

    cpdef ModuleElement _isub_(self, ModuleElement right):
        return self._sub_(right)

    ##################################################
    # Negation
    ##################################################

    def __neg__(self):
        """
        Top-level negation operator for ModuleElements, which
        may choose to implement _neg_ rather than __neg__ for
        consistancy.
        """
        return self._neg_()

    cpdef ModuleElement _neg_(self):
        # default implementation is to try multiplying by -1.
        global coercion_model
        if self._parent._base is None:
            return coercion_model.bin_op(-1, self, mul)
        else:
            return coercion_model.bin_op(self._parent._base(-1), self, mul)

    ##################################################
    # Module element multiplication (scalars, etc.)
    ##################################################
    def __mul__(left, right):
        if have_same_parent(left, right):
            raise arith_error_message(left, right, mul)
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, mul)

    def __imul__(left, right):
        if have_same_parent(left, right):
             raise TypeError
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, imul)

    # rmul -- left * self
    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        Default module left scalar multiplication, which is to try to
        canonically coerce the scalar to the integers and do that
        multiplication, which is always defined.

        The coercion model will catch this error and create the
        appropriate action.
        """
        raise NotImplementedError

    # lmul -- self * right

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Default module left scalar multiplication, which is to try to
        canonically coerce the scalar to the integers and do that
        multiplication, which is always defined.

        The coercion model will catch this error and create the
        appropriate action.
        """
        raise NotImplementedError

    cpdef ModuleElement _ilmul_(self, RingElement right):
        return self._lmul_(right)

    cdef RingElement coerce_to_base_ring(self, x):
        if PY_TYPE_CHECK(x, Element) and (<Element>x)._parent is self._parent._base:
            return x
        try:
            return self._parent._base._coerce_c(x)
        except AttributeError:
            return self._parent._base(x)

    ##################################################
    # Other properties
    ##################################################
    def order(self):              ### DO NOT OVERRIDE THIS!!! Instead, override additive_order.
        """
        Return the additive order of self.
        """
        return self.additive_order()

    def additive_order(self):
        """
        Return the additive order of self.
        """
        raise NotImplementedError

########################################################################
# Monoid
########################################################################

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
    def __mul__(left, right):
        """
        Top-level multiplication operator for monoid elements.
        See extensive documentation at the top of element.pyx.
        """
        global coercion_model
        if have_same_parent(left, right):
            return (<MonoidElement>left)._mul_(<MonoidElement>right)
        try:
            return coercion_model.bin_op(left, right, mul)
        except TypeError, msg:
            if isinstance(left, (int, long)) and left==1:
                return right
            elif isinstance(right, (int, long)) and right==1:
                return left
            raise TypeError, msg


    cpdef MonoidElement _mul_(left, MonoidElement right):
        """
        Pyrex classes should override this function to implement multiplication.
        See extensive documentation at the top of element.pyx.
        """
        raise TypeError

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
        """
        Return the (integral) power of self.
        """
        if dummy is not None:
            raise RuntimeError, "__pow__ dummy argument not used"
        return generic_power_c(self,n,None)

    def __nonzero__(self):
        return True

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

    cpdef ModuleElement _rmul_(self, RingElement left):
        return self._lmul_(left)

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Default module left scalar multiplication, which is to try to
        canonically coerce the scalar to the integers and do that
        multiplication, which is always defined.

        The coercion model will catch this error and create the
        appropriate action.
        """
        raise NotImplementedError

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

    def __div__(left, right):
        if have_same_parent(left, right):
            return left._div_(right)
        global coercion_model
        return coercion_model.bin_op(left, right, div)

    cpdef MultiplicativeGroupElement _div_(self, MultiplicativeGroupElement right):
        """
        Pyrex classes should override this function to implement division.
        See extensive documentation at the top of element.pyx.
        """
        return self * ~right

    def __invert__(self):
        if self.is_one():
            return self
        return 1/self


def is_RingElement(x):
    """
    Return True if x is of type RingElement.
    """
    return IS_INSTANCE(x, RingElement)

cdef class RingElement(ModuleElement):
    ##################################################
    def is_one(self):
        return self == self._parent(1)

    ##################################
    # Multiplication
    ##################################

    cpdef ModuleElement _lmul_(self, RingElement right):
        # We raise an error to invoke the default action of coercing into self
        raise NotImplementedError, "parents %s %s %s" % (parent_c(self), parent_c(right), parent_c(self) is parent_c(right))

    cpdef ModuleElement _rmul_(self, RingElement left):
        # We raise an error to invoke the default action of coercing into self
        raise NotImplementedError

    def __mul__(self, right):
        """
        Top-level multiplication operator for ring elements.
        See extensive documentation at the top of element.pyx.

        AUTHOR:

            Gonzalo Tornaria (2007-06-25) - write base-extending test cases and fix them

        TEST CASES:

            (scalar * vector)

            sage: x, y = var('x, y')

            sage: parent(ZZ(1)*vector(ZZ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(QQ(1)*vector(ZZ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(ZZ(1)*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(QQ(1)*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field

            sage: parent(QQ(1)*vector(ZZ[x],[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x](1)*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ(1)*vector(ZZ[x][y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x][y](1)*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ[x](1)*vector(ZZ[x][y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x][y](1)*vector(QQ[x],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ[y](1)*vector(ZZ[x][y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x][y](1)*vector(QQ[y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(ZZ[x](1)*vector(ZZ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(ZZ[x](1)*vector(QQ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(QQ[x](1)*vector(ZZ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(QQ[x](1)*vector(QQ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'

            (scalar * matrix)

            sage: parent(ZZ(1)*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: parent(QQ(1)*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(ZZ(1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(QQ(1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: parent(QQ(1)*matrix(ZZ[x],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x](1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ(1)*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x][y](1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ[x](1)*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x][y](1)*matrix(QQ[x],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ[y](1)*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ[x][y](1)*matrix(QQ[y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(ZZ[x](1)*matrix(ZZ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(ZZ[x](1)*matrix(QQ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'
            sage: parent(QQ[x](1)*matrix(ZZ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(QQ[x](1)*matrix(QQ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'

        """
        # Try fast pathway if they are both RingElements and the parents match.
        # (We know at least one of the arguments is a RingElement. So if their
        # types are *equal* (fast to check) then they are both RingElements.
        # Otherwise use the slower test via PY_TYPE_CHECK.)
        if have_same_parent(self, right):
            if  (<RefPyObject *>self).ob_refcnt < inplace_threshold:
                return (<RingElement>self)._imul_(<RingElement>right)
            else:
                return (<RingElement>self)._mul_(<RingElement>right)
        global coercion_model
        return coercion_model.bin_op(self, right, mul)

    cpdef RingElement _mul_(self, RingElement right):
        """
        Pyrex classes should override this function to implement multiplication.
        See extensive documentation at the top of element.pyx.
        """
        raise TypeError, arith_error_message(self, right, mul)

    def __imul__(left, right):
        if have_same_parent(left, right):
            if  (<RefPyObject *>left).ob_refcnt <= inplace_threshold:
                return (<RingElement>left)._imul_(<RingElement>right)
            else:
                return (<RingElement>left)._mul_(<RingElement>right)

        global coercion_model
        return coercion_model.bin_op(left, right, imul)

    cpdef RingElement _imul_(RingElement self, RingElement right):
        return self._mul_(right)

    def __pow__(self, n, dummy):
        """
        Retern the (integral) power of self.

        EXAMPLE:
            sage: a = Integers(389)['x']['y'](37)
            sage: a^2
            202
            sage: a^388
            1
            sage: a^(2^120)
            81
            sage: a^0
            1
            sage: a^1 == a
            True
            sage: a^2 * a^3 == a^5
            True
            sage: (a^3)^2 == a^6
            True
            sage: a^57 * a^43 == a^100
            True
            sage: a^(-1) == 1/a
            True
            sage: a^200 * a^(-64) == a^136
            True

        TESTS:
            sage: 2r**(SR(2)-1-1r)
            1
            sage: 2r**(SR(1/2))
            Traceback (most recent call last):
            ...
            NotImplementedError: non-integral exponents not supported
        """
        if dummy is not None:
            raise RuntimeError, "__pow__ dummy argument not used"
        return generic_power_c(self,n,None)

    ##################################
    # Division
    ##################################

    def __truediv__(self, right):
        # in sage all divs are true
        if not PY_TYPE_CHECK(self, Element):
            return coercion_model.bin_op(self, right, div)
        return self.__div__(right)

    def __div__(self, right):
        """
        Top-level multiplication operator for ring elements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent(self, right):
            if  (<RefPyObject *>self).ob_refcnt < inplace_threshold:
                return (<RingElement>self)._idiv_(<RingElement>right)
            else:
                return (<RingElement>self)._div_(<RingElement>right)
        global coercion_model
        return coercion_model.bin_op(self, right, div)

    cpdef RingElement _div_(self, RingElement right):
        """
        Pyrex classes should override this function to implement division.
        See extensive documentation at the top of element.pyx.
        """
        try:
            return self._parent.fraction_field()(self, right)
        except AttributeError:
            if not right:
                raise ZeroDivisionError, "Cannot divide by zero"
            else:
                raise TypeError, arith_error_message(self, right, div)

    def __idiv__(self, right):
        """
        Top-level division operator for ring elements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent(self, right):
            if  (<RefPyObject *>self).ob_refcnt <= inplace_threshold:
                return (<RingElement>self)._idiv_(<RingElement>right)
            else:
                return (<RingElement>self)._div_(<RingElement>right)
        global coercion_model
        return coercion_model.bin_op(self, right, idiv)

    cpdef RingElement _idiv_(RingElement self, RingElement right):
        return self._div_(right)

    def __invert__(self):
        if self.is_one():
            return self
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

    def is_nilpotent(self):
        """
        Return True if self is nilpotent, i.e., some power of self
        is 0.
        """
        if self.is_unit():
            return False
        if self.is_zero():
            return True
        raise NotImplementedError

    def abs(self):
        """
        Return the absolute value of self.  (This just calls the __abs__
        method, so it is equivalent to the abs() built-in function.)

        EXAMPLES:
            sage: RR(-1).abs()
            1.00000000000000
            sage: ZZ(-1).abs()
            1
            sage: CC(I).abs()
            1.00000000000000
            sage: Mod(-15, 37).abs()
            Traceback (most recent call last):
            ...
            ArithmeticError: absolute valued not defined on integers modulo n.
        """
        return self.__abs__()



def is_CommutativeRingElement(x):
    """
    Return True if x is of type CommutativeRingElement.
    """
    return IS_INSTANCE(x, CommutativeRingElement)

cdef class CommutativeRingElement(RingElement):
    def _im_gens_(self, codomain, im_gens):
        if len(im_gens) == 1 and self._parent.gen(0) == 1:
            return codomain(self)
        raise NotImplementedError

    def inverse_mod(self, I):
        r"""
        Return an inverse of self modulo the ideal $I$, if defined,
        i.e., if $I$ and self together generate the unit ideal.
        """
        raise NotImplementedError

    def divides(self, x):
        """
        Return True if self divides x.

        EXAMPLES:
            sage: P.<x> = PolynomialRing(QQ)
            sage: x.divides(x^2)
            True
            sage: x.divides(x^2+2)
            False
            sage: (x^2+2).divides(x)
            False
            sage: P.<x> = PolynomialRing(ZZ)
            sage: x.divides(x^2)
            True
            sage: x.divides(x^2+2)
            False
            sage: (x^2+2).divides(x)
            False
        """
        return (x % self) == 0

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
            2*y^2 + 2*y*z + 2*z^2

        Notice above that $x$ is eliminated.  In the next example,
        both $y$ and $z$ are eliminated.

            sage: (x^2 + y^2 + z^2).mod( (x - y, y - z) )
            3*z^2
            sage: f = (x^2 + y^2 + z^2)^2; f
            x^4 + 2*x^2*y^2 + y^4 + 2*x^2*z^2 + 2*y^2*z^2 + z^4
            sage: f.mod( (x - y, y - z) )
            9*z^4

        In this example $y$ is eliminated.
            sage: (x^2 + y^2 + z^2).mod( (x^3, y - z) )
            x^2 + 2*z^2
        """
        from sage.rings.all import is_Ideal
        if not is_Ideal(I) or not I.ring() is self._parent:
            I = self._parent.ideal(I)
            #raise TypeError, "I = %s must be an ideal in %s"%(I, self.parent())
        return I.reduce(self)

cdef class Vector(ModuleElement):
    cdef bint is_sparse_c(self):
        raise NotImplementedError

    cdef bint is_dense_c(self):
        raise NotImplementedError

    def __imul__(left, right):
        if have_same_parent(left, right):
            return (<Vector>left)._dot_product_(<Vector>right)
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, imul)


    def __mul__(left, right):
        """
        Multiplication of vector by vector, matrix, or scalar

        AUTHOR:

            Gonzalo Tornaria (2007-06-21) - write test cases and fix them

        NOTE:

            scalar * vector is implemented (and tested) in class RingElement
            matrix * vector is implemented (and tested) in class Matrix

        TEST CASES:

            (vector * vector)

            sage: x, y = var('x, y')

            sage: parent(vector(ZZ,[1,2])*vector(ZZ,[1,2]))
            Integer Ring
            sage: parent(vector(ZZ,[1,2])*vector(QQ,[1,2]))
            Rational Field
            sage: parent(vector(QQ,[1,2])*vector(ZZ,[1,2]))
            Rational Field
            sage: parent(vector(QQ,[1,2])*vector(QQ,[1,2]))
            Rational Field

            sage: parent(vector(QQ,[1,2,3,4])*vector(ZZ[x],[1,2,3,4]))
            Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x],[1,2,3,4])*vector(QQ,[1,2,3,4]))
            Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ,[1,2,3,4])*vector(ZZ[x][y],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2,3,4])*vector(QQ,[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[x],[1,2,3,4])*vector(ZZ[x][y],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2,3,4])*vector(QQ[x],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[y],[1,2,3,4])*vector(ZZ[x][y],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2,3,4])*vector(QQ[y],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(ZZ[x],[1,2,3,4])*vector(ZZ[y],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ[x],[1,2,3,4])*vector(QQ[y],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ[x],[1,2,3,4])*vector(ZZ[y],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ[x],[1,2,3,4])*vector(QQ[y],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'

            (vector * matrix)

            sage: parent(vector(ZZ,[1,2])*matrix(ZZ,2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(vector(QQ,[1,2])*matrix(ZZ,2,2,[1,2,3,4]))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(ZZ,[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Vector space of dimension 2 over Rational Field

            sage: parent(vector(QQ,[1,2])*matrix(ZZ[x],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x],[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ,[1,2])*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[x],[1,2])*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2])*matrix(QQ[x],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[y],[1,2])*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2])*matrix(QQ[y],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(ZZ[x],[1,2])*matrix(ZZ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ[x],[1,2])*matrix(QQ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ[x],[1,2])*matrix(ZZ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ[x],[1,2])*matrix(QQ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'

            (vector * scalar)

            sage: parent(vector(ZZ,[1,2])*ZZ(1))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(vector(QQ,[1,2])*ZZ(1))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(ZZ,[1,2])*QQ(1))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2])*QQ(1))
            Vector space of dimension 2 over Rational Field

            sage: parent(vector(QQ,[1,2])*ZZ[x](1))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x],[1,2])*QQ(1))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ,[1,2])*ZZ[x][y](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2])*QQ(1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[x],[1,2])*ZZ[x][y](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2])*QQ[x](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[y],[1,2])*ZZ[x][y](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2])*QQ[y](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(ZZ[x],[1,2])*ZZ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ[x],[1,2])*QQ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ[x],[1,2])*ZZ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ[x],[1,2])*QQ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

        """
        if have_same_parent(left, right):
            return (<Vector>left)._dot_product_(<Vector>right)
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, mul)

    cpdef Element _dot_product_(Vector left, Vector right):
        raise TypeError, arith_error_message(left, right, mul)

    cpdef Vector _pairwise_product_(Vector left, Vector right):
        raise TypeError, "unsupported operation for '%s' and '%s'"%(parent_c(left), parent_c(right))

    def __div__(self, right):
        if PY_IS_NUMERIC(right):
            right = py_scalar_to_element(right)
        if PY_TYPE_CHECK(right, RingElement):
            # Let __mul__ do the job
            return self.__mul__(~right)
        if PY_TYPE_CHECK(right, Vector):
            try:
                W = right.parent().submodule([right])
                return W.coordinates(self)[0] / W.coordinates(right)[0]
            except ArithmeticError:
                if right.is_zero():
                    raise ZeroDivisionError, "division by zero vector"
                else:
                    raise ArithmeticError, "vector is not in free module"
        raise TypeError, arith_error_message(self, right, div)

    def _magma_convert_(self, magma):
        """
        EXAMPLES:
            sage: v = vector([1,2,3])
            sage: mv = magma(v); mv                     # optional - magma
            (1 2 3)
            sage: mv.Type()                             # optional - magma
            ModTupRngElt
            sage: mv.Parent()                           # optional - magma
            Full RSpace of degree 3 over Integer Ring

            sage: v = vector(QQ, [1/2, 3/4, 5/6])
            sage: mv = magma(v); mv                     # optional - magma
            (1/2 3/4 5/6)
            sage: mv.Type()                             # optional - magma
            ModTupFldElt
            sage: mv.Parent()                           # optional - magma
            Full Vector space of degree 3 over Rational Field

        A more demanding example:
            sage: R.<x,y,z> = QQ[]
            sage: v = vector([x^3, y, 2/3*z + x/y])
            sage: magma(v)                              # optional - magma
            (            x^3               y (2/3*y*z + x)/y)
            sage: magma(v).Parent()                     # optional - magma
            Full Vector space of degree 3 over Multivariate rational function field of rank 3 over Rational Field

        """
        V = magma(self._parent)
        v = []
        for x in self.list():
            v.append(x._magma_init_(magma))
        s = '%s![%s]'%(V.name(), ','.join(v))
        return magma(s)


def is_Vector(x):
    return IS_INSTANCE(x, Vector)

cdef class Matrix(AlgebraElement):

    cdef bint is_sparse_c(self):
        raise NotImplementedError

    cdef bint is_dense_c(self):
        raise NotImplementedError

    def __imul__(left, right):
        if have_same_parent(left, right):
            return (<Matrix>left)._matrix_times_matrix_(<Matrix>right)
        else:
            global coercion_model
            return coercion_model.bin_op(left, right, imul)

    def __mul__(left, right):
        """
        Multiplication of matrix by matrix, vector, or scalar

        AUTHOR:

            Gonzalo Tornaria (2007-06-25) - write test cases and fix them

        NOTE:

            scalar * matrix is implemented (and tested) in class RingElement
            vector * matrix is implemented (and tested) in class Vector

        TEST CASES:

            (matrix * matrix)

            sage: x, y = var('x, y')

            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(ZZ[x],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*matrix(QQ[x],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ[y],2,2,[1,2,3,4])*matrix(ZZ[x][y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*matrix(QQ[y],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*matrix(ZZ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*matrix(QQ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'
            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*matrix(ZZ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*matrix(QQ[y],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'

            (matrix * vector)

            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*vector(ZZ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(ZZ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(ZZ[x],[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(ZZ[x][y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*vector(ZZ[x][y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*vector(QQ[x],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ[y],2,2,[1,2,3,4])*vector(ZZ[x][y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*vector(QQ[y],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*vector(ZZ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*vector(QQ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*vector(ZZ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*vector(QQ[y],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'

            (matrix * scalar)

            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*ZZ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*ZZ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*ZZ[x](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*ZZ[x][y](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*ZZ[x][y](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*QQ[x](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ[y],2,2,[1,2,3,4])*ZZ[x][y](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ[x][y],2,2,[1,2,3,4])*QQ[y](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*ZZ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(ZZ[x],2,2,[1,2,3,4])*QQ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Rational Field'
            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*ZZ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(QQ[x],2,2,[1,2,3,4])*QQ[y](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

        """
        if have_same_parent(left, right):
            return (<Matrix>left)._matrix_times_matrix_(<Matrix>right)
        else:
            global coercion_model
            return coercion_model.bin_op(left, right, mul)

    cdef Vector _vector_times_matrix_(matrix_right, Vector vector_left):
        raise TypeError

    cdef Vector _matrix_times_vector_(matrix_left, Vector vector_right):
        raise TypeError

    cdef Matrix _matrix_times_matrix_(left, Matrix right):
        raise TypeError



def is_Matrix(x):
    return IS_INSTANCE(x, Matrix)

def is_IntegralDomainElement(x):
    """
    Return True if x is of type IntegralDomainElement.
    """
    return IS_INSTANCE(x, IntegralDomainElement)

cdef class IntegralDomainElement(CommutativeRingElement):
    def is_nilpotent(self):
        return self.is_zero()


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
        if not PY_TYPE_CHECK(right, Element) or not ((<Element>right)._parent is self._parent):
            return coercion_model.bin_op(self, right, lcm)
        return self._lcm(right)

    def gcd(self, right):
        """
        Returns the gcd of self and right, or 0 if both are 0.
        """
        if not PY_TYPE_CHECK(right, Element) or not ((<Element>right)._parent is self._parent):
            return coercion_model.bin_op(self, right, gcd)
        return self._gcd(right)

    def xgcd(self, right):
        r"""
        Return the extended gcd of self and other, i.e., elements $r, s, t$ such that
        $$
           r = s \cdot self + t \cdot other.
        $$

        Note: there is no guarantee on minimality of the cofactors.
        In the integer case, see documentation for Integer._xgcd() to
        obtain minimal cofactors.
        """
        if not PY_TYPE_CHECK(right, Element) or not ((<Element>right)._parent is self._parent):
            return coercion_model.bin_op(self, right, xgcd)
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

    def __divmod__(self, other):
        """
        Return the quotient and remainder of self divided by other.

        EXAMPLES:
            sage: divmod(5,3)
            (1, 2)
        """
        return self.quo_rem(other)

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
            sage: R.<x> = ZZ[]
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

    def __floordiv__(self, other):
        return self / other

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
        return not not self

    def _gcd(self, FieldElement other):
        """
        Return the greatest common divisor of self and other.
        """
        if self.is_zero() and other.is_zero():
            return self
        else:
            return self._parent(1)

    def _lcm(self, FieldElement other):
        """
        Return the least common multiple of self and other.
        """
        if self.is_zero() and other.is_zero():
            return self
        else:
            return self._parent(1)

    def _xgcd(self, FieldElement other):
        R = self._parent
        if not self.is_zero():
            return R(1), ~self, R(0)
        elif not other.is_zero():
            return R(1), R(0), ~self
        else: # both are 0
            return self, self, self


    def quo_rem(self, right):
        if not isinstance(right, FieldElement) or not (right._parent is self._parent):
            right = self.parent()(right)
        return self/right, 0

    def divides(self, FieldElement other):
        r"""
        Check whether self divides other, for field elements.

        Since this is a field, all values divide all other values,
        except that zero does not divide any non-zero values.

        EXAMPLES:
            sage: K.<rt3> = QQ[sqrt(3)]
            sage: K(0).divides(rt3)
            False
            sage: rt3.divides(K(17))
            True
            sage: K(0).divides(K(0))
            True
            sage: rt3.divides(K(0))
            True
        """
        if not (other._parent is self._parent):
            other = self.parent()(other)
        return bool(self) or other.is_zero()

## def is_FiniteFieldElement(x):
##     """
##     Return True if x is of type FiniteFieldElement.
##     """
##     return IS_INSTANCE(x, FiniteFieldElement)

cdef class FiniteFieldElement(FieldElement):

    def _im_gens_(self, codomain, im_gens):
        """
        Used for applying homomorphisms of finite fields.

        EXAMPLES:
            sage: k.<a> = FiniteField(73^2, 'a')
            sage: K.<b> = FiniteField(73^4, 'b')
            sage: phi = k.hom([ b^(73*73+1) ])
            sage: phi(0)
            0
            sage: phi(a)
            7*b^3 + 13*b^2 + 65*b + 71

            sage: phi(a+3)
            7*b^3 + 13*b^2 + 65*b + 1
        """
        ## NOTE: see the note in sage/rings/number_field_element.pyx,
        ## in the comments for _im_gens_ there -- something analogous
        ## applies here.
        return codomain(self.polynomial()(im_gens[0]))

    def minpoly(self,var='x'):
        """
        Returns the minimal polynomial of this element
        (over the corresponding prime subfield).
        EXAMPLES:
            sage: k.<a> = FiniteField(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: b=a**20;p=b.charpoly("x");p
            x^2 + 15*x + 4
            sage: factor(p)
            (x + 17)^2
            sage: b.minpoly('x')
            x + 17
        """
        p=self.charpoly(var);
        for q in p.factor():
            if q[0](self)==0:
                return q[0]
        # This shouldn't be reached, but you never know!
        raise ArithmeticError("Could not find the minimal polynomial")

        ## We have two names for the same method
        ## for compatibility with sage.matrix
    def minimal_polynomial(self,var='x'):
        """
        Returns the minimal polynomial of this element
        (over the corresponding prime subfield).
        EXAMPLES:
            sage: k.<a> = FiniteField(3^4)
            sage: parent(a)
            Finite Field in a of size 3^4
            sage: b=a**20;p=charpoly(b,"y");p
            y^4 + 2*y^2 + 1
            sage: factor(p)
            (y^2 + 1)^2
            sage: b.minimal_polynomial('y')
            y^2 + 1
        """
        return self.minpoly(var)

    def vector(self, reverse=False):
        r"""
        See \code{_vector_}.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 + 1
            sage: e.vector() # random-ish error message
            doctest:1: DeprecationWarning:The function vector is replaced by _vector_.
            (1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        from sage.misc.misc import deprecation
        deprecation("The function vector is replaced by _vector_.")
        return self._vector_()

    def _vector_(self, reverse=False):
        """
        Return a vector in self.parent().vector_space() matching
        self. The most significant bit is to the right.

        INPUT:
            reverse -- reverse the order of the bits
                       from little endian to big endian.

        EXAMPLES:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: k(v)
            a^2 + 1

            sage: k.<a> = GF(3^16)
            sage: e = 2*a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: k(v)
            2*a^2 + 1

        You can also compute the vector in the other order:

            sage: e._vector_(reverse=True)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1)
        """
        #vector(foo) might pass in ZZ
        if PY_TYPE_CHECK(reverse, Parent):
            raise TypeError, "Base field is fixed to prime subfield."

        k = self.parent()

        v = self.polynomial().list()

        ret = [v[i] for i in range(len(v))]

        for i in range(k.degree() - len(ret)):
            ret.append(0)

        if reverse:
            ret = list(reversed(ret))
        return k.vector_space()(ret)

    def matrix(self, reverse=False):
        r"""
        See \code{_matrix_}.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 + 1
            sage: e.matrix() # random-ish error message
            doctest:1: DeprecationWarning:The function matrix is replaced by _matrix_.
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
            [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
            [1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0]
            [0 1 0 1 0 0 0 0 0 0 0 0 0 0 1 1]
            [0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 1]
            [0 0 0 1 0 1 0 0 0 0 0 0 0 0 1 0]
            [0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1]
        """
        from sage.misc.misc import deprecation
        deprecation("The function matrix is replaced by _matrix_.")
        return self._matrix_()


    def _matrix_(self, reverse=False):
        """
        Return the matrix of right multiplication by the element on
        the power basis $1, x, x^2, \ldots, x^{d-1}$ for the field
        extension.  Thus the \emph{rows} of this matrix give the images
        of each of the $x^i$.

        INPUT:
            reverse -- if True act on vectors in reversed order

        EXAMPLE:
            sage: k.<a> = GF(2^4)
            sage: a._vector_(reverse=True), a._matrix_(reverse=True) * a._vector_(reverse=True)
            ((0, 0, 1, 0), (0, 1, 0, 0))
            sage: vector(a), matrix(a) * vector(a)
            ((0, 1, 0, 0), (0, 0, 1, 0))
        """
        import sage.matrix.matrix_space

        K = self.parent()
        a = K.gen()
        x = K(1)
        d = K.degree()

        columns = []

        if not reverse:
            l = xrange(d)
        else:
            l = reversed(range(d))

        for i in l:
            columns.append( (self * x)._vector_() )
            x *= a

        k = K.base_ring()
        M = sage.matrix.matrix_space.MatrixSpace(k, d)

        if reverse:
            return M(columns)
        else:
            return M(columns).transpose()
    def _latex_(self):
        r"""
        Return the latex representation of self, which is just the
        latex representation of the polynomial representation of self.

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: b._latex_()
            'b'
            sage: (b^2+1)._latex_()
            'b + 4'
        """
        if self.parent().degree()>1:
            return self.polynomial()._latex_()
        else:
            return str(self)

    def _pari_init_(self, var=None):
        """
        Return string that when evaluated in PARI defines this element.

        EXAMPLES:
            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b._pari_init_()
            'Mod(b, Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'
            sage: (2*b+3)._pari_init_()
            'Mod(2*b + 3, Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'
        """
        g = self.parent()._finite_field_ext_pari_modulus_as_str()
        f = str(self.polynomial())
        s = 'Mod(%s, %s)'%(f, g)
        if var is None:
            return s
        return s.replace(self.parent().variable_name(), var)

    def charpoly(self, var='x', algorithm='matrix'):
        """
        Return the characteristic polynomial of self as a polynomial with given variable.

        INPUT:
            var -- string (default: 'x')
            algorithm -- string (default: 'matrix')
                         'matrix' -- return the charpoly computed from the
                                     matrix of left multiplication by self
                         'pari' -- use pari's charpoly routine on polymods,
                                   which is not very good except in small cases

        The result is not cached.

        EXAMPLES:
            sage: k.<a> = GF(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: a.charpoly('X')
            X^2 + 18*X + 2
            sage: a^2 + 18*a + 2
            0
            sage: a.charpoly('X', algorithm='pari')
            X^2 + 18*X + 2
        """
        if algorithm == 'matrix':
            return self._matrix_().charpoly(var)
        elif algorithm == 'pari':
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self.parent().prime_subfield(), var)
            return R(self._pari_().charpoly('x').lift())
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm

    def norm(self):
        """
        Return the norm of self down to the prime subfield.

        This is the product of the Galois conjugates of self.

        EXAMPLES:
            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b.norm()
            2
            sage: b.charpoly('t')
            t^2 + 4*t + 2

        Next we consider a cubic extension:
            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: a.norm()
            2
            sage: a.charpoly('t')
            t^3 + 3*t + 3
            sage: a * a^5 * (a^25)
            2
        """
        f = self.charpoly('x')
        n = f[0]
        if f.degree() % 2 != 0:
            return -n
        else:
            return n

    def trace(self):
        """
        Return the trace of this element, which is the sum of the
        Galois conjugates.

        EXAMPLES:
            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: a.trace()
            0
            sage: a.charpoly('t')
            t^3 + 3*t + 3
            sage: a + a^5 + a^25
            0
            sage: z = a^2 + a + 1
            sage: z.trace()
            2
            sage: z.charpoly('t')
            t^3 + 3*t^2 + 2*t + 2
            sage: z + z^5 + z^25
            2
        """
        return self.parent().prime_subfield()(self._pari_().trace().lift())

    def multiplicative_order(self):
        """
        Return the multiplicative order of this field element.

        """
        import sage.rings.arith

        if self.is_zero():
            return ArithmeticError, "Multiplicative order of 0 not defined."
        n = self._parent.order() - 1
        order = 1
        for p, e in sage.rings.arith.factor(n):
            # Determine the power of p that divides the order.
            a = self**(n/(p**e))
            while a != 1:
                order = order * p
                a = a**p

        return order

    def nth_root(self, int n, extend = False, all = False):
        r"""
        Returns an nth root of self.

        INPUT:
            n -- integer >= 1 (must fit in C int type)
            extend -- bool (default: True); if True, return an nth
                 root in an extension ring, if necessary. Otherwise,
                 raise a ValueError if the root is not in the base
                 ring.  Warning: this option is not implemented!
            all -- bool (default: False); if True, return all nth
                 roots of self, instead of just one.

        OUTPUT:
            If self has an nth root, returns one (if all == False) or a list of
            all of them (if all == True).  Otherwise, raises a ValueError (if
            extend = False) or a NotImplementedError (if extend = True).

        WARNING:
            The 'extend' option is not implemented (yet).

        AUTHOR:
            -- David Roe (2007-10-3)

        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.integer import Integer
        if extend:
            raise NotImplementedError
        R = PolynomialRing(self.parent(), "x")
        f = R([-self] + [self.parent()(Integer(0))] * (n - 1) + [self.parent()(1)])
        L = f.roots()
        if all:
            return [x[0] for x in L]
        else:
            if len(L) == 0:
                raise ValueError, "no nth root"
            else:
                return L[0][0]

    def additive_order(self):
        """
        Return the additive order of this finite field element.

        EXAMPLES:
            sage: k.<a> = FiniteField(2^12, 'a')
            sage: b = a^3 + a + 1
            sage: b.additive_order()
            2
            sage: k(0).additive_order()
            1
        """
        if self.is_zero():
            from sage.rings.integer import Integer
            return Integer(1)
        return self.parent().characteristic()

    def pth_power(self, int k = 1):
        """
        Return the $p^k$th power of self, where $p$ is the characteristic
        of the field.

        INPUT:
            k -- integer (default: 1, must fit in C int type)

        OUTPUT:
	    The $p^k$th power of self.

        Note that if $k$ is negative, then this computes the appropriate root.

	EXAMPLES:
            sage: F.<a> = GF(29^2)
	    sage: z = a^2 + 5*a + 1
            sage: z.pth_power()
            19*a + 20
            sage: z.pth_power(10)
            10*a + 28
            sage: z.pth_power(-10) == z
            True
            sage: F.<b> = GF(2^12)
            sage: y = b^3 + b + 1
            sage: y == (y.pth_power(-3))^(2^3)
            True
            sage: y.pth_power(2)
            b^7 + b^6 + b^5 + b^4 + b^3 + b
        """
        p = self.additive_order()
        n = self.parent().degree()
        return self**(p**(k % n))

    frobenius = pth_power

    def pth_root(self, int k = 1):
        """
        Return the $p^k$th root of self, where $p$ is the characteristic
        of the field.

        INPUT:
            k -- integer (default: 1, must fit in C int type)

        OUTPUT:
	    The $p^k$th root of self.

        Note that if $k$ is negative, then this computes the appropriate power.

	EXAMPLES:
            sage: F.<b> = GF(2^12)
            sage: y = b^3 + b + 1
            sage: y == (y.pth_root(3))^(2^3)
            True
            sage: y.pth_root(2)
            b^11 + b^10 + b^9 + b^7 + b^5 + b^4 + b^2 + b
        """
        return self.pth_power(-k)


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
    def __invert__(self):
        from sage.rings.all import ZZ
        return ZZ(0)

cdef class PlusInfinityElement(InfinityElement):
    pass

cdef class MinusInfinityElement(InfinityElement):
    pass

include "coerce.pxi"


#################################################################################
#
#  Coercion of elements
#
#################################################################################

cpdef canonical_coercion(x, y):
    """
    canonical_coercion(x,y) is what is called before doing an
    arithmetic operation between x and y.  It returns a pair (z,w)
    such that z is got from x and w from y via canonical coercion and
    the parents of z and w are identical.

    EXAMPLES:
        sage: A = Matrix([[0,1],[1,0]])
        sage: canonical_coercion(A,1)
        ([0 1]
        [1 0], [1 0]
        [0 1])
    """
    global coercion_model
    return coercion_model.canonical_coercion(x,y)

cpdef bin_op(x, y, op):
    global coercion_model
    return coercion_model.bin_op(x,y,op)


def coerce(Parent p, x):
    try:
        return p._coerce_c(x)
    except AttributeError:
        return p(x)

def coerce_cmp(x,y):
    global coercion_model
    cdef int c
    try:
        x, y = coercion_model.canonical_coercion(x, y)
        return cmp(x,y)
    except TypeError:
        c = cmp(type(x), type(y))
        if c == 0: c = -1
        return c

# We define this base class here to avoid circular cimports.
cdef class CoercionModel:
    """
    Most basic coersion scheme. If it doesn't already match, throw an error.
    """
    cpdef canonical_coercion(self, x, y):
        if parent_c(x) is parent_c(y):
            return x,y
        raise TypeError, "no common canonical parent for objects with parents: '%s' and '%s'"%(parent_c(x), parent_c(y))

    cpdef bin_op(self, x, y, op):
        if parent_c(x) is parent_c(y):
            return op(x,y)
        raise TypeError, arith_error_message(x,y,op)

import coerce
cdef CoercionModel coercion_model = coerce.CoercionModel_cache_maps()

def get_coercion_model():
    """
    Return the global coercion model.

    EXAMPLES:
       sage: import sage.structure.element as e
       sage: cm = e.get_coercion_model()
       sage: cm
       <sage.structure.coerce.CoercionModel_cache_maps object at ...>
    """
    return coercion_model

def set_coercion_model(cm):
    global coercion_model
    coercion_model = cm

def coercion_traceback(dump=True):
    r"""
    This function is very helpful in debugging coercion errors. It prints
    the tracebacks of all the errors caught in the coercion detection. Note
    that failure is cached, so some errors may be omitted the second time
    around (as it remembers not to retry failed paths for speed reasons.

    EXAMPLES:
        sage: 1 + 1/5
        6/5
        sage: coercion_traceback()  # Should be empty, as all went well.
        sage: 1/5 + GF(5).gen()
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for '+': 'Rational Field' and 'Finite Field of size 5'
        sage: coercion_traceback()
        Traceback (most recent call last):
        ...
        TypeError: no common canonical parent for objects with parents: 'Rational Field' and 'Finite Field of size 5'
    """
    if dump:
        for exc_info in coercion_model.exception_stack():
            print ''.join(traceback.format_exception(*exc_info))
    else:
        return coercion_model.exception_stack()

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




######################

def generic_power(a, n, one=None):
    """
    Computes $a^n$, where $n$ is an integer, and $a$ is an object which
    supports multiplication.  Optionally an additional argument,
    which is used in the case that n == 0:

        one: the "unit" element, returned directly (can be anything)

    If this is not supplied, int(1) is returned.

    EXAMPLES:
	sage: from sage.structure.element import generic_power
        sage: generic_power(int(12),int(0))
        1
        sage: generic_power(int(0),int(100))
        0
        sage: generic_power(Integer(10),Integer(0))
        1
        sage: generic_power(Integer(0),Integer(23))
        0
        sage: sum([generic_power(2,i) for i in range(17)]) #test all 4-bit combinations
        131071
        sage: F = Zmod(5)
        sage: a = generic_power(F(2), 5); a
        2
        sage: a.parent() is F
        True
        sage: a = generic_power(F(1), 2)
        sage: a.parent() is F
        True

        sage: generic_power(int(5), 0)
        1
    """

    return generic_power_c(a,n,one)

cdef generic_power_c(a, nn, one):
    try:
        n = PyNumber_Index(nn)
    except TypeError:
        raise NotImplementedError, "non-integral exponents not supported"

    if not a:
        # a is zero, return it, or raise an exception if n is zero
        if not n:
            raise ArithmeticError, "0^0 is undefined."
        elif n < 0:
            raise ZeroDivisionError
        else:
            return a
    elif not n:
        # n is zero, let's try our best to return a one that doesn't suck
        if one is None:
            if PY_TYPE_CHECK(a, Element):
                return (<Element>a)._parent(1)
            try:
                try:
                    return a.parent()(1)
                except AttributeError:
                    return type(a)(1)
            except:
                return 1 #oops, the one sucks
        else:
            return one

    if n < 0:
        # negative exponent! flip it all around.
        a = ~a
        n = -n

    if n < 4:
        # These cases will probably be called often
        # and don't benefit from the code below
        if n == 1:
            return a
        elif n == 2:
            return a*a
        elif n == 3:
            return a*a*a

    # check for idempotence, and store the result otherwise
    aa = a*a
    if aa == a:
        return a

    # since we've computed a^2, let's start squaring there
    # so, let's keep the least-significant bit around, just
    # in case.
    m = n & 1
    n = n >> 1

    # One multiplication can be saved by starting with
    # the second-smallest power needed rather than with 1
    # we've already squared a, so let's start there.
    apow = aa
    while n&1 == 0:
        apow = apow*apow
        n = n >> 1
    power = apow
    n = n >> 1

    # now multiply that least-significant bit in...
    if m:
        power = power * a

    # and this is straight from the book.
    while n != 0:
        apow = apow*apow
        if n&1 != 0:
            power = power*apow
        n = n >> 1

    return power
