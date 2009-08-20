"""
Direct Access to Singular's Functions via libSingular.

This file implements a wrapper around the Singular interpreter which
allows to call any function directly from Sage without string parsing
or interprocess communication overhead.

AUTHORS:

- Michael Brickenstein (2009-07): initial implementation, overall design
- Martin Albrecht (2009-07): clean up, enhancements, etc.

EXAMPLES::

    sage: from sage.libs.singular.function import singular_function
    sage: P.<a,b,c,d> = PolynomialRing(GF(7))
    sage: std = singular_function('std')
    sage: I = sage.rings.ideal.Cyclic(P)
    sage: std(I)
    [a + b + c + d,
     b^2 + 2*b*d + d^2,
     b*c^2 + c^2*d - b*d^2 - d^3,
     b*c*d^2 + c^2*d^2 - b*d^3 + c*d^3 - d^4 - 1,
     b*d^4 + d^5 - b - d,
     c^3*d^2 + c^2*d^3 - c - d,
     c^2*d^4 + b*c - b*d + c*d - 2*d^2]

TESTS::

    sage: from sage.libs.singular.function import singular_function
    sage: std = singular_function('std')
    sage: loads(dumps(std)) == std
    True
"""
#*****************************************************************************
#       Copyright (C) 2009 Michael Brickenstein <brickenstein@mfo.de>
#       Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/stdsage.pxi"
include "../../ext/interrupt.pxi"

from sage.structure.sage_object cimport SageObject

from sage.rings.integer cimport Integer

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal

from sage.rings.polynomial.multi_polynomial_ideal_libsingular cimport sage_ideal_to_singular_ideal, singular_ideal_to_sage_sequence

from sage.libs.singular.decl cimport leftv, idhdl, poly, ideal, ring, number, intvec
from sage.libs.singular.decl cimport sleftv_bin, omAllocBin, omFreeBin, omStrDup
from sage.libs.singular.decl cimport iiMake_proc, iiExprArith1, iiExprArith2, iiExprArith3, iiExprArithM, iiLibCmd
from sage.libs.singular.decl cimport ggetid, IDEAL_CMD, CMD_M, POLY_CMD, PROC_CMD, RING_CMD, NUMBER_CMD, INT_CMD, INTVEC_CMD
from sage.libs.singular.decl cimport MODUL_CMD, LIST_CMD, MATRIX_CMD, VECTOR_CMD, STRING_CMD, V_LOAD_LIB, V_REDEFINE
from sage.libs.singular.decl cimport IsCmd, rChangeCurrRing, currRing, p_Copy
from sage.libs.singular.decl cimport IDROOT, enterid, currRingHdl
from sage.libs.singular.decl cimport errorreported, verbose, Sy_bit

from sage.libs.singular.singular cimport sa2si, si2sa, si2sa_intvec

from sage.interfaces.singular import get_docstring

from sage.misc.misc import get_verbose

cdef leftv* new_leftv(void *data, res_type):
    """
    INPUT:

    - ``data`` - some Singular data this interpreter object points to

    - ``res_type`` - the type of that data
    """
    cdef leftv* res
    res = <leftv*>omAllocBin(sleftv_bin)
    res.Init()
    res.data = data
    res.rtyp = res_type
    return res

cdef free_leftv(leftv * args):
    """
    Kills this ``leftv`` and all ``leftv``s in the tail.

    INPUT:

    - ``args`` - a list of Singular arguments
    """
    args.CleanUp()
    omFreeBin(args, sleftv_bin)

cdef class Converter(SageObject):
    """
    A :class:`Converter` interfaces between Sage objects and Singular
    interpreter objects.
    """
    cdef leftv *args
    cdef MPolynomialRing_libsingular _ring

    def __init__(self, args, ring):
        """
        Create a new argument list.

        INPUT:

        - ``args`` - a list of Python objects

        - ``ring`` - a multivariate polynomial ring

        EXAMPLE::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: Converter([a,b,c],ring=P)
            Singular Converter in Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
        """
        self.args = NULL
        self._ring = ring

        for a in args:
            if PY_TYPE_CHECK(a, MPolynomial_libsingular):
                self.append_polynomial(<MPolynomial_libsingular> a)
            elif PY_TYPE_CHECK(a, MPolynomialIdeal):
                self.append_ideal(a)
            elif PY_TYPE_CHECK(a, int) or PY_TYPE_CHECK(a, long):
                self.append_int(a)
            elif PY_TYPE_CHECK(a, basestring):
                self.append_str(a)
            elif a.parent() is self._ring.base_ring():
                self.append_number(a)
            elif PY_TYPE_CHECK(a, Integer):
                self.append_int(a)
            # ring
            # string
            # vector
            # matrix
            # list
            # modul
            else:
                raise ValueError("unknown argument type '%s'"%(type(a),))

    def ring(self):
        """
        Return the ring in which the arguments of this list live.

        EXAMPLE::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: Converter([a,b,c],ring=P).ring()
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
        """
        return self._ring

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: Converter([a,b,c],ring=P) # indirect doctest
            Singular Converter in Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
        """
        return "Singular Converter in %s"%(self._ring)

    def __dealloc__(self):
        if self.args:
            free_leftv(self.args)

    def __len__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: len(Converter([a,b,c],ring=P))
            3
        """
        cdef leftv * v
        v=self.args
        cdef int l
        l=0
        while v != NULL:
            l=l+1
            v=v.next
        return l

    cdef leftv* pop_front(self) except NULL:
        """
        Pop a Singular element from the front of the list.
        """
        assert(self.args != NULL)
        cdef leftv *res = self.args
        self.args = self.args.next
        res.next = NULL
        return res

    cdef void _append_leftv(self, leftv *v):
        """
        Append a new Singular element to the list.
        """
        cdef leftv* last
        if not self.args == NULL:
            last = self.args
            while not last.next == NULL:
                last=last.next
            last.next=v
        else:
            self.args = v

    cdef void _append(self, void* data, int res_type):
        """
        Create a new ``leftv`` and append it to the list.

        INPUT:

        - ``data`` - the raw data

        - ``res_type`` - the type of the data
        """
        self._append_leftv( new_leftv(data, res_type) )

    cdef void append_polynomial(self, MPolynomial_libsingular p):
        """
        Append the polynomial ``p`` to the list.
        """
        cdef poly* _p = p_Copy(p._poly, <ring*>(<MPolynomialRing_libsingular>p._parent)._ring)
        self._append(_p, POLY_CMD)

    cdef void append_ideal(self,  i):
        """
        Append the ideal ``i`` to the list.
        """
        cdef ideal* singular_ideal = sage_ideal_to_singular_ideal(i)
        self._append(singular_ideal, IDEAL_CMD)

    cdef void append_number(self, n):
        """
        Append the number ``n`` to the list.
        """
        cdef number *_n =  sa2si(n, self._ring._ring)
        self._append(<void *>_n, NUMBER_CMD)

    cdef void append_int(self, n):
        """
        Append the integer ``n`` to the list.
        """
        cdef long _n =  n
        self._append(<void*>_n, INT_CMD)

    cdef void append_str(self, n):
        """
        Append the string ``n`` to the list.
        """
        cdef char *_n = <char *>n
        self._append(omStrDup(_n), STRING_CMD)

    cdef to_python(self, leftv* to_convert):
        """
        Convert the ``leftv`` to a Python object.

        INPUT:

        - ``to_convert`` - a Singular ``leftv``
        """
        cdef MPolynomial_libsingular res_poly
        cdef int rtyp = to_convert.rtyp
        if rtyp == IDEAL_CMD:
            return singular_ideal_to_sage_sequence(<ideal*>to_convert.data, self._ring._ring, self._ring)

        elif rtyp == POLY_CMD:
            res_poly = MPolynomial_libsingular(self._ring)
            res_poly._poly = <poly*>to_convert.data
            return res_poly

        elif rtyp == INT_CMD:
            return <long>to_convert.data

        elif rtyp == NUMBER_CMD:
            return si2sa(<number *>to_convert.data, self._ring._ring, self._ring.base_ring())

        elif rtyp == INTVEC_CMD:
            return si2sa_intvec(<intvec *>to_convert.data)

        elif rtyp == STRING_CMD:
            ret = <char *>to_convert.data
            return ret

        elif rtyp == VECTOR_CMD:
            raise NotImplementedError("rtyp 'VECTOR_CMD' not implement yet.")

        elif rtyp == RING_CMD:
            raise NotImplementedError("rtyp 'RING_CMD' not implement yet.")

        elif rtyp == MATRIX_CMD:
            raise NotImplementedError("rtyp 'MATRIX_CMD' not implement yet.")

        elif rtyp == LIST_CMD:
            raise NotImplementedError("rtyp 'LIST_CMD' not implement yet.")

        elif rtyp == MODUL_CMD:
            raise NotImplementedError("rtyp 'MODUL_CMD' not implement yet.")

        else:
            raise NotImplementedError("rtyp %d not implemented."%(rtyp))

cdef class BaseCallHandler:
    """
    A call handler is an abstraction which hides the details of the
    implementation differences between kernel and library functions.
    """
    cdef leftv* handle_call(self, Converter argument_list) except NULL:
        """
        Actual function call.
        """
        raise NotImplementedError

    cdef bint free_res(self):
        """
        Do we need to free the result object.
        """
        return False

cdef class LibraryCallHandler(BaseCallHandler):
    """
    A call handler is an abstraction which hides the details of the
    implementation differences between kernel and library functions.

    This class implements calling a library function.

    .. note::

        Do not construct this class directly, use
        :func:`singular_function` instead.
    """
    def __init__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import LibraryCallHandler
            sage: LibraryCallHandler()
            <sage.libs.singular.function.LibraryCallHandler object at 0x...>
        """
        super(LibraryCallHandler, self).__init__()

    cdef leftv* handle_call(self, Converter argument_list):
        return iiMake_proc(self.proc_idhdl, NULL, argument_list.args)

    cdef bint free_res(self):
        """
        We do not need to free the result object for library
        functions.
        """
        return False

cdef class KernelCallHandler(BaseCallHandler):
    """
    A call handler is an abstraction which hides the details of the
    implementation differences between kernel and library functions.

    This class implements calling a kernel function.

    .. note::

        Do not construct this class directly, use
        :func:`singular_function` instead.
    """
    def __init__(self, cmd_n, arity):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import KernelCallHandler
            sage: KernelCallHandler(0,0)
            <sage.libs.singular.function.KernelCallHandler object at 0x...>
        """
        super(KernelCallHandler, self).__init__()
        self.cmd_n = cmd_n
        self.arity = arity

    cdef leftv* handle_call(self, Converter argument_list):
        cdef leftv * res
        res = <leftv*> omAllocBin(sleftv_bin)
        res.Init()
        cdef leftv *arg1
        cdef leftv *arg2
        cdef leftv *arg3
        if self.arity != CMD_M:
            number_of_arguments=len(argument_list)

            if number_of_arguments == 1:
                arg1 = argument_list.pop_front()
                iiExprArith1(res, arg1, self.cmd_n)
                free_leftv(arg1)

            elif number_of_arguments == 2:
                arg1 = argument_list.pop_front()
                arg2 = argument_list.pop_front()
                iiExprArith2(res, arg1, self.cmd_n, arg2, self.cmd_n>255)
                free_leftv(arg1)
                free_leftv(arg2)

            elif number_of_arguments == 3:
                arg1 = argument_list.pop_front()
                arg2 = argument_list.pop_front()
                arg3 = argument_list.pop_front()
                iiExprArith3(res, self.cmd_n, arg1, arg2, arg3)
                free_leftv(arg1)
                free_leftv(arg2)
                free_leftv(arg3)
        else:
          iiExprArithM(res, argument_list.args, self.cmd_n)

        return res

    cdef bint free_res(self):
        """
        We need to free the result object for kernel functions.
        """
        return True

cdef class SingularFunction(SageObject):
    """
    The base class for Singular functions either from the kernel or
    from the library.
    """
    def __init__(self, name):
        """
        INPUT:

        - ``name`` - the name of the function

        EXAMPLE::

            sage: from sage.libs.singular.function import SingularFunction
            sage: SingularFunction('foobar')
            foobar (singular function)
        """
        self._name = name

        global currRingHdl
        if currRingHdl == NULL:
            currRingHdl = enterid("my_awesome_sage_ring", 0, RING_CMD, &IDROOT, 1)
            currRingHdl.data.uring.ref += 1

    cdef BaseCallHandler get_call_handler(self):
        """
        Return a call handler which does the actual work.
        """
        raise NotImplementedError

    cdef bint function_exists(self):
        """
        Return ``True`` if the function exists in this interface.
        """
        raise NotImplementedError

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import SingularFunction
            sage: SingularFunction('foobar') # indirect doctest
            foobar (singular function)
        """
        return "%s (singular function)" %(self._name)

    def __call__(self, *args, MPolynomialRing_libsingular ring=None, bint interruptible=True):
        """
        Call this function with the provided arguments ``args`` in the
        ring ``R``.

        INPUT:

        - ``args`` - a list of arguments

        - ``ring`` - a multivariate polynomial ring

        - ``interruptible`` - if ``True`` pressing Ctrl-C during the
                              execution of this function will
                              interrupt the computation (default: ``True``)

        EXAMPLE::

            sage: from sage.libs.singular.function import singular_function
            sage: size = singular_function('size')
            sage: P.<a,b,c> = PolynomialRing(QQ)
            sage: size(a, ring=P)
            1
            sage: size(2r,ring=P)
            1
            sage: size(2, ring=P)
            1
            sage: size(2)
            Traceback (most recent call last):
            ...
            ValueError: Could not detect ring.
            sage: size(Ideal([a*b + c, a + 1]), ring=P)
            2
            sage: size(Ideal([a*b + c, a + 1]))
            2
            sage: size(1,2, ring=P)
            Traceback (most recent call last):
            ...
            RuntimeError

            sage: size('foobar', ring=P)
            6
        """
        ring = self.common_ring(args, ring)
        return call_function(self, args, ring, interruptible)

    def _sage_doc_(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import singular_function
            sage: groebner = singular_function('groebner')
            sage: groebner._sage_doc_()
             '\ngroebner\n--------\n\nProcedure from library `standard.lib\'...
        """
        return get_docstring(self._name)

    cdef MPolynomialRing_libsingular common_ring(self, tuple args, ring=None):
        """
        Return the common ring for the argument list ``args``.

        If ``ring`` is not ``None`` this routine checks whether it is
        the parent/ring of all members of ``args`` instead.

        If no common ring was found a ``ValueError`` is raised.

        INPUT:

        - ``args`` - a list of Python objects

        - ``ring`` - an optional ring to check
        """
        for a in args:
            if PY_TYPE_CHECK(a, MPolynomialIdeal):
                ring2 = a.ring()
            elif PY_TYPE_CHECK(a, MPolynomial_libsingular):
                ring2 = a.parent()
            elif PY_TYPE_CHECK(a, int) or PY_TYPE_CHECK(a, long) or PY_TYPE_CHECK(a, basestring):
                continue
            elif ring is not None:
                a.parent() is ring
                continue

            if ring is None:
                ring = ring2
            elif ring is not ring2:
                raise ValueError("Rings do not match up.")
        if ring is None:
            raise ValueError("Could not detect ring.")
        return <MPolynomialRing_libsingular>ring

    def __reduce__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import singular_function
            sage: groebner = singular_function('groebner')
            sage: groebner == loads(dumps(groebner))
            True
        """
        return singular_function, (self._name,)

    def __cmp__(self, other):
        """
        EXAMPLE::

            sage: from sage.libs.singular.function import singular_function
            sage: groebner = singular_function('groebner')
            sage: groebner == singular_function('groebner')
            True
            sage: groebner == singular_function('std')
            False
            sage: groebner == 1
            False
        """
        if not PY_TYPE_CHECK(other, SingularFunction):
            return cmp(type(self),type(other))
        else:
            return cmp(self._name, (<SingularFunction>other)._name)

cdef inline call_function(SingularFunction self, tuple args, MPolynomialRing_libsingular R, bint signal_handler=True):
    global currRingHdl
    global errorreported

    cdef ring *si_ring = R._ring

    if si_ring != currRing:
        rChangeCurrRing(si_ring)

    if currRingHdl.data.uring!= currRing:
        currRingHdl.data.uring.ref -= 1
        currRingHdl.data.uring = currRing # ref counting?
        currRingHdl.data.uring.ref += 1

    cdef Converter argument_list = Converter(args, R)

    cdef leftv * _res

    if signal_handler:
        _sig_on
        _res = self.call_handler.handle_call(argument_list)
        _sig_off
    else:
        _res = self.call_handler.handle_call(argument_list)

    if errorreported:
        errorreported = 0
        raise RuntimeError

    res = argument_list.to_python(_res)

    if self.call_handler.free_res():
        free_leftv(_res)
    else:
        _res.CleanUp()

    return res

cdef class SingularLibraryFunction(SingularFunction):
    """
    EXAMPLES::

        sage: from sage.libs.singular.function import SingularLibraryFunction
        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: I = R.ideal(x, x+1)
        sage: f = SingularLibraryFunction("groebner")
        sage: f(I)
        [1]
    """
    def __init__(self, name):
        """
        Construct a new Singular kernel function.

        EXAMPLES::

            sage: from sage.libs.singular.function import SingularLibraryFunction
            sage: R.<x,y> = PolynomialRing(QQ, order='lex')
            sage: I = R.ideal(x + 1, x*y + 1)
            sage: f = SingularLibraryFunction("groebner")
            sage: f(I)
            [y - 1, x + 1]
        """
        super(SingularLibraryFunction,self).__init__(name)
        self.call_handler = self.get_call_handler()

    cdef BaseCallHandler get_call_handler(self):
        cdef idhdl* singular_idhdl = ggetid(self._name, 0)
        if singular_idhdl==NULL:
            raise NameError("Function '%s' is not defined."%self._name)
        if singular_idhdl.typ!=PROC_CMD:
            raise ValueError("Not a procedure")

        cdef LibraryCallHandler res = LibraryCallHandler()
        res.proc_idhdl = singular_idhdl
        return res

    cdef bint function_exists(self):
        cdef idhdl* singular_idhdl = ggetid(self._name, 0)
        return singular_idhdl!=NULL

cdef class SingularKernelFunction(SingularFunction):
    """
    EXAMPLES::

        sage: from sage.libs.singular.function import SingularKernelFunction
        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: I = R.ideal(x, x+1)
        sage: f = SingularKernelFunction("std")
        sage: f(I)
        [1]
    """
    def __init__(self, name):
        """
        Construct a new Singular kernel function.

        EXAMPLES::

            sage: from sage.libs.singular.function import SingularKernelFunction
            sage: R.<x,y> = PolynomialRing(QQ, order='lex')
            sage: I = R.ideal(x + 1, x*y + 1)
            sage: f = SingularKernelFunction("std")
            sage: f(I)
            [y - 1, x + 1]
        """
        super(SingularKernelFunction,self).__init__(name)
        self.call_handler = self.get_call_handler()

    cdef BaseCallHandler get_call_handler(self):
        cdef int cmd_n = -1
        arity = IsCmd(self._name, cmd_n) # call by reverence for CMD_n
        if cmd_n == -1:
            raise NameError("Function '%s' is not defined."%self._name)

        return KernelCallHandler(cmd_n, arity)

    cdef bint function_exists(self):
        cdef int cmd_n = -1
        arity = IsCmd(self._name, cmd_n) # call by reverence for CMD_n
        return cmd_n != -1


def singular_function(name):
    """
    Construct a new libSingular function object for the given
    ``name``.

    This function works both for interpreter and built-in functions.

    INPUT:

    - ``name`` - the name of the function

    EXAMPLE::

        sage: from sage.libs.singular.function import singular_function
        sage: P.<x,y,z> = PolynomialRing(QQ)
        sage: f = 3*x*y + 2*z + 1
        sage: g = 2*x + 1/2
        sage: I = Ideal([f,g])

        sage: number_e = singular_function('number_e')
        Traceback (most recent call last):
        ...
        NameError: Function 'number_e' is not defined.

        sage: from sage.libs.singular.function import lib as singular_lib
        sage: singular_lib('general.lib')
        sage: number_e = singular_function('number_e')
        sage: number_e(10r,ring=P)
        67957045707/25000000000
        sage: RR(number_e(10r,ring=P))
        2.71828182828000

        sage: std = singular_function("std")
        sage: std(I)
        [3*y - 8*z - 4, 4*x + 1]
    """

    cdef SingularFunction fnc
    try:
        return SingularKernelFunction(name)
    except NameError:
        return SingularLibraryFunction(name)

def lib(name):
    """
    Load the Singular library ``name``.

    INPUT:

    - ``name`` - a Singular library name

    EXAMPLE::

        sage: from sage.libs.singular.function import singular_function
        sage: from sage.libs.singular.function import lib as singular_lib
        sage: singular_lib('general.lib')
        sage: primes = singular_function('primes')
        sage: primes(2,10, ring=GF(127)['x,y,z'])
        (2, 3, 5, 7)
    """
    global verbose
    cdef int vv = verbose

    if get_verbose() <= 0:
        verbose &= ~Sy_bit(V_LOAD_LIB)

    if get_verbose() <= 0:
        verbose &= ~Sy_bit(V_REDEFINE)

    cdef bint failure = iiLibCmd(omStrDup(name), 1, 1, 1)
    verbose = vv

    if failure:
        raise NameError("Library '%s' not found."%(name,))
