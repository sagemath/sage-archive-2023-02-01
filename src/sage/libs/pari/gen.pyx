"""
PARI C-library interface

AUTHORS:
    -- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta
             (this involved changing almost every doc strings, among other
             things; the precision behavior of PARI seems to change
             from any version to the next...).
    -- William Stein (2006-03-06): added newtonpoly
    -- Justin Walker: contributed some of the function definitions
    -- Gonzalo Tornaria: improvements to conversions; much better error handling.

EXAMPLES:
    sage: pari('5! + 10/x')
    (120*x + 10)/x
    sage: pari('intnum(x=0,13,sin(x)+sin(x^2) + x)')
    85.18856819515268446242866615                # 32-bit
    85.188568195152684462428666150825866897      # 64-bit
    sage: f = pari('x^3-1')
    sage: v = f.factor(); v
    [x - 1, 1; x^2 + x + 1, 1]
    sage: v[0]   # indexing is 0-based unlike in GP.
    [x - 1, x^2 + x + 1]~
    sage: v[1]
    [1, 1]~

Arithmetic obeys the usual coercion rules.
    sage: type(pari(1) + 1)
    <type 'sage.libs.pari.gen.gen'>
    sage: type(1 + pari(1))
    <type 'sage.libs.pari.gen.gen'>
"""

import math
import types
from sage.misc.misc import xsrange
import operator
import sage.structure.element
from sage.structure.element cimport ModuleElement, RingElement, Element
from sage.structure.parent cimport Parent

#include '../../ext/interrupt.pxi'
include 'pari_err.pxi'
include 'setlvalue.pxi'
include '../../ext/stdsage.pxi'

# The unique running Pari instance.
cdef PariInstance pari_instance, P
#pari_instance = PariInstance(200000000, 500000)
pari_instance = PariInstance(100000000, 500000)
#pari_instance = PariInstance(75000000, 500000)
#pari_instance = PariInstance(50000000, 500000)
P = pari_instance   # shorthand notation

# so Galois groups are represented in a sane way
# See the polgalois section of the PARI users manual.
new_galois_format = 1

cdef pari_sp mytop

# keep track of the stack
cdef pari_sp stack_avma

# real precision
cdef unsigned long prec
prec = GP_DATA.fmt.sigd

# Also a copy of PARI accessible from external pure python code.
pari = pari_instance

# temp variables
cdef GEN t0,t1,t2,t3,t4

cdef t0GEN(x):
    global t0
    t0 = P.toGEN(x)

cdef t1GEN(x):
    global t1
    t1 = P.toGEN(x)

cdef t2GEN(x):
    global t2
    t2 = P.toGEN(x)

cdef t3GEN(x):
    global t3
    t3 = P.toGEN(x)

cdef t4GEN(x):
    global t4
    t4 = P.toGEN(x)

cdef class gen(sage.structure.element.RingElement):
    """
    Python extension class that models the PARI GEN type.
    """
    def __init__(self):
        self.b = 0
        self._parent = P

    def parent(self):
        return pari_instance

    cdef void init(self, GEN g, pari_sp b):
        """
            g -- PARI GEN
            b -- pointer to memory chunk where PARI gen lives
                 (if nonzero then this memory is freed when the object
                  goes out of scope)
        """
        self.g = g
        self.b = b
        self._parent = P

    def __dealloc__(self):
        if self.b:
            sage_free(<void*> self.b)

    def __repr__(self):
        return P.GEN_to_str(self.g)

    def __hash__(self):
        return hash(P.GEN_to_str(self.g))

    def _testclass(self):
        import test
        T = test.testclass()
        T._init(self)
        return T

    cdef GEN _gen(self):
        return self.g

    def list(self):
        if typ(self.g) == t_POL:
            raise NotImplementedError, \
                "please report, t_POL.list() is broken, should not be used!"
        if typ(self.g) == t_SER:
            raise NotImplementedError, \
                "please report, t_SER.list() is broken, should not be used!"
        #return list(self.Vecrev())
        return list(self.Vec())

    def __reduce__(self):
        """
        EXAMPLES:
            sage: f = pari('x^3 - 3')
            sage: loads(dumps(f)) == f
            True
        """
        s = str(self)
        import sage.libs.pari.gen_py
        return sage.libs.pari.gen_py.pari, (s,)

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        _sig_on
        return P.new_gen(gadd(self.g, (<gen>right).g))

    def _add_unsafe(gen self, gen right):
        """
        VERY FAST addition of self and right on stack (and leave on
        stack) without any type checking.

        Basically, this is often about 10 times faster than just typing "self + right".
        The drawback is that (1) if self + right would give an error in PARI, it will
        totally crash SAGE, and (2) the memory used by self + right is *never*
        returned -- it gets allocated on the PARI stack and will never be freed.

        EXAMPLES:
            sage: pari(2)._add_unsafe(pari(3))
            5
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gadd(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z,0)
        mytop = avma
        return w

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        _sig_on
        return P.new_gen(gsub(self.g, (<gen> right).g))

    def _sub_unsafe(gen self, gen right):
        """
        VERY FAST subtraction of self and right on stack (and leave on
        stack) without any type checking.

        Basically, this is often about 10 times faster than just typing "self - right".
        The drawback is that (1) if self - right would give an error in PARI, it will
        totally crash SAGE, and (2) the memory used by self + right is *never*
        returned -- it gets allocated on the PARI stack and will never be freed.

        EXAMPLES:
            sage: pari(2)._sub_unsafe(pari(3))
            -1
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gsub(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z, 0)
        mytop = avma
        return w

    cdef RingElement _mul_c_impl(self, RingElement right):
        _sig_on
        return P.new_gen(gmul(self.g, (<gen>right).g))

    def _mul_unsafe(gen self, gen right):
        """
        VERY FAST multiplication of self and right on stack (and leave on
        stack) without any type checking.

        Basically, this is often about 10 times faster than just typing "self * right".
        The drawback is that (1) if self - right would give an error in PARI, it will
        totally crash SAGE, and (2) the memory used by self + right is *never*
        returned -- it gets allocated on the PARI stack and will never be freed.

        EXAMPLES:
            sage: pari(2)._mul_unsafe(pari(3))
            6
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gmul(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z, 0)
        mytop = avma
        return w

    cdef RingElement _div_c_impl(self, RingElement right):
        _sig_on
        return P.new_gen(gdiv(self.g, (<gen>right).g))

    def _div_unsafe(gen self, gen right):
        """
        VERY FAST division of self and right on stack (and leave on
        stack) without any type checking.

        Basically, this is often about 10 times faster than just typing "self / right".
        The drawback is that (1) if self - right would give an error in PARI, it will
        totally crash SAGE, and (2) the memory used by self + right is *never*
        returned -- it gets allocated on the PARI stack and will never be freed.

        EXAMPLES:
            sage: pari(2)._div_unsafe(pari(3))
            2/3
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gdiv(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z, 0)
        mytop = avma
        return w

    #################################################################

    def _mod(gen self, gen other):
        _sig_on
        return P.new_gen(gmod(self.g, other.g))

    def __mod__(self, other):
        if isinstance(self, gen) and isinstance(other, gen):
            return self._mod(other)
        return sage.structure.element.bin_op(self, other, operator.mod)

    def __pow__(self, n, m):
        t0GEN(self)
        t1GEN(n)
        _sig_on
        return P.new_gen(gpow(t0, t1, prec))

    def __neg__(gen self):
        _sig_on
        return P.new_gen(gneg(self.g))

    def __xor__(gen self, n):
        raise RuntimeError, "Use ** for exponentiation, not '^', which means xor\n"+\
              "in Python, and has the wrong precedence."

    def __rshift__(gen self, long n):
        _sig_on
        return P.new_gen(gshift(self.g, -n))

    def __lshift__(gen self, long n):
        _sig_on
        return P.new_gen(gshift(self.g, n))

    def __invert__(gen self):
        _sig_on
        return P.new_gen(ginv(self.g))

    ###########################################
    # ACCESS
    ###########################################
    #def __getattr__(self, attr):
    def getattr(self, attr):
        t0GEN(str(self) + '.' + str(attr))
        return self.new_gen(t0)

    def __getitem__(gen self, n):
        """
        Return a *copy* of the nth entry.

        The indexing is 0-based, like in Python.  However, *copying*
        the nth entry instead of returning reference is different
        than what one usually does in Python.  However, we do it
        this way for consistency with the PARI/GP interpreter, which
        does make copies.

        EXAMPLES:
            sage: p = pari('1 + 2*x + 3*x^2')
            sage: p[0]
            1
            sage: p[2]
            3
            sage: p[100]
            0
            sage: p[-1]
            0
            sage: q = pari('x^2 + 3*x^3 + O(x^6)')
            sage: q[3]
            3
            sage: q[5]
            0
            sage: q[6]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: m = pari('[1,2;3,4]')
            sage: m[0]
            [1, 3]~
            sage: m[1,0]
            3
            sage: l = pari('List([1,2,3])')
            sage: l[1]
            2
            sage: s = pari('"hello, world!"')
            sage: s[0]
            'h'
            sage: s[4]
            'o'
            sage: s[12]
            '!'
            sage: s[13]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: v = pari('[1,2,3]')
            sage: v[0]
            1
            sage: c = pari('Col([1,2,3])')
            sage: c[1]
            2
            sage: sv = pari('Vecsmall([1,2,3])')
            sage: sv[2]
            3
            sage: type(sv[2])
            <type 'int'>
            sage: tuple(pari(3/5))
            (3, 5)
            sage: tuple(pari('1 + 5*I'))
            (1, 5)
            sage: tuple(pari('Qfb(1, 2, 3)'))
            (1, 2, 3)
            sage: pari(57)[0]
            Traceback (most recent call last):
            ...
            TypeError: unindexable object

        """
        if typ(self.g) in (t_INT, t_REAL, t_PADIC, t_QUAD):
            # these are definitely scalar!
            raise TypeError, "unindexable object"

        if typ(self.g) in (t_INTMOD, t_POLMOD):
            # if we keep going we would get:
            #   [0] = modulus
            #   [1] = lift to t_INT or t_POL
            # do we want this? maybe the other way around?
            raise TypeError, "unindexable object"

        #if typ(self.g) in (t_FRAC, t_RFRAC):
            # generic code gives us:
            #   [0] = numerator
            #   [1] = denominator

        #if typ(self.g) == t_COMPLEX:
            # generic code gives us
            #   [0] = real part
            #   [1] = imag part

        #if type(self.g) in (t_QFR, t_QFI):
            # generic code works ok

        if typ(self.g) == t_POL:
            return self.polcoeff(n)

        if typ(self.g) == t_SER:
            bound = valp(self.g) + lg(self.g) - 2
            if n >= bound:
                raise IndexError, "index out of bounds"
            return self.polcoeff(n)

        if isinstance(n, tuple):
            if typ(self.g) != t_MAT:
                raise TypeError, "an integer is required"
            if len(n) != 2:
                raise TypeError, "index must be an integer or a 2-tuple (i,j)"
            i, j = n[0], n[1]
            if i < 0 or i >= self.nrows():
                raise IndexError, "row index out of bounds"
            if j < 0 or j >= self.ncols():
                raise IndexError, "column index out of bounds"
            return P.new_ref(gmael(self.g,j+1,i+1), self)

        if n < 0 or n >= glength(self.g):
            raise IndexError, "index out of bounds"
        if typ(self.g) == t_LIST:
            return P.new_ref(gel(self.g,n+2), self)
        if typ(self.g) == t_STR:
            return chr( (<char *>(self.g+1))[n] )
        if typ(self.g) == t_VECSMALL:
            return self.g[n+1]
        return P.new_ref(gel(self.g,n+1), self)


    def __getslice__(self,  Py_ssize_t i,  Py_ssize_t j):
        """
        EXAMPLES:
            sage: v = pari(xrange(20))
            sage: v[2:5]
            [2, 3, 4]
            sage: v[:]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[10:]
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[:5]
            [0, 1, 2, 3, 4]
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[-1]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: v[:-3]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
            sage: v[5:]
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        """
        cdef long l, k
        l = glength(self.g)
        if j >= l: j = l
        if i < 0: i = 0
        if j <= i: return P.vector(0)
        v = P.vector(j-i)
        for k from i <= k < j:
            v[k-i] = self[k]
        return v

    def gen_length(gen self):
        return lg(self.g)

    def __setitem__(gen self, n, y):
        r"""
        Set the nth entry to a reference to y.

        \begin{notice}
        \begin{itemize}
            \item There is a known BUG: If v is a vector and entry i of v is a vector,
                       \code{v[i][j] = x}
               should set entry j of v[i] to x.  Instead it sets it to nonsense.
               I do not understand why this occurs.  The following is a safe way
               to do the same thing:
                        \code{tmp = v[i]; tmp[j] = x    }

            \item The indexing is 0-based, like everywhere else in Python, but
               {\em unlike} in GP/PARI.

            \item Assignment sets the nth entry to a reference to y,
               assuming y is an object of type gen.  This is the same
               as in Python, but {\em different} than what happens in the
               gp interpreter, where assignment makes a copy of y.

            \item Because setting creates references it is {\em possible} to
               make circular references, unlike in GP.  Do {\em not} do
               this (see the example below).  If you need circular
               references, work at the Python level (where they work
               well), not the PARI object level.
        \end{itemize}
        \end{notice}

        EXAMPLES:
            sage: v = pari(range(10))
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[0] = 10
            sage: w = pari([5,8,-20])
            sage: v
            [10, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[1] = w
            sage: v
            [10, [5, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w[0] = -30
            sage: v
            [10, [-30, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: t = v[1]; t[1] = 10   # don't do v[1][1] !!! (because of mysterious BUG...)
            sage: v
            [10, [-30, 10, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w
            [-30, 10, -20]

        Finally, we create a circular reference:
            sage: v = pari([0])
            sage: w = pari([v])
            sage: v
            [0]
            sage: w
            [[0]]
            sage: v[0] = w
            sage: # Now there is a circular reference.  Trying to access v[0] will crash SAGE.

        """
        cdef gen x
        _sig_on
        x = pari(y)
        if isinstance(n, tuple):
            try:
                self.refers_to[n] = x
            except TypeError:
                self.refers_to = {n:x}
            i = n[0]
            j = n[1]
            if i < 0 or i >= self.nrows():
                raise IndexError, "row i(=%s) must be between 0 and %s"%(i,self.nrows())
            if j < 0 or j >= self.ncols():
                raise IndexError, "column j(=%s) must be between 0 and %s"%(j,self.ncols())
            (<GEN>(self.g)[j+1])[i+1] = <long> x.g
            return

        if n < 0 or n >= glength(self.g):
            raise IndexError, "index (%s) must be between 0 and %s"%(n,glength(self.g)-1)

        # so python memory manager will work correctly
        # and not free x if PARI part of self is the
        # only thing pointing to it.
        try:
            self.refers_to[n] = x
        except TypeError:
            self.refers_to = {n:x}

        if typ(self.g) == t_POL:
            n = n + 1
        (self.g)[n+1] = <long>(x.g)

    def __len__(gen self):
        return glength(self.g)



    ###########################################
    # comparisons
    # I had to put the call to gcmp in another
    # function since otherwise I can't trap
    # the PariError it will sometimes raise.
    # (This might be a bug/shortcoming to SageX.)
    # Annoyingly the _cmp method always has
    # to be not cdef'd.
    ###########################################

    def _cmp(gen self, gen other):
        cdef int result
        _sig_on
        result = gcmp(self.g, other.g)
        _sig_off
        return result

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Comparisons

        First uses PARI's cmp routine; if it decides the objects are
        not comparable, it then compares the underlying strings (since
        in Python all objects are supposed to be comparable).

        EXAMPLES:
            sage: a = pari(5)
            sage: b = 10
            sage: a < b
            True
            sage: a <= b
            True
            sage: a <= 5
            True
            sage: a > b
            False
            sage: a >= b
            False
            sage: a >= pari(10)
            False
            sage: a == 5
            True
            sage: a is 5
            False

            sage: pari(2.5) > None
            False
            sage: pari(3) == pari(3)
            True
            sage: pari('x^2 + 1') == pari('I-1')
            False
            sage: pari(I) == pari(I)
            True
        """
        try:
            return left._cmp(right)
        except PariError:
            pass
        return cmp(str(left),str(right))

    def copy(gen self):
        return P.new_gen(forcecopy(self.g))

    ###########################################
    # Conversion --> Python
    # Try to convert to a meaningful python object
    # in various ways
    ###########################################

    def list_str(gen self):
        """
        Return str that might correctly evaluate to a Python-list.
        """
        s = str(self)
        if s[:4] == "Mat(":
            s = "[" + s[4:-1] + "]"
        s = s.replace("~","")
        if s.find(";") != -1:
            s = s.replace(";","], [")
            s = "[" + s + "]"
            return eval(s)
        else:
            return eval(s)

    def __hex__(gen self):
        """
        Return the hexadecimal digits of self in lower case.

        EXAMPLES:
            sage: print hex(pari(0))
            0
            sage: print hex(pari(15))
            f
            sage: print hex(pari(16))
            10
            sage: print hex(pari(16938402384092843092843098243))
            36bb1e3929d1a8fe2802f083
            sage: print hex(long(16938402384092843092843098243))
            0x36bb1e3929d1a8fe2802f083L
            sage: print hex(pari(-16938402384092843092843098243))
            -36bb1e3929d1a8fe2802f083
        """
        cdef GEN x
        cdef long lx, *xp
        cdef long w
        cdef char *s, *sp
        cdef char *hexdigits
        hexdigits = "0123456789abcdef"
        cdef int i, j
        cdef int size
        x = self.g
        if typ(x) != t_INT:
            raise TypeError, "gen must be of PARI type t_INT"
        if not signe(x):
            return "0"
        lx = lgefint(x)-2  # number of words
        size = lx*2*BYTES_IN_LONG
        s = <char *>sage_malloc(size+2) # 1 char for sign, 1 char for '\0'
        sp = s + size+1
        sp[0] = 0
        xp = int_LSW(x)
        for i from 0 <= i < lx:
            w = xp[0]
            for j from 0 <= j < 2*BYTES_IN_LONG:
                sp = sp-1
                sp[0] = hexdigits[w & 15]
                w = w>>4
            xp = int_nextW(xp)
        # remove leading zeros!
        while sp[0] == c'0':
            sp = sp+1
        if signe(x) < 0:
            sp = sp-1
            sp[0] = c'-'
        k = sp
        sage_free(s)
        return k

    def __int__(gen self):
        """
        Return Python int.  Very fast, and if the number is too large to
        fit into a C int, a Python long is returned instead.

        EXAMPLES:

            sage: int(pari(0))
            0
            sage: int(pari(10))
            10
            sage: int(pari(-10))
            -10
            sage: int(pari(123456789012345678901234567890))
            123456789012345678901234567890L
            sage: int(pari(-123456789012345678901234567890))
            -123456789012345678901234567890L
            sage: int(pari(2^31-1))
            2147483647
            sage: int(pari(-2^31))
            -2147483648
        """
        cdef GEN x
        cdef long lx, *xp
        if  typ(self.g)==t_POL and self.poldegree()<=0:
            # this is a work around for the case that e.g.
            # GF(q).random_element() returns zero, which is of type
            # t_POL(MOD)
            x = polcoeff0(self.g,0,self.get_var(-1))
        else:
            x = self.g
        if typ(x) != t_INT:
            raise TypeError, "gen must be of PARI type t_INT or t_POL of degree 0"
        if not signe(x):
            return 0
        lx = lgefint(x)-3   # take 1 to account for the MSW
        xp = int_MSW(x)
        # special case 1 word so we return int if it fits
        if not lx:
            if   signe(x) |  xp[0] > 0:     # both positive
                return xp[0]
            elif signe(x) & -xp[0] < 0:     # both negative
                return -xp[0]
        i = <ulong>xp[0]
        while lx:
            xp = int_precW(xp)
            i = i << BITS_IN_LONG | <ulong>xp[0]
            lx = lx-1
        if signe(x) > 0:
            return i
        else:
            return -i
        # NOTE: Could use int_unsafe below, which would be much faster, but
        # the default PARI prints annoying stuff to the screen when
        # the number is large.

    def int_unsafe(gen self):
        """
        Returns int form of self, but raises an exception if int does
        not fit into a long integer.

        This is about 5 times faster than the usual int conversion.
        """
        return gtolong(self.g)

    def intvec_unsafe(self):
        """
        Returns Python int list form of entries of self, but raises an
        exception if int does not fit into a long integer.  Here self
        must be a vector.
        """
        cdef int n, L
        cdef object v
        cdef GEN g
        g = self.g
        if typ(g) != t_VEC:
            raise TypeError, "gen must be of PARI type t_VEC"

        L = glength(g)
        v = []
        for n from 0 <= n < L:
            v.append(gtolong(<GEN> (g[n+1])))
        return v

    def python_list_small(gen self):
        """
        Return a Python list of the PARI gens.  This object must be of
        type t_VECSMALL, and the resulting list contains python 'int's

        EXAMPLES:
            sage: v=pari([1,2,3,10,102,10]).Vecsmall()
            sage: w = v.python_list_small()
            sage: w
            [1, 2, 3, 10, 102, 10]
            sage: type(w[0])
            <type 'int'>
        """
        cdef long n
        if typ(self.g) != t_VECSMALL:
            raise TypeError, "Object (=%s) must be of type t_VECSMALL."%self
        V = []
        m = glength(self.g)
        for n from 0 <= n < m:
            V.append(self.g[n+1])
        return V

    def python_list(gen self):
        """
        Return a Python list of the PARI gens.  This object must be of
        type t_VEC

        INPUT: None
        OUTPUT:
           list -- Python list whose elements are the
                   elements of the input gen.
        EXAMPLES:
            sage: v=pari([1,2,3,10,102,10])
            sage: w = v.python_list()
            sage: w
            [1, 2, 3, 10, 102, 10]
            sage: type(w[0])
            <type 'sage.libs.pari.gen.gen'>
        """
        if typ(self.g) != t_VEC:
            raise TypeError, "Object (=%s) must be of type t_VEC."%self
        cdef long n, m
        cdef gen t
        m = glength(self.g)
        V = []
        for n from 0 <= n < m:
            t = P.new_ref(<GEN> (self.g[n+1]), V)
            V.append(t)
        return V

    def python(self, precision=0, bits_prec=None):
        """
        Return Python eval of self.
        """
        if not bits_prec is None:
            precision = int(bits_prec * 3.4) + 1
        import sage.libs.pari.gen_py
        cdef long orig_prec
        if precision:
            orig_prec = P.get_real_precision()
            P.set_real_precision(precision)
            x = sage.libs.pari.gen_py.python(self)
            P.set_real_precision(orig_prec)
            return x
        else:
            return sage.libs.pari.gen_py.python(self)

    def _sage_(self, precision=0):
        """
        Return SAGE version of self.
        """
        import sage.libs.pari.gen_py
        cdef long orig_prec
        if precision:
            orig_prec = P.get_real_precision()
            P.set_real_precision(precision)
            x = sage.libs.pari.gen_py.python(self)
            P.set_real_precision(orig_prec)
            return x
        else:
            return sage.libs.pari.gen_py.python(self)

    def _eval_(gen self):
        return self.python()

    def __long__(gen self):
        """
        Return Python long.
        """
        return long(int(self))

    def __float__(gen self):
        """
        Return Python float.
        """
        cdef double d
        _sig_on
        d = gtodouble(self.g)
        _sig_off
        return d

    def __nonzero__(self):
        """
        EXAMPLES:
            sage: pari('1').__nonzero__()
            True
            sage: pari('x').__nonzero__()
            True
            sage: bool(pari(0))
            False
        """
        return not gcmp0(self.g)


    ###########################################
    # arith1.c
    ###########################################
    def isprime(gen self, flag=0):
        """
        isprime(x, flag=0): Returns True if x is a PROVEN prime
        number, and False otherwise.

        INPUT:
            flag -- int
                    0 (default): use a combination of algorithms.
                    1: certify primality using the Pocklington-Lehmer Test.
                    2: certify primality using the APRCL test.
        OUTPUT:
            bool -- True or False

        EXAMPLES:
            sage: pari(9).isprime()
            False
            sage: pari(17).isprime()
            True
            sage: n = pari(561)    # smallest Carmichael number
            sage: n.isprime()      # not just a pseudo-primality test!
            False
            sage: n.isprime(1)
            False
            sage: n.isprime(2)
            False
        """
        _sig_on
        t = bool(gisprime(self.g, flag) != stoi(0))
        _sig_off
        return t

    def hclassno(gen n):
        """
        Computes the Hurwitz-Kronecker class number of n.

	EXAMPLES:
            sage: pari(-10007).hclassno()
            77
            sage: pari(-3).hclassno()
	    1/3
        """
        _sig_on
        return P.new_gen(hclassno(n.g))

    def ispseudoprime(gen self, flag=0):
        """
        ispseudoprime(x, flag=0): Returns True if x is a pseudo-prime
        number, and False otherwise.

        INPUT:
            flag -- int
                    0 (default): checks whether x is a Baillie-Pomerance-Selfridge-Wagstaff pseudo prime (strong Rabin-Miller pseudo prime for base 2, followed by strong Lucas test for the sequence (P,-1), P smallest positive integer such that P^2 - 4 is not a square mod x).
                    > 0: checks whether x is a strong Miller-Rabin pseudo prime for flag randomly chosen bases (with end-matching to catch square roots of -1).

        OUTPUT:
            bool -- True or False

        EXAMPLES:
            sage: pari(9).ispseudoprime()
            False
            sage: pari(17).ispseudoprime()
            True
            sage: n = pari(561)     # smallest Carmichael number
            sage: n.ispseudoprime() # not just any old pseudo-primality test!
            False
            sage: n.ispseudoprime(2)
            False
        """
        cdef GEN z
        _sig_on
        z = gispseudoprime(self.g, flag)
        _sig_off
        return bool(gcmp(z, stoi(0)))

    def ispower(gen self, k=None):
        r"""
        Determine whether or not self is a perfect k-th power.
        If k is not specified, find the largest k so that self
        is a k-th power.

        NOTE: There is a BUG in the PARI C-library function (at least
        in PARI 2.2.12-beta) that is used to implement this function!
        This is in GP:
        \begin{verbatim}
           ? p=nextprime(10^100); n=p^100; m=p^2; m^50==n; ispower(n,50)
        \end{verbatim}

        INPUT:
            k -- int (optional)

        OUTPUT:
            power -- int, what power it is
            g -- what it is a power of

        EXAMPLES:
            sage: pari(9).ispower()
            (2, 3)
            sage: pari(17).ispower()
            (1, 17)
            sage: pari(17).ispower(2)
            (False, None)
            sage: pari(17).ispower(1)
            (1, 17)
            sage: pari(2).ispower()
            (1, 2)
        """
        cdef int n
        cdef GEN x

        _sig_on
        if k is None:
            _sig_on
            n = isanypower(self.g, &x)
            if n == 0:
                _sig_off
                return 1, self
            else:
                return n, P.new_gen(x)
        else:
            k = int(k)
            t0GEN(k)
            _sig_on
            n = ispower(self.g, t0, &x)
            if n == 0:
                _sig_off
                return False, None
            else:
                return k, P.new_gen(x)

    ###########################################
    # 1: Standard monadic or dyadic OPERATORS
    ###########################################
    def divrem(gen x, y, var=-1):
        """
        divrem(x, y, {v}): Euclidean division of x by y giving as a
            2-dimensional column vector the quotient and the
            remainder, with respect to v (to main variable if v is
            omitted).
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(divrem(x.g, t0, P.get_var(var)))

    def lex(gen x, y):
        """

        lex(x,y): Compare x and y lexicographically (1 if x>y, 0 if
            x==y, -1 if x<y)

        """
        t0GEN(y)
        _sig_on
        return lexcmp(x.g, t0)

    def max(gen x, y):
        """
        max(x,y): Return the maximum of x and y.
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gmax(x.g, t0))

    def min(gen x, y):
        """
        min(x,y): Return the minumum of x and y.
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gmin(x.g, t0))

    def shift(gen x, long n):
        """
        shift(x,n): shift x left n bits if n>=0, right -n bits if n<0.
        """
        return P.new_gen(gshift(x.g, n))

    def shiftmul(gen x, long n):
        """
        shiftmul(x,n): Return the product of x by $2^n$.
        """
        return P.new_gen(gmul2n(x.g, n))

    def moebius(gen x):
        """
        moebius(x): Moebius function of x.
        """
        return P.new_gen(gmu(x.g))

    def sign(gen x):
        """
        sign(x): Return the sign of x, where x is of type integer, real
            or fraction.
        """
        # Pari throws an error if you attempt to take the sign of
        # a complex number.
        _sig_on
        return gsigne(x.g)

    def vecmax(gen x):
        """
        vecmax(x): Return the maximum of the elements of the vector/matrix x,
        """
        _sig_on
        return P.new_gen(vecmax(x.g))


    def vecmin(gen x):
        """
        vecmin(x): Return the maximum of the elements of the vector/matrix x,
        """
        _sig_on
        return P.new_gen(vecmin(x.g))



    ###########################################
    # 2: CONVERSIONS and similar elementary functions
    ###########################################


    def Col(gen x):
        """
        Col(x): Transforms the object x into a column vector.

        The vector will have only one component, except in the
        following cases:

           * When x is a vector or a quadratic form, the resulting
             vector is the initial object considered as a column
             vector.

           * When x is a matrix, the resulting vector is the column of
             row vectors comprising the matrix.

           * When x is a character string, the result is a column of
             individual characters.

           * When x is a polynomial, the coefficients of the vector
             start with the leading coefficient of the polynomial.

           * When x is a power series, only the significant
             coefficients are taken into account, but this time by
             increasing order of degree.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('1.5').Col()
            [1.500000000000000000000000000]~               # 32-bit
            [1.5000000000000000000000000000000000000]~     # 64-bit
            sage: pari([1,2,3,4]).Col()
            [1, 2, 3, 4]~
            sage: pari('[1,2; 3,4]').Col()
            [[1, 2], [3, 4]]~
            sage: pari('"SAGE"').Col()
            ["S", "A", "G", "E"]~
            sage: pari('3*x^3 + x').Col()
            [3, 0, 1, 0]~
            sage: pari('x + 3*x^3 + O(x^5)').Col()
            [1, 0, 3, 0]~
        """
        _sig_on
        return P.new_gen(gtocol(x.g))

    def List(gen x):
        """
        List(x): transforms the PARI vector or list x into a list.

        EXAMPLES:
            sage: v = pari([1,2,3])
            sage: v
            [1, 2, 3]
            sage: v.type()
            't_VEC'
            sage: w = v.List()
            sage: w
            List([1, 2, 3])
            sage: w.type()
            't_LIST'
        """
        _sig_on
        return P.new_gen(gtolist(x.g))

    def Mat(gen x):
        """
        Mat(x): Returns the matrix defined by x.

           * If x is already a matrix, a copy of x is created and
             returned.

           * If x is not a vector or a matrix, this function returns a
             1x1 matrix.

           * If x is a row (resp. column) vector, this functions
             returns a 1-row (resp. 1-column) matrix, *unless* all
             elements are column (resp. row) vectors of the same
             length, in which case the vectors are concatenated
             sideways and the associated big matrix is returned.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- a PARI matrix

        EXAMPLES:
            sage: x = pari(5)
            sage: x.type()
            't_INT'
            sage: y = x.Mat()
            sage: y
            Mat(5)
            sage: y.type()
            't_MAT'
            sage: x = pari('[1,2;3,4]')
            sage: x.type()
            't_MAT'
            sage: x = pari('[1,2,3,4]')
            sage: x.type()
            't_VEC'
            sage: y = x.Mat()
            sage: y
            Mat([1, 2, 3, 4])
            sage: y.type()
            't_MAT'

            sage: v = pari('[1,2;3,4]').Vec(); v
            [[1, 3]~, [2, 4]~]
            sage: v.Mat()
            [1, 2; 3, 4]
            sage: v = pari('[1,2;3,4]').Col(); v
            [[1, 2], [3, 4]]~
            sage: v.Mat()
            [1, 2; 3, 4]
        """
        _sig_on
        return P.new_gen(gtomat(x.g))

    def Mod(gen x, y):
        """
        Mod(x, y): Returns the object x modulo y, denoted Mod(x, y).

        The input y must be a an integer or a polynomial:

           * If y is an INTEGER, x must also be an integer, a rational
             number, or a p-adic number compatible with the modulus y.

           * If y is a POLYNOMIAL, x must be a scalar (which is not a polmod),
             a polynomial, a rational function, or a power series.

        WARNING: This function is not the same as x % y, the result of
        which is an integer or a polynomial.

        INPUT:
            x -- gen
            y -- integer or polynomial

        OUTPUT:
            gen -- intmod or polmod

        EXAMPLES:
            sage: z = pari(3)
            sage: x = z.Mod(pari(7))
            sage: x
            Mod(3, 7)
            sage: x^2
            Mod(2, 7)
            sage: x^100
            Mod(4, 7)
            sage: x.type()
            't_INTMOD'

            sage: f = pari("x^2 + x + 1")
            sage: g = pari("x")
            sage: a = g.Mod(f)
            sage: a
            Mod(x, x^2 + x + 1)
            sage: a*a
            Mod(-x - 1, x^2 + x + 1)
            sage: a.type()
            't_POLMOD'
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gmodulcp(x.g,t0))

    def Pol(self, v=-1):
        """
        Pol(x, {v}): convert x into a polynomial with main variable v
        and return the result.

           * If x is a scalar, returns a constant polynomial.

           * If x is a power series, the effect is identical
              to \kbd{truncate}, i.e.~it chops off the $O(X^k)$.

           * If x is a vector, this function creates the polynomial
             whose coefficients are given in x, with x[0]
             being the leading coefficient (which can be zero).

        WARNING: This is *not* a substitution function. It will not
        transform an object containing variables of higher priority
        than v:
            sage: pari('x+y').Pol('y')
            Traceback (most recent call last):
            ...
            PariError:  (8)

        INPUT:
            x -- gen
            v -- (optional) which variable, defaults to 'x'
        OUTPUT:
            gen -- a polynomial
        EXAMPLES:
            sage: v = pari("[1,2,3,4]")
            sage: f = v.Pol()
            sage: f
            x^3 + 2*x^2 + 3*x + 4
            sage: f*f
            x^6 + 4*x^5 + 10*x^4 + 20*x^3 + 25*x^2 + 24*x + 16

            sage: v = pari("[1,2;3,4]")
            sage: v.Pol()
            [1, 3]~*x + [2, 4]~
        """
        _sig_on
        return P.new_gen(gtopoly(self.g, P.get_var(v)))

    def Polrev(self, v=-1):
        """
        Polrev(x, {v}): Convert x into a polynomial with main variable
        v and return the result.  This is the reverse of Pol if x is a
        vector, otherwise it is identical to Pol.   By "reverse" we mean
        that the coefficients are reversed.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- a polynomial
        EXAMPLES:
            sage: v = pari("[1,2,3,4]")
            sage: f = v.Polrev()
            sage: f
            4*x^3 + 3*x^2 + 2*x + 1
            sage: v.Pol()
            x^3 + 2*x^2 + 3*x + 4
            sage: v.Polrev('y')
            4*y^3 + 3*y^2 + 2*y + 1

        Note that Polrev does *not* reverse the coefficients of a polynomial!
            sage: f
            4*x^3 + 3*x^2 + 2*x + 1
            sage: f.Polrev()
            4*x^3 + 3*x^2 + 2*x + 1
            sage: v = pari("[1,2;3,4]")
            sage: v.Polrev()
            [2, 4]~*x + [1, 3]~
        """
        _sig_on
        return P.new_gen(gtopolyrev(self.g, P.get_var(v)))

    def Qfb(gen a, b, c, D=0):
        """
        Qfb(a,b,c,{D=0.}): Returns the binary quadratic form
        $$
                   ax^2 + bxy + cy^2.
        $$
        The optional D is 0 by default and initializes Shanks's
        distance if $b^2 - 4ac > 0$.

        NOTE: Negative definite forms are not implemented, so use their
        positive definitine counterparts instead.  (I.e., if f is a
        negative definite quadratic form, then -f is positive
        definite.)

        INPUT:
            a -- gen
            b -- gen
            c -- gen
            D -- gen (optional, defaults to 0)
        OUTPUT:
            gen -- binary quadratic form
        EXAMPLES:
            sage: pari(3).Qfb(7, 2)
            Qfb(3, 7, 2, 0.E-250)               # 32-bit
            Qfb(3, 7, 2, 0.E-693)               # 64-bit
        """
        t0GEN(b); t1GEN(c); t2GEN(D)
        _sig_on
        return P.new_gen(Qfb0(a.g, t0, t1, t2, prec))


    def Ser(gen x, v=-1):
        """
        Ser(x,{v=x}): Create a power series from x with main variable v
        and return the result.

           * If x is a scalar, this gives a constant power series with
             precision given by the default series precision, as returned
             by get_series_precision().

           * If x is a polynomial, the precision is the greatest of
             get_series_precision() and the degree of the polynomial.

           * If x is a vector, the precision is similarly given, and
             the coefficients of the vector are understood to be the
             coefficients of the power series starting from the
             constant term (i.e.~the reverse of the function Pol).

        WARNING: This is *not* a substitution function. It will not
        transform an object containing variables of higher priority than v.

        INPUT:
            x -- gen
            v -- PARI variable (default: x)
        OUTPUT:
            gen -- PARI object of PARI type t_SER
        EXAMPLES:
            sage: pari(2).Ser()
            2 + O(x^16)
            sage: x = pari([1,2,3,4,5])
            sage: x.Ser()
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + O(x^5)
            sage: f = x.Ser('v'); print f
            1 + 2*v + 3*v^2 + 4*v^3 + 5*v^4 + O(v^5)
            sage: pari(1)/f
            1 - 2*v + v^2 + O(v^5)
            sage: pari(1).Ser()
            1 + O(x^16)
        """
        _sig_on
        return P.new_gen(gtoser(x.g, P.get_var(v)))


    def Set(gen x):
        """
        Set(x): convert x into a set, i.e. a row vector of strings in
        increasing lexicographic order.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- a vector of strings in increasing lexicographic order.
        EXAMPLES:
            sage: pari([1,5,2]).Set()
            ["1", "2", "5"]
            sage: pari([]).Set()     # the empty set
            []
            sage: pari([1,1,-1,-1,3,3]).Set()
            ["-1", "1", "3"]
            sage: pari(1).Set()
            ["1"]
            sage: pari('1/(x*y)').Set()
            ["1/(y*x)"]
            sage: pari('["bc","ab","bc"]').Set()
            ["ab", "bc"]
        """
        _sig_on
        return P.new_gen(gtoset(x.g))


    def Str(self):
        """
        Str(self): Return the print representation of self as a PARI
        object is returned.

        INPUT:
            self -- gen
        OUTPUT:
            gen -- a PARI gen of type t_STR, i.e., a PARI string
        EXAMPLES:
            sage: pari([1,2,['abc',1]]).Str()
            [1, 2, [abc, 1]]
            sage: pari('[1,1, 1.54]').Str()
            [1, 1, 1.540000000000000000000000000]        # 32-bit
            [1, 1, 1.5400000000000000000000000000000000000]        # 64-bit
            sage: pari(1).Str()       # 1 is automatically converted to string rep
            1
            sage: x = pari('x')       # PARI variable "x"
            sage: x.Str()             # is converted to string rep.
            x
            sage: x.Str().type()
            't_STR'
        """
        cdef char* c
        _sig_on
        c = GENtostr(self.g)
        v = self.new_gen(strtoGENstr(c))
        free(c)
        return v


    def Strchr(gen x):
        """
        Strchr(x): converts x to a string, translating each integer
        into a character (in ASCII).

        NOTE: Vecsmall is (essentially) the inverse to Strchr().

        INPUT:
            x -- PARI vector of integers
        OUTPUT:
            gen -- a PARI string
        EXAMPLES:
            sage: pari([65,66,123]).Strchr()
            AB{
            sage: pari('"SAGE"').Vecsmall()   # pari('"SAGE"') --> PARI t_STR
            Vecsmall([83, 65, 71, 69])
            sage: _.Strchr()
            SAGE
            sage: pari([83, 65, 71, 69]).Strchr()
            SAGE
        """
        _sig_on
        return P.new_gen(Strchr(x.g))

    def Strexpand(gen x):
        """
        Strexpand(x): Concatenate the entries of the vector x into a
        single string, performing tilde expansion.

        NOTE: I have no clue what the point of this function is. -- William
        """
        if x.type() != 't_VEC':
            raise TypeError, "x must be of type t_VEC."
        _sig_on
        return P.new_gen(Strexpand(x.g))


    def Strtex(gen x):
        r"""
        Strtex(x): Translates the vector x of PARI gens to TeX format
        and returns the resulting concatenated strings as a PARI t_STR.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- PARI t_STR (string)
        EXAMPLES:
            sage: v=pari('x^2')
            sage: v.Strtex()
            x^2
            sage: v=pari(['1/x^2','x'])
            sage: v.Strtex()
            \frac{1}{x^2}x
            sage: v=pari(['1 + 1/x + 1/(y+1)','x-1'])
            sage: v.Strtex()
            \frac{ \left(y
             + 2\right)  x
             + \left(y
             + 1\right) }{ \left(y
             + 1\right)  x}x
             - 1
        """
        if x.type() != 't_VEC':
            x = P.vector(1, [x])
        _sig_on
        return P.new_gen(Strtex(x.g))

    def printtex(gen x):
        return x.Strtex()

    def Vec(gen x):
        """
        Vec(x): Transforms the object x into a vector.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- of PARI type t_VEC
        EXAMPLES:
            sage: pari(1).Vec()
            [1]
            sage: pari('x^3').Vec()
            [1, 0, 0, 0]
            sage: pari('x^3 + 3*x - 2').Vec()
            [1, 0, 3, -2]
            sage: pari([1,2,3]).Vec()
            [1, 2, 3]
            sage: pari('ab').Vec()
            [1, 0]
        """
        _sig_on
        return P.new_gen(gtovec(x.g))

    def Vecrev(gen x):
        """
        Vecrev(x): Transforms the object x into a vector.
        Identical to Vec(x) except when x is
        -- a polynomial, this is the reverse of Vec.
        -- a power series, this includes low-order zero coefficients.
        -- a laurant series, raises an exception

        INPUT:
            x -- gen
        OUTPUT:
            gen -- of PARI type t_VEC
        EXAMPLES:
            sage: pari(1).Vecrev()
            [1]
            sage: pari('x^3').Vecrev()
            [0, 0, 0, 1]
            sage: pari('x^3 + 3*x - 2').Vecrev()
            [-2, 3, 0, 1]
            sage: pari([1, 2, 3]).Vecrev()
            [1, 2, 3]
            sage: pari('Col([1, 2, 3])').Vecrev()
            [1, 2, 3]
            sage: pari('[1, 2; 3, 4]').Vecrev()
            [[1, 3]~, [2, 4]~]
            sage: pari('ab').Vecrev()
            [0, 1]
            sage: pari('x^2 + 3*x^3 + O(x^5)').Vecrev()
            [0, 0, 1, 3, 0]
            sage: pari('x^-2 + 3*x^3 + O(x^5)').Vecrev()
            Traceback (most recent call last):
            ...
            ValueError: Vecrev() is not defined for Laurent series
        """
        cdef long lx, vx, i
        cdef GEN y
        if typ(x.g) == t_POL:
            lx = lg(x.g)
            y = cgetg(lx-1, t_VEC)
            for i from 1 <= i <= lx-2:
                # no need to copy, since new_gen will deep copy
                __set_lvalue__(gel(y,i), gel(x.g,i+1))
            return P.new_gen(y)
        elif typ(x.g) == t_SER:
            lx = lg(x.g)
            vx = valp(x.g)
            if vx < 0:
                raise ValueError, "Vecrev() is not defined for Laurent series"
            y = cgetg(vx+lx-1, t_VEC)
            for i from 1 <= i <= vx:
                __set_lvalue__(gel(y,i), gen_0)
            for i from 1 <= i <= lx-2:
                # no need to copy, since new_gen will deep copy
                __set_lvalue__(gel(y,vx+i), gel(x.g,i+1))
            return P.new_gen(y)
        else:
            return x.Vec()

    def Vecsmall(gen x):
        """
        Vecsmall(x): transforms the object x into a t_VECSMALL.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- PARI t_VECSMALL
        EXAMPLES:
            sage: pari([1,2,3]).Vecsmall()
            Vecsmall([1, 2, 3])
            sage: pari('"SAGE"').Vecsmall()
            Vecsmall([83, 65, 71, 69])
            sage: pari(1234).Vecsmall()
            Vecsmall([1234])
        """
        _sig_on
        return P.new_gen(gtovecsmall(x.g))

    def binary(gen x):
        """
        binary(x): gives the vector formed by the binary digits of
        abs(x), where x is of type t_INT.

        INPUT:
            x -- gen of type t_INT
        OUTPUT:
            gen -- of type t_VEC
        EXAMPLES:
            sage: pari(0).binary()
            [0]
            sage: pari(-5).binary()
            [1, 0, 1]
            sage: pari(5).binary()
            [1, 0, 1]
            sage: pari(2005).binary()
            [1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]

            sage: pari('"2"').binary()
            Traceback (most recent call last):
            ...
            TypeError: x (=2) must be of type t_INT, but is of type t_STR.
        """
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        _sig_on
        return P.new_gen(binaire(x.g))

    def bitand(gen x, y):
        """
        bitand(x,y): Bitwise and of two integers x and y. Negative
        numbers behave as if modulo some large power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(8).bitand(4)
            0
            sage: pari(8).bitand(8)
            8
            sage: pari(6).binary()
            [1, 1, 0]
            sage: pari(7).binary()
            [1, 1, 1]
            sage: pari(6).bitand(7)
            6
            sage: pari(19).bitand(-1)
            19
            sage: pari(-1).bitand(-1)
            -1
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gbitand(x.g, t0))


    def bitneg(gen x, long n=-1):
        r"""
        bitneg(x,{n=-1}): Bitwise negation of the integer x truncated
        to n bits.  n=-1 (the default) represents an infinite sequence
        of the bit 1.  Negative numbers behave as if modulo some large
        power of 2.

        With n=-1, this function returns -n-1.  With n >= 0, it returns
        a number a such that  $a\cong -n-1 \pmod{2^n}$.

        INPUT:
            x -- gen (t_INT)
            n -- long, default = -1
        OUTPUT:
            gen -- t_INT
        EXAMPLES:
            sage: pari(10).bitneg()
            -11
            sage: pari(1).bitneg()
            -2
            sage: pari(-2).bitneg()
            1
            sage: pari(-1).bitneg()
            0
            sage: pari(569).bitneg()
            -570
            sage: pari(569).bitneg(10)
            454
            sage: 454 % 2^10
            454
            sage: -570 % 2^10
            454
        """
        _sig_on
        return P.new_gen(gbitneg(x.g,n))


    def bitnegimply(gen x, y):
        """
        bitnegimply(x,y): Bitwise negated imply of two integers x and
        y, in other words, x BITAND BITNEG(y). Negative numbers behave
        as if modulo big power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(14).bitnegimply(0)
            14
            sage: pari(8).bitnegimply(8)
            0
            sage: pari(8+4).bitnegimply(8)
            4
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gbitnegimply(x.g, t0))


    def bitor(gen x, y):
        """
        bitor(x,y): Bitwise or of two integers x and y. Negative
        numbers behave as if modulo big power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(14).bitor(0)
            14
            sage: pari(8).bitor(4)
            12
            sage: pari(12).bitor(1)
            13
            sage: pari(13).bitor(1)
            13
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gbitor(x.g, t0))


    def bittest(gen x, long n):
        """
        bittest(x, long n): Returns bit number n (coefficient of $2^n$ in binary)
        of the integer x. Negative numbers behave as if modulo a big power of 2.

        INPUT:
           x -- gen (pari integer)
        OUTPUT:
           bool -- a Python bool
        EXAMPLES:
            sage: x = pari(6)
            sage: x.bittest(0)
            False
            sage: x.bittest(1)
            True
            sage: x.bittest(2)
            True
            sage: x.bittest(3)
            False
            sage: pari(-3).bittest(0)
            True
            sage: pari(-3).bittest(1)
            False
            sage: [pari(-3).bittest(n) for n in range(10)]
            [True, False, True, True, True, True, True, True, True, True]
        """
        _sig_on
        b = bool(bittest(x.g, n))
        _sig_off
        return b

    def bitxor(gen x, y):
        """
        bitxor(x,y): Bitwise exclusive or of two integers x and y.
        Negative numbers behave as if modulo big power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(6).bitxor(4)
            2
            sage: pari(0).bitxor(4)
            4
            sage: pari(6).bitxor(0)
            6
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gbitxor(x.g, t0))


    def ceil(gen x):
        """
        Return the smallest integer >= x.

        INPUT:
           x -- gen
        OUTPUT:
           gen -- integer
        EXAMPLES:
            sage: pari(1.4).ceil()
            2
            sage: pari(-1.4).ceil()
            -1
            sage: pari('x').ceil()
            x
            sage: pari('x^2+5*x+2.2').ceil()
            x^2 + 5*x + 2.200000000000000000000000000             # 32-bit
            x^2 + 5*x + 2.2000000000000000000000000000000000000   # 64-bit
            sage: pari('3/4').ceil()
            1
        """
        _sig_on
        return P.new_gen(gceil(x.g))

    def centerlift(gen x, v=-1):
        """
        centerlift(x,{v}): Centered lift of x.  This function returns
        exactly the same thing as lift, except if x is an integer mod.

        INPUT:
            x -- gen
            v -- var (default: x)
        OUTPUT:
            gen
        EXAMPLES:
            sage: x = pari(-2).Mod(5)
            sage: x.centerlift()
            -2
            sage: x.lift()
            3
            sage: f = pari('x-1').Mod('x^2 + 1')
            sage: f.centerlift()
            x - 1
            sage: f.lift()
            x - 1
            sage: f = pari('x-y').Mod('x^2+1')
            sage: f
            Mod(x - y, x^2 + 1)
            sage: f.centerlift('x')
            x - y
            sage: f.centerlift('y')
            Mod(x - y, x^2 + 1)
        """
        _sig_on
        return P.new_gen(centerlift0(x.g, P.get_var(v)))


    def changevar(gen x, y):
        """
        changevar(gen x, y): change variables of x according to the vector y.

        WARNING: This doesn't seem to work right at all in SAGE (!).
        Use with caution.   *STRANGE*

        INPUT:
            x -- gen
            y -- gen (or coercible to gen)
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('x^3+1').changevar(pari(['y']))
            y^3 + 1
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(changevar(x.g, t0))

    def component(gen x, long n):
        """
        component(x, long n): Return n'th component of the internal
        representation of x.  This this function is 1-based
        instead of 0-based.

        NOTE: For vectors or matrices, it is simpler to use x[n-1]. For
        list objects such as is output by nfinit, it is easier to use
        member functions.

        INPUT:
            x -- gen
            n -- C long (coercible to)
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari([0,1,2,3,4]).component(1)
            0
            sage: pari([0,1,2,3,4]).component(2)
            1
            sage: pari([0,1,2,3,4]).component(4)
            3
            sage: pari('x^3 + 2').component(1)
            2
            sage: pari('x^3 + 2').component(2)
            0
            sage: pari('x^3 + 2').component(4)
            1

            sage: pari('x').component(0)
            Traceback (most recent call last):
            ...
            PariError:  (8)
        """
        _sig_on
        return P.new_gen(compo(x.g, n))

    def conj(gen x):
        """
        conj(x): Return the algebraic conjugate of x.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('x+1').conj()
            x + 1
            sage: pari('x+I').conj()
            x - I
            sage: pari('1/(2*x+3*I)').conj()
            1/(2*x - 3*I)
            sage: pari([1,2,'2-I','Mod(x,x^2+1)', 'Mod(x,x^2-2)']).conj()
            [1, 2, 2 + I, Mod(-x, x^2 + 1), Mod(-x, x^2 - 2)]
            sage: pari('Mod(x,x^2-2)').conj()
            Mod(-x, x^2 - 2)
            sage: pari('Mod(x,x^3-3)').conj()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (20)
        """
        _sig_on
        return P.new_gen(gconj(x.g))

    def conjvec(gen x):
        """
        conjvec(x): Returns the vector of all conjugates of the
        algebraic number x.  An algebraic number is a polynomial over
        Q modulo an irreducible polynomial.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('Mod(1+x,x^2-2)').conjvec()
            [-0.4142135623730950488016887242, 2.414213562373095048801688724]~                         # 32-bit
            [-0.41421356237309504880168872420969807857, 2.4142135623730950488016887242096980786]~     # 64-bit
            sage: pari('Mod(x,x^3-3)').conjvec()
            [1.442249570307408382321638311, -0.7211247851537041911608191554 + 1.249024766483406479413179544*I, -0.7211247851537041911608191554 - 1.249024766483406479413179544*I]~           # 32-bit
            [1.4422495703074083823216383107801095884, -0.72112478515370419116081915539005479419 + 1.2490247664834064794131795437632536350*I, -0.72112478515370419116081915539005479419 - 1.2490247664834064794131795437632536350*I]~       # 64-bit
        """
        _sig_on
        return P.new_gen(conjvec(x.g, prec))

    def denominator(gen x):
        """
        denominator(x): Return the denominator of x.  When x is a
        vector, this is the least common multiple of the denominators
        of the components of x.

        what about poly?
        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('5/9').denominator()
            9
            sage: pari('(x+1)/(x-2)').denominator()
            x - 2
            sage: pari('2/3 + 5/8*x + 7/3*x^2 + 1/5*y').denominator()
            1
            sage: pari('2/3*x').denominator()
            1
            sage: pari('[2/3, 5/8, 7/3, 1/5]').denominator()
            120
        """
        _sig_on
        return P.new_gen(denom(x.g))

    def floor(gen x):
        """
        floor(x): Return the floor of x, which is the largest integer <= x.
        This function also works component-wise on polynomials, vectors, etc.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('5/9').floor()
            0
            sage: pari('11/9').floor()
            1
            sage: pari('1.17').floor()
            1
            sage: pari('x').floor()
            x
            sage: pari('x+1.5').floor()
            x + 1.500000000000000000000000000             # 32-bit
            x + 1.5000000000000000000000000000000000000   # 64-bit
            sage: pari('[1.5,2.3,4.99]').floor()
            [1, 2, 4]
            sage: pari('[[1.1,2.2],[3.3,4.4]]').floor()
            [[1, 2], [3, 4]]

            sage: pari('"hello world"').floor()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (20)
        """
        _sig_on
        return P.new_gen(gfloor(x.g))

    def frac(gen x):
        """
        frac(x): Return the fractional part of x, which is x - floor(x).

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('1.7').frac()
            0.7000000000000000000000000000            # 32-bit
            0.70000000000000000000000000000000000000  # 64-bit
            sage: pari('sqrt(2)').frac()
            0.4142135623730950488016887242            # 32-bit
            0.41421356237309504880168872420969807857  # 64-bit
            sage: pari('sqrt(-2)').frac()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (20)
        """
        _sig_on
        return P.new_gen(gfrac(x.g))

    def imag(gen x):
        """
        imag(x): Return the imaginary part of x.  This function also
        works component-wise.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('1+2*I').imag()
            2
            sage: pari('sqrt(-2)').imag()
            1.414213562373095048801688724              # 32-bit
            1.4142135623730950488016887242096980786    # 64-bit
            sage: pari('x+I').imag()
            1
            sage: pari('x+2*I').imag()
            2
            sage: pari('(1+I)*x^2+2*I').imag()
            x^2 + 2
            sage: pari('[1,2,3] + [4*I,5,6]').imag()
            [4, 0, 0]
        """
        _sig_on
        return P.new_gen(gimag(x.g))

    def length(gen x):
        """
        length(x): Return the number of non-code words in x.  If x
        is a string, this is the number of characters of x.

        ?? terminator ?? carriage return ??
        """
        return glength(x.g)

    def lift(gen x, v=-1):
        """
        lift(x,{v}): Returns the lift of an element of Z/nZ to Z or
        R[x]/(P) to R[x] for a type R if v is omitted.  If v is given,
        lift only polymods with main variable v.  If v does not occur
        in x, lift only intmods.

        INPUT:
            x -- gen
            v -- (optional) variable
        OUTPUT:
            gen
        EXAMPLES:
            sage: x = pari("x")
            sage: a = x.Mod('x^3 + 17*x + 3')
            sage: a
            Mod(x, x^3 + 17*x + 3)
            sage: b = a^4; b
            Mod(-17*x^2 - 3*x, x^3 + 17*x + 3)
            sage: b.lift()
            -17*x^2 - 3*x

        ??? more examples
        """
        if v == -1:
            _sig_on
            return P.new_gen(lift(x.g))
        _sig_on
        return P.new_gen(lift0(x.g, P.get_var(v)))

    def numbpart(gen x):
        """
        numbpart(x): returns the number of partitions of x.

        EXAMPLES:
            sage: pari(20).numbpart()
            627
            sage: pari(100).numbpart()
            190569292
        """
        _sig_on
        return P.new_gen(numbpart(x.g))

    def numerator(gen x):
        """
        numerator(x): Returns the numerator of x.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:

        """
        return P.new_gen(numer(x.g))


    def numtoperm(gen k, long n):
        """
        numtoperm(k, n): Return the permutation number k (mod n!) of n
        letters, where n is an integer.

        INPUT:
            k -- gen, integer
            n -- int
        OUTPUT:
            gen -- vector (permutation of {1,...,n})
        EXAMPLES:
        """
        _sig_on
        return P.new_gen(numtoperm(n, k.g))


    def padicprec(gen x, p):
        """
        padicprec(x,p): Return the absolute p-adic precision of the object x.

        INPUT:
            x -- gen
        OUTPUT:
            int
        EXAMPLES:
        """
        cdef gen _p
        _p = pari(p)
        if typ(_p.g) != t_INT:
            raise TypeError, "p (=%s) must be of type t_INT, but is of type %s."%(
                _p, _p.type())
        return padicprec(x.g, _p.g)

    def permtonum(gen x):
        """
        permtonum(x): Return the ordinal (between 1 and n!) of permutation vector x.
        ??? Huh ???  say more.  what is a perm vector.  0 to n-1 or 1-n.

        INPUT:
            x -- gen (vector of integers)
        OUTPUT:
            gen -- integer
        EXAMPLES:
        """
        if typ(x.g) != t_VEC:
            raise TypeError, "x (=%s) must be of type t_VEC, but is of type %s."%(x,x.type())
        return P.new_gen(permtonum(x.g))

    def precision(gen x, long n=-1):
        """
        precision(x,{n}): Change the precision of x to be n, where n
        is a C-integer). If n is omitted, output the real precision of x.

        INPUT:
            x -- gen
            n -- (optional) int
        OUTPUT:
            nothing
          or
            gen if n is omitted
        EXAMPLES:
        """
        if n <= -1:
            return precision(x.g)
        return P.new_gen(precision0(x.g, n))

    def random(gen N):
        r"""
        \code{random(\{N=$2^31$\})}: Return a pseudo-random integer between 0 and $N-1$.

        INPUT:
            N -- gen, integer
        OUTPUT:
            gen -- integer
        EXAMPLES:
        """
        if typ(N.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(N,N.type())
        _sig_on
        return P.new_gen(genrand(N.g))

    def real(gen x):
        """
        real(x): Return the real part of x.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
        """
        _sig_on
        return P.new_gen(greal(x.g))

    def round(gen x, estimate=False):
        """
        round(x,estimat=False):  If x is a real number, returns x rounded
        to the nearest integer (rounding up).  If the optional argument
        estimate is True, also returns the binary exponent e of the difference
        between the original and the rounded value (the "fractional part")
        (this is the integer ceiling of log_2(error)).

        When x is a general PARI object, this function returns the result
        of rounding every coefficient at every level of PARI object.
        Note that this is different than what the truncate function
        does (see the example below).

        One use of round is to get exact results after a long
        approximate computation, when theory tells you that the
        coefficients must be integers.

        INPUT:
            x -- gen
            estimate -- (optional) bool, False by default
        OUTPUT:
            * if estimate == False, return a single gen.
            * if estimate == True, return rounded verison of x and
              error estimate in bits, both as gens.
        EXAMPLES:
            sage: pari('1.5').round()
            2
            sage: pari('1.5').round(True)
            (2, -1)
            sage: pari('1.5 + 2.1*I').round()
            2 + 2*I
            sage: pari('1.0001').round(True)
            (1, -14)
            sage: pari('(2.4*x^2 - 1.7)/x').round()
            (2*x^2 - 2)/x
            sage: pari('(2.4*x^2 - 1.7)/x').truncate()
            2.400000000000000000000000000*x                # 32-bit
            2.4000000000000000000000000000000000000*x      # 64-bit
        """
        cdef int n
        if not estimate:
            _sig_on
            return P.new_gen(ground(x.g))
        cdef long e
        cdef gen y
        _sig_on
        y = P.new_gen(grndtoi(x.g, &e))
        return y, e

    def simplify(gen x):
        """
        simplify(x): Simplify the object x as much as possible, and return
        the result.

        A complex or quadratic number whose imaginary part is an exact 0
        (i.e., not an approximate one such as O(3) or 0.E-28) is converted
        to its real part, and a a polynomial of degree 0 is converted to
        its constant term.  Simplification occurs recursively.

        This function is useful before using arithmetic functions, which
        expect integer arguments:

        EXAMPLES:
            sage: y = pari('y')
            sage: x = pari('9') + y - y
            sage: x
            9
            sage: x.type()
            't_POL'
            sage: x.factor()
            matrix(0,2)
            sage: pari('9').factor()
            Mat([3, 2])
            sage: x.simplify()
            9
            sage: x.simplify().factor()
            Mat([3, 2])
            sage: x = pari('1.5 + 0*I')
            sage: x.type()
            't_COMPLEX'
            sage: x.simplify()
            1.500000000000000000000000000                  # 32-bit
            1.5000000000000000000000000000000000000        # 64-bit
            sage: y = x.simplify()
            sage: y.type()
            't_REAL'
        """
        _sig_on
        return P.new_gen(simplify(x.g))

    def sizebyte(gen x):
        """
        sizebyte(x): Return the total number of bytes occupied by the
        complete tree of the object x.  Note that this number depends
        on whether the computer is 32-bit or 64-bit (see examples).

        INPUT:
            x -- gen
        OUTPUT:
            int (a Python int)
        EXAMPLES:
            sage: pari('1').sizebyte()
            12           # 32-bit
            24           # 64-bit
            sage: pari('10').sizebyte()
            12           # 32-bit
            24           # 64-bit
            sage: pari('10000000000000').sizebyte()
            16           # 32-bit
            24           # 64-bit
            sage: pari('10^100').sizebyte()
            52           # 32-bit
            64           # 64-bit
            sage: pari('x').sizebyte()
            36           # 32-bit
            72           # 64-bit
            sage: pari('x^20').sizebyte()
            264          # 32-bit
            528          # 64-bit
            sage: pari('[x, 10^100]').sizebyte()
            100          # 32-bit
            160          # 64-bit
        """
        return taille2(x.g)

    def sizedigit(gen x):
        """

        sizedigit(x): Return a quick estimate for the maximal number of
        decimal digits before the decimal point of any component of x.

        INPUT:
            x -- gen
        OUTPUT:
            int -- Python integer
        EXAMPLES:
            sage: x = pari('10^100')
            sage: x.Str().length()
            101
            sage: x.sizedigit()
            101

        Note that digits after the decimal point are ignored.
            sage: x = pari('1.234')
            sage: x
            1.234000000000000000000000000              # 32-bit
            1.2340000000000000000000000000000000000    # 64-bit
            sage: x.sizedigit()
            1

        The estimate can be one too big:
            sage: pari('7234.1').sizedigit()
            4
            sage: pari('9234.1').sizedigit()
            5
        """
        return sizedigit(x.g)

    def truncate(gen x, estimate=False):
        """
        truncate(x,estimate=False):  Return the truncation of x.
        If estimate is True, also return the number of error bits.

        When x is in the real numbers, this means that the part
        after the decimal point is chopped away, e is the binary
        exponent of the difference between the original and truncated
        value (the "fractional part").    If x is a rational
        function, the result is the integer part (Euclidean
        quotient of numerator by denominator) and if requested
        the error estimate is 0.

        When truncate is applied to a power series (in X), it
        transforms it into a polynomial or a rational function with
        denominator a power of X, by chopping away the $O(X^k)$.
        Similarly, when applied to a p-adic number, it transforms it
        into an integer or a rational number by chopping away the
        $O(p^k)$.

        INPUT:
            x -- gen
            estimate -- (optional) bool, which is False by default
        OUTPUT:
            * if estimate == False, return a single gen.
            * if estimate == True, return rounded verison of x and
              error estimate in bits, both as gens.
        OUTPUT:
        EXAMPLES:
            sage: pari('(x^2+1)/x').round()
            (x^2 + 1)/x
            sage: pari('(x^2+1)/x').truncate()
            x
            sage: pari('1.043').truncate()
            1
            sage: pari('1.043').truncate(True)
            (1, -5)
            sage: pari('1.6').truncate()
            1
            sage: pari('1.6').round()
            2
            sage: pari('1/3 + 2 + 3^2 + O(3^3)').truncate()
            34/3
            sage: pari('sin(x+O(x^10))').truncate()
            1/362880*x^9 - 1/5040*x^7 + 1/120*x^5 - 1/6*x^3 + x
            sage: pari('sin(x+O(x^10))').round()   # each coefficient has abs < 1
            x + O(x^10)
        """
        if not estimate:
            _sig_on
            return P.new_gen(gtrunc(x.g))
        cdef long e
        cdef gen y
        _sig_on
        y = P.new_gen(gcvtoi(x.g, &e))
        return y, e

    def valuation(gen x, p):
        """
        valuation(x,p): Return the valuation of x with respect to p.

        The valuation is the highest exponent of p dividing x.

           * If p is an integer, x must be an integer, an intmod whose
             modulus is divisible by p, a rational number, a p-adic
             number, or a polynomial or power series in which case the
             valuation is the minimal of the valuations of the
             coefficients.

           * If p is a polynomial, x must be a polynomial or a
             rational fucntion.  If p is a monomial then x may also be
             a power series.

           * If x is a vector, complex or quadratic number, then the
             valuation is the minimum of the component valuations.

           * If x = 0, the result is $2^31-1$ on 32-bit machines or
             $2^63-1$ on 64-bit machines if x is an exact object.
             If x is a p-adic number or power series, the result
             is the exponent of the zero.

        INPUT:
            x -- gen
            p -- coercible to gen
        OUTPUT:
            gen -- integer
        EXAMPLES:
            sage: pari(9).valuation(3)
            2
            sage: pari(9).valuation(9)
            1
            sage: x = pari(9).Mod(27); x.valuation(3)
            2
            sage: pari('5/3').valuation(3)
            -1
            sage: pari('9 + 3*x + 15*x^2').valuation(3)
            1
            sage: pari([9,3,15]).valuation(3)
            1
            sage: pari('9 + 3*x + 15*x^2 + O(x^5)').valuation(3)
            1

            sage: pari('x^2*(x+1)^3').valuation(pari('x+1'))
            3
            sage: pari('x + O(x^5)').valuation('x')
            1
            sage: pari('2*x^2 + O(x^5)').valuation('x')
            2

            sage: pari(0).valuation(3)
            2147483647            # 32-bit
            9223372036854775807   # 64-bit
        """
        t0GEN(p)
        _sig_on
        v = ggval(x.g, t0)
        _sig_off
        return v

    def variable(gen x):
        """
        variable(x): Return the main variable of the object x, or p
        if x is a p-adic number.

        This function raises a TypeError exception on scalars, i.e.,
        on objects with no variable associated to them.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('x^2 + x -2').variable()
            x
            sage: pari('1+2^3 + O(2^5)').variable()
            2
            sage: pari('x+y0').variable()
            x
            sage: pari('y0+z0').variable()
            y0
        """
        _sig_on
        return P.new_gen(gpolvar(x.g))


    ###########################################
    # 3: TRANSCENDENTAL functions
    # AUTHORS: Pyrex Code, docs -- Justin Walker (justin@mac.com)
    #          Examples, docs   -- William Stein
    ###########################################

    def abs(gen self):
        """
        Returns the absolute value of x (its modulus, if x is complex).
        Rational functions are not allowed.  Contrary to most transcendental
        functions, an exact argument is not converted to a real number before
        applying abs and an exact result is returned if possible.

        EXAMPLES:
            sage: x = pari("-27.1")
            sage: x.abs()
            27.10000000000000000000000000               # 32-bit
            27.100000000000000000000000000000000000     # 64-bit

        If x is a polynomial, returns -x if the leading coefficient is real
        and negative else returns x.  For a power series, the constant
        coefficient is considered instead.

        EXAMPLES:
            sage: pari('x-1.2*x^2').abs()
            1.200000000000000000000000000*x^2 - x              # 32-bit
            1.2000000000000000000000000000000000000*x^2 - x    # 64-bit
        """
        _sig_on
        return P.new_gen(gabs(self.g, prec))

    def acos(gen x):
        r"""
        The principal branch of $\cos^{-1}(x)$, so that $\Re(\acos(x))$
        belongs to $[0,Pi]$. If $x$ is real and $|x| > 1$,
        then $\acos(x)$ is complex.

        EXAMPLES:
            sage: pari('0.5').acos()
            1.047197551196597746154214461               # 32-bit
            1.0471975511965977461542144610931676281     # 64-bit
            sage: pari('1.1').acos()
            -0.4435682543851151891329110664*I             # 32-bit
            -0.44356825438511518913291106635249808665*I   # 64-bit
            sage: pari('1.1+I').acos()
            0.8493430542452523259630143655 - 1.097709866825328614547942343*I                         # 32-bit
            0.84934305424525232596301436546298780187 - 1.0977098668253286145479423425784723444*I     # 64-bit
        """
        _sig_on
        return P.new_gen(gacos(x.g, prec))

    def acosh(gen x):
        r"""
        The principal branch of $\cosh^{-1}(x)$, so that
        $\Im(\acosh(x))$ belongs to $[0,Pi]$. If $x$ is real and $x <
        1$, then $\acosh(x)$ is complex.

        EXAMPLES:
            sage: pari(2).acosh()
            1.316957896924816708625046347              # 32-bit
            1.3169578969248167086250463473079684440    # 64-bit
            sage: pari(0).acosh()
            1.570796326794896619231321692*I            # 32-bit
            1.5707963267948966192313216916397514421*I  # 64-bit
            sage: pari('I').acosh()
            0.8813735870195430252326093250 + 1.570796326794896619231321692*I    # 32-bit
            0.88137358701954302523260932497979230902 + 1.5707963267948966192313216916397514421*I   # 64-bit
        """
        _sig_on
        return P.new_gen(gach(x.g, prec))

    def agm(gen x, y):
        r"""
        The arithmetic-geometric mean of x and y.  In the case of complex
        or negative numbers, the principal square root is always chosen.
        p-adic or power series arguments are also allowed.  Note that a p-adic
        AGM exists only if x/y is congruent to 1 modulo p (modulo 16 for p=2).
        x and y cannot both be vectors or matrices.

        EXAMPLES:
            sage: pari('2').agm(2)
            2.000000000000000000000000000             # 32-bit
            2.0000000000000000000000000000000000000   # 64-bit
            sage: pari('0').agm(1)
            0
            sage: pari('1').agm(2)
            1.456791031046906869186432383             # 32-bit
            1.4567910310469068691864323832650819750   # 64-bit
            sage: pari('1+I').agm(-3)
            -0.9647317222908759112270275374 + 1.157002829526317260939086020*I                        # 32-bit
            -0.96473172229087591122702753739366831917 + 1.1570028295263172609390860195427517825*I    # 64-bit
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(agm(x.g, t0, prec))

    def arg(gen x):
        r"""
        arg(x): argument of x,such that $-\pi < \arg(x) \leq \pi$.

        EXAMPLES:
            sage: pari('2+I').arg()
            0.4636476090008061162142562315              # 32-bit
	    0.46364760900080611621425623146121440203    # 64-bit
        """
        _sig_on
        return P.new_gen(garg(x.g,prec))

    def asin(gen x):
        r"""
        The principal branch of $\sin^{-1}(x)$, so that
        $\Re(\asin(x))$ belongs to $[-\pi/2,\pi/2]$. If $x$ is real
        and $|x| > 1$ then $\asin(x)$ is complex.

        EXAMPLES:
            sage: pari(pari('0.5').sin()).asin()
            0.5000000000000000000000000000               # 32-bit
            0.50000000000000000000000000000000000000     # 64-bit
            sage: pari(2).asin()
            1.570796326794896619231321692 + 1.316957896924816708625046347*I                      # 32-bit
            1.5707963267948966192313216916397514421 + 1.3169578969248167086250463473079684440*I  # 64-bit
        """
        _sig_on
        return P.new_gen(gasin(x.g, prec))

    def asinh(gen x):
        r"""
        The principal branch of $\sinh^{-1}(x)$, so that $\Im(\asinh(x))$ belongs
        to $[-\pi/2,\pi/2]$.

        EXAMPLES:
            sage: pari(2).asinh()
            1.443635475178810342493276740                # 32-bit
            1.4436354751788103424932767402731052694      # 64-bit
            sage: pari('2+I').asinh()
            1.528570919480998161272456185 + 0.4270785863924761254806468833*I      # 32-bit
            1.5285709194809981612724561847936733933 + 0.42707858639247612548064688331895685930*I      # 64-bit
        """
        _sig_on
        return P.new_gen(gash(x.g, prec))

    def atan(gen x):
        r"""
        The principal branch of $\tan^{-1}(x)$, so that $\Re(\atan(x))$ belongs
        to $]-\pi/2, \pi/2[$.

        EXAMPLES:
            sage: pari(1).atan()
            0.7853981633974483096156608458              # 32-bit
            0.78539816339744830961566084581987572104    # 64-bit
            sage: pari('1.5+I').atan()
            1.107148717794090503017065460 + 0.2554128118829953416027570482*I                         # 32-bit
            1.1071487177940905030170654601785370401 + 0.25541281188299534160275704815183096744*I     # 64-bit
        """
        _sig_on
        return P.new_gen(gatan(x.g, prec))

    def atanh(gen x):
        r"""
        The principal branch of $\tanh^{-1}(x)$, so that
        $\Im(\atanh(x))$ belongs to $]-\pi/2,\pi/2]$.  If $x$ is real
        and $|x| > 1$ then $\atanh(x)$ is complex.

        EXAMPLES:
            sage: pari(0).atanh()
            0.E-250   # 32-bit
            0.E-693   # 64-bit
            sage: pari(2).atanh()
            0.5493061443340548456976226185 + 1.570796326794896619231321692*I           # 32-bit
            0.54930614433405484569762261846126285232 + 1.5707963267948966192313216916397514421*I     # 64-bit
        """
        _sig_on
        return P.new_gen(gath(x.g, prec))

    def bernfrac(gen x):
        r"""
        The Bernoulli number $B_x$, where $B_0 = 1$, $B_1 = -1/2$,
        $B_2 = 1/6,\ldots,$ expressed as a rational number. The
        argument $x$ should be of type integer.

        EXAMPLES:
            sage: pari(18).bernfrac()
            43867/798
            sage: [pari(n).bernfrac() for n in range(10)]
            [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0]
        """
        _sig_on
        return P.new_gen(bernfrac(x))

    def fibonacci(gen x):
        r"""
        Return the fibonacci number of index x.

        EXAMPLES:
            sage: pari(18).bernfrac()
            43867/798
            sage: [pari(n).bernfrac() for n in range(10)]
            [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0]
        """
        _sig_on
        return P.new_gen(fibo(long(x)))

    def bernreal(gen x):
        r"""
        The Bernoulli number $B_x$, as for the function bernfrac, but
        $B_x$ is returned as a real number (with the current
        precision).

        EXAMPLES:
            sage: pari(18).bernreal()
            54.97117794486215538847117794                  # 32-bit
            54.971177944862155388471177944862155388        # 64-bit
        """
        _sig_on
        return P.new_gen(bernreal(x, prec))

    def bernvec(gen x):
        r"""
        Creates a vector containing, as rational numbers, the
        Bernoulli numbers $B_0, B_2,\ldots, B_{2x}$.  This routine is
        obsolete.  Use bernfrac instead each time you need a Bernoulli
        number in exact form.

        Note:  this  routine  is  implemented  using  repeated  independent
        calls to bernfrac, which is faster than the standard recursion in
        exact arithmetic.

        EXAMPLES:
            sage: pari(8).bernvec()
            [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
            sage: [pari(2*n).bernfrac() for n in range(9)]
            [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
        """
        _sig_on
        return P.new_gen(bernvec(x))

    def besselh1(gen nu, x):
        r"""
        The $H^1$-Bessel function of index $\nu$ and argument $x$.

        EXAMPLES:
            sage: pari(2).besselh1(3)
            0.4860912605858910769078310941 - 0.1604003934849237296757682995*I                        # 32-bit
            0.48609126058589107690783109411498403480 - 0.16040039348492372967576829953798091810*I    # 64-bit
        """
        t0GEN(x)
        _sig_on
        return P.new_gen(hbessel1(nu.g, t0, prec))

    def besselh2(gen nu, x):
        r"""
        The $H^2$-Bessel function of index $\nu$ and argument $x$.

        EXAMPLES:
            sage: pari(2).besselh2(3)
            0.4860912605858910769078310941 + 0.1604003934849237296757682995*I                         # 32-bit
            0.48609126058589107690783109411498403480 + 0.16040039348492372967576829953798091810*I     # 64-bit
        """
        t0GEN(x)
        _sig_on
        return P.new_gen(hbessel2(nu.g, t0, prec))

    def besselj(gen nu, x):
        r"""
        Bessel J function (Bessel function of the first kind), with
        index $\nu$ and argument $x$.  If $x$ converts to a power
        series, the initial factor $(x/2)^{\nu}/\Gamma(\nu+1)$ is
        omitted (since it cannot be represented in PARI when $\nu$ is not
        integral).

        EXAMPLES:
            sage: pari(2).besselj(3)
            0.4860912605858910769078310941            # 32-bit
            0.48609126058589107690783109411498403480  # 64-bit
        """
        t0GEN(x)
        _sig_on
        return P.new_gen(jbessel(nu.g, t0, prec))

    def besseljh(gen nu, x):
        """
        J-Bessel function of half integral index (Speherical Bessel
        function of the first kind).  More precisely, besseljh(n,x)
        computes $J_{n+1/2}(x)$ where n must an integer, and x is any
        complex value.  In the current implementation (PARI, version
        2.2.11), this function is not very accurate when $x$ is small.

        EXAMPLES:
            sage: pari(2).besseljh(3)
            0.4127100322097159934374967959      # 32-bit
            0.41271003220971599343749679594186271499    # 64-bit
        """
        t0GEN(x)
        _sig_on
        return P.new_gen(jbesselh(nu.g, t0, prec))

    def besseli(gen nu, x):
        """
        Bessel I function (Bessel function of the second kind), with
        index $\nu$ and argument $x$.  If $x$ converts to a power
        series, the initial factor $(x/2)^{\nu}/\Gamma(\nu+1)$ is
        omitted (since it cannot be represented in PARI when $\nu$ is not
        integral).

        EXAMPLES:
            sage: pari(2).besseli(3)
            2.245212440929951154625478386              # 32-bit
            2.2452124409299511546254783856342650577    # 64-bit
            sage: pari(2).besseli('3+I')
            1.125394076139128134398383103 + 2.083138226706609118835787255*I      # 32-bit
            1.1253940761391281343983831028963896470 + 2.0831382267066091188357872547036161842*I    # 64-bit
        """
        t0GEN(x)
        _sig_on
        return P.new_gen(ibessel(nu.g, t0, prec))

    def besselk(gen nu, x, long flag=0):
        """
        nu.besselk(x, flag=0): K-Bessel function (modified Bessel
        function of the second kind) of index nu, which can be
        complex, and argument x.

        INPUT:
            nu -- a complex number
            x -- real number (positive or negative)
            flag -- default: 0 or 1: use hyperu  (hyperu is much slower for
                    small x, and doesn't work for negative x).

        WARNING/TODO -- with flag = 1 this function is incredibly slow
        (on 64-bit Linux) as it is implemented using it from the C
        library, but it shouldn't be (e.g., it's not slow via the GP
        interface.)  Why?

        EXAMPLES:
            sage: pari('2+I').besselk(3)
            0.04559077184075505871203211094 + 0.02891929465820812820828883526*I     # 32-bit
            0.045590771840755058712032110938791854704 + 0.028919294658208128208288835257608789842*I     # 64-bit

            sage: pari('2+I').besselk(-3)
            -4.348708749867516799575863067 - 5.387448826971091267230878827*I        # 32-bit
            -4.3487087498675167995758630674661864255 - 5.3874488269710912672308788273655523057*I  # 64-bit

            sage.: pari('2+I').besselk(300, flag=1)
            3.742246033197275082909500148 E-132 + 2.490710626415252262644383350 E-134*I      # 32-bit
            3.7422460331972750829095001475885825717 E-132 + 2.4907106264152522626443833495225745762 E-134*I   # 64-bit

        """
        t0GEN(x)
        _sig_on
        return P.new_gen(kbessel0(nu.g, t0, flag, prec))

    def besseln(gen nu, x):
        """
        nu.besseln(x): Bessel N function (Spherical Bessel function of
        the second kind) of index nu and argument x.

        EXAMPLES:
            sage: pari('2+I').besseln(3)
            -0.2807755669582439141676487005 - 0.4867085332237257928816949747*I     # 32-bit
            -0.28077556695824391416764870046044688871 - 0.48670853322372579288169497466916637395*I    # 64-bit
        """
        t0GEN(x)
        _sig_on
        return P.new_gen(nbessel(nu.g, t0, prec))

    def cos(gen self):
        """
        The cosine function.

        EXAMPLES:
            sage: x = pari('1.5')
            sage: x.cos()
            0.07073720166770291008818985142     # 32-bit
	    0.070737201667702910088189851434268709084    # 64-bit
            sage: pari('1+I').cos()
            0.8337300251311490488838853943 - 0.9888977057628650963821295409*I   # 32-bit
            0.83373002513114904888388539433509447980 - 0.98889770576286509638212954089268618864*I   # 64-bit
            sage: pari('x+O(x^8)').cos()
            1 - 1/2*x^2 + 1/24*x^4 - 1/720*x^6 + 1/40320*x^8 + O(x^9)
        """
        _sig_on
        return P.new_gen(gcos(self.g, prec))

    def cosh(gen self):
        """
        The hyperbolic cosine function.

        EXAMPLES:
            sage: x = pari('1.5')
            sage: x.cosh()
            2.352409615243247325767667965               # 32-bit
            2.3524096152432473257676679654416441702     # 64-bit
            sage: pari('1+I').cosh()
            0.8337300251311490488838853943 + 0.9888977057628650963821295409*I                       # 32-bit
            0.83373002513114904888388539433509447980 + 0.98889770576286509638212954089268618864*I   # 64-bit
            sage: pari('x+O(x^8)').cosh()
            1 + 1/2*x^2 + 1/24*x^4 + 1/720*x^6 + O(x^8)
        """
        _sig_on
        return P.new_gen(gch(self.g, prec))

    def cotan(gen x):
        """
        The cotangent of x.

        EXAMPLES:
            sage: pari(5).cotan()
            -0.2958129155327455404277671681     # 32-bit
            -0.29581291553274554042776716808248528607    # 64-bit

        On a 32-bit computer computing the cotangent of $\pi$ doesn't
        raise an error, but instead just returns a very large number.
        On a 64-bit computer it raises a RuntimeError.

            sage: pari('Pi').cotan()
            1.980704062 E28                      # 32-bit
	    Traceback (most recent call last):   # 64-bit
            ...	                                 # 64-bit
            PariError: division by zero (46)     # 64-bit
        """
        _sig_on
        return P.new_gen(gcotan(x.g, prec))

    def dilog(gen x):
        r"""
        The principal branch of the dilogarithm of $x$, i.e. the analytic
        continuation of the power series $\log_2(x) = \sum_{n>=1} x^n/n^2$.

        EXAMPLES:
            sage: pari(1).dilog()
            1.644934066848226436472415167              # 32-bit
            1.6449340668482264364724151666460251892    # 64-bit
            sage: pari('1+I').dilog()
            0.6168502750680849136771556875 + 1.460362116753119547679775739*I    # 32-bit
            0.61685027506808491367715568749225944595 + 1.4603621167531195476797757394917875976*I   # 64-bit
        """
        _sig_on
        return P.new_gen(dilog(x.g, prec))

    def eint1(gen x, long n=0):
        r"""
        x.eint1({n}): exponential integral E1(x):
        $$
            \int_{x}^{\infty} \frac{e^{-t}}{t} dt
        $$
        If n is present, output the vector
            [eint1(x), eint1(2*x), ..., eint1(n*x)].
        This is faster than repeatedly calling eint1(i*x).

        REFERENCE: See page 262, Prop 5.6.12, of Cohen's book
        "A Course in Computational Algebraic Number Theory".

        EXAMPLES:

        """
        if n <= 0:
            _sig_on
            return P.new_gen(eint1(x.g, prec))
        else:
            _sig_on
            return P.new_gen(veceint1(x.g, stoi(n), prec))

    def erfc(gen x):
        r"""
        Return the complementary error function:
             $$(2/\sqrt{\pi}) \int_{x}^{\infty} e^{-t^2} dt.$$

        EXAMPLES:
            sage: pari(1).erfc()
            0.1572992070502851306587793649                # 32-bit
            0.15729920705028513065877936491739074070      # 64-bit
        """
        _sig_on
        return P.new_gen(gerfc(x.g, prec))

    def eta(gen x, flag=0):
        r"""
        x.eta({flag=0}): if flag=0, $\eta$ function without the $q^{1/24}$;
        otherwise $\eta$ of the complex number $x$ in the upper half plane
        intelligently computed using $\SL(2,\Z)$ transformations.

        DETAILS: This functions computes the following.  If the input
        $x$ is a complex number with positive imaginary part, the
        result is $\prod_{n=1}^{\infty} (q-1^n)$, where $q=e^{2 i \pi
        x}$.  If $x$ is a power series (or can be converted to a power
        series) with positive valuation, the result it
        $\prod_{n=1}^{\infty} (1-x^n)$.

        EXAMPLES:
            sage: pari('I').eta()
            0.9981290699259585132799623222             # 32-bit
            0.99812906992595851327996232224527387813   # 64-bit
        """
        if flag == 1:
            _sig_on
            return P.new_gen(trueeta(x.g, prec))
        _sig_on
        return P.new_gen(eta(x.g, prec))

    def exp(gen self):
        """
        x.exp(): exponential of x.

        EXAMPLES:
            sage: pari(0).exp()
            1.000000000000000000000000000               # 32-bit
            1.0000000000000000000000000000000000000     # 64-bit
            sage: pari(1).exp()
            2.718281828459045235360287471               # 32-bit
            2.7182818284590452353602874713526624978     # 64-bit
            sage: pari('x+O(x^8)').exp()
            1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + 1/720*x^6 + 1/5040*x^7 + O(x^8)
        """
        _sig_on
        return P.new_gen(gexp(self.g, prec))

    def gamma(gen s, precision=0):
        """
        s.gamma({precision}): Gamma function at s.

        INPUT:
            s -- gen (real or complex number
            precision -- optional precisiion

        OUTPUT:
            gen -- value of the Gamma function at s.

        EXAMPLES:
            sage: pari(2).gamma()
            1.000000000000000000000000000              # 32-bit
            1.0000000000000000000000000000000000000    # 64-bit
            sage: pari(5).gamma()
            24.00000000000000000000000000              # 32-bit
            24.000000000000000000000000000000000000    # 64-bit
            sage: pari('1+I').gamma()
            0.4980156681183560427136911175 - 0.1549498283018106851249551305*I    # 32-bit
            0.49801566811835604271369111746219809195 - 0.15494982830181068512495513048388660520*I    # 64-bit
        """
        if not precision: precision = prec
        _sig_on
        return P.new_gen(ggamma(s.g, precision))

    def gammah(gen s):
        """
        x.gammah(): Gamma function evaluated at the argument x+1/2,
        for x an integer.

        EXAMPLES:
            sage: pari(2).gammah()
            1.329340388179137020473625613               # 32-bit
            1.3293403881791370204736256125058588871     # 64-bit
            sage: pari(5).gammah()
            52.34277778455352018114900849               # 32-bit
            52.342777784553520181149008492418193679     # 64-bit
            sage: pari('1+I').gammah()
            0.5753151880634517207185443722 + 0.08821067754409390991246464371*I     # 32-bit
            0.57531518806345172071854437217501119058 + 0.088210677544093909912464643706507454993*I     # 64-bit
        """
        _sig_on
        return P.new_gen(ggamd(s.g, prec))

    def hyperu(gen a, b, x):
        r"""
        a.hyperu(b,x): U-confluent hypergeometric function.

	WARNING/TODO: This function is \emph{extremely slow} as
        implemented when used from the C library.  If you use the GP
        interpreter inteface it is vastly faster, so clearly this
        issue could be fixed with a better understanding of GP/PARI.
        Volunteers?

        EXAMPLES:
            sage.: pari(1).hyperu(2,3)
            0.3333333333333333333333333333              # 32-bit
            0.33333333333333333333333333333333333333    # 64-bit
        """
        t0GEN(b)
        t1GEN(x)
        _sig_on
        return P.new_gen(hyperu(a.g, t0, t1, prec))


    def incgam(gen s, x, y=None, precision=0):
        r"""
        s.incgam(x, {y}, {precision}): incomplete gamma function. y
        is optional and is the precomputed value of gamma(s).

        NOTE: This function works for any complex input (unlike in
        older version of PARI).

        INPUT:
            s, x, y -- gens
            precision -- option precision

        OUTPUT:
            gen -- value of the incomplete Gamma function at s.

        EXAMPLES:
            sage.: pari('1+I').incgam('3-I')
            -0.04582978599199457259586742326 + 0.04336968187266766812050474478*I        # 32-bit
            -0.045829785991994572595867423261490338705 + 0.043369681872667668120504744775954724733*I    # 64-bit
        """
        if not precision:
            precision = prec
        t0GEN(x)
        if y is None:
            _sig_on
            return P.new_gen(incgam(s.g, t0, precision))
        else:
            t1GEN(y)
            _sig_on
            return P.new_gen(incgam0(s.g, t0, t1, precision))

    def incgamc(gen s, x):
        r"""
        s.incgamc(x): complementary incomplete gamma function.

        The arguments $x$ and $s$ are complex numbers such that $s$ is
        not a pole of $\Gamma$ and $|x|/(|s|+1)$ is not much larger
        than $1$ (otherwise, the convergence is very slow).  The
        function returns the value of the integral $\int_{0}^{x}
        e^{-t} t^{s-1} dt.$

        EXAMPLES:
            sage: pari(1).incgamc(2)
            0.8646647167633873081060005050               # 32-bit
            0.86466471676338730810600050502751559659     # 64-bit
        """
        t0GEN(x)
        _sig_on
        return P.new_gen(incgamc(s.g, t0, prec))


    def log(gen self):
        r"""
        x.log(): natural logarithm of x.

        This function returns the principal branch of the natural
        logarithm of $x$, i.e., the branch such that $\Im(\log(x)) \in
        ]-\pi, \pi].$ The result is complex (with imaginary part equal
        to $\pi$) if $x\in \R$ and $x<0$.  In general, the algorithm
        uses the formula
        $$
            \log(x) \simeq \frac{\pi}{2{\rm agm}(1,4/s)} - m\log(2),
        $$
        if $s=x 2^m$ is large enough.  (The result is exact to $B$
        bits provided that $s>2^{B/2}$.)  At low accuracies, this
        function computes $\log$ using the series expansion near $1$.

        Note that $p$-adic arguments can also be given as input,
        with the convention that $\log(p)=0$.  Hence, in
        particular, $\exp(\log(x))/x$ is not in general
        equal to $1$ but instead to a $(p-1)th$ root of
        unity (or $\pm 1$ if $p=2$) times a power of $p$.

        EXAMPLES:
            sage: pari(5).log()
            1.609437912434100374600759333                 # 32-bit
            1.6094379124341003746007593332261876395       # 64-bit
            sage: pari('I').log()
            0.E-250 + 1.570796326794896619231321692*I             # 32-bit
            0.E-693 + 1.5707963267948966192313216916397514421*I   # 64-bit
        """
        _sig_on
        return P.new_gen(glog(self.g, prec))

    def lngamma(gen x):
        r"""
        x.lngamma(): logarithm of the gamma function of x.

        This function returns the principal branch of the logarithm of
        the gamma function of $x$.  The function $\log(\Gamma(x))$ is
        analytic on the complex plane with non-positive integers
        removed.  This function can have much larger inputs than
        $\Gamma$ itself.

        The $p$-adic analogue of this function is unfortunately not
        implemented.

        EXAMPLES:
            sage: pari(100).lngamma()
            359.1342053695753987760440105                # 32-bit
            359.13420536957539877604401046028690961      # 64-bit
        """
        _sig_on
        return P.new_gen(glngamma(x.g,prec))

    def polylog(gen x, long m, flag=0):
        """
        x.polylog(m,{flag=0}): m-th polylogarithm of x. flag is
        optional, and can be 0: default, 1: D_m~-modified m-th polylog
        of x, 2: D_m-modified m-th polylog of x, 3: P_m-modified m-th
        polylog of x.

        TODO: Add more explanation, copied from the PARI manual.

        EXAMPLES:
            sage: pari(10).polylog(3)
            5.641811414751341250600725771 - 8.328202076980270580884185850*I                          # 32-bit
            5.6418114147513412506007257705287671110 - 8.3282020769802705808841858505904310076*I      # 64-bit
            sage: pari(10).polylog(3,0)
            5.641811414751341250600725771 - 8.328202076980270580884185850*I                          # 32-bit
            5.6418114147513412506007257705287671110 - 8.3282020769802705808841858505904310076*I      # 64-bit
            sage: pari(10).polylog(3,1)
            0.5237784535024110488342571116              # 32-bit
            0.52377845350241104883425711161605950842    # 64-bit
            sage: pari(10).polylog(3,2)
            -0.4004590561634505605364328952             # 32-bit
            -0.40045905616345056053643289522452400363   # 64-bit
        """
        _sig_on
        return P.new_gen(polylog0(m, x.g, flag, prec))

    def psi(gen x):
        r"""
        x.psi(): psi-function at x.

        Return the $\psi$-function of $x$, i.e., the logarithmic
        derivative $\Gamma'(x)/\Gamma(x)$.

        EXAMPLES:
            sage: pari(1).psi()
            -0.5772156649015328606065120901              # 32-bit
            -0.57721566490153286060651209008240243104    # 64-bit
        """
        _sig_on
        return P.new_gen(gpsi(x.g, prec))


    def sin(gen x):
        """
        x.sin(): The sine of x.

        EXAMPLES:
            sage: pari(1).sin()
            0.8414709848078965066525023216             # 32-bit
            0.84147098480789650665250232163029899962   # 64-bit
            sage: pari('1+I').sin()
            1.298457581415977294826042366 + 0.6349639147847361082550822030*I                       # 32-bit
            1.2984575814159772948260423658078156203 + 0.63496391478473610825508220299150978151*I   # 64-bit
        """
        _sig_on
        return P.new_gen(gsin(x.g, prec))

    def sinh(gen self):
        """
        The hyperbolic sine function.

        EXAMPLES:
            sage: pari(0).sinh()
            0.E-250   # 32-bit
            0.E-693   # 64-bit
            sage: pari('1+I').sinh()
            0.6349639147847361082550822030 + 1.298457581415977294826042366*I                         # 32-bit
            0.63496391478473610825508220299150978151 + 1.2984575814159772948260423658078156203*I     # 64-bit
        """
        _sig_on
        return P.new_gen(gsh(self.g, prec))

    def sqr(gen x):
        """
        x.sqr(): square of x. NOT identical to x*x.

        TODO: copy extensive notes about this function
        from PARI manual.  Put examples below.

        EXAMPLES:
            sage: pari(2).sqr()
            4
        """
        _sig_on
        return P.new_gen(gsqr(x.g))


    def sqrt(gen x, precision=0):
        """
        x.sqrt({precision}): The square root of x.

        EXAMPLES:
            sage: pari(2).sqrt()
            1.414213562373095048801688724               # 32-bit
            1.4142135623730950488016887242096980786     # 64-bit
        """
        if not precision: precision = prec
        _sig_on
        return P.new_gen(gsqrt(x.g, precision))

    def sqrtn(gen x, n):
        r"""
        x.sqrtn(n): return the principal branch of the n-th root
        of x, i.e., the one such that
              $\arg(\sqrt(x)) \in ]-\pi/n, \pi/n]$.
        Also returns a second argument which is a suitable root
        of unity allowing one to recover all the other roots.
        If it was not possible to find such a number, then this
        second return value is 0.  If the argument is present and
        no square root exists, return 0 instead of raising an error.

        NOTE: intmods (modulo a prime) and $p$-adic numbers are
        allowed as arguments.

        INPUT:
            x -- gen
            n -- integer
        OUTPUT:
            gen -- principal n-th root of x
            gen -- z that gives the other roots

        EXAMPLES:
            sage: s, z = pari(2).sqrtn(5)
            sage: z
            0.3090169943749474241022934172 + 0.9510565162951535721164393334*I                         # 32-bit
            0.30901699437494742410229341718281905886 + 0.95105651629515357211643933337938214340*I     # 64-bit
            sage: s
            1.148698354997035006798626947               # 32-bit
            1.1486983549970350067986269467779275894     # 64-bit
            sage: s^5
            2.000000000000000000000000000               # 32-bit
            2.0000000000000000000000000000000000000     # 64-bit
            sage: (s*z)^5
            2.000000000000000000000000000 - 1.396701498 E-250*I                      # 32-bit
            2.0000000000000000000000000000000000000 - 1.0689317613194482765 E-693*I  # 64-bit
        """
	# TODO: ???  lots of good examples in the PARI docs ???
        cdef GEN zetan
        t0GEN(n)
        _sig_on
        ans = P.new_gen_noclear(gsqrtn(x.g, t0, &zetan, prec))
        return ans, P.new_gen(zetan)

    def tan(gen x):
        """
        x.tan() -- tangent of x

        EXAMPLES:
            sage: pari(2).tan()
            -2.185039863261518991643306102                   # 32-bit
            -2.1850398632615189916433061023136825434         # 64-bit
            sage: pari('I').tan()
            0.E-250 + 0.7615941559557648881194582826*I            # 32-bit
            0.E-693 + 0.76159415595576488811945828260479359041*I  # 64-bit
        """
        _sig_on
        return P.new_gen(gtan(x.g, prec))

    def tanh(gen x):
        """
        x.tanh() -- hyperbolic tangent of x

        EXAMPLES:
            sage: pari(1).tanh()
            0.7615941559557648881194582826             # 32-bit
            0.76159415595576488811945828260479359041   # 64-bit
            sage: pari('I').tanh()
            0.E-250 + 1.557407724654902230506974807*I              # 32-bit
            -5.344658806597241382 E-694 + 1.5574077246549022305069748074583601731*I  # 64-bit
        """
        _sig_on
        return P.new_gen(gth(x.g, prec))

    def teichmuller(gen x):
        r"""
        teichmuller(x): teichmuller character of p-adic number x.

        This is the unique $(p-1)$th root of unity congruent to
        $x/p^{v_p(x)}$ modulo $p$.

        EXAMPLES:
            sage: pari('2+O(7^5)').teichmuller()
            2 + 4*7 + 6*7^2 + 3*7^3 + O(7^5)
        """
        _sig_on
        return P.new_gen(teich(x.g))

    def theta(gen q, z):
        """
        q.theta(z): Jacobi sine theta-function.

        EXAMPLES:
            sage: pari('0.5').theta(2)
            1.632025902952598833772353216               # 32-bit
            1.6320259029525988337723532162682089972     # 64-bit
        """
        t0GEN(z)
        _sig_on
        return P.new_gen(theta(q.g, t0, prec))

    def thetanullk(gen q, long k):
        """
        q.thetanullk(k): return the k-th derivative at z=0 of theta(q,z)
        EXAMPLES:
            sage: pari('0.5').thetanullk(1)
            0.5489785325603405618549383537             # 32-bit
            0.54897853256034056185493835370857284861   # 64-bit
        """
        _sig_on
        return P.new_gen(thetanullk(q.g, k, prec))

    def weber(gen x, flag=0):
        r"""
        x.weber({flag=0}): One of Weber's f function of x.
        flag is optional, and can be
           0: default, function f(x)=exp(-i*Pi/24)*eta((x+1)/2)/eta(x)
              such that $j=(f^{24}-16)^3/f^{24}$,
           1: function f1(x)=eta(x/2)/eta(x) such that
              $j=(f1^24+16)^3/f2^{24}$,
           2: function f2(x)=sqrt(2)*eta(2*x)/eta(x) such that
              $j=(f2^{24}+16)^3/f2^{24}$.

        TODO: Add further explanation from PARI manual.

        EXAMPLES:
            sage: pari('I').weber()
            1.189207115002721066717499971 - 6.98350749 E-251*I     # 32-bit
            1.1892071150027210667174999705604759153 + 0.E-693*I    # 64-bit
            sage: pari('I').weber(1)
            1.090507732665257659207010656             # 32-bit
            1.0905077326652576592070106557607079790   # 64-bit
            sage: pari('I').weber(2)
            1.090507732665257659207010656             # 32-bit
            1.0905077326652576592070106557607079790   # 64-bit
        """
        _sig_on
        return P.new_gen(weber0(x.g, flag, prec))


    def zeta(gen s, precision=0):
        """
        zeta(s): Riemann zeta function at s with s a complex
                 or a p-adic number.

        TODO: Add extensive explanation from PARI user's manual.

        INPUT:
            s -- gen (real or complex number)

        OUTPUT:
            gen -- value of zeta at s.

        EXAMPLES:
            sage: pari(2).zeta()
            1.644934066848226436472415167             # 32-bit
            1.6449340668482264364724151666460251892   # 64-bit
            sage: pari('Pi^2/6')
            1.644934066848226436472415167             # 32-bit
            1.6449340668482264364724151666460251892   # 64-bit
            sage: pari(3).zeta()
            1.202056903159594285399738162             # 32-bit
            1.2020569031595942853997381615114499908   # 64-bit
        """
        if not precision: precision = prec
        _sig_on
        return P.new_gen(gzeta(s.g, precision))

    ###########################################
    # 4: NUMBER THEORETICAL functions
    ###########################################

    def bezout(gen x, y):
        cdef gen u, v, g
        cdef GEN U, V, G
        t0GEN(y)
        _sig_on
        G = gbezout(x.g, t0, &U, &V)
        _sig_off
        g = P.new_gen_noclear(G)
        u = P.new_gen_noclear(U)
        v = P.new_gen(V)
        return g, u, v

    def binomial(gen x, long k):
        _sig_on
        return P.new_gen(binome(x.g, k))

    def contfrac(gen x, b=0, long lmax=0):
        """
        contfrac(x,{b},{lmax}): continued fraction expansion of x (x
        rational, real or rational function). b and lmax are both
        optional, where b is the vector of numerators of the continued
        fraction, and lmax is a bound for the number of terms in the
        continued fraction expansion.
        """
        t0GEN(b)
        _sig_on
        return P.new_gen(contfrac0(x.g, t0, lmax))

    def contfracpnqn(gen x, b=0, long lmax=0):
        """
        contfracpnqn(x): [p_n,p_{n-1}; q_n,q_{n-1}] corresponding to the continued
        fraction x.
        """
        _sig_on
        return P.new_gen(pnqn(x.g))


    def gcd(gen x, y, long flag=0):
        """
        gcd(x,{y},{flag=0}): greatest common divisor of x and y. flag
        is optional, and can be 0: default, 1: use the modular gcd
        algorithm (x and y must be polynomials), 2 use the
        subresultant algorithm (x and y must be polynomials)
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(gcd0(x.g, t0, flag))

    def issquare(gen x, find_root=False):
        """
        issquare(x,{&n}): true(1) if x is a square, false(0) if not.
        If find_root is given, also returns the exact square root if
        it was computed.
        """
        cdef GEN G, t
        cdef gen g
        if find_root:
            _sig_on
            t = gcarrecomplet(x.g, &G)
            _sig_off
            v = bool(P.new_gen_noclear(t))
            if v:
                return v, P.new_gen(G)
            else:
                return v, None
        else:
            _sig_on
            return P.new_gen(gcarreparfait(x.g))


    def issquarefree(gen self):
        """
        EXAMPLES:
            sage: pari(10).issquarefree()
            True
            sage: pari(20).issquarefree()
            False
        """
        _sig_on
        t = bool(issquarefree(self.g))
        _sig_off
        return t

    def lcm(gen x, y):
        """
        Return the least common multiple of x and y.
        EXAMPLES:
            sage: pari(10).lcm(15)
            30
        """
        t0GEN(y)
        _sig_on
        return P.new_gen(glcm(x.g, t0))

    def numdiv(gen n):
        """
        Return the number of divisors of the integer n.

        EXAMPLES:
            sage: pari(10).numdiv()
            4
        """
        _sig_on
        return P.new_gen(gnumbdiv(n.g))

    def phi(gen n):
        """
        Return the Euler phi function of n.
        EXAMPLES:
            sage: pari(10).phi()
            4
        """
        _sig_on
        return P.new_gen(phi(n.g))

    def primepi(gen x):
        """
        Return the number of primes $\leq x$.

        EXAMPLES:
            sage: pari(7).primepi()
            4
            sage: pari(100).primepi()
            25
            sage: pari(1000).primepi()
            168
            sage: pari(100000).primepi()
            9592
        """
        _sig_on
        return P.new_gen(primepi(x.g))

    def sumdiv(gen n):
        """
        Return the sum of the divisors of $n$.

        EXAMPLES:
            sage: pari(10).sumdiv()
            18
        """
        _sig_on
        return P.new_gen(sumdiv(n.g))

    def sumdivk(gen n, long k):
        """
        Return the sum of the k-th powers of the divisors of n.

        EXAMPLES:
            sage: pari(10).sumdivk(2)
            130
        """
        _sig_on
        return P.new_gen(sumdivk(n.g, k))

    def xgcd(gen x, y):
        """
        Returns u,v,d such that d=gcd(x,y) and u*x+v*y=d.

        EXAMPLES:
            sage: pari(10).xgcd(15)
            (5, -1, 1)
        """
        return x.bezout(y)


    ##################################################
    # 5: Elliptic curve functions
    ##################################################

    def ellinit(self, int flag=0, precision=0):
        if not precision:
            precision = prec
        _sig_on
        return P.new_gen(ellinit0(self.g, flag, precision))

    def ellglobalred(self):
        _sig_on
        return self.new_gen(globalreduction(self.g))

    def elladd(self, z0, z1):
        """
        elladd(self, z0, z1)

        Sum of the points z0 and z1 on this elliptic curve.

        INPUT:
            self -- elliptic curve E
            z0 -- point on E
            z1 -- point on E

        OUTPUT:
            point on E

        EXAMPLES:
        First we create an elliptic curve:

            sage: e = pari([0, 1, 1, -2, 0]).ellinit()
            sage: str(e)[:65]   # first part of output
            '[0, 1, 1, -2, 0, 4, -4, 1, -3, 112, -856, 389, 1404928/389, [0.90'

        Next we add two points on the elliptic curve.  Notice that
        the Python lists are automatically converted to PARI objects so
        you don't have to do that explicitly in your code.

            sage: e.elladd([1,0,1], [-1,1,1])
            [-3/4, -15/8]
        """
        t0GEN(z0); t1GEN(z1)
        _sig_on
        return self.new_gen(addell(self.g, t0, t1))

    def ellak(self, n):
        r"""
        e.ellak(n): Returns the coefficient $a_n$ of the $L$-function of
        the elliptic curve e, i.e. the coefficient of a newform of
        weight 2 newform.

        \begin{notice}
        The curve $e$ {\em must} be a medium or long vector of the type given
        by ellinit. For this function to work for every n and not just
        those prime to the conductor, e must be a minimal Weierstrass
        equation. If this is not the case, use the function
        ellminimalmodel first before using ellak (or you will get
        INCORRECT RESULTS!)
        \end{notice}

        INPUT:
            e -- a PARI elliptic curve.
            n -- integer ..

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellak(6)
            2
            sage: e.ellak(2005)
            2
            sage: e.ellak(-1)
            0
            sage: e.ellak(0)
            0
        """
        t0GEN(n)
        _sig_on
        return self.new_gen(akell(self.g, t0))


    def ellan(self, long n, python_ints=False):
        """
        Return the Fourier coefficients of the modular form attached
        to this elliptic curve.

        INPUT:
            n -- a long integer
            python_ints -- bool (default is False); if True, return a
                           list of Python ints instead of a PARI gen
                           wrapper.

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellan(3)
            [1, -2, -1]
            sage: e.ellan(20)
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2]
            sage: e.ellan(-1)
            []
            sage: v = e.ellan(10, python_ints=True); v
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2]
            sage: type(v)
            <type 'list'>
            sage: type(v[0])
            <type 'int'>
        """
        _sig_on
        cdef GEN g
        if python_ints:
            g = anell(self.g, n)
            v = [gtolong(<GEN> g[i+1]) for i in range(glength(g))]
            (<PariInstance>pari).clear_stack()
            return v
        else:
            return self.new_gen(anell(self.g, n))

    def ellap(self, p):
        r"""
        e.ellap(p): Returns the prime-indexed coefficient $a_p$ of the
        $L$-function of the elliptic curve $e$, i.e. the coefficient of a
        newform of weight 2 newform.

        \begin{notice}
        If p is not prime, this function will return an {\bf incorrect}
        answer.

        The curve e must be a medium or long vector of the type given
        by ellinit. For this function to work for every n and not just
        those prime to the conductor, e must be a minimal Weierstrass
        equation. If this is not the case, use the function
        ellminimalmodel first before using ellap (or you will get
        INCORRECT RESULTS!)
        \end{notice}

        INPUT:
            e -- a PARI elliptic curve.
            p -- prime integer ..

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellap(2)
            -2
            sage: e.ellap(2003)
            4
            sage: e.ellak(-1)
            0
        """
        t0GEN(p)
        _sig_on
        return self.new_gen(apell(self.g, t0))


    def ellaplist(self, long n, python_ints=False):
        r"""
        e.ellaplist(n): Returns a PARI list of all the prime-indexed
        coefficient $a_p$ of the $L$-function of the elliptic curve
        $e$, i.e. the coefficient of a newform of weight 2 newform.

        INPUT:
            n -- a long integer
            python_ints -- bool (default is False); if True, return a
                           list of Python ints instead of a PARI gen
                           wrapper.

        \begin{notice}
        The curve e must be a medium or long vector of the type given
        by ellinit. For this function to work for every n and not just
        those prime to the conductor, e must be a minimal Weierstrass
        equation. If this is not the case, use the function
        ellminimalmodel first before using ellanlist (or you will get
        INCORRECT RESULTS!)
        \end{notice}

        INPUT:
            e -- a PARI elliptic curve.
            n -- an integer

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: v = e.ellaplist(10); v
            [-2, -1, 1, -2]
            sage: type(v)
            <type 'sage.libs.pari.gen.gen'>
            sage: v.type()
            't_VEC'
            sage: e.ellan(10)
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2]
            sage: v = e.ellaplist(10, python_ints=True); v
            [-2, -1, 1, -2]
            sage: type(v)
            <type 'list'>
            sage: type(v[0])
            <type 'int'>
        """
        # 1. make a table of primes up to n.
        if n < 2:
            return self.new_gen(zerovec(0))
        cdef GEN g
        pari.init_primes(n+1)
        t0GEN(n)
        _sig_on
        g = primes(gtolong(primepi(t0)))

        # 2. Replace each prime in the table by apell of it.
        cdef long i

        if python_ints:
            v = [gtolong(apell(self.g, <GEN> g[i+1])) \
                        for i in range(glength(g))]
            (<PariInstance>pari).clear_stack()
            return v
        else:
            for i from 0 <= i < glength(g):
                g[i+1] = <long> apell(self.g, <GEN> g[i+1])
            return self.new_gen(g)


    def ellbil(self, z0, z1):
        """
        EXAMPLES:
            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: e.ellbil([1, 0, 1], [-1, 1, 1])
            0.4181889844988605856298894582              # 32-bit
             0.41818898449886058562988945821587638238   # 64-bit
        """
##         Increasing the precision does not increase the precision
##         result, since quantities related to the elliptic curve were
##         computed to low precision.
##             sage: set_real_precision(10)
##             sage: e.ellbil([1, 0, 1], [-1, 1, 1])
##             0.4181889844988605856298894585
##         However, if we recompute the elliptic curve after increasing
##         the precision, then the bilinear pairing will be computed to
##         higher precision as well.
##             sage: e = pari([0,1,1,-2,0]).ellinit()
##             sage: e.ellbil([1, 0, 1], [-1, 1, 1])
##             0.4181889844988605856298894582
##             sage: set_real_precision(5)
        t0GEN(z0); t1GEN(z1)
        _sig_on
        return self.new_gen(bilhell(self.g, t0, t1, prec))

    def ellchangecurve(self, ch):
        """
        EXAMPLES:
            sage: e = pari([1,2,3,4,5]).ellinit()
            sage: e.ellglobalred()
            [10351, [1, -1, 0, -1], 1]
            sage: f = e.ellchangecurve([1,-1,0,-1])
            sage: f[:5]
            [1, -1, 0, 4, 3]
        """
        t0GEN(ch)
        _sig_on
        return self.new_gen(coordch(self.g, t0))

    def elleta(self):
        _sig_on
        return self.new_gen(elleta(self.g, prec))

    def ellheight(self, a, flag=0):
        t0GEN(a)
        _sig_on
        return self.new_gen(ellheight0(self.g, t0, flag, prec))

    def ellheightmatrix(self, x):
        """
        ellheightmatrix(e,x)

        Returns the height matrix for vector of points x on elliptic curve e using
        theta functions.
        """
        t0GEN(x)
        _sig_on
        return self.new_gen(mathell(self.g, t0, prec))

    def ellisoncurve(self, x):
        t0GEN(x)
        _sig_on
        t = bool(oncurve(self.g, t0) == 1)
        _sig_off
        return t

    def elllocalred(self, p):
        t0GEN(p)
        _sig_on
        return self.new_gen(elllocalred(self.g, t0))

    def elllseries(self, s, A=1):
        t0GEN(s); t1GEN(A)
        _sig_on
        return self.new_gen(lseriesell(self.g, t0, t1, prec))

    def ellminimalmodel(self):
        """
        ellminimalmodel(e): return the standard minimal integral model
        of the rational elliptic curve e and the corresponding change
        of variables.
        INPUT:
            e -- gen (that defines an elliptic curve)
        OUTPUT:
            gen -- minimal model
            gen -- change of coordinates
        EXAMPLES:
            sage: e = pari([1,2,3,4,5]).ellinit()
            sage: F, ch = e.ellminimalmodel()
            sage: F[:5]
            [1, -1, 0, 4, 3]
            sage: ch
            [1, -1, 0, -1]
            sage: e.ellchangecurve(ch)[:5]
            [1, -1, 0, 4, 3]
        """
        cdef GEN x, y
        cdef gen model, change
        cdef pari_sp t
        _sig_on
        x = ellminimalmodel(self.g, &y)
        change = self.new_gen_noclear(y)
        model = self.new_gen(x)
        return model, change

    def ellorder(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(orderell(self.g, t0))

    def ellordinate(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(ordell(self.g, t0, prec))

    def ellpointtoz(self, P):
        t0GEN(P)
        _sig_on
        return self.new_gen(zell(self.g, t0, prec))

    def ellpow(self, z, n):
        t0GEN(z); t1GEN(n)
        _sig_on
        return self.new_gen(powell(self.g, t0, t1))

    def ellrootno(self, p=1):
        t0GEN(p)
        _sig_on
        return ellrootno(self.g, t0)

    def ellsigma(self, z, flag=0):
        t0GEN(z)
        _sig_on
        return self.new_gen(ellsigma(self.g, t0, flag, prec))

    def ellsub(self, z1, z2):
        t0GEN(z1); t1GEN(z2)
        _sig_on
        return self.new_gen(subell(self.g, t0, t1))

    def elltaniyama(self):
        _sig_on
        return self.new_gen(taniyama(self.g))

    def elltors(self, flag=0):
        _sig_on
        return self.new_gen(elltors0(self.g, flag))


    def ellzeta(self, z):
        t0GEN(z)
        _sig_on
        return self.new_gen(ellzeta(self.g, t0, prec))

    def ellztopoint(self, z):
        t0GEN(z)
        _sig_on
        return self.new_gen(pointell(self.g, t0, prec))

    def omega(self):
        """
        e.omega(): return basis for the period lattice of the elliptic curve e.

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.omega()
            [1.269209304279553421688794617, 0.6346046521397767108443973084 + 1.458816616938495229330889613*I]   # 32-bit
            [1.2692093042795534216887946167545473052, 0.63460465213977671084439730837727365260 + 1.4588166169384952293308896129036752572*I]   # 64-bit
        """
        return self[14:16]

    def disc(self):
        """
        e.disc(): return the discriminant of the elliptic curve e.

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.disc()
            -161051
            sage: _.factor()
            [-1, 1; 11, 5]
        """
        return self[11]

    def j(self):
        """
        e.j(): return the j-invariant of the elliptic curve e.

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.j()
            -122023936/161051
            sage: _.factor()
            [-1, 1; 2, 12; 11, -5; 31, 3]
        """
        return self[12]




    def ellj(self):
        _sig_on
        return P.new_gen(jell(self.g, prec))


    ###########################################
    # 6: Functions related to general NUMBER FIELDS
    ###########################################
    def bnfcertify(self):
        r"""
        \code{bnf} being as output by \code{bnfinit}, checks whether
        the result is correct, i.e. whether the calculation of the
        contents of self are correct without assuming the Generalized
        Riemann Hypothesis. If it is correct, the answer is 1. If not,
        the program may output some error message, but more probably
        will loop indefinitely. In \emph{no} occasion can the program
        give a wrong answer (barring bugs of course): if the program
        answers 1, the answer is certified.

        \note{WARNING! By default, most of the bnf routines depend on
        the correctness of a heuristic assumption which is stronger
        than GRH.  In order to obtain a provably-correct result you
        \emph{must} specify $c=c_2=12$ for the technical optional
        parameters to the function. There are known counterexamples
        for smaller $c$ (which is the default).}

        """
        cdef long n
        _sig_on
        n = certifybuchall(self.g)
        _sig_off
        return n

    def bnfinit(self, long flag=0, tech=None):
        if tech is None:
            _sig_on
            return P.new_gen(bnfinit0(self.g, flag, <GEN>0, prec))
        else:
            t0GEN(tech)
            _sig_on
            return P.new_gen(bnfinit0(self.g, flag, t0, prec))

    def bnfisintnorm(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(bnfisintnorm(self.g, t0))

    def bnfisprincipal(self, x, long flag=1):
        t0GEN(x)
        _sig_on
        return self.new_gen(isprincipalall(self.g, t0, flag))

    def bnfnarrow(self):
        _sig_on
        return self.new_gen(buchnarrow(self.g))

    def bnfunit(self):
        _sig_on
        return self.new_gen(buchfu(self.g))

    def dirzetak(self, n):
        t0GEN(n)
        _sig_on
        return self.new_gen(dirzetak(self.g, t0))

    def idealadd(self, x, y):
        t0GEN(x); t1GEN(y)
        _sig_on
        return self.new_gen(idealadd(self.g, t0, t1))

    def idealdiv(self, x, y, long flag=0):
        t0GEN(x); t1GEN(y)
        _sig_on
        return self.new_gen(idealdiv0(self.g, t0, t1, flag))

    def idealfactor(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(idealfactor(self.g, t0))

    def idealhnf(self, a, b=None):
        t0GEN(a)
        _sig_on
        if b is None:
            return self.new_gen(idealhermite(self.g, t0))
        else:
            t1GEN(b)
            return self.new_gen(idealhnf0(self.g, t0, t1))

    def idealmul(self, x, y, long flag=0):
        t0GEN(x); t1GEN(y)
        _sig_on
        if flag == 0:
            return self.new_gen(idealmul(self.g, t0, t1))
        else:
            return self.new_gen(idealmulred(self.g, t0, t1, prec))

    def idealnorm(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(idealnorm(self.g, t0))

    def idealtwoelt(self, x, a=None):
        t0GEN(x)
        _sig_on
        if a is None:
            return self.new_gen(ideal_two_elt0(self.g, t0, NULL))
        else:
            t1GEN(a)
            return self.new_gen(ideal_two_elt0(self.g, t0, t1))

    def idealval(self, x, p):
        t0GEN(x); t1GEN(p)
        _sig_on
        return idealval(self.g, t0, t1)

    def modreverse(self):
        """
        modreverse(x): reverse polymod of the polymod x, if it exists.

        EXAMPLES:
        """
        _sig_on
        return self.new_gen(polymodrecip(self.g))

    def nfbasis(self, long flag=0, p=0):
        cdef gen _p
        cdef GEN g
        if p != 0:
            _p = self.pari(p)
            g = _p.g
        else:
            g = <GEN>NULL
        _sig_on
        return self.new_gen(nfbasis0(self.g, flag, g))

    def nffactor(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(nffactor(self.g, t0))

    def nfgenerator(self):
        f = self[0]
        x = f.variable()
        return x.Mod(f)

    def nfinit(self, long flag=0):
        _sig_on
        return P.new_gen(nfinit0(self.g, flag, prec))

    def rnfcharpoly(self, T, a, v='x'):
        t0GEN(T); t1GEN(a); t2GEN(v)
        _sig_on
        return self.new_gen(rnfcharpoly(self.g, t0, t1, gvar(t2)))

    def rnfdisc(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(rnfdiscf(self.g, t0))

    def rnfeltabstorel(self, x):
        t0GEN(x)
        _sig_on
        polymodmod = self.new_gen(rnfelementabstorel(self.g, t0))
        return polymodmod.centerlift().centerlift()

    def rnfeltreltoabs(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(rnfelementreltoabs(self.g, t0))

    def rnfequation(self, poly, long flag=0):
        t0GEN(poly)
        _sig_on
        return self.new_gen(rnfequation0(self.g, t0, flag))

    def rnfidealabstorel(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(rnfidealabstorel(self.g, t0))

    def rnfidealhnf(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(rnfidealhermite(self.g, t0))

    def rnfidealnormrel(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(rnfidealnormrel(self.g, t0))

    def rnfidealreltoabs(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(rnfidealreltoabs(self.g, t0))

    def rnfidealtwoelt(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(rnfidealtwoelement(self.g, t0))

    def rnfinit(self, poly):
        """
        EXAMPLES:
        We construct a relative number field.
            sage: f = pari('y^3+y+1')
            sage: K = f.nfinit()
            sage: x = pari('x'); y = pari('y')
            sage: g = x^5 - x^2 + y
            sage: L = K.rnfinit(g)
        """
        t0GEN(poly)
        _sig_on
        return P.new_gen(rnfinitalg(self.g, t0, prec))

    def rnfisfree(self, poly):
        t0GEN(poly)
        _sig_on
        return rnfisfree(self.g, t0)





    ##################################################
    # 7: POLYNOMIALS and power series
    ##################################################
    def reverse(self):
        """
        Return the polynomial obtained by reversing the coefficients
        of this polynomial.
        """
        return self.Vec().Polrev()

    def deriv(self, v=-1):
        _sig_on
        return self.new_gen(deriv(self.g, self.get_var(v)))

    def eval(self, x):
        t0GEN(x)
        _sig_on
        return self.new_gen(poleval(self.g, t0))

    def __call__(self, x):
        return self.eval(x)

    def factorpadic(self, p, long r=20, long flag=0):
        """
        self.factorpadic(p,{r=20},{flag=0}): p-adic factorization of the
        polynomial x to precision r. flag is optional and may be set
        to 0 (use round 4) or 1 (use Buchmann-Lenstra)
        """
        t0GEN(p)
        _sig_on
        return self.new_gen(factorpadic0(self.g, t0, r, flag))

    def factormod(self, p, long flag=0):
        """
        x.factormod(p,{flag=0}): factorization mod p of the polynomial
        x using Berlekamp. flag is optional, and can be 0: default or
        1: simple factormod, same except that only the degrees of the
        irreducible factors are given.
        """
        t0GEN(p)
        _sig_on
        return self.new_gen(factormod0(self.g, t0, flag))

    def intformal(self, y=-1):
        """
        x.intformal({y}): formal integration of x with respect to the main
        variable of y, or to the main variable of x if y is omitted
        """
        _sig_on
        return self.new_gen(integ(self.g, self.get_var(y)))

    def padicappr(self, a):
        """
        x.padicappr(a): p-adic roots of the polynomial x congruent to a mod p
        """
        t0GEN(a)
        _sig_on
        return self.new_gen(padicappr(self.g, t0))

    def newtonpoly(self, p):
        """
        self.newtonpoly(p): Newton polygon of polynomial x with respect
        to the prime p.
        """
        t0GEN(p)
        _sig_on
        return self.new_gen(newtonpoly(self.g, t0))

    def polcoeff(self, long n, var=-1):
        """
        EXAMPLES:
            sage: f = pari("x^2 + y^3 + x*y")
            sage: f
            x^2 + y*x + y^3
            sage: f.polcoeff(1)
            y
            sage: f.polcoeff(3)
            0
            sage: f.polcoeff(3, "y")
            1
            sage: f.polcoeff(1, "y")
            x
        """
        _sig_on
        return self.new_gen(polcoeff0(self.g, n, self.get_var(var)))

    def polcompositum(self, pol2, long flag=0):
        t0GEN(pol2)
        _sig_on
        return self.new_gen(polcompositum0(self.g, t0, flag))

    def poldegree(self, var=-1):
        """
        f.poldegree(var={x}): Return the degree of this polynomial.
        """
        _sig_on
        n = poldegree(self.g, self.get_var(var))
        _sig_off
        return n

    def poldisc(self, var=-1):
        """
        f.poldist(var={x}):  Return the discriminant of this polynomial.
        """
        _sig_on
        return self.new_gen(poldisc0(self.g, self.get_var(var)))

    def poldiscreduced(self):
        _sig_on
        return self.new_gen(reduceddiscsmith(self.g))

    def polgalois(self):
        """
        f.polgalois(): Galois group of the polynomial f
        """
        _sig_on
        return self.new_gen(galois(self.g, prec))

    def polhensellift(self, y, p, long e):
        """
        self.polhensellift(y, p, e): lift the factorization y of
        self modulo p to a factorization modulo $p^e$ using Hensel
        lift. The factors in y must be pairwise relatively prime
        modulo p.
        """
        t0GEN(y)
        t1GEN(p)
        _sig_on
        return self.new_gen(polhensellift(self.g, t0, t1, e))

    def polisirreducible(self):
        """
        f.polisirreducible(): Returns True if f is an irreducible
        non-constant polynomial, or False if f is reducible or
        constant.
        """
        _sig_on
        return bool(self.new_gen(gisirreducible(self.g)))


    def pollead(self, v=-1):
        """
        self.pollead({v}): leading coefficient of polynomial or series
        self, or self itself if self is a scalar. Error
        otherwise. With respect to the main variable of self if v is
        omitted, with respect to the variable v otherwise
        """
        _sig_on
        return self.new_gen(pollead(self.g, self.get_var(v)))

    def polrecip(self):
        _sig_on
        return self.new_gen(polrecip(self.g))

    def polred(self, flag=0, fa=None):
        _sig_on
        if fa is None:
            return self.new_gen(polred0(self.g, flag, NULL))
        else:
            t0GEN(fa)
            return self.new_gen(polred0(self.g, flag, t0))

    def polresultant(self, y, var=-1, flag=0):
        t0GEN(y)
        _sig_on
        return self.new_gen(polresultant0(self.g, t0, self.get_var(var), flag))

    def polroots(self, flag=0):
        """
        polroots(x,{flag=0}): complex roots of the polynomial x. flag
        is optional, and can be 0: default, uses Schonhage's method modified
        by Gourdon, or 1: uses a modified Newton method.
        """
        _sig_on
        return self.new_gen(roots0(self.g, flag, prec))

    def polrootsmod(self, p, flag=0):
        t0GEN(p)
        _sig_on
        return self.new_gen(rootmod0(self.g, t0, flag))

    def polrootspadic(self, p, r=20):
        t0GEN(p)
        _sig_on
        return self.new_gen(rootpadic(self.g, t0, r))

    def polrootspadicfast(self, p, r=20):
        t0GEN(p)
        _sig_on
        return self.new_gen(rootpadicfast(self.g, t0, r))

    def polsturm(self, a, b):
        t0GEN(a)
        t1GEN(b)
        _sig_on
        n = sturmpart(self.g, t0, t1)
        _sig_off
        return n

    def polsylvestermatrix(self, g):
        t0GEN(g)
        _sig_on
        return self.new_gen(sylvestermatrix(self.g, t0))

    def polsym(self, long n):
        _sig_on
        return self.new_gen(polsym(self.g, n))

    def serconvol(self, g):
        t0GEN(g)
        _sig_on
        return self.new_gen(convol(self.g, t0))

    def serlaplace(self):
        _sig_on
        return self.new_gen(laplace(self.g))

    def serreverse(self):
        """
        serreverse(f): reversion of the power series f.

        If f(t) is a series in t with valuation 1, find the
        series g(t) such that g(f(t)) = t.

        EXAMPLES:
            sage: f = pari('x+x^2+x^3+O(x^4)'); f
            x + x^2 + x^3 + O(x^4)
            sage: g = f.serreverse(); g
            x - x^2 + x^3 + O(x^4)
            sage: f.subst('x',g)
            x + O(x^4)
            sage: g.subst('x',f)
            x + O(x^4)
        """
        _sig_on
        return self.new_gen(recip(self.g))

    def thueinit(self, flag=0):
        _sig_on
        return self.new_gen(thueinit(self.g, flag, prec))





    ###########################################
    # 8: Vectors, matrices, LINEAR ALGEBRA and sets
    ###########################################

    def vecextract(self, y, z=None):
        r"""
        self.vecextract(y,{z}): extraction of the components of the
        matrix or vector x according to y and z. If z is omitted, y
        designates columns, otherwise y corresponds to rows and z to
        columns. y and z can be vectors (of indices), strings
        (indicating ranges as in"1..10") or masks
        (integers whose binary representation indicates the indices
        to extract, from left to right 1, 2, 4, 8, etc.)

        \note{This function uses the PARI row and column indexing,
        so the first row or column is indexed by 1 instead of 0.}
        """
        t0GEN(y)
        if z is None:
            _sig_on
            return P.new_gen(extract(self.g, t0))
        else:
            t1GEN(z)
            _sig_on
            return P.new_gen(extract0(self.g, t0, t1))

    def ncols(self):
        """
        Return the number of rows of self.

        EXAMPLES:
            sage: pari('matrix(19,8)').ncols()
            8
        """
        cdef long n
        _sig_on
        n = glength(self.g)
        _sig_off
        return n

    def nrows(self):
        """
        Return the number of rows of self.

        EXAMPLES:
            sage: pari('matrix(19,8)').nrows()
            19
        """
        cdef long n
        _sig_on
        n = glength(<GEN>(self.g[1]))
        _sig_off
        return n

    def mattranspose(self):
        """
        Transpose of the matrix self.

        EXAMPLES:
            sage: pari('[1,2,3; 4,5,6;  7,8,9]').mattranspose()
            [1, 4, 7; 2, 5, 8; 3, 6, 9]
        """
        _sig_on
        return self.new_gen(gtrans(self.g)).Mat()

    def matadjoint(self):
        """
        matadjoint(x): adjoint matrix of x.

        EXAMPLES:
            sage: pari('[1,2,3; 4,5,6;  7,8,9]').matadjoint()
            [-3, 6, -3; 6, -12, 6; -3, 6, -3]
            sage: pari('[a,b,c; d,e,f; g,h,i]').matadjoint()
            [(i*e - h*f), (-i*b + h*c), (f*b - e*c); (-i*d + g*f), i*a - g*c, -f*a + d*c; (h*d - g*e), -h*a + g*b, e*a - d*b]
        """
        _sig_on
        return self.new_gen(adj(self.g)).Mat()

    def qflllgram(self, long flag=0):
        """
        qflllgram(x,{flag=0}): LLL reduction of the lattice whose gram
        matrix is x (gives the unimodular transformation matrix). flag
        is optional and can be 0: default,1: lllgramint algorithm for
        integer matrices, 4: lllgramkerim giving the kernel and the
        LLL reduced image, 5: lllgramkerimgen same when the matrix has
        polynomial coefficients, 8: lllgramgen, same as qflllgram when
        the coefficients are polynomials.
        """
        _sig_on
        return self.new_gen(qflllgram0(self.g,flag,prec)).Mat()

    def lllgram(self):
        return self.qflllgram(0)

    def lllgramint(self):
        return self.qflllgram(1)

    def qfminim(self, B, max, long flag=0):
        """
        qfminim(x,{bound},{maxnum},{flag=0}): number of vectors of
        square norm <= bound, maximum norm and list of vectors for the
        integral and definite quadratic form x; minimal non-zero
        vectors if bound=0. flag is optional, and can be 0: default;
        1: returns the first minimal vector found (ignore maxnum); 2:
        as 0 but uses a more robust, slower implementation, valid for
        non integral quadratic forms.
        """
        _sig_on
        t0GEN(B)
        t1GEN(max)
        return self.new_gen(qfminim0(self.g,t0,t1,flag,precdl))

    def qfrep(self, B, long flag=0):
        """
        qfrep(x,B,{flag=0}): vector of (half) the number of vectors of
        norms from 1 to B for the integral and definite quadratic form
        x. Binary digits of flag mean 1: count vectors of even norm
        from 1 to 2B, 2: return a t_VECSMALL instead of a t_VEC.
        """
        _sig_on
        t0GEN(B)
        return self.new_gen(qfrep0(self.g,t0,flag))

    def matsolve(self, B):
        """
        matsolve(B): Solve the linear system Mx=B for an invertible matrix M

        matsolve(B) uses gaussian elimination to solve Mx=B, where M is
        invertible and B is a column vector.

        The corresponding pari library routine is gauss. The gp-interface
        name matsolve has been given preference here.

        INPUT:
            B -- a column vector of the same dimension as the square matrix self

        EXAMPLES:
            sage: pari('[1,1;1,-1]').matsolve(pari('[1;0]'))
            [1/2; 1/2]
        """
        _sig_on
        t0GEN(B)
        return self.new_gen(gauss(self.g,t0))

    def matker(self, long flag=0):
        """
        Return a basis of the kernel of this matrix.

        INPUT:
            flag -- optional; may be set to
                0: default;
                non-zero: x is known to have integral entries.

        EXAMPLES:
            sage: pari('[1,2,3;4,5,6;7,8,9]').matker()
            [1; -2; 1]

        With algorithm 1, even if the matrix has integer entries the
        kernel need nto be saturated (which is weird):
            sage: pari('[1,2,3;4,5,6;7,8,9]').matker(1)
            [3; -6; 3]
            sage: pari('matrix(3,3,i,j,i)').matker()
            [-1, -1; 1, 0; 0, 1]
            sage: pari('[1,2,3;4,5,6;7,8,9]*Mod(1,2)').matker()
            [Mod(1, 2); Mod(0, 2); 1]
        """
        _sig_on
        return self.new_gen(matker0(self.g, flag))

    def matkerint(self, long flag=0):
        """
        Return the integer kernel of a matrix.

        This is the LLL-reduced Z-basis of the kernel of the matrix x
        with integral entries.

        INPUT:
            flag -- optional, and may be set to
                   0: default, uses a modified LLL,
                   1: uses matrixqz.

        EXAMPLES:
            sage: pari('[2,1;2,1]').matker()
            [-1/2; 1]
            sage: pari('[2,1;2,1]').matkerint()
            [-1; 2]

        This is worrisome (so be careful!):
            sage: pari('[2,1;2,1]').matkerint(1)
            Mat(1)
        """
        _sig_on
        return self.new_gen(matkerint0(self.g, flag))

    def matdet(self, long flag=0):
        """
        Return the determinant of this matrix.

        INPUT:
            flag  -- (optional) flag
                     0: using Gauss-Bareiss.
                     1: use classical gaussian elimination (slightly better for integer entries)

        EXAMPLES:
            sage: pari('[1,2; 3,4]').matdet(0)
            -2
            sage: pari('[1,2; 3,4]').matdet(1)
            -2
        """
        _sig_on
        return self.new_gen(det0(self.g, flag))

    def trace(self):
        """
        Return the trace of this PARI object.

        EXAMPLES:
            sage: pari('[1,2; 3,4]').trace()
            5
        """
        _sig_on
        return self.new_gen(gtrace(self.g))

    def mathnf(self, flag=0):
        """
        A.mathnf({flag=0}): (upper triangular) Hermite normal form of
        A, basis for the lattice formed by the columns of A.

        INPUT:
            flag -- optional, value range from 0 to 4 (0 if
            omitted), meaning :
                0: naive algorithm
                1: Use Batut's algorithm -- output 2-component vector
                   [H,U] such that H is the  HNF of A, and U is a
                   unimodular matrix such that xU=H.
                3: Use Batut's algorithm. Output [H,U,P] where P is
                   a permutation matrix such that P A U = H.
                4: As 1, using a heuristic variant of LLL reduction
                   along the way.

        EXAMPLES:
            sage: pari('[1,2,3; 4,5,6;  7,8,9]').mathnf()
            [6, 1; 3, 1; 0, 1]
        """
        _sig_on
        return self.new_gen(mathnf0(self.g, flag))

    def mathnfmod(self, d):
        """
        Returns the Hermite normal form if d is a multiple of the determinant

        Beware that PARI's concept of a hermite normal form is an upper
        triangular matrix with the same column space as the input matrix.

        INPUT:
            d -- multiple of the determinant of self

        EXAMPLES:
            sage: M=matrix([[1,2,3],[4,5,6],[7,8,11]])
	    sage: d=M.det()
	    sage: pari(M).mathnfmod(d)
            [6, 4, 3; 0, 1, 0; 0, 0, 1]

	Note that d really needs to be a multiple of the discriminant, not
	just of the exponent of the cokernel:

            sage: M=matrix([[1,0,0],[0,2,0],[0,0,6]])
	    sage: pari(M).mathnfmod(6)
	    [1, 0, 0; 0, 1, 0; 0, 0, 6]
	    sage: pari(M).mathnfmod(12)
	    [1, 0, 0; 0, 2, 0; 0, 0, 6]

        """
        _sig_on
        t0GEN(d)
        return self.new_gen(hnfmod(self.g, t0))

    def mathnfmodid(self, d):
        """
        Returns the Hermite Normal Form of M concatenated with d*Identity

        Beware that PARI's concept of a Hermite normal form is a maximal rank
        upper triangular matrix with the same column space as the input matrix.

        INPUT:
            d -- Determines

        EXAMPLES:
            sage: M=matrix([[1,0,0],[0,2,0],[0,0,6]])
	    sage: pari(M).mathnfmodid(6)
            [1, 0, 0; 0, 2, 0; 0, 0, 6]

	This routine is not completely equivalent to mathnfmod:

	    sage: pari(M).mathnfmod(6)
	    [1, 0, 0; 0, 1, 0; 0, 0, 6]
        """
        _sig_on
        t0GEN(d)
        return self.new_gen(hnfmodid(self.g, t0))

    def matsnf(self, flag=0):
        """
        x.matsnf({flag=0}): Smith normal form (i.e. elementary
        divisors) of the matrix x, expressed as a vector d. Binary
        digits of flag mean 1: returns [u,v,d] where d=u*x*v,
        otherwise only the diagonal d is returned, 2: allow polynomial
        entries, otherwise assume x is integral, 4: removes all
        information corresponding to entries equal to 1 in d.

        EXAMPLES:
            sage: pari('[1,2,3; 4,5,6;  7,8,9]').matsnf()
            [0, 3, 1]
        """
        _sig_on
        return self.new_gen(matsnf0(self.g, flag))

    def matfrobenius(self, flag=0):
        r"""
        M.matfrobenius({flag=0}): Return the Frobenius form of the
        square matrix M. If flag is 1, return only the elementary
        divisors (a list of polynomials). If flag is 2, return a
        two-components vector [F,B] where F is the Frobenius form and
        B is the basis change so that $M=B^{-1} F B$.

        EXAMPLES:
            sage: a = pari('[1,2;3,4]')
            sage: a.matfrobenius()
            [0, 2; 1, 5]
            sage: a.matfrobenius(flag=1)
            [x^2 - 5*x - 2]
            sage: a.matfrobenius(2)
            [[0, 2; 1, 5], [1, -1/3; 0, 1/3]]
            sage: v = a.matfrobenius(2)
            sage: v[0]
            [0, 2; 1, 5]
            sage: v[1]^(-1)*v[0]*v[1]
            [1, 2; 3, 4]

        We let t be the matrix of $T_2$ acting on modular symbols of level 43,
        which was computed using \code{ModularSymbols(43,sign=1).T(2).matrix()}:

            sage: t = pari('[3, -2, 0, 0; 0, -2, 0, 1; 0, -1, -2, 2; 0, -2, 0, 2]')
            sage: t.matfrobenius()
            [0, 0, 0, -12; 1, 0, 0, -2; 0, 1, 0, 8; 0, 0, 1, 1]
            sage: t.charpoly('x')
            x^4 - x^3 - 8*x^2 + 2*x + 12
            sage: t.matfrobenius(1)
            [x^4 - x^3 - 8*x^2 + 2*x + 12]

        AUTHOR:
           -- 2006-04-02: Martin Albrecht
        """
        _sig_on
        return self.new_gen(matfrobenius(self.g, flag, 0))


    ###########################################
    # 9: SUMS, products, integrals and similar functions
    ###########################################


    ###########################################
    # polarit2.c
    ###########################################
    def factor(gen self, limit=-1):
        """
        Return the factorization of x.

        lim is optional and can be set whenever x is of (possibly
        recursive) rational type. If lim is set return partial
        factorization, using primes up to lim (up to primelimit if
        lim=0).

        EXAMPLES:
            sage: pari('x^10-1').factor()
            [x - 1, 1; x + 1, 1; x^4 - x^3 + x^2 - x + 1, 1; x^4 + x^3 + x^2 + x + 1, 1]
            sage: pari(2^100-1).factor()
            [3, 1; 5, 3; 11, 1; 31, 1; 41, 1; 101, 1; 251, 1; 601, 1; 1801, 1; 4051, 1; 8101, 1; 268501, 1]

        PARI doesn't have an algorithm for factoring multivariate polynomials:

            sage: pari('x^3 - y^3').factor()
            Traceback (most recent call last):
            ...
            PariError: sorry, (15)
        """
        _sig_on
        return P.new_gen(factor0(self.g, limit))


    ###########################################
    # misc (classify when I know where they go)
    ###########################################

    def hilbert(x, y, p):
        t0GEN(y)
        t1GEN(p)
        _sig_on
        return hil0(x.g, t0, t1)

    def chinese(self, y):
        t0GEN(y)
        _sig_on
        return P.new_gen(chinese(self.g, t0))

    def order(self):
        _sig_on
        return P.new_gen(order(self.g))

    def znprimroot(self):
        _sig_on
        return P.new_gen(ggener(self.g))

    def __abs__(self):
        return self.abs()

    def norm(gen self):
        _sig_on
        return P.new_gen(gnorm(self.g))

    def nextprime(gen self):
        """
        nextprime(x): smallest pseudoprime >= x
        """
        #NOTE: This is much faster than MAGMA's NextPrime with Proof := False.
        _sig_on
        return P.new_gen(gnextprime(self.g))

    def subst(self, var, y):
        """
        EXAMPLES:
           sage: x = pari("x"); y = pari("y")
           sage: f = pari('x^3 + 17*x + 3')
           sage: f.subst(x, y)
           y^3 + 17*y + 3
           sage: f.subst(x, "z")
           z^3 + 17*z + 3
           sage: f.subst(x, "z")^2
           z^6 + 34*z^4 + 6*z^3 + 289*z^2 + 102*z + 9
           sage: f.subst(x, "x+1")
           x^3 + 3*x^2 + 20*x + 21
           sage: f.subst(x, "xyz")
           xyz^3 + 17*xyz + 3
           sage: f.subst(x, "xyz")^2
           xyz^6 + 34*xyz^4 + 6*xyz^3 + 289*xyz^2 + 102*xyz + 9
        """
        cdef long n
        n = P.get_var(var)
        t0GEN(y)
        _sig_on
        return P.new_gen(gsubst(self.g, n, t0))

    def substpol(self, y, z):
        t0GEN(y)
        t1GEN(z)
        _sig_on
        return self.new_gen(gsubstpol(self.g, t0, t1))

    def taylor(self, v=-1):
        _sig_on
        return self.new_gen(tayl(self.g, self.get_var(v), precdl))

    def thue(self, rhs, ne):
        t0GEN(rhs)
        t1GEN(ne)
        _sig_on
        return self.new_gen(thue(self.g, t0, t1))

    def charpoly(self, var=-1, flag=0):
        """
        charpoly(A,{v=x},{flag=0}): det(v*Id-A) = characteristic
        polynomial of A using the comatrix. flag is optional and may
        be set to 1 (use Lagrange interpolation) or 2 (use Hessenberg
        form), 0 being the default.
        """
        _sig_on
        return P.new_gen(charpoly0(self.g, P.get_var(var), flag))


    def kronecker(gen self, y):
        t0GEN(y)
        _sig_on
        return P.new_gen(gkronecker(self.g, t0))


    def type(gen self):
        return str(type_name(typ(self.g)))


    def polinterpolate(self, ya, x):
        """
        self.polinterpolate({ya},{x},{&e}): polynomial interpolation
        at x according to data vectors self, ya (i.e. return P such
        that P(self[i]) = ya[i] for all i).  Also return an error
        estimate on the returned value.
        """
        t0GEN(ya)
        t1GEN(x)
        cdef GEN dy, g
        _sig_on
        g = polint(self.g, t0, t1, &dy)
        _sig_off
        dif = self.new_gen_noclear(dy)
        return self.new_gen(g), dif

    def algdep(self, long n, bit=0):
        """
        EXAMPLES:
            sage: n = pari.set_real_precision (200)
            sage: w1 = pari('z1=2-sqrt(26); (z1+I)/(z1-I)')
            sage: f = w1.algdep(12); f
            545*x^11 - 297*x^10 - 281*x^9 + 48*x^8 - 168*x^7 + 690*x^6 - 168*x^5 + 48*x^4 - 281*x^3 - 297*x^2 + 545*x
            sage: f(w1)
            7.75513996 E-200 + 5.70672991 E-200*I     # 32-bit
            3.780069700150794274 E-209 - 9.362977321012524836 E-211*I   # 64-bit
            sage: f.factor()
            [x, 1; x + 1, 2; x^2 + 1, 1; x^2 + x + 1, 1; 545*x^4 - 1932*x^3 + 2790*x^2 - 1932*x + 545, 1]
            sage: pari.set_real_precision(n)
            200
        """
        cdef long b
        b = bit
        _sig_on
        return self.new_gen(algdep0(self.g, n, b, prec))

    def concat(self, y):
        t0GEN(y)
        _sig_on
        return self.new_gen(concat(self.g, t0))

    def lindep(self, flag=0):
        _sig_on
        return self.new_gen(lindep0(self.g, flag, prec))

    def listinsert(self, obj, long n):
        t0GEN(obj)
        _sig_on
        return self.new_gen(listinsert(self.g, t0, n))

    def listput(self, obj, long n):
        t0GEN(obj)
        _sig_on
        return self.new_gen(listput(self.g, t0, n))



    def elleisnum(self, long k, int flag=0):
        """
        om.elleisnum(k, {flag=0}, {prec}): om=[om1,om2] being a
            2-component vector giving a basis of a lattice L and k an
            even positive integer, computes the numerical value of the
            Eisenstein series of weight k. When flag is non-zero and
            k=4 or 6, this gives g2 or g3 with the correct
            normalization.

        INPUT:
            om -- gen, 2-component vector giving a basis of a lattice L
            k  -- int (even positive)
            flag -- int (default 0)
            pref -- precision

        OUTPUT:
            gen -- numerical value of E_k

        EXAMPLES:
            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: om = e.omega()
            sage: om
            [2.490212560855055075321357792, 1.971737701551648204422407698*I]                       # 32-bit
            [2.4902125608550550753213577919423024602, 1.9717377015516482044224076981513423349*I]   # 64-bit
            sage: om.elleisnum(2)
            -5.288649339654257621791534695              # 32-bit
            -5.2886493396542576217915346952045657616    # 64-bit
            sage: om.elleisnum(4)
            112.0000000000000000000000000               # 32-bit
            112.00000000000000000000000000000000000     # 64-bit
            sage: om.elleisnum(100)
            2.153142485760776361127070349 E50              # 32-bit
            2.1531424857607763611270703492586704424 E50    # 64-bit
        """
        _sig_on
        return self.new_gen(elleisnum(self.g, k, flag, prec))

    def ellwp(self, z='z', long n=20, long flag=0):
        """
        ellwp(E, z,{flag=0}): Return the complex value of the Weierstrass
        P-function at z on the lattice defined by e.

        INPUT:
            E -- list OR elliptic curve
                  list -- [om1, om2], which are Z-generators for a lattice
                  elliptic curve -- created using ellinit

            z -- (optional) complex number  OR string (default = "z")
                   complex number -- any number in the complex plane
                   string (or PARI variable) -- name of a variable.

            n -- int (optional: default 20) if z is a variable, compute up to at least o(z^n).

            flag -- int: 0 (default): compute only P(z)
                         1 compute [P(z),P'(z)]
                         2 consider om or E as an elliptic curve and use P-function to
                           compute the point on E (with the Weierstrass equation for E)
                           P(z)
                           for that curve (identical to ellztopoint in this case).

        OUTPUT:
            gen -- complex number or list of two complex numbers

        EXAMPLES:

        We first define the elliptic curve X_0(11):
            sage: E = pari([0,-1,1,-10,-20]).ellinit()

        Compute P(1).
            sage: E.ellwp(1)
            13.96586952574849779469497770 + 1.465441833 E-249*I                       # 32-bit
            13.965869525748497794694977695009324221 + 5.607702583566084181 E-693*I    # 64-bit

        Compute P(1+I), where I = sqrt(-1).
            sage: E.ellwp(pari("1+I"))
            -1.115106825655550879209066492 + 2.334190523074698836184798794*I                       # 32-bit
            -1.1151068256555508792090664916218413986 + 2.3341905230746988361847987936140321964*I   # 64-bit
            sage: E.ellwp("1+I")
            -1.115106825655550879209066492 + 2.334190523074698836184798794*I                       # 32-bit
            -1.1151068256555508792090664916218413986 + 2.3341905230746988361847987936140321964*I   # 64-bit

        The series expansion, to the default 20 precision:
            sage: E.ellwp()
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)

        Compute the series for wp to lower precision:
            sage: E.ellwp(n=4)
            z^-2 + 31/15*z^2 + O(z^4)

        Next we use the version where the input is generators for a lattice:
            sage: pari(['1.2692', '0.63 + 1.45*I']).ellwp(1)
            13.96561469366894364802003736 + 0.0006448292728105361474541633318*I                        # 32-bit
            13.965614693668943648020037358850990554 + 0.00064482927281053614745416280316868200698*I    # 64-bit

        With flag 1 compute the pair P(z) and P'(z):
            sage: E.ellwp(1, flag=1)
            [13.96586952574849779469497770 + 1.465441833 E-249*I, 50.56193008800727525558465690 + 4.46944479 E-249*I]    # 32-bit
            [13.965869525748497794694977695009324221 + 5.607702583566084181 E-693*I, 50.561930088007275255584656898892400699 + 1.7102908181111172423 E-692*I]   # 64-bit

        With flag=2, the computed pair is (x,y) on the curve instead of [P(z),P'(z)]:
            sage: E.ellwp(1, flag=2)
            [14.29920285908183112802831103 + 1.465441833 E-249*I, 50.06193008800727525558465690 + 4.46944479 E-249*I]    # 32-bit
            [14.299202859081831128028311028342657555 + 5.607702583566084181 E-693*I, 50.061930088007275255584656898892400699 + 1.7102908181111172423 E-692*I]   # 64-bit
        """
        t0GEN(z)
        if n < 0:
            n = 0
        if n%2 == 1:
            n = n + 1
        _sig_on
        return self.new_gen(ellwp0(self.g, t0, flag, prec, (n+2)/2))

    def ellchangepoint(self, y):
        """
        self.ellchangepoint(y): change data on point or vector of points self
                             on an elliptic curve according to y=[u,r,s,t]

        EXAMPLES:
            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: x = pari([1,0,1])
            sage: e.ellisoncurve([1,4,4])
            False
            sage: e.ellisoncurve(x)
            True
            sage: f = e.ellchangecurve([1,2,3,-1])
            sage: f[:5]   # show only first five entries
            [6, -2, -1, 17, 8]
            sage: x.ellchangepoint([1,2,3,-1])
            [-1, 4]
            sage: f.ellisoncurve([-1,4])
            True
        """
        t0GEN(y)
        _sig_on
        return self.new_gen(pointch(self.g, t0))

    ##################################################
    # Technical functions that can be used by other
    # classes that derive from gen.
    ##################################################
    cdef gen pari(self, object x):
        return pari(x)

    cdef gen new_gen(self, GEN x):
        return P.new_gen(x)

    cdef gen new_gen_noclear(self, GEN x):
        return P.new_gen_noclear(x)

    cdef GEN _deepcopy_to_python_heap(self, GEN x, pari_sp* address):
        return P.deepcopy_to_python_heap(x, address)

    cdef int get_var(self, v):
        return P.get_var(v)



cdef unsigned long num_primes

cdef class PariInstance(sage.structure.parent_base.ParentWithBase):
    def __init__(self, long size=16000000, unsigned long maxprime=500000):
        """
        Initialize the PARI system.

        INPUT:
            size -- long, the number of bytes for the initial PARI stack
                    (see note below)
            maxprime -- unsigned long, upper limit on a precomputed prime
                        number table  (default: 500000)

        NOTES:

            * In py_pari, the PARI stack is different than in gp or the
              PARI C library.  In Python, instead of the PARI stack
              holding the results of all computations, it *only* holds the
              results of an individual computation.  Each time a new
              Python/PARI object is computed, it it copied to its own
              space in the Python heap, and the memory it occupied on the
              PARI stack is freed.  Thus it is not necessary to make the
              stack very large.  Also, unlike in PARI, if the stack does
              overflow, in most cases the PARI stack is automatically
              increased and the relevant step of the computation rerun.

              This design obviously involves some performance penalties
              over the way PARI works, but it scales much better and is
              far more robus for large projects.

            * If you do not want prime numbers, put maxprime=2, but be
              careful because many PARI functions require this table.  If
              you get the error message "not enough precomputed primes",
              increase this parameter.

        """
        if bot:
            return  # pari already initialized.

        global initialized, num_primes, ZERO, ONE, TWO, avma, top, bot

        #print "Initializing PARI (size=%s, maxprime=%s)"%(size,maxprime)
        pari_init(1024, maxprime)

        _sig_on
        init_stack(size)
        _sig_off

        GP_DATA.fmt.prettyp = 0

        # Take control of SIGSEGV back from PARI.
        import signal
        signal.signal(signal.SIGSEGV, signal.SIG_DFL)

        # We do the following, since otherwise the IPython pager
        # causes sage to crash when it is exited early.  This is
        # because when PARI was initialized it set a trap for this
        # signal.
        import signal
        signal.signal(signal.SIGPIPE, _my_sigpipe)
        initialized = 1
        stack_avma = avma
        num_primes = maxprime
        self.ZERO = self(0)    # todo: gen_0
        self.ONE = self(1)
        self.TWO = self(2)

    def __dealloc__(self):
        # TODO -- add pari free here
        pass

    def __repr__(self):
        return "Interface to the PARI C library"

    def __hash__(self):
        return 907629390   # hash('pari')

    cdef has_coerce_map_from_c_impl(self, x):
        return True

    def __richcmp__(left, right, int op):
        """
        EXAMPLES:
            sage: pari == pari
            True
            sage: pari == gp
            False
            sage: pari == 5
            False
        """
        return (<Parent>left)._richcmp(right, op)

    def default(self, variable, value=None):
        if not value is None:
            return self('default(%s, %s)'%(variable, value))
        return self('default(%s)'%variable)

    def set_debug_level(self, level):
        """
        Set the debug PARI C library variable.
        """
        self.default('debug', int(level))

    def get_debug_level(self):
        """
        Set the debug PARI C library variable.
        """
        return int(self.default('debug'))

    cdef GEN toGEN(self, x) except NULL:
        cdef gen _x
        if isinstance(x, gen):
            _x = x
            return _x.g
        s = str(x)
        cdef GEN g
        _sig_on
        g = flisseq(s)
        _sig_off
        return g


    def set_real_precision(self, long n):
        """
        Sets the PARI default real precision, both for creation of
        new objects and for printing.  This is the number of digits
        *IN DECIMAL* that real numbers are printed or computed to
        by default.

        Returns the previous PARI real precision.
        """
        cdef unsigned long k
        global prec

        k = GP_DATA.fmt.sigd
        s = str(n)
        _sig_on
        sd_realprecision(s, 2)
        prec = GP_DATA.fmt.sigd
        _sig_off
        return int(k)  # Python int

    def get_real_precision(self):
        return GP_DATA.fmt.sigd

    def set_series_precision(self, long n):
        global precdl
        precdl = n

    def get_series_precision(self):
        return precdl


    ###########################################
    # Create a gen from a GEN object.
    # This *steals* a reference to the GEN, as it
    # frees the memory the GEN occupied.
    ###########################################

    cdef gen new_gen(self, GEN x):
        """
        Create a new gen, then free the *entire* stack and call _sig_off.
        """
        cdef gen g
        g = _new_gen(x)

        global mytop, avma
        avma = mytop
        _sig_off
        return g

    cdef void clear_stack(self):
        """
        Clear the entire PARI stack and turn off call _sig_off.
        """
        global mytop, avma
        avma = mytop
        _sig_off

    cdef gen new_gen_noclear(self, GEN x):
        """
        Create a new gen, but don't free any memory on the stack and
        don't call _sig_off.
        """
        return _new_gen(x)

    def double_to_gen(self, x):
        cdef double dx
        dx = float(x)
        return self.double_to_gen_c(dx)

    cdef gen double_to_gen_c(self, double x):
        """
        Create a new gen with the value of the double x, using Pari's
        dbltor.

        EXAMPLES:
            sage: pari.double_to_gen(1)
            1.0000000000000000000
            sage: pari.double_to_gen(1e30)
            1.0000000000000000199 E30
            sage: pari.double_to_gen(0)
            0.E-15
            sage: pari.double_to_gen(-sqrt(RDF(2)))
            -1.4142135623730951455
        """
        # Pari has an odd concept where it attempts to track the accuracy
        # of floating-point 0; a floating-point zero might be 0.0e-20
        # (meaning roughly that it might represent any number in the
        # range -1e-20 <= x <= 1e20).

        # Pari's dbltor converts a floating-point 0 into the Pari real
        # 0.0e-307; Pari treats this as an extremely precise 0.  This
        # can cause problems; for instance, the Pari incgam() function can
        # be very slow if the first argument is very precise.

        # So we translate 0 into a floating-point 0 with 53 bits
        # of precision (that's the number of mantissa bits in an IEEE
        # double).

        if x == 0:
            return self.new_gen(real_0_bit(-53))
        else:
            return self.new_gen(dbltor(x))

    cdef GEN double_to_GEN(self, double x):
        if x == 0:
            return real_0_bit(-53)
        else:
            return dbltor(x)

    def complex(self, re, im):
        """
        Create a new complex number, initialized from re and im.
        """
        t0GEN(re)
        t1GEN(im)
        cdef GEN cp
        cp = cgetg(3, t_COMPLEX)
        __set_lvalue__(gel(cp, 1), t0)
        __set_lvalue__(gel(cp, 2), t1)
        return self.new_gen(cp)

    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address):
        return deepcopy_to_python_heap(x, address)

    cdef gen new_ref(self, GEN x, g):
        cdef gen p
        p = gen()
        p.init(x, 0)
        try:
            p.refers_to[-1] = g  # so underlying memory won't get deleted
                                 # out from under us.
        except TypeError:
            p.refers_to = {-1:g}
        return p

    cdef gen adapt(self, s):
        if isinstance(s, gen):
            return s
        return pari(s)

    def __call__(self, s):
        """
        Create the PARI object got by evaluating s using PARI.
        """
        cdef pari_sp prev_avma
        global avma

        prev_avma = avma

        if isinstance(s, gen):
            return s
        try:
            return s._pari_()
        except AttributeError:
            pass
        if isinstance(s, (types.ListType, types.XRangeType,
                            types.TupleType, xsrange)):
            v = self.vector(len(s))
            for i, x in enumerate(s):
                v[i] = self(x)
            return v
        elif isinstance(s, bool):
            if s:
                return self.ONE
            return self.ZERO

        cdef GEN g
        t = str(s)
        _sig_str('evaluating PARI string')
        g = gp_read_str(t)
        _sig_off
        return self.new_gen(g)

    cdef _coerce_c_impl(self, x):
        """
        Implicit canonical coercion into a PARI object.
        """
        try:
            return self(x)
        except (TypeError, AttributeError):
            raise TypeError, "no canonical coercion of %s into PARI"%x
        if isinstance(x, gen):
            return x
        raise TypeError, "x must be a PARI object"

    cdef _an_element_c_impl(self):  # override this in SageX
        return self.ZERO

    def new_with_prec(self, s, long precision=0):
        r"""
        pari.new_with_bits_prec(self, s, precision) creates s as a PARI gen
        with at lest precision decimal \emph{digits} of precision.
        """
        global prec
        cdef unsigned long old_prec
        old_prec = prec
        if not precision:
            precision = prec
        self.set_real_precision(precision)
        x = self(s)
        self.set_real_precision(old_prec)
        return x

    def new_with_bits_prec(self, s, long precision=0):
        r"""
        pari.new_with_bits_prec(self, s, precision) creates s as a PARI gen
        with (at most) precision \emph{bits} of precision.
        """
        global prec

        cdef unsigned long old_prec
        old_prec = prec
        precision = long(precision / 3.4)+1     # be safe, since log_2(10) = 3.3219280948873626
        if not precision:
            precision = prec
        self.set_real_precision(precision)
        x = self(s)
        self.set_real_precision(old_prec)
        return x



    cdef int get_var(self, v):
        """
        Converts a Python string into a PARI variable reference number.
        Or if v = -1, returns -1.
        """
        if v != -1:
            s = str(v)
            return fetch_user_var(s)
        return -1


    ############################################################
    # conversion between GEN and string types
    # Note that string --> GEN evaluates the string in PARI,
    # where GEN_to_str returns a Python string.
    ############################################################
    cdef object GEN_to_str(self, GEN g):
        cdef char* c
        cdef int n
        _sig_on
        c = GENtostr(g)
        _sig_off
        s = str(c)
        free(c)
        return s


    ############################################################
    # Initialization
    ############################################################

    def allocatemem(self, silent=False):
        r"""
        Double the \emph{PARI} stack.
        """
        if not silent:
            print "Doubling the PARI stack."
        init_stack(0)

    def pari_version(self):
        return str(PARIVERSION)

    def init_primes(self, _M):
        """
        Recompute the primes table including at least all primes up to
        M (but possibly more).

        EXAMPLES:
            sage: pari.init_primes(200000)
        """
        cdef unsigned long M
        M = _M
        global diffptr, num_primes
        if M <= num_primes:
            return
        #if not silent:
        #    print "Extending PARI prime table up to %s"%M
        free(<void*> diffptr)
        num_primes = M
        diffptr = initprimes(M)


    ##############################################
    ## Support for GP Scripts
    ##############################################


    def read(self, filename):
        r"""
        Read a script from the named filename into the interpreter, where
        s is a string.  The functions defined in the script are then
        available for use from SAGE/PARI.

        EXAMPLE:

            If foo.gp is a script that contains
            \begin{verbatim}
                {foo(n) =
                    n^2
                }
            \end{verbatim}
            and you type \code{read("foo.gp")}, then the command
            \code{pari("foo(12)")} will create the Python/PARI gen which
            is the integer 144.

        CONSTRAINTS:
            The PARI script must *not* contain the following function calls:

                 print, default, ???    (please report any others that cause trouble)
        """
        F = open(filename).read()
        while True:
            i = F.find("}")
            if i == -1:
                _read_script(F)
                break
            s = F[:i+1]
            _read_script(s)
            F = F[i+1:]
        return

    ##############################################

    def prime_list(self, long n):
        """
        prime_list(n): returns list of the first n primes

        To extend the table of primes use pari.init_primes(M).

        INPUT:
            n -- C long
        OUTPUT:
            gen -- PARI list of first n primes

        EXAMPLES:
            sage: pari.prime_list(0)
            []
            sage: pari.prime_list(-1)
            []
            sage: pari.prime_list(3)
            [2, 3, 5]
            sage: pari.prime_list(10)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
            sage: pari.prime_list(20)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
            sage: len(pari.prime_list(1000))
            1000
        """
        if n >= 2:
            self.nth_prime(n)
        _sig_on
        return self.new_gen(primes(n))

    def primes_up_to_n(self, long n):
        """
        Return the primes <= n as a pari list.
        """
        if n <= 1:
            return pari([])
        self.init_primes(n+1)
        return self.prime_list(pari(n).primepi())

##         cdef long k
##         k = (n+10)/math.log(n)
##         p = 2
##         while p <= n:
##             p = self.nth_prime(k)
##             k = 2
##         v = self.prime_list(k)
##         return v[:pari(n).primepi()]

    def __nth_prime(self, long n):
        """
        nth_prime(n): returns the n-th prime, where n is a C-int
        """
        global num_primes

        if n <= 0:
            raise ArithmeticError, "nth prime meaningless for negative n (=%s)"%n
        cdef GEN g
        _sig_on
        g = prime(n)
        return self.new_gen(g)


    def nth_prime(self, long n):
        try:
            return self.__nth_prime(n)
        except PariError:
            self.init_primes(max(2*num_primes,20*n))
            return self.nth_prime(n)

    def euler(self):
        """
        Return Euler's constant to the current real precision.

        EXAMPLES:
            sage: pari.euler()
            0.5772156649015328606065120901             # 32-bit
            0.57721566490153286060651209008240243104   # 64-bit
            sage: orig = pari.get_real_precision(); orig
            28         # 32-bit
            38         # 64-bit
            sage: _ = pari.set_real_precision(50)
            sage: pari.euler()
            0.57721566490153286060651209008240243104215933593992

        We restore precision to the default.
            sage: pari.set_real_precision(orig)
            50

        Euler is returned to the original precision:
            sage: pari.euler()
            0.5772156649015328606065120901              # 32-bit
            0.57721566490153286060651209008240243104    # 64-bit
        """
        if not precision:
            precision = prec
        _sig_on
        consteuler(precision)
        return self.new_gen(geuler)

    def pi(self):
        """
        Return the value of the constant pi = 3.1415 to the current
        real precision.

        EXAMPLES:
            sage: pari.pi()
            3.141592653589793238462643383             # 32-bit
            3.1415926535897932384626433832795028842   # 64-bit
            sage: orig_prec = pari.set_real_precision(5)
            sage: pari.pi()
            3.1416

            sage: pari.get_real_precision()
            5
            sage: _ = pari.set_real_precision(40)
            sage: pari.pi()
            3.141592653589793238462643383279502884197


        We restore precision to the default.
            sage: _ = pari.set_real_precision(orig_prec)

        Note that pi is now computed to the original precision:

            sage: pari.pi()
            3.141592653589793238462643383              # 32-bit
            3.1415926535897932384626433832795028842    # 64-bit
        """
        if not precision:
            precision = prec
        _sig_on
        constpi(precision)
        return self.new_gen(gpi)

    def pollegendre(self, long n, v=-1):
        """
        pollegendre(n,{v=x}): legendre polynomial of degree n (n
        C-integer), in variable v
        """
        _sig_on
        return self.new_gen(legendre(n, self.get_var(v)))

    def poltchebi(self, long n, v=-1):
        _sig_on
        return self.new_gen(tchebi(n, self.get_var(v)))

    def factorial(self, long n):
        """
        Return the factorial of the integer n as a PARI gen.
        """
        _sig_on
        return self.new_gen(mpfact(n))

    def polcyclo(self, long n, v=-1):
        _sig_on
        return self.new_gen(cyclo(n, self.get_var(v)))

    def polsubcyclo(self, long n, long d, v=-1):
        _sig_on
        return self.new_gen(polsubcyclo(n, d, self.get_var(v)))

    def polzagier(self, long n, long m):
        _sig_on
        return self.new_gen(polzag(n, m))

    def listcreate(self, long n):
        _sig_on
        return self.new_gen(listcreate(n))

    def vector(self, long n, entries=None):
        """
        vector(long n, entries=None):
        Create and return the length n PARI vector with given list of entries.
        """
        cdef gen v
        _sig_on
        v = self.new_gen(zerovec(n))
        if entries != None:
            if len(entries) != n:
                raise IndexError, "length of entries (=%s) must equal n (=%s)"%\
                      (len(entries), n)
            for i, x in enumerate(entries):
                v[i] = x
        return v

    def matrix(self, long m, long n, entries=None):
        """
        matrix(long m, long n, entries=None):
        Create and return the m x n PARI matrix with given list of entries.
        """
        cdef int i, j, k
        cdef gen A
        _sig_on
        A = self.new_gen(zeromat(m,n)).Mat()
        if entries != None:
            if len(entries) != m*n:
                raise IndexError, "len of entries (=%s) must be %s*%s=%s"%(len(entries),m,n,m*n)
            k = 0
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    A[i,j] = entries[k]
                    k = k + 1
        return A

    def finitefield_init(self, p, long n, var=-1):
        """
        finitefield_init(p, long n, var="x"): Return a polynomial f(x) so
        that the extension of F_p of degree n is k = F_p[x]/(f).  Moreover,
        the element x (mod f) of k is a generator for the multiplicative
        group of k.

        INPUT:
            p -- int, a prime number
            n -- int, positive integer
            var -- str, default to "x", but could be any pari variable.
        OUTPUT:
            pari polynomial mod p -- defines field
        EXAMPLES:
            sage: pari.finitefield_init(97,1)
            Mod(1, 97)*x + Mod(92, 97)

        The last entry in each of the following two lines is
        determined by a random algorithm.
            sage: pari.finitefield_init(7,2)       # random
            Mod(1, 7)*x^2 + Mod(6, 7)*x + Mod(3, 7)
            sage: pari.finitefield_init(2,3)       # random
            Mod(1, 2)*x^3 + Mod(1, 2)*x^2 + Mod(1, 2)
        """
        cdef gen _p, _f2, s
        cdef int err
        cdef long x
        cdef GEN v, g
        _p = self(int(p))
        if n < 1:
            raise ArithmeticError, "Degree n (=%s) must be at least 1."%n
        if _p < 2 or not _p.isprime():
            raise ArithmeticError, "Prime p (=%s) must be prime."%_p
        x = self.get_var(var)
        if n == 1:
            _sig_on
            return self.new_gen(ffinit(_p.g, n, x)) - _p.znprimroot()
        _sig_on
        f = self.new_gen(ffinit(_p.g, n, x))
        _f2 = f.lift()
        _sig_on
        g = FpXQ_gener(_f2.g, _p.g)
        s = self.new_gen(g)*self.ONE.Mod(p)
        return s.Mod(f).charpoly(var)


    ##############################################


cdef int init_stack(size_t size) except -1:
    cdef size_t s

    global top, bot, avma, stack_avma, mytop

    # delete this if get core dumps and change the 2* to a 1* below.
    if bot:
        sage_free(<void*>bot)

    if size == 0:
        size = 2*(top-bot)

    # if size == -1, then allocate the biggest chunk possible
    if size == -1:
        s = 4294967295
        while True:
            s = fix_size(s)
            bot = <pari_sp> sage_malloc(s)
            if bot:
                break
            s = s/2
    else:
        # Decide on size
        s = fix_size(size)
        # Alocate memory for new stack using Python's memory allocator.
        # As explained in the python/C api reference, using this instead
        # of malloc is much better (and more platform independent, etc.)
        bot = <pari_sp> sage_malloc(s)
        if not bot:
            raise MemoryError, "Unable to allocate %s bytes memory for PARI."%(<long>size)
    #endif
    top = bot + s
    mytop = top
    avma = top
    stack_avma = avma


def _my_sigpipe(signum, frame):
    # If I do anything, it messes up Ipython's pager.
    pass

cdef size_t fix_size(size_t a):
    cdef size_t b
    b = a - (a & (BYTES_IN_LONG-1))     # sizeof(long) | b <= a
    if b < 1024:
        b = 1024
    return b

cdef GEN deepcopy_to_python_heap(GEN x, pari_sp* address):
    cdef size_t s
    cdef pari_sp tmp_bot, tmp_top, tmp_avma
    global avma, bot, top, mytop
    cdef GEN h

    tmp_top = top
    tmp_bot = bot
    tmp_avma = avma

    h = forcecopy(x)
    s = <size_t> (tmp_avma - avma)

    #print "Allocating %s bytes for PARI/Python object"%(<long> s)
    bot = <pari_sp> sage_malloc(s)
    top = bot + s
    avma = top
    h = forcecopy(x)
    address[0] = bot

    # Restore the stack to how it was before x was created.
    top = tmp_top
    bot = tmp_bot
    avma = tmp_avma
    return h

# The first one makes a separate copy on the heap, so the stack
# won't overflow -- but this costs more time...
cdef gen _new_gen(GEN x):
    cdef GEN h
    cdef pari_sp address
    cdef gen y
    _sig_on
    h = deepcopy_to_python_heap(x, &address)
    y = PY_NEW(gen)
    y.init(h, address)
    _sig_off
    return y

cdef gen xxx_new_gen(GEN x):
    cdef gen y
    y = PY_NEW(gen)
    y.init(x, 0)
    _sig_off
    return y

def min(x,y):
    """
    min(x,y): Return the minimum of x and y.
    """
    if x <= y:
        return x
    return y

def max(x,y):
    """
    max(x,y): Return the maximum of x and y.
    """
    if x >= y:
        return x
    return y

cdef int _read_script(char* s) except -1:
    _sig_on
    gp_read_str(s)
    _sig_off

    # We set top to avma, so the script's affects won't be undone
    # when we call new_gen again.
    global mytop, top, avma
    mytop = avma
    return 0


#######################
# Base gen class
#######################


cdef extern from "pari/pari.h":
    char *errmessage[]
    int user
    int errpile
    int noer

def __errmessage(d):
    if d <= 0 or d > noer:
        return "unknown"
    return errmessage[d]

# FIXME: we derive PariError from RuntimeError, for backward
# compatibility with code that catches the latter. Once this is
# in production, we should change the base class to StandardError.
from exceptions import RuntimeError

# can we have "cdef class" ?
# because of the inheritance, need to somehow "import" the builtin
# exception class...
class PariError (RuntimeError):

    errmessage = staticmethod(__errmessage)

    def __repr__(self):
        return "PariError(%d)"%self.args[0]

    def __str__(self):
        return "%s (%d)"%(self.errmessage(self.args[0]),self.args[0])

# We expose a trap function to C.
# If this function returns without raising an exception,
# the code is retried.
# This is a proof-of-concept example.
# THE TRY CODE IS NOT REENTRANT -- NO CALLS TO PARI FROM HERE !!!
#              - Gonzalo Tornario

cdef void _pari_trap "_pari_trap" (long errno, long retries) except *:
    """
    TESTS:
        sage: v = pari.listcreate(10^9)
        Traceback (most recent call last):
        ...
        RuntimeError: The PARI stack overflowed.  It has automatically been doubled using pari.allocatemem().  Please retry your computation, possibly after you manually call pari.allocatemem() a few times.
    """
    if retries > 100:
        raise RuntimeError, "_pari_trap recursion too deep"
    if errno == errpile:
        P.allocatemem()
        raise RuntimeError, "The PARI stack overflowed.  It has automatically been doubled using pari.allocatemem().  Please retry your computation, possibly after you manually call pari.allocatemem() a few times."

        #print "Stack overflow! (%d retries so far)"%retries
        #print " enlarge the stack."
        P.allocatemem(silent=True)
    elif errno == user:
        raise Exception, "PARI user exception"
    else:
        raise PariError, errno



def vecsmall_to_intlist(gen v):
    """
    INPUT:
        v -- a gen of type Vecsmall
    OUTPUT:
        a Python list of Python ints
    """
    if typ(v.g) != t_VECSMALL:
        raise TypeError, "input v must be of type vecsmall (use v.Vecsmall())"
    return [v.g[k+1] for k in range(glength(v.g))]
