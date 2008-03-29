"""
Multivariate polynomials via libSINGULAR.

AUTHORS:
    -- Martin Albrecht <malb@informatik.uni-bremen.de>

TODO:
   -- implement Real, Complex, NumberFields

TESTS:
    sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
    sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = MPolynomialRing_libsingular(GF(2^8,'a'),3)
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = MPolynomialRing_libsingular(GF(127),3)
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True

"""

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"

import os
import sage.rings.memory

# conversion routines
from sage.libs.singular.singular cimport Conversion
cdef Conversion co
co = Conversion()

# polynomial imports
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.polynomial_ring import PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.polydict import ETuple
from sage.rings.polynomial.pbori import BooleanPolynomial


# base ring imports
from sage.rings.rational_field import RationalField
from sage.rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.ring import FiniteField as FiniteField_generic
from sage.rings.finite_field_givaro cimport FiniteField_givaroElement
from sage.rings.finite_field_givaro cimport FiniteField_givaro
from sage.rings.number_field.number_field import NumberField_generic
from sage.rings.rational cimport Rational
from sage.rings.complex_field import is_ComplexField
from sage.rings.real_mpfr import is_RealField
from sage.rings.integer_ring import is_IntegerRing
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer cimport Integer

from sage.structure.parent cimport Parent
from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens

from sage.structure.element cimport EuclideanDomainElement
from sage.structure.element cimport RingElement
from sage.structure.element cimport ModuleElement
from sage.structure.element cimport Element
from sage.structure.element cimport CommutativeRingElement

from sage.structure.factorization import Factorization
from sage.structure.sequence import Sequence

from sage.interfaces.all import macaulay2
from sage.interfaces.singular import singular as singular_default, is_SingularElement, SingularElement
from sage.interfaces.macaulay2 import macaulay2 as macaulay2_default, is_Macaulay2Element

from sage.misc.misc import mul
from sage.misc.sage_eval import sage_eval
from sage.misc.latex import latex

import sage.libs.pari.gen
import polynomial_element


# shared library loading
cdef extern from "dlfcn.h":
    void *dlopen(char *, long)
    char *dlerror()
    void dlclose(void *handle)

cdef init_singular():
    """
    This initializes the Singular library. Right now, this is a hack.

    SINGULAR has a concept of compiled extension modules similar to
    SAGE. For this, the compiled modules need to see the symbols from
    the main programm. However, SINGULAR is a shared library in this
    context these symbols are not known globally. The work around so
    far is to load the library again and to specifiy RTLD_GLOBAL.
    """
    cdef void *handle


    for extension in ["so", "dylib", "dll"]:
        lib = os.environ['SAGE_LOCAL']+"/lib/libsingular."+extension
        if os.path.exists(lib):
            handle = dlopen(lib, 256+1)
            break

    if handle == NULL:
        print dlerror()
        raise ImportError, "cannot load libSINGULAR library"

    # load SINGULAR
    siInit(lib)

    # steal memory manager back or weird things may happen
    sage.rings.memory.pmem_malloc()

    dlclose(handle)

    singular_options[0] = singular_options[0] | Sy_bit(OPT_REDSB)

 # call it
init_singular()


# mapping str --> SINGULAR representation
order_dict = {"dp":ringorder_dp,
              "Dp":ringorder_Dp,
              "lp":ringorder_lp,
              "rp":ringorder_rp,
              "ds":ringorder_ds,
              "Ds":ringorder_Ds,
              "ls":ringorder_ls,
              }

cdef class MPolynomialRing_libsingular(MPolynomialRing_generic):
    """
    A multivariate polynomial ring over QQ or GF(p) implemented using SINGULAR.

    """
    def __init__(self, base_ring, n, names, order='degrevlex'):
        """

        Construct a multivariate polynomial ring subject to the following conditions:

        INPUT:
            base_ring -- base ring (must be either GF(p) (p prime) or QQ)
            n -- number of variables (must be at least 1)
            names -- names of ring variables, may be string of list/tuple
            order -- term order (default: degrevlex)

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: P
            Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: f = 27/113 * x^2 + y*z + 1/2; f
            27/113*x^2 + y*z + 1/2

            sage: P.term_order()
            Degree reverse lexicographic term order

            sage: P = MPolynomialRing(GF(127),3,names='abc', order='lex')
            sage: P
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1; f
            57*a^2*b + 43*c + 1

            sage: P.term_order()
            Lexicographic term order

        """
        cdef char **_names
        cdef char *_name
        cdef int i
        cdef int nblcks
        cdef int offset
        cdef int characteristic
        cdef MPolynomialRing_libsingular k
        cdef MPolynomial_libsingular minpoly
        cdef lnumber *nmp

        is_extension = False

        n = int(n)
        if n<1:
            raise ArithmeticError, "number of variables must be at least 1"

        self.__ngens = n

        order = TermOrder(order, n)

        MPolynomialRing_generic.__init__(self, base_ring, n, names, order)

        self._has_singular = True

        assert(n == len(self._names))

        _names = <char**>omAlloc0(sizeof(char*)*(len(self._names)))

        for i from 0 <= i < n:
            _name = self._names[i]
            _names[i] = omStrDup(_name)

        # from the SINGULAR source code documentation for the rInit function
        ##  characteristic --------------------------------------------------
        ##  input: 0 ch=0 : Q     parameter=NULL    ffChar=FALSE   float_len (done)
        ##         0    1 : Q(a,...)        *names         FALSE             (todo)
        ##         0   -1 : R               NULL           FALSE  0
        ##         0   -1 : R               NULL           FALSE  prec. >6
        ##         0   -1 : C               *names         FALSE  prec. 0..?
        ##         p    p : Fp              NULL           FALSE             (done)
        ##         p   -p : Fp(a)           *names         FALSE             (done)
        ##         q    q : GF(q=p^n)       *names         TRUE              (todo)

        if PY_TYPE_CHECK(base_ring, FiniteField_prime_modn):
            if base_ring.characteristic() <= 2147483629:
                characteristic = base_ring.characteristic()
            else:
                raise TypeError, "p must be <= 2147483629."

        elif PY_TYPE_CHECK(base_ring, RationalField):
            characteristic = 0

        elif PY_TYPE_CHECK(base_ring, FiniteField_generic):
            if base_ring.characteristic() <= 2147483629:
                characteristic = -base_ring.characteristic() # note the negative characteristic
            else:
                raise TypeError, "characteristic must be <= 2147483629."
            k = MPolynomialRing_libsingular(base_ring.prime_subfield(), 1, base_ring.variable_name(), 'lex')
            minpoly = base_ring.polynomial()(k.gen())
            is_extension = True

        elif PY_TYPE_CHECK(base_ring, NumberField_generic):
            raise NotImplementedError, "Number fields are not fully supported yet."
            characteristic = 1
            k = MPolynomialRing_libsingular(RationalField(), 1, base_ring.variable_name(), 'lex')
            minpoly = base_ring.polynomial()(k.gen())
            is_extension = True
        else:
            raise NotImplementedError, "Only GF(q) and QQ are supported."

        self._ring = <ring*>omAlloc0Bin(sip_sring_bin)
        self._ring.ch = characteristic
        self._ring.N = n
        self._ring.names  = _names

        if is_extension:
            rChangeCurrRing(k._ring)
            self._ring.algring = rCopy0(k._ring)
            rComplete(self._ring.algring,1)
            self._ring.P = self._ring.algring.N
            #self._ring.parameter = self._ring.algring.names
            self._ring.parameter = <char**>omAlloc0(sizeof(char*)*2)
            self._ring.parameter[0] = omStrDup(self._ring.algring.names[0])

            nmp = <lnumber*>omAlloc0Bin(rnumber_bin)
            nmp.z= <napoly*>p_Copy(minpoly._poly, self._ring.algring) # fragile?
            nmp.s=2

            self._ring.minpoly=<number*>nmp

        nblcks = len(order.blocks)
        offset = 0

        self._ring.wvhdl  = <int **>omAlloc0((nblcks + 2) * sizeof(int*))
        self._ring.order  = <int *>omAlloc0((nblcks + 2) * sizeof(int *))
        self._ring.block0 = <int *>omAlloc0((nblcks + 2) * sizeof(int *))
        self._ring.block1 = <int *>omAlloc0((nblcks + 2) * sizeof(int *))
        self._ring.OrdSgn = 1


        for i from 0 <= i < nblcks:
            self._ring.order[i] = order_dict.get(order.blocks[i][0], ringorder_lp)
            self._ring.block0[i] = offset + 1
            if order.blocks[i][1] == 0: # may be zero in some cases
                self._ring.block1[i] = offset + n
            else:
                self._ring.block1[i] = offset + order.blocks[i][1]
            offset = self._ring.block1[i]

        self._ring.order[nblcks] = ringorder_C

        rComplete(self._ring, 1)
        self._ring.ShortOut = 0

        rChangeCurrRing(self._ring)
        self._one_element = <MPolynomial_libsingular>co.new_MP(self,p_ISet(1, self._ring))
        self._zero_element = <MPolynomial_libsingular>co.new_MP(self,NULL)

    def __dealloc__(self):
        """
        """
        cdef ring *oldRing = NULL
        if currRing != self._ring:
            oldRing = currRing
            rChangeCurrRing(self._ring)
            rDelete(self._ring)
        else:
            (&currRing)[0] = NULL
            rDelete(self._ring)
        if oldRing != NULL:
            rChangeCurrRing(oldRing)

    cdef _coerce_c_impl(self, element):
        """

        Coerces elements to self.

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]

            We can coerce elements of self to self

            sage: P._coerce_(x*y + 1/2)
            x*y + 1/2

        We can coerce elements for a ring with the same algebraic properties

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P == R
            True

            sage: P is R
            False

            sage: P._coerce_(x*y + 1)
            x*y + 1

            We can coerce base ring elements

            sage: P._coerce_(3/2)
            3/2

            sage: P._coerce_(ZZ(1))
            1

            sage: P._coerce_(int(1))
            1

            sage: k.<a> = GF(2^8)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: P._coerce_(a)
            (a)

        """
        cdef poly *_p
        cdef ring *_ring
        cdef number *_n
        cdef poly *mon
        cdef int i

        _ring = self._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        if PY_TYPE_CHECK(element, MPolynomial_libsingular):
            if element.parent() is <object>self:
                return element
            elif element.parent() == self:
                # is this safe?
                _p = p_Copy((<MPolynomial_libsingular>element)._poly, _ring)
            elif self.base_ring().has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), self.base_ring()))
            else:
                raise TypeError, "parents do not match"

        elif PY_TYPE_CHECK(element, MPolynomial_polydict):
            if element.parent() == self:
                _p = p_ISet(0, _ring)
                for (m,c) in element.element().dict().iteritems():
                    mon = p_Init(_ring)
                    p_SetCoeff(mon, co.sa2si(c, _ring), _ring)
                    for pos in m.nonzero_positions():
                        p_SetExp(mon, pos+1, m[pos], _ring)
                    p_Setm(mon, _ring)
                    _p = p_Add_q(_p, mon, _ring)
            elif self.base_ring().has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), self.base_ring()))
            else:
                raise TypeError, "parents do not match"

        elif isinstance(element, polynomial_element.Polynomial):
            if self.base_ring().has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), self.base_ring()))
            else:
                raise TypeError, "incompatable parents"

        elif PY_TYPE_CHECK(element, CommutativeRingElement):
            # base ring elements
            if  <Parent>element.parent() is self._base:
                # shortcut for GF(p)
                if PY_TYPE_CHECK(self._base, FiniteField_prime_modn):
                    _p = p_ISet(int(element), _ring)
                else:
                    _n = co.sa2si(element,_ring)
                    _p = p_NSet(_n, _ring)

            # also accepting ZZ
            elif element.parent() is IntegerRing():
                if PY_TYPE_CHECK(self._base, RationalField):
                    _n = co.sa2si_ZZ(element,_ring)
                    _p = p_NSet(_n, _ring)
                else: # GF(p)
                    _p = p_ISet(int(element),_ring)
            else:
                # fall back to base ring
                return self._base._coerce_c(element)

        # Accepting int
        elif PY_TYPE_CHECK(element, int):
            _p = p_ISet(int(element), _ring)
        else:
            raise TypeError, "Cannot coerce element"

        return co.new_MP(self,_p)

    def __call__(self, element):
        """
        Construct a new element in self.

        INPUT:
            element -- several types are supported, see below

        EXAMPLE:
            Call supports all conversions _coerce_ supports, plus:

        Coercion from strings:
            sage: P.<x,y,z> = QQ[]
            sage: P('x+y + 1/4')
            x + y + 1/4

        Coercion from SINGULAR elements:
            sage: P._singular_()
            //   characteristic : 0
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_().set_ring()
            sage: P(singular('x + 3/4'))
            x + 3/4

        Coercion from symbolic variables:
            sage: x,y,z = var('x,y,z')
            sage: R = QQ[x,y,z]
            sage: R(x)
            x

        Coercion from 'similar' rings:
            sage: P.<x,y,z> = QQ[]
            sage: R.<a,b,c> = MPolynomialRing(ZZ,3)
            sage: P(a)
            x

        Coercion from PARI objects:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: P(pari('x^2 + y'))
            x^2 + y
            sage: P(pari('x*y'))
            x*y

        Coercion from boolean polynomials:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: P(B.gen(0))
            x

        If everything else fails, we try to coerce to the base ring:
            sage: R.<x,y,z> = GF(3)[]
            sage: R(1/2)
            -1

        TESTS:
        Coerce in a polydict where a coefficient reduces to 0 but isn't 0.
            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S( (5*x*y + x + 17*y)._mpoly_dict_recursive() )
            xx + 2*yy

        Coerce in a polynomial one of whose coefficients reduces to 0.
            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S(5*x*y + x + 17*y)
            xx + 2*yy

        Some other examples that illustrate the same coercion idea:
            sage: R.<x,y> = ZZ[]
            sage: S.<xx,yy> = GF(25,'a')[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

            sage: S.<xx,yy> = Integers(5)[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy
        """
        cdef poly *_p, *mon
        cdef ring *_ring = self._ring
        cdef unsigned int pos
        rChangeCurrRing(_ring)

        # try to coerce first
        try:
            return self._coerce_c_impl(element)
        except TypeError:
            pass

        if PY_TYPE_CHECK(element, SingularElement) or \
           PY_TYPE_CHECK(element, sage.libs.pari.gen.gen):
            element = str(element)

        if PY_TYPE_CHECK(element, basestring):
            # let python do the parsing
            d = self.gens_dict()
            if PY_TYPE_CHECK(self._base, FiniteField_givaro):
                d[str(self._base.gen())]=self._base.gen()
            try:
                if '/' in element:
                    element = sage_eval(element,d)
                else:
                    element = element.replace("^","**")
                    element = eval(element, d, {})
            except SyntaxError:
                raise TypeError

            # we need to do this, to make sure that we actually get an
            # element in self.
            return self._coerce_c(element)

        if PY_TYPE_CHECK(element, MPolynomial_libsingular):
            if element.parent() is not self and element.parent() != self and  element.parent().ngens() == self.ngens():
                # Map the variables in some crazy way (but in order,
                # of course).  This is here since R(blah) is supposed
                # to be "make an element of R if at all possible with
                # no guarantees that this is mathematically solid."
                # TODO: We can do this faster without the dict
                _p = p_ISet(0, _ring)
                K = self.base_ring()
                for (m,c) in element.dict().iteritems():
                    try:
                        c = K(c)
                        if not c: continue
                    except TypeError, msg:
                        p_Delete(&_p, _ring)
                        raise TypeError, msg
                    mon = p_Init(_ring)
                    p_SetCoeff(mon, co.sa2si(c , _ring), _ring)
                    for pos in m.nonzero_positions():
                        p_SetExp(mon, pos+1, m[pos], _ring)
                    p_Setm(mon, _ring)
                    _p = p_Add_q(_p, mon, _ring)
                return co.new_MP(self, _p)

        if PY_TYPE_CHECK(element, MPolynomial_polydict):
            if element.parent().ngens() == self.ngens():
                # Map the variables in some crazy way (but in order,
                # of course).  This is here since R(blah) is supposed
                # to be "make an element of R if at all possible with
                # no guarantees that this is mathematically solid."
                _p = p_ISet(0, _ring)
                K = self.base_ring()
                for (m,c) in element.element().dict().iteritems():
                    try:
                        c = K(c)
                        if not c: continue
                    except TypeError, msg:
                        p_Delete(&_p, _ring)
                        raise TypeError, msg
                    mon = p_Init(_ring)
                    p_SetCoeff(mon, co.sa2si(c , _ring), _ring)
                    for pos in m.nonzero_positions():
                        p_SetExp(mon, pos+1, m[pos], _ring)
                    p_Setm(mon, _ring)
                    _p = p_Add_q(_p, mon, _ring)
                return co.new_MP(self, _p)

        if PY_TYPE_CHECK(element, BooleanPolynomial) and \
               element.parent().ngens() == _ring.N and \
               element.parent().variable_names() == self.variable_names():
            if element.constant():
                if element:
                    return self._one_element
                else:
                    return self._zero_element
            return eval(str(element),self.gens_dict())

        if PY_TYPE_CHECK(element, dict):
            _p = p_ISet(0, _ring)
            K = self.base_ring()
            for (m,c) in element.iteritems():
                try:
                    c = K(c)
                    if not c: continue
                except TypeError, msg:
                    p_Delete(&_p, _ring)
                    raise TypeError, msg
                mon = p_Init(_ring)
                p_SetCoeff(mon, co.sa2si(c , _ring), _ring)
                if len(m) != self.ngens():
                    raise TypeError, "tuple key must have same length as ngens"
                for pos from 0 <= pos < len(m):
                    if m[pos]:
                        p_SetExp(mon, pos+1, m[pos], _ring)
                p_Setm(mon, _ring)
                _p = p_Add_q(_p, mon, _ring)

            return co.new_MP(self, _p)

        if hasattr(element,'_polynomial_'):
            # SymbolicVariable
            return element._polynomial_(self)

        if is_Macaulay2Element(element):
            return self(repr(element))

        try:
            return self(str(element))
        except TypeError:
            pass

        # now try calling the base ring's __call__ methods
        element = self.base_ring()(element)
        _p = p_NSet(co.sa2si(element,_ring), _ring)
        return co.new_MP(self,_p)

    def _repr_(self):
        """
        EXAMPLE:
            sage: P.<x,y> = QQ[]
            sage: P
            Multivariate Polynomial Ring in x, y over Rational Field

        """
        varstr = ", ".join([ rRingVar(i,self._ring)  for i in range(self.__ngens) ])
        return "Multivariate Polynomial Ring in %s over %s"%(varstr,self._base)

    def ngens(self):
        """
        Returns the number of variables in self.

        EXAMPLES:
            sage: P.<x,y> = QQ[]
            sage: P.ngens()
            2

            sage: k.<a> = GF(2^16)
            sage: P = PolynomialRing(k,1000,'x')
            sage: P.ngens()
            1000

        """
        return int(self.__ngens)

    #def gens(self):
        #"""
        #Return the tuple of variables in self.

        #EXAMPLES:
            #sage: P.<x,y,z> = QQ[]
            #sage: P.gens()
            #(x, y, z)

            #sage: P = MPolynomialRing(QQ,10,'x')
            #sage: P.gens()
            #(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)

            #sage: P.<SAGE,SINGULAR> = MPolynomialRing(QQ,2) # weird names
            #sage: P.gens()
            #(SAGE, SINGULAR)

        #"""
        #return tuple([self.gen(i) for i in range(self.__ngens)  ])

    def gen(self, int n=0):
        """
        Returns the n-th generator of self.

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: P.gen(),P.gen(1)
            (x, y)

            sage: P = MPolynomialRing(GF(127),1000,'x')
            sage: P.gen(500)
            x500

            sage: P.<SAGE,SINGULAR> = MPolynomialRing(QQ,2) # weird names
            sage: P.gen(1)
            SINGULAR

        """
        cdef poly *_p
        cdef ring *_ring = self._ring

        if n < 0 or n >= self.__ngens:
            raise ValueError, "Generator not defined."

        rChangeCurrRing(_ring)
        _p = p_ISet(1,_ring)
        p_SetExp(_p, n+1, 1, _ring)
        p_Setm(_p, _ring);

        return co.new_MP(self,_p)

    def ideal(self, *gens, **kwds):
        """
        Create an ideal in this polynomial ring.

        INPUT:
            *gens -- list or tuple of generators (or several input
                  arguments)
            coerce -- bool (default: True); this must be a keyword
                  argument. Only set it to False if you are certain
                  that each generator is already in the ring.


        EXAMPLE:
            sage: P.<x,y,z> = QQ[]
            sage: sage.rings.ideal.Katsura(P)
            Ideal (x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y + 2*y*z - y) of Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: P.ideal([x + 2*y + 2*z-1, 2*x*y + 2*y*z-y, x^2 + 2*y^2 + 2*z^2-x])
            Ideal (x + 2*y + 2*z - 1, 2*x*y + 2*y*z - y, x^2 + 2*y^2 + 2*z^2 - x) of Multivariate Polynomial Ring in x, y, z over Rational Field

        """
        coerce = kwds.get('coerce', True)
        if len(gens) == 1:
            gens = gens[0]
        if is_SingularElement(gens):
            gens = list(gens)
            coerce = True
        elif is_Macaulay2Element(gens):
            gens = list(gens)
            coerce = True
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return MPolynomialIdeal(self, gens, coerce=False)

    def _macaulay2_(self, macaulay2=macaulay2_default):
        """
        Create a M2 representation of self if Macaulay2 is installed.

        INPUT:
            macaulay2 -- M2 interpreter (default: macaulay2_default)

        EXAMPLES:
            sage: R.<x,y> = ZZ[]
            sage: macaulay2(R)        # optional
            ZZ [x, y, MonomialOrder => GRevLex, MonomialSize => 16]

            sage: R.<x,y> = QQ[]
            sage: macaulay2(R)        # optional
            QQ [x, y, MonomialOrder => GRevLex, MonomialSize => 16]

            sage: R.<x,y> = GF(17)[]
            sage: print macaulay2(R)        # optional
            ZZ
            -- [x, y, MonomialOrder => GRevLex, MonomialSize => 16]
            17
        """
        try:
            R = self.__macaulay2
            if R is None or not (R.parent() is macaulay2):
                raise ValueError
            R._check_valid()
            return R
        except (AttributeError, ValueError):
            self.__macaulay2 = self._macaulay2_set_ring(macaulay2)
        return self.__macaulay2

    def _macaulay2_set_ring(self, macaulay2):
        if not self.__m2_set_ring_cache is None:
            base_str, gens, order = self.__m2_set_ring_cache
        else:
            if self.base_ring().is_prime_field():
                if self.characteristic() == 0:
                    base_str = "QQ"
                else:
                    base_str = "ZZ/" + str(self.characteristic())
            elif is_IntegerRing(self.base_ring()):
                base_str = "ZZ"
            else:
                raise TypeError, "no conversion of to a Macaulay2 ring defined"
            gens = str(self.gens())
            order = self.term_order().macaulay2_str()
            self.__m2_set_ring_cache = (base_str, gens, order)
        return macaulay2.ring(base_str, gens, order)

    def _can_convert_to_singular(self):
        """
        Returns True
        """
        return True

    def _singular_(self, singular=singular_default):
        """
        Create a SINGULAR (as in the CAS) representation of self. The
        result is cached.

        INPUT:
            singular -- SINGULAR interpreter (default: singular_default)

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: P._singular_()
            //   characteristic : 0
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_() is P._singular_()
            True

            sage: P._singular_().name() == P._singular_().name()
            True

            sage: k.<a> = GF(3^3)
            sage: P.<x,y,z> = PolynomialRing(k,3)
            sage: P._singular_()
            //   characteristic : 3
            //   1 parameter    : a
            //   minpoly        : (a^3-a+1)
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_() is P._singular_()
            True

            sage: P._singular_().name() == P._singular_().name()
            True


        TESTS:
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x> = MPolynomialRing_libsingular(QQ,1)
            sage: P._singular_()
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

        """
        try:
            R = self.__singular
            if R is None or not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if self.base_ring().is_prime_field():
                return R
            if self.base_ring().is_finite():
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                    singular.eval("minpoly=%s"%(self.__minpoly))
                    self.__minpoly = singular.eval('minpoly')[1:-1] # store in correct format
            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self, singular=singular_default):
        """
        Create a SINGULAR (as in the CAS) representation of self. The
        result is NOT cached.

        INPUT:
            singular -- SINGULAR interpreter (default: singular_default)

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: P._singular_init_()
            //   characteristic : 0
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C
            sage: P._singular_init_() is P._singular_init_()
            False

            sage: P._singular_init_().name() == P._singular_init_().name()
            False

        TESTS:
            sage: P.<x> = QQ[]
            sage: P._singular_init_()
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

        """
        if self.ngens()==1:
            _vars = str(self.gen())
            if "*" in _vars: # 1.000...000*x
                _vars = _vars.split("*")[1]
            order = 'lp'
        else:
            _vars = str(self.gens())
            order = self.term_order().singular_str()

        if is_RealField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = self.base_ring().precision()
            digits = sage.rings.arith.ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(real,%d,0)"%digits, _vars, order=order)

        elif is_ComplexField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = self.base_ring().precision()
            digits = sage.rings.arith.ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(complex,%d,0,I)"%digits, _vars,  order=order)

        elif self.base_ring().is_prime_field():
            self.__singular = singular.ring(self.characteristic(), _vars, order=order)

        elif self.base_ring().is_finite(): #must be extension field
            gen = str(self.base_ring().gen())
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order)
            self.__minpoly = (str(self.base_ring().modulus()).replace("x",gen)).replace(" ","")
            if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                singular.eval("minpoly=%s"%(self.__minpoly) )
                self.__minpoly = singular.eval('minpoly')[1:-1]
            self.__singular = r
        else:
            raise TypeError, "no conversion to a Singular ring defined"

        return self.__singular

    def __hash__(self):
        """
        Return a hash for self, that is, a hash of the string representation of self

        EXAMPLE:
            sage: P.<x,y,z> = QQ[]
            sage: hash(P)      # somewhat random output
            967902441410893180 # 64-bit
            -1767675994        # 32-bit
        """
        return hash(self.__repr__())

    def __richcmp__(left, right, int op):
        return (<Parent>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
        """
        Multivariate polynomial rings are said to be equal if:
         * their base rings match
         * their generator names match
         * their term orderings match

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: R.<x,y,z> = QQ[]
            sage: P == R
            True

            sage: R.<x,y,z> = MPolynomialRing(GF(127),3)
            sage: P == R
            False

            sage: R.<x,y> = MPolynomialRing(QQ,2)
            sage: P == R
            False

            sage: R.<x,y,z> = MPolynomialRing(QQ,3,order='invlex')
            sage: P == R
            False


        """
        if PY_TYPE_CHECK(right, MPolynomialRing_libsingular) or PY_TYPE_CHECK(right, MPolynomialRing_polydict_domain):
            return cmp( (left.base_ring(), map(str, left.gens()), left.term_order()),
                        (right.base_ring(), map(str, right.gens()), right.term_order())
                        )
        else:
            return cmp(type(left),type(right))

    def __reduce__(self):
        """
        Serializes self.

        EXAMPLES:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P = MPolynomialRing(GF(127),3,names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^8,'F'),3,names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^16,'B'),3,names='abc')
            sage: P == loads(dumps(P))
            True

        """
        return sage.rings.polynomial.multi_polynomial_libsingular.unpickle_MPolynomialRing_libsingular, ( self.base_ring(),
                                                                                               map(str, self.gens()),
                                                                                               self.term_order() )


    def __temporarily_change_names(self, names, latex_names):
        """
        This is used by the variable names context manager.
        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self)._ring
        cdef char **_names, **_orig_names
        cdef char *_name
        cdef int i

        if len(names) != _ring.N:
            raise TypeError, "len(names) doesn't equal self.ngens()"

        old = self._names, self._latex_names
        (self._names, self._latex_names) = names, latex_names

        _names = <char**>omAlloc0(sizeof(char*)*_ring.N)
        for i from 0 <= i < _ring.N:
            _name = names[i]
            _names[i] = omStrDup(_name)

        _orig_names = _ring.names
        _ring.names = _names

        for i from 0 <= i < _ring.N:
            omFree(_orig_names[i])
        omFree(_orig_names)

        return old

    ### The following methods are handy for implementing Groebner
    ### basis algorithms. They do only superficial type/sanity checks
    ### and should be called carefully.

    def monomial_quotient(self, MPolynomial_libsingular f, MPolynomial_libsingular g, coeff=False):
        """
        Return f/g, where both f and g are treated as
        monomials. Coefficients are ignored by default.

        INPUT:
            f -- monomial
            g -- monomial
            coeff -- divide coefficents as well (default: False)

        EXAMPLE:
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_quotient(3/2*x*y,x)
            y

            sage: P.monomial_quotient(3/2*x*y,x,coeff=True)
            3/2*y

        TESTS:
            sage: R.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_quotient(x*y,x)
            y

            sage: P.monomial_quotient(x*y,R.gen())
            y

            sage: P.monomial_quotient(P(0),P(1))
            0

            sage: P.monomial_quotient(P(1),P(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

            sage: P.monomial_quotient(P(3/2),P(2/3), coeff=True)
            9/4

            sage: P.monomial_quotient(x,y) # Note the wrong result
            x*y^1048575*z^1048575 # 64-bit
            x*y^65535*z^65535 # 32-bit

            sage: P.monomial_quotient(x,P(1))
            x

        NOTE: Assumes that the head term of f is a multiple of the
        head term of g and return the multiplicant m. If this rule is
        violated, funny things may happen.

        """
        cdef poly *res
        cdef ring *r = self._ring

        if not <ParentWithBase>self is f._parent:
            f = self._coerce_c(f)
        if not <ParentWithBase>self is g._parent:
            g = self._coerce_c(g)

        if(r != currRing): rChangeCurrRing(r)

        if not f._poly:
            return self._zero_element
        if not g._poly:
            raise ZeroDivisionError

        res = pDivide(f._poly,g._poly)
        if coeff:
            p_SetCoeff(res, n_Div( p_GetCoeff(f._poly, r) , p_GetCoeff(g._poly, r), r), r)
        else:
            p_SetCoeff(res, n_Init(1, r), r)
        return co.new_MP(self, res)

    def monomial_divides(self, MPolynomial_libsingular a, MPolynomial_libsingular b):
        """
        Return False if a does not divide b and True otherwise.

        INPUT:
            a -- monomial
            b -- monomial

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_divides(x*y*z, x^3*y^2*z^4)
            True
            sage: P.monomial_divides(x^3*y^2*z^4, x*y*z)
            False

        TESTS:
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_divides(P(1), P(0))
            True
            sage: P.monomial_divides(P(1), x)
            True

        """
        cdef poly *_a
        cdef poly *_b
        cdef ring *_r
        if a._parent is not b._parent:
            b = (<MPolynomialRing_libsingular>a._parent)._coerce_c(b)

        _a = a._poly
        _b = b._poly
        _r = (<MPolynomialRing_libsingular>a._parent)._ring

        if _a == NULL:
            raise ZeroDivisionError
        if _b == NULL:
            return True

        if not p_DivisibleBy(_a, _b, _r):
            return False
        else:
            return True


    def monomial_lcm(self, MPolynomial_libsingular f, MPolynomial_libsingular g):
        """
        LCM for monomials. Coefficients are ignored.

        INPUT:
            f -- monomial
            g -- monomial

        EXAMPLE:
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_lcm(3/2*x*y,x)
            x*y

        TESTS:
            sage: R.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_lcm(x*y,R.gen())
            x*y

            sage: P.monomial_lcm(P(3/2),P(2/3))
            1

            sage: P.monomial_lcm(x,P(1))
            x

        """
        cdef poly *m = p_ISet(1,self._ring)

        if not <ParentWithBase>self is f._parent:
            f = self._coerce_c(f)
        if not <ParentWithBase>self is g._parent:
            g = self._coerce_c(g)

        if f._poly == NULL:
            if g._poly == NULL:
                return self._zero_element
            else:
                raise ArithmeticError, "cannot compute lcm of zero and nonzero element"
        if g._poly == NULL:
            raise ArithmeticError, "cannot compute lcm of zero and nonzero element"

        if(self._ring != currRing): rChangeCurrRing(self._ring)

        pLcm(f._poly, g._poly, m)
        p_Setm(m, self._ring)
        return co.new_MP(self,m)

    def monomial_reduce(self, MPolynomial_libsingular f, G):
        """
        Try to find a g in G where g.lm() divides f. If found (flt,g)
        is returned, (0,0) otherwise, where flt is f/g.lm().

        It is assumed that G is iterable and contains ONLY elements in
        self.

        INPUT:
            f -- monomial
            G -- list/set of mpolynomials

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, 1/2  ]
            sage: P.monomial_reduce(f,G)
            (y, 1/4*x*y + 2/7)

        TESTS:
            sage: P.<x,y,z> = QQ[]
            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, 1/2  ]

            sage: P.monomial_reduce(P(0),G)
            (0, 0)

            sage: P.monomial_reduce(f,[P(0)])
            (0, 0)

        """
        cdef poly *m = f._poly
        cdef ring *r = self._ring
        cdef poly *flt

        if not m:
            return f,f

        for g in G:
            if PY_TYPE_CHECK(g, MPolynomial_libsingular) \
                   and (<MPolynomial_libsingular>g) \
                   and p_LmDivisibleBy((<MPolynomial_libsingular>g)._poly, m, r):
                flt = pDivide(f._poly, (<MPolynomial_libsingular>g)._poly)
                #p_SetCoeff(flt, n_Div( p_GetCoeff(f._poly, r) , p_GetCoeff((<MPolynomial_libsingular>g)._poly, r), r), r)
                p_SetCoeff(flt, n_Init(1, r), r)
                return co.new_MP(self,flt), g
        return self._zero_element,self._zero_element

    def monomial_pairwise_prime(self, MPolynomial_libsingular g, MPolynomial_libsingular h):
        """
        Return True if h and g are pairwise prime. Both are treated as monomials.

        INPUT:
            h -- monomial
            g -- monomial

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_pairwise_prime(x^2*z^3, y^4)
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, 3/4*y^3)
            False

        TESTS:
            sage: Q.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_pairwise_prime(x^2*z^3, Q('y^4'))
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, Q(0))
            True

            sage: P.monomial_pairwise_prime(P(1/2),x)
            False


        """
        cdef int i
        cdef ring *r
        cdef poly *p, *q

        if h._parent is not g._parent:
            g = (<MPolynomialRing_libsingular>h._parent)._coerce_c(g)

        r = (<MPolynomialRing_libsingular>h._parent)._ring
        p = g._poly
        q = h._poly

        if p == NULL:
            if q == NULL:
                return False #GCD(0,0) = 0
            else:
                return True #GCD(x,0) = 1

        elif q == NULL:
            return True # GCD(0,x) = 1

        elif p_IsConstant(p,r) or p_IsConstant(q,r): # assuming a base field
            return False

        for i from 1 <= i <= r.N:
            if p_GetExp(p,i,r) and p_GetExp(q,i,r):
                return False
        return True



    def monomial_all_divisors(self, MPolynomial_libsingular t):
        """
        Return a list of all monomials that divide t, coefficients are
        ignored.

        INPUT:
            t -- a monomial

        OUTPUT:
            a list of monomials


        EXAMPLE:
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_all_divisors(x^2*z^3)
            [x, x^2, z, x*z, x^2*z, z^2, x*z^2, x^2*z^2, z^3, x*z^3, x^2*z^3]

        ALGORITHM: addwithcarry idea by Toon Segers
        """

        M = list()

        cdef ring *_ring = self._ring
        cdef poly *maxvector = t._poly
        cdef poly *tempvector = p_ISet(1, _ring)

        pos = 1

        while not p_ExpVectorEqual(tempvector, maxvector, _ring):
          tempvector = addwithcarry(tempvector, maxvector, pos, _ring)
          M.append(co.new_MP(self, p_Copy(tempvector,_ring)))
        return M


cdef inline int polyLengthBounded(poly *p, int bound):
    cdef poly *n = p
    cdef int count = 0
    while n and count < bound:
        n = pNext(n)
        count += 1
    return count

def unpickle_MPolynomialRing_libsingular(base_ring, names, term_order):
    """
    inverse function for MPolynomialRing_libsingular.__reduce__

    """
    from sage.rings.polynomial.polynomial_ring_constructor import _get_from_cache
    key = (base_ring,  tuple(names), len(names), False, term_order)
    R = _get_from_cache(key)
    if not R is None:
        return R

    return MPolynomialRing_libsingular(base_ring, len(names), names, term_order)

cdef inline MPolynomial_libsingular new_MP(MPolynomialRing_libsingular parent, poly *juice):
    """
    Construct a new MPolynomial_libsingular element
    """
    cdef MPolynomial_libsingular p
    p = PY_NEW(MPolynomial_libsingular)
    p._parent = <ParentWithBase>parent
    p._poly = juice
    return p

cdef class MPolynomial_libsingular(sage.rings.polynomial.multi_polynomial.MPolynomial):
    """
    A multivariate polynomial implemented using libSINGULAR.
    """
    def __init__(self, MPolynomialRing_libsingular parent):
        """
        Construct a zero element in parent.
        """
        self._poly = NULL
        self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        # for some mysterious reason, various things may be NULL in some cases
        if self._parent is not <ParentWithBase>None and (<MPolynomialRing_libsingular>self._parent)._ring != NULL and self._poly != NULL:
            p_Delete(&self._poly, (<MPolynomialRing_libsingular>self._parent)._ring)

    def __call__(self, *x, **kwds):
        r"""
        Evaluate this multi-variate polynomial at $x$, where $x$ is
        either the tuple of values to substitute in, or one can use
        functional notation $f(a_0,a_1,a_2, \ldots)$ to evaluate $f$
        with the ith variable replaced by $a_i$.

        INPUT:
            x -- a list of elements in self.parent()
            or **kwds -- a dictionary of variable-name:value pairs.

        EXAMPLES:
            sage: P.<x,y,z> = QQ[]
            sage: f = 3/2*x^2*y + 1/7 * y^2 + 13/27
            sage: f(0,0,0)
            13/27

            sage: f(1,1,1)
            803/378
            sage: 3/2 + 1/7 + 13/27
            803/378

            sage: f(45/2,19/3,1)
            7281167/1512

            sage: f(1,2,3).parent()
            Rational Field

        TESTS:
            sage: P.<x,y,z> = QQ[]
            sage: P(0)(1,2,3)
            0
            sage: P(3/2)(1,2,3)
            3/2

            sage: R.<a,b,y> = QQ[]
            sage: f = a*y^3 + b*y - 3*a*b*y
            sage: f(a=5,b=3,y=10)
            4580
            sage: f(5,3,10)
            4580
        """
        if len(kwds) > 0:
            f = self.subs(**kwds)
            if len(x) > 0:
                return f(*x)
            else:
                return f

        cdef int l = len(x)
        cdef MPolynomialRing_libsingular parent = (<MPolynomialRing_libsingular>self._parent)
        cdef ring *_ring = parent._ring

        cdef poly *_p

        if l == 1 and (PY_TYPE_CHECK(x, tuple) or PY_TYPE_CHECK(x, list)):
            x = x[0]
            l = len(x)

        if l != parent._ring.N:
            raise TypeError, "number of arguments does not match number of variables in parent"

        try:
            x = [parent._coerce_c(e) for e in x]
        except TypeError:
            # give up, evaluate functional
            y = parent.base_ring()(0)
            for (m,c) in self.dict().iteritems():
                y += c*mul([ x[i]**m[i] for i in m.nonzero_positions()])
            return y

        cdef ideal *to_id = idInit(l,1)
        for i from 0 <= i < l:
            to_id.m[i]= p_Copy( (<MPolynomial_libsingular>x[i])._poly, _ring)

        cdef ideal *from_id=idInit(1,1)
        from_id.m[0] = self._poly

        cdef ideal *res_id = fast_map(from_id, _ring, to_id, _ring)
        cdef poly *res = res_id.m[0]

        from_id.m[0] = NULL
        res_id.m[0] = NULL

        id_Delete(&to_id, _ring)
        id_Delete(&from_id, _ring)
        id_Delete(&res_id, _ring)

        if p_IsConstant(res, _ring) and all([e in parent._base for e in x]):
            # I am sure there must be a better way to do this...
            return parent._base(co.new_MP(parent, res))
        else:
            return co.new_MP(parent, res)

    # you may have to replicate this boilerplate code in derived classes if you override
    # __richcmp__.  The python documentation at  http://docs.python.org/api/type-structs.html
    # explains how __richcmp__, __hash__, and __cmp__ are tied together.
    def __hash__(self):
        return self._hash_c()

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Compare left and right and return -1, 0, and 1 for <,==, and > respectively.

        EXAMPLES:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3, order='degrevlex')
            sage: x == x
            True

            sage: x > y
            True
            sage: y^2 > x
            True

            sage: (2/3*x^2 + 1/2*y + 3) > (2/3*x^2 + 1/4*y + 10)
            True

        TESTS:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3, order='degrevlex')
            sage: x > P(0)
            True

            sage: P(0) == P(0)
            True

            sage: P(0) < P(1)
            True

            sage: x > P(1)
            True

            sage: 1/2*x < 3/4*x
            True

            sage: (x+1) > x
            True

            sage: f = 3/4*x^2*y + 1/2*x + 2/7
            sage: f > f
            False
            sage: f < f
            False
            sage: f == f
            True

            sage: P.<x,y,z> = MPolynomialRing(GF(127),3, order='degrevlex')
            sage: (66*x^2 + 23) > (66*x^2 + 2)
            True


        """
        cdef ring *r
        cdef poly *p, *q
        cdef number *h
        cdef int ret = 0

        if left is right:
            return ret

        r = (<MPolynomialRing_libsingular>left._parent)._ring
        if(r != currRing): rChangeCurrRing(r)
        p = (<MPolynomial_libsingular>left)._poly
        q = (<MPolynomial_libsingular>right)._poly

        # handle special cases first (slight slowdown, as special
        # cases are - well - special
        if p==NULL:
            if q==NULL:
                # compare 0, 0
                return 0
            elif p_IsConstant(q,r):
                # compare 0, const
                return 1-2*n_GreaterZero(p_GetCoeff(q,r), r) # -1: <, 1: > #
        elif q==NULL:
            if p_IsConstant(p,r):
                # compare const, 0
                return -1+2*n_GreaterZero(p_GetCoeff(p,r), r) # -1: <, 1: >
        #else

        while ret==0 and p!=NULL and q!=NULL:
            ret = p_Cmp( p, q, r)

            if ret==0:
                h = n_Sub(p_GetCoeff(p, r),p_GetCoeff(q, r), r)
                # compare coeffs
                ret = -1+n_IsZero(h, r)+2*n_GreaterZero(h, r) # -1: <, 0:==, 1: >
                n_Delete(&h, r)
            p = pNext(p)
            q = pNext(q)

        if ret==0:
            if p==NULL and q != NULL:
                ret = -1
            elif p!=NULL and q==NULL:
                ret = 1

        return ret

    cdef ModuleElement _add_c_impl( left, ModuleElement right):
        """
        Add left and right.

        EXAMPLE:
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: 3/2*x + 1/2*y + 1
            3/2*x + 1/2*y + 1

        """
        cdef MPolynomial_libsingular res

        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        _l = p_Copy(left._poly, _ring)
        _r = p_Copy((<MPolynomial_libsingular>right)._poly, _ring)

        _p= p_Add_q(_l, _r, _ring)

        return co.new_MP((<MPolynomialRing_libsingular>left._parent),_p)

    cdef ModuleElement _iadd_c_impl( left, ModuleElement right):
        """
        Add left and right inplace.

        EXAMPLE:
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: 3/2*x + 1/2*y + 1
            3/2*x + 1/2*y + 1

        """
        cdef MPolynomial_libsingular res

        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        _l = left._poly
        _r = p_Copy((<MPolynomial_libsingular>right)._poly, _ring)

        _p= p_Add_q(_l, _r, _ring)

        left._poly = _p
        return left

    cdef ModuleElement _sub_c_impl( left, ModuleElement right):
        """
        Subtract left and right.

        EXAMPLE:
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: 3/2*x - 1/2*y - 1
            3/2*x - 1/2*y - 1

        """
        cdef MPolynomial_libsingular res

        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        _l = p_Copy(left._poly, _ring)
        _r = p_Copy((<MPolynomial_libsingular>right)._poly, _ring)

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p= p_Add_q(_l, p_Neg(_r, _ring), _ring)

        return co.new_MP((<MPolynomialRing_libsingular>left._parent),_p)

    cdef ModuleElement _isub_c_impl( left, ModuleElement right):
        """
        Subtract left and right inplace.

        EXAMPLE:
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: 3/2*x - 1/2*y - 1
            3/2*x - 1/2*y - 1

        """
        cdef MPolynomial_libsingular res

        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        _l = left._poly
        _r = p_Copy((<MPolynomial_libsingular>right)._poly, _ring)

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p= p_Add_q(_l, p_Neg(_r, _ring), _ring)
        left._poly = _p
        return left

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        """
        Multiply self with a base ring element.

        EXAMPLE:
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: 3/2*x
            3/2*x
        """

        cdef number *_n
        cdef ring *_ring
        cdef poly *_p

        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        if not left:
            return (<MPolynomialRing_libsingular>self._parent)._zero_element

        _n = co.sa2si(left,_ring)

        _p = pp_Mult_nn(self._poly,_n,_ring)
        n_Delete(&_n, _ring)
        return co.new_MP((<MPolynomialRing_libsingular>self._parent),_p)

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        # all currently implemented rings are commutative
        return self._rmul_c_impl(right)

    cdef RingElement  _mul_c_impl(left, RingElement right):
        # all currently implemented rings are commutative
        """
        Multiply left and right.

        EXAMPLE:
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2*x + 1/2*y + 1)
            9/4*x^2 - 1/4*y^2 - y - 1
        """
        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p = pp_Mult_qq(left._poly, (<MPolynomial_libsingular>right)._poly, _ring)
        return co.new_MP(left._parent,_p)

    cdef RingElement  _imul_c_impl(left, RingElement right):
        # all currently implemented rings are commutative
        """
        Multiply left and right inplace.

        EXAMPLE:
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2*x + 1/2*y + 1)
            9/4*x^2 - 1/4*y^2 - y - 1
        """
        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p = pp_Mult_qq(left._poly, (<MPolynomial_libsingular>right)._poly, _ring)
        p_Delete(&left._poly, _ring)
        left._poly = _p
        return left

    cdef RingElement  _div_c_impl(left, RingElement right):
        """
        Divide left by right

        EXAMPLES:
            sage: R.<x,y>=MPolynomialRing(QQ,2)
            sage: f = (x + y)/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        Note that / is still a constructor for elements of the
        fraction field in all cases as long as both arguments have the
        same parent.

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: R.<x,y>=MPolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: g = x
            sage: h = f/g; h
            (x^3 + y)/x
            sage: h.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field

        TESTS:
            sage: R.<x,y>=MPolynomialRing(QQ,2)
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        """
        cdef poly *p
        cdef ring *r
        cdef number *n
        if (<MPolynomial_libsingular>right).is_constant_c():

            p = (<MPolynomial_libsingular>right)._poly
            if p == NULL:
                raise ZeroDivisionError
            r = (<MPolynomialRing_libsingular>(<MPolynomial_libsingular>left)._parent)._ring
            n = p_GetCoeff(p, r)
            n = nInvers(n)
            p = pp_Mult_nn(left._poly, n,  r)
            n_Delete(&n,r)
            return co.new_MP(left._parent, p)
        else:
            return (<MPolynomialRing_libsingular>left._parent).fraction_field()(left,right)

    def __pow__(MPolynomial_libsingular self,int exp,ignored):
        """
        Return self^(exp).

        EXAMPLE:
            sage: R.<x,y>=MPolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f^2
            x^6 + 2*x^3*y + y^2
            sage: g = f^(-1); g
            1/(x^3 + y)
            sage: type(g)
            <class 'sage.rings.fraction_field_element.FractionFieldElement'>
        """
        cdef ring *_ring
        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        cdef poly *_p
        cdef int _exp

        _exp = exp

        if _exp < 0:
            return 1/(self**(-_exp))
        if _exp > 65535:
            raise TypeError,  "exponent is too large, max. is 65535"

        if(_ring != currRing): rChangeCurrRing(_ring)
        cdef int count = polyLengthBounded(self._poly,15)
        if count >= 15 or _exp > 15:
            _sig_on
        _p = pPower( p_Copy(self._poly,_ring),_exp)
        if count >= 15 or _exp > 15:
            _sig_off
        return co.new_MP((<MPolynomialRing_libsingular>self._parent),_p)

    def __neg__(self):
        """
        Return -self.

        EXAMPLE:
            sage: R.<x,y>=MPolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: -f
            -x^3 - y
        """
        cdef ring *_ring
        _ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if(_ring != currRing): rChangeCurrRing(_ring)

        return co.new_MP((<MPolynomialRing_libsingular>self._parent),\
                                           p_Neg(p_Copy(self._poly,_ring),_ring))

    def _repr_(self):
        s =  self._repr_short_c()
        s = s.replace("+"," + ").replace("-"," - ")
        if s.startswith(" - "):
            return "-" + s[3:]
        else:
            return s

    def _repr_short(self):
        """
        This is a faster but less pretty way to print polynomials. If available
        it uses the short SINGULAR notation.
        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if _ring.CanShortOut:
            _ring.ShortOut = 1
            s = self._repr_short_c()
            _ring.ShortOut = 0
        else:
            s = self._repr_short_c()
        return s

    cdef _repr_short_c(self):
        """
        Raw SINGULAR printing.
        """
        rChangeCurrRing((<MPolynomialRing_libsingular>self._parent)._ring)
        s = p_String(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring, (<MPolynomialRing_libsingular>self._parent)._ring)
        return s

    def _latex_(self):
        r"""
        Return a polynomial latex representation of self.

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: latex(f)
            - x^{2} y - \frac{25}{27} y^{3} - z^{2}

        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef int n = _ring.N
        cdef int j, e
        cdef poly *p = self._poly
        poly = ""
        gens = self.parent().latex_variable_names()
        base = self.parent().base()

        while p:
            sign_switch = False

            # First determine the multinomial:
            multi = ""
            for j from 1 <= j <= n:
                e = p_GetExp(p, j, _ring)
                if e > 0:
                    multi += " "+gens[j-1]
                if e > 1:
                    multi += "^{%d}"%e
            multi = multi.lstrip().rstrip()

            # Next determine coefficient of multinomial
            c =  co.si2sa( p_GetCoeff(p, _ring), _ring, base)
            if len(multi) == 0:
                multi = latex(c)
            elif c != 1:
                if  c == -1:
                    if len(poly) > 0:
                        sign_switch = True
                    else:
                        multi = "- %s"%(multi)
                else:
                    multi = "%s %s"%(latex(c),multi)

            # Now add on coefficiented multinomials
            if len(poly) > 0:
                if sign_switch:
                    poly = poly + " - "
                else:
                    poly = poly + " + "
            poly = poly + multi

            p = pNext(p)

        poly = poly.lstrip().rstrip()
        poly = poly.replace("+ -","- ")

        if len(poly) == 0:
            return "0"
        return poly

    def _macaulay2_(self, macaulay2=macaulay2):
        """
        Return corresponding Macaulay2 polynomial.

        WARNING: Two identical rings are not canonically isomorphic in
        M2, so we require the user to explicitly set the ring, since
        there is no way to know if the ring has been set or not, and
        setting it twice screws everything up.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(GF(7), 2)   # optional
            sage: f = (x^3 + 2*y^2*x)^7; f          # optional
            x^21 + 2*x^7*y^14

        Always call the Macaulay2 ring conversion on the parent polynomial
        ring before converting a copy of elements to Macaulay2:
            sage: macaulay2(R)                      # optional
            ZZ/7 [x, y, MonomialOrder => GRevLex, MonomialSize => 16]
            sage: h = f._macaulay2_(); h            # optional
            x^21+2*x^7*y^14
            sage: k = (x+y)._macaulay2_()           # optional
            sage: k + h                             # optional
            x^21+2*x^7*y^14+x+y
            sage: R(h)                              # optional
            x^21 + 2*x^7*y^14
            sage: R(h^20) == f^20                   # optional
            True
        """
        try:
            if self.__macaulay2[macaulay2].parent() is macaulay2:
                return self.__macaulay2[macaulay2]
        except (TypeError, AttributeError):
            self.__macaulay2 = {}
        except KeyError:
            pass
        #self.parent()._macaulay2_set_ring(macaulay2)
        z = macaulay2(repr(self))
        self.__macaulay2[macaulay2] = z
        return z

    def _repr_with_changed_varnames(self, varnames):
        """
        Return string representing self but change the variable names
        to varnames.

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: print f._repr_with_changed_varnames(['FOO', 'BAR', 'FOOBAR'])
            -FOO^2*BAR - 25/27*BAR^3 - FOOBAR^2

        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef char **_names
        cdef char **_orig_names
        cdef char *_name
        cdef int i

        if len(varnames) != _ring.N:
            raise TypeError, "len(varnames) doesn't equal self.parent().ngens()"


        _names = <char**>omAlloc0(sizeof(char*)*_ring.N)
        for i from 0 <= i < _ring.N:
            _name = varnames[i]
            _names[i] = omStrDup(_name)

        _orig_names = _ring.names
        _ring.names = _names
        s = str(self)
        _ring.names = _orig_names

        for i from 0 <= i < _ring.N:
            omFree(_names[i])
        omFree(_names)

        return s

    def degree(self, MPolynomial_libsingular x=None):
        """
        Return the maximal degree of self in x, where x must be one of the
        generators for the parent of self.

        INPUT:
            x -- multivariate polynmial (a generator of the parent of self)
                 If x is not specified (or is None), return the total degree,
                 which is the maximum degree of any monomial.

        OUTPUT:
            integer

        EXAMPLE:
            sage: R.<x, y> = MPolynomialRing(QQ, 2)
            sage: f = y^2 - x^9 - x
            sage: f.degree(x)
            9
            sage: f.degree(y)
            2
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(x)
            3
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(y)
            10

        TESTS:
            sage: P.<x, y> = MPolynomialRing(QQ, 2)
            sage: P(0).degree(x)
            0
            sage: P(1).degree(x)
            0

        """
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef poly *p = self._poly
        cdef int deg, _deg

        deg = 0

        if not x:
            return self.total_degree()

        # TODO: we can do this faster
        if not x in self._parent.gens():
            raise TypeError, "x must be one of the generators of the parent."
        for i from 1 <= i <= r.N:
            if p_GetExp(x._poly, i, r):
                break
        while p:
            _deg =  p_GetExp(p,i,r)
            if _deg > deg:
                deg = _deg
            p = pNext(p)

        return deg

    def newton_polytope(self):
        """
        Return the Newton polytope of this polynomial.

        You should have the optional polymake package installed.

        EXAMPLES:
            sage: R.<x,y> = MPolynomialRing(QQ,2)
            sage: f = 1 + x*y + x^3 + y^3
            sage: P = f.newton_polytope()
            sage: P
            Convex hull of points [[1, 0, 0], [1, 0, 3], [1, 1, 1], [1, 3, 0]]
            sage: P.facets()
            [(0, 1, 0), (3, -1, -1), (0, 0, 1)]
            sage: P.is_simple()
            True

        TESTS:
            sage: R.<x,y> = MPolynomialRing(QQ,2)
            sage: R(0).newton_polytope()
            Convex hull of points []
            sage: R(1).newton_polytope()
            Convex hull of points [[1, 0, 0]]

        """
        from sage.geometry.all import polymake
        e = self.exponents()
        a = [[1] + list(v) for v in e]
        P = polymake.convex_hull(a)
        return P

    def total_degree(self):
        """
        Return the total degree of self, which is the maximum degree
        of all monomials in self.

        EXAMPLES:
            sage: R.<x,y,z> = MPolynomialRing(QQ, 3)
            sage: f=2*x*y^3*z^2
            sage: f.total_degree()
            6
            sage: f=4*x^2*y^2*z^3
            sage: f.total_degree()
            7
            sage: f=99*x^6*y^3*z^9
            sage: f.total_degree()
            18
            sage: f=x*y^3*z^6+3*x^2
            sage: f.total_degree()
            10
            sage: f=z^3+8*x^4*y^5*z
            sage: f.total_degree()
            10
            sage: f=z^9+10*x^4+y^8*x^2
            sage: f.total_degree()
            10

        TESTS:
            sage: R.<x,y,z> = MPolynomialRing(QQ, 3)
            sage: R(0).total_degree()
            0
            sage: R(1).total_degree()
            0
        """
        cdef poly *p = self._poly
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef int l
        if self._poly == NULL:
            return 0
        if(r != currRing): rChangeCurrRing(r)
        return pLDeg(p,&l,r)

    def coefficient(self, degrees):
        """
        Return the coefficient of the variables with the degrees
        specified in the python dictionary \code{degrees}.  Mathematically,
        this is the coefficient in the base ring adjoined by the variables
        of this ring not listed in \code{degrees}.  However, the result
        has the same parent as this polynomial.

        This function contrasts with the function \code{monomial_coefficient}
        which returns the coefficient in the base ring of a monomial.

        INPUT:
            degrees -- Can be any of:
                -- a dictionary of degree restrictions
                -- a list of degree restrictions (with None in the unrestricted variables)
                -- a monomial (very fast, but not as flexible)

        OUTPUT:
            element of the parent of self

        SEE ALSO:
            For coefficients of specific monomials, look at \ref{monomial_coefficient}.

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: f=x*y+y+5
            sage: f.coefficient({x:0,y:1})
            1
            sage: f.coefficient({x:0})
            y + 5
            sage: f=(1+y+y^2)*(1+x+x^2)
            sage: f.coefficient({x:0})
            y^2 + y + 1
            sage: f.coefficient([0,None])
            y^2 + y + 1
            sage: f.coefficient(x)
            y^2 + y + 1
            sage: # Be aware that this may not be what you think!
            sage: # The physical appearance of the variable x is deceiving -- particularly if the exponent would be a variable.
            sage: f.coefficient(x^0) # outputs the full polynomial
            x^2*y^2 + x^2*y + x*y^2 + x^2 + x*y + y^2 + x + y + 1
            sage: R.<x,y> = GF(389)[]
            sage: f=x*y+5
            sage: c=f.coefficient({x:0,y:0}); c
            5
            sage: parent(c)
            Multivariate Polynomial Ring in x, y over Finite Field of size 389

        AUTHOR:
            -- Joel B. Mohler (2007.10.31)
        """
        cdef poly *_degrees = <poly*>0
        cdef poly *p = self._poly
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef poly *newp = p_ISet(0,r)
        cdef poly *newptemp
        cdef int i
        cdef int flag
        cdef int gens = self._parent.ngens()
        cdef int *exps = <int*>sage_malloc(sizeof(int)*gens)
        for i from 0<=i<gens:
            exps[i] = -1

        if PY_TYPE_CHECK(degrees, MPolynomial_libsingular) and self._parent is (<MPolynomial_libsingular>degrees)._parent:
            _degrees = (<MPolynomial_libsingular>degrees)._poly
            if pLength(_degrees) != 1:
                raise TypeError, "degrees must be a monomial"
            for i from 0<=i<gens:
                if p_GetExp(_degrees,i+1,r)!=0:
                    exps[i] = p_GetExp(_degrees,i+1,r)
        elif type(degrees) is list:
            for i from 0<=i<gens:
                if degrees[i] is None:
                    exps[i] = -1
                else:
                    exps[i] = int(degrees[i])
        elif type(degrees) is dict:
            # Extract the ordered list of degree specifications from the dictionary
            poly_vars = self.parent().gens()
            for i from 0<=i<gens:
                try:
                    exps[i] = degrees[poly_vars[i]]
                except KeyError:
                    pass
        else:
            raise TypeError, "The input degrees must be a dictionary of variables to exponents."

        # Extract the monomials that match the specifications
        while(p):
            flag = 0
            for i from 0<=i<gens:
                if exps[i] != -1 and p_GetExp(p,i+1,r)!=exps[i]:
                    #print i, p_GetExp(p,i+1,r), exps[i]
                    flag = 1
            if flag == 0:
                newptemp = p_LmInit(p,r)
                p_SetCoeff(newptemp,n_Copy(p_GetCoeff(p,r),r),r)
                for i from 0<=i<gens:
                    if exps[i] != -1:
                        p_SetExp(newptemp,i+1,0,r)
                p_Setm(newptemp,r)
                newp = p_Add_q(newp,newptemp,r)
            p = pNext(p)

        sage_free(exps)

        return co.new_MP(self.parent(),newp)

    def monomial_coefficient(self, MPolynomial_libsingular mon):
        """
        Return the coefficient in the base ring of the monomial mon in self, where mon
        must have the same parent as self.

        This function contrasts with the function \code{coefficient}
        which returns the coefficient of a monomial viewing this polynomial in a
        polynomial ring over a base ring having fewer variables.

        INPUT:
            mon -- a monomial

        OUTPUT:
            coefficient in base ring

        SEE ALSO:
            For coefficients in a base ring of fewer variables, look at \ref{coefficient}.

        EXAMPLES:
            sage: P.<x,y> = MPolynomialRing(QQ)

            The parent of the return is a member of the base ring.
            sage: f = 2 * x * y
            sage: c = f.monomial_coefficient(x*y); c
            2
            sage: c.parent()
            Rational Field

            sage: f = y^2 + y^2*x - x^9 - 7*x + 5*x*y
            sage: f.monomial_coefficient(y^2)
            1
            sage: f.monomial_coefficient(x*y)
            5
            sage: f.monomial_coefficient(x^9)
            -1
            sage: f.monomial_coefficient(x^10)
            0
        """
        cdef poly *p = self._poly
        cdef poly *m = mon._poly
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring

        if not mon._parent is self._parent:
            raise TypeError, "mon must have same parent as self"

        while(p):
            if p_ExpVectorEqual(p, m, r) == 1:
                return co.si2sa(p_GetCoeff(p, r), r, (<MPolynomialRing_libsingular>self._parent)._base)
            p = pNext(p)

        return (<MPolynomialRing_libsingular>self._parent)._base._zero_element

    def dict(self):
        """
        Return a dictionary representing self. This dictionary is in
        the same format as the generic MPolynomial: The dictionary
        consists of ETuple:coefficient pairs.

        EXAMPLE:
            sage: R.<x,y,z> = MPolynomialRing(QQ, 3)
            sage: f=2*x*y^3*z^2 + 1/7*x^2 + 2/3
            sage: f.dict()
            {(1, 3, 2): 2, (0, 0, 0): 2/3, (2, 0, 0): 1/7}
        """
        cdef poly *p
        cdef ring *r
        cdef int n
        cdef int v
        r = (<MPolynomialRing_libsingular>self._parent)._ring
        if r!=currRing: rChangeCurrRing(r)
        base = (<MPolynomialRing_libsingular>self._parent)._base
        p = self._poly
        pd = dict()
        while p:
            d = dict()
            for v from 1 <= v <= r.N:
                n = p_GetExp(p,v,r)
                if n!=0:
                    d[v-1] = n

            pd[ETuple(d,r.N)] = co.si2sa(p_GetCoeff(p, r), r, base)

            p = pNext(p)
        return pd

    cdef long _hash_c(self):
        """
        This hash incorporates the variable name in an effort to respect the obvious inclusions
        into multi-variable polynomial rings.

        The tuple algorithm is borrowed from http://effbot.org/zone/python-hash.htm.

        EXAMPLES:
            sage: R.<x>=QQ[]
            sage: S.<x,y>=QQ[]
            sage: hash(S(1/2))==hash(1/2)  # respect inclusions of the rationals
            True
            sage: hash(S.0)==hash(R.0)  # respect inclusions into mpoly rings
            True
            sage: # the point is to make for more flexible dictionary look ups
            sage: d={S.0:12}
            sage: d[R.0]
            12
        """
        cdef poly *p
        cdef ring *r
        cdef int n
        cdef int v
        r = (<MPolynomialRing_libsingular>self._parent)._ring
        if r!=currRing: rChangeCurrRing(r)
        base = (<MPolynomialRing_libsingular>self._parent)._base
        p = self._poly
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        var_name_hash = [hash(vn) for vn in self._parent.variable_names()]
        cdef long c_hash
        while p:
            c_hash = hash(co.si2sa(p_GetCoeff(p, r), r, base))
            if c_hash != 0: # this is always going to be true, because we are sparse (correct?)
                # Hash (self[i], gen_a, exp_a, gen_b, exp_b, gen_c, exp_c, ...) as a tuple according to the algorithm.
                # I omit gen,exp pairs where the exponent is zero.
                result_mon = c_hash
                for v from 1 <= v <= r.N:
                    n = p_GetExp(p,v,r)
                    if n!=0:
                        result_mon = (1000003 * result_mon) ^ var_name_hash[v-1]
                        result_mon = (1000003 * result_mon) ^ n
                result += result_mon

            p = pNext(p)
        if result == -1:
            return -2
        return result

    def __iter__(self):
        """
        Facilitates iterating over the monomials of self,
        returning tuples of the form (coeff, mon) for each
        non-zero monomial.

        NOTE: This function creates the entire list upfront because
              Cython doesn't (yet) support iterators.

        EXAMPLES:
            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: f = 3*x^3*y + 16*x + 7
            sage: [(c,m) for c,m in f]
            [(3, x^3*y), (16, x), (7, 1)]
            sage: f = P.random_element(12,14)
            sage: sum(c*m for c,m in f) == f
            True
        """
        L = zip(self.coefficients(), self.monomials())
        return iter(L)

    def __getitem__(self,x):
        """
        same as self.monomial_coefficent but for exponent vectors.

        INPUT:
            x -- a tuple or, in case of a single-variable MPolynomial
                 ring x can also be an integer.

        EXAMPLES:
            sage: R.<x, y> = MPolynomialRing(QQ, 2)
            sage: f = -10*x^3*y + 17*x*y
            sage: f[3,1]
            -10
            sage: f[1,1]
            17
            sage: f[0,1]
            0

            sage: R.<x> = MPolynomialRing(GF(7),1); R
            Multivariate Polynomial Ring in x over Finite Field of size 7
            sage: f = 5*x^2 + 3; f
            -2*x^2 + 3
            sage: f[2]
            5
        """

        cdef poly *m
        cdef poly *p = self._poly
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef int i

        if PY_TYPE_CHECK(x, MPolynomial_libsingular):
            return self.monomial_coefficient(x)
        if not PY_TYPE_CHECK(x, tuple):
            try:
                x = tuple(x)
            except TypeError:
                x = (x,)

        if len(x) != (<MPolynomialRing_libsingular>self._parent).__ngens:
            raise TypeError, "x must have length self.ngens()"

        m = p_ISet(1,r)
        i = 1
        for e in x:
            p_SetExp(m, i, int(e), r)
            i += 1
        p_Setm(m, r)

        while(p):
            if p_ExpVectorEqual(p, m, r) == 1:
                p_Delete(&m,r)
                return co.si2sa(p_GetCoeff(p, r), r, (<MPolynomialRing_libsingular>self._parent)._base)
            p = pNext(p)

        p_Delete(&m,r)
        return (<MPolynomialRing_libsingular>self._parent)._base._zero_element

    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES:
            sage: R.<a,b,c> = PolynomialRing(QQ, 3)
            sage: R.<a,b,c> = MPolynomialRing(QQ, 3)
            sage: f = a^3 + b + 2*b^2
            sage: f.exponents()
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
        """
        cdef poly *p
        cdef ring *r
        cdef int v
        r = (<MPolynomialRing_libsingular>self._parent)._ring

        p = self._poly

        pl = list()
        while p:
            ml = list()
            for v from 1 <= v <= r.N:
                ml.append(p_GetExp(p,v,r))
            pl.append(ETuple(ml))

            p = pNext(p)
        return pl

    def is_unit(self):
        """
        Return True if self is a unit.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: R.<x,y> = MPolynomialRing(QQ, 2)
            sage: (x+y).is_unit()
            False
            sage: R(0).is_unit()
            False
            sage: R(-1).is_unit()
            True
            sage: R(-1 + x).is_unit()
            False
            sage: R(2).is_unit()
            True
        """
        return bool(p_IsUnit(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring))

    def inverse_of_unit(self):
        """
        Return the inverse of self if self is a unit.

        EXAMPLES:

            sage: R.<x,y> = MPolynomialRing(QQ, 2)
        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if(_ring != currRing): rChangeCurrRing(_ring)

        if not p_IsUnit(self._poly, _ring):
            raise ArithmeticError, "is not a unit"
        else:
            return co.new_MP(self._parent,pInvers(0,self._poly,NULL))

    def is_homogeneous(self):
        """
        Return True if self is a homogeneous polynomial.

        EXAMPLES:
            sage: P.<x,y> = PolynomialRing(RationalField(), 2)
            sage: (x+y).is_homogeneous()
            True
            sage: (x.parent()(0)).is_homogeneous()
            True
            sage: (x+y^2).is_homogeneous()
            False
            sage: (x^2 + y^2).is_homogeneous()
            True
            sage: (x^2 + y^2*x).is_homogeneous()
            False
            sage: (x^2*y + y^2*x).is_homogeneous()
            True

        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if(_ring != currRing): rChangeCurrRing(_ring)
        return bool(pIsHomogeneous(self._poly))

    cpdef _homogenize(self, int var):
        r"""
        Return \code{self} if \code{self} is homogeneous.  Otherwise
        return a homogenized polynomial constructed by modifying the
        degree of the variable with index \code{var}.

        INPUT:
            var -- an integer indicating which variable to use to
                    homogenize (0 <= var < parent(self).ngens())

        OUTPUT:
            a multivariate polynomial

        EXAMPLES:
            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: P.<x,y> = MPolynomialRing(QQ,2)
            sage: f = x^2 + y + 1 + 5*x*y^10
            sage: g = f.homogenize('z'); g # indirect doctest
            5*x*y^10 + x^2*z^9 + y*z^10 + z^11
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: f._homogenize(0)
            2*x^11 + x^10*y + 5*x*y^10

        SEE: \code{self.homogenize}
        """
        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>self._parent
        cdef MPolynomial_libsingular f

        if self.is_homogeneous():
            return self

        if var < parent._ring.N:
            return co.new_MP(parent, pHomogen(p_Copy(self._poly, parent._ring), var+1))
        else:
            raise TypeError, "var must be < self.parent().ngens()"

    def is_monomial(self):
        return not self._poly.next

    def subs(self, fixed=None, **kw):
        """
        Fixes some given variables in a given multivariate polynomial and
        returns the changed multivariate polynomials. The polynomial
        itself is not affected.  The variable,value pairs for fixing are
        to be provided as dictionary of the form {variable:value}.

        This is a special case of evaluating the polynomial with some of
        the variables constants and the others the original variables, but
        should be much faster if only few variables are to be fixed.

        INPUT:
            fixed -- (optional) dict with variable:value pairs
            **kw -- names parameters

        OUTPUT:
            new MPolynomial

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f(5,y)
            25*y^2 + y + 30
            sage: f.subs({x:5})
            25*y^2 + y + 30
            sage: f.subs(x=5)
            25*y^2 + y + 30

            sage: P.<x,y,z> = PolynomialRing(GF(2),3)
            sage: f = x + y + 1
            sage: f.subs({x:y+1})
            0
            sage: f.subs(x=y)
            1
            sage: f.subs(x=x)
            x + y + 1
            sage: f.subs({x:z})
            y + z + 1
            sage: f.subs(x=z+1)
            y + z

            sage: f.subs(x=1/y)
            (y^2 + y + 1)/y
            sage: f.subs({x:1/y})
            (y^2 + y + 1)/y

        TESTS:
            sage: P.<x,y,z> = QQ[]
            sage: f = y
            sage: f.subs({y:x}).subs({x:z})
            z

        NOTE: The evaluation is performed by evalutating every
        variable:value pair separately.  This has side effects if
        e.g. x=y, y=z is provided. If x=y is evaluated first, all x
        variables will be replaced by z eventually.

        """
        cdef int mi, i, need_map, try_symbolic

        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>self._parent
        cdef ring *_ring = parent._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        cdef poly *_p = p_Copy(self._poly, _ring)
        cdef poly *_f

        cdef ideal *to_id = idInit(_ring.N,1)
        cdef ideal *from_id
        cdef ideal *res_id
        need_map = 0
        try_symbolic = 0

        if fixed is not None:
            for m,v in fixed.iteritems():
                if PY_TYPE_CHECK(m,int) or PY_TYPE_CHECK(m,Integer):
                    mi = m+1
                elif PY_TYPE_CHECK(m,MPolynomial_libsingular) and <MPolynomialRing_libsingular>m.parent() is parent:
                    for i from 0 < i <= _ring.N:
                        if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                            mi = i
                            break
                    if i > _ring.N:
                        raise TypeError, "key does not match"
                else:
                    raise TypeError, "keys do not match self's parent"
                try:
                    v = parent._coerce_c(v)
                except TypeError:
                    try_symbolic = 1
                    break
                _f = (<MPolynomial_libsingular>v)._poly
                if _f == NULL or pNext(_f) == NULL:
                    _p = pSubst(_p, mi, _f)
                else:
                    need_map = 1
                    to_id.m[mi-1] = p_Copy(_f, _ring)

        if not try_symbolic:
            gd = parent.gens_dict()
            for m,v in kw.iteritems():
                m = gd[m]
                for i from 0 < i <= _ring.N:
                    if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                        mi = i
                        break
                if i > _ring.N:
                    raise TypeError, "key does not match"
                try:
                    v = parent._coerce_c(v)
                except TypeError:
                    try_symbolic = 1
                    break
                _f = (<MPolynomial_libsingular>v)._poly
                if _f == NULL or pNext(_f) == NULL:
                    _p = pSubst(_p, mi, _f)
                else:
                    if to_id.m[mi-1] != NULL:
                        p_Delete(&to_id.m[mi-1],_ring)
                    to_id.m[mi-1] = p_Copy(_f, _ring)
                    need_map = 1

            if need_map:
                for mi from 0 <= mi < _ring.N:
                    if to_id.m[mi] == NULL:
                        to_id.m[mi] = p_ISet(1,_ring)
                        p_SetExp(to_id.m[mi], mi+1, 1, _ring)
                        p_Setm(to_id.m[mi], _ring)

                from_id=idInit(1,1)
                from_id.m[0] = _p

                res_id = fast_map(from_id, _ring, to_id, _ring)
                _p = res_id.m[0]

                from_id.m[0] = NULL
                res_id.m[0] = NULL

                id_Delete(&from_id, _ring)
                id_Delete(&res_id, _ring)

        id_Delete(&to_id, _ring)

        if not try_symbolic:
            return co.new_MP(parent,_p)

        # now as everything else failed, try to do it symbolically as in call

        g = list(parent.gens())

        if fixed is not None:
            for m,v in fixed.iteritems():
                if PY_TYPE_CHECK(m,int) or PY_TYPE_CHECK(m,Integer):
                    mi = m+1
                elif PY_TYPE_CHECK(m,MPolynomial_libsingular) and <MPolynomialRing_libsingular>m.parent() is parent:
                    for i from 0 < i <= _ring.N:
                        if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                            mi = i
                            break
                    if i > _ring.N:
                        raise TypeError, "key does not match"
                else:
                    raise TypeError, "keys do not match self's parent"

                g[mi-1] = v

        for m,v in kw.iteritems():
            m = gd[m]
            for i from 0 < i <= _ring.N:
                if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                    mi = i
                    break
            if i > _ring.N:
                raise TypeError, "key does not match"

            g[mi-1] = v

        return self(*g)

    def monomials(self):
        """
        Return the list of monomials in self. The returned list is
        decreasingly ordered by the term ordering of self.parent().

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: f = x + 3/2*y*z^2 + 2/3
            sage: f.monomials()
            [y*z^2, x, 1]
            sage: f = P(3/2)
            sage: f.monomials()
            [1]

        TESTS:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: f = x
            sage: f.monomials()
            [x]
            sage: f = P(0)
            sage: f.monomials()
            [0]


        """
        l = list()
        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>self._parent
        cdef ring *_ring = parent._ring
        cdef poly *p = p_Copy(self._poly, _ring)
        cdef poly *t

        if p == NULL:
            return [parent._zero_element]

        while p:
            t = pNext(p)
            p.next = NULL
            p_SetCoeff(p, n_Init(1,_ring), _ring)
            p_Setm(p, _ring)
            l.append( co.new_MP(parent,p) )
            p = t

        return l

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate polynomial.

        EXAMPLES:
            sage: P.<x, y> = MPolynomialRing(QQ,2)
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x^2
            sage: f.constant_coefficient()
            0
        """
        cdef poly *p = self._poly
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        if p == NULL:
            return (<MPolynomialRing_libsingular>self._parent)._base._zero_element

        while p.next:
            p = pNext(p)

        if p_LmIsConstant(p, r):
            return co.si2sa( p_GetCoeff(p, r), r, (<MPolynomialRing_libsingular>self._parent)._base )
        else:
            return (<MPolynomialRing_libsingular>self._parent)._base._zero_element

    def is_univariate(self):
        """
        Return True if self is a univariate polynomial, that is if
        self contains only one variable.

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: f = x^2 + 1
            sage: f.is_univariate()
            True
            sage: f = y*x^2 + 1
            sage: f.is_univariate()
            False
            sage: f = P(0)
            sage: f.is_univariate()
            True
        """
        return bool(len(self._variable_indices_(sort=False))<2)

    def univariate_polynomial(self, R=None):
        """
        Returns a univariate polynomial associated to this
        multivariate polynomial.

        INPUT:
            R -- (default: None) PolynomialRing

        If this polynomial is not in at most one variable, then a
        ValueError exception is raised.  This is checked using the
        is_univariate() method.  The new Polynomial is over the same
        base ring as the given MPolynomial and in the variable 'x' if
        no ring 'ring' is provided.

        EXAMPLES:
            sage: R.<x, y> = MPolynomialRing(QQ,2)
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: polynomial must involve at most one variable
            sage: g = f.subs({x:10}); g
            700*y^2 - 2*y + 305
            sage: g.univariate_polynomial ()
            700*x^2 - 2*x + 305
            sage: g.univariate_polynomial(PolynomialRing(QQ,'z'))
            700*z^2 - 2*z + 305
        """
        cdef poly *p = self._poly
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        k = self.base_ring()

        if not self.is_univariate():
            raise TypeError, "polynomial must involve at most one variable"

        #construct ring if none
        if R == None:
            R =  PolynomialRing(k,'x')

        zero = k(0)
        coefficients = [zero] * (self.degree() + 1)

        while p:
            coefficients[pTotaldegree(p, r)] = co.si2sa(p_GetCoeff(p, r), r, k)
            p = pNext(p)

        return R(coefficients)


    def _variable_indices_(self, sort=True):
        """
        Return the indices of all variables occuring in self.
        This index is the index as SAGE uses them (starting at zero), not
        as SINGULAR uses them (starting at one).

        INPUT:
            sort -- specifies whether the indices shall be sorted

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: f = x*z^2 + z + 1
            sage: f._variable_indices_()
            [0, 2]

        """
        cdef poly *p
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef int i
        s = set()
        p = self._poly
        while p:
            for i from 1 <= i <= r.N:
                if p_GetExp(p,i,r):
                    s.add(i-1)
            p = pNext(p)
        if sort:
            return sorted(s)
        else:
            return list(s)

    def variables(self, sort=True):
        """
        Return a list of all variables occuring in self.

        INPUT:
            sort -- specifies whether the indices shall be sorted

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: f = x*z^2 + z + 1
            sage: f.variables()
            [z, x]
            sage: f.variables(sort=False)
            [x, z]

        """
        cdef poly *p, *v
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        if(r != currRing): rChangeCurrRing(r)
        cdef int i
        l = list()
        si = set()
        p = self._poly
        while p:
            for i from 1 <= i <= r.N:
                if i not in si and p_GetExp(p,i,r):
                    v = p_ISet(1,r)
                    p_SetExp(v, i, 1, r)
                    p_Setm(v, r)
                    l.append(co.new_MP(self._parent, v))
                    si.add(i)
            p = pNext(p)
        if sort:
            return sorted(l)
        else:
            return l

    def variable(self, i=0):
        """
        Return the i-th variable occuring in self. The index i is the
        index in self.variables().

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: f = x*z^2 + z + 1
            sage: f.variables()
            [z, x]
            sage: f.variable(1)
            x
        """
        return self.variables()[i]

    def nvariables(self):
        """
        """
        return self._variable_indices_(sort=False)

    def is_constant(self):
        """
        """
        return bool(p_IsConstant(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring))

    cdef int is_constant_c(self):
        return p_IsConstant(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring)

    def lm(MPolynomial_libsingular self):
        """
        Returns the lead monomial of self with respect to the term
        order of self.parent(). In SAGE a monomial is a product of
        variables in some power without a coefficient.

        EXAMPLES:
            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: R.<x,y,z>=MPolynomialRing(GF(7),3,order='lex')
            sage: f = x^1*y^2 + y^3*z^4
            sage: f.lm()
            x*y^2
            sage: f = x^3*y^2*z^4 + x^3*y^2*z^1
            sage: f.lm()
            x^3*y^2*z^4

            sage: R.<x,y,z>=MPolynomialRing(QQ,3,order='deglex')
            sage: f = x^1*y^2*z^3 + x^3*y^2*z^0
            sage: f.lm()
            x*y^2*z^3
            sage: f = x^1*y^2*z^4 + x^1*y^1*z^5
            sage: f.lm()
            x*y^2*z^4

            sage: R.<x,y,z>=PolynomialRing(GF(127),3,order='degrevlex')
            sage: f = x^1*y^5*z^2 + x^4*y^1*z^3
            sage: f.lm()
            x*y^5*z^2
            sage: f = x^4*y^7*z^1 + x^4*y^2*z^3
            sage: f.lm()
            x^4*y^7*z

        """
        cdef poly *_p
        cdef ring *_ring
        _ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if self._poly == NULL:
            return (<MPolynomialRing_libsingular>self._parent)._zero_element
        _p = p_Head(self._poly, _ring)
        p_SetCoeff(_p, n_Init(1,_ring), _ring)
        p_Setm(_p,_ring)
        return co.new_MP((<MPolynomialRing_libsingular>self._parent), _p)


    def lc(MPolynomial_libsingular self):
        """
        Leading coefficient of self. See self.lm() for details.
        """

        cdef poly *_p
        cdef ring *_ring
        cdef number *_n
        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        if self._poly == NULL:
            return (<MPolynomialRing_libsingular>self._parent)._base._zero_element

        if(_ring != currRing): rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        return co.si2sa(_n, _ring, (<MPolynomialRing_libsingular>self._parent)._base)

    def lt(MPolynomial_libsingular self):
        """
        Leading term of self. In SAGE a term is a product of variables
        in some power AND a coefficient.

        See self.lm() for details
        """
        if self._poly == NULL:
            return (<MPolynomialRing_libsingular>self._parent)._zero_element

        return co.new_MP((<MPolynomialRing_libsingular>self._parent),
                                           p_Head(self._poly,(<MPolynomialRing_libsingular>self._parent)._ring))

    def is_zero(self):
        if self._poly is NULL:
            return True
        else:
            return False

    def __nonzero__(self):
        if self._poly:
            return True
        else:
            return False

    def __floordiv__(self, right):
        """
        Perform division with remainder and return the quotient.

        INPUT:
            right -- something coercable to an MPolynomial_libsingular in self.parent()

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(GF(32003),3)
            sage: R.<x,y,z> = MPolynomialRing(GF(32003),3)
            sage: f = y*x^2 + x + 1
            sage: f//x
            x*y + 1
            sage: f//y
            x^2
        """
        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>(<MPolynomial_libsingular>self)._parent
        cdef ring *r = parent._ring
        if(r != currRing): rChangeCurrRing(r)
        cdef MPolynomial_libsingular _self, _right
        cdef poly *quo

        _self = self

        if not PY_TYPE_CHECK(right, MPolynomial_libsingular) or (<ParentWithBase>parent is not (<MPolynomial_libsingular>right)._parent):
            _right = parent._coerce_c(right)
        else:
            _right = right

        if right.is_zero():
            raise ZeroDivisionError

        cdef int count = polyLengthBounded(_self._poly,15)
        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
            _sig_on
        quo = singclap_pdivide( _self._poly, _right._poly )
        if count >= 15:
            _sig_off
        return co.new_MP(parent, quo)

    def factor(self):
        """
        Return the factorization of self.

        EXAMPLE:
            sage: R.<x,y,z> = PolynomialRing(GF(32003),3)
            sage: R.<x,y,z> = MPolynomialRing(GF(32003),3)
            sage: f = 9*(x-1)^2*(y+z)
            sage: f.factor()
            (9) * (y + z) * (x - 1)^2

            sage: R.<x,w,v,u> = QQ['x','w','v','u']
            sage: p = (4*v^4*u^2 - 16*v^2*u^4 + 16*u^6 - 4*v^4*u + 8*v^2*u^3 + v^4)
            sage: p.factor()
            (-2*v^2*u + 4*u^3 + v^2)^2
            sage: R.<a,b,c,d> = QQ[]
            sage: f =  (-2) * (a - d) * (-a + b) * (b - d) * (a - c) * (b - c) * (c - d)
            sage: F = f.factor(); F
            (-2) * (c - d) * (b - d) * (b - c) * (-a + b) * (a - d) * (a - c)
            sage: F[0][0]
            c - d
            sage: F.unit_part()
            -2

        Factorization of multivariate polynomials over non-prime
        finite fields is only implemented in Singular, and
        unfortunately Singular is currently very buggy at this
        computation.  So we disable it in Sage:

            sage: k.<a> = GF(9)
            sage: R.<x,y> = PolynomialRing(k)
            sage: f = (x-a)*(y-a)
            sage: f.factor()
            Traceback (most recent call last):
            ...
            NotImplementedError: factorization of multivariate polynomials over non-prime fields explicitly disabled due to bugs in Singular
        """
        cdef ring *_ring
        cdef poly *ptemp
        cdef intvec *iv
        cdef int *ivv
        cdef ideal *I
        cdef MPolynomialRing_libsingular parent
        cdef int i

        if self.base_ring().is_finite() and not self.base_ring().is_prime_field():
            raise NotImplementedError, "factorization of multivariate polynomials over non-prime fields explicitly disabled due to bugs in Singular"

        parent = self._parent
        _ring = parent._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        # I make a temporary copy of the poly in self because singclap_factorize appears to modify it's parameter
        ptemp = p_Copy(self._poly,_ring)
        iv = NULL
        cdef int count = polyLengthBounded(self._poly,5)
        if count >= 5:
            _sig_on
        I = singclap_factorize ( ptemp, &iv , 0) #delete iv at some point
        if count >= 5:
            _sig_off

        ivv = iv.ivGetVec()
        v = [(co.new_MP(parent, p_Copy(I.m[i],_ring)) , ivv[i])   for i in range(1,I.ncols)]

        unit = co.new_MP(parent, p_Copy(I.m[0],_ring))

        F = Factorization(v,unit)
        F.sort()

        delete(iv)
        id_Delete(&I,_ring)
        p_Delete(&ptemp,_ring)

        return F

    def lift(self, I):
        """
        given an ideal I = (f_1,...,f_r) and some g (== self) in I,
        find s_1,...,s_r such that g = s_1 f_1 + ... + s_r f_r

        EXAMPLE:
            sage: A.<x,y> = PolynomialRing(QQ,2,order='degrevlex')
            sage: I = A.ideal([x^10 + x^9*y^2, y^8 - x^2*y^7 ])
            sage: f = x*y^13 + y^12
            sage: M = f.lift(I)
            sage: M
            [y^7, x^7*y^2 + x^8 + x^5*y^3 + x^6*y + x^3*y^4 + x^4*y^2 + x*y^5 + x^2*y^3 + y^4]
            sage: sum( map( mul , zip( M, I.gens() ) ) ) == f
            True
        """

        cdef ideal *fI = idInit(1,1)
        cdef ideal *_I
        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>self._parent
        cdef int i = 0
        cdef int j
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef ideal *res

        if PY_TYPE_CHECK(I, MPolynomialIdeal):
            I = I.gens()

        _I = idInit(len(I),1)

        for f in I:
            if not (PY_TYPE_CHECK(f,MPolynomial_libsingular) \
                    and <MPolynomialRing_libsingular>(<MPolynomial_libsingular>f)._parent is parent):
                try:
                    f = parent._coerce_c(f)
                except TypeError, msg:
                    id_Delete(&fI,r)
                    id_Delete(&_I,r)
                    raise TypeError, msg

            _I.m[i] = p_Copy((<MPolynomial_libsingular>f)._poly, r)
            i+=1

        fI.m[0]= p_Copy(self._poly, r)

        res = idLift(_I, fI, NULL, 0, 0, 0)
        l = []
        for i from 0 <= i < IDELEMS(res):
            for j from 1 <= j <= IDELEMS(_I):
                l.append( co.new_MP(parent, pTakeOutComp1(&res.m[i], j)) )

        id_Delete(&fI, r)
        id_Delete(&_I, r)
        id_Delete(&res, r)
        return Sequence(l, check=False, immutable=True)

    def reduce(self,I):
        """
        Return the normal form of self w.r.t. I, i.e. return the
        remainder of self with respect to the polynomials in I. If the
        polynomial set/list I is not a Groebner basis the result is
        not canonical.

        INPUT:
            I -- a list/set of polynomials in self.parent(). If I is an ideal,
                 the generators are used.

        EXAMPLE:

            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: f1 = -2 * x^2 + x^3
            sage: f2 = -2 * y + x* y
            sage: f3 = -x^2 + y^2
            sage: F = Ideal([f1,f2,f3])
            sage: g = x*y - 3*x*y^2
            sage: g.reduce(F)
            -6*y^2 + 2*y
            sage: g.reduce(F.gens())
            -6*y^2 + 2*y

        """
        cdef ideal *_I
        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>self._parent
        cdef int i = 0
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef poly *res

        if(r != currRing): rChangeCurrRing(r)

        if PY_TYPE_CHECK(I, MPolynomialIdeal):
            I = I.gens()

        _I = idInit(len(I),1)
        for f in I:
            if not (PY_TYPE_CHECK(f,MPolynomial_libsingular) \
                   and <MPolynomialRing_libsingular>(<MPolynomial_libsingular>f)._parent is parent):
                try:
                    f = parent._coerce_c(f)
                except TypeError, msg:
                    id_Delete(&_I,r)
                    raise TypeError, msg

            _I.m[i] = p_Copy((<MPolynomial_libsingular>f)._poly, r)
            i+=1

        #the second parameter would be qring!
        res = kNF(_I, NULL, p_Copy(self._poly, r))
        return co.new_MP(parent,res)

    def gcd(self, right, algorithm=None):
        """
        Return the greates common divisor of self and right.

        INPUT:
            right -- polynomial
            algorithm -- 'ezgcd' -- EZGCD algorithm
                         'modular' -- multi-modular algorithm (default)

        EXAMPLES:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: f = (x*y*z)^6 - 1
            sage: g = (x*y*z)^4 - 1
            sage: f.gcd(g)
            x^2*y^2*z^2 - 1
            sage: GCD([x^3 - 3*x + 2, x^4 - 1, x^6 -1])
            x - 1

        TESTS:
            sage: Q.<x,y,z> = MPolynomialRing(QQ,3)
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: P(0).gcd(Q(0))
            0
            sage: x.gcd(1)
            1

        """
        cdef MPolynomial_libsingular _right
        cdef poly *_res
        cdef ring *_ring

        if algorithm is None:
            algorithm = "modular"

        if algorithm == "ezgcd":
            Off(SW_USE_CHINREM_GCD)
            On(SW_USE_EZGCD)
        elif algorithm == "modular":
            On(SW_USE_CHINREM_GCD)
            Off(SW_USE_EZGCD)
        else:
            raise TypeError, "algorithm %s not supported"%(algorithm)

        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        if not PY_TYPE_CHECK(right, MPolynomial_libsingular):
            _right = (<MPolynomialRing_libsingular>self._parent)._coerce_c(right)
        else:
            _right = (<MPolynomial_libsingular>right)

        cdef int count = polyLengthBounded(self._poly,20)+polyLengthBounded(_right._poly,20)
        if count >= 20:
            _sig_on
        _res = singclap_gcd(p_Copy(self._poly, _ring), p_Copy(_right._poly, _ring))
        if count >= 20:
            _sig_off

        return co.new_MP((<MPolynomialRing_libsingular>self._parent), _res)

    def lcm(self, MPolynomial_libsingular g):
        """
        Return the least common multiple of self and g.

        INPUT:
            g -- polynomial

        OUTPUT:
            polynomial

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(QQ,3)
            sage: p = (x+y)*(y+z)
            sage: q = (z^4+2)*(y+z)
            sage: lcm(p,q)
            x*y*z^4 + y^2*z^4 + x*z^5 + y*z^5 + 2*x*y + 2*y^2 + 2*x*z + 2*y*z

        NOTE: This only works for GF(p) and QQ as base rings
        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef poly *ret, *prod, *gcd
        cdef MPolynomial_libsingular _g
        if(_ring != currRing): rChangeCurrRing(_ring)

        if self._parent is not g._parent:
            _g = (<MPolynomialRing_libsingular>self._parent)._coerce_c(g)
        else:
            _g = <MPolynomial_libsingular>g

        cdef int count = polyLengthBounded(self._poly,20)+polyLengthBounded(_g._poly,20)
        if count >= 20:
            _sig_on
        gcd = singclap_gcd(p_Copy(self._poly, _ring), p_Copy(_g._poly, _ring))
        prod = pp_Mult_qq(self._poly, _g._poly, _ring)
        ret = singclap_pdivide(prod , gcd )
        p_Delete(&prod, _ring)
        p_Delete(&gcd, _ring)
        if count >= 20:
            _sig_off
        return co.new_MP(self._parent, ret)

    def is_squarefree(self):
        """
        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if(_ring != currRing): rChangeCurrRing(_ring)
        return bool(singclap_isSqrFree(self._poly))

    def quo_rem(self, MPolynomial_libsingular right):
        """
        Returns quotient and remainder of self and right.

        EXAMPLES:
            sage: R.<x,y> = MPolynomialRing(QQ,2)
            sage: f = y*x^2 + x + 1
            sage: f.quo_rem(x)
            (x*y + 1, 1)
            sage: f.quo_rem(y)
            (x^2, x + 1)

        TESTS:
            sage: R.<x,y> = MPolynomialRing(QQ,2)
            sage: R(0).quo_rem(R(1))
            (0, 0)
            sage: R(1).quo_rem(R(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        """
        #cdef ideal *selfI, *rightI, *R, *res
        cdef poly *quo, *rem
        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>self._parent
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring
        if(r != currRing): rChangeCurrRing(r)

        if self._parent is not right._parent:
            right = self._parent._coerce_c(right)

        if right.is_zero():
            raise ZeroDivisionError

        cdef int count = polyLengthBounded(self._poly,15)
        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
            _sig_on
        quo = singclap_pdivide( self._poly, right._poly )
        rem = p_Add_q(p_Copy(self._poly, r), p_Neg(pp_Mult_qq(right._poly, quo, r), r), r)
        if count >= 15:
            _sig_off
        return co.new_MP(parent, quo), co.new_MP(parent, rem)

    def _magma_(self, magma=None):
        """
        Returns the MAGMA representation of self.

        EXAMPLES:
            sage: R.<x,y> = MPolynomialRing(GF(2),2)
            sage: f = y*x^2 + x +1
            sage: f._magma_() #optional
            x^2*y + x + 1
        """
        if magma is None:
            # TODO: import this globally
            import sage.interfaces.magma
            magma = sage.interfaces.magma.magma

        magma_gens = [e.name() for e in self.parent()._magma_().gens()]
        f = self._repr_with_changed_varnames(magma_gens)
        return magma(f)

    def _singular_(self, singular=singular_default, have_ring=False):
        """
        Return a SINGULAR (as in the CAS) element for this
        element. The result is cached.

        INPUT:
            singular -- interpreter (default: singular_default)
            have_ring -- should the correct ring not be set in SINGULAR first (default:False)

        EXAMPLES:
            sage: P.<x,y,z> = PolynomialRing(GF(127),3)
            sage: x._singular_()
            x
            sage: f =(x^2 + 35*y + 128); f
            x^2 + 35*y + 1
            sage: x._singular_().name() == x._singular_().name()
            True


        TESTS:
            sage: P.<x,y,z> = MPolynomialRing(GF(127),3)
            sage: P.<x,y,z> = PolynomialRing(GF(127),3)
            sage: P(0)._singular_()
            0

        """
        if not have_ring:
            self.parent()._singular_(singular).set_ring() #this is expensive

        try:
            if self.__singular is None:
                return self._singular_init_c(singular, True)

            self.__singular._check_valid()

            if self.__singular.parent() is singular:
                return self.__singular

        except (AttributeError, ValueError):
            pass

        return self._singular_init_c(singular, True)

    def _singular_init_(self,singular=singular_default, have_ring=False):
        """
        Return a new SINGULAR (as in the CAS) element for this element.

        INPUT:
            singular -- interpreter (default: singular_default)
            have_ring -- should the correct ring not be set in SINGULAR first (default:False)

        EXAMPLES:
            sage: P.<x,y,z> = PolynomialRing(GF(127),3)
            sage: x._singular_init_()
            x
            sage: (x^2+37*y+128)._singular_init_()
            x^2+37*y+1
            sage: x._singular_init_().name() == x._singular_init_().name()
            False

        TESTS:
            sage: P(0)._singular_init_()
            0
        """
        return self._singular_init_c(singular, have_ring)

    cdef _singular_init_c(self,singular, have_ring):
        """
        See MPolynomial_libsingular._singular_init_

        """
        if not have_ring:
            self.parent()._singular_(singular).set_ring() #this is expensive

        self.__singular = singular(str(self))
        return self.__singular

    def sub_m_mul_q(self, MPolynomial_libsingular m, MPolynomial_libsingular q):
        """
        Return self - m*q, where m must be a monomial and q a
        polynomial.

        INPUT:
            m -- a monomial
            q -- a polynomial

        EXAMPLE:
            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: x.sub_m_mul_q(y,z)
            -y*z + x

        TESTS:
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: Q.<x,y,z>=MPolynomialRing(QQ,3)
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: P(0).sub_m_mul_q(P(0),P(1))
            0
            sage: x.sub_m_mul_q(Q.gen(1),Q.gen(2))
            -y*z + x

         """
        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring

        if not self._parent is m._parent:
            m = self._parent._coerce_c(m)
        if not self._parent is q._parent:
            q = self._parent._coerce_c(q)

        if m._poly and m._poly.next:
            raise ArithmeticError, "m must be a monomial"
        elif not m._poly:
            return self

        return co.new_MP(self._parent, p_Minus_mm_Mult_qq(p_Copy(self._poly, r), m._poly, q._poly, r))

    def add_m_mul_q(self, MPolynomial_libsingular m, MPolynomial_libsingular q):
        """
        Return self + m*q, where m must be a monomial and q a
        polynomial.

       INPUT:
            m -- a monomial
            q -- a polynomial

        EXAMPLE:
            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: x.add_m_mul_q(y,z)
            y*z + x

        TESTS:
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R.<x,y,z>=MPolynomialRing(QQ,3)
            sage: P.<x,y,z>=MPolynomialRing(QQ,3)
            sage: P(0).add_m_mul_q(P(0),P(1))
            0
            sage: x.add_m_mul_q(R.gen(),R.gen(1))
            x*y + x
         """

        cdef ring *r = (<MPolynomialRing_libsingular>self._parent)._ring

        if not self._parent is m._parent:
            m = self._parent._coerce_c(m)
        if not self._parent is q._parent:
            q = self._parent._coerce_c(q)

        if m._poly and m._poly.next:
            raise ArithmeticError, "m must be a monomial"
        elif not m._poly:
            return self

        return co.new_MP(self._parent, p_Plus_mm_Mult_qq(p_Copy(self._poly, r), m._poly, q._poly, r))


    def __reduce__(self):
        """

        Serialize self.

        EXAMPLES:
            sage: P.<x,y,z> = PolynomialRing(QQ,3, order='degrevlex')
            sage: f = 27/113 * x^2 + y*z + 1/2
            sage: f == loads(dumps(f))
            True

            sage: P = PolynomialRing(GF(127),3,names='abc')
            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1
            sage: f == loads(dumps(f))
            True

        """
        return sage.rings.polynomial.multi_polynomial_libsingular.unpickle_MPolynomial_libsingular, ( self._parent, self.dict() )

    def _im_gens_(self, codomain, im_gens):
        """

        INPUT:
            codomain
            im_gens

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: f = R.hom([y,x], R)
            sage: f(x^2 + 3*y^5)
            3*x^5 + y^2

            sage: R.<a,b,c,d> = QQ[]
            sage: S.<u> = QQ[]
            sage: h = R.hom([0,0,0,u], S)
            sage: h((a+d)^3)
            u^3
        """
        #TODO: very slow
        n = self.parent().ngens()
        if n == 0:
            return codomain._coerce_(self)
        y = codomain(0)
        for (m,c) in self.dict().iteritems():
            y += codomain(c)*mul([ im_gens[i]**m[i] for i in range(n) if m[i]])
        return y


    def _derivative(self, MPolynomial_libsingular var):
        """
        Differentiates self with respect to the provided variable. This
        is completely symbolic so it is also defined over e.g. finite
        fields.

        INPUT:
            variable -- the derivative is taken with respect to variable
            have_ring -- ignored, accepted for compatibility reasons

        SEE ALSO:
            self.derivative()

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
            sage: f._derivative(x)
            9*x^2*y^2 + 3
            sage: f._derivative(y)
            6*x^3*y + 10*y

        The derivative is also defined over finite fields:

            sage: R.<x,y> = PolynomialRing(GF(2**8, 'a'),2)
            sage: f = x^3*y^2 + y^2 + x + 2
            sage: f._derivative(x)
            x^2*y^2 + 1

        """
        if var is None:
            raise ValueError, "you must specify which variable with respect to which to differentiate"

        cdef int i, var_i

        cdef poly *p
        if var._parent is not self._parent:
            raise TypeError, "provided variable is not in same ring as self"
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if _ring != currRing:
            rChangeCurrRing(_ring)

        var_i = -1
        for i from 0 <= i <= _ring.N:
            if p_GetExp(var._poly, i, _ring):
                if var_i == -1:
                    var_i = i
                else:
                    raise TypeError, "provided variable is not univariate"

        if var_i == -1:
            raise TypeError, "provided variable is constant"

        p = pDiff(self._poly, var_i)
        return co.new_MP(self._parent,p)


    def resultant(self, MPolynomial_libsingular other, variable=None):
        """
        computes the resultant of self and the first argument with
        respect to the variable given as the second argument.

        If a second argument is not provide the first variable of
        self.parent() is chosen.

        INPUT:
            other -- polynomial in self.parent()
            variable -- optional variable (of type polynomial) in self.parent() (default: None)

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: a = x+y
            sage: b = x^3-y^3
            sage: c = a.resultant(b); c
            -2*y^3
            sage: d = a.resultant(b,y); d
            2*x^3

            The SINGULAR example:
            sage: R.<x,y,z> = PolynomialRing(GF(32003),3)
            sage: f = 3 * (x+2)^3 + y
            sage: g = x+y+z
            sage: f.resultant(g,x)
            3*y^3 + 9*y^2*z + 9*y*z^2 + 3*z^3 - 18*y^2 - 36*y*z - 18*z^2 + 35*y + 36*z - 24

        TESTS:
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x,y> = MPolynomialRing_libsingular(QQ,2,order='degrevlex')
            sage: a = x+y
            sage: b = x^3-y^3
            sage: c = a.resultant(b); c
            -2*y^3
            sage: d = a.resultant(b,y); d
            2*x^3

        """
        cdef ring *_ring = (<MPolynomialRing_libsingular>self._parent)._ring
        cdef poly *rt

        if variable is None:
            variable = self.parent().gen(0)

        if not self._parent is other._parent:
            raise TypeError, "first parameter needs to be an element of self.parent()"

        if not variable.parent() is self.parent():
            raise TypeError, "second parameter needs to be an element of self.parent() or None"

        cdef int count = polyLengthBounded(self._poly,20)+polyLengthBounded(other._poly,20)
        if count >= 20:
            _sig_on
        rt =  singclap_resultant(self._poly, other._poly, (<MPolynomial_libsingular>variable)._poly )
        if count >= 20:
            _sig_off
        return co.new_MP(self._parent, rt)

    def coefficients(self):
        """
        Return the nonzero coefficients of this polynomial in a list.
        The returned list is decreasingly ordered by the term ordering
        of self.parent().

        EXAMPLES:
            sage: R.<x,y,z> = MPolynomialRing(QQ,3,order='degrevlex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]

            sage: R.<x,y,z> = MPolynomialRing(QQ,3,order='lex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        AUTHOR:
            -- didier deshommes
        """
        cdef poly *p
        cdef ring *r
        r = (<MPolynomialRing_libsingular>self._parent)._ring
        if r!=currRing: rChangeCurrRing(r)
        base = (<MPolynomialRing_libsingular>self._parent)._base
        p = self._poly
        coeffs = list()
        while p:
            coeffs.append(co.si2sa(p_GetCoeff(p, r), r, base))
            p = pNext(p)
        return coeffs

    def gradient(self):
        """
        Return a list of partial derivatives of self, ordered by the
        variables of self.parent().

        EXAMPLE:
           sage: P.<x,y,z> = PolynomialRing(QQ,3)
           sage: f= x*y + 1
           sage: f.gradient()
           [y, x, 0]
        """
        cdef ring *r
        cdef int k

        r = (<MPolynomialRing_libsingular>self._parent)._ring
        if r!=currRing: rChangeCurrRing(r)
        i = []
        for k from 0 < k <= r.N:
            i.append( co.new_MP(self._parent, pDiff(self._poly, k)))

        return i

def unpickle_MPolynomial_libsingular(MPolynomialRing_libsingular R, d):
    """
    Deserialize a MPolynomial_libsingular object

    INPUT:
        R -- the base ring
        d -- a Python dictionary as returned by MPolynomial_libsingular.dict

    """
    cdef ring *r = R._ring
    cdef poly *m, *p
    cdef int _i, _e
    p = p_ISet(0,r)
    rChangeCurrRing(r)
    for mon,c in d.iteritems():
        m = p_Init(r)
        for i,e in mon.sparse_iter():
            _i = i
            if _i >= r.N:
                p_Delete(&p,r)
                p_Delete(&m,r)
                raise TypeError, "variable index too big"
            _e = e
            if _e <= 0:
                p_Delete(&p,r)
                p_Delete(&m,r)
                raise TypeError, "exponent too small"
            p_SetExp(m, _i+1,_e, r)
        p_SetCoeff(m, co.sa2si(c, r), r)
        p_Setm(m,r)
        p = p_Add_q(p,m,r)
    return co.new_MP(R,p)


cdef poly *addwithcarry(poly *tempvector, poly *maxvector, int pos, ring *_ring):
    if p_GetExp(tempvector, pos, _ring) < p_GetExp(maxvector, pos, _ring):
      p_SetExp(tempvector, pos, p_GetExp(tempvector, pos, _ring)+1, _ring)
    else:
      p_SetExp(tempvector, pos, 0, _ring)
      tempvector = addwithcarry(tempvector, maxvector, pos + 1, _ring)
    p_Setm(tempvector, _ring)
    return tempvector


