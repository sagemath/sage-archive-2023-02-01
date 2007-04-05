# __nodoctest__
"""
Multivariate polynomials over QQ and GF(p) implemented using SINGULAR as backend.

AUTHORS:
    Martin Albrecht <malb@informatik.uni-bremen.de>


TESTS:
    sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
    sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3, order='degrevlex')
    sage: P == loads(dumps(P))
    True

    sage: f = 27/113 * x^2 + y*z + 1/2
    sage: f == loads(dumps(P))
    True

    sage: P = MPolynomialRing_libsingular(GF(127),3,names='abc')
    sage: P == loads(dumps(P))
    True

    sage: a,b,c = P.gens()
    sage: f = 57 * a^2*b + 43 * c + 1
    sage: f == loads(dumps(f))
    True


"""
# We do this as we get a link error for init_csage(). However, on
# obscure plattforms (Windows) we might need to link to csage anyway.

cdef extern from "stdsage.h":
    ctypedef void PyObject
    object PY_NEW(object t)
    int PY_TYPE_CHECK(object o, object t)
    PyObject** FAST_SEQ_UNSAFE(object o)
    void init_csage()

    void  sage_free(void *p)
    void* sage_realloc(void *p, size_t n)
    void* sage_malloc(size_t)

import os
import sage.rings.memory

from sage.libs.singular.singular import Conversion
from sage.libs.singular.singular cimport Conversion

cdef Conversion co
co = Conversion()

from sage.rings.multi_polynomial_ring import singular_name_mapping, TermOrder
from sage.rings.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polydict import ETuple

from sage.rings.rational_field import RationalField
from sage.rings.finite_field import FiniteField_prime_modn

from  sage.rings.rational cimport Rational

from sage.interfaces.singular import singular as singular_default, is_SingularElement, SingularElement
from sage.interfaces.macaulay2 import macaulay2 as macaulay2_default, is_Macaulay2Element
from sage.structure.factorization import Factorization

from complex_field import is_ComplexField
from real_mpfr import is_RealField


from sage.rings.integer_ring import IntegerRing
from sage.structure.element cimport EuclideanDomainElement, \
     RingElement, \
     ModuleElement, \
     Element, \
     CommutativeRingElement

from sage.rings.integer cimport Integer

from sage.structure.parent cimport Parent
from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens

from sage.misc.sage_eval import sage_eval

# shared library loading
cdef extern from "dlfcn.h":
    void *dlopen(char *, long)
    char *dlerror()

cdef extern from "string.h":
    char *strdup(char *s)


cdef init_singular():
    """
    This initializes the Singular library. Right now, this is a hack.

    SINGULAR has a concept of compiled extension modules similar to
    SAGE. For this, the compiled modules need to see the symbols from
    the main programm. However, SINGULAR is a shared library in this
    context these symbols are not known globally. The work around so
    far is to load the library again and to specifiy RTLD_GLOBAL.
    """

    # This is a work around until we found a way to export those
    # symbols without loading the lib again

    cdef void *handle

    lib = os.environ['SAGE_LOCAL']+"/lib/libsingular.so"

    handle = dlopen(lib, 256+1)

    if handle == NULL:
        print dlerror()

    # Load Singular
    siInit(lib)

    # Steal Memory Manager back or weird things may happen
    sage.rings.memory.pmem_malloc()

# call it
init_singular()

cdef class MPolynomialRing_libsingular(MPolynomialRing_generic):
    """
    A multivariate polynomial ring over QQ or GF(p) implemented using SINGULAR.

    """
    def __init__(self, base_ring, n, names, order='degrevlex'):
        """

        Constructs a multivariate polynomial ring subject to the following conditions.

        INPUT:
            base_ring -- base ring (must be either GF(p) (p prime) or QQ)
            n -- number of variables (must be at least 1)
            names -- names of ring variables, may be string of list/tuple
            order -- term order (default: degrevlex)

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P
            Polynomial Ring in x, y, z over Rational Field

            sage: f = 27/113 * x^2 + y*z + 1/2; f
            27/113*x^2 + y*z + 1/2

            sage: P.term_order()
            Degree reverse lexicographic term order

            sage: P = MPolynomialRing_libsingular(GF(127),3,names='abc', order='lex')
            sage: P
            Polynomial Ring in a, b, c over Finite Field of size 127

            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1; f
            57*a^2*b + 43*c + 1

            sage: P.term_order()
            Lexicographic term order

        """
        cdef char **_names
        cdef char *_name
        cdef int i
        cdef int characteristic

        n = int(n)
        if n<1:
            raise ArithmeticError, "number of variables must be at least 1"

        self.__ngens = n

        MPolynomialRing_generic.__init__(self, base_ring, n, names, TermOrder(order))

        self._has_singular = True

        _names = <char**>sage_malloc(sizeof(char*)*len(self._names))

        for i from 0 <= i < n:
            _name = self._names[i]
            _names[i] = strdup(_name)

        if PY_TYPE_CHECK(base_ring, FiniteField_prime_modn):
            characteristic = base_ring.characteristic()

        elif PY_TYPE_CHECK(base_ring, RationalField):
            characteristic = 0

        else:
            raise NotImplementedError, "Only GF(p) and QQ are supported right now, sorry"

        try:
            order = singular_name_mapping[order]
        except KeyError:
            pass

        self._ring = rDefault(characteristic, n, _names)
        if(self._ring != currRing): rChangeCurrRing(self._ring)

        rUnComplete(self._ring)

        omFree(self._ring.wvhdl)
        omFree(self._ring.order)
        omFree(self._ring.block0)
        omFree(self._ring.block1)

        self._ring.wvhdl  = <int **>omAlloc0(3 * sizeof(int*))
        self._ring.order  = <int *>omAlloc0(3* sizeof(int *))
        self._ring.block0 = <int *>omAlloc0(3 * sizeof(int *))
        self._ring.block1 = <int *>omAlloc0(3 * sizeof(int *))

        if order == "dp":
            self._ring.order[0] = ringorder_dp
        elif order == "Dp":
            self._ring.order[0] = ringorder_Dp
        elif order == "lp":
            self._ring.order[0] = ringorder_lp
        elif order == "rp":
            self._ring.order[0] = ringorder_rp
        else:
            self._ring.order[0] = ringorder_lp

        self._ring.order[1] = ringorder_C

        self._ring.block0[0] = 1
        self._ring.block1[0] = n

        rComplete(self._ring, 1)
        self._ring.ShortOut = 0

        sage_free(_names)

    def __dealloc__(self):
        """
        """
        rDelete(self._ring)

    cdef _coerce_c_impl(self, element):
        """

        Coerces elements to self.

        EXAMPLES:
             sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
             sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3)

             We can coerce elements of self to self

             sage: P._coerce_(x*y + 1/2)
             x*y + 1/2

             We can coerce elements for a ring with the same algebraic properties

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

        """
        cdef poly *_p
        cdef ring *_ring
        cdef number *_n
        cdef poly *rm
        cdef int ok
        cdef int i
        cdef ideal *from_id, *to_id, *res_id

        _ring = self._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        if PY_TYPE_CHECK(element, MPolynomial_libsingular):
            if element.parent() is <object>self:
                return element
            elif element.parent() == self:
                # is this safe?
                _p = p_Copy((<MPolynomial_libsingular>element)._poly, _ring)
            else:
                TypeError, "parents do not match"

        elif PY_TYPE_CHECK(element, CommutativeRingElement):
            # Accepting ZZ
            if element.parent() is IntegerRing():
                _p = p_ISet(int(element), _ring)

            elif  <Parent>element.parent() is self._base:
                # Accepting GF(p)
                if PY_TYPE_CHECK(self._base, FiniteField_prime_modn):
                    _p = p_ISet(int(element), _ring)

                # Accepting QQ
                elif PY_TYPE_CHECK(self._base, RationalField):
                    _n = co.sa2si_QQ(element,_ring)
                    _p = p_NSet(_n, _ring)
                else:
                    raise NotImplementedError
            else:
                raise TypeError, "base rings must be identical"

        # Accepting int
        elif PY_TYPE_CHECK(element, int):
            _p = p_ISet(int(element), _ring)
        else:
            raise TypeError, "Cannot coerce element"

        return new_MP(self,_p)

    def __call__(self, element):
        """
        EXAMPLE:
             Call supports all conversions _coerce_ supports, plus:

             sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
             sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
             sage: P('x+y + 1/4')
             x + y + 1/4
             sage: P._singular_()
             //   characteristic : 0
             //   number of vars : 3
             //        block   1 : ordering dp
             //                  : names    x y z
             //        block   2 : ordering C
             sage: P(singular('x + 3/4'))
             x + 3/4
        """
        cdef poly *_m, *_p, *_tmp
        cdef ring *_r

        if PY_TYPE_CHECK(element, SingularElement):
            element = str(element)

        if PY_TYPE_CHECK(element,str):
            # let python do the the parsing
            return sage_eval(element,self.gens_dict())

               # this almost does what I want, besides variables with 0-9 in their names
##             _r = self._ring
##             if(_r != currRing): rChangeCurrRing(_r)
##             # improve this
##             element = element.replace('^','').replace('**','').replace(' ','').strip()
##             monomials = element.split('+')
##             _p = p_ISet(0, _r)
##             # wrong for e.g. x0,..,x7
##             for m in monomials:
##                 _m = p_ISet(1,_r)
##                 for var in m.split('*'):
##                     p_Read(var, _tmp, _r)
##                     _m = p_Mult_q(_m,_tmp, _r)

##                 _p = p_Add_q(_p,_m,_r)

##             return new_MP(self,_p)

        return self._coerce_c_impl(element)

    def _repr_(self):
        varstr = ", ".join([ rRingVar(i,self._ring)  for i in range(self.__ngens) ])
        return "Polynomial Ring in %s over %s"%(varstr,self._base)

    def ngens(self):
        """
        Returns the number of variables in self.

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x,y> = MPolynomialRing_libsingular(QQ, 2)
            sage: P.ngens()
            2
            sage: P = MPolynomialRing_libsingular(GF(127),1000,'x')
            sage: P.ngens()
            1000

        """
        return int(self.__ngens)

    def gens(self):
        """
        Return the tuple of variables in self.

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P.gens()
            (x, y, z)

            sage: P = MPolynomialRing_libsingular(QQ,10,'x')
            sage: P.gens()
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)

            sage: P.<SAGE,SINGULAR> = MPolynomialRing_libsingular(QQ,2) # weird names
            sage: P.gens()
            (SAGE, SINGULAR)

        """
        return tuple([self.gen(i) for i in range(self.__ngens)  ])

    def gen(self, int n=0):
        """
        Returns the n-th generator of self.

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P.gen(),P.gen(1)
            (x, y)

            sage: P = MPolynomialRing_libsingular(GF(127),1000,'x')
            sage: P.gen(500)
            x500

            sage: P.<SAGE,SINGULAR> = MPolynomialRing_libsingular(QQ,2) # weird names
            sage: P.gen(1)
            SINGULAR

        """
        cdef poly *_p
        cdef ring *_ring

        if n < 0 or n >= self.__ngens:
            raise ValueError, "Generator not defined."

        _ring = self._ring
        _p = p_ISet(1,_ring)

        # oddly enough, Singular starts counting a 1!!!
        p_SetExp(_p, n+1, 1, _ring)
        p_Setm(_p, _ring);

        return new_MP(self,_p)

    def ideal(self, gens, coerce=True):
        """
        Create an ideal in this polynomial ring.

        """
        if is_SingularElement(gens):
            gens = list(gens)
            coerce = True
        if is_Macaulay2Element(gens):
            gens = list(gens)
            coerce = True
        elif not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return MPolynomialIdeal(self, gens, coerce=False)

    def _macaulay2_(self, macaulay2=macaulay2_default):
        try:
            R = self.__macaulay2
            if R is None or not (R.parent() is macaulay2):
                raise ValueError
            R._check_valid()
            return R
        except (AttributeError, ValueError):
            if self.base_ring().is_prime_field():
                if self.characteristic() == 0:
                    base_str = "QQ"
                else:
                    base_str = "ZZ/" + str(self.characteristic())
            elif is_IntegerRing(self.base_ring()):
                base_str = "ZZ"
            else:
                raise TypeError, "no conversion of to a Macaulay2 ring defined"
            self.__macaulay2 = macaulay2.ring(base_str, str(self.gens()), \
                                              self.term_order().macaulay2_str())
        return self.__macaulay2

    def _singular_(self, singular=singular_default):
        try:
            R = self.__singular
            if R is None or not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if self.base_ring().is_prime_field():
                return R
            if self.base_ring().is_finite():
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != self.__minpoly:
                    singular.eval("minpoly=%s"%(self.__minpoly))
            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self, singular=singular_default):
        """
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
            self.__minpoly = "("+(str(self.base_ring().modulus()).replace("x",gen)).replace(" ","")+")"
            singular.eval("minpoly=%s"%(self.__minpoly) )

            self.__singular = r
        else:
            raise TypeError, "no conversion to a Singular ring defined"

        return self.__singular

    def __hash__(self):
        return hash(self.__repr__())

    def __richcmp__(left, right, int op):
        return (<Parent>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
        """
        Multivariate polynomial rings are said to be equal if:
         * their base rings match
         * their generator names match
         * their term orderings match
        """
        if PY_TYPE_CHECK(right, MPolynomialRing_libsingular):
            return cmp( (left.base_ring(), map(str, left.gens()), left.term_order()),
                        (right.base_ring(), map(str, right.gens()), right.term_order())
                        )
        else:
            return cmp(type(left),type(right))

    def __reduce__(self):
        """
        Serializes self.

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x,y,z> = MPolynomialRing_libsingular(QQ,3, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P = MPolynomialRing_libsingular(GF(127),3,names='abc')
            sage: P == loads(dumps(P))
            True

        """
        return sage.rings.multi_polynomial_libsingular.unpickle_MPolynomialRing_libsingular, ( self.base_ring(),
                                                                                               map(str, self.gens()),
                                                                                               self.term_order() )


def unpickle_MPolynomialRing_libsingular(base_ring, names, term_order):
    """
    inverse function for MPolynomialRing_libsingular.__reduce__

    """
    return MPolynomialRing_libsingular(base_ring, len(names), names, term_order)

cdef new_MP(MPolynomialRing_libsingular parent, poly *juice):
    """
    Construct a new MPolynomial_libsingular element
    """
    cdef MPolynomial_libsingular p
    p = PY_NEW(MPolynomial_libsingular)
    p._parent = <ParentWithBase>parent
    p._poly = juice
    return p

cdef class MPolynomial_libsingular(sage.rings.multi_polynomial.MPolynomial):
    """

    """
    def __init__(self, MPolynomialRing_libsingular parent):
        self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        if self._poly:
            p_Delete(&self._poly, (<MPolynomialRing_libsingular>self._parent)._ring)

    def __call__(self, *x):
        cdef int l = len(x)
        cdef MPolynomialRing_libsingular parent = (<MPolynomialRing_libsingular>self._parent)
        cdef ring *_ring = parent._ring

        cdef poly *_p

        if l != parent._ring.N:
            raise TypeError, "number of arguments does not match number of variables in parent"

        cdef ideal *to_id = idInit(l,1)

        try:
            for i from 0 <= i < l:
                e = x[i] # TODO: optimize this line
                to_id.m[i]= p_Copy( (<MPolynomial_libsingular>(<MPolynomialRing_libsingular>parent._coerce_c(x[i])))._poly, _ring)

        except TypeError:
            id_Delete(&to_id, _ring)
            raise TypeError, "cannot coerce in arguments"

        cdef ideal *from_id=idInit(1,1)
        from_id.m[0] = p_Copy(self._poly, _ring)

        cdef ideal *res_id = fast_map(from_id, _ring, to_id, _ring)
        cdef poly *res = p_Copy(res_id.m[0], _ring)

        id_Delete(&to_id, _ring)
        id_Delete(&from_id, _ring)
        id_Delete(&res_id, _ring)
        return new_MP(parent, res)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        cdef ring *_r

        _r = (<MPolynomialRing_libsingular>left._parent)._ring
        if(_r != currRing): rChangeCurrRing(_r)
        ret =  p_Cmp( (<MPolynomial_libsingular>left)._poly, (<MPolynomial_libsingular>right)._poly, _r)
        return ret

    cdef ModuleElement _add_c_impl( left, ModuleElement right):
        cdef MPolynomial_libsingular res

        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        _l = p_Copy(left._poly, _ring)
        _r = p_Copy((<MPolynomial_libsingular>right)._poly, _ring)

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p= p_Add_q(_l, _r, _ring)

        return new_MP((<MPolynomialRing_libsingular>left._parent),_p)

    cdef ModuleElement _sub_c_impl( left, ModuleElement right):
        cdef MPolynomial_libsingular res

        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        _l = p_Copy(left._poly, _ring)
        _r = p_Copy((<MPolynomial_libsingular>right)._poly, _ring)

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p= p_Add_q(_l, p_Neg(_r, _ring), _ring)

        return new_MP((<MPolynomialRing_libsingular>left._parent),_p)


    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        cdef number *_n
        cdef ring *_ring
        cdef poly *_p

        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        if PY_TYPE_CHECK((<MPolynomialRing_libsingular>self._parent)._base, FiniteField_prime_modn):
            _n = n_Init(int(left),_ring)

        elif PY_TYPE_CHECK((<MPolynomialRing_libsingular>self._parent)._base, RationalField):
            _n = co.sa2si_QQ(left,_ring)

        _p = pp_Mult_nn(self._poly,_n,_ring)
        n_Delete(&_n, _ring)
        return new_MP((<MPolynomialRing_libsingular>self._parent),_p)

    cdef RingElement  _mul_c_impl(left, RingElement right):
        cdef poly *_l, *_r, *_p
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>left._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p = pp_Mult_qq(left._poly, (<MPolynomial_libsingular>right)._poly, _ring)

        return new_MP(left._parent,_p)

    cdef RingElement  _div_c_impl(left, RingElement right):
        raise NotImplementedError


    def __pow__(MPolynomial_libsingular self,int exp,ignored):
        """
        """
        cdef ring *_ring
        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        cdef poly *_p


        if exp < 0:
            raise ArithmeticError, "Cannot comput negative powers of polynomials"

        if(_ring != currRing): rChangeCurrRing(_ring)
        _p = pPower( p_Copy(self._poly,(<MPolynomialRing_libsingular>self._parent)._ring),exp)
        return new_MP((<MPolynomialRing_libsingular>self._parent),_p)


    def __neg__(self):
        cdef ring *_ring
        _ring = (<MPolynomialRing_libsingular>self._parent)._ring
        if(_ring != currRing): rChangeCurrRing(_ring)

        return new_MP((<MPolynomialRing_libsingular>self._parent),\
                                           p_Neg(p_Copy(self._poly,_ring),_ring))

    def _repr_(self):
        s =  self._repr_short_c()
        return s.replace("+"," + ")

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
        s = p_String(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring, (<MPolynomialRing_libsingular>self._parent)._ring)
        return s

    def _latex(self):
        raise NotImplementedError

    def _repr_with_changed_varnames(self, varnames):
        raise NotImplementedError

    def degree(self, x=None):
        return pDeg(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring)

    def newton_polytope(self):
        raise NotImplementedError

    def total_degree(self):
        """
        Return the total degree of self, which is the
        maximum degree of any monomial in self.

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R.<x,y,z> = MPolynomialRing_libsingular(QQ, 3)
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
        """
        return pTotaldegree(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring)

    def monomial_coefficient(self, mon):
        raise NotImplementedError

    def dict(self):
        cdef poly *p
        cdef ring *r
        cdef int n
        cdef int v
        r = (<MPolynomialRing_libsingular>self._parent)._ring
        base = (<MPolynomialRing_libsingular>self._parent)._base
        p = self._poly
        pd = dict()
        while p:
            d = dict()
            # assuming that SINGULAR stores monomial dense
            for v from 1 <= v <= r.N:
                n = p_GetExp(p,v,r)
                if n!=0:
                    d[v-1] = n

            pd[ETuple(d,r.N)] = co.si2sa(p_GetCoeff(p, r), r, base)

            p = p.next
        return pd

    def __getitem__(self,x):
        raise NotImplementedError

    def coefficient(self, mon):
        raise NotImplementedError

    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES:
           sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
           sage: R.<a,b,c> = MPolynomialRing_libsingular(QQ, 3)
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
            pl.append(tuple(ml))

            p = p.next
        return pl

    def is_unit(self):
        """
        Return True if self is a unit.

        EXAMPLES:
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
        raise NotImplementedError

    def is_homogeneous(self):
        """
        Return True if self is a homogeneous polynomial.

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: P.<x,y> = MPolynomialRing_libsingular(RationalField(), 2)
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

    def homogenize(self, var='h'):
        raise NotImplementedError

    def is_monomial(self):
        return not self._poly.next

    def fix(self, fixed):
        """
        Fixes some given variables in a given multivariate polynomial and
        returns the changed multivariate polynomials. The polynomial
        itself is not affected.  The variable,value pairs for fixing are
        to be provided as dictionary of the form {variable:value}.

        This is a special case of evaluating the polynomial with some of
        the variables constants and the others the original variables, but
        should be much faster if only few variables are to be fixed.

        INPUT:
            fixed -- dict with variable:value pairs

        OUTPUT:
            new MPolynomial

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: x, y = MPolynomialRing_libsingular(QQ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f(5,y)
            25*y^2 + y + 30
            sage: f.fix({x:5})
            25*y^2 + y + 30

        """
        cdef int mi, i

        cdef MPolynomialRing_libsingular parent = <MPolynomialRing_libsingular>self._parent
        cdef ring *_ring = parent._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        cdef poly *_p = p_Copy(self._poly, _ring)

        for m,v in fixed.iteritems():
            if PY_TYPE_CHECK(m,int) or PY_TYPE_CHECK(m,Integer):
                mi = m+1
            elif PY_TYPE_CHECK(m,MPolynomial_libsingular) and <MPolynomialRing_libsingular>m.parent() is parent:
                for i from 0 < i <= _ring.N:
                    mi = p_GetExp((<MPolynomial_libsingular>m)._poly,i, _ring)
                    if mi!=0:
                        break
            else:
                raise TypeError, "keys do not match self's parent"

            _p = pSubst(p_Copy(_p, _ring), mi, (<MPolynomial_libsingular>parent._coerce_c(v))._poly)

        return new_MP(parent,_p)

    def monomials(self):
        raise NotImplementedError

    def constant_coefficent(self):
        raise NotImplementedError

    def is_univariate(self):
        raise NotImplementedError

    def univariate_polynomial(self, R=None):
        raise NotImplementedError

    def _variable_indices_(self):
        raise NotImplementedError

    def variables(self):
        raise NotImplementedError

    def variable(self):
        raise NotImplementedError

    def nvariables(self):
        raise NotImplementedError

    def is_constant(self):
        return bool(p_IsConstant(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring))

    def __hash__(self):
        s = p_String(self._poly, (<MPolynomialRing_libsingular>self._parent)._ring, (<MPolynomialRing_libsingular>self._parent)._ring)
        return hash(s)


    def lm(MPolynomial_libsingular self):
        """
        Returns the lead monomial of self with respect to the term order of
        self.parent().

        EXAMPLES:
             sage: from sage.rings.multi_polynomial_libsingular import MPolynomialRing_libsingular

             sage: R.<x,y,z>=MPolynomialRing_libsingular(GF(7),3,order='lex')
             sage: f = x^1*y^2 + y^3*z^4
             sage: f.lm()
             x*y^2
             sage: f = x^3*y^2*z^4 + x^3*y^2*z^1
             sage: f.lm()
             x^3*y^2*z^4

             sage: R.<x,y,z>=MPolynomialRing_libsingular(QQ,3,order='deglex')
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
        _p = p_Head(self._poly, _ring)
        p_SetCoeff(_p, n_Init(int(1),_ring), _ring)
        p_Setm(_p,_ring)
        return new_MP((<MPolynomialRing_libsingular>self._parent), _p)


    def lc(MPolynomial_libsingular self):
        """
        Leading coefficient of self. See self.lm() for details.
        """

        cdef poly *_p
        cdef ring *_ring
        cdef number *_n
        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        return co.si2sa(_n, _ring, (<MPolynomialRing_libsingular>self._parent)._base)

    def lt(MPolynomial_libsingular self):
        """
        Leading term of self. See self.lm() for details
        """
        return new_MP((<MPolynomialRing_libsingular>self._parent),
                                           p_Head(self._poly,(<MPolynomialRing_libsingular>self._parent)._ring))

    def is_zero(self):
        raise NotImplementedError

    def __floordiv__(self,right):
        raise NotImplementedError

    def factor(self, param=0):

        cdef ring *_ring
        cdef intvec *iv
        cdef int *ivv
        cdef ideal *I
        cdef MPolynomialRing_libsingular parent
        cdef int i

        parent = self._parent
        _ring = parent._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        iv = NULL
        I = singclap_factorize ( self._poly, &iv , int(param)) #delete iv at some point

        if param!=1:
            ivv = iv.ivGetVec()
            v = [(new_MP(parent, p_Copy(I.m[i],_ring)) , ivv[i])   for i in range(I.ncols)]
        else:
            v = [(new_MP(parent, p_Copy(I.m[i],_ring)) , 1)   for i in range(I.ncols)]

        # TODO: peel of 1

        F = Factorization(v)
        F.sort()

        # delete intvec
        # delete ideal

        return F

    def lift(self):
        raise NotImplementedError

    def gcd(self, right):

        cdef MPolynomial_libsingular _right
        cdef poly *_res
        cdef ring *_ring

        _ring = (<MPolynomialRing_libsingular>self._parent)._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        if not PY_TYPE_CHECK(right, MPolynomial_libsingular):
            _right = self.parent()._coerce_c_impl(right)
        else:
            _right = (<MPolynomial_libsingular>right)

        _res = singclap_gcd(p_Copy(self._poly, _ring), p_Copy(_right._poly, _ring))

        return new_MP((<MPolynomialRing_libsingular>self._parent), _res)

    def quo_rem(self, right):
        raise NotImplementedError

    def _magma_(self):
        raise NotImplementedError

    def _singular_(self, singular=singular_default, have_ring=False):

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
        return self._singular_init_c(singular, have_ring)

    cdef _singular_init_c(self,singular, have_ring):
        if not have_ring:
            self.parent()._singular_(singular).set_ring() #this is expensive

        self.__singular = singular(str(self))
        return self.__singular


