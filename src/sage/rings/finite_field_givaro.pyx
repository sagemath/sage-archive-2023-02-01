r"""
Finite Extension Fields of cardinality up to $2^{16}$

SAGE includes the Givaro finite field library, for highly optimized
arithmetic in finite fields.

NOTES: The arithmetic is performed by the Givaro C++ library which
uses Zech logs internally to represent finite field elements. This
implementation is the default finite extension field implementation in
SAGE for the cardinality $< 2^{16}$, as it is vastly faster than the
PARI implementation which uses polynomials to represent finite field
elements. Some functionality in this class however is implemented
using the PARI implementation.

EXAMPLES:
    sage: k = GF(5); type(k)
    <class 'sage.rings.finite_field_prime_modn.FiniteField_prime_modn'>
    sage: k = GF(5^2,'c'); type(k)
    <type 'sage.rings.finite_field_givaro.FiniteField_givaro'>
    sage: k = GF(2^16,'c'); type(k)
    <type 'sage.rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e'>
    sage: k = GF(3^16,'c'); type(k)
    <class 'sage.rings.finite_field_ext_pari.FiniteField_ext_pari'>

    sage: n = previous_prime_power(2^16 - 1)
    sage: while is_prime(n):
    ...    n = previous_prime_power(n)
    sage: factor(n)
    251^2
    sage: k = GF(n,'c'); type(k)
    <type 'sage.rings.finite_field_givaro.FiniteField_givaro'>

AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-06-05)
     -- William Stein (2006-12-07): editing, lots of docs, etc.
     -- Robert Bradshaw (2007-05-23): is_square/sqrt, pow.
"""


#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../ext/interrupt.pxi"

# this fails with a scoping error:
#   AttributeError: CVoidType instance has no attribute 'scope'
#include "../ext/stdsage.pxi"

cdef extern from "stdsage.h":
    ctypedef void PyObject
    object PY_NEW(object t)
    int PY_TYPE_CHECK(object o, object t)
    void init_csage()
#init_csage()

from sage.misc.randstate cimport randstate, current_randstate
from sage.rings.ring cimport FiniteField
from sage.rings.ring cimport Ring
from sage.structure.element cimport FiniteFieldElement, Element, RingElement, ModuleElement
from sage.rings.finite_field_element import FiniteField_ext_pariElement
from sage.structure.sage_object cimport SageObject
import operator
import sage.rings.arith
import finite_field
import finite_field_ext_pari

import sage.interfaces.gap
from sage.libs.pari.all import pari
from sage.libs.pari.gen import gen

from sage.structure.parent  cimport Parent
from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens

cdef object is_IntegerMod
cdef object IntegerModRing_generic
cdef object Integer
cdef object Rational
cdef object is_Polynomial
cdef object ConwayPolynomials
cdef object conway_polynomial
cdef object MPolynomial
cdef object Polynomial
cdef object FreeModuleElement

cdef void late_import():
    """
    """
    global is_IntegerMod, \
           IntegerModRing_generic, \
           Integer, \
           Rational, \
           is_Polynomial, \
           ConwayPolynomials, \
           conway_polynomial, \
           MPolynomial, \
           Polynomial, \
           FreeModuleElement

    if is_IntegerMod is not None:
        return

    import sage.rings.integer_mod
    is_IntegerMod = sage.rings.integer_mod.is_IntegerMod

    import sage.rings.integer_mod_ring
    IntegerModRing_generic = sage.rings.integer_mod_ring.IntegerModRing_generic

    import sage.rings.integer
    Integer = sage.rings.integer.Integer

    import sage.rings.rational
    Rational = sage.rings.rational.Rational

    import sage.rings.polynomial.polynomial_element
    is_Polynomial = sage.rings.polynomial.polynomial_element.is_Polynomial

    import sage.databases.conway
    ConwayPolynomials = sage.databases.conway.ConwayPolynomials

    import sage.rings.finite_field
    conway_polynomial = sage.rings.finite_field.conway_polynomial

    import sage.rings.polynomial.multi_polynomial_element
    MPolynomial = sage.rings.polynomial.multi_polynomial_element.MPolynomial

    import sage.rings.polynomial.polynomial_element
    Polynomial = sage.rings.polynomial.polynomial_element.Polynomial

    import sage.modules.free_module_element
    FreeModuleElement = sage.modules.free_module_element.FreeModuleElement

cdef FiniteField_givaro parent_object(Element o):
    return <FiniteField_givaro>(o._parent)

cdef class FiniteField_givaro(FiniteField):
    def __init__(FiniteField_givaro self, q, name="a",  modulus=None, repr="poly", cache=False):
        """
        Finite Field. These are implemented using Zech logs and the
        cardinality must be < 2^16. By default conway polynomials are
        used as minimal polynomial.

        INPUT:
            q     -- p^n (must be prime power)
            name  -- variable used for poly_repr (default: 'a')
            modulus -- you may provide a minimal polynomial to use for
                     reduction or 'random' to force a random
                     irreducible polynomial. (default: None, a conway
                     polynomial is used if found. Otherwise a random
                     polynomial is used)
            repr  -- controls the way elements are printed to the user:
                     (default: 'poly')
                     'log': repr is element.log_repr()
                     'int': repr is element.int_repr()
                     'poly': repr is element.poly_repr()
            cache -- if True a cache of all elements of this field is
                     created. Thus, arithmetic does not create new
                     elements which speeds calculations up. Also, if
                     many elements are needed during a calculation
                     this cache reduces the memory requirement as at
                     most self.order() elements are created. (default: False)

        OUTPUT:
            Givaro finite field with characteristic p and cardinality p^n.

        EXAMPLES:

            By default conway polynomials are used:

            sage: k.<a> = GF(2**8)
            sage: -a ^ k.degree()
            a^4 + a^3 + a^2 + 1
            sage: f = k.modulus(); f
            x^8 + x^4 + x^3 + x^2 + 1

            You may enforce a modulus:

            sage: P.<x> = PolynomialRing(GF(2))
            sage: f = x^8 + x^4 + x^3 + x + 1 # Rijndael Polynomial
            sage: k.<a> = GF(2^8, modulus=f)
            sage: k.modulus()
            x^8 + x^4 + x^3 + x + 1
            sage: a^(2^8)
            a

            You may enforce a random modulus:

            sage: k = GF(3**5, 'a', modulus='random')
            sage: k.modulus() # random polynomial
            x^5 + 2*x^4 + 2*x^3 + x^2 + 2

            Three different representations are possible:

            sage: sage.rings.finite_field_givaro.FiniteField_givaro(9,repr='poly').gen()
            a
            sage: sage.rings.finite_field_givaro.FiniteField_givaro(9,repr='int').gen()
            3
            sage: sage.rings.finite_field_givaro.FiniteField_givaro(9,repr='log').gen()
            5
        """

        # we are calling late_import here because this constructor is
        # called at least once before any arithmetic is performed.
        late_import()

        cdef intvec cPoly

        self._kwargs = {}

        if repr=='poly':
            self.repr = 0
        elif repr=='log':
            self.repr = 1
        elif repr=='int':
            self.repr = 2
        else:
            raise ValueError, "Unknown representation %s"%repr

        self._kwargs['repr'] = repr
        self._kwargs['cache'] = cache

        if q >= 1<<16:
            raise ValueError, "q must be < 2^16"

        q = Integer(q)
        if q < 2:
            raise ValueError, "q  must be a prime power"
        F = q.factor()
        if len(F) > 1:
            raise ValueError, "q must be a prime power"
        p = F[0][0]
        k = F[0][1]

        ParentWithGens.__init__(self, finite_field.FiniteField(p), name, normalize=False)

        self._is_conway = False
        if modulus is None or modulus=="random":
            if k>1 and ConwayPolynomials().has_polynomial(p, k) and modulus!="random":
                modulus = conway_polynomial(p, k)
                self._is_conway = True
            else:
                _sig_on
                self.objectptr = gfq_factorypk(p,k)
                self._zero_element = make_FiniteField_givaroElement(self,self.objectptr.zero)
                self._one_element = make_FiniteField_givaroElement(self,self.objectptr.one)
                _sig_off
                if cache:
                    self._array = self.gen_array()
                return

        if is_Polynomial(modulus):
            modulus = modulus.list()

        if PY_TYPE_CHECK(modulus, list) or PY_TYPE_CHECK(modulus, tuple):

            for i from 0 <= i < len(modulus):
                cPoly.push_back(int( modulus[i] % p ))

            _sig_on
            self.objectptr = gfq_factorypkp(p, k,cPoly)
            _sig_off
            self._zero_element = make_FiniteField_givaroElement(self,self.objectptr.zero)
            self._one_element = make_FiniteField_givaroElement(self,self.objectptr.one)

            if cache:
                self._array = self.gen_array()
            return

        raise TypeError, "Cannot understand modulus"

    cdef gen_array(FiniteField_givaro self):
        """
        Generates an array/list/tuple containing all elements of self indexed by their
        power with respect to the internal generator.
        """
        cdef int i

        array = list()
        for i from 0 <= i < self.order_c():
            array.append( make_FiniteField_givaroElement(self,i) )
        return tuple(array)

    def __dealloc__(FiniteField_givaro self):
        """
        Free the memory occupied by this Givaro finite field.
        """
        delete(self.objectptr)

    def characteristic(FiniteField_givaro self):
        """
        Return the characteristic of this field.

        EXAMPLES:
            sage: p = GF(19^5,'a').characteristic(); p
            19
            sage: type(p)
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(self.objectptr.characteristic())

    def order(FiniteField_givaro self):
        """
        Return the cardinality of this field.

        OUTPUT:
            Integer -- the number of elements in self.

        EXAMPLES:
            sage: n = GF(19^5,'a').order(); n
            2476099
            sage: type(n)
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(self.order_c())

    cdef order_c(FiniteField_givaro self):
        return self.objectptr.cardinality()

    def degree(FiniteField_givaro self):
        r"""
        If \code{self.cardinality() == p^n} this method returns $n$.

        OUTPUT:
            Integer -- the degree

        EXAMPLES:
            sage: GF(3^4,'a').degree()
            4
        """
        return Integer(self.objectptr.exponent())

    def is_atomic_repr(FiniteField_givaro self):
        """
        Return whether elements of self are printed using an atomic
        representation.

        EXAMPLES:
            sage: GF(3^4,'a').is_atomic_repr()
            False
        """
        if self.repr==0: #modulus
            return False
        else:
            return True

    def random_element(FiniteField_givaro self):
        """
        Return a random element of self.

        EXAMPLES:
            sage: k = GF(23**3, 'a')
            sage: e = k.random_element()
            sage: type(e)
            <type 'sage.rings.finite_field_givaro.FiniteField_givaroElement'>
        """
        cdef int seed = current_randstate().c_random()
        cdef int res
        cdef GivRandom generator = GivRandomSeeded(seed)
        res = self.objectptr.random(generator,res)
        return make_FiniteField_givaroElement(self,res)

    def __call__(FiniteField_givaro self, e):
        """
        Coerces several data types to self.

        INPUT:
            e -- data to coerce

        EXAMPLES:

            FiniteField_givaroElement are accepted where the parent
            is either self, equals self or is the prime subfield

            sage: k = GF(2**8, 'a')
            sage: k.gen() == k(k.gen())
            True


            Floats, ints, longs, Integer are interpreted modulo characteristic

            sage: k(2)
            0

            Floats coerce in:
            sage: k(float(2.0))
            0

            Rational are interpreted as
                             self(numerator)/self(denominator).
            Both may not be >= self.characteristic().

            sage: k = GF(3**8, 'a')
            sage: k(1/2) == k(1)/k(2)
            True

            Free modulo elements over self.prime_subfield() are interpreted 'little endian'

            sage: k = GF(2**8, 'a')
            sage: e = k.vector_space().gen(1); e
            (0, 1, 0, 0, 0, 0, 0, 0)
            sage: k(e)
            a

            Strings are evaluated as polynomial representation of elements in self

            sage: k('a^2+1')
            a^2 + 1

        Univariate polynomials coerce into finite fields by evaluating
        the polynomial at the field's generator:
            sage: from sage.rings.finite_field_givaro import FiniteField_givaro
            sage: R.<x> = QQ[]
            sage: k, a = FiniteField_givaro(5^2, 'a').objgen()
            sage: k(R(2/3))
            4
	    sage: k(x^2)
	    a + 3
            sage: R.<x> = GF(5)[]
            sage: k(x^3-2*x+1)
            2*a + 4

            sage: x = polygen(QQ)
            sage: k(x^25)
            a

            sage: Q, q = FiniteField_givaro(5^3,'q').objgen()
            sage: L = GF(5)
            sage: LL.<xx> = L[]
            sage: Q(xx^2 + 2*xx + 4)
            q^2 + 2*q + 4


        Multivariate polynomials only coerce if constant:
            sage: R = k['x,y,z']; R
            Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 5^2
            sage: k(R(2))
            2
            sage: R = QQ['x,y,z']
            sage: k(R(1/5))
            Traceback (most recent call last):
            ...
	    ZeroDivisionError: division by zero in finite field.


            PARI elements are interpreted as finite field elements; this PARI flexibility
            is (absurdly!) liberal:

            sage: k = GF(2**8, 'a')
            sage: k(pari('Mod(1,2)'))
            1
            sage: k(pari('Mod(2,3)'))
            0
            sage: k(pari('Mod(1,3)*a^20'))
            a^7 + a^5 + a^4 + a^2

            GAP elements need to be finite field elements:

            sage: from sage.rings.finite_field_givaro import FiniteField_givaro
            sage: x = gap('Z(13)')
            sage: F = FiniteField_givaro(13)
            sage: F(x)
            2
            sage: F(gap('0*Z(13)'))
            0
            sage: F = FiniteField_givaro(13^2)
            sage: x = gap('Z(13)')
            sage: F(x)
            2
            sage: x = gap('Z(13^2)^3')
            sage: F(x)
            12*a + 11
            sage: F.multiplicative_generator()^3
            12*a + 11

            sage: k.<a> = GF(29^3)
            sage: k(48771/1225)
            28

            sage: from sage.rings.finite_field_givaro import FiniteField_givaro
            sage: F9 = FiniteField_givaro(9)
            sage: F81 = FiniteField_givaro(81)
            sage: F81(F9.gen())
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce from a finite field other than the prime subfield


        """

        cdef int res
        cdef int g
        cdef int x
        cdef int e_int


        ########

        if PY_TYPE_CHECK(e, FiniteField_givaroElement):
            if e.parent() is self:
                return e
            if e.parent() == self:
                return make_FiniteField_givaroElement(self,(<FiniteField_givaroElement>e).element)
            if e.parent() is self.prime_subfield_C() or e.parent() == self.prime_subfield_C():
                res = self.int_to_log(int(e))
            else:
                raise TypeError, "unable to coerce from a finite field other than the prime subfield"


        elif PY_TYPE_CHECK(e, int) or \
             PY_TYPE_CHECK(e, Integer) or \
             PY_TYPE_CHECK(e, long) or is_IntegerMod(e):
            try:
                e_int = e
                if e != e_int:       # overflow in Pyrex is often not detected correctly... but this is bullet proof.
                                     # sometimes it is detected correctly, so we do have to use exceptions though.
                                     # todo -- be more eloquent here!!
                    raise OverflowError
                res = self.objectptr.initi(res,e_int)
            except OverflowError:
                res = self.objectptr.initi(res,int(e)%int(self.objectptr.characteristic()))

        elif PY_TYPE_CHECK(e, float):
            res = self.objectptr.initd(res,e)

        elif PY_TYPE_CHECK(e, str):
            return self(eval(e.replace("^","**"),self.gens_dict()))

        elif PY_TYPE_CHECK(e, FreeModuleElement):
            if self.vector_space() != e.parent():
                raise TypeError, "e.parent must match self.vector_space"
            ret = self._zero_element
            for i in range(len(e)):
                ret = ret + self(int(e[i]))*self.gen()**i
            return ret

        elif PY_TYPE_CHECK(e, MPolynomial):
            if e.is_constant():
                return self(e.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

        elif PY_TYPE_CHECK(e, Polynomial):
            if e.is_constant():
                return self(e.constant_coefficient())
            else:
                return e(self.gen())

        elif PY_TYPE_CHECK(e, Rational):
            num = e.numer()
            den = e.denom()
            return self(num)/self(den)

        elif PY_TYPE_CHECK(e, gen):
            pass # handle this in next if clause

        elif PY_TYPE_CHECK(e,FiniteField_ext_pariElement):
            # reduce FiniteFieldElements to pari
            e = e._pari_()

        elif sage.interfaces.gap.is_GapElement(e):
            from sage.interfaces.gap import gfq_gap_to_sage
            return gfq_gap_to_sage(e, self)

        else:
            raise TypeError, "unable to coerce"

        if PY_TYPE_CHECK(e, gen):
            e = e.lift().lift()
            try:
                res = self.int_to_log(e[0])
            except TypeError:
                res = self.int_to_log(e)

            g = self.objectptr.sage_generator()
            x = self.objectptr.one

            for i from 0 < i <= e.poldegree():
                x = self.objectptr.mul(x,x,g)
                res = self.objectptr.axpyin( res, self.int_to_log(e[i]) , x)

        return make_FiniteField_givaroElement(self,res)

    cdef _coerce_c_impl(self, x):
        """
        Coercion accepts elements of self.parent(), ints, and prime subfield elements.
        """
        cdef int r, cx

        if PY_TYPE_CHECK(x, int) \
               or PY_TYPE_CHECK(x, long) or PY_TYPE_CHECK(x, Integer):
            cx = x % self.characteristic()
            r = (<FiniteField_givaro>self).objectptr.initi(r, cx%self.objectptr.characteristic())
            return make_FiniteField_givaroElement(<FiniteField_givaro>self,r)

        if PY_TYPE_CHECK(x, FiniteFieldElement) or \
               PY_TYPE_CHECK(x, FiniteField_givaroElement) or is_IntegerMod(x):
            K = x.parent()
            if K is <object>self:
                return x
            if PY_TYPE_CHECK(K, IntegerModRing_generic) \
                   and K.characteristic() % self.characteristic() == 0:
                return self(int(x))
            if K.characteristic() == self.characteristic():
                if K.degree() == 1:
                    return self(int(x))
                elif self.degree() % K.degree() == 0:
                    # This is where we *would* do coercion from one nontrivial finite field to another...
                    raise TypeError, 'no canonical coercion defined'
        raise TypeError, 'no canonical coercion defined'

    def gen(FiniteField_givaro self, int n=0):
        r"""
        Return a generator of self. All elements x of self are
        expressed as $\log_{self.gen()}(p)$ internally. If self is
        a prime field this method returns a generator of the
        multiplicative group (i.e., a  primitive root).

        EXAMPLES:
            sage: k = GF(3^4, 'b'); k.gen()
            b
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: only one generator
	    sage: F=sage.rings.finite_field_givaro.FiniteField_givaro(31)
	    sage: F.gen()
	    3
	    sage: F.gen().multiplicative_order()
	    30
        """
        if n > 0:
            raise IndexError, "only one generator"
        cdef int r
        from sage.rings.arith import primitive_root

        if self.degree() == 1:
            return self(primitive_root(self.order_c()))
        else:
            return make_FiniteField_givaroElement(self,self.objectptr.sage_generator())

    cdef prime_subfield_C(FiniteField_givaro self):
        if self._prime_subfield is None:
            self._prime_subfield = finite_field.FiniteField(self.characteristic())
        return self._prime_subfield

    def prime_subfield(FiniteField_givaro self):
        r"""
        Return the prime subfield $\FF_p$ of self if self is $\FF_{p^n}$.

        EXAMPLES:
            sage: GF(3^4, 'b').prime_subfield()
            Finite Field of size 3

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: S.prime_subfield()
            Finite Field of size 5
            sage: type(S.prime_subfield())
            <class 'sage.rings.finite_field_prime_modn.FiniteField_prime_modn'>
        """
        return self.prime_subfield_C()


    def log_to_int(FiniteField_givaro self, int n):
        r"""
        Given an integer $n$ this method returns $i$ where $i$
        satisfies \code{self.gen()^n == i}, if the result is
        interpreted as an integer.

        INPUT:
            n -- log representation of a finite field element

        OUTPUT:
            integer representation of a finite field element.

        EXAMPLE:
            sage: k = GF(2**8, 'a')
            sage: k.log_to_int(4)
            16
            sage: k.log_to_int(20)
            180
        """
        cdef int ret

        if n<0:
            raise ArithmeticError, "Cannot serve negative exponent %d"%n
        elif n>=self.order_c():
            raise IndexError, "n=%d must be < self.order()"%n
        _sig_on
        ret = int(self.objectptr.convert(ret, n))
        _sig_off
        return ret

    def int_to_log(FiniteField_givaro self, int n):
        r"""
        Given an integer $n$ this method returns $i$ where $i$ satisfies
        \code{self.gen()^i==(n\%self.characteristic())}.

        INPUT:
            n -- integer representation of an finite field element

        OUTPUT:
            log representation of n

        EXAMPLE:
            sage: k = GF(7**3, 'a')
            sage: k.int_to_log(4)
            228
            sage: k.int_to_log(3)
            57
            sage: k.gen()^57
            3
        """
        cdef int r
        _sig_on
        ret =  int(self.objectptr.initi(r,n))
        _sig_off
        return ret

    def fetch_int(FiniteField_givaro self, int n):
        r"""
        Given an integer $n$ return a finite field element in self
        which equals $n$ under the condition that  self.gen() is set to
        self.characteristic().

        EXAMPLE:
            sage: k.<a> = GF(2^8)
            sage: k.fetch_int(8)
            a^3
            sage: e = k.fetch_int(151); e
            a^7 + a^4 + a^2 + a + 1
            sage: 2^7 + 2^4 + 2^2 + 2 + 1
            151
        """
        cdef GivaroGfq *k = self.objectptr
        cdef int ret = k.zero
        cdef int a = k.sage_generator()
        cdef int at = k.one
        cdef unsigned int ch = k.characteristic()
        cdef int _n, t, i

        if n<0 or n>k.cardinality():
            raise TypeError, "n must be between 0 and self.order()"

        _n = n

        for i from 0 <= i < k.exponent():
            t = k.initi(t, _n%ch)
            ret = k.axpy(ret, t, at, ret)
            at = k.mul(at,at,a)
            _n = _n/ch
        return make_FiniteField_givaroElement(self, ret)

    def polynomial(self, name=None):
        """
        Return the defining polynomial of this field as an element of
        self.polynomial_ring().

        This is the same as the characteristic polynomial of the
        generator of self.

        INPUT:
            name -- optional name of the generator

        EXAMPLES:
            sage: k = GF(3^4, 'a')
            sage: k.polynomial()
            a^4 + 2*a^3 + 2
        """
        if self._polynomial is None or name is not None:
            quo = int(-(self.gen()**(self.degree())))
            b   = int(self.characteristic())

            ret = []
            for i in range(self.degree()):
                ret.append(quo%b)
                quo = quo/b
            ret = ret + [1]
            R = self.polynomial_ring(name)
            if name is None:
                self._polynomial = R(ret)
                return self._polynomial
            else:
                return R(ret)
        else:
            return self._polynomial

    def _pari_modulus(self):
        """
        EXAMPLES:
            sage: GF(3^4,'a')._pari_modulus()
            Mod(1, 3)*a^4 + Mod(2, 3)*a^3 + Mod(2, 3)
        """
        f = pari(str(self.modulus()))
        return f.subst('x', 'a') * pari("Mod(1,%s)"%self.characteristic())

    def _finite_field_ext_pari_(self):  # todo -- cache
        """
        Return a FiniteField_ext_pari isomorphic to self with the same
        defining polynomial.

        EXAMPLES:
            sage: GF(3^4,'z')._finite_field_ext_pari_()
            Finite Field in z of size 3^4
        """
        f = self.polynomial()
        return finite_field_ext_pari.FiniteField_ext_pari(self.order_c(),
                                                          self.variable_name(), f)

    def _finite_field_ext_pari_modulus_as_str(self):  # todo -- cache
        """
        Return a string representing the modulus as a PARI string,
        this is mainly of internal interest and might go away
        eventually.

        EXAMPLE:
            sage: GF(3^4,'z')._finite_field_ext_pari_modulus_as_str()
            'Mod(1, 3)*z^4 + Mod(2, 3)*z^3 + Mod(2, 3)'
            sage: pari('Mod(1, 3)*z^4 + Mod(2, 3)*z^3 + Mod(2, 3)')
            Mod(1, 3)*z^4 + Mod(2, 3)*z^3 + Mod(2, 3)
        """
        return self._finite_field_ext_pari_().modulus()._pari_init_()

    def __iter__(FiniteField_givaro self):
        """
        Finite fields may be iterated over:

        EXAMPLE:
            sage: list(GF(2**2, 'a'))
            [0, a, a + 1, 1]
        """
        return FiniteField_givaro_iterator(self)

    def __richcmp__(left, right, int op):
        r"""
        Compare \code{self} with \code{right}.

        EXAMPLE:
            sage: k.<a> = GF(2^3)
            sage: j.<b> = GF(3^4)
            sage: k == j
            False

            sage: GF(2^3,'a') == copy(GF(2^3,'a'))
            True
        """
        return (<Parent>left)._richcmp(right, op)

    def __hash__(FiniteField_givaro self):
        """
        The hash of a Givaro finite field is a hash over it's
        characterstic polynomial and the string 'givaro'

        EXAMPLES:
            sage: hash(GF(3^4, 'a'))
            -417021630                # 32-bit
            1006006598732398914       # 64-bit
        """
        if self._hash is None:
            pass
        if self.degree()>1:
            self._hash = hash((self.characteristic(),self.polynomial(),self.variable_name(),"givaro"))
        else:
            self._hash = hash((self.characteristic(),self.variable_name(),"givaro"))
        return self._hash

    def _element_repr(FiniteField_givaro self, FiniteField_givaroElement e):
        """
        Wrapper for log, int, and poly representations.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._element_repr(a^20)
            '2*a^3 + 2*a^2 + 2'

            sage: k = sage.rings.finite_field_givaro.FiniteField_givaro(3^4,'a', repr='int')
            sage: a = k.gen()
            sage: k._element_repr(a^20)
            '74'

            sage: k = sage.rings.finite_field_givaro.FiniteField_givaro(3^4,'a', repr='log')
            sage: a = k.gen()
            sage: k._element_repr(a^20)
            '20'
        """
        if self.repr==0:
            return self._element_poly_repr(e)
        elif self.repr==1:
            return self._element_log_repr(e)
        else:
            return self._element_int_repr(e)

    def _element_log_repr(FiniteField_givaro self, FiniteField_givaroElement e):
        r"""
        Return str(i) where \code{base.gen()^i==self}.

        TODO -- this isn't what it does -- see below.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._element_log_repr(a^20)
            '20'
            sage: k._element_log_repr(a)    # TODO -- why 41 ?? why not 1?
            '41'
        """
        return str(int(e.element))

    def _element_int_repr(FiniteField_givaro self, FiniteField_givaroElement e):
        """
        Return integer representation of e.

	Elements of this field are represented as ints in as follows:
        for $e \in \FF_p[x]$ with $e = a_0 + a_1x + a_2x^2 + \cdots $, $e$ is
        represented as: $n= a_0 + a_1  p + a_2  p^2 + \cdots$.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._element_int_repr(a^20)
            '74'
        """
        return str(int(e))

    def _element_poly_repr(FiniteField_givaro self, FiniteField_givaroElement e, varname = None):
        """
        Return a polynomial expression in base.gen() of self.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._element_poly_repr(a^20)
            '2*a^3 + 2*a^2 + 2'
        """
        if varname is None:
            variable = self.variable_name()
        else:
            variable = varname

        quo = self.log_to_int(e.element)
        b   = int(self.characteristic())

        ret = ""
        for i in range(self.degree()):
            coeff = quo%b
            if coeff != 0:
                if i>0:
                    if coeff==1:
                        coeff=""
                    else:
                        coeff=str(coeff)+"*"
                    if i>1:
                        ret = coeff + variable + "^" + str(i) + " + " + ret
                    else:
                        ret = coeff + variable + " + " + ret
                else:
                    ret = str(coeff) + " + " + ret
            quo = quo/b
        if ret == '':
            return "0"
        return ret[:-3]

    def a_times_b_plus_c(FiniteField_givaro self,FiniteField_givaroElement a, FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return r = a*b + c. This is faster than multiplying a and b
        first and adding c to the result.

        INPUT:
            a -- FiniteField_givaroElement
            b -- FiniteField_givaroElement
            c -- FiniteField_givaroElement

        EXAMPLE:
            sage: k.<a> = GF(2**8)
            sage: k.a_times_b_plus_c(a,a,k(1))
            a^2 + 1
        """
        cdef int r

        r = self.objectptr.axpy(r, a.element, b.element, c.element)
        return make_FiniteField_givaroElement(self,r)

    def a_times_b_minus_c(FiniteField_givaro self,FiniteField_givaroElement a, FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return r = a*b - c.

        INPUT:
            a -- FiniteField_givaroElement
            b -- FiniteField_givaroElement
            c -- FiniteField_givaroElement

        EXAMPLE:
            sage: k.<a> = GF(3**3)
            sage: k.a_times_b_minus_c(a,a,k(1))
            a^2 + 2
        """
        cdef int r

        r = self.objectptr.axmy(r, a.element, b.element, c.element, )
        return make_FiniteField_givaroElement(self,r)

    def c_minus_a_times_b(FiniteField_givaro self,FiniteField_givaroElement a,
                          FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return r = c - a*b.

        INPUT:
            a -- FiniteField_givaroElement
            b -- FiniteField_givaroElement
            c -- FiniteField_givaroElement

        EXAMPLE:
            sage: k.<a> = GF(3**3)
            sage: k.c_minus_a_times_b(a,a,k(1))
            2*a^2 + 1
        """
        cdef int r

        r = self.objectptr.amxy(r , a.element, b.element, c.element, )
        return make_FiniteField_givaroElement(self,r)

    def __reduce__(FiniteField_givaro self):
        """
        Pickle self:

        EXAMPLE:
            sage: k.<a> = GF(2**8)
            sage: loads(dumps(k)) == k
            True
        """
        if self._array is None:
            cache = 0
        else:
            cache = 1

        return sage.rings.finite_field_givaro.unpickle_FiniteField_givaro, \
               (self.order_c(),(self.variable_name(),),
                map(int,list(self.modulus())),int(self.repr),int(self._array is not None))

def unpickle_FiniteField_givaro(order,variable_name,modulus,rep,cache):
    """
    EXAMPLE:
       sage: k = GF(3**7, 'a')
       sage: loads(dumps(k)) == k # indirect doctest
       True
    """
    from sage.rings.arith import is_prime

    if rep == 0:
        rep = 'poly'
    elif rep == 1:
        rep = 'log'
    elif rep == 2:
        rep = 'int'

    if not is_prime(order):
        return FiniteField_givaro(order,variable_name,modulus,rep,cache=cache)
    else:
        return FiniteField_givaro(order,cache=cache)

cdef class FiniteField_givaro_iterator:
    """
    Iterator over FiniteField_givaro elements of degree 1. We iterate
    over fields of higher degree using the VectorSpace iterator.

    EXAMPLES:
        sage: for x in GF(2^2,'a'): print x
        0
        a
        a + 1
        1
    """

    def __init__(self, FiniteField_givaro parent):
        """
            sage: k.<a> = GF(3^4)
            sage: i = iter(k) # indirect doctest
            sage: i
            Iterator over Finite Field in a of size 3^4
        """
        self._parent = parent
        self.iterator = -1

    def __next__(self):
        """
            sage: k.<a> = GF(3^4)
            sage: i = iter(k) # indirect doctest
            sage: i.next()
            0
            sage: i.next()
            2*a
        """

        self.iterator=self.iterator+1

        if self.iterator==self._parent.order_c():
            self.iterator = -1
            raise StopIteration

        return make_FiniteField_givaroElement(self._parent,self.iterator)

    def __repr__(self):
        """
        EXAMPLE:
            sage: k.<a> = GF(3^4)
            sage: i = iter(k)
            sage: i # indirect doctest
            Iterator over Finite Field in a of size 3^4
        """
        return "Iterator over %s"%self._parent

cdef FiniteField_givaro_copy(FiniteField_givaro orig):
    cdef FiniteField_givaro copy
    copy = FiniteField_givaro(orig.characteristic()**orig.degree())
    delete(copy.objectptr)
    copy.objectptr = gfq_factorycopy(gfq_deref(orig.objectptr))
    copy._zero_element = make_FiniteField_givaroElement(copy,copy.objectptr.zero)
    copy._one_element = make_FiniteField_givaroElement(copy,copy.objectptr.one)

    return copy

cdef class FiniteField_givaroElement(FiniteFieldElement):
    """
    An element of a (Givaro) finite field.
    """

    def __init__(FiniteField_givaroElement self, parent ):
        """
        Initializes an element in parent. It's much better to use
        parent(<value>) or any specialized method of parent
        like gen() instead. In general do not call this
        constructor directly.

        Alternatively you may provide a value which is directly
        assigned to this element. So the value must represent the
        log_g of the value you wish to assign.

        INPUT:
            parent -- base field

        OUTPUT:
            finite field element.

        EXAMPLE:
            sage: k.<a> = GF(5^2)
            sage: from sage.rings.finite_field_givaro import FiniteField_givaroElement
            sage: FiniteField_givaroElement(k)
            0

        """
        self._parent = <ParentWithBase> parent  # explicit case required for C++
        self.element = 0

    cdef FiniteField_givaroElement _new_c(self, int value):
        return make_FiniteField_givaroElement(parent_object(self), value)

    def __dealloc__(FiniteField_givaroElement self):
        pass

    def _repr_(FiniteField_givaroElement self):
        """
        EXAMPLE:
            sage: k.<FOOBAR> = GF(2^4)
            sage: FOOBAR #indirect doctest
            FOOBAR

            sage: k.<FOOBAR> = GF(2^4,repr='log')
            sage: FOOBAR
            1

            sage: k.<FOOBAR> = GF(2^4,repr='int')
            sage: FOOBAR
            2
        """
        return (<FiniteField_givaro>self._parent)._element_repr(self)

    def parent(self):
        """
        Return parent finite field.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: (a*a).parent()
            Finite Field in a of size 3^4
        """
        return (<FiniteField_givaro>self._parent)

    def __nonzero__(FiniteField_givaroElement self):
        r"""
        Return True if \code{self != k(0)}.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: a.is_zero()
            False
            sage: k(0).is_zero()
            True
        """
        return not (<FiniteField_givaro>self._parent).objectptr.isZero(self.element)

    def is_one(FiniteField_givaroElement self):
        r"""
        Return True if \code{self == k(1)}.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: a.is_one()
            False
            sage: k(1).is_one()
            True
        """
        return (<FiniteField_givaro>self._parent).objectptr.isOne(self.element)

    def is_unit(FiniteField_givaroElement self):
        """
        Return True if self is nonzero, so it is a unit as an element of the
        finite field.

        EXAMPLES:
            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: a.is_unit()
            True
            sage: k(0).is_unit()
            False
        """
        return not (<FiniteField_givaro>self._parent).objectptr.isZero(self.element)
        # **WARNING** Givaro seems to define unit to mean in the prime field,
        # which is totally wrong!  It's a confusion with the underlying polynomial
        # representation maybe??  That's why the following is commented out.
        # return (<FiniteField_givaro>self._parent).objectptr.isunit(self.element)


    def is_square(FiniteField_givaroElement self):
        """
        Return True if self is a square in self.parent()

        EXAMPLES:
            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: a.is_square()
            False
            sage: v = set([x^2 for x in k])
            sage: [x.is_square() for x in v]
            [True, True, True, True, True]
            sage: [x.is_square() for x in k if not x in v]
            [False, False, False, False]

        ALGORITHM:
            Elements are stored as powers of generators, so we simply check
            to see if it is an even power of a generator.

        TESTS:
            sage: K = GF(27, 'a')
            sage: set([a*a for a in K]) == set([a for a in K if a.is_square()])
            True
            sage: K = GF(25, 'a')
            sage: set([a*a for a in K]) == set([a for a in K if a.is_square()])
            True
            sage: K = GF(16, 'a')
            sage: set([a*a for a in K]) == set([a for a in K if a.is_square()])
            True
        """
        cdef FiniteField_givaro field = <FiniteField_givaro>self._parent
        if field.objectptr.characteristic() == 2:
            return True
        elif self.element == field.objectptr.one:
            return True
        else:
            return self.element % 2 == 0

    def sqrt(FiniteField_givaroElement self, all=False, extend=False):
        """
        Return a square root of this finite field element in its
        parent, if there is one.  Otherwise, raise a ValueError.

        EXAMPLES:
            sage: k.<a> = GF(7^2)
            sage: k(2).sqrt()
            3
            sage: k(3).sqrt()
            2*a + 6
            sage: k(3).sqrt()**2
            3
            sage: k(4).sqrt()
            2
            sage: k.<a> = GF(7^3)
            sage: k(3).sqrt()
            Traceback (most recent call last):
            ...
            ValueError: must be a perfect square.

        ALGORITHM:
            Self is stored as $a^k$ for some generator $a$.
            Return $a^(k/2)$ for even $k$.

        TESTS:
            sage: K = GF(49, 'a')
            sage: all([a.sqrt()*a.sqrt() == a for a in K if a.is_square()])
            True
            sage: K = GF(27, 'a')
            sage: all([a.sqrt()*a.sqrt() == a for a in K if a.is_square()])
            True
            sage: K = GF(8, 'a')
            sage: all([a.sqrt()*a.sqrt() == a for a in K if a.is_square()])
            True
        """
        if all:
            a = self.sqrt()
            return [a, -a] if -a != a else [a]
        cdef FiniteField_givaro field = <FiniteField_givaro>self._parent
        if self.element == field.objectptr.one:
            return make_FiniteField_givaroElement(field, field.objectptr.one)
        elif self.element % 2 == 0:
            return make_FiniteField_givaroElement(field, self.element/2)
        elif field.objectptr.characteristic() == 2:
            return make_FiniteField_givaroElement(field, (field.objectptr.cardinality() - 1 + self.element)/2)
        elif extend:
            raise NotImplementedError # TODO: fix this once we have nested embeddings of finite fields
        else:
            raise ValueError, "must be a perfect square."

    def nth_root(self, int n, extend = False, all = False):
        r"""
        Returns an nth root of self.

        INPUT:
            n -- integer >= 1 (must fit in C int type)
            extend -- bool (default: True); if True, return a square
                 root in an extension ring, if necessary. Otherwise,
                 raise a ValueError if the square is not in the base
                 ring.
            all -- bool (default: False); if True, return all square
                 roots of self, instead of just one.

        OUTPUT:
           If self has an nth root, returns one (if all == False) or a list of
           all of them (if all == True).  Otherwise, raises a ValueError (if
           extend = False) or a NotImplementedError (if extend = True).

        AUTHOR:
           -- David Roe (2007-10-3)

        EXAMPLES:
        sage: k.<a> = GF(29^2)
        sage: b = a^2 + 5*a + 1
        sage: b.nth_root(11)
        3*a + 20
        sage: b.nth_root(5)
        Traceback (most recent call last):
        ...
        ValueError: no nth root
        sage: b.nth_root(5, all = True)
        []
        sage: b.nth_root(3, all = True)
        [14*a + 18, 10*a + 13, 5*a + 27]
        """
        if extend:
            raise NotImplementedError
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
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

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two elements.

        EXAMPLE:
            sage: k.<b> = GF(9**2)
            sage: b^10 + 2*b
            2*b^3 + 2*b^2 + 2*b + 1
        """
        cdef int r
        r = parent_object(self).objectptr.add(r, self.element ,
                                              (<FiniteField_givaroElement>right).element )
        return make_FiniteField_givaroElement(parent_object(self),r)

    cdef ModuleElement _iadd_c_impl(self, ModuleElement right):
        """
        Add two elements inplace.

        EXAMPLE:
            sage: k.<b> = GF(9**2)
            sage: b^10 + 2*b
            2*b^3 + 2*b^2 + 2*b + 1
        """
        cdef int r
        self.element = parent_object(self).objectptr.add(r, self.element ,
                                                         (<FiniteField_givaroElement>right).element )
        return self


    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Multiply two elements:

        EXAMPLE:
            sage: k.<c> = GF(7**4)
            sage: 3*c
            3*c
            sage: c*c
            c^2
        """
        cdef int r
        r = parent_object(self).objectptr.mul(r, self.element,
                                              (<FiniteField_givaroElement>right).element)
        return make_FiniteField_givaroElement(parent_object(self),r)


    cdef RingElement _imul_c_impl(self, RingElement right):
        """
        Multiply two elements inplace.

        EXAMPLE:
            sage: k.<c> = GF(7**4)
            sage: 3*c
            3*c
            sage: c*c
            c^2
        """
        cdef int r
        self.element = parent_object(self).objectptr.mul(r, self.element,
                                                         (<FiniteField_givaroElement>right).element)
        return self

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Divide two elements

        EXAMPLE:
            sage: k.<g> = GF(2**8)
            sage: g/g
            1

            sage: k(1) / k(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field.
        """
        cdef int r
        if (<FiniteField_givaroElement>right).element == 0:
            raise ZeroDivisionError, 'division by zero in finite field.'
        r = parent_object(self).objectptr.div(r, self.element,
                                              (<FiniteField_givaroElement>right).element)
        return make_FiniteField_givaroElement(parent_object(self),r)

    cdef RingElement _idiv_c_impl(self, RingElement right):
        """
        Divide two elements inplace

        EXAMPLE:
            sage: k.<g> = GF(2**8)
            sage: g/g
            1

            sage: k(1) / k(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field.
        """

        cdef int r
        if (<FiniteField_givaroElement>right).element == 0:
            raise ZeroDivisionError, 'division by zero in finite field.'
        self.element = parent_object(self).objectptr.div(r, self.element,
                                                         (<FiniteField_givaroElement>right).element)
        return self

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two elements

        EXAMPLE:
            sage: k.<a> = GF(3**4)
            sage: k(3) - k(1)
            2
            sage: 2*a - a^2
            2*a^2 + 2*a
        """
        cdef int r
        r = parent_object(self).objectptr.sub(r, self.element,
                                              (<FiniteField_givaroElement>right).element)
        return make_FiniteField_givaroElement(parent_object(self),r)

    cdef ModuleElement _isub_c_impl(self, ModuleElement right):
        """
        Subtract two elements inplace

        EXAMPLE:
            sage: k.<a> = GF(3**4)
            sage: k(3) - k(1)
            2
            sage: 2*a - a^2
            2*a^2 + 2*a
        """
        cdef int r
        self.element = parent_object(self).objectptr.sub(r, self.element,
                                                         (<FiniteField_givaroElement>right).element)
        return self

    def __neg__(FiniteField_givaroElement self):
        """
        Negative of an element.

        EXAMPLES:
            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: -a
            2*a
        """
        cdef int r

        r = (<FiniteField_givaro>self._parent).objectptr.neg(r, self.element)
        return make_FiniteField_givaroElement((<FiniteField_givaro>self._parent),r)

    def __invert__(FiniteField_givaroElement self):
        """
        Return the multiplicative inverse of an element.

        EXAMPLES:
            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: ~a
            a + 2
            sage: ~a*a
            1
        """
        cdef int r

        (<FiniteField_givaro>self._parent).objectptr.inv(r, self.element)
        return make_FiniteField_givaroElement((<FiniteField_givaro>self._parent),r)

    def __pow__(FiniteField_givaroElement self, exp, other):
        """
        EXAMPLE:
            sage: K.<a> = GF(3^3, 'a')
            sage: a^3 == a*a*a
            True
            sage: b = a+1
            sage: b^5 == b^2 * b^3
            True
            sage: b^(-1) == 1/b
            True
            sage: b = K(-1)
            sage: b^2 == 1
            True

        ALGORITHM:
            Givaro objects are stored as integers $i$ such that $self=a^i$, where
            $a$ is a generator of $K$ (though not necessarily the one returned by K.gens()).
            Now it is trivial to compute $(a^i)^exp = a^(i*exp)$, and reducing the exponent
            mod the multiplicative order of $K$.

        AUTHOR:
            Robert Bradshaw
        """
        _exp = int(exp)
        if _exp != exp:
            raise ValueError, "exponent must be an integer"
        exp = _exp

        cdef int r
        cdef int order
        cdef FiniteField_givaro field
        field = <FiniteField_givaro>self._parent

        if (field.objectptr).isOne(self.element):
            return self

        elif exp == 0:
            if (field.objectptr).isZero(self.element):
                raise ArithmeticError, "0^0 is undefined."
            return make_FiniteField_givaroElement(field, field.objectptr.one)

        elif (field.objectptr).isZero(self.element):
            if exp < 0:
                raise ZeroDivisionError
            return make_FiniteField_givaroElement(field, field.objectptr.zero)

        order = (field.order_c()-1)
        exp = exp % order

        r = exp
        if r == 0:
            return make_FiniteField_givaroElement(field, field.objectptr.one)

        r = (r * self.element) % order
        if r < 0:
            r = r + order
        if r == 0:
            return make_FiniteField_givaroElement(field, field.objectptr.one)

        return make_FiniteField_givaroElement(field, r)

    def __richcmp__(left, right, int op):
        """
        EXAMPLE:
            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: a == k('a') # indirect doctest
            True
            sage: a == a + 1
            False
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Comparison of finite field elements is correct or equality
        tests and somewhat random for < and > type of
        comparisons. This implementation performs these tests by
        comparing the underlying int representations.

        EXAMPLES:
            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: a == k('a')
            True
            sage: a == a + 1
            False

            Even though inequality tests do return answers, they
            really make no sense as finite fields are unordered. Thus,
            you cannot rely on the result as it is implementation specific.

            sage: a < a^2
            True
            sage: a^2 < a
            False
        """
        cdef FiniteField_givaro F
        F = left._parent

        return cmp(F.log_to_int(left.element), F.log_to_int((<FiniteField_givaroElement>right).element) )

    def __int__(FiniteField_givaroElement self):
        """
        Return the int representation of self.  When self is in the
        prime subfield, the integer returned is equal to self and not
        to \code{log_repr}.

        Elements of this field are represented as ints in as follows:
        for $e \in \FF_p[x]$ with $e = a_0 + a_1x + a_2x^2 + \cdots $, $e$ is
        represented as: $n= a_0 + a_1  p + a_2  p^2 + \cdots$.

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: int(k(4))
            4
            sage: int(b)
            5
            sage: type(int(b))
            <type 'int'>
        """
        return (<FiniteField_givaro>self._parent).log_to_int(self.element)


    def log_to_int(FiniteField_givaroElement self):
        r"""
        Returns the int representation of self, as a SAGE integer.   Use
        int(self) to directly get a Python int.

        Elements of this field are represented as ints in as follows:
        for $e \in \FF_p[x]$ with $e = a_0 + a_1x + a_2x^2 + \cdots $, $e$ is
        represented as: $n= a_0 + a_1  p + a_2  p^2 + \cdots$.

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: k(4).log_to_int()
            4
            sage: b.log_to_int()
            5
            sage: type(b.log_to_int())
            <type 'sage.rings.integer.Integer'>
        """
        return Integer((<FiniteField_givaro>self._parent).log_to_int(self.element))

    def log(FiniteField_givaroElement self, base):
        """
        Return the log to the base b of self, i.e., an integer n
        such that b^n = self.

        WARNING:  TODO -- This is currently implemented by solving the discrete
        log problem -- which shouldn't be needed because of how finit field
        elements are represented.

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: a = b^7
            sage: a.log(b)
            7
        """
        q = (<FiniteField_givaro> self.parent()).order()
        b = self.parent()(base)
        return sage.groups.generic.discrete_log(self, b, q-1)

    def int_repr(FiniteField_givaroElement self):
        r"""
        Return the sring representation of self as an int (as returned
        by \code{log_to_int}).

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+1).int_repr()
            '6'
        """
        return (<FiniteField_givaro>self._parent)._element_int_repr(self)

    def log_repr(FiniteField_givaroElement self):
        r"""
        Return the log representation of self as a string.  See the
        documentation of the \code{_element_log_repr} function of the
        parent field.

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+2).log_repr()
            '3'
        """
        return (<FiniteField_givaro>self._parent)._element_log_repr(self)

    def poly_repr(FiniteField_givaroElement self):
        r"""
        Return representation of this finite field element as a polynomial
        in the generator.

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+2).poly_repr()
            'b + 2'
        """
        return (<FiniteField_givaro>self._parent)._element_poly_repr(self)

    def polynomial(FiniteField_givaroElement self, name=None):
        """
        Return self viewed as a polynomial over self.parent().prime_subfield().

        EXAMPLES:
            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: f = (b^2+1).polynomial(); f
            b + 4
            sage: type(f)
            <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_mod_p'>
            sage: parent(f)
            Univariate Polynomial Ring in b over Finite Field of size 5
        """
        cdef FiniteField_givaro F
        F = self._parent
        quo = F.log_to_int(self.element)
        b   = int(F.characteristic())
        ret = []
        for i in range(F.degree()):
            coeff = quo%b
            ret.append(coeff)
            quo = quo/b
        if not name is None and F.variable_name() != name:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            return PolynomialRing(F.prime_subfield_C(), name)(ret)
        else:
            return F.polynomial_ring()(ret)

    def _finite_field_ext_pari_element(FiniteField_givaroElement self, k=None):
        """
        Return an element of k supposed to match this element.  No
        checks if k equals self.parent() are performed.

        INPUT:
            k -- (optional) FiniteField_ext_pari

        OUTPUT:
            k.gen()^(self.log_repr())

        EXAMPLES:
            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b.charpoly('x')
            x^2 + 4*x + 2
            sage: P = S._finite_field_ext_pari_(); type(P)
            <class 'sage.rings.finite_field_ext_pari.FiniteField_ext_pari'>
            sage: c = b._finite_field_ext_pari_element(P); c
            b
            sage: type(c)
            <class 'sage.rings.finite_field_element.FiniteField_ext_pariElement'>
            sage: c.charpoly('x')
            x^2 + 4*x + 2

        The PARI field is automatically determined if it is not given:
            sage: d = b._finite_field_ext_pari_element(); d
            b
            sage: type(d)
            <class 'sage.rings.finite_field_element.FiniteField_ext_pariElement'>
        """
        if k is None:
            k = (<FiniteField_givaro>self._parent)._finite_field_ext_pari_()
        elif not isinstance(k, finite_field_ext_pari.FiniteField_ext_pari):
            raise TypeError, "k must be a pari finite field."

        variable = k.gen()._pari_()

        quo = int(self)
        b   = int((<FiniteField_givaro>self._parent).characteristic())

        ret = k._pari_one() - k._pari_one()    # TODO -- weird
        i = 0
        while quo!=0:
            coeff = quo%b
            if coeff != 0:
                ret = coeff * variable ** i + ret
            quo = quo/b
            i = i+1
        return k(ret)

    def _pari_(FiniteField_givaroElement self, var=None):
        """
        Return PARI representation of this finite field element

        INPUT:
            var -- optional variable string (default: none)

        EXAMPLES:
            sage: k.<a> = GF(5^3)
            sage: a._pari_()
            Mod(a, Mod(1, 5)*a^3 + Mod(3, 5)*a + Mod(3, 5))

            sage: a._pari_('b')
            Mod(b, Mod(1, 5)*b^3 + Mod(3, 5)*b + Mod(3, 5))
        """
        return pari(self._pari_init_(var))

    def _magma_init_(self):
        """
        Return a string representation of self that MAGMA can
        understand.

        EXAMPLE:
            sage: k.<a> = GF(3^5)
            sage: k._magma_init_() # requires optional MAGMA installation
            'ext< GF(3) | Polynomial(GF(3), [GF(3)!1,GF(3)!2,GF(3)!0,GF(3)!0,GF(3)!0,GF(3)!1]) >'
        """
        km = self.parent()._magma_()
        vn = km.gen(1).name()
        return self.parent()._element_poly_repr(self,vn)

    def multiplicative_order(FiniteField_givaroElement self):
        """
        Return the multiplicative order of this field element.

        EXAMPLES:
            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b.multiplicative_order()
            24
            sage: (b^6).multiplicative_order()
            4
        """
        # TODO -- I'm sure this can be made vastly faster
        # using how elements are represented as a power of the generator ??

        # code copy'n'pasted from finite_field_element.py
        import sage.rings.arith

        if self.__multiplicative_order is not None:
            return self.__multiplicative_order
        else:
            if self.is_zero():
                return ArithmeticError, "Multiplicative order of 0 not defined."
            n = (parent_object(self)).order_c() - 1
            order = 1
            for p, e in sage.rings.arith.factor(n):
                # Determine the power of p that divides the order.
                a = self**(n/(p**e))
                while a != 1:
                    order = order * p
                    a = a**p

            self.__multiplicative_order = order
            return order

    def __copy__(self):
        """
        Return a copy of this element.  Actually just returns self, since
        finite field elements are immutable.

        EXAMPLES:
            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: c = copy(b); c
            b
            sage: c is b
            True
            sage: copy(5r) is 5r
            True
        """
        return self

    def _gap_init_(FiniteField_givaroElement self):
        """
        Return a string that evaluates to the GAP representation of
        this element.

        A NotImplementedError is raised if self.parent().modulus() is
        not a Conway polynomial, as the isomorphism of finite fields is
        not implemented yet.

        EXAMPLES:
            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: (4*b+3)._gap_init_()
            'Z(25)^3'
            sage: S(gap('Z(25)^3'))
            4*b + 3
        """
        #copied from finite_field_element.py
        cdef FiniteField_givaro F
        F = parent_object(self)
        if not F._is_conway:
            raise NotImplementedError, "conversion of (Givaro) finite field element to GAP not implemented except for fields defined by Conway polynomials."
        if F.order_c() > 65536:
            raise TypeError, "order (=%s) must be at most 65536."%F.order_c()
        if self == 0:
            return '0*Z(%s)'%F.order_c()
        assert F.degree() > 1
        g = F.multiplicative_generator()
        n = self.log(g)
        return 'Z(%s)^%s'%(F.order_c(), n)

    def __hash__(FiniteField_givaroElement self):
        """
        Return the hash of this finite field element.  We hash the parent
        and the underlying integer representation of this element.

        EXAMPLES:
            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: hash(a)
            5
        """
        return hash(self.log_to_int())

    def vector(FiniteField_givaroElement self, reverse=False):
        """
        Return a vector in self.parent().vector_space() matching
        self. The most significant bit is to the right.

        INPUT:
            reverse -- reverse the order of the bits
                       from little endian to big endian.

        EXAMPLES:
            sage: k.<a> = GF(2^4)
            sage: e = a^2 + 1
            sage: v = e.vector()
            sage: v
            (1, 0, 1, 0)
            sage: k(v)
            a^2 + 1

            sage: k.<a> = GF(3^4)
            sage: e = 2*a^2 + 1
            sage: v = e.vector()
            sage: v
            (1, 0, 2, 0)
            sage: k(v)
            2*a^2 + 1

            You can also compute the vector in the other order:

            sage: e.vector(reverse=True)
            (0, 2, 0, 1)
        """
        cdef FiniteField_givaro k = <FiniteField_givaro>self._parent

        quo = k.log_to_int(self.element)
        b   = int(k.characteristic())

        ret = []
        for i in range(k.degree()):
            coeff = quo%b
            ret.append(coeff)
            quo = quo/b
        if reverse:
            ret = list(reversed(ret))
        return k.vector_space()(ret)

    def __reduce__(FiniteField_givaroElement self):
        """
        Used for supporting pickling of finite field elements.

        EXAMPLE:
            sage: k = GF(2**8, 'a')
            sage: e = k.random_element()
            sage: loads(dumps(e)) == e
            True
        """
        return unpickle_FiniteField_givaroElement,(parent_object(self),self.element)



def unpickle_FiniteField_givaroElement(FiniteField_givaro parent, int x):
    """
    EXAMPLE:
        sage: k = GF(3**4, 'a')
        sage: e = k.random_element()
        sage: loads(dumps(e)) == e # indirect doctest
        True
    """
    return make_FiniteField_givaroElement(parent, x)

cdef inline FiniteField_givaroElement make_FiniteField_givaroElement(FiniteField_givaro parent, int x):
    cdef FiniteField_givaroElement y

    if parent._array is None:
        #y = FiniteField_givaroElement(parent)
        y = PY_NEW(FiniteField_givaroElement)
        y._parent = <ParentWithBase> parent
        y.element = x
        return y
    else:
        return <FiniteField_givaroElement>parent._array[x]

