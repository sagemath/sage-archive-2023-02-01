r"""
 Finite (Extension) Fields of cardinality < 2**16.

 The arithmetic is performed by the Givaro C++ library which uses Zech
 logs internally to represent finite field elements. This
 implementation is the default finite extension field implementation
 in SAGE for the cardinality < 2**16 as it is way faster than the PARI
 implementation which uses polynomials to represent finite field
 elements. Some functionality in this class however is implemented
 using the PARI implementation.

 AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-06-05)


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


from sage.rings.ring cimport FiniteField
from sage.rings.coerce import bin_op
from sage.structure.element cimport FiniteFieldElement, Element
from sage.rings.finite_field_element import FiniteField_ext_pariElement
from sage.structure.sage_object cimport SageObject
import operator
import sage.rings.arith

import sage.interfaces.gap
from sage.libs.pari.all import pari
from sage.libs.pari.gen import gen

cdef extern from "stdlib.h":
    ctypedef int size_t
    void free(void *ptr)
    void *malloc(size_t size)
    void *realloc(void *ptr, size_t size)
    size_t strlen(char *s)
    char *strcpy(char *dest, char *src)

## cdef extern from 'interrupt.h':
##     int _sig_on, _sig_off, _sig_check
##     void _sig_str(char*)
cdef int _sig_on
cdef int _sig_off
cdef int _sig_check

cdef class FiniteField_givaro(FiniteField) #forward declaration



cdef extern from "Python.h":
    ctypedef struct PyTypeObject
    ctypedef struct PyObject
    int PyObject_TypeCheck(object o, PyTypeObject *t)

cdef extern from "givaro/givrandom.h":
    ctypedef struct GivRandom "GivRandom":
        pass

cdef extern from "givaro/givgfq.h":
    ctypedef struct intvec "std::vector<unsigned int>":
        void (* push_back)(int elem)

    ctypedef struct constintvec "const std::vector<unsigned int>"

    intvec intvec_factory "std::vector<unsigned int>"(int len)

cdef extern from "givaro/givgfq.h":

    ctypedef struct GivaroGfq "GFqDom<int>":
        #attributes
        unsigned int one
        unsigned int zero

        # methods
        int (* mul)(int r, int a, int b)
        int (* add)(int r, int a, int b)
        int (* sub)(int r, int a, int b)
        int (* div)(int r, int a, int b)
        int (* inv)(int r, int x)
        int (* neg)(int r, int x)
        int (* mulin)(int a, int b)
        unsigned int (* characteristic)()
        unsigned int (* cardinality)()
        int (* exponent)()
        int (* random)(GivRandom gen, int res)
        int (* initi "init")(int res, int e)
        int (* initd "init")(int res, double e)
        int (* axpyin)(int r, int a, int x)
        int (* sage_generator)() # SAGE specific method, not found upstream
        int (* write)(int r, int p)
        int (* read)(int r, int p)
        int (* axpy)(int r, int a, int b, int c)
        int (* axmy)(int r, int a, int b, int c)
        int (* amxy)(int r, int a, int b, int c)
        int (* isZero)(int e)
        int (* isOne)(int e)
        int (* isunit)(int e)

    GivaroGfq *gfq_factorypk "new GFqDom<int>" (unsigned int p, unsigned int k)
    # SAGE specific method, not found upstream
    GivaroGfq *gfq_factorypkp "new GFqDom<int>" (unsigned int p, unsigned int k, intvec poly)
    GivaroGfq *gfq_factorycopy "new GFqDom<int>"(GivaroGfq orig)
    GivaroGfq  gfq_deref "*"(GivaroGfq *orig)
    void delete "delete "(void *o)
    int gfq_element_factory "GFqDom<int>::Element"()

cdef class FiniteField_givaroElement(FiniteFieldElement) # forward declaration

cdef FiniteField_givaro parent_object(Element o):
    return <FiniteField_givaro>(o._parent)

cdef PyTypeObject *type_object(object o):
    return <PyTypeObject*><PyObject*>o

cdef class FiniteField_givaro(FiniteField):
    """
    Fnite Field. These are implemented using Zech logs and the
    cardinality must be < 2^16. See FiniteField_ext_pari for larger
    cardinalities.
    """

    cdef GivaroGfq *objectptr # C++ object
    cdef object _polynomial_ring
    cdef object _prime_subfield
    cdef object _array
    cdef int repr

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
            sage: -a^k.degree()
            a^4 + a^3 + a^2 + 1
            sage: f = k.modulus(); f
            x^8 + x^4 + x^3 + x^2 + 1


            You may enforce a modulus:

            sage: k = GF(2)
            sage: P.<x> = PolynomialRing(k)
            sage: f = x^8+x^4+x^3+x+1 # Rijndael Polynomial
            sage: k.<a> = GF(2**8,modulus=f)
            sage: k.modulus()
            x^8 + x^4 + x^3 + x + 1
            sage: -a^k.degree()
            a^4 + a^3 + a + 1

            You may enforce a random modulus:

            sage: k = GF(3**5,modulus='random')
            sage: k.modulus() # random polynomial
            x^5 + 2*x^4 + 2*x^3 + x^2 + 2

            Three different representations are possible:

            sage: k.<a> = sage.rings.finite_field_givaro.FiniteField_givaro(9,repr='poly')
            sage: a
            a
            sage: k.<a> = sage.rings.finite_field_givaro.FiniteField_givaro(9,repr='int')
            sage: a
            3
            sage: k.<a> = sage.rings.finite_field_givaro.FiniteField_givaro(9,repr='log')
            sage: a
            5
        """

        from sage.rings.polynomial_element import is_Polynomial
        import sage.databases.conway
        from sage.rings.finite_field import conway_polynomial
        from sage.rings.integer import Integer

        cdef intvec cPoly

        if repr=='poly':
            self.repr = 0
        elif repr=='log':
            self.repr = 1
        elif repr=='int':
            self.repr = 2
        else:
            raise ValueError, "Unknown representation %s"%repr

        if q >= 1<<16:
            raise ArithmeticError, "q must be < 2^16"

        q = Integer(q)
        if q < 2:
            raise ArithmeticError, "q  must be a prime power"
        F = q.factor()
        if len(F) > 1:
            raise ArithmeticError, "q must be a prime power"
        p = F[0][0]
        k = F[0][1]

        self.assign_names(name)

        if modulus is None or modulus=="random":
            if k>1 and sage.databases.conway.ConwayPolynomials().has_polynomial(p, k) and modulus!="random":
                modulus = conway_polynomial(p, k)
            else:
                _sig_on
                self.objectptr = gfq_factorypk(p,k)
                _sig_off
                if cache:
                    self._array = self.gen_array()
                return

        if is_Polynomial(modulus):
            modulus = modulus.list()

        if PyObject_TypeCheck(modulus,type_object(list)) or PyObject_TypeCheck(modulus, type_object(tuple)):

            for i from 0 <= i < len(modulus):
                cPoly.push_back(int(modulus[i]))

            _sig_on
            self.objectptr = gfq_factorypkp(p, k,cPoly)
            _sig_off
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
        """
        delete(self.objectptr)

    def __repr__(FiniteField_givaro self):
        if self.degree()>1:
            return "Finite Field in %s of size %d^%d"%(self.variable_name(),self.characteristic(),self.degree())
        else:
            return "Finite Field of size %d"%(self.characteristic())

    def characteristic(FiniteField_givaro self):
        """
        Return integer representing characteristic of the domain.

        """
        return int(self.objectptr.characteristic())

    def order(FiniteField_givaro self):
        """
        Return integer representing cardinality of the domain.

        """
        return self.cardinality()

    cdef int order_c(FiniteField_givaro self):
        """
        C equivalent of self.order()
        """
        return self.objectptr.cardinality()

    def __len__(self):
        """

        len(k) is defined as k.cardinality()

        EXAMPLE:
            sage: k = GF(23**3)
            sage: len(k)
            12167
            sage: k = GF(2)
            sage: len(k)
            2
        """
        return self.cardinality()

    def cardinality(FiniteField_givaro self):
        """
        Return integer representing cardinality of the field.

        """
        return int(self.objectptr.cardinality())

    def degree(FiniteField_givaro self):
        """
        If self.cardinality() == p^n this method returns n.

        """
        return int(self.objectptr.exponent())

    def is_atomic_repr(FiniteField_givaro self):
        if self.repr==0: #modulus
            return False
        else:
            return True

    def is_prime_field(FiniteField_givaro self):
        """
        Return True if self is a prime field

        EXAMPLES:
            sage: GF(3**7).is_prime_field()
            False

        """
        return bool(self.degree()==1)

    def is_prime(FiniteField_givaro self):
        """
        Return True if self has prime cardinality.
        """
        return bool(self.degree()==1)

    def random_element(FiniteField_givaro self):
        """
        Return a random element of self.

        EXAMPLES:
            sage: k = GF(23**3)
            sage: e = k.random_element()
            sage: type(e)
            <type 'finite_field_givaro.FiniteField_givaroElement'>

        """
        cdef int res
        cdef GivRandom generator
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

            sage: k = GF(2**8)
            sage: k.gen() == k(k.gen())
            True


            Floats, ints, longs, Integer are interpreted modulo characteristic

            sage: k(2)
            0
            sage: k(float(2.0))
            0

            Rational are interpreted as
            self(numerator)/self(denominator). Both may not be >=
            self.characteristic().

            sage: k = GF(3**8)
            sage: k(1/2) == k(1)/k(2)
            True

            Free modulo elements over self.prime_subfield() are interpreted 'little endian'

            sage: k = GF(2**8)
            sage: e = k.vector_space().gen(1); e
            (0, 1, 0, 0, 0, 0, 0, 0)
            sage: k(e)
            a

            Strings are evaluated as polynomial representation of elements in self

            sage: k('a^2+1')
            a^2 + 1

            PARI elements are interpreted as finite field elements, the strategy is pretty liberal:

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

        """

        from sage.rings.multi_polynomial_element import MPolynomial
        from sage.rings.polynomial_element import Polynomial
        from sage.modules.free_module_element import FreeModuleElement
        from sage.rings.integer_mod import is_IntegerMod
        from sage.rings.rational import Rational
        from sage.rings.integer import Integer

        cdef int res
        cdef int g
        cdef int x

        ########

        if PyObject_TypeCheck(e, type_object(FiniteField_givaroElement)):
            if e.parent() is self:
                return e
            if e.parent() == self:
                return e
            if e.parent() is self.prime_subfield_C() or e.parent() == self.prime_subfield_C():
                res = self.int2log(int(e))

        elif PyObject_TypeCheck(e, type_object(int)) or \
             PyObject_TypeCheck(e, type_object(Integer)) or \
             PyObject_TypeCheck(e, type_object(long)) or is_IntegerMod(e):
            res = self.objectptr.initi(res,int(e))

        elif PyObject_TypeCheck(e, type_object(float)):
            res = self.objectptr.initd(res,e)

        elif PyObject_TypeCheck(e, type_object(str)):
            return self(eval(e.replace("^","**"),{str(self.variable_name()):self.gen()}))

        elif PyObject_TypeCheck(e, type_object(FreeModuleElement)):
            if self.vector_space() != e.parent():
                raise TypeError, "e.parent must match self.vector_space"
            ret = self.zero()
            for i in range(len(e)):
                ret = ret + self(int(e[i]))*self.gen()**i
            return ret

        elif sage.interfaces.gap.is_GapElement(e):
            return gap_to_givaro(e, self)

        elif PyObject_TypeCheck(e, type_object(MPolynomial)) or PyObject_TypeCheck(e, type_object(Polynomial)):
            if e.is_constant():
                return self(e.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

        elif PyObject_TypeCheck(e, type_object(Rational)):
            num = e.numer()
            den = e.denom()
            if num>=self.characteristic() or den>=self.characteristic():
                raise TypeError, "unable to coerce"
            return self(num)/self(den)

        elif PyObject_TypeCheck(e, type_object(gen)):
            pass # handle this in next if clause

        elif isinstance(e,FiniteField_ext_pariElement):
            # reduce FiniteFieldElements to pari
            e = e._pari_()
        else:
            raise TypeError, "unable to coerce"

        if PyObject_TypeCheck(e, type_object(gen)):
            e = e.lift().lift()
            try:
                res = self.int2log(e[0])
            except TypeError:
                res = self.int2log(e)

            g = self.objectptr.sage_generator()
            x = self.objectptr.one

            for i from 0 < i <= e.poldegree():
                x = self.objectptr.mul(x,x,g)
                res = self.objectptr.axpyin( res, self.int2log(e[i]) , x)

        return make_FiniteField_givaroElement(self,res)

    def _coerce_(self, x):
        """
        Coercion accepts elements of self.parent(), ints, and prime subfield elements.
        """
        #from finite_field.py

        from sage.rings.finite_field_element import FiniteFieldElement
        from sage.rings.integer_mod import is_IntegerMod
        from sage.rings.integer_mod_ring import IntegerModRing_generic
        from sage.rings.integer import Integer

        if PyObject_TypeCheck(x, type_object(int)) \
               or PyObject_TypeCheck(x,type_object(long)) or PyObject_TypeCheck(x, type_object(Integer)):
            return self(x)

        if PyObject_TypeCheck(x, type_object(FiniteFieldElement)) or \
               PyObject_TypeCheck(x, type_object(FiniteField_givaroElement)) or is_IntegerMod(x):
            K = x.parent()
            if K is self:
                return x
            if PyObject_TypeCheck(K, type_object(IntegerModRing_generic)) \
                   and K.characteristic() % self.characteristic() == 0:
                return self(int(x))
            if K.characteristic() == self.characteristic():
                if K.degree() == 1:
                    return self(int(x))
                elif self.degree() % K.degree() == 0:
                    # This is where we *would* do coercion from one nontrivial finite field to another...
                    raise TypeError, 'no canonical coercion defined'
        raise TypeError, 'no canonical coercion defined'


    def one(FiniteField_givaro self):
        """
        Return 1 element in self, which satisfies 1*p=p for every
        element of self != 0.
        """
        return make_FiniteField_givaroElement(self,self.objectptr.one)

    def zero(FiniteField_givaro self):
        """
        Return 0 element in self, which satisfies 0+p=p for every
        element of self.
        """
        return make_FiniteField_givaroElement(self,self.objectptr.zero)


    def gen(FiniteField_givaro self, ignored=None):
        """
        Return a generator of self. All elements x of self are
        expressed as log_{self.gen()}(p) internally. If self is
        a prime field this method returns 1.
        """
        cdef int r
        from sage.rings.arith import primitive_root

        if self.degree() == 1:
            return self(primitive_root(self.order_c()))
        else:
            return make_FiniteField_givaroElement(self,self.objectptr.sage_generator())

##     def multiplicative_generator(FiniteField_givaro self):
##         """
##         Return the multiplicative generator if self.degree() > 1
##         which equals self.gen() and raises NotImplementedError
##         otherwise.
##         """
##         if self.degree() > 1:
##             return self.gen()
##         else:
##             raise NotImplementedError, "multiplicative generator for FiniteField_givaro(p) not implemented"

    cdef prime_subfield_C(FiniteField_givaro self):
        if self._prime_subfield is None:
            self._prime_subfield = FiniteField_givaro(self.characteristic())
        return self._prime_subfield

    def prime_subfield(FiniteField_givaro self):
        """
        Return the prime subfield GF(p) of self if self is GF(p^n)
        """
        return self.prime_subfield_C()

    def base_ring(FiniteField_givaro self):
        """

        Same as prime_subfield.
        """
        return self.prime_subfield_C()

    def log2int(FiniteField_givaro self, int p):
        """
        Given an integer p this method returns i where i satisfies
        self.gen()^p == i. if the result is interpreted as an integer.

        INPUT:
            p -- log representation of a finite field element

        OUTPUT:
            integer representation of a finite field element.

        EXAMPLE:
            sage: k = GF(2**8)
            sage: k.log2int(4)
            16
            sage: k.log2int(20)
            180

        """
        cdef int ret

        if p<0:
            raise ArithmeticError, "Cannot serve negative exponent %d"%p
        elif p>=self.order_c():
            raise IndexError, "p=%d must be < self.order()"%p
        _sig_on
        ret = int(self.objectptr.write(ret, p))
        _sig_off
        return ret

    def int2log(FiniteField_givaro self, int p):
        """
        Given an integer p this method returns i where i satisfies
        self.gen()^i==(p\%self.characteristic())

        INPUT:
            p -- integer representation of an finite field element

        OUTPUT:
            log representation of p

        EXAMPLE:
        sage: k = GF(7**3)
        sage: k.int2log(4)
        228
        sage: k.int2log(3)
        57
        sage: k.gen()^57
        3
        """
        cdef int r
        _sig_on
        ret =  int(self.objectptr.read(r,p))
        _sig_off
        return ret

    def polynomial(self):
        """
        Minimal polynomial of self in self.polynomial_ring().
        """

        quo = int(-(self.gen()**(self.degree())))
        b   = int(self.characteristic())

        ret = []
        for i in range(self.degree()):
            ret.append(quo%b)
            quo = quo/b
        ret = ret + [1]
        from sage.rings.polynomial_ring import PolynomialRing
        R = PolynomialRing(self.prime_subfield_C())
        return R(ret)

    def modulus(self):
        """
        Minimal polynomial of self in self.polynomial_ring().
        """
        return self.polynomial()

    def _pari_modulus(self):
        """
        """
        f = pari(str(self.modulus()))
        return f.subst('x', 'a') * pari("Mod(1,%s)"%self.characteristic())

    cdef polynomial_ring_C(self):
        """
        """
        if self._polynomial_ring is None:
            from sage.rings.polynomial_ring import PolynomialRing
            self._polynomial_ring = PolynomialRing(self.prime_subfield_C(),self.variable_name())
            return self._polynomial_ring
        else:
            return self._polynomial_ring

    def polynomial_ring(self):
        """
        Return the polynomial ring over the prime subfield in the
        same variable as this finite field.
        """
        return self.polynomial_ring_C()

    def _finite_field_ext_pari_(self):
        """
        Return  a FiniteField_ext_pari matching self.
        """

        from sage.rings.finite_field import FiniteField_ext_pari
        from sage.rings.finite_field import FiniteField_prime_modn
        return FiniteField_ext_pari(self.order_c(),self.variable_name(),self.polynomial())

    def vector_space(FiniteField_givaroElement self):
         """

         Returns self interpreted as a VectorSpace over
         self.prime_subfield()

         EXAMPLE:
             sage: k = GF(3**5)
             sage: k.vector_space()
             Vector space of dimension 5 over Finite Field of size 3

         """
         import sage.modules.all
         V = sage.modules.all.VectorSpace(self.prime_subfield(),self.degree())
         return V

    def __iter__(FiniteField_givaro self):
        """
        Finite fields may be iterated over:

        EXAMPLE:
            sage: list(GF(2**2))
            [0, 1, a, a + 1]
        """
        if self.degree()>1:
            return FiniteField.__iter__(self)
        else:
            return FiniteField_givaro_iterator(self)

    def __cmp__(FiniteField_givaro self, other):
        """
        Finite Fields are considered to be equal if
         * their implementation is the same (Givaro)
         * their characteristics match
         * their degrees match
         * their moduli match if degree > 1

        This will probably change in the future so that different
        implementations are equal.

        """
        if not PyObject_TypeCheck(other,type_object(FiniteField_givaro)):
            return 1
        if self.characteristic()!=other.characteristic():
            return 1
        if self.degree()!=other.degree():
            return 1
        if self.degree()>1:
            if self.polynomial()!=other.polynomial():
                return 1
            if self.variable_name()!=other.variable_name():
                return 1
        return 0

    def __hash__(FiniteField_givaro self):
        """
        A hash of a Givaro finite field is a hash over it's
        characterstic, polynomial and the string 'givaro'
        """
        if self.degree()>1:
            return hash((self.characteristic(),self.polynomial(),"givaro"))
        else:
            return hash((self.characteristic(),"givaro"))

    def _element_repr(FiniteField_givaro self, FiniteField_givaroElement e):
        """
        Wrapper for log, int, and poly representations.
        """
        if self.repr==0:
            return self._element_poly_repr(e)
        elif self.repr==1:
            return self._element_log_repr(e)
        else:
            return self._element_int_repr(e)

    def _element_log_repr(FiniteField_givaro self, FiniteField_givaroElement e):
        """
        Return str(i) where base.gen()^i=self
        """
        return str(int(e.object))

    def _element_int_repr(FiniteField_givaro self, FiniteField_givaroElement e):
        """
	elements of this field will be written in the following
        manner: for e in ZZp[x] with e = a0 + a1x + a2x^2 + ..., e is
        represented as: 'n' where n = a0 + a1 * p + a2 * p^2 + ...

        """
        return str(int(e))

    def _element_poly_repr(FiniteField_givaro self, FiniteField_givaroElement e):
        """
        Return a polynomial expression in base.gen() of self.
        """
        variable = self.variable_name()

        quo = self.log2int(e.object)
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

    def abpc(FiniteField_givaro self,FiniteField_givaroElement a, FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return r = ab + c. This is faster than multiplying a and b
        first and adding c to the result.

        INPUT:
            a -- FiniteField_givaroElement
            b -- FiniteField_givaroElement
            c -- FiniteField_givaroElement

        EXAMPLE:
            sage: k.<a> = GF(2**8)
            sage: k.abpc(a,a,k(1))
            a^2 + 1

        """
        cdef int r

        r = self.objectptr.axpy(r, a.object, b.object, c.object)
        return make_FiniteField_givaroElement(self,r)

    def abmc(FiniteField_givaro self,FiniteField_givaroElement a, FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return r = ab - c.

        INPUT:
            a -- FiniteField_givaroElement
            b -- FiniteField_givaroElement
            c -- FiniteField_givaroElement

        EXAMPLE:
            sage: k.<a> = GF(3**3)
            sage: k.abmc(a,a,k(1))
            a^2 + 2

        """

        cdef int r

        r = self.objectptr.axmy(r, a.object, b.object, c.object, )
        return make_FiniteField_givaroElement(self,r)

    def cmab(FiniteField_givaro self,FiniteField_givaroElement a, FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return r = c - ab.


        INPUT:
            a -- FiniteField_givaroElement
            b -- FiniteField_givaroElement
            c -- FiniteField_givaroElement

        EXAMPLE:
            sage: k.<a> = GF(3**3)
            sage: k.cmab(a,a,k(1))
            2*a^2 + 1
        """
        cdef int r

        r = self.objectptr.amxy(r , a.object, b.object, c.object, )
        return make_FiniteField_givaroElement(self,r)

##     def _add(FiniteField_givaro self, int r, int l):
##         """
##         This is the fastest way to add two Givaro finite field
##         elements using SAGE. Given r and l this method calculates s
##         such that self.gen()^s = self.gen()^r + self.gen()^l.

##         INPUT:
##             r -- int representing an exponent of self.gen()
##             l -- int representing an exponent of self.gen()

##         EXAMPLE:
##             sage: k.<a> = GF(2**8)
##             sage: k._add(int(10),int(20))
##             31
##             sage: (a^10+a^20).log_repr()
##             '31'
##         """
##         cdef int res
##         return self.objectptr.add(res, r , l )

##     def _mul(FiniteField_givaro self, int r, int l):
##         """
##         This is the fastest way to multiply two Givaro finite field
##         elements using SAGE. Given r and l this method calculates s
##         such that self.gen()^s = self.gen()^r * self.gen()^l.

##         INPUT:
##             r -- int representing an exponent of self.gen()
##             l -- int representing an exponent of self.gen()

##         EXAMPLE:
##             sage: k.<a> = GF(2**8)
##             sage: k._mul(int(10),int(20))
##             30
##             sage: (a^10*a^20).log_repr()
##             '30'

##         """
##         cdef int res
##         return self.objectptr.mul(res, r , l )

##     def _div(FiniteField_givaro self, int r, int l):
##         """
##         This is the fastest way to divide two Givaro finite field
##         elements using SAGE. Given r and l this method calculates s
##         such that self.gen()^s = self.gen()^r / self.gen()^l.

##         INPUT:
##             r -- int representing an exponent of self.gen()
##             l -- int representing an exponent of self.gen()

##         EXAMPLE:
##             sage: k.<a> = GF(2**8)
##             sage: k._div(int(10),int(20))
##             245
##             sage: (a^10/a^20).log_repr()
##             '245'

##         """
##         cdef int res
##         return self.objectptr.div(res, r , l )

##     def _sub(FiniteField_givaro self, int r, int l):
##         """
##         This is the fastest way to subtract two Givaro finite field
##         elements using SAGE. Given r and l this method calculates s
##         such that self.gen()^s = self.gen()^r + self.gen()^l.

##         INPUT:
##             r -- int representing an exponent of self.gen()
##             l -- int representing an exponent of self.gen()

##         EXAMPLE:
##             sage: k.<a> = GF(2**8)
##             sage: k._sub(int(10),int(20))
##             31
##             sage: (a^10-a^20).log_repr()
##             '31'
##         """
##         cdef int res
##         return self.objectptr.sub(res, r , l )

    def __reduce__(FiniteField_givaro self):
        """

        Pickle self:

        EXAMPLE:
            sage: k = GF(2**8)
            sage: loads(dumps(k)) == k
            True

        """
        if self._array is None:
            cache = 0
        else:
            cache = 1

        return sage.rings.finite_field_givaro.unpickle_FiniteField_givaro, \
               (self.order_c(),self.variable_name(),map(int,list(self.modulus())),int(self.repr),cache)

def unpickle_FiniteField_givaro(order,variable_name,modulus,rep,cache):
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
    """
    cdef int iterator
    cdef FiniteField_givaro _parent

    def __init__(self, FiniteField_givaro parent):
        self._parent = parent
        self.iterator = -1

    def __next__(self):
        """
        """

        self.iterator=self.iterator+1

        if self.iterator==self._parent.characteristic():
            self.iterator = -1
            raise StopIteration

        return make_FiniteField_givaroElement(self._parent,self._parent.int2log(self.iterator))

    def __repr__(self):
        return "Iterator over %s"%self._parent

cdef FiniteField_givaro_copy(FiniteField_givaro orig):
    cdef FiniteField_givaro copy
    copy = FiniteField_givaro(orig.characteristic()**orig.degree())
    delete(copy.objectptr)
    copy.objectptr = gfq_factorycopy(gfq_deref(orig.objectptr))
    return copy

cdef class FiniteField_givaroElement(FiniteFieldElement):
    """
    Element in FiniteField_givaro.
    """
    cdef int object
    cdef object __multiplicative_order

    def __init__(FiniteField_givaroElement self, SageObject parent ):
        """
        Initializes an element in parent. It's much better to use
        parent(<value>) or any specialized method of parent
        (one,zero,gen) instead.

        Alternatively you may provide a value which is directly
        assigned to this element. So the value must represent the
        log_g of the value you wish to assign.

        INPUT:
            parent -- base field

        OUTPUT:
            finite field element.
        """
        self._parent = parent
        self.object = 0

    def __dealloc__(FiniteField_givaroElement self):
        pass

    def __repr__(FiniteField_givaroElement self):
        return (<FiniteField_givaro>self._parent)._element_repr(self)

    def parent(self):
        """
        Return parent finite field.
        """
        return (<FiniteField_givaro>self._parent)

    def is_zero(FiniteField_givaroElement self):
        """
        Return True if self == k(0).
        """
        return bool((<FiniteField_givaro>self._parent).objectptr.isZero(self.object))

    def is_one(FiniteField_givaroElement self):
        """
        Return True if self == k(1)
        """
        return bool((<FiniteField_givaro>self._parent).objectptr.isOne(self.object))

    def is_unit(FiniteField_givaroElement self):
        """
        Return True if self is an element of the prime subfield.
        """
        return bool((<FiniteField_givaro>self._parent).objectptr.isunit(self.object))


    def is_square(FiniteField_givaroElement self):
        """
        Return True if self is a square in self.parent()
        """
        #copied from finite_field_element.py
        cdef FiniteField_givaro K
        K = (<FiniteField_givaro>self._parent)
        if K.characteristic() == 2:
            return True
        n = K.order_c() - 1
        a = self**(n / 2)
        return bool(a == 1)


    def __add__(self, other):
        """
        Add two elements.

        EXAMPLE:
            sage: k.<b> = GF(9**2)
            sage: b^10 + 2*b
            2*b^3 + 2*b^2 + 2*b + 1

        """
        cdef int r

        if not PyObject_TypeCheck(self, type_object(FiniteField_givaroElement)):
            other,self = self,other

        if not PyObject_TypeCheck(other, type_object(FiniteField_givaroElement)):
            return bin_op(self,other,operator.add)

        else:
            r = parent_object(self).objectptr.add(r,
                                                  (<FiniteField_givaroElement>self).object ,
                                                  (<FiniteField_givaroElement>other).object )


            return make_FiniteField_givaroElement(parent_object(self),r)

    def __mul__(self, other):
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


        if not PyObject_TypeCheck(self, type_object(FiniteField_givaroElement)):
            other,self = self,other

        if not PyObject_TypeCheck(other,type_object(FiniteField_givaroElement)):
            return bin_op(self,other,operator.mul)
        else:
            r = parent_object(self).objectptr.mul(r,
                                                  (<FiniteField_givaroElement>self).object,
                                                  (<FiniteField_givaroElement>other).object)
            return make_FiniteField_givaroElement(parent_object(self),r)

    def __div__(self, other):
        """
        Divide two elements

        EXAMPLE:
            sage: k.<g> = GF(2**8)
            sage: g/g
            1
        """
        cdef int r

        if not PyObject_TypeCheck(self, type_object(FiniteField_givaroElement)):
            other,self = self,other

        if not PyObject_TypeCheck(other, type_object(FiniteField_givaroElement)):
            return bin_op(self,other,operator.div)
        else:
            r = parent_object(self).objectptr.div(r,
                                                  (<FiniteField_givaroElement>self).object,
                                                  (<FiniteField_givaroElement>other).object)
            return make_FiniteField_givaroElement(parent_object(self),r)

    def __sub__(self, other):
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

        if not PyObject_TypeCheck(self, type_object(FiniteField_givaroElement)):
            other,self = self,other

        if not PyObject_TypeCheck(other, type_object(FiniteField_givaroElement)):
            return bin_op(self,other,operator.sub)
        else:
            r = parent_object(self).objectptr.sub(r,
                                                  (<FiniteField_givaroElement>self).object,
                                                  (<FiniteField_givaroElement>other).object)
            return make_FiniteField_givaroElement(parent_object(self),r)

    def __neg__(FiniteField_givaroElement self):
        cdef int r

        r = (<FiniteField_givaro>self._parent).objectptr.neg(r, self.object)
        return make_FiniteField_givaroElement((<FiniteField_givaro>self._parent),r)

    def __invert__(FiniteField_givaroElement self):
        cdef int r

        (<FiniteField_givaro>self._parent).objectptr.inv(r, self.object)
        return make_FiniteField_givaroElement((<FiniteField_givaro>self._parent),r)


    def __pow__(FiniteField_givaroElement self, int exp, other):
        #There doesn't seem to exist a power function for FiniteField_givaro. So we
        #had to write one. It is pretty clumbsy (read: slow) right now

        cdef int power
        cdef int i
        cdef int epow2
        cdef GivaroGfq *field

        field = (<FiniteField_givaro>self._parent).objectptr

        exp = exp % int(field.cardinality()-1)

        if field.isOne(self.object):
            return self

        if exp==0:
            return make_FiniteField_givaroElement((<FiniteField_givaro>self._parent),field.one)

        power = field.one
        i = 0;
        epow2 = self.object;
        while (exp>>i) > 0:
            if (exp>>i) & 1:
                field.mulin(power,epow2)
            field.mulin(epow2,epow2)
            i = i + 1

        return make_FiniteField_givaroElement((<FiniteField_givaro>self._parent),power)

    def __cmp__(self, other):
        """
        Comparison of finite field elements is performed by comparing
        their int representation.
        """
        if PyObject_TypeCheck(self,type_object(FiniteField_givaroElement)) and \
               PyObject_TypeCheck(other,type_object(FiniteField_givaroElement)):
            if not parent_object(self) is parent_object(other):
                if not parent_object(self) == parent_object(other):
                    return 1
            return cmp(int(self),int(other))
        return cmp(int(self),int(other))

    def __richcmp__(self, right, int op):
        cdef int res

        if PyObject_TypeCheck(self,type_object(FiniteField_givaroElement)) and \
           PyObject_TypeCheck(right,type_object(FiniteField_givaroElement)):
            return (<FiniteField_givaroElement>self).richcmp_C(<FiniteField_givaroElement>right,op)
        else:
            res = cmp(int(self),int(right))
            if op == 0:  #<
                return bool(res < 0)
            if op == 2: #==
                return bool(res == 0)
            if op == 4: #>
                return bool(res > 0)
            if op == 1: #<=
                return bool(res <= 0)
            if op == 3: #!=
                return bool(res != 0)
            if op == 5: #>=
                return bool(res >= 0)

    cdef richcmp_C(FiniteField_givaroElement self, FiniteField_givaroElement right, int op):

        if not parent_object(self) is parent_object(right):
            if not parent_object(self) == parent_object(right):
                if not self.is_zero() and right.is_zero():
                    if op == 3:
                        return True
                    else:
                        return False

        if op == 0:  #<
            return bool(int(self) < int(right))
        if op == 2: #==
            return bool(int(self) == int(right))
        if op == 4: #>
            return bool(int(self) > int(right))
        if op == 1: #<=
            return bool(int(self) <= int(right))
        if op == 3: #!=
            return bool(int(self) != int(right))
        if op == 5: #>=
            return bool(int(self) >= int(right))

    def __int__(FiniteField_givaroElement self):
        """
        Return self coerced to an int. The integer returned is
        equivalent to the representation of self and not to log_repr.
        """
        return (<FiniteField_givaro>self._parent).log2int(self.object)


    def logint(FiniteField_givaroElement self):
        """
        Return i where base.gen()^i=self
        """
        return int(self.object)

    def log(FiniteField_givaroElement self, a):
        #copied from finite_field_element.py
        q = (self.parent()).order_c() - 1
        return sage.rings.arith.discrete_log_generic(self, a, q)

    def int_repr(FiniteField_givaroElement self):
        return (<FiniteField_givaro>self._parent)._element_int_repr(self)

    def log_repr(FiniteField_givaroElement self):
        return (<FiniteField_givaro>self._parent)._element_log_repr(self)

    def poly_repr(FiniteField_givaroElement self):
        return (<FiniteField_givaro>self._parent)._element_poly_repr(self)

    def polynomial(FiniteField_givaroElement self):
        """
        Return self as polynomial over self.parent().prime_subfield()
        """
        quo = (<FiniteField_givaro>self._parent).log2int(self.object)
        b   = int((<FiniteField_givaro>self._parent).characteristic())

        ret = []
        for i in range((<FiniteField_givaro>self._parent).degree()):
            coeff = quo%b
            ret.append(coeff)
            quo = quo/b
        return (<FiniteField_givaro>self._parent).polynomial_ring()(ret)

    def _latex_(FiniteField_givaroElement self):
        if (<FiniteField_givaro>self._parent).degree()>1:
            return self.polynomial()._latex_()
        else:
            return str(self)

    def _finite_field_ext_pari_(FiniteField_givaroElement self, k=None):
        """
        Return an element of k supposed to match this element.  No
        checks if k equals self.parent() are performed.

        INPUT:
            k -- FiniteField_ext_pari

        OUTPUT:
            k.gen()^(self.log_repr())

        """
        if k is None:
            k=(<FiniteField_givaro>self._parent)._sage_()

        variable = k.gen()._pari_()

        quo = int(self)
        b   = (<FiniteField_givaro>self._parent).characteristic()

        ret = k._pari_one() - k._pari_one()
        i = 0
        while quo!=0:
            coeff = quo%b
            if coeff != 0:
                ret = coeff * variable ** i + ret
            quo = quo/b
            i = i+1
        return k(ret)

    def _pari_init_(FiniteField_givaroElement self):
        """
        """
        k=(parent_object(self))._finite_field_ext_pari_()

        variable = k.gen()._pari_()

        quo = int(self)
        b   = (parent_object(self)).characteristic()

        ret = k._pari_one() - k._pari_one() # there is no pari_zero
        i = 0
        while quo!=0:
            coeff = quo%b
            if coeff != 0:
                ret = coeff * variable ** i + ret
            quo = quo/b
            i = i+1
        return ret

    def multiplicative_order(FiniteField_givaroElement self):
        """
        """
        # code copy'n'pasted from finite_field_element.py
        import sage.rings.arith
        from sage.rings.integer import Integer

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

    def copy(self):
        """
        """
        return make_FiniteField_givaroElement((parent_object(self)),self.object)

    def _gap_init_(FiniteField_givaroElement self):
        """
        This is only correct if self.parent().modulus() is a conway
        polynomial as the isomorphism of finite fields is not
        implemented yet.

        """
        #copied from finite_field_element.py
        cdef FiniteField_givaro F
        F = parent_object(self)
        if F.order_c() > 65536:
            raise TypeError, "order (=%s) must be at most 65536."%F.order_c()
        if self == 0:
            return '0*Z(%s)'%F.order_c()
        assert F.degree() > 1
        g = F.multiplicative_generator()
        n = g.log(self)
        return 'Z(%s)^%s'%(F.order_c(), n)

    def charpoly(FiniteField_givaroElement self):
        """
        Return characteristic polynomial of self.

        """
        from sage.rings.polynomial_ring import PolynomialRing
        R = PolynomialRing(parent_object(self).prime_subfield_C())
        return R(self._pari_().charpoly().lift())


    def norm(FiniteField_givaroElement self):
        """
        Return norm of self.

        """
        return self.charpoly()[0]

    def trace(FiniteField_givaroElement self):
        """
        Return norm of self.

        """
        return parent_object(self).prime_subfield()(self._pari_().trace().lift())

    def __hash__(FiniteField_givaroElement self):
        # GF elements are hashed by hashing their string
        # representation but string representations are slow. So we
        # hash the log and the int representation which should provide
        # the same level of distinction.
        return hash((self.object,(parent_object(self)).log2int(self.object),"givaro"))

    def square_root(FiniteField_givaroElement self):
        """
        Return a square root of this finite field element in its
        finite field, if there is one.  Otherwise, raise a ValueError.

        EXAMPLES:
          sage: k.<a> = GF(7^2)
          sage: k(2).square_root()
          4
          sage: k(3).square_root()
          5*a + 1
          sage: k(3).square_root()**2
          3
          sage: k(4).square_root()
          5
          sage: k = GF(7^3)
          sage: k(3).square_root()
          Traceback (most recent call last):
          ...
          ValueError: must be a perfect square.

        """
        from sage.rings.polynomial_ring import PolynomialRing
        R = PolynomialRing(parent_object(self))
        f = R([-self, 0, 1])
        g = f.factor()
        if len(g) == 2 or g[0][1] == 2:
            return -g[0][0][0]
        raise ValueError, "must be a perfect square."

    def __reduce__(FiniteField_givaroElement self):
        """

        EXAMPLE:
            sage: k.<a> = GF(2**8)
            sage: e = k.random_element()
            sage: loads(dumps(e)) == e
            True
        """
        return unpickle_FiniteField_givaroElement,(parent_object(self),self.object)


def unpickle_FiniteField_givaroElement(FiniteField_givaro parent, int x):
    return make_FiniteField_givaroElement(parent, x)

cdef make_FiniteField_givaroElement(FiniteField_givaro parent, int x):
    """
    """
    cdef FiniteField_givaroElement y
    if parent._array is None:
        y = FiniteField_givaroElement(parent)
        y.object = x
        return y
    else:
        return parent._array[x]

cdef gap_to_givaro(x, FiniteField_givaro F):
    """
    INPUT:
        x -- gap finite field element
        F -- Givaro finite field
    OUTPUT:
        element of F

    EXAMPLES:
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

    AUTHOR:
        -- David Joyner and William Stein
        -- Martin Albrecht (copied from gap_to_sage)
    """
    import sage.interfaces.gap
    s = str(x)
    if s[:2] == '0*':
        return F(0)
    i1 = s.index("(")
    i2 = s.index(")")
    q  = eval(s[i1+1:i2].replace('^','**'))
    if q == F.order_c():
        K = F
    else:
        K = FiniteField_givaro(q)
    if s.find(')^') == -1:
        e = 1
    else:
        e = int(s[i2+2:])

    if F.degree() == 1:
        g = int(sage.interfaces.gap.gap.eval('Int(Z(%s))'%q))
    else:
        g = K.multiplicative_generator()
    return F(K(g**e))

