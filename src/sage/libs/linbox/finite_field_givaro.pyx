r"""nodoctest
 AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-06-05)

"""


#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
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


from sage.ext.ring cimport FiniteField
from sage.rings.coerce import bin_op
from sage.ext.sage_object cimport SageObject
import operator
import sage.rings.arith

import sage.interfaces.gap
from sage.libs.pari.all import pari
from sage.libs.pari.gen import gen

## cdef extern from 'interrupt.h':
##     int _sig_on, _sig_off, _sig_check
##     void _sig_str(char*)
cdef int _sig_on
cdef int _sig_off
cdef int _sig_check

cdef extern from "givaro/givrandom.h":
    ctypedef struct GivRandom "GivRandom":
        pass

cdef extern from "linbox/integer.h":
    ctypedef struct intvec "std::vector<LinBox::integer>":
        void (* push_back)(int elem)

    intvec intvec_factory "std::vector<LinBox::integer>"(int len)

cdef extern from "linbox/field/givaro-gfq.h":

    ctypedef struct GivaroGfq "LinBox::GivaroGfq":
        #attributes
        int one
        int zero

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
        int (* gen)(int r)
        int (* write)(int r, int p)
        int (* read)(int r, int p)
        int (* axpy)(int r, int a, int b, int c)
        int (* axmy)(int r, int a, int b, int c)
        int (* amxy)(int r, int a, int b, int c)
        int (* isZero)(int e)
        int (* isOne)(int e)
        int (* isunit)(int e)

    GivaroGfq *gfq_factorypk "new LinBox::GivaroGfq" (int p, int k)
    GivaroGfq *gfq_factorypkp "new LinBox::GivaroGfq" (int p, int k, intvec poly)
    GivaroGfq *gfq_factorycopy "new LinBox::GivaroGfq"(GivaroGfq orig)
    GivaroGfq  gfq_deref "*"(GivaroGfq *orig)
    void delete "delete "(void *o)
    int gfq_element_factory "LinBox::GivaroGfq::Element"()

cdef class GFq_element(SageObject) # forward declaration

cdef class GFq(FiniteField):
    """
    Givaro Finite Field. These are implemented using Zech logs and
    therefor the cardinality must be < 2^20.
    """

    cdef GivaroGfq *objectptr
    cdef int repr

    def __init__(GFq self, q, name="a",  poly=None, repr="poly"):
        """
        Givaro Finite Field. These are implemented using Zech logs and
        therefor the cardinality must be < 2^20.

        INPUT:
            q     -- p^n (must be prime power)
            name  -- variable used for poly_repr (default: 'a')
            poly  -- you may provide a minimal polynomial to use for
                     reduction. (default: None, so a random polynomial will be
                     used by the Givaro library)
            repr  -- controls the way elements are printed to the user:
                     (default: 'poly')
                     'log': repr is element.log_repr()
                     'int': repr is element.int_repr()
                     'poly': repr is element.poly_repr()


        OUTPUT:
            Givaro finite field with characteristic p and cardinality p^n.
        """

        from sage.rings.polynomial_element import is_Polynomial
        import sage.databases.conway
        from sage.rings.finite_field import conway_polynomial
        from sage.ext.integer import Integer

        cdef intvec cPoly

        if repr=='poly':
            self.repr = 0
        elif repr=='log':
            self.repr = 1
        elif repr=='int':
            self.repr = 2
        else:
            raise ValueError, "Unknown representation %s"%repr

        q = Integer(q)
        if q < 2:
            raise ArithmeticError, "q (=%s) must be a prime power"%q
        F = q.factor()
        if len(F) > 1:
            raise ArithmeticError, "q (=%s) must be a prime power"%q
        p = F[0][0]
        k = F[0][1]

        self.assign_names(name)

        if poly==None or poly=="random":
            if k>1 and sage.databases.conway.ConwayPolynomials().has_polynomial(p, k) and poly!="random":
                poly = conway_polynomial(p, k)
            else:
                _sig_on
                self.objectptr = gfq_factorypk(p,k)
                _sig_off
                return

        if is_Polynomial(poly):
            poly = poly.list()

        if isinstance(poly,list) or isinstance(poly,tuple):
            #cPoly = intvec_factory(len(poly))

            for i from 0 <= i < len(poly):
                cPoly.push_back(int(poly[i]))

            _sig_on
            self.objectptr = gfq_factorypkp(p,k,cPoly)
            _sig_off
            return

        raise TypeError, "Cannot understand %s"%(poly)



    def __dealloc__(GFq self):
        """
        """

        delete(self.objectptr)

    def __repr__(GFq self):
        if self.degree()>1:
            return "Finite Field in %s of size %d^%d"%(self.variable_name(),self.characteristic(),self.degree())
        else:
            return "Finite Field of size %d"%(self.characteristic())

    def characteristic(GFq self):
        """
        Return integer representing characteristic of the domain.
        Returns a positive integer to all domains with finite
        characteristic, and returns 0 to signify a domain of infinite
        characteristic.

        OUTPUT:
            integer representing characteristic of the domain.
        """
        return int(self.objectptr.characteristic())

    def order(GFq self):
        return self.cardinality()

    def __len__(self):
        return self.cardinality()

    def cardinality(GFq self):
        """
        Return integer representing cardinality of the domain.
        Returns a non-negative integer for all domains with finite
        cardinality, and returns -1 to signify a domain of infinite
        cardinality.

        OUTPUT:
            integer representing cardinality of the domain
        """
        return int(self.objectptr.cardinality())

    def degree(GFq self):
        """
        If self.cardinality() == p^n this method returns n.

        OUTPUT:
             log_{self.characteristic()}(self.cardinality())
        """
        return int(self.objectptr.exponent())

    def is_atomic_repr(GFq self):
        if self.repr==0: #poly
            return False
        else:
            return True

    def is_prime_field(GFq self):
        """
        """
        return bool(self.degree()==1)

    def random_element(GFq self):
        """
        """

        cdef int res
        cdef GivRandom generator
        res = self.objectptr.random(generator,res)
        return make_GFq_element(self,res)

    def __call__(GFq self, e):
        """
        Coerces several data types to self.
        """

        from sage.rings.multi_polynomial_element import MPolynomial
        from sage.rings.polynomial_element import Polynomial
        from sage.modules.free_module_element import FreeModuleElement
        from sage.rings.finite_field_element import FiniteFieldElement
        from sage.rings.integer_mod import IntegerMod
        from sage.rings.rational import Rational
        from sage.ext.integer import Integer


        cdef int res
        cdef int g
        cdef int x

        res = gfq_element_factory()

        if isinstance(e,FiniteFieldElement):
            # reduce FiniteFieldElements to pari
            e = e._pari_()

        ########

        if isinstance(e, (Integer,int,long,IntegerMod)):
            res = self.objectptr.initi(res,int(e))

        elif isinstance(e,float):
            res = self.objectptr.initd(res,e)

        elif isinstance(e,gen):
            e = e.lift().lift()
            try:
                res = self.int2log(e[0])
            except TypeError:
                res = self.int2log(e)

            g = self.objectptr.gen(g)

            for i from 0 < i <= e.poldegree():
                #x = GFq_element_pow(self.objectptr, g, i);
                res = self.objectptr.axpyin( res, self.int2log(e[i]) , x)

        elif isinstance(e,GFq_element):
            if e.parent() == self:
                return e
            if e.parent() == self.prime_subfield():
                res = self.int2log(int(e))

        elif isinstance(e,str):
            return self(eval(e.replace("^","**"),{str(self.variable_name()):self.gen()}))

        elif isinstance(e,FreeModuleElement):
            ret = self.zero()
            for i in range(len(e)):
                ret = ret + self(int(e[i]))*self.gen()**i
            return ret

        elif sage.interfaces.gap.is_GapElement(e):
            return gap_to_givaro(e, self)

        elif isinstance(e, (MPolynomial, Polynomial)):
            if e.is_constant():
                return self(e.constant_coefficient())
            else:
                raise TypeError, "no coercion of non-constant polynomial %s into %s defined."%(e,self)

        elif isinstance(e, Rational):
            num = e.numer()
            den = e.denom()
            if num>=self.characteristic() or den>=self.characteristic():
                raise TypeError, "Unable to coerce %s into %s."%(e,self)
            return self(num)/self(den)
        else:
            raise TypeError, "Unable to coerce %s into %s."%(e,self)

        return make_GFq_element(self,res)

    def _coerce_(self, x):
        #from finite_field.py

        from sage.rings.finite_field_element import FiniteFieldElement
        from sage.rings.integer_mod import IntegerMod
        from sage.rings.integer_mod_ring import IntegerModRing_generic
        from sage.ext.integer import Integer

        if isinstance(x, (int, long, Integer)):
            return self(x)

        if isinstance(x, (FiniteFieldElement,IntegerMod,GFq_element)):
            K = x.parent()
            if K is self:
                return x
            if isinstance(K, IntegerModRing_generic) and K.characteristic() % self.characteristic() == 0:
                return self(int(x))
            if K.characteristic() == self.characteristic():
                if K.degree() == 1:
                    return self(int(x))
                elif self.degree() % K.degree() == 0:
                    # This is where we *would* do coercion from one nontrivial finite field to another...
                    raise TypeError, 'no canonical coercion of %s to %s defined; adding this is planned (see rings/finite_field.py to help!)'%(x, self)
        raise TypeError, 'no canonical coercion of %s to %s'%(x, self)


    def one(GFq self):
        """
        Returns 1 element in self, which satisfies 1*p=p for
        every element of self != 0.
        """
        return make_GFq_element(self,self.objectptr.one)

    def zero(GFq self):
        """
        Returns 0 element in self, which satisfies 0+p=p for
        every element of self.
        """
        return make_GFq_element(self,self.objectptr.zero)


    def gen(GFq self, ignored=None):
        """
        Returns a generator of self. All elements x of self are
        expressed as log_{self.gen()}(p) internally.
        """
        cdef int r
        from sage.rings.arith import primitive_root

        if self.degree() == 1:
            return self(primitive_root(self.order()))
        else:
            return make_GFq_element(self,self.objectptr.gen(r))

    def multiplicative_generator(GFq self):
        """
        """
        return self.gen()

    def prime_subfield(GFq self):
        return self.base_ring()

    def base_ring(GFq self):
        return GFq(self.characteristic())

    def log2int(GFq self, int p):
        """
        Given an integer p this method returns i where i satisfies
        self.gen()^p == i.

        INPUT:
            p -- log representation of a finite field element

        OUTPUT:
            integer representation of a finite field element.
        """
        cdef int ret

        if p<0:
            raise ArithmeticError, "Cannot serve negative exponent %d"%p
        elif p>=self.order():
            raise IndexError, "p=%d must be < self.order()"%p
        _sig_on
        ret = int(self.objectptr.write(ret, p))
        _sig_off
        return ret

    def int2log(GFq self, int p):
        """
        Given an integer p this method returns i where i satisfies
        self.gen()^i==(p\%self.characteristic())

        INPUT:
            p -- integer representation of an finite field element

        OUTPUT:
            log representation of p
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
        R = PolynomialRing(GFq(self.characteristic()))

        return R(ret)

    def modulus(self):
        """
        Minimal polynomial of self in self.polynomial_ring().
        """
        return self.polynomial()

    def _pari_modulus(self):
        """
        """
        return self._sage_()._pari_modulus()

    def polynomial_ring(self):
        """
        Returns the polynomial ring over the prime subfield in the
        same variable as this finite field.
        """
        from sage.rings.polynomial_ring import PolynomialRing
        return PolynomialRing(GFq(self.characteristic()),self.variable_name())

    def _sage_(self):
        """
        Returns a SAGE Finite Field (which is a FiniteField_ext_pari)
        matching self.
        """

        from sage.rings.finite_field import FiniteField_ext_pari
        from sage.rings.finite_field import FiniteField_prime_modn
        if self.degree()>1:
            return FiniteField_ext_pari(self.cardinality(),self.variable_name(),self.polynomial())
        else:
            return FiniteField_prime_modn(self.cardinality(),self.variable_name())

    def vector_space(GFq_element self):
         """
         """
         import sage.modules.all
         V = sage.modules.all.VectorSpace(self.prime_subfield(),self.degree())
         return V

    def __iter__(GFq self):
        if self.degree()>1:
            return FiniteField.__iter__(self)
        else:
            return GFq_iterator(self)

    def __cmp__(GFq self, other):
        if not isinstance(other,GFq):
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

    def __hash__(GFq self):
        """
        """
        if self.degree()>1:
            return hash((self.characteristic(),self.polynomial(),"givaro"))
        else:
            return hash((self.characteristic(),"givaro"))

    def _element_repr(GFq self, GFq_element e):
        """
        Wrapper for log, int, and poly representations.
        """
        if self.repr==0:
            return self._element_poly_repr(e)
        elif self.repr==1:
            return self._element_log_repr(e)
        else:
            return self._element_int_repr(e)

    def _element_log_repr(GFq self, GFq_element e):
        """
        Returns str(i) where base.gen()^i=self
        """
        return str(int(e.object))

    def _element_int_repr(GFq self, GFq_element e):
        """
	elements of this field will be written in the following
        manner: for e in ZZp[x] with e = a0 + a1x + a2x^2 + ..., e is
        represented as: 'n' where n = a0 + a1 * p + a2 * p^2 + ...

        """
        return str(int(e))

    def _element_poly_repr(GFq self, GFq_element e):
        """
        Returns a polynomial expression in base.gen() of self.
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

    def axpy(GFq self,GFq_element a, GFq_element b, GFq_element c):
        """
        r <-  c + a * b mod p
        """
        cdef int r

        r = self.objectptr.axpy(r, a.object, b.object, c.object)
        return make_GFq_element(self,r)

    def axmy(GFq self,GFq_element a, GFq_element b, GFq_element c):
        """
        r <- a * b - c mod p
        """

        cdef int r

        r = self.objectptr.axmy(r, a.object, b.object, c.object, )
        return make_GFq_element(self,r)

    def amxy(GFq self,GFq_element a, GFq_element b, GFq_element c):
        """
        r <- c - a * b mod p
        """
        cdef int r

        r = self.objectptr.amxy(r , a.object, b.object, c.object, )
        return make_GFq_element(self,r)

    def _add(GFq self, int r, int l):
        cdef int res
        return int(self.objectptr.add(res, r , l ))

    def _mul(GFq self, int r, int l):
        cdef int res
        return int(self.objectptr.mul(res, r , l ))

    def _div(GFq self, int r, int l):
        cdef int res
        return int(self.objectptr.div(res, r , l ))

    def _sub(GFq self, int r, int l):
        cdef int res
        return int(self.objectptr.sub(res, r , l ))


cdef class GFq_iterator:
    """
    Iterator over GFq elements of degree 1. We iterate over fields of
    higher degree using the VectorSpace iterator.
    """
    cdef int iterator
    cdef GFq _parent

    def __init__(self, GFq parent):
        self._parent = parent
        self.iterator = -1

    def __next__(self):
        """
        """

        self.iterator=self.iterator+1

        if self.iterator==self._parent.characteristic():
            self.iterator = -1
            raise StopIteration

        return make_GFq_element(self._parent,self._parent.int2log(self.iterator))

    def __repr__(self):
        return "Iterator over %s"%self._parent

cdef GFq_copy(GFq orig):
    cdef GFq copy
    copy = GFq(orig.characteristic()**orig.degree())
    delete(copy.objectptr)
    copy.objectptr = gfq_factorycopy(gfq_deref(orig.objectptr))
    return copy

cdef class GFq_element(SageObject):
    cdef int object
    cdef GFq _parent
    cdef object __multiplicative_order

    def __init__(GFq_element self, GFq parent ):
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

    def __dealloc__(GFq_element self):
        pass

    def __repr__(GFq_element self):
        return self._parent._element_repr(self)

    def parent(self):
        """
        Returns parent finite field.
        """
        return self._parent

    def is_zero(GFq_element self):
        """
        Returns True if self == k(0).
        """
        return bool(self._parent.objectptr.isZero(self.object))

    def is_one(GFq_element self):
        """
        Returns True if self == k(1)
        """
        return bool(self._parent.objectptr.isOne(self.object))

    def is_unit(GFq_element self):
        """
        Returns True if self is an element of the prime subfield.
        """
        return bool(self._parent.objectptr.isunit(self.object))


    def is_square(GFq_element self):
        """
        """
        #copied from finite_field_element.py
        K = self._parent
        if K.characteristic() == 2:
            return True
        n = K.order() - 1
        a = self**(n / 2)
        return bool(a == 1)

    def __add__(self, other):
        cdef int r

        if not isinstance(self,GFq_element):
            other,self = self,other

        if not isinstance(other,GFq_element):
            return bin_op(self,other,operator.add)
        else:
            r = (<GFq_element>self)._parent.objectptr.add(r, (<GFq_element>self).object , (<GFq_element>other).object )
            return make_GFq_element((<GFq_element>self)._parent,r)

    def __mul__(self, other):
        cdef int r

        if not isinstance(self,GFq_element):
            other,self = self,other

        if not isinstance(other,GFq_element):
            return bin_op(self,other,operator.mul)
        else:
            r = (<GFq_element>self)._parent.objectptr.mul(r, (<GFq_element>self).object, (<GFq_element>other).object)
            return make_GFq_element((<GFq_element>self)._parent,r)

    def __div__(self, other):

        cdef int r

        if not isinstance(self,GFq_element):
            other,self = self,other

        if not isinstance(other,GFq_element):
            return bin_op(self,other,operator.div)
        else:
            r = (<GFq_element>self)._parent.objectptr.div(r, (<GFq_element>self).object, (<GFq_element>other).object)
            return make_GFq_element((<GFq_element>self)._parent,r)

    def __sub__(self, other):

        cdef int r

        if not isinstance(self,GFq_element):
            other,self = self,other

        if not isinstance(other,GFq_element):
            return bin_op(self,other,operator.sub)
        else:
            r = (<GFq_element>self)._parent.objectptr.sub(r, (<GFq_element>self).object, (<GFq_element>other).object)
            return make_GFq_element((<GFq_element>self)._parent,r)

    def __neg__(GFq_element self):
        cdef int r

        r = self._parent.objectptr.neg(r, self.object)
        return make_GFq_element(self._parent,r)

    def __invert__(GFq_element self):
        cdef int r

        self._parent.objectptr.inv(r, self.object)
        return make_GFq_element(self._parent,r)


    def __pow__(GFq_element self, int exp, other):
        #There doesn't seem to exist a power function for GFq. So we
        #had to write one. It is pretty clumbsy (read: slow) right now

        cdef int power
        cdef int i
        cdef int epow2
        cdef GivaroGfq *field

        field = self._parent.objectptr

        exp = exp % (self._parent.order()-1)

        if field.isOne(self.object):
            return self

        if exp==0:
            return make_GFq_element(self._parent,field.one)

        power = field.one
        i = 0;
        epow2 = self.object;
        while (exp>>i) > 0:
            if (exp>>i) & 1:
                field.mulin(power,epow2)
            field.mulin(epow2,epow2)
            i = i + 1

        return make_GFq_element(self._parent,power)

    def add(GFq_element self,GFq_element other):
        cdef int r
        r = self._parent.objectptr.add(r, self.object , other.object )
        return make_GFq_element(self._parent,r)

    def mul(GFq_element self,GFq_element other):
        cdef int r
        r = self._parent.objectptr.mul(r, self.object , other.object )
        return make_GFq_element(self._parent,r)


    def div(GFq_element self,GFq_element  other):
        cdef int r
        r = self._parent.objectptr.div(r, self.object , other.object )
        return make_GFq_element(self._parent,r)

    def sub(GFq_element self,GFq_element other):
        cdef int r
        r = self._parent.objectptr.sub(r, self.object , other.object )
        return make_GFq_element(self._parent,r)

    def __cmp__(self, other):
        return cmp(int(self),int(other))

    def __richcmp__(self, right, int op):
        if op == 0:  #<
            try:
                return bool(int(self) < int(right))
            except:
                return False
        if op == 2: #==
            try:
                return bool(int(self) == int(right))
            except:
                return False
        if op == 4: #>
            try:
                return bool(int(self) > int(right))
            except:
                return False
        if op == 1: #<=
            try:
                return bool(int(self) <= int(right))
            except:
                return False
        if op == 3: #!=
            try:
                return bool(int(self) != int(right))
            except:
                return True
        if op == 5: #>=
            try:
                return bool(int(self) >= int(right))
            except:
                return False

    def __int__(GFq_element self):
        """
        Returns self coerced to an int. The integer returned is
        equivalent to the representation of self and not to log_repr.
        """
        return self._parent.log2int(self.object)


    def logint(GFq_element self):
        """
        Returns i where base.gen()^i=self
        """
        return int(self.object)

    def log(GFq_element self, a):
        #copied from finite_field_element.py
        q = (self.parent()).order() - 1
        return sage.rings.arith.discrete_log_generic(self, a, q)



    def int_repr(GFq_element self):
        return self._parent._element_int_repr(self)

    def log_repr(GFq_element self):
        return self._parent._element_log_repr(self)

    def poly_repr(GFq_element self):
        return self._parent._element_poly_repr(self)

    def polynomial(GFq_element self):
        quo = self._parent.log2int(self.object)
        b   = int(self._parent.characteristic())

        ret = []
        for i in range(self._parent.degree()):
            coeff = quo%b
            ret.append(coeff)
            quo = quo/b
        return self._parent.polynomial_ring()(ret)

    def _latex_(GFq_element self):
        if self._parent.degree()>1:
            return self.polynomial()._latex_()
        else:
            return str(self)

    def _sage_(GFq_element self, k=None):
        """
        Returns an element of k supposed to match this element.  No
        checks if k equals self.parent() are performed.

        INPUT:
            k -- SAGE finite field

        OUTPUT:
            k.gen()^(self.log_repr())

        """
        if k==None:
            k=self._parent._sage_()

        variable = k.gen()._pari_()

        quo = int(self)
        b   = self._parent.characteristic()

        ret = k._pari_one() - k._pari_one()
        i = 0
        while quo!=0:
            coeff = quo%b
            if coeff != 0:
                ret = coeff * variable ** i + ret
            quo = quo/b
            i = i+1
        return k(ret)

    def _pari_init_(GFq_element self):
        k=self._parent._sage_()

        variable = k.gen()._pari_()

        quo = int(self)
        b   = self._parent.characteristic()

        ret = k._pari_one() - k._pari_one()
        i = 0
        while quo!=0:
            coeff = quo%b
            if coeff != 0:
                ret = coeff * variable ** i + ret
            quo = quo/b
            i = i+1
        return ret

    def multiplicative_order(GFq_element self):
        """
        """
        # code copy'n'pasted from finite_field_element.py
        import sage.rings.arith
        from sage.ext.integer import Integer

        if self.__multiplicative_order!=None:
            return self.__multiplicative_order
        else:
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
            self.__multiplicative_order = order
            return order

    def copy(self):
        return make_GFq_element(self._parent,self.object)

    def _gap_init_(GFq_element self):
        #copied from finite_field_element.py
        F = self._parent
        if F.order() > 65536:
            raise TypeError, "order (=%s) must be at most 65536."%F.order()
        if self == 0:
            return '0*Z(%s)'%F.order()
        assert F.degree() > 1
        g = F.multiplicative_generator()
        n = g.log(self)
        return 'Z(%s)^%s'%(F.order(), n)

    def charpoly(GFq_element self):
        return self._sage_().charpoly()

    def norm(GFq_element self):
        return self._sage_().norm()

    def trace(GFq_element self):
        return self._sage_().trace()

    def __hash__(GFq_element self):
        # GF elements are hashed by hashing their string
        # representation but string representations are slow. So we
        # hash the log and the int representation which should provide
        # the same level of distinction.
        return hash((self.object,self._parent.log2int(self.object),"givaro"))

cdef make_GFq_element(GFq parent, int x):
    """
    """
    cdef GFq_element y
    y = GFq_element(parent)
    y.object = x
    return y

cdef gap_to_givaro(x, F):
    """
    INPUT:
        x -- gap finite field element
        F -- Givaro finite field
    OUTPUT:
        element of F

    EXAMPLES:
        sage: x = gap('Z(13)')
        sage: F = linbox.GFq(13)
        sage: F(x)
        2
        sage: F(gap('0*Z(13)'))
        0
        sage: F = linbox.GFq(13^2)
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
    if q == F.order():
        K = F
    else:
        K = GFq(q)
    if s.find(')^') == -1:
        e = 1
    else:
        e = int(s[i2+2:])

    if F.degree() == 1:
        g = int(sage.interfaces.gap.gap.eval('Int(Z(%s))'%q))
    else:
        g = K.multiplicative_generator()
    return F(K(g**e))
