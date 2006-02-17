"""
Number Fields
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
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

# There will be one running instance of GP for all
# number field calculations that use the interpreter.
from sage.interfaces.gp import Gp

import sage.interfaces.all
import sage.misc.preparser

_gp = None
def gp():
    global _gp
    if not _gp is None:
        return _gp
    else:
        _gp = Gp()
        return _gp

import operator

import weakref

from sage.misc.latex import latex

import sage.rings.arith as arith
import sage.rings.field as field
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.infinity as infinity
import sage.rings.rational as rational
import sage.rings.integer as integer
import sage.rings.polynomial_ring as polynomial_ring
import sage.rings.polynomial_element as polynomial_element
import sage.rings.ideal as ideal
import sage.rings.complex_field

import sage.groups.abelian


import number_field_element


from sage.libs.all import pari, gen

Q = rational_field.RationalField()
Z = integer_ring.IntegerRing()
R = polynomial_ring.PolynomialRing(Q)

_objsNumberField = {}
def NumberField(polynomial, name='a', check=False):
    global _objsNumberField
    key = (polynomial, name)
    if _objsNumberField.has_key(key):
        K = _objsNumberField[key]()
        if K != None:
            return K
    if polynomial.degree() == 2:
        K = NumberField_quadratic(polynomial, name, check)
    else:
        K = NumberField_generic(polynomial, name, check)

    _objsNumberField[key] = weakref.ref(K)
    return K

def QuadraticField(D, name='a', check=False):
    x = polynomial_ring.PolynomialRing(Q, 'x').gen()
    return NumberField(x**2 - D, name, check)

def is_QuadraticField(x):
    return isinstance(x, NumberField_quadratic)

def is_NumberField(x):
    return isinstance(x, NumberField_generic)

def CyclotomicField(n):
    return NumberField_cyclotomic(n)

def is_CyclotomicField(x):
    return isinstance(x, NumberField_cyclotomic)


class NumberField_generic(field.Field):
    """
    EXAMPLES:
        sage: K = NumberField(x^3 - 2, 'a'); K
        Number Field in a with defining polynomial x^3 - 2
        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, polynomial, name=None, check=True):
        if not isinstance(polynomial, polynomial_element.Polynomial):
            raise TypeError, "polynomial (=%s) must be a polynomial"%polynomial

        if check:
            polynomial = R(polynomial)

        if check and not polynomial.is_irreducible():
            raise ValueError, "defining polynomial (%s) must be irreducible"%polynomial

        if not polynomial.is_monic():
            raise NotImplementedError, "number fields for non-monic polynomials not yet implemented."

        self.assign_names(name)
        self.__polynomial = polynomial
        self.__degree = polynomial.degree()

    def __repr__(self):
        return "Number Field in %s with defining polynomial %s"%(
                   self.variable_name(), self.polynomial())

    def _latex_(self):
        return "%s[%s]/(%s)"%(latex(Q), self.variable_name(), self.polynomial()._latex_(self.variable_name()))

    def __call__(self, x):
        """
        Coerce x into this number field.
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            if x.parent() == self:
                return x
            # todo: more general coercision if embedding have been asserted

        if not isinstance(x, (int, long, rational.Rational, integer.Integer, gen, list)):
            raise TypeError, "Cannot coerce %s into %s"%(x,self)

        return number_field_element.NumberFieldElement(self, x)

    def _coerce_(self, x):
        if isinstance(x, number_field_element.NumberFieldElement):
            if x.parent() == self:
                return x
        if isinstance(x, (rational.Rational, integer.Integer, int, long)):
            return number_field_element.NumberFieldElement(self, x)
        raise TypeError

    def __cmp__(self, other):
        if not isinstance(other, NumberField_generic):
            return -1
        if self.variable_name() != other.variable_name():
            return -1
        return self.__polynomial.__cmp__(other.__polynomial)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            # We need that elements of the base ring of the polynomial
            # ring map canonically into codomain.
            codomain._coerce_(rational.Rational(1))
            f = self.defining_polynomial()
            return codomain(f(im_gens[0])) == 0
        except TypeError, ValueError:
            return False

    def base_ring(self):
        return rational_field.Q

    def pari_polynomial(self):
        """
        PARI polynomial corresponding to polynomial that defines
        this field.
        """
        try:
            return self.__pari_polynomial
        except AttributeError:
            self.__pari_polynomial = self.polynomial()._pari_()
            return self.__pari_polynomial

    def pari_nf(self):
        """
        PARI number field corresponding to this field.
        """
        try:
            return self.__pari_nf
        except AttributeError:
            f = self.pari_polynomial()
            self.__pari_nf = f.nfinit()
            return self.__pari_nf

    def pari_bnf(self):
        """
        PARI big number field corresponding to this field.
        """
        try:
            return self.__pari_bnf
        except AttributeError:
            f = self.pari_polynomial()
            self.__pari_bnf = f.bnfinit()
            return self.__pari_bnf

    def characteristic(self):
        return 0

    def class_group(self):
        """
        WARNING: Assume GRH, etc. !!
          TODO: Change to use bnf_certify, unless user requests not to.
        """
        try:
            return self.__class_group
        except AttributeError:
            k = self.pari_bnf()
            s = str(k[7][0])  # it's the [8][1] entry in pari, but the python interface is 0 based.
            s = s.replace(";",",")
            s = eval(s)
            self.__class_group = sage.groups.abelian.AbelianGroup(s[1])
        return self.__class_group

    def class_number(self):
        return self.class_group().order()

    def composite_fields(self, other):
        """
        List of all possible composite fields formed from self and other.
        """
        if not isinstance(other, NumberField_generic):
            raise TypeError, "other must be a number field."
        f = self.pari_polynomial()
        g = other.pari_polynomial()
        C = f.polcompositum(g)
        R = self.polynomial().parent()
        C = [R(h) for h in C]
        return [NumberField(h) for h in C]

    def degree(self):
        return self.__degree

    def discriminant(self, v=None):
        """
        Returns the discriminant of the ring of integers of the number field,
        or if v is specified, the determinant of the trace pairing
        on the elements of the list v.

        INPUT:
            v (optional) -- list of element of this number field
        OUTPUT:
            Integer if v is omitted, and Rational otherwise.
        """
        if v == None:
            try:
                return self.__disc
            except AttributeError:
                self.__disc = Z(str(self.pari_nf()[2]))
                return self.__disc
        else:
            return Q(self.trace_pairing(v).det())

    disc = discriminant

    def factor_integer(self, n):
        """
        Ideal factorization of the principal ideal of the ring
        of integers generated by n.
        """
        F = list(self.pari_nf().idealfactor(n))
        P, exps = F[0], F[1]
        A = []
        for i, p in enumerate(P):
            B = [Z(x) for x in p[1]]
            I = FractionalIdeal(self, T="generators", data=(p[0], B))
            I._ramification = p[2]
            I._residue_class_degree = p[3]
            A.append((I,exps[i]))

        return A

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "Only one generator."
        try:
            return self.__gen
        except AttributeError:
            if self.__polynomial != None:
                X = self.__polynomial.parent().gen()
            else:
                X = PolynomialRing(rational_field.RationalField()).gen()
            self.__gen = number_field_element.NumberFieldElement(self, X)
            return self.__gen

    def is_field(self):
        return True

    def galois_group(self, pari_group = False, use_kash=False):
        r"""
        Return the Galois group of the Galois closure of this number
        field as an abstract group.

        For more (important!) documentation, so the documentation
        for Galois groups of polynomials over $\Q$, e.g., by
        typing \code{K.polynomial().galois_group?}, where $K$
        is a number field.

        EXAMPLES:
            sage: NumberField(x^3-2).galois_group(pari_group=True)
            PARI group [6, -1, 2, "S3"] of degree 3

            sage: NumberField(x-1).galois_group()    # optional database_gap package
            Transitive group number 1 of degree 1
            sage: NumberField(x^2+2).galois_group()  # optional database_gap package
            Transitive group number 1 of degree 2
            sage: NumberField(x^3-2).galois_group()  # optional database_gap package
            Transitive group number 2 of degree 3
        """
        return self.polynomial().galois_group(pari_group = pari_group, use_kash = use_kash)


    def integral_basis(self):
        """
        Return a list of elements of this number field that are a basis
        for the full ring of integers.

        EXAMPLES:
            sage: x = PolynomialRing(QQ).gen()
            sage: K = NumberField(x^5+10*x+1, 'a')
            sage: K.integral_basis()
            [1, a, a^2, a^3, a^4]

        Next we compute the ring of integers of a cubic field in which 2
        is an "essential discriminant divisor", so the ring of integers
        is not generated by a single element.
            sage: K = NumberField(x^3 + x^2 - 2*x + 8, 'a')
            sage: K.integral_basis()
            [1, a, 1/2*a^2 + 1/2*a]
        """
        try:
            return self.__integral_basis
        except AttributeError:
            f = self.pari_polynomial()
            B = f.nfbasis()
            R = self.polynomial().parent()
            self.__integral_basis = [self(R(g).list()) for g in B]
        return self.__integral_basis

    def ngens(self):
        return 1

    def order(self):
        return infinity.infinity

    def order_table(self):
        return []

    def polynomial(self):
        return self.__polynomial

    def defining_polynomial(self):
        return self.__polynomial

    def polynomial_ring(self):
        return R

    def polynomial_quotient_ring(self):
        """
        Return the polynomial quotient ring isomorphic to this number field.

        EXAMPLES:
            sage: R = PolynomialRing(RationalField(), 'x'); x = R.gen()
            sage: K = NumberField(x^3 + 2*x - 5, 'alpha')
            sage: K.polynomial_quotient_ring()
            Univariate Quotient Polynomial Ring in alpha over Rational Field with modulus x^3 + 2*x - 5
        """
        return self.polynomial_ring().quotient(self.polynomial(), self.variable_name())

    def trace_pairing(self, v):
        """
        Returns the trace pairing on the elements of the list v.
        """
        import sage.matrix.matrix_space
        A = sage.matrix.matrix_space.MatrixSpace(Q,len(v))(0)
        for i in range(len(v)):
            for j in range(i,len(v)):
                t = (v[i]*v[j]).trace()
                A[i,j] = t
                A[j,i] = t
        return A

    def units(self):
        try:
            return self.__units
        except AttributeError:
            B = self.pari_bnf().bnfunit()
            R = self.polynomial().parent()
            self.__units = [self(R(g)) for g in B]
            return self.__units

    def zeta(self, n=2):
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        else:
            raise NotImplementedError


class NumberField_cyclotomic(NumberField_generic):
    """
    Create a cyclotomic extension of the rational field.

    The command CyclotomicField(n) creates the n-th cyclotomic
    field, got by adjoing an n-th root of unity to the rational
    field.

    EXAMPLES:
        sage: CyclotomicField(3)
        Cyclotomic Field of order 3 and degree 2
        sage: CyclotomicField(18)
        Cyclotomic Field of order 18 and degree 6
        sage: z = CyclotomicField(6).gen(); z
        zeta_6
        sage: z^3
        -1
        sage: (1+z)^3
        6*zeta_6 - 3

        sage: K = CyclotomicField(197)
        sage: loads(K.dumps()) == K
        True
        sage: loads((z^2).dumps()) == z^2
        True
    """
    def __init__(self, n):
        f = R.cyclotomic_polynomial(n)
        NumberField_generic.__init__(self, f, name="zeta_%s"%n, check=False)
        n = integer.Integer(n)
        zeta = self.gen()
        zeta._set_multiplicative_order(n)
        self.__zeta_order = n

    def __repr__(self):
        return "Cyclotomic Field of order %s and degree %s"%(
                self.zeta_order(), self.degree())

    def _latex_(self):
        return "%s(\\zeta_{%s})"%(latex(Q), self.__zeta_order)

    def __call__(self, x):
        """
        Create an element of this cyclotomic field from $x$.

        EXAMPLES:
        The following example illustrates coercion from the cyclotomic
        field Q(zeta_42) to the cyclotomic field Q(zeta_6), in a case
        where such coercion is defined:

            sage: k42 = CyclotomicField(42)
            sage: k6 = CyclotomicField(6)
            sage: a = k42.gen(0)
            sage: b = a^7
            sage: b
            zeta_42^7
            sage: k6(b)
            zeta_6
            sage: b^2
            zeta_42^7 - 1
            sage: k6(b^2)
            zeta_6 - 1

        Coercion of GAP cyclotomic elements is also fully supported.


        """
        if isinstance(x, number_field_element.NumberFieldElement) and \
                isinstance(x.parent(), NumberField_cyclotomic):
            K = x.parent()
            n = K.zeta_order()
            m = self.zeta_order()
            if m % n == 0:   # easy case
                e = m/n
                f = x.polynomial()
                X = f.parent().gen()
                g = f(X**e)
            else:
                n = x.multiplicative_order()
                if m % n == 0:
                    # Harder case.  E.g., x = (zeta_42)^7 and
                    # self.__zeta = zeta_6, so it is possible to
                    # coerce x in, but not zeta_42 in.
                    # Algorithm:
                    #    1. Compute self.__zeta as an element
                    #       of K = parent of x.  Call this y.
                    #    2. Write x as a power r of y.
                    #       TODO: we do step two STUPIDLY.
                    #    3. Return self.__zeta to the power r.
                    y = K(self.zeta())
                    z = y
                    for r in xrange(y.multiplicative_order()):
                        if z == x:
                            return self.zeta()**(r+1)
                        z *= y
                raise TypeError, "Cannot coerce %s into %s"%(x,self)
            return number_field_element.NumberFieldElement(self, g)

        elif sage.interfaces.all.is_GapElement(x):
            s = str(x)
            i = s.find('E(')
            if i == -1:
                return self(rational.Rational(s))
            j = i + s[i:].find(')')
            n = int(s[i+2:j])
            if n == self.zeta_order():
                K = self
            else:
                K = CyclotomicField(n)
            zeta = K.gen()
            s = s.replace('E(%s)'%n,'zeta')
            s = sage.misc.all.sage_eval(s, {'zeta':K.gen()})
            if K is self:
                return s
            else:
                return self(s)

        else:
            return NumberField_generic.__call__(self, x)

    def complex_embedding(self, prec=53):
        r"""
        Return the embedding of this cyclotomic field into the
        approximate complex field with precision prec obtained by
        sending the generator $\zeta$ of self to exp(2*pi*i/n), where
        $n$ is the multiplicative order of $\zeta$.

        EXAMPLES:
            sage: C = CyclotomicField(4)
            sage: C.complex_embedding()
            Ring morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Complex Field with 53 bits of precision
              Defn: zeta_4 |--> 0.000000000000000061232339957367660 + 1.0000000000000000*I

        Note in the example above that the way zeta is computed (using
        sin and cosine in MPFR) means that only the prec bits of the
        number after the decimal point are valid.

            sage: K = CyclotomicField(3)
            sage: phi = K.complex_embedding (10)
            sage: phi(K.0)
            -0.49951 + 0.86621*I
            sage: phi(K.0^3)
            1.0000
            sage: phi(K.0^3 - 1)
            0
            sage: phi(K.0^3 + 7)
            8.0000
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        return self.hom([CC.zeta(self.zeta_order())], check=False)

    def integral_basis(self):
        """
        Return a list of elements of this number field that are a basis
        for the full ring of integers.
        """
        try:
            return self.__integral_basis
        except AttributeError:
            z = self.gen()
            a = self(1)
            B = []
            for n in xrange(self.degree()):
                B.append(a)
                a *= z
            self.__integral_basis = B
        return self.__integral_basis


    def zeta_order(self):
        return self.__zeta_order

    def multiplicative_order_table(self):
        try:
            return self.__multiplicative_order_table
        except AttributeError:
            t = {}
            x = self(1)
            n = self.zeta_order()
            m = 0
            zeta = self.zeta()
            # todo: this desperately needs to be optimized!!!
            for i in range(n):
                t[x.polynomial()] = n//arith.GCD(m,n)   # multiplicative_order of (zeta_n)**m
                x *= zeta
                m += 1
            self.__multiplicative_order_table = t
            return t

    def zeta(self, n=None):
        """
        Returns an element of multiplicative order $n$ in this this
        number field, if there is one.  Raises a ValueError if there
        is not.

        INPUT:
            n -- integer (default: None, returns element of maximal order)

        OUTPUT:
            root of unity

        EXAMPLES:
            sage: k = CyclotomicField(7)
            sage: k.zeta()
            zeta_7
            sage: k.zeta().multiplicative_order()
            7
            sage: k = CyclotomicField(49)
            sage: k.zeta().multiplicative_order()
            49
            sage: k.zeta(7).multiplicative_order()
            7
            sage: k.zeta()
            zeta_49
            sage: k.zeta(7)
            zeta_49^7
        """
        if n is None:
            return self.gen()
        else:
            n = integer.Integer(n)
            z = self.gen()
            m = z.multiplicative_order()
            if m % n != 0:
                raise ValueError, "No %sth root of unity in self"%n
            return z**(m//n)

class NumberField_quadratic(NumberField_generic):
    """
    Create a quadratic extension of the rational field.

    The command QuadraticExtension(a) creates the field Q(sqrt(a)).

    EXAMPLES:
        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3
        sage: QuadraticField(-4)
        Number Field in a with defining polynomial x^2 + 4
    """
    def __init__(self, polynomial, name=None, check=True):
        NumberField_generic.__init__(self, polynomial, name=name, check=check)


    def class_number(self, proof = True):
        """
        Return the size of the class group of self.

        If proof = False (not the default) and the discriminant of the
        field is negative, then the following warning from the PARI
        manual applies: IMPORTANT WARNING: For D<0, this function may
        give incorrect results when the class group has a low exponent
        (has many cyclic factors), because implementing Shank's method
        in full generality slows it down immensely.
        """
        try:
            return self.__class_number
        except AttributeError:
            D = self.discriminant()
            if D < 0 and proof:
                self.__class_number = pari("qfbclassno(%s,1)"%D).python()
            else:
                self.__class_number = pari("qfbclassno(%s)"%D).python()
            return self.__class_number

    def hilbert_class_polynomial(self):
        r"""
        Returns a polynomial over $\Q$ whose roots generate the
        Hilbert class field of this quadratic field.

        \note{Computed using PARI via Schertz's method.  This implementation
        is amazingly fast.}

        EXAMPLES:
            sage: K = NumberField(x^2 + 23)
            sage: K.hilbert_class_polynomial()
            x^3 + x^2 - 1

            sage: K = NumberField(x^2 + 431)
            sage: K.hilbert_class_polynomial()
            x^21 + x^20 - 13*x^19 - 50*x^18 + 592*x^17 - 2403*x^16 + 5969*x^15 - 10327*x^14 + 13253*x^13 - 12977*x^12 + 9066*x^11 - 2248*x^10 - 5523*x^9 + 11541*x^8 - 13570*x^7 + 11315*x^6 - 6750*x^5 + 2688*x^4 - 577*x^3 + 9*x^2 + 15*x + 1
        """
        f = pari('quadhilbert(%s))'%self.discriminant())
        g = R(list(reversed(f.list())))
        return g

    def hilbert_class_field(self):
        r"""
        Returns the Hilbert class field of this quadratic
        field as an absolute extension of $\Q$.  For a polynomial
        that defines a relative extension see the
        \code{hilbert_class_polynomial} command.

        \note{Computed using PARI via Schertz's method.  This implementation
        is amazingly fast.}

        EXAMPLES:
            sage: K = NumberField(x^2 + 23)
            sage: K.hilbert_class_polynomial()
            x^3 + x^2 - 1
            sage: K.hilbert_class_field()
            Number Field in a with defining polynomial x^6 + 2*x^5 + 70*x^4 + 90*x^3 + 1631*x^2 + 1196*x + 12743
        """
        f = self.hilbert_class_polynomial()
        C = self.composite_fields(NumberField(f))
        assert len(C) == 1
        return C[0]

def is_fundamental_discriminant(D):
    d = D % 4
    if not (d in [0,1]):
        return False
    return D != 1 and  D != 0 and \
           (arith.is_squarefree(D) or \
            (d == 0 and (D/4)%4 in [2,3] and arith.is_squarefree(D/4)))


# TODO: To be removed or related to contents of sage/rings/ideal.
class FractionalIdeal:
    def __init__(self, number_field, T, data):
        self._number_field = number_field
        self._T = T
        if T == "generators":
            self._generators = data
            self._characteristic = data[0]
        else:
            raise RuntimeError, "Defining data type %s not known."%T

    def __repr__(self):
        return "Fractional ideal of %s"%self.number_field()

    def characteristic(self):
        return self._characteristic

    def gens(self):
        """
        Two elements that generate the ideal.  The first is an integer and
        the second is a Z-linear combination of the chosen basis for O_K.
        OUTPUT:
            tuple -- (p, v)
        """
        return self._generators

    def hermite_basis(self):
        return self._hermite_basis

    def idele(self):
        return self._idele

    def integral_basis(self):
        """
        Integral basis (over Z), which need not be in Hermite normal form.
        """
        return self._integral_basis

    def is_principal(self):
        raise NotImplementedError

    def number_field(self):
        return self._number_field

    def ramification(self):
        return self._ramification

    def residue_class_degree(self):
        return self._residue_class_degree

    def residue_class_field(self):
        """
        """
        raise NotImplementedError


