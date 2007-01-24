"""
Number Fields

AUTHORS:
   -- William Stein (2004, 2005): initial version
   -- Steven Sivek (2006-05-12): added support for relative extensions
"""

#*****************************************************************************
#       Copyright (C) 2004, 2005, 2006 William Stein <wstein@gmail.com>
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

import sage.interfaces.gap
import sage.misc.preparser
import sage.rings.arith
import sage.rings.complex_field
import sage.rings.ring

import sage.structure.parent_gens

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
import sage.groups.abelian_gps.abelian_group

from sage.structure.parent_gens import ParentWithGens
import number_field_element
from number_field_ideal import convert_from_zk_basis

from sage.libs.all import pari, pari_gen

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()

_nf_cache = {}
def NumberField(polynomial, name=None, check=True, names=None):
    r"""
    Return {\em the} number field defined by the given irreducible
    polynomial and with variable with the given name.  If check is
    True (the default), also verify that the defining polynomial is
    irreducible and over Q.

    INPUT:
        polynomial -- a polynomial over Q  (for now)
        name -- a string (default: 'a'), the name of the generator
        check -- bool (default: True); do type checking and
                 irreducibility checking.

    EXAMPLES:
        sage: z = QQ['z'].0
        sage: K = NumberField(z^2 - 2,'s'); K
        Number Field in s with defining polynomial z^2 - 2
        sage: s = K.0; s
        s
        sage: s*s
        2
        sage: s^2
        2

    EXAMPLES: Constructing a relative number field
        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(x^2 - 2)
        sage: R.<t> = K['t']
        sage: L = K.extension(t^3+t+a, 'b'); L
        Extension by t^3 + t + a of the Number Field in a with defining polynomial x^2 - 2
        sage: L.absolute_field()
        Number Field in b with defining polynomial x^6 + 2*x^4 + x^2 - 2
        sage: b = L.gen()
        sage: a*b
        -b^4 - b^2
        sage: L.lift_to_base(-3*b^3 - 3*b + 1)
        3*a + 1

    Number fields are globally unique.
        sage: K.<a>= NumberField(x^3-5)
        sage: a^3
        5
        sage: L.<a>= NumberField(x^3-5)
        sage: K is L
        True
    """
    if name is None and names is None:
        raise TypeError, "You must specify the name of the generator."
    if not names is None:
        name = names

    name = sage.structure.parent_gens.normalize_names(1, name)
    key = (polynomial, name)
    if _nf_cache.has_key(key):
        K = _nf_cache[key]()
        if not K is None: return K

    R = polynomial.base_ring()
    if R == ZZ:
        polynomial = QQ['x'](polynomial)
    elif isinstance(R, NumberField_generic):
        S = R.extension(polynomial, name)
        _nf_cache[key] = weakref.ref(S)
        return S

    if polynomial.degree() == 2:
        K = NumberField_quadratic(polynomial, name, check)
    else:
        K = NumberField_generic(polynomial, name, check)

    _nf_cache[key] = weakref.ref(K)
    return K

def QuadraticField(D, names, check=False):
    x = polynomial_ring.PolynomialRing(QQ, 'x').gen()
    return NumberField(x**2 - D, names, check)

def is_QuadraticField(x):
    return isinstance(x, NumberField_quadratic)

def is_NumberField(x):
    return isinstance(x, NumberField_generic)

def is_NumberFieldExtension(x):
    return isinstance(x, NumberField_extension)

_cyclo_cache = {}
def CyclotomicField(n, names=None):
    if names is None:
        names = "zeta%s"%n
    names = sage.structure.parent_gens.normalize_names(1, names)
    key = (n, names)
    if _cyclo_cache.has_key(key):
        K = _cyclo_cache[key]()
        if not K is None: return K
    K = NumberField_cyclotomic(n, names)
    _cyclo_cache[key] = weakref.ref(K)
    return K

def is_CyclotomicField(x):
    return isinstance(x, NumberField_cyclotomic)


class NumberField_generic(field.Field):
    """
    EXAMPLES:
        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(x^3 - 2); K
        Number Field in a with defining polynomial x^3 - 2
        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, polynomial, name,
                 latex_name=None, check=True):
        ParentWithGens.__init__(self, QQ, name)
        if not isinstance(polynomial, polynomial_element.Polynomial):
            raise TypeError, "polynomial (=%s) must be a polynomial"%polynomial

        if check:
            if not polynomial.is_irreducible():
                raise ValueError, "defining polynomial (%s) must be irreducible"%polynomial
            if not polynomial.parent().base_ring() == QQ:
                raise TypeError, "polynomial must be defined over rational field"
        if not polynomial.is_monic():
            raise NotImplementedError, "number fields for non-monic polynomials not yet implemented."

        self._assign_names(name)
        if latex_name is None:
            self.__latex_variable_name = self.variable_name()
        else:
            self.__latex_variable_name = latex_name
        self.__polynomial = polynomial
        self.__pari_bnf_certified = False
        self.__absolute_field = self

    def complex_embeddings(self, prec=53):
        r"""
        Return all homomorphisms of this ring into the approximate
        complex field with precision prec.

        EXAMPLES:
            sage: f = x^5 + x + 17
            sage: k.<a> = NumberField(f)
            sage: v = k.complex_embeddings()
            sage: [phi(k.0^2) for phi in v]
            [2.97572074037667, 0.921039066973047 - 3.07553311884577*I, 0.921039066973047 + 3.07553311884577*I, -2.40889943716138 + 1.90254105303505*I, -2.40889943716138 - 1.90254105303505*I]
        """
        try:
            return self.__complex_embeddings[prec]
        except AttributeError:
            self.__complex_embeddings = {}
        except KeyError:
            pass
        CC = sage.rings.complex_field.ComplexField(prec)
        f = self.defining_polynomial().base_extend(CC)
        v = f.roots()
        e = [self.hom([a], check=False) for a in v]
        self.__complex_embeddings[prec] = e
        return e

    def latex_variable_name(self, name=None):
        if name is None:
            return self.__latex_variable_name
        else:
            self.__latex_variable_name = name

    def __repr__(self):
        return "Number Field in %s with defining polynomial %s"%(
                   self.variable_name(), self.polynomial())

    def _latex_(self):
        return "%s[%s]/(%s)"%(latex(QQ), self.variable_name(),
                              self.polynomial()._latex_(self.variable_name()))

    def __call__(self, x):
        """
        Coerce x into this number field.
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return number_field_element.NumberFieldElement(self, x.polynomial())
            return self._coerce_from_other_number_field(x)
        return self._coerce_non_number_field_element_in(x)

    def _coerce_from_other_number_field(self, x):
        f = x.polynomial()
        if f.degree() <= 0:
            return number_field_element.NumberFieldElement(self, f[0])
        # todo: more general coercion if embedding have been asserted
        raise TypeError, "Cannot coerce %s into %s"%(x,self)

    def _coerce_non_number_field_element_in(self, x):
        if isinstance(x, (int, long, rational.Rational,
                              integer.Integer, pari_gen,
                              polynomial_element.Polynomial,
                              list)):
            return number_field_element.NumberFieldElement(self, x)
        raise TypeError, "Cannot coerce %s into %s"%(x,self)

    def _coerce_impl(self, x):
        if isinstance(x, (rational.Rational, integer.Integer, int, long)):
            return number_field_element.NumberFieldElement(self, x)
        raise TypeError

    def category(self):
        from sage.categories.all import NumberFields
        return NumberFields()


    def category(self):
        from sage.categories.all import NumberFields
        return NumberFields()

    def __cmp__(self, other):
        if not isinstance(other, NumberField_generic):
            return -1
        if self.variable_name() != other.variable_name():
            return -1
        return self.__polynomial.__cmp__(other.__polynomial)

    def _ideal_class_(self):
        return sage.rings.number_field.number_field_ideal.NumberFieldIdeal

    def ideal(self, gens):
        r"""
        Return the ideal in $\mathcal{O}_K$ generated by gens.  This
        overrides the \code{sage.rings.ring.Field} method to use the
        \code{sage.rings.ring.Ring} one instead, since we're not really
        concerned with ideals in a field but in its ring of integers.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3-2)
            sage: K.ideal([a])
            Fractional ideal (a) of Number Field in a with defining polynomial x^3 - 2
        """
        return sage.rings.ring.Ring.ideal(self, gens)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            # We need that elements of the base ring of the polynomial
            # ring map canonically into codomain.
            codomain._coerce_(rational.Rational(1))
            f = self.defining_polynomial()
            return codomain(f(im_gens[0])) == 0
        except TypeError, ValueError:
            return False

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

    def pari_bnf(self, certify=False):
        """
        PARI big number field corresponding to this field.
        """
        try:
            if certify:
                self.pari_bnf_certify()
            return self.__pari_bnf
        except AttributeError:
            f = self.pari_polynomial()
            self.__pari_bnf = f.bnfinit()
            if certify:
                self.pari_bnf_certify()
            return self.__pari_bnf

    def pari_bnf_certify(self):
        """
        Run the PARI bnfcertify function to ensure the correctness of answers.
        """
        if not self.__pari_bnf_certified:
            if self.pari_bnf(certify=False).bnfcertify() != 1:
                raise ValueError, "The result is not correct according to bnfcertify"
            self.__pari_bnf_certified = True

    def characteristic(self):
        return 0

    def class_group(self, certify=True):
        r"""
        Return the class group of this field.
        """
        try:
            return self.__class_group
        except AttributeError:
            k = self.pari_bnf(certify)
            s = str(k.getattr('clgp'))
            s = s.replace(";",",")
            s = eval(s)
            self.__class_group = \
               sage.groups.abelian_gps.abelian_group.AbelianGroup(s[1])
        return self.__class_group

    def class_number(self, certify=True):
        return self.class_group(certify).order()

    def composite_fields(self, other, names):
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
        return [NumberField(h, names) for h in C]

    def degree(self):
        return self.polynomial().degree()

    def different(self):
        """
        Compute the different ideal of this number field.
        """
        try:
            return self.__different
        except AttributeError:
            diff = self.pari_nf().getattr('diff')
            zk_basis = self.pari_nf().getattr('zk')
            basis_elts = zk_basis * diff
            R = self.polynomial().parent()
            self.__different = self.ideal([ self(R(x)) for x in basis_elts ])
            return self.__different

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
                self.__disc = ZZ(str(self.pari_nf()[2]))
                return self.__disc
        else:
            return Q(self.trace_pairing(v).det())

    disc = discriminant

    def elements_of_norm(self, n, certify=True):
        r"""
        Return a list of solutions modulo units of positive norm to
        $Norm(a) = n$, where a can be any integer in this number field.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2+1)
            sage: K.elements_of_norm(3)
            []
            sage: K.elements_of_norm(50)
            [7*a - 1, -5*a + 5, a - 7]           # 32-bit
            [7*a - 1, -5*a + 5, -7*a - 1]        # 64-bit
        """
        B = self.pari_bnf(certify).bnfisintnorm(n)
        R = self.polynomial().parent()
        return [self(QQ['x'](R(g))) for g in B]

    def extension(self, poly, name=None, names=None):
        """
        Return the relative extension of this field by a given polynomial.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: t = K['x'].gen()
            sage: L.<b> = K.extension(t^2 + a); L
            Extension by x^2 + a of the Number Field in a with defining polynomial x^3 - 2
        """
        if not names is None: name = names
        if name is None:
            raise TypeError, "the variable name must be specified."
        return NumberField_extension(self, poly, name)

    def factor_integer(self, n):
        r"""
        Ideal factorization of the principal ideal of the ring
        of integers generated by $n$.

	EXAMPLE:
        Here we show how to factor gaussian integers.
        First we form a number field defined by $x^2 + 1$:

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<I> = NumberField(x^2 + 1); K
            Number Field in I with defining polynomial x^2 + 1

        Here are the factors:

	    sage: fi, fj = K.factor_integer(13);fi,fj
            ((Fractional ideal (3*I - 2) of Number Field in I with defining polynomial x^2 + 1, 1),
            (Fractional ideal (-3*I - 2) of Number Field in I with defining polynomial x^2 + 1, 1))

        Now we extract the reduced form of the generators:

	    sage: zi = fi[0].gens_reduced()[0]; zi
            3*I - 2
	    sage: zj = fj[0].gens_reduced()[0]; zj
            -3*I - 2

        We recover the integer that was factor in $\Z[i]$

	    sage: zi*zj
            13

	    AUTHOR:
                -- Alex Clemesha (2006-05-20): examples

        """
        return self.ideal(n).factor()

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
            sage: R.<x> = PolynomialRing(QQ)

            sage: NumberField(x^3-2, 'a').galois_group(pari_group=True)
            PARI group [6, -1, 2, "S3"] of degree 3

            sage: NumberField(x-1, 'a').galois_group()    # optional database_gap package
            Transitive group number 1 of degree 1
            sage: NumberField(x^2+2, 'a').galois_group()  # optional database_gap package
            Transitive group number 1 of degree 2
            sage: NumberField(x^3-2, 'a').galois_group()  # optional database_gap package
            Transitive group number 2 of degree 3
        """
        return self.polynomial().galois_group(pari_group = pari_group, use_kash = use_kash)


    def integral_basis(self):
        """
        Return a list of elements of this number field that are a basis
        for the full ring of integers.

        EXAMPLES:
            sage: x = polygen(QQ,'x')
            sage: K.<a> = NumberField(x^5+10*x+1)
            sage: K.integral_basis()
            [1, a, a^2, a^3, a^4]

        Next we compute the ring of integers of a cubic field in which 2
        is an "essential discriminant divisor", so the ring of integers
        is not generated by a single element.
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
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

    def narrow_class_group(self, certify = True):
        r"""
        Return the narrow class group of this field.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: NumberField(x^3+x+9, 'a').narrow_class_group()
            Multiplicative Abelian Group isomorphic to C2
        """
        try:
            return self.__narrow_class_group
        except AttributeError:
            k = self.pari_bnf(certify)
            s = str(k.bnfnarrow())
            s = s.replace(";",",")
            s = eval(s)
            self.__narrow_class_group = sage.groups.abelian_gps.abelian_group.AbelianGroup(s[1])
        return self.__narrow_class_group

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
        return self.polynomial().parent()

    def polynomial_quotient_ring(self):
        """
        Return the polynomial quotient ring isomorphic to this number field.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K = NumberField(x^3 + 2*x - 5, 'alpha')
            sage: K.polynomial_quotient_ring()
            Univariate Quotient Polynomial Ring in alpha over Rational Field with modulus x^3 + 2*x - 5
        """
        return self.polynomial_ring().quotient(self.polynomial(), self.variable_name())

    def regulator(self, certify=True):
        """
        Return the regulator of this number field.

        Note that PARI computes the regulator to higher precision than
        the SAGE default.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: NumberField(x^2-2, 'a').regulator()
            0.88137358701954305
            sage: NumberField(x^4+x^3+x^2+x+1, 'a').regulator()
            0.96242365011920694
        """
        try:
            return self.__regulator
        except AttributeError:
            k = self.pari_bnf(certify)
            s = str(k.getattr('reg'))
            self.__regulator = eval(s)
        return self.__regulator

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: NumberField(x^2+1, 'a').signature()
            (0, 1)
            sage: NumberField(x^3-2, 'a').signature()
            (1, 1)
            sage: CyclotomicField(7).signature()
            (0, 3)
        """
        r1, r2 = self.pari_nf().getattr('sign')
        return (ZZ(r1), ZZ(r2))

    def trace_pairing(self, v):
        """
        Return the matrix of the trace pairing on the elements of the
        list $v$.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<zeta3> = NumberField(x^2 + 3)
            sage: K.trace_pairing([1,zeta3])
            [ 2  0]
            [ 0 -6]
        """
        import sage.matrix.matrix_space
        A = sage.matrix.matrix_space.MatrixSpace(self.base_ring(), len(v))(0)
        for i in range(len(v)):
            for j in range(i,len(v)):
                t = (self(v[i]*v[j])).trace()
                A[i,j] = t
                A[j,i] = t
        return A

    def units(self, certify = True):
        """
        Return generators for the unit group modulo torsion.

        ALGORITHM: Uses PARI's bnfunit command.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: A = x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3
            sage: K = NumberField(A, 'a')
            sage: K.units()
            [8/275*a^3 - 12/55*a^2 + 15/11*a - 2]
        """
        try:
            return self.__units
        except AttributeError:
            B = self.pari_bnf(certify).bnfunit()
            R = self.polynomial().parent()
            self.__units = [self(R(g)) for g in B]
            return self.__units


    def zeta(self, n=2, all=False):
        """
        Return an n-th root of unity in this field.  If all is True,
        return all of them.

        INPUT:
            n -- positive integer
            all -- bool, default: False.  If True, return a list
                   of all n-th roots of 1)

        If there are no n-th roots of unity in self, this function
        raises a ValueError exception.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^2 + 3, 'zeta3')
            sage: K.zeta(1)
            1
            sage: K.zeta(2)
            -1
            sage: K.zeta(2, all=True)
            [-1]
            sage: K.zeta(3)
            1/2*zeta3 - 1/2
            sage: K.zeta(3, all=True)
            [1/2*zeta3 - 1/2, -1/2*zeta3 - 1/2]
            sage: K.zeta(4)
            Traceback (most recent call last):
            ...
            ValueError: There are no 4-th roots of unity self.
        """
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n (=%s) must be positive"%n
        if n == 1:
            if all:
                return [self(1)]
            else:
                return self(1)
        elif n == 2:
            if all:
                return [self(-1)]
            else:
                return self(-1)
        else:
            field = self.__absolute_field
            f = field.polynomial_ring().cyclotomic_polynomial(n)
            F = polynomial_ring.PolynomialRing(field, 'x')(f)
            R = F.roots()
            if len(R) == 0:
                raise ValueError, "There are no %s-th roots of unity self."%n
            if all:
                return [r[0] for r in R]
            else:
                return R[0][0]

    def zeta_coefficients(self, n):
        """
        Compute the first n coefficients of the Dedekind zeta function
        of this field as a Dirichlet series.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: NumberField(x^2+1, 'a').zeta_coefficients(10)
            [1, 1, 0, 1, 2, 0, 0, 1, 1, 2]
        """
        return self.pari_nf().dirzetak(n)



class NumberField_extension(NumberField_generic):
    """
    EXAMPLES:
        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(x^3 - 2)
        sage: t = K['x'].gen()
        sage: L.<b> = K.extension(t^2+t+a); L
        Extension by x^2 + x + a of the Number Field in a with defining polynomial x^3 - 2
    """
    def __init__(self, base, polynomial, name, latex_name=None, names=None):
        """
        Note: polynomial must be defined in the ring \code{K['x']}, where
        K is the base field.
        """
        if not names is None: name = names
        if not is_NumberField(base):
            raise TypeError, "base (=%s) must be a number field"%base
        if not isinstance(polynomial, polynomial_element.Polynomial):
            raise TypeError, "polynomial (=%s) must be a polynomial"%polynomial
        if name == base.variable_name():
            raise ValueError, "Base field and extension cannot have the same name"
        if polynomial.parent().base_ring() != base:
            raise ValueError, "The polynomial must be defined over the base field"

        # Generate the nf and bnf corresponding to the base field
        # defined as polynomials in y, e.g. for rnfisfree

        # Convert the polynomial defining the base field into a
        # polynomial in y to satisfy PARI's ordering requirements.
        # NOTE: This might not work properly if the base field is not
        #       defined by a polynomial in one variable.  But currently
        #       they are all defined in one variable, so no problem!

        Qx = base.polynomial().parent()
        Qy = (base.polynomial().base_ring())['y']
        phi = Qx.hom([Qy.gen()])
        base_polynomial_y = phi(base.polynomial())

        self.__base_nf = pari(base_polynomial_y).nfinit()
        self.__base_bnf = pari(base_polynomial_y).bnfinit()

        # Use similar methods to convert the polynomial defining the
        # relative extension into a polynomial in x, with y denoting
        # the generator of the base field.
        # NOTE: This should be rewritten if there is a way to extend
        #       homomorphisms K -> K' to homomorphisms K[x] -> K'[x].
        base_field_y = NumberField(base.polynomial(), 'y')
        Kx = base_field_y['x']
        i = base.hom([base_field_y.gen()]) # inclusion K -> K' with a -> y
        rel_coeffs = [i(c) for c in polynomial.coeffs()]
        polynomial_y = Kx(rel_coeffs)

        self.__pari_relative_polynomial = pari(str(polynomial_y))
        self.__rnf = self.__base_nf.rnfinit(self.__pari_relative_polynomial)

        self.__base_field = base
        NumberField_generic.__init__(self, self.absolute_polynomial(), name=name, latex_name=latex_name, check=False)

        self._assign_names(name)
        self.__relative_polynomial = polynomial
        self.__pari_bnf_certified = False

    def __repr__(self):
        return "Extension by %s of the Number Field in %s with defining polynomial %s"%(
            self.polynomial(), self.base_field().variable_name(),
            self.base_field().polynomial())

    def _latex_(self):
        r"""
        Return a \LaTeX representation of the extension.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^3 - 2)
            sage: t = K['x'].gen()
            sage: K.extension(t^2+t+a, 'b')._latex_()
            '\\mathbf{Q}[b,a]/(b^{2} + b + a, a^{3} - 2)'
        """
        return "%s[%s,%s]/(%s, %s)"%(latex(QQ), self.variable_name(), self.base_field().variable_name(), self.polynomial()._latex_(self.variable_name()), self.base_field().polynomial()._latex_(self.base_field().variable_name()))

    def __call__(self, x):
        """
        Coerce x into this number field.
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return number_field_element.NumberFieldElement(self, x.polynomial())
            if x.parent() == self.base_field():
                return self.__base_inclusion(x)

        if not isinstance(x, (int, long, rational.Rational,
                              integer.Integer, pari_gen,
                              polynomial_element.Polynomial,
                              list)):
            raise TypeError, "Cannot coerce %s into %s"%(x,self)

        return number_field_element.NumberFieldElement(self, x)

    def _coerce_impl(self, x):
        if isinstance(x, number_field_element.NumberFieldElement):
            if x.parent() == self:
                return x
            if x.parent() == self.base_field():
                return self.__base_inclusion(x)
        elif isinstance(x, (rational.Rational, integer.Integer, int, long)):
            return number_field_element.NumberFieldElement(self, x)
        raise TypeError

    def __base_inclusion(self, element):
        """
        Given an element of the base field, give its inclusion into this
        extension (according to PARI's rnfeltreltoabs) in terms of the
        generator of this field.
        """
        if not number_field_element.is_NumberFieldElement(element):
            raise TypeError, "element must be a NumberFieldElement"
        if element.parent() != self.base_field():
            raise TypeError, "element must belong to the base field"
        base_field_y = NumberField(self.base_field().polynomial(), 'y')
        phi = self.base_field().hom([base_field_y.gen()])
        expr_x = self.pari_rnf().rnfeltreltoabs(str(phi(element)))

        # Convert to a polynomial in x, then to one in gen(), and return it
        return self(QQ['x'](str(expr_x).replace('^','**')))

    def _ideal_class_(self):
        return sage.rings.number_field.number_field_ideal.NumberFieldIdeal_rel

    def _pari_base_bnf(self, certify=False):
        # No need to certify the same field twice, so we'll just check
        # that the base field is certified.
        if certify:
            self.base_field().pari_bnf_certify()
        return self.__base_bnf

    def _pari_base_nf(self):
        return self.__base_nf

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "Only one generator."
        try:
            return self.__gen
        except AttributeError:
            X = rational_field.RationalField()['x'].gen()
            self.__gen = number_field_element.NumberFieldElement(self, X)
            return self.__gen

    def gen_relative(self):
        """
        Return root of defining polynomial, which is a generator of
        the relative number field over the base.
        """
        try:
            return self.__gen_relative
        except AttributeError:
            rnf = self.pari_rnf()
            f = (pari('x') - rnf[10][2]*rnf[10][1]).lift()
            self.__gen_relative = number_field_element.NumberFieldElement(self, f)
            return self.__gen_relative

            if self.__polynomial != None:
                X = self.__polynomial.parent().gen()
            else:
                X = PolynomialRing(rational_field.RationalField()).gen()
            self.__gen_relative = number_field_element.NumberFieldElement(self, X)
            return self.__gen_relative

    def pari_polynomial(self):
        """
        PARI polynomial corresponding to polynomial that defines
        this field.
        """
        try:
            return self.__pari_polynomial
        except AttributeError:
            self.__pari_polynomial = self.absolute_polynomial()._pari_()
            return self.__pari_polynomial

    def pari_rnf(self):
        return self.__rnf

    def pari_relative_polynomial(self):
        return self.__pari_relative_polynomial

    def absolute_field(self, name=None):
        r"""
        Return this field as an extension of $\Q$ rather than an
        extension of the base field.
        """
        try:
            return self.__absolute_field
        except AttributeError:
            if name is None:
                name = self.variable_name()
            self.__absolute_field = NumberField(self.absolute_polynomial(), name)
            return self.__absolute_field

    def absolute_polynomial(self):
        r"""
        Return the polynomial over $\Q$ which defines this field as an
        extension of the rational numbers.
        """
        try:
            return self.__absolute_polynomial
        except AttributeError:
            pbn = self._pari_base_nf()
            prp = self.pari_relative_polynomial()
            pari_poly = str(pbn.rnfequation(prp)).replace('^', '**')
            R = self.base_field().polynomial().parent()
            self.__absolute_polynomial = R(pari_poly)
            return self.__absolute_polynomial

    def base_field(self):
        return self.__base_field

    def base_ring(self):
        return self.base_field()

    def discriminant(self, certify=True):
        """
        Return the relative discriminant of this extension $L/K$ as
        an ideal of $K$.  If you want the (rational) discriminant of
        $L/Q$, use e.g. \code{L.absolute_field().discriminant()}.

        Note that this uses PARI's \code{rnfdisc} function, which
        according to the documentation takes an \code{nf} parameter in
        GP but a \code{bnf} parameter in the C library.  If the C
        library actually accepts an \code{nf}, then this function
        should be fixed and the \code{certify} parameter removed.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: K.<i> = NumberField(x^2+1)
            sage: t = K['x'].gen()
            sage: L.<b> = K.extension(t^4-i)
            sage: L.discriminant()
            Fractional ideal (256) of Number Field in i with defining polynomial x^2 + 1
        """
        bnf = self._pari_base_bnf(certify)
        K = self.base_field()
        R = K.polynomial().parent()
        D, d = bnf.rnfdisc(self.pari_relative_polynomial())
        return K.ideal([ K(R(x)) for x in convert_from_zk_basis(K, D) ])

    disc = discriminant

    def extension(self, poly, name='b'):
        """
        Raise a NotImplemented error, since relative extensions of relative
        extensions are not yet supported.
        """
        raise NotImplementedError, "relative extensions of relative extensions are not supported"

    def galois_group(self, pari_group = False, use_kash=False):
        r"""
        Return the Galois group of the Galois closure of this number
        field as an abstract group.  Note that even though this is an
        extension $L/K$, the group will be computed as if it were $L/\Q$.

        For more (important!) documentation, so the documentation
        for Galois groups of polynomials over $\Q$, e.g., by
        typing \code{K.polynomial().galois_group?}, where $K$
        is a number field.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^2 + 1)
            sage: R.<t> = PolynomialRing(K)
            sage: L = K.extension(t^5-t+a, 'b')
            sage: L.galois_group()                     # optional
            Transitive group number 22 of degree 10
        """
        return self.absolute_polynomial().galois_group(pari_group = pari_group, use_kash = use_kash)

    def is_free(self, certify=True):
        r"""
        Determine whether or not $L/K$ is free (i.e. if $\mathcal{O}_L$ is
        a free $\mathcal{O}_K$-module).

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^2+6)
            sage: L.<b> = K.extension(K['x'].gen()^2 + 3)    ## extend by x^2+3
            sage: L.is_free()
            False
        """
        base_bnf = self._pari_base_bnf(certify)
        if base_bnf.rnfisfree(self.pari_relative_polynomial()) == 1:
            return True
        return False

    def lift_to_base(self, element):
        """
        Lift an element of this extension into the base field if possible,
        or raise a ValueError if it is not possible.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^3 - 2, 'a')
            sage: R = K['x']
            sage: L = K.extension(R.gen()^2 - K.gen(), 'b')
            sage: b = L.gen()
            sage: L.lift_to_base(b^4)
            a^2
            sage: L.lift_to_base(b)
            Traceback (most recent call last):
            ...
            ValueError: The element b is not in the base field
        """
        poly_xy = self.pari_rnf().rnfeltabstorel( self(element)._pari_() )
        if str(poly_xy).find('x') >= 0:
            raise ValueError, "The element %s is not in the base field"%element
        return self.base_field()( QQ['y'](poly_xy) )

    def polynomial(self):
        return self.__relative_polynomial



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
        zeta6
        sage: z^3
        -1
        sage: (1+z)^3
        6*zeta6 - 3

        sage: K = CyclotomicField(197)
        sage: loads(K.dumps()) == K
        True
        sage: loads((z^2).dumps()) == z^2
        True

        sage: cf12 = CyclotomicField( 12 )
        sage: z12 = cf12.0
        sage: cf6 = CyclotomicField( 6 )
        sage: z6 = cf6.0
        sage: FF = Frac( cf12['x'] )
        sage: x = FF.0
        sage: print z6*x^3/(z6 + x)
        zeta12^2*x^3/(x + zeta12^2)
    """
    def __init__(self, n, names):
        f = QQ['x'].cyclotomic_polynomial(n)
        if names[0][:4] == 'zeta':
            latex_name = "\\zeta_{%s}"%n
        else:
            latex_name = None
        NumberField_generic.__init__(self, f,
                                     name= names,
                                     latex_name=latex_name,
                                     check=False)
        n = integer.Integer(n)
        zeta = self.gen()
        zeta._set_multiplicative_order(n)
        self.__zeta_order = n

    def __repr__(self):
        return "Cyclotomic Field of order %s and degree %s"%(
                self.zeta_order(), self.degree())

    def _latex_(self):
        return "%s(\\zeta_{%s})"%(latex(QQ), self.__zeta_order)

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
            zeta42^7
            sage: k6(b)
            zeta6
            sage: b^2
            zeta42^7 - 1
            sage: k6(b^2)
            zeta6 - 1

        Coercion of GAP cyclotomic elements is also fully supported.


        """
        if isinstance(x, number_field_element.NumberFieldElement):
            if isinstance(x.parent(), NumberField_cyclotomic):
                return self._coerce_from_other_cyclotomic_field(x)
            else:
                return self._coerce_from_other_number_field(x)
        elif sage.interfaces.gap.is_GapElement(x):
            return self._coerce_from_gap(x)
        else:
            return self._coerce_non_number_field_element_in(x)

    def _coerce_from_other_cyclotomic_field(self, x, only_canonical=False):
        """
        Coerce an element x of a cyclotomic field into self, if at all possible.

        INPUT:
            x -- number field element
            only_canonical -- bool (default: False); Attempt to work, even in some
                   cases when x is not in a subfield of the cyclotomics (as long as x is
                   a root of unity).
        """
        K = x.parent()
        if K is self:
            return x
        elif K == self:
            return number_field_element.NumberFieldElement(self, x.polynomial())
        n = K.zeta_order()
        m = self.zeta_order()
        if m % n == 0:   # easy case
            e = m // n
            f = x.polynomial()
            X = f.parent().gen()
            g = f(X**e)
        else:
            if only_canonical:
                raise TypeError
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

    def _coerce_from_gap(self, x):
        """
        Attempt to coerce a GAP number field element into this cyclotomic field.
        """
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
        s = sage.misc.all.sage_eval(s, locals={'zeta':K.gen()})
        if K is self:
            return s
        else:
            return self(s)

    def _coerce_impl(self, x):
        """
        Canonical coercion of x into self.

        Elements of other compatible cyclotomic fields coerce in, as do elements
        of the rings that coerce to all number fields (e.g., integers, rationals).
        """
        if isinstance(x, number_field_element.NumberFieldElement) and \
                isinstance(x.parent(), NumberField_cyclotomic):
            return self._coerce_from_other_cyclotomic_field(x, only_canonical=True)
        return NumberField_generic._coerce_impl(self, x)

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
              Defn: zeta4 |--> 0.0000000000000000612323399573676 + 1.00000000000000*I

        Note in the example above that the way zeta is computed (using
        sin and cosine in MPFR) means that only the prec bits of the
        number after the decimal point are valid.

            sage: K = CyclotomicField(3)
            sage: phi = K.complex_embedding (10)
            sage: phi(K.0)
            -0.49 + 0.86*I
            sage: phi(K.0^3)
            1.0
            sage: phi(K.0^3 - 1)
            0
            sage: phi(K.0^3 + 7)
            8.0
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        return self.hom([CC.zeta(self.zeta_order())], check=False)

    def complex_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the
        approximate complex field with precision prec.

        EXAMPLES:
            sage: C = CyclotomicField(4)
            sage: C.complex_embeddings()
            [Ring morphism:
                  From: Cyclotomic Field of order 4 and degree 2
                  To:   Complex Field with 53 bits of precision
                  Defn: zeta4 |--> 0.0000000000000000612323399573676 + 1.00000000000000*I, Ring morphism:
                  From: Cyclotomic Field of order 4 and degree 2
                  To:   Complex Field with 53 bits of precision
                  Defn: zeta4 |--> -0.000000000000000183697019872102 - 1.00000000000000*I]
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        n = self.zeta_order()
        z = CC.zeta(self.zeta_order())
        X = [m for m in range(n) if sage.rings.arith.gcd(m,n) == 1]
        return [self.hom([z**n], check=False) for n in X]

    def next_split_prime(self, p=2):
        """
        Return the next prime integer $p$ that splits completely in
        this cyclotomic field (and does not ramify).

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: K.next_split_prime(7)
            13
        """
        n = self.zeta_order()
        while True:
            p = sage.rings.arith.next_prime(p)
            if p % n == 1:
                return p

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
            zeta7
            sage: k.zeta().multiplicative_order()
            7
            sage: k = CyclotomicField(49)
            sage: k.zeta().multiplicative_order()
            49
            sage: k.zeta(7).multiplicative_order()
            7
            sage: k.zeta()
            zeta49
            sage: k.zeta(7)
            zeta49^7
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
        sage: QuadraticField(-4, 'b')
        Number Field in b with defining polynomial x^2 + 4
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

        \note{Computed using PARI via Schertz's method.  This
        implementation is quite fast.}

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^2 + 23, 'a')
            sage: K.hilbert_class_polynomial()
            x^3 + x^2 - 1

            sage: K = NumberField(x^2 + 431, 'a')
            sage: K.hilbert_class_polynomial()
            x^21 + x^20 - 13*x^19 - 50*x^18 + 592*x^17 - 2403*x^16 + 5969*x^15 - 10327*x^14 + 13253*x^13 - 12977*x^12 + 9066*x^11 - 2248*x^10 - 5523*x^9 + 11541*x^8 - 13570*x^7 + 11315*x^6 - 6750*x^5 + 2688*x^4 - 577*x^3 + 9*x^2 + 15*x + 1
        """
        f = pari('quadhilbert(%s))'%self.discriminant())
        g = QQ['x'](f)
        return g

    def hilbert_class_field(self, names):
        r"""
        Returns the Hilbert class field of this quadratic
        field as an absolute extension of $\Q$.  For a polynomial
        that defines a relative extension see the
        \code{hilbert_class_polynomial} command.

        \note{Computed using PARI via Schertz's method.  This implementation
        is amazingly fast.}

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^2 + 23, 'a')
            sage: K.hilbert_class_polynomial()
            x^3 + x^2 - 1
            sage: K.hilbert_class_field('h')
            Number Field in h with defining polynomial x^6 + 2*x^5 + 70*x^4 + 90*x^3 + 1631*x^2 + 1196*x + 12743
        """
        f = self.hilbert_class_polynomial()
        C = self.composite_fields(NumberField(f,'x'),names)
        assert len(C) == 1
        return C[0]

def is_fundamental_discriminant(D):
    d = D % 4
    if not (d in [0,1]):
        return False
    return D != 1 and  D != 0 and \
           (arith.is_squarefree(D) or \
            (d == 0 and (D//4)%4 in [2,3] and arith.is_squarefree(D//4)))
