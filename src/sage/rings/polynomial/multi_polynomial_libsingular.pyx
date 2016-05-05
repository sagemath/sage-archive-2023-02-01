r"""
Multivariate Polynomials via libSINGULAR

This module implements specialized and optimized implementations for
multivariate polynomials over many coefficient rings, via a shared library
interface to SINGULAR. In particular, the following coefficient rings are
supported by this implementation:

- the rational numbers `\QQ`,

- the ring of integers `\ZZ`,

- `\ZZ/n\ZZ` for any integer `n`,

- finite fields `\GF{p^n}` for `p` prime and `n > 0`,

- and absolute number fields `\QQ(a)`.

AUTHORS:

The libSINGULAR interface was implemented by

- Martin Albrecht (2007-01): initial implementation

- Joel Mohler (2008-01): misc improvements, polishing

- Martin Albrecht (2008-08): added `\QQ(a)` and `\ZZ` support

- Simon King (2009-04): improved coercion

- Martin Albrecht (2009-05): added `\ZZ/n\ZZ` support, refactoring

- Martin Albrecht (2009-06): refactored the code to allow better
  re-use

- Simon King (2011-03): Use a faster way of conversion from the base
  ring.

- Volker Braun (2011-06): Major cleanup, refcount singular rings, bugfixes.

TODO:

- implement Real, Complex coefficient rings via libSINGULAR

EXAMPLES:

We show how to construct various multivariate polynomial rings::

    sage: P.<x,y,z> = QQ[]
    sage: P
    Multivariate Polynomial Ring in x, y, z over Rational Field

    sage: f = 27/113 * x^2 + y*z + 1/2; f
    27/113*x^2 + y*z + 1/2

    sage: P.term_order()
    Degree reverse lexicographic term order

    sage: P = PolynomialRing(GF(127),3,names='abc', order='lex')
    sage: P
    Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

    sage: a,b,c = P.gens()
    sage: f = 57 * a^2*b + 43 * c + 1; f
    57*a^2*b + 43*c + 1

    sage: P.term_order()
    Lexicographic term order

    sage: z = QQ['z'].0
    sage: K.<s> = NumberField(z^2 - 2)
    sage: P.<x,y> = PolynomialRing(K, 2)
    sage: 1/2*s*x^2 + 3/4*s
    (1/2*s)*x^2 + (3/4*s)

    sage: P.<x,y,z> = ZZ[]; P
    Multivariate Polynomial Ring in x, y, z over Integer Ring

    sage: P.<x,y,z> = Zmod(2^10)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

    sage: P.<x,y,z> = Zmod(3^10)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

    sage: P.<x,y,z> = Zmod(2^100)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1267650600228229401496703205376

    sage: P.<x,y,z> = Zmod(2521352)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352
    sage: type(P)
    <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

    sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 25213521351515232
    sage: type(P)
    <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_with_category'>

We construct the Frobenius morphism on `\GF{5}[x,y,z]` over `\GF{5}`::

    sage: R.<x,y,z> = PolynomialRing(GF(5), 3)
    sage: frob = R.hom([x^5, y^5, z^5])
    sage: frob(x^2 + 2*y - z^4)
    -z^20 + x^10 + 2*y^5
    sage: frob((x + 2*y)^3)
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15
    sage: (x^5 + 2*y^5)^3
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15

We make a polynomial ring in one variable over a polynomial ring in
two variables::

    sage: R.<x, y> = PolynomialRing(QQ, 2)
    sage: S.<t> = PowerSeriesRing(R)
    sage: t*(x+y)
    (x + y)*t

TESTS::

    sage: P.<x,y,z> = QQ[]
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = GF(2^8,'a')[]
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = GF(127)[]
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = GF(127)[]
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True

    sage: Rt.<t> = PolynomialRing(QQ,1)
    sage: p = 1+t
    sage: R.<u,v> = PolynomialRing(QQ, 2)
    sage: p(u/v)
    (u + v)/v

Check if :trac:`6160` is fixed::

    sage: x=var('x')
    sage: K.<j> = NumberField(x-1728)
    sage: R.<b,c> = K[]
    sage: b-j*c
    b - 1728*c
"""

# The Singular API is as follows:
#
#   pXXX does assume the currRing to be set
#   p_XXX does not.
#
# However, sometimes there are bugs, i.e. you call p_XXX and it
# crashes unless currRing is set.
#
# Notable exceptions:
#   * pNext and pIter don't need currRing
#   * p_Normalize apparently needs currRing

include "cysignals/memory.pxi"
include "cysignals/signals.pxi"

# singular types
from sage.libs.singular.decl cimport ring, poly, ideal, intvec, number, currRing

# singular functions
from sage.libs.singular.decl cimport (
    errorreported, libfac_interruptflag,
    p_ISet, rChangeCurrRing, p_Copy, p_Init, p_SetCoeff, p_Setm, p_SetExp, p_Add_q,
    p_NSet, p_GetCoeff, p_Delete, p_GetExp, pNext, rRingVar, omAlloc0, omStrDup,
    omFree, pDivide, p_SetCoeff0, n_Init, p_DivisibleBy, pLcm, p_LmDivisibleBy,
    pDivide, p_IsConstant, p_ExpVectorEqual, p_String, p_LmInit, n_Copy,
    p_IsUnit, pInvers, p_Head, idInit, fast_map, id_Delete,
    pIsHomogeneous, pHomogen, p_Totaldegree, singclap_pdivide, singclap_factorize,
    idLift, IDELEMS, On, Off, SW_USE_CHINREM_GCD, SW_USE_EZGCD,
    p_LmIsConstant, pTakeOutComp1, singclap_gcd, pp_Mult_qq, p_GetMaxExp,
    pLength, kNF, singclap_isSqrFree, p_Neg, p_Minus_mm_Mult_qq, p_Plus_mm_Mult_qq,
    pDiff, singclap_resultant, p_Normalize,
    prCopyR, prCopyR_NoSort )

# singular conversion routines
from sage.libs.singular.singular cimport si2sa, sa2si, overflow_check

# singular poly arith
from sage.libs.singular.polynomial cimport (
    singular_polynomial_call, singular_polynomial_cmp, singular_polynomial_add,
    singular_polynomial_sub, singular_polynomial_neg, singular_polynomial_rmul,
    singular_polynomial_mul, singular_polynomial_div_coeff, singular_polynomial_pow,
    singular_polynomial_str, singular_polynomial_latex,
    singular_polynomial_str_with_changed_varnames, singular_polynomial_deg,
    singular_polynomial_length_bounded, singular_polynomial_subst )

# singular rings
from sage.libs.singular.ring cimport singular_ring_new, singular_ring_reference, singular_ring_delete

# polynomial imports
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict, MPolynomialRing_polydict_domain
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.polydict import ETuple
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing

# base ring imports
from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.rational cimport Rational
from sage.rings.rational_field import RationalField
from sage.rings.complex_field import is_ComplexField
from sage.rings.real_mpfr import is_RealField
from sage.rings.integer_ring import is_IntegerRing, ZZ
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
from sage.rings.number_field.number_field_base cimport NumberField

from sage.arith.all import gcd
from sage.structure.element import coerce_binop

from sage.structure.parent cimport Parent
from sage.structure.category_object cimport CategoryObject

from sage.structure.element cimport EuclideanDomainElement
from sage.structure.element cimport RingElement
from sage.structure.element cimport ModuleElement
from sage.structure.element cimport Element
from sage.structure.element cimport CommutativeRingElement
from sage.structure.element cimport coercion_model

from sage.structure.factorization import Factorization
from sage.structure.sequence import Sequence

from sage.interfaces.all import macaulay2
from sage.interfaces.singular import singular as singular_default, is_SingularElement, SingularElement
from sage.interfaces.macaulay2 import macaulay2 as macaulay2_default, is_Macaulay2Element

from sage.misc.all import prod as mul
from sage.misc.sage_eval import sage_eval

import sage.libs.pari.gen
import polynomial_element

permstore=[]
cdef class MPolynomialRing_libsingular(MPolynomialRing_generic):

    def __cinit__(self):
        """
        The Cython constructor.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: MPolynomialRing_libsingular(QQ, 3, ('x', 'y', 'z'), TermOrder('degrevlex', 3))
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: type(_)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

            sage: P.<x,y,z> = QQ[]; P
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: type(P)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
        """
        self._ring = NULL

    def __init__(self, base_ring, n, names, order='degrevlex'):
        """
        Construct a multivariate polynomial ring subject to the
        following conditions:

        INPUT:

        - ``base_ring`` - base ring (must be either GF(q), ZZ, ZZ/nZZ,
                          QQ or absolute number field)

        - ``n`` - number of variables (must be at least 1)

        - ``names`` - names of ring variables, may be string of list/tuple

        - ``order`` - term order (default: ``degrevlex``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P
            Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: f = 27/113 * x^2 + y*z + 1/2; f
            27/113*x^2 + y*z + 1/2

            sage: P.term_order()
            Degree reverse lexicographic term order

            sage: P = PolynomialRing(GF(127),3,names='abc', order='lex')
            sage: P
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1; f
            57*a^2*b + 43*c + 1

            sage: P.term_order()
            Lexicographic term order

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)
            sage: P.<x,y> = PolynomialRing(K, 2)
            sage: 1/2*s*x^2 + 3/4*s
            (1/2*s)*x^2 + (3/4*s)

            sage: P.<x,y,z> = ZZ[]; P
            Multivariate Polynomial Ring in x, y, z over Integer Ring

            sage: P.<x,y,z> = Zmod(2^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

            sage: P.<x,y,z> = Zmod(3^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

            sage: P.<x,y,z> = Zmod(2^100)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1267650600228229401496703205376

            sage: P.<x,y,z> = Zmod(2521352)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352
            sage: type(P)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

            sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 25213521351515232
            sage: type(P)
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_with_category'>

            sage: P.<x,y,z> = PolynomialRing(Integers(2^32),order='lex')
            sage: P(2^32-1)
            4294967295

        TEST:

        Make sure that a faster coercion map from the base ring is used;
        see :trac:`9944`::

            sage: R.<x,y> = PolynomialRing(ZZ)
            sage: R.coerce_map_from(R.base_ring())
            Polynomial base injection morphism:
              From: Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Integer Ring
        """
        MPolynomialRing_generic.__init__(self, base_ring, n, names, order)
        self._has_singular = True
        assert(n == len(self._names))
        self.__ngens = n
        self._ring = singular_ring_new(base_ring, n, self._names, order)
        self._zero_element = new_MP(self, NULL)
        cdef MPolynomial_libsingular one = new_MP(self, p_ISet(1, self._ring))
        self._one_element = one
        self._one_element_poly = one._poly
        # This polynomial ring should belong to Algebras(base_ring).
        # Algebras(...).parent_class, which was called from MPolynomialRing_generic.__init__,
        # tries to provide a conversion from the base ring, if it does not exist.
        # This is for algebras that only do the generic stuff in their initialisation.
        # But here, we want to use PolynomialBaseringInjection. Hence, we need to
        # wipe the memory and construct the conversion from scratch.
        from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
        base_inject = PolynomialBaseringInjection(base_ring, self)
        self.register_coercion(base_inject)
        #permanently store a reference to this ring until deallocation works reliably
        permstore.append(self)

    def __dealloc__(self):
        r"""
        Deallocate the ring without changing ``currRing``

        TESTS:

        This example caused a segmentation fault with a previous version
        of this method::

            sage: import gc
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R1 = MPolynomialRing_libsingular(GF(5), 2, ('x', 'y'), TermOrder('degrevlex', 2))
            sage: R2 = MPolynomialRing_libsingular(GF(11), 2, ('x', 'y'), TermOrder('degrevlex', 2))
            sage: R3 = MPolynomialRing_libsingular(GF(13), 2, ('x', 'y'), TermOrder('degrevlex', 2))
            sage: _ = gc.collect()
            sage: foo = R1.gen(0)
            sage: del foo
            sage: del R1
            sage: _ = gc.collect()
            sage: del R2
            sage: _ = gc.collect()
            sage: del R3
            sage: _ = gc.collect()
        """
        if self._ring != NULL:  # the constructor did not raise an exception
            singular_ring_delete(self._ring)

    def __copy__(self):
        """
        Copy ``self``.

        The ring is unique and immutable, so we do not copy.

        TESTS::

            sage: import gc
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: from sage.libs.singular.ring import ring_refcount_dict
            sage: gc.collect()  # random output
            sage: n = len(ring_refcount_dict)
            sage: R = MPolynomialRing_libsingular(GF(547), 2, ('x', 'y'), TermOrder('degrevlex', 2))
            sage: len(ring_refcount_dict) == n + 1
            True

            sage: Q = copy(R)   # indirect doctest
            sage: p = R.gen(0) ^2+R.gen(1)^2
            sage: q = copy(p)
            sage: del R
            sage: del Q
            sage: del p
            sage: del q
            sage: gc.collect() # random output
            sage: len(ring_refcount_dict) == n
            False
        """
        return self

    def __deepcopy__(self, memo):
        """
        Deep copy ``self``.

        The ring should be immutable, so we do not copy.

        TESTS::

            sage: R.<x,y> = GF(547)[]
            sage: R is deepcopy(R)   # indirect doctest
            True
        """
        memo[id(self)] = self
        return self

    cpdef _coerce_map_from_(self, other):
        """
        Return True if and only if there exists a coercion map from
        ``other`` to ``self``.

        TEST::

            sage: R.<x,y> = QQ[]
            sage: type(R)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
            sage: R.has_coerce_map_from(ZZ['t'])
            False
            sage: R.coerce_map_from(ZZ['x'])
            Conversion map:
            From: Univariate Polynomial Ring in x over Integer Ring
            To:   Multivariate Polynomial Ring in x, y over Rational Field

        """
        base_ring = self.base_ring()

        if isinstance(other, MPolynomialRing_libsingular):
            if self is other:
                return True
            n = other.ngens()
            if(other.base_ring is base_ring and self.ngens() >= n and
               self.variable_names()[:n] == other.variable_names()):
                return True
            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
                return True
        elif isinstance(other, MPolynomialRing_polydict):
            if self == other:
                return True
            elif other.ngens() == 0:
                return True
            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
                return True
        elif is_PolynomialRing(other):
            if base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
                return True
        elif base_ring.has_coerce_map_from(other):
            return True

    Element = MPolynomial_libsingular

    def _element_constructor_(self, element):
        """
        Construct a new element in this polynomial ring by converting
        ``element`` into ``self`` if possible.

        INPUT:

        - ``element`` -- several types are supported, see below

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]

        We can coerce elements of self to self::

            sage: P._coerce_(x*y + 1/2)
            x*y + 1/2

        We can coerce elements for a ring with the same algebraic properties::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P == R
            True

            sage: P is R
            False

            sage: P._coerce_(x*y + 1)
            x*y + 1

        We can coerce base ring elements::

            sage: P._coerce_(3/2)
            3/2

        and all kinds of integers::

            sage: P._coerce_(ZZ(1))
            1

            sage: P._coerce_(int(1))
            1

            sage: k.<a> = GF(2^8)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: P._coerce_(a)
            (a)

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)
            sage: P.<x,y> = PolynomialRing(K, 2)
            sage: P._coerce_(1/2*s)
            (1/2*s)

        TESTS::

            sage: P.<x,y> = PolynomialRing(GF(127))
            sage: P("111111111111111111111111111111111111111111111111111111111")
            21
            sage: P.<x,y> = PolynomialRing(QQ)
            sage: P("111111111111111111111111111111111111111111111111111111111")
            111111111111111111111111111111111111111111111111111111111
            sage: P("31367566080")
            31367566080

        Check if :trac:`7582` is fixed::

            sage: R.<x,y,z> = PolynomialRing(CyclotomicField(2),3)
            sage: R.coerce(1)
            1

        Check if :trac:`6160` is fixed::

            sage: x=var('x')
            sage: K.<j> = NumberField(x-1728)
            sage: R.<b,c> = K[]
            sage: R.coerce(1)
            1

        Check if coercion from zero variable polynomial rings work (#7951)::

            sage: P = PolynomialRing(QQ,0,'')
            sage: R.<x,y> = QQ[]
            sage: P(5)*x
            5*x
            sage: P = PolynomialRing(ZZ,0,'')
            sage: R.<x,y> = GF(127)[]
            sage: R.coerce(P(5))
            5

        Conversion from strings::

            sage: P.<x,y,z> = QQ[]
            sage: P('x+y + 1/4')
            x + y + 1/4

        Coercion from SINGULAR elements::

            sage: P._singular_()
            //   characteristic : 0
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_().set_ring()
            sage: P(singular('x + 3/4'))
            x + 3/4

        Coercion from symbolic variables::

            sage: R = QQ['x,y,z']
            sage: var('x')
            x
            sage: R(x)
            x

        Coercion from 'similar' rings, which maps by index::

            sage: P.<x,y,z> = QQ[]
            sage: R.<a,b,c> = ZZ[]
            sage: P(a)
            x

        ::

            sage: P.<x,y> = QQ[]
            sage: R.<a,b,c> = QQ[]
            sage: R(x)
            a

        Coercion from PARI objects::

            sage: P.<x,y,z> = QQ[]
            sage: P(pari('x^2 + y'))
            x^2 + y
            sage: P(pari('x*y'))
            x*y

        Coercion from boolean polynomials, also by index::

            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.<x,y,z> = QQ[]
            sage: P(B.gen(0))
            x

        If everything else fails, we try to coerce to the base ring::

            sage: R.<x,y,z> = GF(3)[]
            sage: R(1/2)
            -1

        Finally, conversions from other polynomial rings which are not
        coercions are provided. Variables are mapped as follows. Say,
        we are mapping an element from `P` to `Q` (this ring). If the
        variables of `P` are a subset of `Q`, we perform a name
        preserving conversion::

            sage: P.<y_2, y_1, z_3, z_2, z_1> = GF(3)[]
            sage: Q = GF(3)['y_4', 'y_3', 'y_2', 'y_1', 'z_5', 'z_4', 'z_3', 'z_2', 'z_1']
            sage: Q(y_1*z_2^2*z_1)
            y_1*z_2^2*z_1

        Otherwise, if `P` has less than or equal the number of
        variables as `Q`, we perform a conversion by index::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a), f(b), f(c)
            (c, b, d)

        ::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d']
            sage: f = Q.convert_map_from(P)
            sage: f(a),f(b),f(c)
            (c, b, d)

        In all other cases, we fail::

            sage: P.<a,b,c,f> = GF(2)[]
            sage: Q = GF(2)['c','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a)
            Traceback (most recent call last):
            ...
            TypeError: Could not find a mapping of the passed element to this ring.

        Coerce in a polydict where a coefficient reduces to 0 but isn't 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S( (5*x*y + x + 17*y)._mpoly_dict_recursive() )
            xx + 2*yy

        Coerce in a polynomial one of whose coefficients reduces to 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S(5*x*y + x + 17*y)
            xx + 2*yy

        Some other examples that illustrate the same coercion idea::

            sage: R.<x,y> = ZZ[]
            sage: S.<xx,yy> = GF(25,'a')[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

            sage: S.<xx,yy> = Integers(5)[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

        See :trac:`5292`::

            sage: R.<x> = QQ[]; S.<q,t> = QQ[]; F = FractionField(S);
            sage: x in S
            False
            sage: x in F
            False

        Check if :trac:`8228` is fixed::

            sage: P.<x,y> = Zmod(10)[]; P(0)
            0
            sage: P.<x,y> = Zmod(2^10)[]; P(0)
            0

        And :trac:`7597` is fixed if this does not segfault::

            sage: F2 = GF(2)
            sage: F.<x> = GF(2^8)
            sage: R4.<a,b> = PolynomialRing(F)
            sage: R.<u,v> = PolynomialRing(F2)
            sage: P = a
            sage: (P(0,0).polynomial()[0])*u
            0
            sage: P(a,b)
            a

        Check that :trac:`15746` is fixed::

            sage: R.<x,y> = GF(7)[]
            sage: R(2^31)
            2

        Check that :trac:`17964` is fixed::

            sage: K.<a> = QuadraticField(17)
            sage: Q.<x,y> = K[]
            sage: f = (-3*a)*y + (5*a)
            sage: p = K.primes_above(5)[0]
            sage: R = K.residue_field(p)
            sage: S = R['x','y']
            sage: S(f)
            (2*abar)*y

        """
        cdef poly *_p
        cdef poly *mon
        cdef poly *El_poly
        cdef ring *_ring = self._ring
        cdef number *_n
        cdef ring *El_ring
        cdef long mpos
        cdef MPolynomial_libsingular Element
        cdef MPolynomialRing_libsingular El_parent
        cdef unsigned int i, j
        cdef list ind_map = []
        cdef int e
        if _ring!=currRing: rChangeCurrRing(_ring)

        base_ring = self.base_ring()

        if isinstance(element, MPolynomial_libsingular):
            n = (<MPolynomial_libsingular>element)._parent.ngens()
            if element.parent() is self:
                return element
            elif(base_ring is element.base_ring() and
                 self.ngens() >= n and
                 self.variable_names()[:n] == (<MPolynomial_libsingular>element)._parent.variable_names()):
                if self.term_order() == (<MPolynomial_libsingular>element)._parent.term_order():
                    _p = prCopyR_NoSort((<MPolynomial_libsingular>element)._poly,
                                        (<MPolynomial_libsingular>element)._parent_ring,
                                        _ring)
                else:
                    _p = prCopyR((<MPolynomial_libsingular>element)._poly,
                                 (<MPolynomial_libsingular>element)._parent_ring, _ring)
                return new_MP(self, _p)
            elif base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))

        elif isinstance(element, MPolynomial_polydict):
            if element.parent() == self:
                _p = p_ISet(0, _ring)
                for (m,c) in element.element().dict().iteritems():
                    mon = p_Init(_ring)
                    p_SetCoeff(mon, sa2si(c, _ring), _ring)
                    for pos in m.nonzero_positions():
                        overflow_check(m[pos], _ring)
                        p_SetExp(mon, pos+1, m[pos], _ring)
                    p_Setm(mon, _ring)
                    _p = p_Add_q(_p, mon, _ring)
                return new_MP(self, _p)
            elif element.parent().ngens() == 0:
                # zero variable polynomials
                _p = p_NSet(sa2si(base_ring(element[tuple()]), _ring),
                        _ring)
                return new_MP(self, _p)
            elif base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))

        elif isinstance(element, polynomial_element.Polynomial):
            if base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))

        if isinstance(element, CommutativeRingElement):
            # base ring elements
            if element.parent() is base_ring:
                # shortcut for GF(p)
                if isinstance(base_ring, FiniteField_prime_modn):
                    _p = p_ISet(int(element) % _ring.ch, _ring)
                else:
                    _n = sa2si(element,_ring)
                    _p = p_NSet(_n, _ring)
                return new_MP(self, _p)
            # also accepting ZZ
            elif is_IntegerRing(element.parent()):
                if isinstance(base_ring, FiniteField_prime_modn):
                    _p = p_ISet(int(element),_ring)
                else:
                    _n = sa2si(base_ring(element),_ring)
                    _p = p_NSet(_n, _ring)
                return new_MP(self, _p)
            # fall back to base ring
            try:
                element = base_ring._coerce_c(element)
                _n = sa2si(element,_ring)
                _p = p_NSet(_n, _ring)
                return new_MP(self, _p)
            except TypeError:
                pass

        elif isinstance(element, int) or isinstance(element, long):
            if isinstance(base_ring, FiniteField_prime_modn):
                _p = p_ISet(element % _ring.ch, _ring)
            else:
                _n = sa2si(base_ring(element), _ring)
                _p = p_NSet(_n, _ring)
            return new_MP(self, _p)

        if isinstance(element, (SingularElement, sage.libs.pari.gen.gen)):
            element = str(element)

        if isinstance(element, MPolynomial_libsingular) and element.parent() is not self and element.parent() != self:
            variable_names_s = element.parent().variable_names()
            variable_names_t = self.variable_names()

            if set(variable_names_s).issubset(variable_names_t):
                for v in variable_names_s:
                    ind_map.append(variable_names_t.index(v)+1)
            else:
                ind_map = [i+1 for i in range(_ring.N)]

            if element.parent().ngens() <= self.ngens():
                # Map the variables by indices
                _p = p_ISet(0, _ring)
                K = self.base_ring()
                base = self._base
                Element = <MPolynomial_libsingular>element
                El_poly = Element._poly
                El_parent = Element._parent
                El_ring = Element._parent_ring
                El_base = El_parent._base

                while El_poly:
                    c = si2sa(p_GetCoeff(El_poly, El_ring), El_ring, El_base)
                    try:
                        c = K(c)
                    except TypeError:
                        p_Delete(&_p, _ring)
                        raise
                    if c:
                        mon = p_Init(_ring)
                        p_SetCoeff(mon, sa2si(c, _ring), _ring)
                        for j from 1 <= j <= El_ring.N:
                            e = p_GetExp(El_poly, j, El_ring)
                            if e:
                                p_SetExp(mon, ind_map[j-1], e, _ring)
                        p_Setm(mon, _ring)
                        _p = p_Add_q(_p, mon, _ring)
                    El_poly = pNext(El_poly)
                return new_MP(self, _p)

        if isinstance(element, MPolynomial_polydict):
            variable_names_s = element.parent().variable_names()
            variable_names_t = self.variable_names()

            if set(variable_names_s).issubset(variable_names_t):
                for v in variable_names_s:
                    ind_map.append(variable_names_t.index(v)+1)
            else:
                ind_map = [i+1 for i in range(_ring.N)]

            if element.parent().ngens() <= self.ngens():
                # Map variables by indices
                _p = p_ISet(0, _ring)
                K = self.base_ring()
                for (m,c) in element.element().dict().iteritems():
                    try:
                        c = K(c)
                        if not c: continue
                    except TypeError:
                        p_Delete(&_p, _ring)
                        raise
                    mon = p_Init(_ring)
                    p_SetCoeff(mon, sa2si(c , _ring), _ring)
                    for pos in m.nonzero_positions():
                        overflow_check(m[pos], _ring)
                        p_SetExp(mon, ind_map[pos], m[pos], _ring)
                    p_Setm(mon, _ring)
                    _p = p_Add_q(_p, mon, _ring)
                return new_MP(self, _p)

        from sage.rings.polynomial.pbori import BooleanPolynomial
        if isinstance(element, BooleanPolynomial):
            if element.constant():
                if element:
                    return self._one_element
                else:
                    return self._zero_element

            variable_names_s = set(element.parent().variable_names())
            variable_names_t = self.variable_names()

            if variable_names_s.issubset(variable_names_t):
                return eval(str(element),self.gens_dict())

            elif element.parent().ngens() <= self.ngens():
                Q = element.parent()
                gens_map = dict(zip(Q.variable_names(),self.gens()[:Q.ngens()]))
                return eval(str(element),gens_map)

        if isinstance(element, basestring):
            # let python do the parsing
            d = self.gens_dict()
            if self.base_ring().gen() != 1:
                d[str(self.base_ring().gen())]=self.base_ring().gen()
            try:
                if '/' in element:
                    element = sage_eval(element,d)
                else:
                    element = element.replace("^","**")
                    element = eval(element, d, {})
            except (SyntaxError, NameError):
                raise TypeError("Could not find a mapping of the passed element to this ring.")

            # we need to do this, to make sure that we actually get an
            # element in self.
            return self._coerce_c(element)

        if isinstance(element, dict):
            _p = p_ISet(0, _ring)
            K = self.base_ring()
            for (m,c) in element.iteritems():
                try:
                    c = K(c)
                    if not c: continue
                except TypeError:
                    p_Delete(&_p, _ring)
                    raise
                mon = p_Init(_ring)
                p_SetCoeff(mon, sa2si(c , _ring), _ring)
                if len(m) != self.ngens():
                    raise TypeError("tuple key must have same length as ngens")
                for pos from 0 <= pos < len(m):
                    if m[pos]:
                        overflow_check(m[pos], _ring)
                        p_SetExp(mon, pos+1, m[pos], _ring)
                p_Setm(mon, _ring)
                _p = p_Add_q(_p, mon, _ring)

            return new_MP(self, _p)

        try: #if hasattr(element,'_polynomial_'):
            # SymbolicVariable
            return element._polynomial_(self)
        except AttributeError:
            pass

        if is_Macaulay2Element(element):
            return self(element.external_string())

        try:
            return self(str(element))
        except TypeError:
            pass

        try:
            # now try calling the base ring's __call__ methods
            element = self.base_ring()(element)
            _p = p_NSet(sa2si(element,_ring), _ring)
            return new_MP(self,_p)
        except (TypeError,ValueError):
            raise TypeError("Could not find a mapping of the passed element to this ring.")

    def _repr_(self):
        """
        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: P # indirect doctest
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        varstr = ", ".join([ rRingVar(i,self._ring)  for i in range(self.__ngens) ])
        return "Multivariate Polynomial Ring in %s over %s"%(varstr,self.base_ring())

    def ngens(self):
        """
        Returns the number of variables in this multivariate polynomial ring.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: P.ngens()
            2

            sage: k.<a> = GF(2^16)
            sage: P = PolynomialRing(k,1000,'x')
            sage: P.ngens()
            1000
        """
        return int(self.__ngens)

    def gen(self, int n=0):
        """
        Returns the ``n``-th generator of this multivariate polynomial
        ring.

        INPUT:

        - ``n`` -- an integer ``>= 0``

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.gen(),P.gen(1)
            (x, y)

            sage: P = PolynomialRing(GF(127),1000,'x')
            sage: P.gen(500)
            x500

            sage: P.<SAGE,SINGULAR> = QQ[] # weird names
            sage: P.gen(1)
            SINGULAR
        """
        cdef poly *_p
        cdef ring *_ring = self._ring

        if n < 0 or n >= self.__ngens:
            raise ValueError, "Generator not defined."

        rChangeCurrRing(_ring)
        _p = p_ISet(1, _ring)
        p_SetExp(_p, n+1, 1, _ring)
        p_Setm(_p, _ring);

        return new_MP(self, _p)

    def ideal(self, *gens, **kwds):
        """
        Create an ideal in this polynomial ring.

        INPUT:

        - ``*gens`` - list or tuple of generators (or several input arguments)

        - ``coerce`` - bool (default: ``True``); this must be a
          keyword argument. Only set it to ``False`` if you are certain
          that each generator is already in the ring.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: sage.rings.ideal.Katsura(P)
            Ideal (x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y + 2*y*z - y) of Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: P.ideal([x + 2*y + 2*z-1, 2*x*y + 2*y*z-y, x^2 + 2*y^2 + 2*z^2-x])
            Ideal (x + 2*y + 2*z - 1, 2*x*y + 2*y*z - y, x^2 + 2*y^2 + 2*z^2 - x) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        coerce = kwds.get('coerce', True)
        if len(gens) == 1:
            gens = gens[0]
        from sage.rings.ideal import is_Ideal
        if is_Ideal(gens):
            if gens.ring() is self:
                return gens
            gens = gens.gens()
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
        Create an M2 representation of this polynomial ring if
        Macaulay2 is installed.

        INPUT:

        - ``macaulay2`` - M2 interpreter (default: ``macaulay2_default``)

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: macaulay2(R)        # optional - macaulay2
            ZZ [x, y, MonomialOrder => GRevLex, MonomialSize => 16]

            sage: R.<x,y> = QQ[]
            sage: macaulay2(R)        # optional - macaulay2, indirect doctest
            QQ [x, y, MonomialOrder => GRevLex, MonomialSize => 16]

            sage: R.<x,y> = GF(17)[]
            sage: print macaulay2(R)        # optional - macaulay2
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

    def _macaulay2_set_ring(self, macaulay2=macaulay2_default):
        """
        Set the associated M2 ring.

        INPUT:

        - ``macaulay2`` - M2 instance

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: M2 = P._macaulay2_set_ring() # optional - macaulay2
        """
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

    def _singular_(self, singular=singular_default):
        """
        Create a SINGULAR (as in the computer algebra system)
        representation of this polynomial ring. The result is cached.

        INPUT:

        - ``singular`` - SINGULAR interpreter (default: ``singular_default``)

        EXAMPLES::

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
            sage: P.<x> = QQ[]
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
            if not self.base_ring().is_field():
                return R
            if self.base_ring().is_prime_field():
                return R
            if self.base_ring().is_finite() \
                    or (isinstance(self.base_ring(), NumberField) and self.base_ring().is_absolute()):
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                    singular.eval("minpoly=%s"%(self.__minpoly))
                    self.__minpoly = singular.eval('minpoly')[1:-1] # store in correct format
            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self, singular=singular_default):
        """
        Create a SINGULAR (as in the computer algebra system)
        representation of this polynomial ring. The result is NOT
        cached.

        INPUT:

        - ``singular`` - SINGULAR interpreter (default: ``singular_default``)

        EXAMPLES::

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

            sage: w = var('w')
            sage: R.<x,y> = PolynomialRing(NumberField(w^2+1,'s'))
            sage: singular(R)
            //   characteristic : 0
            //   1 parameter    : s
            //   minpoly        : (s^2+1)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(2**8,'a'),10,'x', order='invlex')
            sage: singular(R)
            //   characteristic : 2
            //   1 parameter    : a
            //   minpoly        : (a^8+a^4+a^3+a^2+1)
            //   number of vars : 10
            //        block   1 : ordering rp
            //                  : names    x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127),2,'x', order='invlex')
            sage: singular(R)
            //   characteristic : 127
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ,2,'x', order='invlex')
            sage: singular(R)
            //   characteristic : 0
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ,2,'x', order='degneglex')
            sage: singular(R)
            //   characteristic : 0
            //   number of vars : 2
            //        block   1 : ordering a
            //                  : names    x0 x1
            //                  : weights   1  1
            //        block   2 : ordering ls
            //                  : names    x0 x1
            //        block   3 : ordering C

            sage: R = PolynomialRing(QQ,'x')
            sage: singular(R)
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127),'x')
            sage: singular(R)
            //   characteristic : 127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = ZZ['x,y']
            sage: singular(R)
            //   coeff. ring is : Integers
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = IntegerModRing(1024)['x,y']
            sage: singular(R)
            //   coeff. ring is : Z/2^10
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = IntegerModRing(15)['x,y']
            sage: singular(R)
            //   coeff. ring is : Z/15
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

        TESTS::

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
            order = self.term_order().singular_str()%dict(ngens=self.ngens())

        base_ring = self.base_ring()

        if is_RealField(base_ring):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = base_ring.precision()
            digits = sage.arith.all.ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(real,%d,0)"%digits, _vars, order=order)

        elif is_ComplexField(base_ring):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = base_ring.precision()
            digits = sage.arith.all.ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(complex,%d,0,I)"%digits, _vars,  order=order)

        elif base_ring.is_prime_field():
            self.__singular = singular.ring(self.characteristic(), _vars, order=order)

        elif base_ring.is_field() and base_ring.is_finite(): #must be extension field
            gen = str(base_ring.gen())
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order)
            self.__minpoly = (str(base_ring.modulus()).replace("x",gen)).replace(" ","")
            if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                singular.eval("minpoly=%s"%(self.__minpoly) )
                self.__minpoly = singular.eval('minpoly')[1:-1]
            self.__singular = r

        elif isinstance(base_ring, NumberField) and base_ring.is_absolute():
            gen = str(base_ring.gen())
            poly = base_ring.polynomial()
            poly_gen = str(poly.parent().gen())
            poly_str = str(poly).replace(poly_gen,gen)
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order, check=False)
            self.__minpoly = (poly_str).replace(" ","")
            if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                singular.eval("minpoly=%s"%(self.__minpoly) )
                self.__minpoly = singular.eval('minpoly')[1:-1]
            self.__singular = r

        elif is_IntegerRing(base_ring):
            self.__singular = singular.ring("(integer)", _vars, order=order)

        elif is_IntegerModRing(base_ring):
            ch = base_ring.characteristic()
            if ch.is_power_of(2):
                exp = ch.nbits() -1
                self.__singular = singular.ring("(integer,2,%d)"%(exp,), _vars, order=order, check=False)
            else:
                self.__singular = singular.ring("(integer,%d)"%(ch,), _vars, order=order, check=False)

        else:
            raise TypeError, "no conversion to a Singular ring defined"

        return self.__singular

    # It is required in cython to redefine __hash__ when __richcmp__ is
    # overloaded. Also just writing
    #         __hash__ = CategoryObject.__hash__
    # doesn't work.
    def __hash__(self):
        """
        Return a hash for this ring, that is, a hash of the string
        representation of this polynomial ring.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: hash(P)      # somewhat random output
            967902441410893180 # 64-bit
            -1767675994        # 32-bit
        """
        return CategoryObject.__hash__(self)

    def __richcmp__(left, right, int op):
        r"""
        Multivariate polynomial rings are said to be equal if:

        - their base rings match,
        - their generator names match and
        - their term orderings match.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: R.<x,y,z> = QQ[]
            sage: P == R
            True

            sage: R.<x,y,z> = GF(127)[]
            sage: P == R
            False

            sage: R.<x,y> = QQ[]
            sage: P == R
            False

            sage: R.<x,y,z> = PolynomialRing(QQ,order='invlex')
            sage: P == R
            False
        """
        return (<Parent>left)._richcmp(right, op)

    cpdef int _cmp_(left, right) except -2:
        """
        Compare ``left`` with ``right``.

        TEST::

            sage: R = QQ['x', 'y']; R
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: R == R
            True
            sage: R == QQ['x','z']
            False

        """
        if isinstance(right, (MPolynomialRing_libsingular, MPolynomialRing_polydict_domain)):
            return cmp((left.base_ring(), map(str, left.gens()), left.term_order()),
                       (right.base_ring(), map(str, right.gens()), right.term_order()))
        else:
            return cmp(type(left),type(right))

    def __reduce__(self):
        """
        Serializes self.

        EXAMPLES:
            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P.<x,y,z> = PolynomialRing(ZZ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(127), names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^8,'F'), names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^16,'B'), names='abc')
            sage: P == loads(dumps(P))
            True
            sage: z = QQ['z'].0
            sage: P = PolynomialRing(NumberField(z^2 + 3,'B'), names='abc')
            sage: P == loads(dumps(P))
            True
        """
        return sage.rings.polynomial.multi_polynomial_libsingular.unpickle_MPolynomialRing_libsingular, \
            (self.base_ring(), self.variable_names(), self.term_order())

    def __temporarily_change_names(self, names, latex_names):
        """
        This is used by the variable names context manager.

        EXAMPLES::

            sage: R.<x,y> = QQ[] # indirect doctest
            sage: with localvars(R, 'z,w'):
            ...       print x^3 + y^3 - x*y
            ...
            z^3 + w^3 - z*w
        """
        cdef ring *_ring = self._ring
        cdef char **_names
        cdef char **_orig_names
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
        r"""
        Return ``f/g``, where both ``f`` and`` ``g`` are treated as
        monomials.

        Coefficients are ignored by default.

        INPUT:

        - ``f`` - monomial
        - ``g`` - monomial
        - ``coeff`` - divide coefficients as well (default: ``False``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_quotient(3/2*x*y,x)
            y

            sage: P.monomial_quotient(3/2*x*y,x,coeff=True)
            3/2*y

        Note, that `\ZZ` behaves different if ``coeff=True``::

            sage: P.monomial_quotient(2*x,3*x)
            1

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: P.monomial_quotient(2*x,3*x,coeff=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot divide these coefficients.

        TESTS::

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

        .. warning::

           Assumes that the head term of f is a multiple of the head
           term of g and return the multiplicant m. If this rule is
           violated, funny things may happen.
        """
        cdef poly *res
        cdef ring *r = self._ring
        cdef number *n
        cdef number *denom

        if self is not f._parent:
            f = self._coerce_c(f)
        if self is not g._parent:
            g = self._coerce_c(g)

        if not f._poly:
            return self._zero_element
        if not g._poly:
            raise ZeroDivisionError

        if r!=currRing: rChangeCurrRing(r)  # pDivide
        res = pDivide(f._poly, g._poly)
        if coeff:
            if r.ringtype == 0 or r.cf.nDivBy(p_GetCoeff(f._poly, r), p_GetCoeff(g._poly, r)):
                n = r.cf.nDiv( p_GetCoeff(f._poly, r) , p_GetCoeff(g._poly, r))
                p_SetCoeff0(res, n, r)
            else:
                raise ArithmeticError("Cannot divide these coefficients.")
        else:
            p_SetCoeff0(res, n_Init(1, r), r)
        return new_MP(self, res)

    def monomial_divides(self, MPolynomial_libsingular a, MPolynomial_libsingular b):
        """
        Return ``False`` if a does not divide b and ``True``
        otherwise.

        Coefficients are ignored.

        INPUT:

        - ``a`` -- monomial

        - ``b`` -- monomial

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_divides(x*y*z, x^3*y^2*z^4)
            True
            sage: P.monomial_divides(x^3*y^2*z^4, x*y*z)
            False

        TESTS::

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
            b = a._parent._coerce_c(b)

        _a = a._poly
        _b = b._poly
        _r = a._parent_ring

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

        - ``f`` - monomial

        - ``g`` - monomial

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_lcm(3/2*x*y,x)
            x*y

        TESTS::

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

        if self is not f._parent:
            f = self._coerce_c(f)
        if self is not g._parent:
            g = self._coerce_c(g)

        if f._poly == NULL:
            if g._poly == NULL:
                return self._zero_element
            else:
                raise ArithmeticError, "Cannot compute LCM of zero and nonzero element."
        if g._poly == NULL:
            raise ArithmeticError, "Cannot compute LCM of zero and nonzero element."

        if(self._ring != currRing): rChangeCurrRing(self._ring)

        pLcm(f._poly, g._poly, m)
        p_Setm(m, self._ring)
        return new_MP(self,m)

    def monomial_reduce(self, MPolynomial_libsingular f, G):
        """
        Try to find a ``g`` in ``G`` where ``g.lm()`` divides
        ``f``. If found ``(flt,g)`` is returned, ``(0,0)`` otherwise,
        where ``flt`` is ``f/g.lm()``.

        It is assumed that ``G`` is iterable and contains *only*
        elements in this polynomial ring.

        Coefficients are ignored.

        INPUT:

        - ``f`` - monomial
        - ``G`` - list/set of mpolynomials

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, 1/2  ]
            sage: P.monomial_reduce(f,G)
            (y, 1/4*x*y + 2/7)

        TESTS::

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
            if isinstance(g, MPolynomial_libsingular) \
                   and (<MPolynomial_libsingular>g) \
                   and g.parent() is self \
                   and p_LmDivisibleBy((<MPolynomial_libsingular>g)._poly, m, r):
                if r!=currRing: rChangeCurrRing(r)  # pDivide
                flt = pDivide(f._poly, (<MPolynomial_libsingular>g)._poly)
                #p_SetCoeff(flt, n_Div( p_GetCoeff(f._poly, r) , p_GetCoeff((<MPolynomial_libsingular>g)._poly, r), r), r)
                p_SetCoeff(flt, n_Init(1, r), r)
                return new_MP(self,flt), g
        return self._zero_element,self._zero_element

    def monomial_pairwise_prime(self, MPolynomial_libsingular g, MPolynomial_libsingular h):
        """
        Return ``True`` if ``h`` and ``g`` are pairwise prime. Both
        are treated as monomials.

        Coefficients are ignored.

        INPUT:

        - ``h`` - monomial
        - ``g`` - monomial

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_pairwise_prime(x^2*z^3, y^4)
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, 3/4*y^3)
            False

        TESTS::

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
        cdef poly *p
        cdef poly *q

        if h._parent is not g._parent:
            g = h._parent._coerce_c(g)

        r = h._parent_ring
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
        Return a list of all monomials that divide ``t``.

        Coefficients are ignored.

        INPUT:

        - ``t`` - a monomial

        OUTPUT:
            a list of monomials


        EXAMPLES::

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
          M.append(new_MP(self, p_Copy(tempvector, _ring)))
        return M

def unpickle_MPolynomialRing_libsingular(base_ring, names, term_order):
    """
    inverse function for ``MPolynomialRing_libsingular.__reduce__``

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(P)) is P # indirect doctest
        True
    """
    from sage.rings.polynomial.polynomial_ring_constructor import _multi_variate
    # If libsingular would be replaced by a different implementation in future
    # sage version, the unpickled ring will belong the new implementation.
    return _multi_variate(base_ring, tuple(names), len(names), False, term_order, None)

cdef class MPolynomial_libsingular(sage.rings.polynomial.multi_polynomial.MPolynomial):
    """
    A multivariate polynomial implemented using libSINGULAR.
    """
    def __init__(self, MPolynomialRing_libsingular parent):
        """
        Construct a zero element in parent.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomial_libsingular
            sage: P = PolynomialRing(GF(32003),3,'x')
            sage: MPolynomial_libsingular(P)
            0
        """
        self._poly = NULL
        self._parent = parent
        self._parent_ring = singular_ring_reference(parent._ring)

    def __dealloc__(self):
        # WARNING: the Cython class self._parent is now no longer accessible!
        if self._poly==NULL:
            # e.g. MPolynomialRing_libsingular._zero_element
            singular_ring_delete(self._parent_ring)
            return
        assert self._parent_ring != NULL # the constructor has no way to raise an exception
        p_Delete(&self._poly, self._parent_ring)
        singular_ring_delete(self._parent_ring)

    def __copy__(self):
        """
        Copy ``self``.

        OUTPUT:

        A copy.

        EXAMPLES::

            sage: F.<a> = GF(7^2)
            sage: R.<x,y> = F[]
            sage: p = a*x^2 + y + a^3; p
            (a)*x^2 + y + (-2*a - 3)
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
            sage: lst = [p,q];
            sage: matrix(ZZ, 2, 2, lambda i,j: bool(lst[i]==lst[j]))
            [1 1]
            [1 1]
            sage: matrix(ZZ, 2, 2, lambda i,j: bool(lst[i] is lst[j]))
            [1 0]
            [0 1]
        """
        return new_MP(self._parent, p_Copy(self._poly, self._parent_ring))

    def __deepcopy__(self, memo={}):
        """
        Deep copy ``self``

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: p = x^2 + y^2
            sage: p is deepcopy(p)
            False
            sage: p == deepcopy(p)
            True
            sage: p.parent() is deepcopy(p).parent()
            True
        """
        cpy = self.__copy__()
        memo[id(self)] = cpy
        return cpy

    cpdef MPolynomial_libsingular _new_constant_poly(self, x, MPolynomialRing_libsingular P):
        r"""
        Quickly create a new constant polynomial with value x in the parent P.

        ASSUMPTION:

        The value x must be an element of the base ring. That assumption is
        not verified.

        EXAMPLE::

            sage: R.<x,y> = QQ[]
            sage: x._new_constant_poly(2/1,R)
            2

        """
        if not x:
            return new_MP(P, NULL)
        cdef poly *_p
        singular_polynomial_rmul(&_p, P._one_element_poly, x, P._ring)
        return new_MP(P, _p)

    def __call__(self, *x, **kwds):
        """
        Evaluate this multi-variate polynomial at ``x``, where ``x``
        is either the tuple of values to substitute in, or one can use
        functional notation ``f(a_0,a_1,a_2, \ldots)`` to evaluate
        ``f`` with the ith variable replaced by ``a_i``.

        INPUT:

        - ``x`` - a list of elements in ``self.parent()``
        - or ``**kwds`` - a dictionary of ``variable-name:value`` pairs.

        EXAMPLES::

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

        TESTS::

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

        See :trac:`8502`::

            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^2+47)
            sage: R.<X,Y,Z> = K[]
            sage: f = X+Y+Z
            sage: a = f(t,t,t); a
            3*t
            sage: a.parent() is K
            True

            sage: R.<X,Y,Z> = QQ[]
            sage: f = X+Y+Z
            sage: a = f(2,3,4); a
            9
            sage: a.parent() is QQ
            True
        """
        if len(kwds) > 0:
            f = self.subs(**kwds)
            if len(x) > 0:
                return f(*x)
            else:
                return f

        cdef int l = len(x)
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *_ring = parent._ring

        if l == 1 and isinstance(x[0], (list, tuple)):
            x = x[0]
            l = len(x)

        if l != parent._ring.N:
            raise TypeError, "number of arguments does not match number of variables in parent"

        try:
            # Attempt evaluation via singular.
            coerced_x = [parent.coerce(e) for e in x]
        except TypeError:
            # give up, evaluate functional
            y = parent.base_ring()(0)
            for (m,c) in self.dict().iteritems():
                y += c*mul([ x[i]**m[i] for i in m.nonzero_positions()])
            return y

        cdef poly *res    # ownership will be transferred to us in the next line
        singular_polynomial_call(&res, self._poly, _ring, coerced_x, MPolynomial_libsingular_get_element)
        res_parent = coercion_model.common_parent(parent._base, *x)

        if res == NULL:
            return res_parent(0)
        if p_LmIsConstant(res, _ring):
            sage_res = si2sa( p_GetCoeff(res, _ring), _ring, parent._base )
            p_Delete(&res, _ring)            # sage_res contains copy
        else:
            sage_res = new_MP(parent, res)   # pass on ownership of res to sage_res

        if parent(sage_res) is not res_parent:
            sage_res = res_parent(sage_res)
        return sage_res

    def __hash__(self):
        """
        This hash incorporates the variable name in an effort to
        respect the obvious inclusions into multi-variable polynomial
        rings.

        The tuple algorithm is borrowed from http://effbot.org/zone/python-hash.htm.

        EXAMPLES::

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
        return self._hash_c()

    cpdef int _cmp_(left, Element right) except -2:
        """
        Compare left and right and return -1, 0, and 1 for <,==, and >
        respectively.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,order='degrevlex')
            sage: x == x
            True

            sage: x > y
            True
            sage: y^2 > x
            True

            sage: (2/3*x^2 + 1/2*y + 3) > (2/3*x^2 + 1/4*y + 10)
            True

        TESTS::

            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
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

            sage: P.<x,y,z> = PolynomialRing(GF(127), order='degrevlex')
            sage: (66*x^2 + 23) > (66*x^2 + 2)
            True
        """
        if left is right:
            return 0
        cdef poly *p = (<MPolynomial_libsingular>left)._poly
        cdef poly *q = (<MPolynomial_libsingular>right)._poly
        cdef ring *r = (<MPolynomial_libsingular>left)._parent_ring
        return singular_polynomial_cmp(p, q, r)

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        Add left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x + 1/2*y + 1 #indirect doctest
            3/2*x + 1/2*y + 1
        """
        cdef ring *r = (<MPolynomial_libsingular>left)._parent_ring
        cdef poly *_p
        singular_polynomial_add(&_p, left._poly,
                                 (<MPolynomial_libsingular>right)._poly, r)
        return new_MP((<MPolynomial_libsingular>left)._parent, _p)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        """
        Subtract left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x - 1/2*y - 1 #indirect doctest
            3/2*x - 1/2*y - 1
        """
        cdef ring *_ring = (<MPolynomial_libsingular>left)._parent_ring
        cdef poly *_p
        singular_polynomial_sub(&_p, left._poly,
                                (<MPolynomial_libsingular>right)._poly,
                                _ring)
        return new_MP((<MPolynomial_libsingular>left)._parent, _p)

    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        Multiply self with a base ring element.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x # indirect doctest
            3/2*x
        """

        cdef ring *_ring = (<MPolynomial_libsingular>self)._parent_ring
        if not left:
            return (<MPolynomial_libsingular>self)._parent._zero_element
        cdef poly *_p
        singular_polynomial_rmul(&_p, self._poly, left, _ring)
        return new_MP((<MPolynomial_libsingular>self)._parent, _p)

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2) # indirect doctest
            9/4*x - 3/4*y - 3/2
        """
        # Note that the argument to _rmul_ and _lmul_ is an
        # element of the base ring.
        # All currently implemented base rings are commutative,
        # So, calling _rmul_ is the correct thing to do.
        return self._rmul_(right)

    cpdef RingElement _mul_(left, RingElement right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2*x + 1/2*y + 1) # indirect doctest
            9/4*x^2 - 1/4*y^2 - y - 1

            sage: P.<x,y> = PolynomialRing(QQ,order='lex')
            sage: (x^2^30) * x^2^30
            Traceback (most recent call last):
            ...
            OverflowError: Exponent overflow (...).
        """
        # all currently implemented rings are commutative
        cdef poly *_p
        singular_polynomial_mul(&_p, left._poly,
                                 (<MPolynomial_libsingular>right)._poly,
                                 (<MPolynomial_libsingular>left)._parent_ring)
        return new_MP((<MPolynomial_libsingular>left)._parent,_p)

    cpdef RingElement _div_(left, RingElement right_ringelement):
        """
        Divide left by right

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = (x + y)/3 # indirect doctest
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        Note that / is still a constructor for elements of the
        fraction field in all cases as long as both arguments have the
        same parent and right is not constant. ::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: g = x
            sage: h = f/g; h
            (x^3 + y)/x
            sage: h.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field

        If we divide over `\ZZ` the result is the same as multiplying
        by 1/3 (i.e. base extension). ::

            sage: R.<x,y> = ZZ[]
            sage: f = (x + y)/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: f = (x + y) * 1/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        But we get a true fraction field if the denominator is not in
        the fraction field of the base ring.""

            sage: f = x/y
            sage: f.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring

        Division will fail for non-integral domains::

            sage: P.<x,y> = Zmod(1024)[]
            sage: x/P(3)
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

            sage: x/3
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

        TESTS::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        cdef poly *p
        cdef MPolynomial_libsingular right = <MPolynomial_libsingular>right_ringelement
        cdef bint is_field = left._parent._base.is_field()
        if p_IsConstant(right._poly, right._parent_ring):
            if is_field:
                singular_polynomial_div_coeff(&p, left._poly, right._poly, right._parent_ring)
                return new_MP(left._parent, p)
            else:
                return left.change_ring(left.base_ring().fraction_field())/right_ringelement
        else:
            return (left._parent).fraction_field()(left,right_ringelement)

    def __pow__(MPolynomial_libsingular self, exp, ignored):
        """
        Return ``self**(exp)``.

        The exponent must be an integer.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f^2
            x^6 + 2*x^3*y + y^2
            sage: g = f^(-1); g
            1/(x^3 + y)
            sage: type(g)
            <type 'sage.rings.fraction_field_element.FractionFieldElement'>

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: P(2)**(-1)
            1/2

            sage: P.<u,v> = PolynomialRing(QQ, 2)
            sage: u^(1/2)
            Traceback (most recent call last):
            ...
            TypeError: non-integral exponents not supported

            sage: P.<x,y> = PolynomialRing(QQ,order='lex')
            sage: (x+y^2^30)^10
            Traceback (most recent call last):
            ....
            OverflowError: Exponent overflow (...).
        """
        if type(exp) is not Integer:
            try:
                exp = Integer(exp)
            except TypeError:
                raise TypeError("non-integral exponents not supported")

        if exp < 0:
            return 1/(self**(-exp))
        elif exp == 0:
            return self._parent._one_element

        cdef ring *_ring = self._parent_ring
        cdef poly *_p
        singular_polynomial_pow(&_p, self._poly, exp, _ring)
        return new_MP(self._parent, _p)

    def __neg__(self):
        """
        Return ``-self``.

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: -f
            -x^3 - y
        """
        cdef ring *_ring = self._parent_ring

        cdef poly *p
        singular_polynomial_neg(&p, self._poly, _ring)
        return new_MP(self._parent, p)

    def _repr_(self):
        """
        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f # indirect doctest
            x^3 + y
        """
        cdef ring *_ring = self._parent_ring
        s = singular_polynomial_str(self._poly, _ring)
        return s

    cpdef _repr_short_(self):
        """
        This is a faster but less pretty way to print polynomials. If
        available it uses the short SINGULAR notation.

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f._repr_short_()
            'x3+y'
        """
        cdef ring *_ring = self._parent_ring
        rChangeCurrRing(_ring)
        if _ring.CanShortOut:
            _ring.ShortOut = 1
            s = p_String(self._poly, _ring, _ring)
            _ring.ShortOut = 0
        else:
            s = p_String(self._poly, _ring, _ring)
        return s

    def _latex_(self):
        """
        Return a polynomial LaTeX representation of this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: latex(f)  # indirect doctest
            - x^{2} y - \frac{25}{27} y^{3} - z^{2}
        """
        cdef ring *_ring = self._parent_ring
        gens = self.parent().latex_variable_names()
        base = self.parent().base()
        return singular_polynomial_latex(self._poly, _ring, base, gens)

    def _repr_with_changed_varnames(self, varnames):
        """
        Return string representing this polynomial but change the
        variable names to ``varnames``.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: print f._repr_with_changed_varnames(['FOO', 'BAR', 'FOOBAR'])
            -FOO^2*BAR - 25/27*BAR^3 - FOOBAR^2
        """
        return  singular_polynomial_str_with_changed_varnames(self._poly, self._parent_ring, varnames)

    def degree(self, MPolynomial_libsingular x=None, int std_grading=False):
        """
        Return the maximal degree of this polynomial in ``x``, where
        ``x`` must be one of the generators for the parent of this
        polynomial.

        INPUT:

        - ``x`` - (default: ``None``) a multivariate polynomial which is (or
          coerces to) a generator of the parent of self. If ``x`` is ``None``,
          return the total degree, which is the maximum degree of any monomial.
          Note that a matrix term ordering alters the grading of the generators
          of the ring; see the tests below.  To avoid this behavior, use either
          ``exponents()`` for the exponents themselves, or the optional
          argument ``std_grading=False``.

        OUTPUT:
            integer

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = y^2 - x^9 - x
            sage: f.degree(x)
            9
            sage: f.degree(y)
            2
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(x)
            3
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(y)
            10

        TESTS::

            sage: P.<x, y> = QQ[]
            sage: P(0).degree(x)
            -1
            sage: P(1).degree(x)
            0

        With a matrix term ordering, the grading of the generators is
        determined by the first row of the matrix.  This affects the behavior
        of ``degree()`` when no variable is specified.
        To evaluate the degree with a standard grading, use the optional
        argument ``std_grading=True``.

            sage: tord = TermOrder(matrix([3,0,1,1,1,0,1,0,0]))
            sage: R.<x,y,z> = PolynomialRing(QQ,'x',3,order=tord)
            sage: (x^3*y+x*z^4).degree()
            9
            sage: (x^3*y+x*z^4).degree(std_grading=True)
            5
            sage: x.degree(x), y.degree(y), z.degree(z)
            (1, 1, 1)

        The following example is inspired by :trac:`11652`::

            sage: R.<p,q,t> = ZZ[]
            sage: poly = p+q^2+t^3
            sage: poly = poly.polynomial(t)[0]
            sage: poly
            q^2 + p

        There is no canonical coercion from ``R`` to the parent of ``poly``, so
        this doesn't work::

            sage: poly.degree(q)
            Traceback (most recent call last):
            ...
            TypeError: argument must canonically coerce to parent

        Using a non-canonical coercion does work, but we require this
        to be done explicitly, since it can lead to confusing results
        if done automatically::

            sage: poly.degree(poly.parent()(q))
            2
            sage: poly.degree(poly.parent()(p))
            1
            sage: T.<x,y> = ZZ[]
            sage: poly.degree(poly.parent()(x))  # noncanonical coercions can be confusing
            1

        The argument to degree has to be a generator::

            sage: pp = poly.parent().gen(0)
            sage: poly.degree(pp)
            1
            sage: poly.degree(pp+1)
            Traceback (most recent call last):
            ...
            TypeError: argument must be a generator

        Canonical coercions are used::

            sage: S = ZZ['p,q']
            sage: poly.degree(S.0)
            1
            sage: poly.degree(S.1)
            2
        """
        cdef ring *r = self._parent_ring
        cdef poly *p = self._poly
        if not x:
            if not std_grading:
                return singular_polynomial_deg(p,NULL,r)
            else:
                return self.total_degree(std_grading=True)

        if not x.parent() is self.parent():
            try:
                x = self.parent().coerce(x)
            except TypeError:
                raise TypeError("argument must canonically coerce to parent")
        if not x.is_generator():
            raise TypeError("argument must be a generator")

        return singular_polynomial_deg(p, x._poly, r)

    def total_degree(self, int std_grading=False):
        """
        Return the total degree of ``self``, which is the maximum degree
        of all monomials in ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
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

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: R(0).total_degree()
            -1
            sage: R(1).total_degree()
            0

        With a matrix term ordering, the grading changes.
        To evaluate the total degree using the standard grading,
        use the optional argument``std_grading=True``.

            sage: tord=TermOrder(matrix([3,0,1,1,1,0,1,0,0]))
            sage: R.<x,y,z> = PolynomialRing(QQ,'x',3,order=tord)
            sage: (x^2*y).total_degree()
            6
            sage: (x^2*y).total_degree(std_grading=True)
            3
        """
        cdef int i, result
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring

        if not std_grading:
            return singular_polynomial_deg(p,NULL,r)
        else:
            result = 0
            while p:
                result = max(result, sum([p_GetExp(p,i,r) for i in xrange(1,r.N+1)]))
                p = pNext(p)
            return result

    def degrees(self):
        """
        Returns a tuple with the maximal degree of each variable in
        this polynomial.  The list of degrees is ordered by the order
        of the generators.

        EXAMPLES::

            sage: R.<y0,y1,y2> = PolynomialRing(QQ,3)
            sage: q = 3*y0*y1*y1*y2; q
            3*y0*y1^2*y2
            sage: q.degrees()
            (1, 2, 1)
            sage: (q + y0^5).degrees()
            (5, 2, 1)
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int i
        cdef list d = [0 for _ in range(r.N)]
        while p:
            for i from 0 <= i < r.N:
                d[i] = max(d[i],p_GetExp(p, i+1, r))
            p = pNext(p)
        return tuple(d)

    def coefficient(self, degrees):
        """
        Return the coefficient of the variables with the degrees
        specified in the python dictionary ``degrees``.
        Mathematically, this is the coefficient in the base ring
        adjoined by the variables of this ring not listed in
        ``degrees``.  However, the result has the same parent as this
        polynomial.

        This function contrasts with the function
        ``monomial_coefficient`` which returns the coefficient in the
        base ring of a monomial.

        INPUT:

        - ``degrees`` - Can be any of:
                - a dictionary of degree restrictions
                - a list of degree restrictions (with None in the unrestricted variables)
                - a monomial (very fast, but not as flexible)

        OUTPUT:
            element of the parent of this element.

        .. note::

           For coefficients of specific monomials, look at :meth:`monomial_coefficient`.

        EXAMPLES::

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


        Be aware that this may not be what you think! The physical
        appearance of the variable x is deceiving -- particularly if
        the exponent would be a variable. ::

            sage: f.coefficient(x^0) # outputs the full polynomial
            x^2*y^2 + x^2*y + x*y^2 + x^2 + x*y + y^2 + x + y + 1
            sage: R.<x,y> = GF(389)[]
            sage: f=x*y+5
            sage: c=f.coefficient({x:0,y:0}); c
            5
            sage: parent(c)
            Multivariate Polynomial Ring in x, y over Finite Field of size 389

        AUTHOR:

        - Joel B. Mohler (2007.10.31)
        """
        cdef poly *_degrees = <poly*>0
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef poly *newp = p_ISet(0,r)
        cdef poly *newptemp
        cdef int i
        cdef int flag
        cdef int gens = self._parent.ngens()
        cdef int *exps = <int*>sig_malloc(sizeof(int)*gens)
        for i from 0<=i<gens:
            exps[i] = -1

        if isinstance(degrees, MPolynomial_libsingular) \
                and self._parent is (<MPolynomial_libsingular>degrees)._parent:
            _degrees = (<MPolynomial_libsingular>degrees)._poly
            if pLength(_degrees) != 1:
                raise TypeError, "degrees must be a monomial"
            for i from 0<=i<gens:
                if p_GetExp(_degrees,i+1,r)!=0:
                    exps[i] = p_GetExp(_degrees,i+1,r)
        elif isinstance(degrees, list):
            for i from 0<=i<gens:
                if degrees[i] is None:
                    exps[i] = -1
                else:
                    exps[i] = int(degrees[i])
        elif isinstance(degrees, dict):
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

        sig_free(exps)
        return new_MP(self._parent, newp)

    def monomial_coefficient(self, MPolynomial_libsingular mon):
        """
        Return the coefficient in the base ring of the monomial mon in
        ``self``, where mon must have the same parent as self.

        This function contrasts with the function ``coefficient``
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` - a monomial

        OUTPUT:
            coefficient in base ring

        SEE ALSO:
            For coefficients in a base ring of fewer variables, look at ``coefficient``.

        EXAMPLES::

            sage: P.<x,y> = QQ[]

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
        cdef ring *r = self._parent_ring

        if not mon._parent is self._parent:
            raise TypeError("mon must have same parent as self.")

        while(p):
            if p_ExpVectorEqual(p, m, r) == 1:
                return si2sa(p_GetCoeff(p, r), r, self._parent._base)
            p = pNext(p)

        return self._parent._base._zero_element

    def dict(self):
        """
        Return a dictionary representing self. This dictionary is in
        the same format as the generic MPolynomial: The dictionary
        consists of ``ETuple:coefficient`` pairs.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f=2*x*y^3*z^2 + 1/7*x^2 + 2/3
            sage: f.dict()
            {(0, 0, 0): 2/3, (1, 3, 2): 2, (2, 0, 0): 1/7}
        """
        cdef poly *p
        cdef ring *r = self._parent_ring
        cdef int n
        cdef int v
        if r!=currRing: rChangeCurrRing(r)
        base = self._parent._base
        p = self._poly
        pd = dict()
        while p:
            d = dict()
            for v from 1 <= v <= r.N:
                n = p_GetExp(p,v,r)
                if n!=0:
                    d[v-1] = n

            pd[ETuple(d,r.N)] = si2sa(p_GetCoeff(p, r), r, base)

            p = pNext(p)
        return pd

    cdef long _hash_c(self) except -1:
        """
        See ``self.__hash__``
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int n
        cdef int v
        if r!=currRing: rChangeCurrRing(r)
        base = self._parent._base
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        var_name_hash = [hash(vn) for vn in self._parent.variable_names()]
        cdef long c_hash
        while p:
            c_hash = hash(si2sa(p_GetCoeff(p, r), r, base))
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

    def __getitem__(self,x):
        """
        Same as ``self.monomial_coefficient`` but for exponent vectors.

        INPUT:

        - ``x`` - a tuple or, in case of a single-variable MPolynomial
          ring x can also be an integer.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = -10*x^3*y + 17*x*y
            sage: f[3,1]
            -10
            sage: f[1,1]
            17
            sage: f[0,1]
            0

            sage: R.<x> = PolynomialRing(GF(7),1); R
            Multivariate Polynomial Ring in x over Finite Field of size 7
            sage: f = 5*x^2 + 3; f
            -2*x^2 + 3
            sage: f[2]
            5
        """
        cdef poly *m
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int i

        if isinstance(x, MPolynomial_libsingular):
            return self.monomial_coefficient(x)
        if not isinstance(x, tuple):
            try:
                x = tuple(x)
            except TypeError:
                x = (x,)

        if len(x) != self._parent.ngens():
            raise TypeError, "x must have length self.ngens()"

        m = p_ISet(1,r)
        i = 1
        for e in x:
            overflow_check(e, r)
            p_SetExp(m, i, int(e), r)
            i += 1
        p_Setm(m, r)

        while(p):
            if p_ExpVectorEqual(p, m, r) == 1:
                p_Delete(&m,r)
                return si2sa(p_GetCoeff(p, r), r, self._parent._base)
            p = pNext(p)

        p_Delete(&m,r)
        return self._parent._base._zero_element

    def exponents(self, as_ETuples=True):
        """
        Return the exponents of the monomials appearing in this
        polynomial.

        INPUT:

        - ``as_ETuples`` - (default: ``True``) if true returns the result as an list of ETuples
                          otherwise returns a list of tuples


        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: f = a^3 + b + 2*b^2
            sage: f.exponents()
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
            sage: f.exponents(as_ETuples=False)
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int v
        cdef list pl, ml

        pl = list()
        ml = range(r.N)
        while p:
            for v from 1 <= v <= r.N:
                ml[v-1] = p_GetExp(p,v,r)

            if as_ETuples:
                pl.append(ETuple(ml))
            else:
                pl.append(tuple(ml))

            p = pNext(p)
        return pl

    def is_unit(self):
        """
        Return ``True`` if self is a unit.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
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

            sage: R.<x,y> = ZZ[]
            sage: R(1).is_unit()
            True
            sage: R(2).is_unit()
            False

        """
        cdef bint is_field = self._parent._base.is_field()
        if is_field:
            if self._parent_ring != currRing: rChangeCurrRing(self._parent_ring) # bug in p_IsUnit
            return bool(p_IsUnit(self._poly, self._parent_ring))
        else:
            return bool(p_IsConstant(self._poly, self._parent_ring) and self.constant_coefficient().is_unit())

    def inverse_of_unit(self):
        """
        Return the inverse of this polynomial if it is a unit.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: x.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: Element is not a unit.

            sage: R(1/2).inverse_of_unit()
            2
        """
        cdef ring *_ring = self._parent_ring
        if(_ring != currRing): rChangeCurrRing(_ring)

        if not p_IsUnit(self._poly, _ring):
            raise ArithmeticError, "Element is not a unit."

        sig_on()
        cdef MPolynomial_libsingular r = new_MP(self._parent, pInvers(0, self._poly, NULL))
        sig_off()
        return r

    def is_homogeneous(self):
        """
        Return ``True`` if this polynomial is homogeneous.

        EXAMPLES::

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
        cdef ring *_ring = self._parent_ring
        if(_ring != currRing): rChangeCurrRing(_ring)
        return bool(pIsHomogeneous(self._poly))

    cpdef _homogenize(self, int var):
        """
        Return ``self`` if ``self`` is homogeneous.  Otherwise return
        a homogenized polynomial constructed by modifying the degree
        of the variable with index ``var``.

        INPUT:

        - ``var`` - an integer indicating which variable to use to
          homogenize (``0 <= var < parent(self).ngens()``)

        OUTPUT:
            a multivariate polynomial

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: f = x^2 + y + 1 + 5*x*y^10
            sage: g = f.homogenize('z'); g # indirect doctest
            5*x*y^10 + x^2*z^9 + y*z^10 + z^11
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: f._homogenize(0)
            2*x^11 + x^10*y + 5*x*y^10

        SEE: ``self.homogenize``
        """
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef MPolynomial_libsingular f

        if self.is_homogeneous():
            return self

        if var < parent._ring.N:
            return new_MP(parent, pHomogen(p_Copy(self._poly, self._parent_ring), var+1))
        else:
            raise TypeError, "var must be < self.parent().ngens()"

    def is_monomial(self):
        """
        Return ``True`` if this polynomial is a monomial.  A monomial
        is defined to be a product of generators with coefficient 1.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: x.is_monomial()
            True
            sage: (2*x).is_monomial()
            False
            sage: (x*y).is_monomial()
            True
            sage: (x*y + x).is_monomial()
            False
        """
        cdef poly *_p
        cdef ring *_ring
        cdef number *_n
        _ring = self._parent_ring

        if self._poly == NULL:
            return True

        if(_ring != currRing): rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        ret = bool((not self._poly.next) and _ring.cf.nIsOne(_n))

        p_Delete(&_p, _ring)
        return ret

    def subs(self, fixed=None, **kw):
        """
        Fixes some given variables in a given multivariate polynomial
        and returns the changed multivariate polynomials. The
        polynomial itself is not affected.  The variable,value pairs
        for fixing are to be provided as dictionary of the form
        ``{variable:value}``.

        This is a special case of evaluating the polynomial with some
        of the variables constants and the others the original
        variables, but should be much faster if only few variables are
        to be fixed.

        INPUT:

        - ``fixed`` - (optional) dict with variable:value pairs
        - ``**kw`` - names parameters

        OUTPUT:
            a new multivariate polynomial

        EXAMPLES::

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

        The parameters are subsituted in order and without side effects::

            sage: R.<x,y>=QQ[]
            sage: g=x+y
            sage: g.subs({x:x+1,y:x*y})
            x*y + x + 1
            sage: g.subs({x:x+1}).subs({y:x*y})
            x*y + x + 1
            sage: g.subs({y:x*y}).subs({x:x+1})
            x*y + x + y + 1

        ::

            sage: R.<x,y> = QQ[]
            sage: f = x + 2*y
            sage: f.subs(x=y,y=x)
            2*x + y

        TESTS::

            sage: P.<x,y,z> = QQ[]
            sage: f = y
            sage: f.subs({y:x}).subs({x:z})
            z

        We test that we change the ring even if there is nothing to do::

            sage: P = QQ['x,y']
            sage: x = var('x')
            sage: parent(P.zero() / x)
            Symbolic Ring

        We are catching overflows::

            sage: R.<x,y> = QQ[]
            sage: n=1000; f = x^n
            sage: try:
            ....:   f.subs(x = x^n)
            ....:   print "no overflow"
            ....: except OverflowError:
            ....:   print "overflow"
            overflow    # 32-bit
            x^1000000   # 64-bit
            no overflow # 64-bit

            sage: n=100000;
            sage: try:
            ....:   f = x^n
            ....:   f.subs(x = x^n)
            ....:   print "no overflow"
            ....: except OverflowError:
            ....:   print "overflow"
            overflow

        Check that there is no more segmentation fault if the polynomial gets 0
        in the middle of a substitution (:trac:`17785`)::

            sage: R.<x,y,z> = QQ[]
            sage: for vx in [0,x,y,z]:
            ....:     for vy in [0,x,y,z]:
            ....:         for vz in [0,x,y,z]:
            ....:             d = {x:vx, y:vy, z:vz}
            ....:             ds = {'x': vx, 'y': vy, 'z': vz}
            ....:             assert x.subs(d) == x.subs(**ds) == vx
            ....:             assert y.subs(d) == y.subs(**ds) == vy
            ....:             assert z.subs(d) == z.subs(**ds) == vz
            ....:             assert (x+y).subs(d) == (x+y).subs(**ds) == vx+vy
        """
        cdef int mi, i, need_map, try_symbolic

        cdef unsigned long degree = 0
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *_ring = parent._ring

        if(_ring != currRing): rChangeCurrRing(_ring)

        cdef poly *_p = p_Copy(self._poly, _ring)
        cdef poly *_f

        cdef ideal *to_id = idInit(_ring.N,1)
        cdef ideal *from_id
        cdef ideal *res_id
        need_map = 0
        try_symbolic = 0

        if _p == NULL:
            # the polynomial is 0. There is nothing to do except to change the
            # ring
            try_symbolic = 1

        if not try_symbolic and fixed is not None:
            for m,v in fixed.iteritems():
                if isinstance(m, (int, Integer)):
                    mi = m+1
                elif isinstance(m,MPolynomial_libsingular) and m.parent() is parent:
                    for i from 0 < i <= _ring.N:
                        if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                            mi = i
                            break
                    if i > _ring.N:
                        id_Delete(&to_id, _ring)
                        p_Delete(&_p, _ring)
                        raise TypeError("key does not match")
                else:
                    id_Delete(&to_id, _ring)
                    p_Delete(&_p, _ring)
                    raise TypeError("keys do not match self's parent")
                try:
                    v = parent._coerce_c(v)
                except TypeError:
                    try_symbolic = 1
                    break
                _f = (<MPolynomial_libsingular>v)._poly
                if p_IsConstant(_f, _ring):
                    singular_polynomial_subst(&_p, mi-1, _f, _ring)
                else:
                    need_map = 1
                    degree = <unsigned long>p_GetExp(_p, mi, _ring) * <unsigned long>p_GetMaxExp(_f, _ring)
                    if  degree > _ring.bitmask:
                        id_Delete(&to_id, _ring)
                        p_Delete(&_p, _ring)
                        raise OverflowError("Exponent overflow (%d)."%(degree))
                    to_id.m[mi-1] = p_Copy(_f, _ring)

                if _p == NULL:
                    # polynomial becomes 0 after some substitution
                    try_symbolic = 1
                    break

        cdef dict gd

        if not try_symbolic:
            gd = parent.gens_dict()
            for m,v in kw.iteritems():
                m = gd[m]
                for i from 0 < i <= _ring.N:
                    if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                        mi = i
                        break
                if i > _ring.N:
                    id_Delete(&to_id, _ring)
                    p_Delete(&_p, _ring)
                    raise TypeError("key does not match")
                try:
                    v = parent._coerce_c(v)
                except TypeError:
                    try_symbolic = 1
                    break
                _f = (<MPolynomial_libsingular>v)._poly
                if p_IsConstant(_f, _ring):
                    singular_polynomial_subst(&_p, mi-1, _f, _ring)
                else:
                    if to_id.m[mi-1] != NULL:
                        p_Delete(&to_id.m[mi-1],_ring)
                    to_id.m[mi-1] = p_Copy(_f, _ring)
                    degree = <unsigned long>p_GetExp(_p, mi, _ring) * <unsigned long>p_GetMaxExp(_f, _ring)
                    if degree > _ring.bitmask:
                        id_Delete(&to_id, _ring)
                        p_Delete(&_p, _ring)
                        raise OverflowError("Exponent overflow (%d)."%(degree))
                    need_map = 1

                if _p == NULL:
                    # the polynomial is 0
                    try_symbolic = 1
                    break

            if need_map:
                for mi from 0 <= mi < _ring.N:
                    if to_id.m[mi] == NULL:
                        to_id.m[mi] = p_ISet(1,_ring)
                        p_SetExp(to_id.m[mi], mi+1, 1, _ring)
                        p_Setm(to_id.m[mi], _ring)

                from_id=idInit(1,1)
                from_id.m[0] = _p

                rChangeCurrRing(_ring)
                res_id = fast_map(from_id, _ring, to_id, _ring)
                _p = res_id.m[0]

                from_id.m[0] = NULL
                res_id.m[0] = NULL

                id_Delete(&from_id, _ring)
                id_Delete(&res_id, _ring)

        id_Delete(&to_id, _ring)

        if not try_symbolic:
            return new_MP(parent,_p)

        # now as everything else failed, try to do it symbolically with call

        cdef list g = list(parent.gens())

        if fixed is not None:
            for m,v in fixed.iteritems():
                if isinstance(m, (int, Integer)):
                    mi = m+1
                elif isinstance(m, MPolynomial_libsingular) and m.parent() is parent:
                    for i from 0 < i <= _ring.N:
                        if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                            mi = i
                            break
                    if i > _ring.N:
                        raise TypeError("key does not match")
                else:
                    raise TypeError("keys do not match self's parent")

                g[mi-1] = v

        gd = parent.gens_dict()
        for m,v in kw.iteritems():
            m = gd[m]
            for i from 0 < i <= _ring.N:
                if p_GetExp((<MPolynomial_libsingular>m)._poly, i, _ring) != 0:
                    mi = i
                    break
            if i > _ring.N:
                raise TypeError("key does not match")

            g[mi-1] = v

        return self(*g)

    def monomials(self):
        """
        Return the list of monomials in self. The returned list is
        decreasingly ordered by the term ordering of
        ``self.parent()``.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = x + 3/2*y*z^2 + 2/3
            sage: f.monomials()
            [y*z^2, x, 1]
            sage: f = P(3/2)
            sage: f.monomials()
            [1]

        TESTS::

            sage: P.<x,y,z> = QQ[]
            sage: f = x
            sage: f.monomials()
            [x]

        Check if :trac:`12706` is fixed::

            sage: f = P(0)
            sage: f.monomials()
            []

        Check if :trac:`7152` is fixed::

            sage: x=var('x')
            sage: K.<rho> = NumberField(x**2 + 1)
            sage: R.<x,y> = QQ[]
            sage: p = rho*x
            sage: q = x
            sage: p.monomials()
            [x]
            sage: q.monomials()
            [x]
            sage: p.monomials()
            [x]
        """
        l = list()
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *_ring = parent._ring
        if(_ring != currRing): rChangeCurrRing(_ring)
        cdef poly *p = p_Copy(self._poly, _ring)
        cdef poly *t

        if p == NULL:
            return []

        while p:
            t = pNext(p)
            p.next = NULL
            p_SetCoeff(p, n_Init(1,_ring), _ring)
            p_Setm(p, _ring)
            l.append( new_MP(parent,p) )
            p = t

        return l

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate
        polynomial.

        EXAMPLES::

            sage: P.<x, y> = QQ[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x^2
            sage: f.constant_coefficient()
            0
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        if p == NULL:
            return self._parent._base._zero_element

        while p.next:
            p = pNext(p)

        if p_LmIsConstant(p, r):
            return si2sa(p_GetCoeff(p, r), r, self._parent._base)
        else:
            return self._parent._base._zero_element

    def univariate_polynomial(self, R=None):
        """
        Returns a univariate polynomial associated to this
        multivariate polynomial.

        INPUT:

        - ``R`` - (default: ``None``) PolynomialRing

        If this polynomial is not in at most one variable, then a
        ``ValueError`` exception is raised.  This is checked using the
        :meth:`is_univariate()` method.  The new Polynomial is over
        the same base ring as the given ``MPolynomial`` and in the
        variable ``x`` if no ring ``R`` is provided.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: polynomial must involve at most one variable
            sage: g = f.subs({x:10}); g
            700*y^2 - 2*y + 305
            sage: g.univariate_polynomial ()
            700*y^2 - 2*y + 305
            sage: g.univariate_polynomial(PolynomialRing(QQ,'z'))
            700*z^2 - 2*z + 305

        Here's an example with a constant multivariate polynomial::

            sage: g = R(1)
            sage: h = g.univariate_polynomial(); h
            1
            sage: h.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        k = self.base_ring()

        if not self.is_univariate():
            raise TypeError, "polynomial must involve at most one variable"

        #construct ring if none
        if R is None:
            if self.is_constant():
                R = self.base_ring()['x']
            else:
                R = self.base_ring()[str(self.variables()[0])]

        zero = k(0)
        coefficients = [zero] * (self.degree() + 1)

        if(r != currRing): rChangeCurrRing(r)

        while p:
            coefficients[p_Totaldegree(p, r)] = si2sa(p_GetCoeff(p, r), r, k)
            p = pNext(p)

        return R(coefficients)

    def is_univariate(self):
        """
        Return ``True`` if self is a univariate polynomial, that is if
        self contains only one variable.

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
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

    def _variable_indices_(self, sort=True):
        """
        Return the indices of all variables occurring in self.  This
        index is the index as Sage uses them (starting at zero), not
        as SINGULAR uses them (starting at one).

        INPUT:

        - ``sort`` - specifies whether the indices shall be sorted

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
            sage: f = x*z^2 + z + 1
            sage: f._variable_indices_()
            [0, 2]

        """
        cdef poly *p
        cdef ring *r = self._parent_ring
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

    def variables(self):
        """
        Return a tuple of all variables occurring in self.

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
            sage: f = x*z^2 + z + 1
            sage: f.variables()
            (x, z)
        """
        cdef poly *p
        cdef poly *v
        cdef ring *r = self._parent_ring
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
                    l.append(new_MP(self._parent, v))
                    si.add(i)
            p = pNext(p)
        return tuple(sorted(l,reverse=True))

    def variable(self, i=0):
        """

        Return the i-th variable occurring in self. The index i is the
        index in ``self.variables()``.

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
            sage: f = x*z^2 + z + 1
            sage: f.variables()
            (x, z)
            sage: f.variable(1)
            z
        """
        return self.variables()[i]

    def nvariables(self):
        """
        Return the number variables in this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(127))
            sage: f = x*y + z
            sage: f.nvariables()
            3
            sage: f = x + y
            sage: f.nvariables()
            2
        """
        return len(self._variable_indices_(sort=False))

    cpdef is_constant(self):
        """
        Return ``True`` if this polynomial is constant.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(127))
            sage: x.is_constant()
            False
            sage: P(1).is_constant()
            True
        """
        return bool(p_IsConstant(self._poly, self._parent_ring))

    def lm(MPolynomial_libsingular self):
        """
        Returns the lead monomial of self with respect to the term
        order of ``self.parent()``. In Sage a monomial is a product of
        variables in some power without a coefficient.

        EXAMPLES::

            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: f = x^1*y^2 + y^3*z^4
            sage: f.lm()
            x*y^2
            sage: f = x^3*y^2*z^4 + x^3*y^2*z^1
            sage: f.lm()
            x^3*y^2*z^4

            sage: R.<x,y,z>=PolynomialRing(QQ,3,order='deglex')
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
        cdef ring *_ring = self._parent_ring
        if self._poly == NULL:
            return self._parent._zero_element
        _p = p_Head(self._poly, _ring)
        p_SetCoeff(_p, n_Init(1,_ring), _ring)
        p_Setm(_p,_ring)
        return new_MP(self._parent, _p)

    def lc(MPolynomial_libsingular self):
        """
        Leading coefficient of this polynomial with respect to the
        term order of ``self.parent()``.

        EXAMPLES::

            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: f = 3*x^1*y^2 + 2*y^3*z^4
            sage: f.lc()
            3

            sage: f = 5*x^3*y^2*z^4 + 4*x^3*y^2*z^1
            sage: f.lc()
            5
        """

        cdef poly *_p
        cdef ring *_ring = self._parent_ring
        cdef number *_n

        if self._poly == NULL:
            return self._parent._base._zero_element

        if(_ring != currRing): rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        ret =  si2sa(_n, _ring, self._parent._base)
        p_Delete(&_p, _ring)
        return ret

    def lt(MPolynomial_libsingular self):
        """
        Leading term of this polynomial. In Sage a term is a product
        of variables in some power and a coefficient.

        EXAMPLES::

            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: f = 3*x^1*y^2 + 2*y^3*z^4
            sage: f.lt()
            3*x*y^2

            sage: f = 5*x^3*y^2*z^4 + 4*x^3*y^2*z^1
            sage: f.lt()
            -2*x^3*y^2*z^4
        """
        if self._poly == NULL:
            return self._parent._zero_element
        return new_MP(self._parent, p_Head(self._poly, self._parent_ring))

    def is_zero(self):
        """
        Return ``True`` if this polynomial is zero.

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: x.is_zero()
            False
            sage: (x-x).is_zero()
            True
        """
        if self._poly is NULL:
            return True
        else:
            return False

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: bool(x) # indirect doctest
            True
            sage: bool(x-x)
            False
        """
        if self._poly:
            return True
        else:
            return False

    cpdef RingElement _floordiv_(self, RingElement right):
        """
        Perform division with remainder and return the quotient.

        INPUT:

        - ``right`` - something coercible to an MPolynomial_libsingular
          in ``self.parent()``

        EXAMPLES::

            sage: R.<x,y,z> = GF(32003)[]
            sage: f = y*x^2 + x + 1
            sage: f//x
            x*y + 1
            sage: f//y
            x^2

            sage: P.<x,y> = ZZ[]
            sage: x//y
            0
            sage: (x+y)//y
            1

            sage: P.<x,y> = QQ[]
            sage: (x+y)//y
            1
            sage: (x)//y
            0

            sage: P.<x,y> = Zmod(1024)[]
            sage: (x+y)//x
            1
            sage: (x+y)//(2*x)
            Traceback (most recent call last):
            ...
            NotImplementedError: Division of multivariate polynomials over non fields by non-monomials not implemented.
        """
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *r = self._parent_ring
        if(r != currRing): rChangeCurrRing(r)
        cdef MPolynomial_libsingular _self, _right
        cdef poly *quo
        cdef poly *temp
        cdef poly *p

        if right.is_zero():
            raise ZeroDivisionError

        if self._parent._base.is_finite() and self._parent._base.characteristic() > 1<<29:
            raise NotImplementedError, "Division of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented."

        _self = <MPolynomial_libsingular>self
        _right = <MPolynomial_libsingular>right

        if r.ringtype != 0:
            if r.ringtype == 4:
                P = parent.change_ring(RationalField())
                f = P(self)//P(right)
                CM = list(f)
                return parent(sum([c.floor()*m for c,m in CM]))
            if _right.is_monomial():
                p = _self._poly
                quo = p_ISet(0,r)
                if r != currRing: rChangeCurrRing(r)   # pDivide
                while p:
                    if p_DivisibleBy(_right._poly, p, r):
                        temp = pDivide(p, _right._poly)
                        p_SetCoeff0(temp, n_Copy(p_GetCoeff(p, r), r), r)
                        quo = p_Add_q(quo, temp, r)
                    p = pNext(p)
                return new_MP(parent, quo)
            raise NotImplementedError("Division of multivariate polynomials over non fields by non-monomials not implemented.")

        cdef int count = singular_polynomial_length_bounded(_self._poly,15)
        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
            sig_on()
        if r!=currRing: rChangeCurrRing(r)   # singclap_pdivide
        quo = singclap_pdivide( _self._poly, _right._poly )
        if count >= 15:
            sig_off()
        f = new_MP(parent, quo)

        return f

    def factor(self, proof=True):
        """
        Return the factorization of this polynomial.

        INPUT:

        - ``proof`` - ignored.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            x^5 + 2*x^3*y^2 + x^4 + 2*x^2*y^2 + x^3 + 2*x*y^2
            sage: F = f.factor()
            sage: F
            x * (x^2 + x + 1) * (x^2 + 2*y^2)

        Next we factor the same polynomial, but over the finite field
        of order 3.::

            sage: R.<x, y> = GF(3)[]
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            x^5 - x^3*y^2 + x^4 - x^2*y^2 + x^3 - x*y^2
            sage: F = f.factor()
            sage: F # order is somewhat random
            (-1) * x * (-x + y) * (x + y) * (x - 1)^2

        Next we factor a polynomial, but over a finite field of order 9.::

            sage: K.<a> = GF(3^2)
            sage: R.<x, y> = K[]
            sage: f = (x^3 + 2*a*y^2*x) * (x^2 + x + 1); f
            x^5 + (-a)*x^3*y^2 + x^4 + (-a)*x^2*y^2 + x^3 + (-a)*x*y^2
            sage: F = f.factor()
            sage: F
            ((-a)) * x * (x - 1)^2 * ((-a + 1)*x^2 + y^2)
            sage: f - F
            0

        Next we factor a polynomial over a number field.::

            sage: p = var('p')
            sage: K.<s> = NumberField(p^3-2)
            sage: KXY.<x,y> = K[]
            sage: factor(x^3 - 2*y^3)
            (x + (-s)*y) * (x^2 + (s)*x*y + (s^2)*y^2)
            sage: k = (x^3-2*y^3)^5*(x+s*y)^2*(2/3 + s^2)
            sage: k.factor()
            ((s^2 + 2/3)) * (x + (s)*y)^2 * (x + (-s)*y)^5 * (x^2 + (s)*x*y + (s^2)*y^2)^5

        This shows that ticket :trac:`2780` is fixed, i.e. that the unit
        part of the factorization is set correctly::

            sage: x = var('x')
            sage: K.<a> = NumberField(x^2 + 1)
            sage: R.<y, z> = PolynomialRing(K)
            sage: f = 2*y^2 + 2*z^2
            sage: F = f.factor(); F.unit()
            2

        Another example::

            sage: R.<x,y,z> = GF(32003)[]
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
            (-2) * (c - d) * (-b + c) * (b - d) * (-a + c) * (-a + b) * (a - d)
            sage: F[0][0]
            c - d
            sage: F.unit()
            -2

        Constant elements are factorized in the base rings. ::

            sage: P.<x,y> = ZZ[]
            sage: P(2^3*7).factor()
            2^3 * 7

        Factorization for finite prime fields with characteristic
        `> 2^{29}` is not supported ::

            sage: q = 1073741789
            sage: T.<aa, bb> = PolynomialRing(GF(q))
            sage: f = aa^2 + 12124343*bb*aa + 32434598*bb^2
            sage: f.factor()
            Traceback (most recent call last):
            ...
            NotImplementedError: Factorization of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.

        Factorization over the integers is now supported, see :trac:`17840`::

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: f = 12 * (3*x*y + 4) * (5*x - 2) * (2*y + 7)^2
            sage: f.factor()
            2^2 * 3 * (2*y + 7)^2 * (5*x - 2) * (3*x*y + 4)
            sage: g = -12 * (x^2 - y^2)
            sage: g.factor()
            (-1) * 2^2 * 3 * (x - y) * (x + y)
            sage: factor(-4*x*y - 2*x + 2*y + 1)
            (-1) * (2*y + 1) * (2*x - 1)

        Factorization over non-integral domains is not supported ::

            sage: R.<x,y> = PolynomialRing(Zmod(4))
            sage: f = (2*x + 1) * (x^2 + x + 1)
            sage: f.factor()
            Traceback (most recent call last):
            ...
            NotImplementedError: Factorization of multivariate polynomials over Ring of integers modulo 4 is not implemented.

        TESTS:

        This shows that :trac:`10270` is fixed::

            sage: R.<x,y,z> = GF(3)[]
            sage: f = x^2*z^2+x*y*z-y^2
            sage: f.factor()
            x^2*z^2 + x*y*z - y^2

        This checks that :trac:`11838` is fixed::

            sage: K = GF(4,'a')
            sage: a = K.gens()[0]
            sage: R.<x,y> = K[]
            sage: p=x^8*y^3 + x^2*y^9 + a*x^9 + a*x*y^4
            sage: q=y^11 + (a)*y^10 + (a + 1)*x*y^3
            sage: f = p*q
            sage: f.factor()
            x * y^3 * (y^8 + (a)*y^7 + (a + 1)*x) * (x^7*y^3 + x*y^9 + (a)*x^8 + (a)*y^4)

        We test several examples which were known to return wrong
        results in the past (see :trac:`10902`)::

            sage: R.<x,y> = GF(2)[]
            sage: p = x^3*y^7 + x^2*y^6 + x^2*y^3
            sage: q = x^3*y^5
            sage: f = p*q
            sage: p.factor()*q.factor()
            x^5 * y^8 * (x*y^4 + y^3 + 1)
            sage: f.factor()
            x^5 * y^8 * (x*y^4 + y^3 + 1)
            sage: f.factor().expand() == f
            True

        ::

            sage: R.<x,y> = GF(2)[]
            sage: p=x^8 + y^8; q=x^2*y^4 + x
            sage: f=p*q
            sage: lf = f.factor()
            sage: f-lf
            0

        ::

            sage: R.<x,y> = GF(3)[]
            sage: p = -x*y^9 + x
            sage: q = -x^8*y^2
            sage: f = p*q
            sage: f
            x^9*y^11 - x^9*y^2
            sage: f.factor()
            y^2 * (y - 1)^9 * x^9
            sage: f - f.factor()
            0

        ::

            sage: R.<x,y> = GF(5)[]
            sage: p=x^27*y^9 + x^32*y^3 + 2*x^20*y^10 - x^4*y^24 - 2*x^17*y
            sage: q=-2*x^10*y^24 + x^9*y^24 - 2*x^3*y^30
            sage: f=p*q; f-f.factor()
            0

        ::

            sage: R.<x,y> = GF(7)[]
            sage: p=-3*x^47*y^24
            sage: q=-3*x^47*y^37 - 3*x^24*y^49 + 2*x^56*y^8 + 3*x^29*y^15 - x^2*y^33
            sage: f=p*q
            sage: f-f.factor()
            0

        The following examples used to give a Segmentation Fault, see
        :trac:`12918` and :trac:`13129`::

            sage: R.<x,y> = GF(2)[]
            sage: f = x^6 + x^5 + y^5 + y^4
            sage: f.factor()
            x^6 + x^5 + y^5 + y^4
            sage: f = x^16*y + x^10*y + x^9*y + x^6*y + x^5 + x*y + y^2
            sage: f.factor()
            x^16*y + x^10*y + x^9*y + x^6*y + x^5 + x*y + y^2

        Test :trac:`12928`::

            sage: R.<x,y> = GF(2)[]
            sage: p = x^2 + y^2 + x + 1
            sage: q = x^4 + x^2*y^2 + y^4 + x*y^2 + x^2 + y^2 + 1
            sage: factor(p*q)
            (x^2 + y^2 + x + 1) * (x^4 + x^2*y^2 + y^4 + x*y^2 + x^2 + y^2 + 1)

        Check that :trac:`13770` is fixed::

            sage: U.<y,t> = GF(2)[]
            sage: f = y*t^8 + y^5*t^2 + y*t^6 + t^7 + y^6 + y^5*t + y^2*t^4 + y^2*t^2 + y^2*t + t^3 + y^2 + t^2
            sage: l = f.factor()
            sage: l[0][0]==t^2 + y + t + 1 or l[1][0]==t^2 + y + t + 1
            True

        The following used to sometimes take a very long time or get
        stuck, see :trac:`12846`. These 100 iterations should take less
        than 1 second::

            sage: K.<a> = GF(4)
            sage: R.<x,y> = K[]
            sage: f = (a + 1)*x^145*y^84 + (a + 1)*x^205*y^17 + x^32*y^112 + x^92*y^45
            sage: for i in range(100):
            ....:     assert len(f.factor()) == 4

        Test for :trac:`20435`::

            sage: x,y = polygen(ZZ,'x,y')
            sage: p = x**2-y**2
            sage: z = factor(p); z
            (x - y) * (x + y)
            sage: z[0][0].parent()
            Multivariate Polynomial Ring in x, y over Integer Ring
        """
        cdef ring *_ring = self._parent_ring
        cdef poly *ptemp
        cdef intvec *iv
        cdef int *ivv
        cdef ideal *I
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef int i

        if _ring!=currRing: rChangeCurrRing(_ring)

        if p_IsConstant(self._poly, _ring):
            return self.constant_coefficient().factor()

        if not self._parent._base.is_field():
            try:
                frac_field = self._parent._base.fraction_field()
                F = self.change_ring(frac_field).factor()
                FF = [(self._parent(f[0]), f[1]) for f in F]
                U = self._parent._base(F.unit()).factor()
                return Factorization(list(U) + FF, unit=U.unit())
            except Exception:
                raise NotImplementedError("Factorization of multivariate polynomials over %s is not implemented."%self._parent._base)

        if self._parent._base.is_finite():
            if self._parent._base.characteristic() > 1<<29:
                raise NotImplementedError, "Factorization of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented."

        # I make a temporary copy of the poly in self because singclap_factorize appears to modify it's parameter
        ptemp = p_Copy(self._poly,_ring)
        iv = NULL
        sig_on()
        if _ring!=currRing: rChangeCurrRing(_ring)   # singclap_factorize
        I = singclap_factorize ( ptemp, &iv , 0)
        sig_off()

        ivv = iv.ivGetVec()
        v = [(new_MP(parent, p_Copy(I.m[i],_ring)) , ivv[i])   for i in range(1,I.ncols)]
        v = [(f,m) for f,m in v if f!=0] # we might have zero in there
        unit = new_MP(parent, p_Copy(I.m[0],_ring))

        F = Factorization(v,unit)
        F.sort()

        del iv
        id_Delete(&I,_ring)

        return F

    def lift(self, I):
        """
        given an ideal ``I = (f_1,...,f_r)`` and some ``g (== self)`` in ``I``,
        find ``s_1,...,s_r`` such that ``g = s_1 f_1 + ... + s_r f_r``.

        A ``ValueError`` exception is raised if ``g (== self)`` does not belong to ``I``.

        EXAMPLES::

            sage: A.<x,y> = PolynomialRing(QQ,2,order='degrevlex')
            sage: I = A.ideal([x^10 + x^9*y^2, y^8 - x^2*y^7 ])
            sage: f = x*y^13 + y^12
            sage: M = f.lift(I)
            sage: M
            [y^7, x^7*y^2 + x^8 + x^5*y^3 + x^6*y + x^3*y^4 + x^4*y^2 + x*y^5 + x^2*y^3 + y^4]
            sage: sum( map( mul , zip( M, I.gens() ) ) ) == f
            True

        Check that :trac:`13671` is fixed::

            sage: R.<x1,x2> = QQ[]
            sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
            sage: f = I.gen(0) + x2*I.gen(1)
            sage: f.lift(I)
            [1, x2]
            sage: (f+1).lift(I)
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not in the ideal
            sage: f.lift(I)
            [1, x2]

        TESTS:

        Check that :trac:`13714` is fixed::

            sage: R.<x1,x2> = QQ[]
            sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
            sage: R.one().lift(I)
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not in the ideal
            sage: foo = I.complete_primary_decomposition() # indirect doctest
            sage: foo[0][0]
            Ideal (x2 - 1, x1 - 1) of Multivariate Polynomial Ring in x1, x2 over Rational Field

        """
        global errorreported, libfac_interruptflag
        if not self._parent._base.is_field():
            raise NotImplementedError, "Lifting of multivariate polynomials over non-fields is not implemented."

        cdef ideal *fI = idInit(1,1)
        cdef ideal *_I
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef int i = 0
        cdef int j
        cdef ring *r = self._parent_ring
        cdef ideal *res

        if isinstance(I, MPolynomialIdeal):
            I = I.gens()

        _I = idInit(len(I),1)

        for f in I:
            if not (isinstance(f,MPolynomial_libsingular) \
                    and (<MPolynomial_libsingular>f)._parent is parent):
                try:
                    f = parent._coerce_c(f)
                except TypeError as msg:
                    id_Delete(&fI,r)
                    id_Delete(&_I,r)
                    raise TypeError, msg

            _I.m[i] = p_Copy((<MPolynomial_libsingular>f)._poly, r)
            i+=1

        fI.m[0]= p_Copy(self._poly, r)

        if r!=currRing: rChangeCurrRing(r)  # idLift
        res = idLift(_I, fI, NULL, 0, 0, 0)
        if errorreported != 0 or libfac_interruptflag != 0:
            errorcode = errorreported
            errorreported = 0
            libfac_interruptflag = 0
            if errorcode == 1:
                raise ValueError("polynomial is not in the ideal")
            raise RuntimeError

        l = []
        for i from 0 <= i < IDELEMS(res):
            for j from 1 <= j <= IDELEMS(_I):
                l.append( new_MP(parent, pTakeOutComp1(&res.m[i], j)) )

        id_Delete(&fI, r)
        id_Delete(&_I, r)
        id_Delete(&res, r)
        return Sequence(l, check=False, immutable=True)

    def reduce(self,I):
        """
        Return the normal form of self w.r.t. ``I``, i.e. return the
        remainder of this polynomial with respect to the polynomials
        in ``I``. If the polynomial set/list ``I`` is not a (strong)
        Groebner basis the result is not canonical.

        A strong Groebner basis ``G`` of ``I`` implies that for every
        leading term ``t`` of ``I`` there exists an element ``g`` of ``G``,
        such that the leading term of ``g`` divides ``t``.

        INPUT:

        - ``I`` - a list/set of polynomials. If ``I`` is an ideal, the
          generators are used.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f1 = -2 * x^2 + x^3
            sage: f2 = -2 * y + x* y
            sage: f3 = -x^2 + y^2
            sage: F = Ideal([f1,f2,f3])
            sage: g = x*y - 3*x*y^2
            sage: g.reduce(F)
            -6*y^2 + 2*y
            sage: g.reduce(F.gens())
            -6*y^2 + 2*y

        `\ZZ` is also supported. ::

            sage: P.<x,y,z> = ZZ[]
            sage: f1 = -2 * x^2 + x^3
            sage: f2 = -2 * y + x* y
            sage: f3 = -x^2 + y^2
            sage: F = Ideal([f1,f2,f3])
            sage: g = x*y - 3*x*y^2
            sage: g.reduce(F)
            -6*y^2 + 2*y
            sage: g.reduce(F.gens())
            -6*y^2 + 2*y

            sage: f = 3*x
            sage: f.reduce([2*x,y])
            3*x
        """
        cdef ideal *_I
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef int i = 0
        cdef ring *r = self._parent_ring
        cdef poly *res

        if r!=currRing: rChangeCurrRing(r)

        if isinstance(I, MPolynomialIdeal):
            try:
                strat = I._groebner_strategy()
                return strat.normal_form(self)
            except (TypeError, NotImplementedError):
                pass
            I = I.gens()

        _I = idInit(len(I),1)
        for f in I:
            if not (isinstance(f,MPolynomial_libsingular) \
                   and (<MPolynomial_libsingular>f)._parent is parent):
                try:
                    f = parent._coerce_c(f)
                except TypeError as msg:
                    id_Delete(&_I,r)
                    raise TypeError, msg

            _I.m[i] = p_Copy((<MPolynomial_libsingular>f)._poly, r)
            i+=1

        #the second parameter would be qring!
        if r!=currRing: rChangeCurrRing(r)  # kNF
        res = kNF(_I, NULL, self._poly)
        id_Delete(&_I,r)
        return new_MP(parent,res)

    @coerce_binop
    def gcd(self, right, algorithm=None, **kwds):
        """
        Return the greatest common divisor of self and right.

        INPUT:

        - ``right`` - polynomial
        - ``algorithm``
          - ``ezgcd`` - EZGCD algorithm
          - ``modular`` - multi-modular algorithm (default)
        - ``**kwds`` - ignored

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = (x*y*z)^6 - 1
            sage: g = (x*y*z)^4 - 1
            sage: f.gcd(g)
            x^2*y^2*z^2 - 1
            sage: GCD([x^3 - 3*x + 2, x^4 - 1, x^6 -1])
            x - 1

            sage: R.<x,y> = QQ[]
            sage: f = (x^3 + 2*y^2*x)^2
            sage: g = x^2*y^2
            sage: f.gcd(g)
            x^2

        We compute a gcd over a finite field::

            sage: F.<u> = GF(31^2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3
            sage: gcd(p,q)  # yes, twice -- tests that singular ring is properly set.
            x^3 + (u + 1)*y^3 + z^3

        We compute a gcd over a number field::

            sage: x = polygen(QQ)
            sage: F.<u> = NumberField(x^3 - 2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3

        TESTS::

            sage: Q.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P(0).gcd(Q(0))
            0
            sage: x.gcd(1)
            1

            sage: k.<a> = GF(9)
            sage: R.<x,y> = PolynomialRing(k)
            sage: f = R.change_ring(GF(3)).gen()
            sage: g = x+y
            sage: g.gcd(f)
            1
            sage: x.gcd(R.change_ring(GF(3)).gen())
            x
        """
        cdef MPolynomial_libsingular _right
        cdef poly *_res
        cdef ring *_ring = self._parent_ring

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

        if _ring.ringtype != 0:
            if _ring.ringtype == 4:
                P = self._parent.change_ring(RationalField())
                res = P(self).gcd(P(right))
                coef = sage.rings.integer.GCD_list(self.coefficients() + right.coefficients())
                return self._parent(coef*res)

            raise NotImplementedError("GCD over rings not implemented.")

        if self._parent._base.is_finite() and self._parent._base.characteristic() > 1<<29:
            raise NotImplementedError, "GCD of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented."

        if(_ring != currRing): rChangeCurrRing(_ring)

        if not (isinstance(right, MPolynomial_libsingular) \
                    and (<MPolynomial_libsingular>right)._parent is self._parent):
            _right = self._parent._coerce_c(right)
        else:
            _right = <MPolynomial_libsingular>right

        cdef int count = singular_polynomial_length_bounded(self._poly,20) \
            + singular_polynomial_length_bounded(_right._poly,20)
        if count >= 20:
            sig_on()
        if _ring!=currRing: rChangeCurrRing(_ring)  # singclap_gcd
        _res = singclap_gcd(p_Copy(self._poly, _ring), p_Copy(_right._poly, _ring))
        if count >= 20:
            sig_off()

        res = new_MP(self._parent, _res)
        return res

    @coerce_binop
    def lcm(self, MPolynomial_libsingular g):
        """
        Return the least common multiple of ``self`` and `g`.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: p = (x+y)*(y+z)
            sage: q = (z^4+2)*(y+z)
            sage: lcm(p,q)
            x*y*z^4 + y^2*z^4 + x*z^5 + y*z^5 + 2*x*y + 2*y^2 + 2*x*z + 2*y*z

            sage: P.<x,y,z> = ZZ[]
            sage: p = 2*(x+y)*(y+z)
            sage: q = 3*(z^4+2)*(y+z)
            sage: lcm(p,q)
            6*x*y*z^4 + 6*y^2*z^4 + 6*x*z^5 + 6*y*z^5 + 12*x*y + 12*y^2 + 12*x*z + 12*y*z

            sage: r.<x,y> = PolynomialRing(GF(2**8, 'a'), 2)
            sage: a = r.base_ring().0
            sage: f = (a^2+a)*x^2*y + (a^4+a^3+a)*y + a^5
            sage: f.lcm(x^4)
            (a^2 + a)*x^6*y + (a^4 + a^3 + a)*x^4*y + (a^5)*x^4

            sage: w = var('w')
            sage: r.<x,y> = PolynomialRing(NumberField(w^4 + 1, 'a'), 2)
            sage: a = r.base_ring().0
            sage: f = (a^2+a)*x^2*y + (a^4+a^3+a)*y + a^5
            sage: f.lcm(x^4)
            (a^2 + a)*x^6*y + (a^3 + a - 1)*x^4*y + (-a)*x^4
        """
        cdef ring *_ring = self._parent_ring
        cdef poly *ret
        cdef poly *prod
        cdef poly *gcd
        cdef MPolynomial_libsingular _g
        if _ring!=currRing: rChangeCurrRing(_ring)

        if _ring.ringtype != 0:
            if _ring.ringtype == 4:
                P = self.parent().change_ring(RationalField())
                py_gcd = P(self).gcd(P(g))
                py_prod = P(self*g)
                return self.parent(py_prod//py_gcd)
            else:
                raise TypeError("LCM over non-integral domains not available.")

        if self._parent is not g._parent:
            _g = self._parent._coerce_c(g)
        else:
            _g = <MPolynomial_libsingular>g

        if self._parent._base.is_finite() and self._parent._base.characteristic() > 1<<29:
            raise NotImplementedError, "LCM of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented."

        cdef int count = singular_polynomial_length_bounded(self._poly,20) \
            + singular_polynomial_length_bounded(_g._poly,20)
        if count >= 20:
            sig_on()
        if _ring!=currRing: rChangeCurrRing(_ring)  # singclap_gcd
        gcd = singclap_gcd(p_Copy(self._poly, _ring), p_Copy(_g._poly, _ring))
        prod = pp_Mult_qq(self._poly, _g._poly, _ring)
        ret = singclap_pdivide(prod , gcd )
        p_Delete(&prod, _ring)
        p_Delete(&gcd, _ring)
        if count >= 20:
            sig_off()
        return new_MP(self._parent, ret)

    def is_squarefree(self):
        """
        Return ``True`` if this polynomial is square free.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: f= x^2 + 2*x*y + 1/2*z
            sage: f.is_squarefree()
            True
            sage: h = f^2
            sage: h.is_squarefree()
            False
        """
        cdef ring *_ring = self._parent_ring

        if self._parent._base.is_finite() and self._parent._base.characteristic() > 1<<29:
            raise NotImplementedError, "is_squarefree of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented."

        if(_ring != currRing): rChangeCurrRing(_ring)
        return bool(singclap_isSqrFree(self._poly))

    @coerce_binop
    def quo_rem(self, MPolynomial_libsingular right):
        """
        Returns quotient and remainder of self and right.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = y*x^2 + x + 1
            sage: f.quo_rem(x)
            (x*y + 1, 1)
            sage: f.quo_rem(y)
            (x^2, x + 1)

            sage: R.<x,y> = ZZ[]
            sage: f = 2*y*x^2 + x + 1
            sage: f.quo_rem(x)
            (2*x*y + 1, 1)
            sage: f.quo_rem(y)
            (2*x^2, x + 1)
            sage: f.quo_rem(3*x)
            (0, 2*x^2*y + x + 1)

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: R(0).quo_rem(R(1))
            (0, 0)
            sage: R(1).quo_rem(R(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        """
        cdef poly *quo
        cdef poly *rem
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *r = self._parent_ring
        if r!=currRing: rChangeCurrRing(r)

        if right.is_zero():
            raise ZeroDivisionError

        if not self._parent._base.is_field():
            py_quo = self//right
            py_rem = self - right*py_quo
            return py_quo, py_rem

        if self._parent._base.is_finite() and self._parent._base.characteristic() > 1<<29:
            raise NotImplementedError, "Division of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented."

        cdef int count = singular_polynomial_length_bounded(self._poly,15)
        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
            sig_on()
        if r!=currRing: rChangeCurrRing(r)   # singclap_pdivide
        quo = singclap_pdivide( self._poly, right._poly )
        rem = p_Add_q(p_Copy(self._poly, r), p_Neg(pp_Mult_qq(right._poly, quo, r), r), r)
        if count >= 15:
            sig_off()
        return new_MP(parent, quo), new_MP(parent, rem)

    def _singular_init_(self, singular=singular_default):
        """
        Return a SINGULAR (as in the computer algebra system) string
        representation for this element.

        INPUT:

        - ``singular`` - interpreter (default: ``singular_default``)

        - ``have_ring`` - should the correct ring not be set in
           SINGULAR first (default: ``False``)

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(127),3)
            sage: x._singular_init_()
            'x'
            sage: (x^2+37*y+128)._singular_init_()
            'x2+37y+1'

        TESTS::

            sage: P(0)._singular_init_()
            '0'
        """
        self._parent._singular_().set_ring()
        return self._repr_short_()

    def sub_m_mul_q(self, MPolynomial_libsingular m, MPolynomial_libsingular q):
        """
        Return ``self - m*q``, where ``m`` must be a monomial and
        ``q`` a polynomial.

        INPUT:

        - ``m`` - a monomial
        - ``q`` - a polynomial

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: x.sub_m_mul_q(y,z)
            -y*z + x

        TESTS::

            sage: Q.<x,y,z>=PolynomialRing(QQ,3)
            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: P(0).sub_m_mul_q(P(0),P(1))
            0
            sage: x.sub_m_mul_q(Q.gen(1),Q.gen(2))
            -y*z + x
         """
        cdef ring *r = self._parent_ring

        if not self._parent is m._parent:
            m = self._parent._coerce_c(m)
        if not self._parent is q._parent:
            q = self._parent._coerce_c(q)

        if m._poly and m._poly.next:
            raise ArithmeticError, "m must be a monomial."
        elif not m._poly:
            return self

        cdef int le = p_GetMaxExp(m._poly, r)
        cdef int lr = p_GetMaxExp(q._poly, r)
        cdef int esum = le + lr

        overflow_check(esum, r)

        return new_MP(self._parent, p_Minus_mm_Mult_qq(p_Copy(self._poly, r), m._poly, q._poly, r))

    def _macaulay2_(self, macaulay2=macaulay2):
        """
        Return a Macaulay2 string representation of this polynomial.

        .. note::

           Two identical rings are not canonically isomorphic in M2,
           so we require the user to explicitly set the ring, since
           there is no way to know if the ring has been set or not,
           and setting it twice screws everything up.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(GF(7), 2)
            sage: f = (x^3 + 2*y^2*x)^7; f          # indirect doctest
            x^21 + 2*x^7*y^14

            sage: h = macaulay2(f); h               # optional - macaulay2
             21     7 14
            x   + 2x y
            sage: k = macaulay2(x+y); k             # optional - macaulay2
            x + y
            sage: k + h                             # optional - macaulay2
             21     7 14
            x   + 2x y   + x + y
            sage: R(h)                              # optional - macaulay2
            x^21 + 2*x^7*y^14
            sage: R(h^20) == f^20                   # optional - macaulay2
            True
        """
        m2_parent = macaulay2(self.parent())
        macaulay2.use(m2_parent)
        return macaulay2(repr(self))

    def add_m_mul_q(self, MPolynomial_libsingular m, MPolynomial_libsingular q):
        """
        Return ``self + m*q``, where ``m`` must be a monomial and
        ``q`` a polynomial.

        INPUT:

        - ``m`` - a monomial
        - ``q``  - a polynomial

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: x.add_m_mul_q(y,z)
            y*z + x

        TESTS::

            sage: R.<x,y,z>=PolynomialRing(QQ,3)
            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: P(0).add_m_mul_q(P(0),P(1))
            0
            sage: x.add_m_mul_q(R.gen(),R.gen(1))
            x*y + x
        """
        cdef ring *r = self._parent_ring

        if not self._parent is m._parent:
            m = self._parent._coerce_c(m)
        if not self._parent is q._parent:
            q = self._parent._coerce_c(q)

        if m._poly and m._poly.next:
            raise ArithmeticError, "m must be a monomial."
        elif not m._poly:
            return self

        cdef int le = p_GetMaxExp(m._poly, r)
        cdef int lr = p_GetMaxExp(q._poly, r)
        cdef int esum = le + lr

        overflow_check(esum, r)

        return new_MP(self._parent, p_Plus_mm_Mult_qq(p_Copy(self._poly, r), m._poly, q._poly, r))

    def __reduce__(self):
        """
        Serialize this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3, order='degrevlex')
            sage: f = 27/113 * x^2 + y*z + 1/2
            sage: f == loads(dumps(f))
            True

            sage: P = PolynomialRing(GF(127),3,names='abc')
            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1
            sage: f == loads(dumps(f))
            True

        TESTS:

        Verify that :trac:`9220` is fixed.

            sage: R=QQ['x']
            sage: S=QQ['x','y']
            sage: h=S.0^2
            sage: parent(h(R.0,0))
            Univariate Polynomial Ring in x over Rational Field
        """
        return sage.rings.polynomial.multi_polynomial_libsingular.unpickle_MPolynomial_libsingular, ( self._parent, self.dict() )

    def _im_gens_(self, codomain, im_gens):
        """
        INPUT:

        - ``codomain``
        - ``im_gens``

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: f = R.hom([y,x], R)
            sage: f(x^2 + 3*y^5)
            3*x^5 + y^2

            sage: R.<a,b,c,d> = QQ[]
            sage: S.<u> = QQ[]
            sage: h = R.hom([0,0,0,u], S) # indirect doctest
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
        Differentiates this polynomial with respect to the provided
        variable. This is completely symbolic so it is also defined
        over finite fields.

        INPUT:

        - ``variable`` - the derivative is taken with respect to variable

        - ``have_ring`` - ignored, accepted for compatibility reasons

        .. note:: See also :meth:`derivative`

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
            sage: f._derivative(x)
            9*x^2*y^2 + 3
            sage: f._derivative(y)
            6*x^3*y + 10*y

        The derivative is also defined over finite fields::

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
        cdef ring *_ring = self._parent_ring
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
        return new_MP(self._parent,p)

    def integral(self, MPolynomial_libsingular var):
        """
        Integrates this polynomial with respect to the provided
        variable.

        One requires that `\QQ` is contained in the ring.

        INPUT:

        - ``variable`` - the integral is taken with respect to variable

        EXAMPLES::

            sage: R.<x, y> = PolynomialRing(QQ, 2)
            sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
            sage: f.integral(x)
            3/4*x^4*y^2 + 5*x*y^2 + 3/2*x^2 + 2*x
            sage: f.integral(y)
            x^3*y^3 + 5/3*y^3 + 3*x*y + 2*y

        Check that :trac:`15896` is solved::

            sage: s = x+y
            sage: s.integral(x)+x
            1/2*x^2 + x*y + x
            sage: s.integral(x)*s
            1/2*x^3 + 3/2*x^2*y + x*y^2

        TESTS::

            sage: z, w = polygen(QQ, 'z, w')
            sage: f.integral(z)
            Traceback (most recent call last):
            ...
            TypeError: the variable is not in the same ring as self

            sage: f.integral(y**2)
            Traceback (most recent call last):
            ...
            TypeError: not a variable in the same ring as self

            sage: x,y = polygen(ZZ,'x,y')
            sage: y.integral(x)
            Traceback (most recent call last):
            ...
            TypeError: the ring must contain the rational numbers
        """
        cdef int index

        ambient_ring = var.parent()
        if ambient_ring is not self._parent:
            raise TypeError("the variable is not in the same ring as self")

        if not ambient_ring.has_coerce_map_from(RationalField()):
            raise TypeError("the ring must contain the rational numbers")

        gens = ambient_ring.gens()
        try:
            index = gens.index(var)
        except ValueError:
            raise TypeError("not a variable in the same ring as self")

        cdef poly *_p
        cdef poly *mon
        cdef ring *_ring = self._parent_ring
        if _ring != currRing:
            rChangeCurrRing(_ring)

        v = ETuple({index: 1}, len(gens))

        _p = p_ISet(0, _ring)
        for (exp, coeff) in self.dict().iteritems():
            nexp = exp.eadd(v)  # new exponent
            mon = p_Init(_ring)
            p_SetCoeff(mon, sa2si(coeff / (1 + exp[index]), _ring), _ring)
            for pos in nexp.nonzero_positions():
                overflow_check(nexp[pos], _ring)
                p_SetExp(mon, pos + 1, nexp[pos], _ring)
            p_Setm(mon, _ring)
            _p = p_Add_q(_p, mon, _ring)
        return new_MP(self._parent, _p)

    def resultant(self, MPolynomial_libsingular other, variable=None):
        """
        Compute the resultant of this polynomial and the first
        argument with respect to the variable given as the second
        argument.

        If a second argument is not provide the first variable of
        the parent is chosen.

        INPUT:

        - ``other`` - polynomial

        - ``variable`` - optional variable (default: ``None``)

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: a = x+y
            sage: b = x^3-y^3
            sage: c = a.resultant(b); c
            -2*y^3
            sage: d = a.resultant(b,y); d
            2*x^3

        The SINGULAR example::

            sage: R.<x,y,z> = PolynomialRing(GF(32003),3)
            sage: f = 3 * (x+2)^3 + y
            sage: g = x+y+z
            sage: f.resultant(g,x)
            3*y^3 + 9*y^2*z + 9*y*z^2 + 3*z^3 - 18*y^2 - 36*y*z - 18*z^2 + 35*y + 36*z - 24

        Resultants are also supported over the Integers::

            sage: R.<x,y,a,b,u>=PolynomialRing(ZZ, 5, order='lex')
            sage: r = (x^4*y^2+x^2*y-y).resultant(x*y-y*a-x*b+a*b+u,x)
            sage: r
            y^6*a^4 - 4*y^5*a^4*b - 4*y^5*a^3*u + y^5*a^2 - y^5 + 6*y^4*a^4*b^2 + 12*y^4*a^3*b*u - 4*y^4*a^2*b + 6*y^4*a^2*u^2 - 2*y^4*a*u + 4*y^4*b - 4*y^3*a^4*b^3 - 12*y^3*a^3*b^2*u + 6*y^3*a^2*b^2 - 12*y^3*a^2*b*u^2 + 6*y^3*a*b*u - 4*y^3*a*u^3 - 6*y^3*b^2 + y^3*u^2 + y^2*a^4*b^4 + 4*y^2*a^3*b^3*u - 4*y^2*a^2*b^3 + 6*y^2*a^2*b^2*u^2 - 6*y^2*a*b^2*u + 4*y^2*a*b*u^3 + 4*y^2*b^3 - 2*y^2*b*u^2 + y^2*u^4 + y*a^2*b^4 + 2*y*a*b^3*u - y*b^4 + y*b^2*u^2

        TESTS::

            sage: P.<x,y> = PolynomialRing(QQ, order='degrevlex')
            sage: a = x+y
            sage: b = x^3-y^3
            sage: c = a.resultant(b); c
            -2*y^3
            sage: d = a.resultant(b,y); d
            2*x^3


            sage: P.<x,y> = PolynomialRing(ZZ,2)
            sage: f = x+y
            sage: g=y^2+x
            sage: f.resultant(g,y)
            x^2 + x
        """
        cdef ring *_ring = self._parent_ring
        cdef poly *rt

        if variable is None:
            variable = self.parent().gen(0)

        if not self._parent is other._parent:
            raise TypeError, "first parameter needs to be an element of self.parent()"

        if not variable.parent() is self.parent():
            raise TypeError, "second parameter needs to be an element of self.parent() or None"


        if self._parent._base.is_finite() and self._parent._base.characteristic() > 1<<29:
            raise NotImplementedError, "Resultants of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented."

        if is_IntegerRing(self._parent._base):
            ret = self.change_ring(RationalField()).resultant(other.change_ring(RationalField()),
                                                              variable.change_ring(RationalField()))
            return ret.change_ring(ZZ)
        elif not self._parent._base.is_field():
            raise ValueError("Resultants require base fields or integer base ring.")

        cdef int count = singular_polynomial_length_bounded(self._poly,20) \
            + singular_polynomial_length_bounded(other._poly,20)
        if count >= 20:
            sig_on()
        if _ring != currRing: rChangeCurrRing(_ring)   # singclap_resultant
        rt =  singclap_resultant(p_Copy(self._poly, _ring), p_Copy(other._poly, _ring),p_Copy((<MPolynomial_libsingular>variable)._poly, _ring))
        if count >= 20:
            sig_off()
        return new_MP(self._parent, rt)

    def discriminant(self,variable):
        """
        Returns the discriminant of self with respect to the given variable.

        INPUT:

          - ``variable`` - The variable with respect to which we compute
              the discriminant

        OUTPUT:

          - An element of the base ring of the polynomial ring.


        EXAMPLES::

            sage: R.<x,y,z>=QQ[]
            sage: f=4*x*y^2 + 1/4*x*y*z + 3/2*x*z^2 - 1/2*z^2
            sage: f.discriminant(x)
            1
            sage: f.discriminant(y)
            -383/16*x^2*z^2 + 8*x*z^2
            sage: f.discriminant(z)
            -383/16*x^2*y^2 + 8*x*y^2

        Note that, unlike the univariate case, the result lives in
        the same ring as the polynomial::

            sage: R.<x,y>=QQ[]
            sage: f=x^5*y+3*x^2*y^2-2*x+y-1
            sage: f.discriminant(y)
            x^10 + 2*x^5 + 24*x^3 + 12*x^2 + 1
            sage: f.polynomial(y).discriminant()
            x^10 + 2*x^5 + 24*x^3 + 12*x^2 + 1
            sage: f.discriminant(y).parent()==f.polynomial(y).discriminant().parent()
            False

        AUTHOR:
            Miguel Marco
        """
        if self.is_zero():
            return self.parent().zero()
        n = self.degree(variable)
        d = self.derivative(variable)
        k = d.degree(variable)

        r = n % 4
        u = -1 # (-1)**(n*(n-1)/2)
        if r == 0 or r == 1:
            u = 1
        an = self.coefficient(variable**n)**(n - k - 2)
        return self.parent()(u * self.resultant(d, variable) * an)

    def coefficients(self):
        """
        Return the nonzero coefficients of this polynomial in a list.
        The returned list is decreasingly ordered by the term ordering
        of the parent.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]

            sage: R.<x,y,z> = PolynomialRing(QQ, order='lex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        AUTHOR:

        - Didier Deshommes
        """
        cdef poly *p
        cdef ring *r = self._parent_ring
        if r!=currRing: rChangeCurrRing(r)
        base = self._parent._base
        p = self._poly
        coeffs = list()
        while p:
            coeffs.append(si2sa(p_GetCoeff(p, r), r, base))
            p = pNext(p)
        return coeffs

    def gradient(self):
        """
        Return a list of partial derivatives of this polynomial,
        ordered by the variables of the parent.

        EXAMPLES::

           sage: P.<x,y,z> = PolynomialRing(QQ,3)
           sage: f= x*y + 1
           sage: f.gradient()
           [y, x, 0]
        """
        cdef ring *r = self._parent_ring
        cdef int k

        if r!=currRing: rChangeCurrRing(r)
        i = []
        for k from 0 < k <= r.N:
            i.append( new_MP(self._parent, pDiff(self._poly, k)))
        return i


    def numerator(self):
        """
        Return a numerator of self computed as self * self.denominator()

        If the base_field of self is the Rational Field then the
        numerator is a polynomial whose base_ring is the Integer Ring,
        this is done for compatibility to the univariate case.

        .. warning::

            This is not the numerator of the rational function
            defined by self, which would always be self since self is a
            polynomial.

        EXAMPLES:

        First we compute the numerator of a polynomial with
        integer coefficients, which is of course self.

        ::

            sage: R.<x, y> = ZZ[]
            sage: f = x^3 + 17*y + 1
            sage: f.numerator()
            x^3 + 17*y + 1
            sage: f == f.numerator()
            True

        Next we compute the numerator of a polynomial with rational
        coefficients.

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = (1/17)*x^19 - (2/3)*y + 1/3; f
            1/17*x^19 - 2/3*y + 1/3
            sage: f.numerator()
            3*x^19 - 34*y + 17
            sage: f == f.numerator()
            False
            sage: f.numerator().base_ring()
            Integer Ring

        We check that the computation of numerator and denominator
        is valid.

        ::

            sage: K=QQ['x,y']
            sage: f=K.random_element()
            sage: f.numerator() / f.denominator() == f
            True

        The following tests against a bug fixed in :trac:`11780`::

            sage: P.<foo,bar> = ZZ[]
            sage: Q.<foo,bar> = QQ[]
            sage: f = Q.random_element()
            sage: f.numerator().parent() is P
            True
        """
        if self.base_ring() == RationalField():
            #This part is for compatibility with the univariate case,
            #where the numerator of a polynomial over RationalField
            #is a polynomial over IntegerRing
            #
            # Trac ticket #11780: Create the polynomial ring over
            # the integers using the (cached) polynomial ring constructor:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            integer_polynomial_ring = PolynomialRing(ZZ,\
            self.parent().ngens(), self.parent().gens(), order =\
            self.parent().term_order())
            return integer_polynomial_ring(self * self.denominator())
        else:
            return self * self.denominator()


def unpickle_MPolynomial_libsingular(MPolynomialRing_libsingular R, d):
    """
    Deserialize an ``MPolynomial_libsingular`` object

    INPUT:

    - ``R`` - the base ring
    - ``d`` - a Python dictionary as returned by :meth:`MPolynomial_libsingular.dict`

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(x)) == x # indirect doctest
        True
    """
    cdef ring *r = R._ring
    cdef poly *m
    cdef poly *p
    cdef int _i, _e
    p = p_ISet(0,r)
    rChangeCurrRing(r)
    for mon,c in d.iteritems():
        m = p_Init(r)
        for i,e in mon.sparse_iter():
            _i = i
            if _i >= r.N:
                p_Delete(&p, r)
                p_Delete(&m, r)
                raise TypeError, "variable index too big"
            _e = e
            if _e <= 0:
                p_Delete(&p, r)
                p_Delete(&m, r)
                raise TypeError, "exponent too small"
            overflow_check(_e, r)
            p_SetExp(m, _i+1, _e, r)
        p_SetCoeff(m, sa2si(c, r), r)
        p_Setm(m, r)
        p = p_Add_q(p, m, r)
    return new_MP(R, p)


cdef inline poly *addwithcarry(poly *tempvector, poly *maxvector, int pos, ring *_ring):
    if p_GetExp(tempvector, pos, _ring) < p_GetExp(maxvector, pos, _ring):
      p_SetExp(tempvector, pos, p_GetExp(tempvector, pos, _ring)+1, _ring)
    else:
      p_SetExp(tempvector, pos, 0, _ring)
      tempvector = addwithcarry(tempvector, maxvector, pos + 1, _ring)
    p_Setm(tempvector, _ring)
    return tempvector


cdef inline MPolynomial_libsingular new_MP(MPolynomialRing_libsingular parent, poly *juice):
    """
    Construct MPolynomial_libsingular from parent and SINGULAR poly.

    INPUT:

    - ``parent`` -- a :class:`MPolynomialRing_libsingular``
      instance. The parent of the polynomial to create.

    - ``juice`` -- a pointer to a Singular ``poly`` C struct.

    OUTPUT:

    A Python object :class:`MPolynomial_libsingular`.

    The ownership of ``juice`` will be transferred to the Python
    object. You must not free it yourself. Singular will modify the
    polynomial, so it is your repsonsiblity to make a copy if the
    Singular data structure is used elsewhere.
    """
    cdef MPolynomial_libsingular p = MPolynomial_libsingular.__new__(MPolynomial_libsingular)
    p._parent = parent
    p._parent_ring = singular_ring_reference(parent._ring)
    p._poly = juice
    p_Normalize(p._poly, p._parent_ring)
    return p


cdef poly *MPolynomial_libsingular_get_element(object self):
    return (<MPolynomial_libsingular>self)._poly
