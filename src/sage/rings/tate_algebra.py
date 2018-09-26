from sage.structure.factory import UniqueFactory
from sage.monoids.monoid import Monoid_class
from sage.rings.ring import CommutativeAlgebra
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.padics.padic_generic import pAdicGeneric

from sage.categories.monoids import Monoids
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from sage.structure.category_object import normalize_names
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.tate_algebra_element import TateAlgebraTerm
from sage.rings.tate_algebra_element import TateAlgebraElement

from sage.rings.polynomial.polydict import ETuple


# Factory
#########

class TateAlgebraFactory(UniqueFactory):
    """
    TODO Doctest here.
    """
    def create_key(self, base, prec=None, log_radii=ZZ(0), names=None, order='degrevlex'):
        if not isinstance(base, pAdicGeneric):
            raise TypeError("The base ring must be a p-adic field")
        # TODO: allow for arbitrary CDVF
        base = base.fraction_field()
        if names is None:
            raise ValueError("You must specify the names of the variables")
        names = normalize_names(-1, names)
        ngens = len(names)
        if not isinstance(log_radii, (list, tuple)):
            log_radii = [ZZ(log_radii)] * ngens
        elif len(log_radii) != ngens:
            raise ValueError("The number of radii does not match the number of variables")
        else:
            log_radii = [ ZZ(r) for r in log_radii ]
        order = TermOrder(order, ngens)
        if prec is None:
            prec = base.precision_cap()
        key = (base, prec, tuple(log_radii), names, order)
        return key

    def create_object(self, version, key):
        (base, prec, log_radii, names, order) = key
        return TateAlgebra_generic(base, prec, log_radii, names, order)

TateAlgebra = TateAlgebraFactory("TateAlgebra")

# Helper function
#################

def _are_parameters_compatible(params1,params2):
    r"""
    Test whether two sets of parameters of Tate algebra are compatible for coercion

    Each set of parameters comprises of a base ring, a list of variable names, a
    tuple of log of convergence radii and a term order. 
    
    INPUT:

    - ``params1`` - a tuple ``(base1,names1,log_radii1,order1)`` where `base1`
      is a discrete valuation ring, `names1` is a list of variable names,
      `log_radii1` is a tuple of integers, and `order1` is a term order

    - ``params2`` - a tuple ``(base2,names2,log_radii2,order2)`` where `base2`
      is a discrete valuation ring, `names2` is a list of variable names,
      `log_radii2` is a tuple of integers, and `order2` is a term order

    OUTPUT:

    - ``True`` if and only it coercion is possible from an algebra with the
      second set of parameters into an algebra with the first set of
      parameters. This is the second ring can be coerced into the first, if
      after this conversion, all convergence radii in the second set are the
      same as in the first set, and if the variable names and term orders are
      equal.
    
    """
    base1,names1,logs1,order1 = params1
    base2,names2,logs2,order2 = params2
    if base1.has_coerce_map_from(base2) and names1 == names2 and order1 == order2:
        ratio = base1.absolute_e() // base2.absolute_e()
        for i in range(len(logs1)) :
            if logs1[i] != ratio * logs2[i]:
                return False
        return True
    return False

# Parent for terms
##################

class TateTermMonoid(Monoid_class):
    r"""
    A base class for Tate algebra terms
    
    A term in `R\{X_1,\dots,X_n\}` is the product of a coefficient in `R` and a
    monomial in the variables `X_1,\dots,X_n`.
    
    Those terms form a pre-ordered monoid, with term multiplication and the
    term order of the parent Tate algebra.
    
    INPUT:
    
    - ``A`` - the Tate series algebra where the terms live
    
    EXAMPLES::
    
        sage: R = pAdicRing(2,prec=10,print_mode='digits')
        sage: A.<x,y> = TateAlgebra(R,log_radii=[1,1],order="lex")
        sage: T1 = A.monoid_of_terms(); T1
        Monoid of terms in x (val >= -1), y (val >= -1) over 2-adic Field with capped relative precision 10

    One may also create a Tate term monoid in the following way, but this is
    not recommended.
        
        sage: from sage.rings.tate_algebra import *
        sage: T2 = TateTermMonoid(A)

    """
    
    def __init__(self, A):
        r"""
        A base class for Tate algebra terms

        A term in `R\{X_1,\dots,X_n\}` is the product of a coefficient in `R` and a
        monomial in the variables `X_1,\dots,X_n`.

        Those terms form a pre-ordered monoid, with term multiplication and the
        term order of the parent Tate algebra.

        INPUT:

        - ``A`` - the Tate series algebra where the terms live

        EXAMPLES::

            sage: R = pAdicRing(2,prec=10,print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R,log_radii=(1,1),order="lex")
            sage: T1 = A.monoid_of_terms(); T1
            Monoid of terms in x (val >= -1), y (val >= -1) over 2-adic Field with capped relative precision 10

        One may also create a Tate term monoid in the following way, but this is
        not recommended.
        
            sage: from sage.rings.tate_algebra import *
            sage: T2 = TateTermMonoid(A)
        

        """
        # This function is not exposed to the user
        # so we do not check the inputs
        # TODO : At the moment it *is* exposed
        self.element_class = TateAlgebraTerm
        names = A.variable_names()
        Monoid_class.__init__(self, names)
        self._base = A.base_ring()
        self._names = names
        self._ngens = len(self._names)
        self._log_radii = ETuple(A.log_radii())
        self._order = A.term_order()
        self._parent_alg = A

    def _repr_(self):
        r"""
        Return a printable representation of a Tate term monoid

        EXAMPLES::

            sage: R = pAdicRing(2,prec=10,print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R,log_radii=[1,1],order="lex")
            sage: A.monoid_of_terms() # indirect doctest
            Monoid of terms in x (val >= -1), y (val >= -1) over 2-adic Field with capped relative precision 10

        """
        if self._ngens == 0:
            return "Monoid of terms over %s" % self._base
        vars = ""
        for i in range(self._ngens):
            vars += ", %s (val >= %s)" % (self._names[i], -self._log_radii[i])
        return "Monoid of terms in %s over %s" % (vars[2:], self._base)

    def _coerce_map_from_(self, R):
        r"""
        Test whether the term monoid `R` can be coerced into the term monoid

        `R` can be coerced into the term monoid if and only if the same coercion
        can be done on their parent algebras.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10

        A ring can be coerced into a monoid of terms if and only if it can be
        coerced into its base ring.
        
            sage: B = ZZ; B
            Integer Ring
            sage: T.has_coerce_map_from(B) # indirect doctest
            True
            sage: B = GF(2); B
            Finite Field of size 2
            sage: T.has_coerce_map_from(B) # indirect doctest
            False

        The base rings must be coercible.
        
            sage: S.<a> = Zq(4); S
            2-adic Unramified Extension Ring in a defined by x^2 + x + 1
            sage: B.<x,y> = TateAlgebra(S); B
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Unramified Extension Field in a defined by x^2 + x + 1
            sage: U = B.monoid_of_terms()
            sage: U.has_coerce_map_from(T) # indirect doctest
            True
            sage: T.has_coerce_map_from(U) # indirect doctest
            False

        Note that a Tate algebra cannot be coerced into a monoid of terms:

            sage: U.has_coerce_map_from(A) # indirect doctest
            False
            sage: T.has_coerce_map_from(B) # indirect doctest
            False

        Variable names must match exactly:
        
            sage: B.<x,z> = TateAlgebra(R)
            sage: U = B.monoid_of_terms()
            sage: T.has_coerce_map_from(U) # indirect doctest
            False
            sage: U.has_coerce_map_from(T) # indirect doctest
            False

        Term orders must match exactly:
        
            sage: B.<x,y> = TateAlgebra(R,order="lex"); B
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: U = B.monoid_of_terms()
            sage: T.has_coerce_map_from(U) # indirect doctest
            False
            sage: U.has_coerce_map_from(T) # indirect doctest
            False

        Variable names must also be in the same order:
        
            sage: B.<y,x> = TateAlgebra(R); B
            Tate Algebra in y (val >= 0), x (val >= 0) over 2-adic Field with capped relative precision 10
            sage: U = B.monoid_of_terms()
            sage: T.has_coerce_map_from(U) # indirect doctest
            False
            sage: U.has_coerce_map_from(T) # indirect doctest
            False

        """
        base = self._base
        if base.has_coerce_map_from(R):
            return True
        if isinstance(R, TateTermMonoid):
            params_self = (base,self._names,self._log_radii,self._order)
            params_R = (R.base_ring(),R.variable_names(),R.log_radii(),R.term_order())
            return _are_parameters_compatible(params_self,params_R)

    def algebra_of_series(self):
        r"""
        Return the Tate algebra corresponding to the Tate term monoid

        EXAMPLES::

            sage: R = Zp(2,10,print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: AA = T.algebra_of_series()
            sage: A = AA
        
        """
        return self._parent_alg    
            
    def base_ring(self):
        r"""
        Return the base ring of the Tate term monoid

        EXAMPLES::

            sage: R = Zp(2,10,print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: T.base_ring()
            2-adic Field with capped relative precision 10

        """
        
        return self._base

    def variable_names(self):
        r"""
        Return the names of the variables of the Tate term monoid

        EXAMPLES::

            sage: R = Zp(2,10,print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: T.variable_names()
            ('x', 'y')

        """
        return self._names

    def log_radii(self):
        r"""
        Return the log radii of convergence of the Tate term monoid

        EXAMPLES::

            sage: R = Zp(2,10,print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: T.log_radii()
            (0, 0)

        """
        return tuple(self._log_radii)

    def term_order(self):
        r"""
        Return the term order on the Tate term monoid

        EXAMPLES::

            sage: R = Zp(2,10,print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: T.term_order()
            Degree reverse lexicographic term order
        
        """
        return self._order

    def ngens(self):
        r"""
        Return the number of variables in the Tate term monoid

        EXAMPLES::

            sage: R = Zp(2,10,print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: T.ngens()
            2

        """
        return self._ngens


# Tate algebras
###############

class TateAlgebra_generic(CommutativeAlgebra):
    r"""
    Create a Tate series ring over a given complete discrete valuation
    ring.

    Given a complete discrete valuation ring `R`, variables `X_1,\dots,X_k`
    and convergence radii `r_,\dots, r_n` in `\mathbb{R}_{>0}`, the Tate
    algebra `R{X_1,\dots,X_k}` is the algebra of power series with
    coefficients `a_{i_1,\dots,i_n}` in `R` and such that
    `|a_{i_1,\dots,i_n}|*r_1^{-i_1}*\dots*r_n^{-i_n}` tends to 0 as
    `i_1,\dots,i_n` go towards infinity.


    INPUT:

    - ``base`` - a complete discrete valuation ring or field

    - ``names`` - names of the indeterminates

    - ``log-radii`` - (default: ``0``) the value(s) `\log(r_i)`. If only
      one number l is given, all `r_i`'s are defined with `\log(r_i)=l`.

    - ``prec`` - the default precision used if an exact object
      must be changed to an approximate object in order to do an
      arithmetic operation. If left as ``None``, it will be set to
      the precision of the base ring, if any. Otherwise,
      it will be set to the cap in precision of the base
      ring, if any. Otherwise, it will be set to the global
      default (20).

    - ``order`` - (default: ``degrevlex``) the monomial ordering 
      used to break ties when comparing terms with the same 
      coefficient valuation

    EXAMPLES::

        sage: R = Zp(2, 10, print_mode='digits'); R
        2-adic Ring with capped relative precision 10

    ::

        sage: A.<x,y> = TateAlgebra(R, order='lex'); A
        Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10

    The term ordering is used to determine how series are displayed. Terms
    are compared first according to the valuation of their coefficient, and
    ties are broken using the monomial ordering.

        sage: A.term_order()
        Lexicographic term order
        sage: f = 2 + y^5 + x^2; f
        (...0000000001)*x^2 + (...0000000001)*y^5 + (...00000000010)
        sage: B.<x,y> = TateAlgebra(R); B
        Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
        sage: B.term_order()
        Degree reverse lexicographic term order
        sage: B(f)
        (...0000000001)*y^5 + (...0000000001)*x^2 + (...00000000010)

    By order of priority, the precision cap is taken to be the given value, the
    precision of the base ring, the precision cap of the base ring, or the
    global default (20).  With a given value:

        sage: A.<x,y> = TateAlgebra(R,prec=5); A
        Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
        sage: A.precision_cap()
        5

    With the precision of the base ring:

        sage: A.<x,y> = TateAlgebra(R); A
        Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
        sage: A.precision_cap()
        10

    """
    def __init__(self, field, prec, log_radii, names, order, integral=False):
        """
        Initialize the Tate algebra

        TESTS::
        
        """
        self.element_class = TateAlgebraElement
        self._field = field
        self._cap = prec
        self._log_radii = ETuple(log_radii)  # TODO: allow log_radii in QQ (but ETuple does not support this)
        self._names = names
        self._ngens = len(names)
        self._order = order
        self._integral = integral
        if integral:
            base = field.integer_ring()
        else:
            base = field
        CommutativeAlgebra.__init__(self, base, names, category=CommutativeAlgebras(base))
        self._polynomial_ring = PolynomialRing(field, names, order=order)
        one = field(1)
        if integral:
            self._gens = [ self((one << log_radii[i].ceil()) * self._polynomial_ring.gen(i)) for i in range(self._ngens) ]
            self._integer_ring = self
        else:
            self._gens = [ self(g) for g in self._polynomial_ring.gens() ]
            self._integer_ring = TateAlgebra_generic(field, prec, log_radii, names, order, integral=True)
        self._parent_terms = TateTermMonoid(self)
        self._oneterm = self._parent_terms(one, ETuple([0]*self._ngens))

    def _an_element_(self):
        r"""
        Return an element of the Tate series algebra

        EXAMPLES::
         
        """
        return self.element_class(0)

    def _coerce_map_from_(self, R):
        r"""
        Test whether the ring `R` can be coerced into the Tate algebra.

        R can be coerced into the Tate algebra if it can be coerced into its
        base ring, or if R is a Tate algebra whose base ring can be coerced into
        the base ring, which has the same variables with the same monomial
        order, and if the series of R are converging on a larger ball than those
        of the algebra.

        INPUT:

        - ``R`` - the ring to be coerced

        EXAMPLES::
        
            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10

        B can be coerced into A if B can be coerced into the base ring of A

            sage: B = ZZ; B
            Integer Ring
            sage: A.has_coerce_map_from(B) # indirect doctest
            True
            sage: B = GF(2); B
            Finite Field of size 2
            sage: A.has_coerce_map_from(B) # indirect doctest
            False

        B can be coerced into A if B is a Tate algebra over a base ring which
        can be coerced into the base ring of A, with the same variables and term
        orders, and a larger convergence radius:

            sage: S.<a> = Zq(4); S
            2-adic Unramified Extension Ring in a defined by x^2 + x + 1
            sage: B.<x,y> = TateAlgebra(S); B
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Unramified Extension Field in a defined by x^2 + x + 1
            sage: B.has_coerce_map_from(A) # indirect doctest
            True

        If the base ring of B cannot be coerced into the base ring of A,
        coercion cannot happen:
        
            sage: A.has_coerce_map_from(B) # indirect doctest
            False

        If B has different variables than A, coercion cannot happen:

            sage: B.<x,z> = TateAlgebra(R)
            sage: A.has_coerce_map_from(B) # indirect doctest
            False

        If B has a different term order than A, coercion cannot happen:

            sage: B.<x,y> = TateAlgebra(R,order="lex"); B
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: B.has_coerce_map_from(A) # indirect doctest
            False
            sage: B.<y,x> = TateAlgebra(R); B
            Tate Algebra in y (val >= 0), x (val >= 0) over 2-adic Field with capped relative precision 10
            sage: B.has_coerce_map_from(A) # indirect doctest
            False

        """
        base = self._base
        if base.has_coerce_map_from(R):
            return True
        if isinstance(R, (TateTermMonoid, TateAlgebra_generic)):
            params_self = (base,self._names,self._log_radii,self._order)
            params_R = (R.base_ring(),R.variable_names(),R.log_radii(),R.term_order())
            return _are_parameters_compatible(params_self,params_R)

    def _ideal_class_(self, n):
        r"""
        Return the class of ideals in the Tate algebra

        INPUT:

        - ``n`` - number of generators

        EXAMPLE::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A._ideal_class_(3)
            <class 'sage.rings.tate_algebra_ideal.TateAlgebraIdeal'>
        
        .. NOTE::

            The argument ``n`` is disregarded in the current implementation.
        """
        from sage.rings.tate_algebra_ideal import TateAlgebraIdeal
        return TateAlgebraIdeal

    #def _pushout_(self, R):
    #    if not isinstance(R, TateAlgebra_generic):  # should we allow PolynomialRing as well?
    #        return None
    #    if self._names != R.variable_names():
    #        return None
    #    if self._order != R.term_order():
    #        return None
    #    Sbase = self._base
    #    Rbase = R.base_ring()
    #    base = pushout(Sbase, Rbase)
    #    Sratio = base.absolute_e() // Sbase.absolute_e()
    #    Rratio = base.absolute_e() // Rbase.absolute_e()
    #    log_radii = tuple([ min(self._log_radii[i] * Sratio, R.log_radii()[i] * Rratio) for i in range(self._ngens) ])
    #    cap = min(self._cap * Sratio, R.precision_cap() * Rratio)
    #    return TateAlgebra_generic(base, self._names, log_radii, cap, self._order)

    def gen(self, n=0):
        r"""
        Returns the ``n``'th generator of the algebra

        INPUT:

        - ``n`` - (default: ``0``) the generator to return

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A.gen()
            (...0000000001)*x
            sage: A.gen(0)
            (...0000000001)*x
            sage: A.gen(1)
            (...0000000001)*y
            sage: A.gen(2)
            Traceback (most recent call last):
            ...
            ValueError: Generator not defined
        
        """
        try:
            return self._gens[n]
        except IndexError:
            raise ValueError("Generator not defined")

    def gens(self):
        r"""
        Returns the list of generators of the algebra.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A.gens()
            ((...0000000001)*x, (...0000000001)*y)
        
        """
        return tuple(self._gens)

    def ngens(self):
        """
        Returns the number of generators of this algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A.ngens()
            2

        """
        return self._ngens

    def _repr_(self):
        """
        Returns a printable representation of this algebra.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            
        """
        vars = ""
        for i in range(self._ngens):
            vars += ", %s (val >= %s)" % (self._names[i], -self._log_radii[i])
        if self._integral:
            return "Integer ring of the Tate Algebra in %s over %s" % (vars[2:], self._field)
        else:
            return "Tate Algebra in %s over %s" % (vars[2:], self._field)

    def variable_names(self):
        """
        Return the list of the names of the variables of the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A.variable_names()
            ('x', 'y')
        
        """
        return self._names

    def log_radii(self):
        """
        Returns the list of the logs of the convergence radii of the series of the
        algebra.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A.log_radii()
            (0, 0)
            sage: B.<x,y> = TateAlgebra(R,log_radii=1); B
            Tate Algebra in x (val >= -1), y (val >= -1) over 2-adic Field with capped relative precision 10
            sage: B.log_radii()
            (1, 1)
            sage: C.<x,y> = TateAlgebra(R,log_radii=(1,-1)); C
            Tate Algebra in x (val >= -1), y (val >= 1) over 2-adic Field with capped relative precision 10
            sage: C.log_radii()
            (1, -1)

        """
        return self._log_radii

    def integer_ring(self):
        return self._integer_ring

    def monoid_of_terms(self):
        return self._parent_terms

    def term_order(self):
        """
        Returns the monomial order used in this algebra.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
        """
        return self._order

    def precision_cap(self):
        """
        Return the precision cap of the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A.precision_cap()
            10
        """
        return self._cap

    def characteristic(self):
        """
        Returns the characteristic of this algebra.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: A.characteristic()
            0
        """
        return self.base_ring().characteristic()

    #def ideal(self, gens):
    #    pass


