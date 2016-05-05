r"""
Noncommutative Polynomials via libSINGULAR/Plural

This module provides specialized and optimized implementations for
noncommutative multivariate polynomials over many coefficient rings, via the
shared library interface to SINGULAR. In particular, the following coefficient
rings are supported by this implementation:

- the rational numbers `\QQ`, and

- finite fields `\GF{p}` for `p` prime

AUTHORS:

The PLURAL wrapper is due to

  - Burcin Erocal (2008-11 and 2010-07): initial implementation and concept

  - Michael Brickenstein (2008-11 and 2010-07): initial implementation and concept

  - Oleksandr Motsak (2010-07): complete overall noncommutative functionality and first release

  - Alexander Dreyer (2010-07): noncommutative ring functionality and documentation

  - Simon King (2011-09): left and two-sided ideals; normal forms; pickling;
    documentation

The underlying libSINGULAR interface was implemented by

  - Martin Albrecht (2007-01): initial implementation

  - Joel Mohler (2008-01): misc improvements, polishing

  - Martin Albrecht (2008-08): added `\QQ(a)` and `\ZZ` support

  - Simon King (2009-04): improved coercion

  - Martin Albrecht (2009-05): added `\ZZ/n\ZZ` support, refactoring

  - Martin Albrecht (2009-06): refactored the code to allow better
    re-use

.. TODO::

    extend functionality towards those of libSINGULARs commutative part

EXAMPLES:

We show how to construct various noncommutative polynomial rings::

    sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
    sage: P.<x,y,z> = A.g_algebra(relations={y*x:-x*y}, order = 'lex')

    sage: P
    Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y}

    sage: y*x + 1/2
    -x*y + 1/2

    sage: A.<x,y,z> = FreeAlgebra(GF(17), 3)
    sage: P.<x,y,z> = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
    sage: P
    Noncommutative Multivariate Polynomial Ring in x, y, z over Finite Field of size 17, nc-relations: {y*x: -x*y}

    sage: y*x + 7
    -x*y + 7


Raw use of this class; *this is not the intended use!*
::

    sage: from sage.matrix.constructor import Matrix
    sage: c = Matrix(3)
    sage: c[0,1] = -2
    sage: c[0,2] = 1
    sage: c[1,2] = 1

    sage: d = Matrix(3)
    sage: d[0, 1] = 17
    sage: P = QQ['x','y','z']
    sage: c = c.change_ring(P)
    sage: d = d.change_ring(P)

    sage: from sage.rings.polynomial.plural import NCPolynomialRing_plural
    sage: R.<x,y,z> = NCPolynomialRing_plural(QQ, c = c, d = d, order=TermOrder('lex',3),category=Algebras(QQ))
    sage: R
    Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -2*x*y + 17}

    sage: R.term_order()
    Lexicographic term order

    sage: a,b,c = R.gens()
    sage: f = 57 * a^2*b + 43 * c + 1; f
    57*x^2*y + 43*z + 1

TESTS::

    sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
    sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
    sage: TestSuite(P).run()
    sage: loads(dumps(P)) is P
    True

"""
include "cysignals/memory.pxi"

from sage.categories.algebras import Algebras

# singular rings

from sage.libs.singular.function cimport RingWrap

from sage.libs.singular.polynomial cimport (singular_polynomial_call, singular_polynomial_cmp, singular_polynomial_add, singular_polynomial_sub, singular_polynomial_neg, singular_polynomial_pow, singular_polynomial_mul, singular_polynomial_rmul, singular_polynomial_deg, singular_polynomial_str_with_changed_varnames, singular_polynomial_latex, singular_polynomial_str, singular_polynomial_div_coeff)

import sage.libs.singular.ring
from sage.libs.singular.ring cimport singular_ring_new, singular_ring_delete, wrap_ring, singular_ring_reference

from sage.libs.singular.singular cimport si2sa, sa2si, overflow_check

from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import is_IntegerRing

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular
from sage.rings.polynomial.multi_polynomial_ring_generic import MPolynomialRing_generic
from sage.rings.polynomial.multi_polynomial_ideal import NCPolynomialIdeal

from sage.rings.polynomial.polydict import ETuple
from sage.rings.ring import check_default_category
from sage.structure.element cimport CommutativeRingElement, Element, ModuleElement
from sage.structure.factory import UniqueFactory
from sage.structure.parent cimport Parent
from sage.structure.parent_gens cimport ParentWithGens
from sage.rings.polynomial.term_order import TermOrder

class G_AlgFactory(UniqueFactory):
    """
    A factory for the creation of g-algebras as unique parents.

    TESTS::

        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: H = A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
        sage: H is A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y}) # indirect doctest
        True

    """
    def create_object(self, version, key, **extra_args):
        """
        Create a g-algebra to a given unique key.

        INPUT:

        - ``key`` - a 6-tuple, formed by a base ring, a tuple of names, two
          matrices over a polynomial ring over the base ring with the given
          variable names, a term order, and a category
        - ``extra_args`` - a dictionary, whose only relevant key is 'check'.

        TEST::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H=A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
            sage: sorted(H.relations().iteritems(),key=str)
            [(y*x, x*y - z), (z*x, x*z + 2*x), (z*y, y*z - 2*y)]

        """
        # key = (base_ring,names, c,d, order, category)
        # extra args: check
        base_ring,names,c,d,order,category = key
        check = extra_args.get('check')
        return NCPolynomialRing_plural(base_ring, names, c,d, order, category, check)
    def create_key_and_extra_args(self, base_ring, c,d, names=None, order=None,
                                 category=None,check=None):
        """
        Create a unique key for g-algebras.

        INPUT:

        - ``base_ring`` - a ring
        - ``c,d`` - two matrices
        - ``names`` - a tuple or list of names
        - ``order`` - (optional) term order
        - ``category`` - (optional) category
        - ``check`` - optional bool

        TEST::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
            sage: H is A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y}) # indirect doctest
            True

        """
        if names is None:
            raise ValueError("The generator names must be provided")

        # Get the number of names:
        names = tuple(names)
        n = len(names)
        if not isinstance(order, TermOrder):
            order = TermOrder(order or 'degrevlex', n)

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        P = PolynomialRing(base_ring, n, names, order=order)
        # The names may have been normalised in P:
        names = P.variable_names()
        c = c.change_ring(P)
        c.set_immutable()
        d = d.change_ring(P)
        d.set_immutable()

        # Get the correct category
        category=check_default_category(Algebras(base_ring),category)

        # Extra arg
        if check is None:
            return (base_ring,names,c,d,order,category),{}
        return (base_ring,names,c,d,order,category),{'check':check}

g_Algebra = G_AlgFactory('sage.rings.polynomial.plural.g_Algebra')

cdef class NCPolynomialRing_plural(Ring):
    """
    A non-commutative polynomial ring.

    EXAMPLES::

        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: H = A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
        sage: H._is_category_initialized()
        True
        sage: H.category()
        Category of algebras over Rational Field
        sage: TestSuite(H).run()

    Note that two variables commute if they are not part of the given
    relations::

        sage: H.<x,y,z> = A.g_algebra({z*x:x*z+2*x, z*y:y*z-2*y})
        sage: x*y == y*x
        True

    """
    def __init__(self, base_ring, names, c, d, order, category, check = True):
        """
        Construct a noncommutative polynomial G-algebra subject to the following conditions:

        INPUT:

        - ``base_ring`` - base ring (must be either `\GF{q}`, `\ZZ`, `\ZZ/n\ZZ`, `\QQ` or absolute number field)
        - ``names`` - a tuple of names of ring variables
        - ``c``, ``d``- upper triangular matrices of coefficients,
          resp. commutative polynomials, satisfying the nondegeneracy
          conditions, which are to be tested if ``check`` is ``True``. These
          matrices describe the noncommutative relations:

            ``self.gen(j)*self.gen(i) == c[i, j] * self.gen(i)*self.gen(j) + d[i, j],``

          where ``0 <= i < j < self.ngens()``. Note that two variables
          commute if they are not part of one of these relations.
        - ``order`` - term order
        - ``check`` - check the noncommutative conditions (default: ``True``)

        TESTS:

        It is strongly recommended to construct a g-algebra using
        :class:`G_AlgFactory`. The following is just for documenting
        the arguments of the ``__init__`` method::

            sage: from sage.matrix.constructor  import Matrix
            sage: c0 = Matrix(3)
            sage: c0[0,1] = -1
            sage: c0[0,2] = 1
            sage: c0[1,2] = 1

            sage: d0 = Matrix(3)
            sage: d0[0, 1] = 17
            sage: P = QQ['x','y','z']
            sage: c = c0.change_ring(P)
            sage: d = d0.change_ring(P)

            sage: from sage.rings.polynomial.plural import NCPolynomialRing_plural
            sage: P.<x,y,z> = NCPolynomialRing_plural(QQ, c = c, d = d, order=TermOrder('lex',3), category=Algebras(QQ))

            sage: P # indirect doctest
            Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y + 17}

            sage: P(x*y)
            x*y

            sage: f = 27/113 * x^2 + y*z + 1/2; f
            27/113*x^2 + y*z + 1/2

            sage: P.term_order()
            Lexicographic term order

            sage: from sage.rings.polynomial.plural import NCPolynomialRing_plural
            sage: P = GF(7)['x','y','z']
            sage: c = c0.change_ring(P)
            sage: d = d0.change_ring(P)
            sage: P.<x,y,z> = NCPolynomialRing_plural(GF(7), c = c, d = d, order=TermOrder('degrevlex',3), category=Algebras(GF(7)))

            sage: P # indirect doctest
            Noncommutative Multivariate Polynomial Ring in x, y, z over Finite Field of size 7, nc-relations: {y*x: -x*y + 3}

            sage: P(x*y)
            x*y

            sage: f = 3 * x^2 + y*z + 5; f
            3*x^2 + y*z - 2

            sage: P.term_order()
            Degree reverse lexicographic term order
        """
        n = len(names)
        self._relations = None

        P = c.base_ring()
        self._c = c
        self._d = d

        from sage.libs.singular.function import singular_function
        ncalgebra = singular_function('nc_algebra')

        cdef RingWrap rw = ncalgebra(self._c, self._d, ring = P)

        #       rw._output()
        self._ring = singular_ring_reference(rw._ring)
        self._ring.ShortOut = 0

        self.__ngens = n
        self.__term_order = order

        Ring.__init__(self, base_ring, names, category=category)
        self._populate_coercion_lists_()

        #MPolynomialRing_generic.__init__(self, base_ring, n, names, order)
        #self._has_singular = True
        assert(n == len(self._names))

        self._one_element = new_NCP(self, p_ISet(1, self._ring))
        self._zero_element = new_NCP(self, NULL)


        if check:
            from sage.libs.singular.function_factory import ff
            test = ff.nctools__lib.ndcond(ring = self)
            if (len(test) != 1) or (test[0] != 0):
                raise ValueError("NDC check failed!")

    def __reduce__(self):
        """
        TESTS::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
            sage: H is A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
            True
            sage: H is loads(dumps(H))  # indirect doctest
            True

        Check that :trac:`17224` is fixed::

            sage: from sage.rings.polynomial.term_order import TermOrder
            sage: F.<x,y> = FreeAlgebra(QQ)
            sage: g = F.g_algebra({y*x:-x*y}, order=TermOrder('wdegrevlex', [1,2]))
            sage: loads(dumps(g)) == g
            True
        """
        return g_Algebra, (self.base_ring(),self._c,self._d,
                            self.variable_names(),
                            self.term_order(),
                            self.category())

    def __dealloc__(self):
        r"""
        Carefully deallocate the ring, without changing "currRing"
        (since this method can be at unpredictable times due to garbage
        collection).

        TESTS:

        This example caused a segmentation fault with a previous version
        of this method::

            sage: import gc
            sage: from sage.rings.polynomial.plural import NCPolynomialRing_plural
            sage: from sage.algebras.free_algebra import FreeAlgebra
            sage: A1.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: R1 = A1.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y}, order=TermOrder('degrevlex', 2))
            sage: A2.<x,y,z> = FreeAlgebra(GF(5), 3)
            sage: R2 = A2.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y}, order=TermOrder('degrevlex', 2))
            sage: A3.<x,y,z> = FreeAlgebra(GF(11), 3)
            sage: R3 = A3.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y}, order=TermOrder('degrevlex', 2))
            sage: A4.<x,y,z> = FreeAlgebra(GF(13), 3)
            sage: R4 = A4.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y}, order=TermOrder('degrevlex', 2))
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
        singular_ring_delete(self._ring)

    def _element_constructor_(self, element):
        """
        Make sure element is a valid member of self, and return the constructed element.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')

        We can construct elements from the base ring::

            sage: P(1/2)
            1/2

        and all kinds of integers::

            sage: P(17)
            17
            sage: P(int(19))
            19
            sage: P(long(19))
            19

        TESTS:

        Check conversion from self::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P.<x,y,z> = A.g_algebra(relations={y*x:-x*y}, order = 'lex')

            sage: P._element_constructor_(1/2)
            1/2

            sage: P._element_constructor_(x*y)
            x*y

            sage: P._element_constructor_(y*x)
            -x*y

        Testing special cases::

            sage: P._element_constructor_(1)
            1

            sage: P._element_constructor_(0)
            0
        """

        if element == 0:
            return self._zero_element
        if element == 1:
            return self._one_element

        cdef poly *_p
        cdef ring *_ring,
        cdef number *_n

        _ring = self._ring

        base_ring = self.base_ring()

        if(_ring != currRing): rChangeCurrRing(_ring)


        if isinstance(element, NCPolynomial_plural):

            if element.parent() is <object>self:
                return element
            elif element.parent() == self:
                # is this safe?
                _p = p_Copy((<NCPolynomial_plural>element)._poly, _ring)

        elif isinstance(element, CommutativeRingElement):
            # base ring elements
            if  <Parent>element.parent() is base_ring:
                # shortcut for GF(p)
                if isinstance(base_ring, FiniteField_prime_modn):
                    _p = p_ISet(int(element) % _ring.ch, _ring)
                else:
                    _n = sa2si(element,_ring)
                    _p = p_NSet(_n, _ring)

            # also accepting ZZ
            elif is_IntegerRing(element.parent()):
                if isinstance(base_ring, FiniteField_prime_modn):
                    _p = p_ISet(int(element),_ring)
                else:
                    _n = sa2si(base_ring(element),_ring)
                    _p = p_NSet(_n, _ring)
            else:
                # fall back to base ring
                element = base_ring._coerce_c(element)
                _n = sa2si(element,_ring)
                _p = p_NSet(_n, _ring)

        # Accepting int
        elif isinstance(element, int):
            if isinstance(base_ring, FiniteField_prime_modn):
                _p = p_ISet(int(element) % _ring.ch,_ring)
            else:
                _n = sa2si(base_ring(element),_ring)
                _p = p_NSet(_n, _ring)

        # and longs
        elif isinstance(element, long):
            if isinstance(base_ring, FiniteField_prime_modn):
                element = element % self.base_ring().characteristic()
                _p = p_ISet(int(element),_ring)
            else:
                _n = sa2si(base_ring(element),_ring)
                _p = p_NSet(_n, _ring)

        else:
            raise NotImplementedError("not able to interprete "+repr(element) +
                                      " of type "+ repr(type(element)) +
                                      " as noncommutative polynomial")  ### ??????
        return new_NCP(self,_p)



    cpdef _coerce_map_from_(self, S):
       """
       The only things that coerce into this ring are:

       - the integer ring
       - other localizations away from fewer primes

       EXAMPLES::

           sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
           sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
           sage: P._coerce_map_from_(ZZ)
           True
       """

       if self.base_ring().has_coerce_map_from(S):
           return True



    def __hash__(self):
        """
        Return a hash for this noncommutative ring, that is, a hash of the string
        representation of this polynomial ring.

        NOTE:

        G-algebras are unique parents, provided that the g-algebra constructor
        is used. Thus, the hash simply is the memory address of the g-algebra
        (so, it is a session hash, but no stable hash). It is possible to
        destroy uniqueness of g-algebras on purpose, but that's your own
        problem if you do those things.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: {P:2}[P]            # indirect doctest
            2

        """
        return <Py_hash_t> <void *> self

    def __pow__(self, n, _):
        """
        Return the free module of rank `n` over this ring.

        NOTE:

        This is not properly implemented yet. Thus, there is
        a warning.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P.<x,y,z> = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: P^3
            d...: UserWarning: You are constructing a free module
            over a noncommutative ring. Sage does not have a concept
            of left/right and both sided modules, so be careful.
            It's also not guaranteed that all multiplications are
            done from the right side.
            d...: UserWarning: You are constructing a free module
            over a noncommutative ring. Sage does not have a concept
            of left/right and both sided modules, so be careful.
            It's also not guaranteed that all multiplications are
            done from the right side.
            Ambient free module of rank 3 over Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y}

        """
        from sage.modules.all import FreeModule
        return FreeModule(self, n)

    def term_order(self):
        """
        Return the term ordering of the noncommutative ring.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: P.term_order()
            Lexicographic term order

            sage: P = A.g_algebra(relations={y*x:-x*y})
            sage: P.term_order()
            Degree reverse lexicographic term order
        """
        return self.__term_order

    def is_commutative(self):
        """
        Return ``False``.

        .. TODO:: Provide a mathematically correct answer.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: P.is_commutative()
            False
        """
        return False

    def is_field(self):
        """
        Return ``False``.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: P.is_field()
            False
        """
        return False

    def _repr_(self):
        """
        EXAMPLE::

            sage: A.<x,y> = FreeAlgebra(QQ, 2)
            sage: H.<x,y> = A.g_algebra({y*x:-x*y})
            sage: H # indirect doctest
            Noncommutative Multivariate Polynomial Ring in x, y over Rational Field, nc-relations: {y*x: -x*y}
            sage: x*y
            x*y
            sage: y*x
            -x*y
        """
#TODO: print the relations
        varstr = ", ".join([ rRingVar(i,self._ring)  for i in range(self.__ngens) ])
        return "Noncommutative Multivariate Polynomial Ring in %s over %s, nc-relations: %s"%(varstr,self.base_ring(), self.relations())


    def _ringlist(self):
        """
        Return an internal list representation of the noncommutative ring.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: P._ringlist()
            [
                                                                       [ 0 -1  1]
                                                                       [ 0  0  1]
            0, ['x', 'y', 'z'], [['lp', (1, 1, 1)], ['C', (0,)]], [0], [ 0  0  0],
            <BLANKLINE>
            [0 0 0]
            [0 0 0]
            [0 0 0]
            ]
        """
        cdef ring* _ring = self._ring
        if(_ring != currRing): rChangeCurrRing(_ring)
        from sage.libs.singular.function import singular_function
        ringlist = singular_function('ringlist')
        result = ringlist(self, ring=self)
        return result


    def relations(self, add_commutative = False):
        """
        Return the relations of this g-algebra.

        INPUT:

        ``add_commutative`` (optional bool, default ``False``)

        OUTPUT:

        The defining relations. There are some implicit relations:
        Two generators commute if they are not part of any given
        relation. The implicit relations are not provided, unless
        ``add_commutative==True``.

        EXAMPLE::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H.<x,y,z> = A.g_algebra({z*x:x*z+2*x, z*y:y*z-2*y})
            sage: x*y == y*x
            True
            sage: H.relations()
            {z*x: x*z + 2*x, z*y: y*z - 2*y}
            sage: H.relations(add_commutative=True)
            {y*x: x*y, z*x: x*z + 2*x, z*y: y*z - 2*y}

        """
        if add_commutative:
            if self._relations_commutative is not None:
                return self._relations_commutative

            from sage.algebras.free_algebra import FreeAlgebra
            A = FreeAlgebra( self.base_ring(), self.ngens(), self.gens() )

            res = {}
            n = self.ngens()
            for r in range(0, n-1, 1):
                for c in range(r+1, n, 1):
                    res[ A.gen(c) * A.gen(r) ] = self.gen(c) * self.gen(r) # C[r, c] * P.gen(r) * P.gen(c) + D[r, c]
            self._relations_commutative = res
            return res

        if self._relations is not None:
            return self._relations

        from sage.algebras.free_algebra import FreeAlgebra
        A = FreeAlgebra( self.base_ring(), self.ngens(), self.gens() )

        res = {}
        n = self.ngens()
        for r in range(0, n-1, 1):
            for c in range(r+1, n, 1):
                if  (self.gen(c) * self.gen(r) != self.gen(r) * self.gen(c)):
                    res[ A.gen(c) * A.gen(r) ] = self.gen(c) * self.gen(r) # C[r, c] * P.gen(r) * P.gen(c) + D[r, c]

        self._relations = res
        return self._relations

    def ngens(self):
        """
        Returns the number of variables in this noncommutative polynomial ring.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P.<x,y,z> = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: P.ngens()
            3
        """
        return int(self.__ngens)

    def gen(self, int n=0):
        """
        Returns the ``n``-th generator of this noncommutative polynomial
        ring.

        INPUT:

        - ``n`` -- an integer ``>= 0``

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
            sage: P.gen(),P.gen(1)
            (x, y)

        Note that the generators are not cached::

            sage: P.gen(1) is P.gen(1)
            False

        """
        cdef poly *_p
        cdef ring *_ring = self._ring

        if n < 0 or n >= self.__ngens:
            raise ValueError("Generator not defined.")

        rChangeCurrRing(_ring)
        _p = p_ISet(1,_ring)
        p_SetExp(_p, n+1, 1, _ring)
        p_Setm(_p, _ring);

        return new_NCP(self,_p)

    def ideal(self, *gens, **kwds):
        """
        Create an ideal in this polynomial ring.

        INPUT:

        - ``*gens`` - list or tuple of generators (or several input arguments)
        - ``coerce`` - bool (default: ``True``); this must be a
          keyword argument. Only set it to ``False`` if you are certain
          that each generator is already in the ring.
        - ``side`` - string (either "left", which is the default, or "twosided")
          Must be a keyword argument. Defines whether the ideal is a left ideal
          or a two-sided ideal. Right ideals are not implemented.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P.<x,y,z> = A.g_algebra(relations={y*x:-x*y}, order = 'lex')

            sage: P.ideal([x + 2*y + 2*z-1, 2*x*y + 2*y*z-y, x^2 + 2*y^2 + 2*z^2-x])
            Left Ideal (x + 2*y + 2*z - 1, 2*x*y + 2*y*z - y, x^2 - x + 2*y^2 + 2*z^2) of Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y}
            sage: P.ideal([x + 2*y + 2*z-1, 2*x*y + 2*y*z-y, x^2 + 2*y^2 + 2*z^2-x], side="twosided")
            Twosided Ideal (x + 2*y + 2*z - 1, 2*x*y + 2*y*z - y, x^2 - x + 2*y^2 + 2*z^2) of Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y}

        """
        coerce = kwds.get('coerce', True)
        if len(gens) == 1:
            gens = gens[0]
        #if is_SingularElement(gens):
        #    gens = list(gens)
        #    coerce = True
        #elif is_Macaulay2Element(gens):
        #    gens = list(gens)
        #    coerce = True
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return NCPolynomialIdeal(self, gens, coerce=False, side=kwds.get('side','left'))

    def _list_to_ring(self, L):
        """
        Convert internal list representation to noncommutative ring.

        EXAMPLES::

           sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
           sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
           sage: rlist = P._ringlist();
           sage: Q = P._list_to_ring(rlist)
           sage: Q # indirect doctest
           <noncommutative RingWrap>
        """

        cdef ring* _ring = self._ring
        if(_ring != currRing): rChangeCurrRing(_ring)

        from sage.libs.singular.function import singular_function
        ring = singular_function('ring')
        return ring(L, ring=self)

# TODO: Implement this properly!
#    def quotient(self, I):
#        """
#        Construct quotient ring of ``self`` and the two-sided Groebner basis of `ideal`
#
#        EXAMPLE::
#
#            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
#            sage: H = A.g_algebra(relations={y*x:-x*y},  order='lex')
#            sage: I = H.ideal([H.gen(i) ^2 for i in [0, 1]]).twostd()
#
#            sage: Q = H.quotient(I); Q
#            Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y}
#
#        TESTS:
#
#        check coercion bug::
#
#            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
#            sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
#            sage: rlist = P._ringlist();
#            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
#            sage: H = A.g_algebra(relations={y*x:-x*y},  order='lex')
#            sage: I = H.ideal([H.gen(i) ^2 for i in [0, 1]]).twostd()
#            sage: Q = H.quotient(I); Q
#            Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y}
#            sage: Q.gen(0)^2
#            0
#            sage: Q.gen(1) * Q.gen(0)
#            -x*y
#        """
#        L = self._ringlist()
#        L[3] = I.twostd()
#        W = self._list_to_ring(L)
#        return new_NRing(W, self.base_ring())


    ### The following methods are handy for implementing Groebner
    ### basis algorithms. They do only superficial type/sanity checks
    ### and should be called carefully.

    def monomial_quotient(self, NCPolynomial_plural f, NCPolynomial_plural g, coeff=False):
        r"""
        Return ``f/g``, where both ``f`` and ``g`` are treated as
        monomials.

        Coefficients are ignored by default.

        INPUT:

        - ``f`` - monomial
        - ``g`` - monomial
        - ``coeff`` - divide coefficients as well (default: ``False``)

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_quotient(3/2*x*y,x,coeff=True)
            3/2*y

        Note that `\ZZ` behaves differently if ``coeff=True``::

            sage: P.monomial_quotient(2*x,3*x)
            1
            sage: P.monomial_quotient(2*x,3*x,coeff=True)
            2/3

        TESTS::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: R.inject_variables()
            Defining x, y, z

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_quotient(x*y,x)
            y

            sage: P.monomial_quotient(x*y,R.gen()) # not tested
            y

            sage: P.monomial_quotient(P(0),P(1))
            0

            sage: P.monomial_quotient(P(1),P(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

            sage: P.monomial_quotient(P(3/2),P(2/3), coeff=True)
            9/4

            sage: P.monomial_quotient(x,P(1))
            x

        TESTS::

            sage: P.monomial_quotient(x,y) # Note the wrong result
            x*y^...

        .. WARNING::

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

        if(r != currRing): rChangeCurrRing(r)

        if not f._poly:
            return self._zero_element
        if not g._poly:
            raise ZeroDivisionError

        res = pDivide(f._poly,g._poly)
        if coeff:
            if r.ringtype == 0 or r.cf.nDivBy(p_GetCoeff(f._poly, r), p_GetCoeff(g._poly, r)):
                n = r.cf.nDiv( p_GetCoeff(f._poly, r) , p_GetCoeff(g._poly, r))
                p_SetCoeff0(res, n, r)
            else:
                raise ArithmeticError("Cannot divide these coefficients.")
        else:
            p_SetCoeff0(res, n_Init(1, r), r)
        return new_NCP(self, res)

    def monomial_divides(self, NCPolynomial_plural a, NCPolynomial_plural b):
        """
        Return ``False`` if ``a`` does not divide ``b`` and ``True``
        otherwise.

        Coefficients are ignored.

        INPUT:

        - ``a`` -- monomial

        - ``b`` -- monomial

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_divides(x*y*z, x^3*y^2*z^4)
            True
            sage: P.monomial_divides(x^3*y^2*z^4, x*y*z)
            False

        TESTS::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: Q = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: Q.inject_variables()
            Defining x, y, z

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_divides(P(1), P(0))
            True
            sage: P.monomial_divides(P(1), x)
            True
        """
        cdef poly *_a
        cdef poly *_b
        cdef ring *_r
        if a._parent is not b._parent:
            b = (<NCPolynomialRing_plural>a._parent)._coerce_c(b)

        _a = a._poly
        _b = b._poly
        _r = (<NCPolynomialRing_plural>a._parent)._ring

        if _a == NULL:
            raise ZeroDivisionError
        if _b == NULL:
            return True

        if not p_DivisibleBy(_a, _b, _r):
            return False
        else:
            return True


    def monomial_lcm(self, NCPolynomial_plural f, NCPolynomial_plural g):
        """
        LCM for monomials. Coefficients are ignored.

        INPUT:

        - ``f`` - monomial

        - ``g`` - monomial

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_lcm(3/2*x*y,x)
            x*y

        TESTS::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: R.inject_variables()
            Defining x, y, z

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_lcm(x*y,R.gen()) # not tested
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
                raise ArithmeticError("Cannot compute LCM of zero and nonzero element.")
        if g._poly == NULL:
            raise ArithmeticError("Cannot compute LCM of zero and nonzero element.")

        if(self._ring != currRing): rChangeCurrRing(self._ring)

        pLcm(f._poly, g._poly, m)
        p_Setm(m, self._ring)
        return new_NCP(self,m)

    def monomial_reduce(self, NCPolynomial_plural f, G):
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

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, 1/2  ]
            sage: P.monomial_reduce(f,G)
            (y, 1/4*x*y + 2/7)

        TESTS::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: Q = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: Q.inject_variables()
            Defining x, y, z

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z
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
            if isinstance(g, NCPolynomial_plural) \
                   and (<NCPolynomial_plural>g) \
                   and p_LmDivisibleBy((<NCPolynomial_plural>g)._poly, m, r):
                flt = pDivide(f._poly, (<NCPolynomial_plural>g)._poly)
                #p_SetCoeff(flt, n_Div( p_GetCoeff(f._poly, r) , p_GetCoeff((<NCPolynomial_plural>g)._poly, r), r), r)
                p_SetCoeff(flt, n_Init(1, r), r)
                return new_NCP(self,flt), g
        return self._zero_element,self._zero_element

    def monomial_pairwise_prime(self, NCPolynomial_plural g, NCPolynomial_plural h):
        """
        Return ``True`` if ``h`` and ``g`` are pairwise prime.

        Both ``h`` and ``g`` are treated as monomials.

        Coefficients are ignored.

        INPUT:

        - ``h`` - monomial
        - ``g`` - monomial

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_pairwise_prime(x^2*z^3, y^4)
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, 3/4*y^3)
            False

        TESTS::

            sage: A.<x1,y1,z1> = FreeAlgebra(QQ, 3)
            sage: Q = A.g_algebra(relations={y1*x1:-x1*y1},  order='lex')
            sage: Q.inject_variables()
            Defining x1, y1, z1

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

            sage: P.monomial_pairwise_prime(x^2*z^3, x1^4) # not tested
            True

            sage: P.monomial_pairwise_prime((2)*x^3*y^2, Q.zero()) # not tested
            True

            sage: P.monomial_pairwise_prime(2*P.one(),x)
            False
        """
        cdef int i
        cdef ring *r
        cdef poly *p
        cdef poly *q

        if h._parent is not g._parent:
            g = (<NCPolynomialRing_plural>h._parent)._coerce_c(g)

        r = (<NCPolynomialRing_plural>h._parent)._ring
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

    def monomial_all_divisors(self, NCPolynomial_plural t):
        """
        Return a list of all monomials that divide ``t``.

        Coefficients are ignored.

        INPUT:

        - ``t`` - a monomial

        OUTPUT:

        a list of monomials

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, y, z

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
          M.append(new_NCP(self, p_Copy(tempvector,_ring)))
        return M

def unpickle_NCPolynomial_plural(NCPolynomialRing_plural R, d):
    """
    Auxiliary function to unpickle a non-commutative polynomial.

    TEST::

        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: H.<x,y,z> = A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
        sage: p = x*y+2*z+4*x*y*z*x
        sage: loads(dumps(p)) == p  # indirect doctest
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
                p_Delete(&p,r)
                p_Delete(&m,r)
                raise TypeError("variable index too big")
            _e = e
            if _e <= 0:
                p_Delete(&p,r)
                p_Delete(&m,r)
                raise TypeError("exponent too small")
            overflow_check(_e, r)
            p_SetExp(m, _i+1,_e, r)
        p_SetCoeff(m, sa2si(c, r), r)
        p_Setm(m,r)
        p = p_Add_q(p,m,r)
    return new_NCP(R,p)


cdef class NCPolynomial_plural(RingElement):
    """
    A noncommutative multivariate polynomial implemented using libSINGULAR.
    """
    def __init__(self, NCPolynomialRing_plural parent):
        """
        Construct a zero element in parent.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: from sage.rings.polynomial.plural import NCPolynomial_plural
            sage: NCPolynomial_plural(H)
            0
        """
        self._poly = NULL
        self._parent = parent

    def __dealloc__(self):
        # TODO: Warn otherwise!
        # for some mysterious reason, various things may be NULL in some cases
        if self._parent is not None and (<NCPolynomialRing_plural>self._parent)._ring != NULL and self._poly != NULL:
            p_Delete(&self._poly, (<NCPolynomialRing_plural>self._parent)._ring)

#    def __call__(self, *x, **kwds): # ?

    def __reduce__(self):
        """
        TEST::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H.<x,y,z> = A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
            sage: loads(dumps(x*y+2*z+4*x*y*z*x))
            4*x^2*y*z + 8*x^2*y - 4*x*z^2 + x*y - 8*x*z + 2*z

        """
        return unpickle_NCPolynomial_plural, (self._parent, self.dict())

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

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y

            sage: x == x
            True

            sage: x > y
            True
            sage: y^2 > x
            False

        TESTS::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y

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
        """
        if left is right:
            return 0
        cdef poly *p = (<NCPolynomial_plural>left)._poly
        cdef poly *q = (<NCPolynomial_plural>right)._poly
        cdef ring *r = (<NCPolynomialRing_plural>left._parent)._ring
        return singular_polynomial_cmp(p, q, r)

    cpdef ModuleElement _add_( left, ModuleElement right):
        """
        Adds left and right.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: 3/2*x + 1/2*y + 1 # indirect doctest
            3/2*x + 1/2*y + 1
        """
        cdef poly *_p
        singular_polynomial_add(&_p, left._poly,
                                 (<NCPolynomial_plural>right)._poly,
                                 (<NCPolynomialRing_plural>left._parent)._ring)
        return new_NCP((<NCPolynomialRing_plural>left._parent), _p)

    cpdef ModuleElement _sub_( left, ModuleElement right):
        """
        Subtract left and right.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: 3/2*x - 1/2*y - 1 # indirect doctest
            3/2*x - 1/2*y - 1

        """
        cdef ring *_ring = (<NCPolynomialRing_plural>left._parent)._ring

        cdef poly *_p
        singular_polynomial_sub(&_p, left._poly,
                                (<NCPolynomial_plural>right)._poly,
                                _ring)
        return new_NCP((<NCPolynomialRing_plural>left._parent), _p)

    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        Multiply ``self`` with a base ring element.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: 3/2*x # indirect doctest
            3/2*x
        """

        cdef ring *_ring = (<NCPolynomialRing_plural>self._parent)._ring
        if not left:
            return (<NCPolynomialRing_plural>self._parent)._zero_element
        cdef poly *_p
        singular_polynomial_rmul(&_p, self._poly, left, _ring)
        return new_NCP((<NCPolynomialRing_plural>self._parent),_p)

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Multiply ``self`` with a base ring element.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: x* (2/3) # indirect doctest
            2/3*x
        """
        return self._rmul_(right)

    cpdef RingElement  _mul_(left, RingElement right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: (3/2*x - 1/2*y - 1) * (3/2*x + 1/2*y + 1) # indirect doctest
            9/4*x^2 + 3/2*x*y - 3/4*z - 1/4*y^2 - y - 1

        TEST::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: (x^2^30) * x^2^30
            Traceback (most recent call last):
            ...
            OverflowError: Exponent overflow (...).
        """
        # all currently implemented rings are commutative
        cdef poly *_p
        singular_polynomial_mul(&_p, left._poly,
                                 (<NCPolynomial_plural>right)._poly,
                                 (<NCPolynomialRing_plural>left._parent)._ring)
        return new_NCP((<NCPolynomialRing_plural>left._parent),_p)

    cpdef RingElement _div_(left, RingElement right):
        """
        Divide left by right

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = (x + y)/3 # indirect doctest
            sage: f.parent()
            Noncommutative Multivariate Polynomial Ring in x, z, y over Rational Field, nc-relations: {y*x: -x*y + z}

        TESTS::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        cdef poly *p
        cdef bint is_field = left._parent._base.is_field()
        if p_IsConstant((<NCPolynomial_plural>right)._poly, (<NCPolynomialRing_plural>right._parent)._ring):
            if is_field:
                singular_polynomial_div_coeff(&p, left._poly, (<NCPolynomial_plural>right)._poly, (<NCPolynomialRing_plural>right._parent)._ring)
                return new_NCP(left._parent, p)
            else:
                return left.change_ring(left.base_ring().fraction_field())/right
        else:
            return (<NCPolynomialRing_plural>left._parent).fraction_field()(left,right)

    def __pow__(NCPolynomial_plural self, exp, ignored):
        """
        Return ``self**(exp)``.

        The exponent must be an integer.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = x^3 + y
            sage: f^2
            x^6 + x^2*z + y^2

        TESTS::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
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
            return (<NCPolynomialRing_plural>self._parent)._one_element

        cdef ring *_ring = (<NCPolynomialRing_plural>self._parent)._ring
        cdef poly *_p
        singular_polynomial_pow(&_p, self._poly, exp, _ring)
        return new_NCP((<NCPolynomialRing_plural>self._parent),_p)

    def __neg__(self):
        """
        Return ``-self``.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = x^3 + y
            sage: -f
            -x^3 - y
        """
        cdef ring *_ring = (<NCPolynomialRing_plural>self._parent)._ring

        cdef poly *p
        singular_polynomial_neg(&p, self._poly, _ring)
        return new_NCP((<NCPolynomialRing_plural>self._parent), p)

    def reduce(self, I):
        """

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H.<x,y,z> = A.g_algebra({y*x:x*y-z, z*x:x*z+2*x, z*y:y*z-2*y})
            sage: I = H.ideal([y^2, x^2, z^2-H.one()],coerce=False)

        The result of reduction is not the normal form, if one reduces
        by a list of polynomials::

            sage: (x*z).reduce(I.gens())
            x*z

        However, if the argument is an ideal, then a normal form (reduction
        with respect to a two-sided Groebner basis) is returned::

            sage: (x*z).reduce(I)
            -x

        The Groebner basis shows that the result is correct::

            sage: I.std() #random
            Left Ideal (z^2 - 1, y*z - y, x*z + x, y^2, 2*x*y - z - 1, x^2) of
            Noncommutative Multivariate Polynomial Ring in x, y, z over Rational
            Field, nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(I.std().gens(),key=str)
            [2*x*y - z - 1, x*z + x, x^2, y*z - y, y^2, z^2 - 1]


        """
        cdef ideal *_I
        cdef NCPolynomialRing_plural parent = <NCPolynomialRing_plural>self._parent
        cdef int i = 0
        cdef ring *r = parent._ring
        cdef poly *res

        if(r != currRing): rChangeCurrRing(r)

        if isinstance(I, NCPolynomialIdeal):
            try:
                strat = I._groebner_strategy()
                return strat.normal_form(self)
            except (TypeError, NotImplementedError) as msg:
                pass
            I = I.gens()

        _I = idInit(len(I),1)
        for f in I:
            if not (isinstance(f, NCPolynomial_plural) \
                   and <NCPolynomialRing_plural>(<NCPolynomial_plural>f)._parent is parent):
                try:
                    f = parent._coerce_c(f)
                except TypeError as msg:
                    id_Delete(&_I,r)
                    raise TypeError(msg)

            _I.m[i] = p_Copy((<NCPolynomial_plural>f)._poly, r)
            i+=1

        #the second parameter would be qring!
        res = kNF(_I, NULL, self._poly)
        id_Delete(&_I,r)
        return new_NCP(parent,res)

    def _repr_(self):
        """
        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = x^3 + y*x*z + z
            sage: f # indirect doctest
            x^3 - x*z*y + z^2 + z
        """
        cdef ring *_ring = (<NCPolynomialRing_plural>self._parent)._ring
        s = singular_polynomial_str(self._poly, _ring)
        return s

    cpdef _repr_short_(self):
        """
        This is a faster but less pretty way to print polynomials. If
        available it uses the short SINGULAR notation.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = x^3 + y
            sage: f._repr_short_()
            'x3+y'
        """
        cdef ring *_ring = (<NCPolynomialRing_plural>self._parent)._ring
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

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: latex(f) # indirect doctest
            - x^{2} y - z^{2} - \frac{25}{27} y^{3}
        """
        cdef ring *_ring = (<NCPolynomialRing_plural>self._parent)._ring
        gens = self.parent().latex_variable_names()
        base = self.parent().base()
        return singular_polynomial_latex(self._poly, _ring, base, gens)

    def _repr_with_changed_varnames(self, varnames):
        """
        Return string representing this polynomial but change the
        variable names to ``varnames``.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: print f._repr_with_changed_varnames(['FOO', 'BAR', 'FOOBAR'])
            -FOO^2*FOOBAR - BAR^2 - 25/27*FOOBAR^3
        """
        return  singular_polynomial_str_with_changed_varnames(self._poly, (<NCPolynomialRing_plural>self._parent)._ring, varnames)

    def degree(self, NCPolynomial_plural x=None):
        """
        Return the maximal degree of this polynomial in ``x``, where
        ``x`` must be one of the generators for the parent of this
        polynomial.

        INPUT:

        - ``x`` - multivariate polynomial (a generator of the parent of
          self) If x is not specified (or is ``None``), return the total
          degree, which is the maximum degree of any monomial.

        OUTPUT:

        integer

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
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

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: P(0).degree(x)
            -1
            sage: P(1).degree(x)
            0

        """
        cdef ring *r = (<NCPolynomialRing_plural>self._parent)._ring
        cdef poly *p = self._poly
        if not x:
            return singular_polynomial_deg(p,NULL,r)

        # TODO: we can do this faster
        if not x in self._parent.gens():
            raise TypeError("x must be one of the generators of the parent.")

        return singular_polynomial_deg(p, (<NCPolynomial_plural>x)._poly, r)

    def total_degree(self):
        """
        Return the total degree of ``self``, which is the maximum degree
        of all monomials in ``self``.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
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

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: R(0).total_degree()
            -1
            sage: R(1).total_degree()
            0
        """
        cdef poly *p = self._poly
        cdef ring *r = (<NCPolynomialRing_plural>self._parent)._ring
        return singular_polynomial_deg(p,NULL,r)

    def degrees(self):
        """
        Returns a tuple with the maximal degree of each variable in
        this polynomial.  The list of degrees is ordered by the order
        of the generators.

        EXAMPLES::

            sage: A.<y0,y1,y2> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y1*y0:-y0*y1 + y2},  order='lex')
            sage: R.inject_variables()
            Defining y0, y1, y2
            sage: q = 3*y0*y1*y1*y2; q
            3*y0*y1^2*y2
            sage: q.degrees()
            (1, 2, 1)
            sage: (q + y0^5).degrees()
            (5, 2, 1)
        """
        cdef poly *p = self._poly
        cdef ring *r = (<NCPolynomialRing_plural>self._parent)._ring
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
        :meth:`monomial_coefficient` which returns the coefficient in the
        base ring of a monomial.

        INPUT:

        - ``degrees`` - Can be any of:
                - a dictionary of degree restrictions
                - a list of degree restrictions (with None in the unrestricted variables)
                - a monomial (very fast, but not as flexible)

        OUTPUT:

        element of the parent of this element.

        .. NOTE::

           For coefficients of specific monomials, look at :meth:`monomial_coefficient`.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f=x*y+y+5
            sage: f.coefficient({x:0,y:1})
            1
            sage: f.coefficient({x:0})
            y + 5
            sage: f=(1+y+y^2)*(1+x+x^2)
            sage: f.coefficient({x:0})
            z + y^2 + y + 1

            sage: f.coefficient(x)
            y^2 - y + 1

            sage: f.coefficient([0,None]) # not tested
            y^2 + y + 1

        Be aware that this may not be what you think! The physical
        appearance of the variable x is deceiving -- particularly if
        the exponent would be a variable. ::

            sage: f.coefficient(x^0) # outputs the full polynomial
            x^2*y^2 + x^2*y + x^2 + x*y^2 - x*y + x + z + y^2 + y + 1

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f=x*y+5
            sage: c=f.coefficient({x:0,y:0}); c
            5
            sage: parent(c)
            Noncommutative Multivariate Polynomial Ring in x, z, y over Finite Field of size 389, nc-relations: {y*x: -x*y + z}

        AUTHOR:

        - Joel B. Mohler (2007-10-31)
        """
        cdef poly *_degrees = <poly*>0
        cdef poly *p = self._poly
        cdef ring *r = (<NCPolynomialRing_plural>self._parent)._ring
        cdef poly *newp = p_ISet(0,r)
        cdef poly *newptemp
        cdef int i
        cdef int flag
        cdef int gens = self._parent.ngens()
        cdef int *exps = <int*>sig_malloc(sizeof(int)*gens)
        for i from 0<=i<gens:
            exps[i] = -1

        if isinstance(degrees, NCPolynomial_plural) and self._parent is (<NCPolynomial_plural>degrees)._parent:
            _degrees = (<NCPolynomial_plural>degrees)._poly
            if pLength(_degrees) != 1:
                raise TypeError("degrees must be a monomial")
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
            raise TypeError("The input degrees must be a dictionary of variables to exponents.")

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

        return new_NCP(self.parent(),newp)

    def monomial_coefficient(self, NCPolynomial_plural mon):
        """
        Return the coefficient in the base ring of the monomial ``mon`` in
        ``self``, where ``mon`` must have the same parent as ``self``.

        This function contrasts with the function :meth:`coefficient`
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` - a monomial

        OUTPUT:

        coefficient in base ring

        .. SEEALSO::

            For coefficients in a base ring of fewer variables, look at :meth:`coefficient`

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y

            The parent of the return is a member of the base ring.
            sage: f = 2 * x * y
            sage: c = f.monomial_coefficient(x*y); c
            2
            sage: c.parent()
            Finite Field of size 389

            sage: f = y^2 + y^2*x - x^9 - 7*x + 5*x*y
            sage: f.monomial_coefficient(y^2)
            1
            sage: f.monomial_coefficient(x*y)
            5
            sage: f.monomial_coefficient(x^9)
            388
            sage: f.monomial_coefficient(x^10)
            0
        """
        cdef poly *p = self._poly
        cdef poly *m = mon._poly
        cdef ring *r = (<NCPolynomialRing_plural>self._parent)._ring

        if not mon._parent is self._parent:
            raise TypeError("mon must have same parent as self.")

        while(p):
            if p_ExpVectorEqual(p, m, r) == 1:
                return si2sa(p_GetCoeff(p, r), r, (<NCPolynomialRing_plural>self._parent)._base)
            p = pNext(p)

        return (<NCPolynomialRing_plural>self._parent)._base._zero_element

    def dict(self):
        """
        Return a dictionary representing ``self``. This dictionary is in
        the same format as the generic MPolynomial: The dictionary
        consists of ``ETuple:coefficient`` pairs.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y

            sage: f = (2*x*y^3*z^2 + (7)*x^2 + (3))
            sage: f.dict()
            {(0, 0, 0): 3, (1, 2, 3): 2, (2, 0, 0): 7}
        """
        cdef poly *p
        cdef ring *r
        cdef int n
        cdef int v
        r = (<NCPolynomialRing_plural>self._parent)._ring
        if r!=currRing: rChangeCurrRing(r)
        base = (<NCPolynomialRing_plural>self._parent)._base
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


    cdef long _hash_c(self):
        """
        See :meth:`__hash__`
        """
        cdef poly *p
        cdef ring *r
        cdef int n
        cdef int v
        r = (<NCPolynomialRing_plural>self._parent)._ring
        if r!=currRing: rChangeCurrRing(r)
        base = (<NCPolynomialRing_plural>self._parent)._base
        p = self._poly
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
        Same as :meth:`monomial_coefficient` but for exponent vectors.

        INPUT:

        - ``x`` - a tuple or, in case of a single-variable MPolynomial
          ring ``x`` can also be an integer.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = (-10*x^3*y + 17*x*y)* ( 15*z^3 + 2*x*y*z - 1); f
            20*x^4*z*y^2 - 150*x^3*z^3*y - 20*x^3*z^2*y + 10*x^3*y - 34*x^2*z*y^2 - 134*x*z^3*y + 34*x*z^2*y - 17*x*y
            sage: f[4,1,2]
            20
            sage: f[1,0,1]
            372
            sage: f[0,0,0]
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
        cdef ring *r = (<NCPolynomialRing_plural>self._parent)._ring
        cdef int i

        if isinstance(x, NCPolynomial_plural):
            return self.monomial_coefficient(x)
        if not isinstance(x, tuple):
            try:
                x = tuple(x)
            except TypeError:
                x = (x,)

        if len(x) != (<NCPolynomialRing_plural>self._parent).__ngens:
            raise TypeError("x must have length self.ngens()")

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
                return si2sa(p_GetCoeff(p, r), r, (<NCPolynomialRing_plural>self._parent)._base)
            p = pNext(p)

        p_Delete(&m,r)
        return (<NCPolynomialRing_plural>self._parent)._base._zero_element

    def exponents(self, as_ETuples=True):
        """
        Return the exponents of the monomials appearing in this polynomial.

        INPUT:

        - ``as_ETuples`` - (default: ``True``) if ``True`` returns the result as an list of ETuples
          otherwise returns a list of tuples


        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y
            sage: f = x^3 + y + 2*z^2
            sage: f.exponents()
            [(3, 0, 0), (0, 2, 0), (0, 0, 1)]
            sage: f.exponents(as_ETuples=False)
            [(3, 0, 0), (0, 2, 0), (0, 0, 1)]
        """
        cdef poly *p
        cdef ring *r
        cdef int v
        cdef list pl, ml

        r = (< NCPolynomialRing_plural>self._parent)._ring
        p = self._poly

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

    def is_homogeneous(self):
        """
        Return ``True`` if this polynomial is homogeneous.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: (x+y+z).is_homogeneous()
            True
            sage: (x.parent()(0)).is_homogeneous()
            True
            sage: (x+y^2+z^3).is_homogeneous()
            False
            sage: (x^2 + y^2).is_homogeneous()
            True
            sage: (x^2 + y^2*x).is_homogeneous()
            False
            sage: (x^2*y + y^2*x).is_homogeneous()
            True
        """
        cdef ring *_ring = (<NCPolynomialRing_plural>self._parent)._ring
        if(_ring != currRing): rChangeCurrRing(_ring)
        return bool(pIsHomogeneous(self._poly))


    def is_monomial(self):
        """
        Return ``True`` if this polynomial is a monomial.

        A monomial is defined to be a product of generators with
        coefficient 1.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
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
        _ring = (<NCPolynomialRing_plural>self._parent)._ring

        if self._poly == NULL:
            return True

        if(_ring != currRing): rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        ret = bool((not self._poly.next) and _ring.cf.nIsOne(_n))

        p_Delete(&_p, _ring)
        return ret

    def monomials(self):
        """
        Return the list of monomials in ``self``

        The returned list is decreasingly ordered by the term ordering
        of ``self.parent()``.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: f = x + (3*2)*y*z^2 + (2+3)
            sage: f.monomials()
            [x, z^2*y, 1]
            sage: f = P(3^2)
            sage: f.monomials()
            [1]

        TESTS::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
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
        cdef NCPolynomialRing_plural parent = <NCPolynomialRing_plural>self._parent
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
            l.append( new_NCP(parent,p) )
            p = t

        return l

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate
        polynomial.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x^2
            sage: f.constant_coefficient()
            0
        """
        cdef poly *p = self._poly
        cdef ring *r = (<NCPolynomialRing_plural>self._parent)._ring
        if p == NULL:
            return (<NCPolynomialRing_plural>self._parent)._base._zero_element

        while p.next:
            p = pNext(p)

        if p_LmIsConstant(p, r):
            return si2sa( p_GetCoeff(p, r), r, (<NCPolynomialRing_plural>self._parent)._base )
        else:
            return (<NCPolynomialRing_plural>self._parent)._base._zero_element

    cpdef is_constant(self):
        """
        Return ``True`` if this polynomial is constant.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(GF(389), 3)
            sage: P = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: P.inject_variables()
            Defining x, z, y
            sage: x.is_constant()
            False
            sage: P(1).is_constant()
            True
        """
        return bool(p_IsConstant(self._poly, (<NCPolynomialRing_plural>self._parent)._ring))

    def lm(NCPolynomial_plural self):
        """
        Returns the lead monomial of ``self`` with respect to the term
        order of ``self.parent()``.

        In Sage a monomial is a product of variables in some power
        without a coefficient.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(GF(7), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, y, z
            sage: f = x^1*y^2 + y^3*z^4
            sage: f.lm()
            x*y^2
            sage: f = x^3*y^2*z^4 + x^3*y^2*z^1
            sage: f.lm()
            x^3*y^2*z^4

            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='deglex')
            sage: R.inject_variables()
            Defining x, y, z
            sage: f = x^1*y^2*z^3 + x^3*y^2*z^0
            sage: f.lm()
            x*y^2*z^3
            sage: f = x^1*y^2*z^4 + x^1*y^1*z^5
            sage: f.lm()
            x*y^2*z^4

            sage: A.<x,y,z> = FreeAlgebra(GF(127), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='degrevlex')
            sage: R.inject_variables()
            Defining x, y, z
            sage: f = x^1*y^5*z^2 + x^4*y^1*z^3
            sage: f.lm()
            x*y^5*z^2
            sage: f = x^4*y^7*z^1 + x^4*y^2*z^3
            sage: f.lm()
            x^4*y^7*z

        """
        cdef poly *_p
        cdef ring *_ring
        _ring = (<NCPolynomialRing_plural>self._parent)._ring
        if self._poly == NULL:
            return (<NCPolynomialRing_plural>self._parent)._zero_element
        _p = p_Head(self._poly, _ring)
        p_SetCoeff(_p, n_Init(1,_ring), _ring)
        p_Setm(_p,_ring)
        return new_NCP((<NCPolynomialRing_plural>self._parent), _p)

    def lc(NCPolynomial_plural self):
        """
        Leading coefficient of this polynomial with respect to the
        term order of ``self.parent()``.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(GF(7), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, y, z

            sage: f = 3*x^1*y^2 + 2*y^3*z^4
            sage: f.lc()
            3

            sage: f = 5*x^3*y^2*z^4 + 4*x^3*y^2*z^1
            sage: f.lc()
            5
        """

        cdef poly *_p
        cdef ring *_ring
        cdef number *_n
        _ring = (<NCPolynomialRing_plural>self._parent)._ring

        if self._poly == NULL:
            return (<NCPolynomialRing_plural>self._parent)._base._zero_element

        if(_ring != currRing): rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        ret =  si2sa(_n, _ring, (<NCPolynomialRing_plural>self._parent)._base)
        p_Delete(&_p, _ring)
        return ret

    def lt(NCPolynomial_plural self):
        """
        Leading term of this polynomial.

        In Sage a term is a product of variables in some power and a
        coefficient.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(GF(7), 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, y, z

            sage: f = 3*x^1*y^2 + 2*y^3*z^4
            sage: f.lt()
            3*x*y^2

            sage: f = 5*x^3*y^2*z^4 + 4*x^3*y^2*z^1
            sage: f.lt()
            -2*x^3*y^2*z^4
        """
        if self._poly == NULL:
            return (<NCPolynomialRing_plural>self._parent)._zero_element

        return new_NCP((<NCPolynomialRing_plural>self._parent),
                                           p_Head(self._poly,(<NCPolynomialRing_plural>self._parent)._ring))

    def is_zero(self):
        """
        Return ``True`` if this polynomial is zero.

        EXAMPLES::

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y

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

            sage: A.<x,z,y> = FreeAlgebra(QQ, 3)
            sage: R = A.g_algebra(relations={y*x:-x*y + z},  order='lex')
            sage: R.inject_variables()
            Defining x, z, y

            sage: bool(x) # indirect doctest
            True
            sage: bool(x-x)
            False
        """
        if self._poly:
            return True
        else:
            return False


#####################################################################


cdef inline NCPolynomial_plural new_NCP(NCPolynomialRing_plural parent,
        poly *juice):
    """
    Construct NCPolynomial_plural from parent and SINGULAR poly.

    EXAMPLES::

        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: H = A.g_algebra({z*x:x*z+2*x, z*y:y*z-2*y})
        sage: H.gen(2)   # indirect doctest
        z

    """
    cdef NCPolynomial_plural p = NCPolynomial_plural.__new__(NCPolynomial_plural)
    p._parent = parent
    p._poly = juice
    p_Normalize(p._poly, parent._ring)
    return p




cpdef MPolynomialRing_libsingular new_CRing(RingWrap rw, base_ring):
    """
    Construct MPolynomialRing_libsingular from ringWrap, assumming the ground field to be base_ring

    EXAMPLES::

        sage: H.<x,y,z> = PolynomialRing(QQ, 3)
        sage: from sage.libs.singular.function import singular_function

        sage: ringlist = singular_function('ringlist')
        sage: ring = singular_function("ring")

        sage: L = ringlist(H, ring=H); L
        [0, ['x', 'y', 'z'], [['dp', (1, 1, 1)], ['C', (0,)]], [0]]

        sage: len(L)
        4

        sage: W = ring(L, ring=H); W
        <RingWrap>

        sage: from sage.rings.polynomial.plural import new_CRing
        sage: R = new_CRing(W, H.base_ring())
        sage: R # indirect doctest
        Multivariate Polynomial Ring in x, y, z over Rational Field

    Check that :trac:`13145` has been resolved::

        sage: h = hash(R.gen() + 1) # sets currRing
        sage: from sage.libs.singular.ring import ring_refcount_dict, currRing_wrapper
        sage: curcnt = ring_refcount_dict[currRing_wrapper()]
        sage: newR = new_CRing(W, H.base_ring())
        sage: ring_refcount_dict[currRing_wrapper()] - curcnt
        1
    """
    assert( rw.is_commutative() )

    cdef MPolynomialRing_libsingular self = <MPolynomialRing_libsingular>MPolynomialRing_libsingular.__new__(MPolynomialRing_libsingular)

    self._ring = rw._ring

    wrapped_ring = wrap_ring(self._ring)
    cur_refcnt = sage.libs.singular.ring.ring_refcount_dict.get(wrapped_ring, 0)
    sage.libs.singular.ring.ring_refcount_dict[wrapped_ring] = cur_refcnt + 1

    self._ring.ShortOut = 0

    self.__ngens = rw.ngens()
    self.__term_order =  TermOrder(rw.ordering_string(), force=True)

    ParentWithGens.__init__(self, base_ring, rw.var_names())
#    self._populate_coercion_lists_()  # ???

    #MPolynomialRing_generic.__init__(self, base_ring, n, names, order)
    self._has_singular = True
#    self._relations = self.relations()

    return self

cpdef NCPolynomialRing_plural new_NRing(RingWrap rw, base_ring):
    """
    Construct NCPolynomialRing_plural from ringWrap, assumming the ground field to be base_ring

    EXAMPLES::

        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: H = A.g_algebra({y*x:x*y-1})
        sage: H.inject_variables()
        Defining x, y, z
        sage: z*x
        x*z
        sage: z*y
        y*z
        sage: y*x
        x*y - 1
        sage: I = H.ideal([y^2, x^2, z^2-1])
        sage: I._groebner_basis_libsingular()
        [1]

        sage: from sage.libs.singular.function import singular_function

        sage: ringlist = singular_function('ringlist')
        sage: ring = singular_function("ring")

        sage: L = ringlist(H, ring=H); L
        [
                                                                   [0 1 1]
                                                                   [0 0 1]
        0, ['x', 'y', 'z'], [['dp', (1, 1, 1)], ['C', (0,)]], [0], [0 0 0],
        <BLANKLINE>
        [ 0 -1  0]
        [ 0  0  0]
        [ 0  0  0]
        ]
        sage: len(L)
        6

        sage: W = ring(L, ring=H); W
        <noncommutative RingWrap>

        sage: from sage.rings.polynomial.plural import new_NRing
        sage: R = new_NRing(W, H.base_ring())
        sage: R # indirect doctest
        Noncommutative Multivariate Polynomial Ring in x, y, z over
        Rational Field, nc-relations: {y*x: x*y - 1}
    """

    assert( not rw.is_commutative() )

    cdef NCPolynomialRing_plural self = <NCPolynomialRing_plural>NCPolynomialRing_plural.__new__(NCPolynomialRing_plural)
    self._ring = rw._ring

    wrapped_ring = wrap_ring(self._ring)
    cur_refcnt = sage.libs.singular.ring.ring_refcount_dict.get(wrapped_ring, 0)
    sage.libs.singular.ring.ring_refcount_dict[wrapped_ring] = cur_refcnt + 1

    self._ring.ShortOut = 0

    self.__ngens = rw.ngens()
    self.__term_order =  TermOrder(rw.ordering_string(), force=True)

    ParentWithGens.__init__(self, base_ring, rw.var_names())
#    self._populate_coercion_lists_()  # ???

    #MPolynomialRing_generic.__init__(self, base_ring, n, names, order)
    self._has_singular = True
    self._relations = self.relations()

    return self


def new_Ring(RingWrap rw, base_ring):
    """
    Constructs a Sage ring out of low level RingWrap, which wraps a pointer to a Singular ring.

    The constructed ring is either commutative or noncommutative depending on the Singular ring.

    EXAMPLES::

        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: H = A.g_algebra({y*x:x*y-1})
        sage: H.inject_variables()
        Defining x, y, z
        sage: z*x
        x*z
        sage: z*y
        y*z
        sage: y*x
        x*y - 1
        sage: I = H.ideal([y^2, x^2, z^2-1])
        sage: I._groebner_basis_libsingular()
        [1]

        sage: from sage.libs.singular.function import singular_function

        sage: ringlist = singular_function('ringlist')
        sage: ring = singular_function("ring")

        sage: L = ringlist(H, ring=H); L
        [
                                                                   [0 1 1]
                                                                   [0 0 1]
        0, ['x', 'y', 'z'], [['dp', (1, 1, 1)], ['C', (0,)]], [0], [0 0 0],
        <BLANKLINE>
        [ 0 -1  0]
        [ 0  0  0]
        [ 0  0  0]
        ]
        sage: len(L)
        6

        sage: W = ring(L, ring=H); W
        <noncommutative RingWrap>

        sage: from sage.rings.polynomial.plural import new_Ring
        sage: R = new_Ring(W, H.base_ring()); R
        Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: x*y - 1}

    """
    #    import warnings
    #    warnings.warn("This is a hack. Please, use it on your own risk...")
    if rw.is_commutative():
        return new_CRing(rw, base_ring)
    return new_NRing(rw, base_ring)

def SCA(base_ring, names, alt_vars, order='degrevlex'):
    """
    Return a free graded-commutative algebra

    This is also known as a free super-commutative algebra.

    INPUT:

    - ``base_ring`` -- the ground field
    - ``names`` -- a list of variable names
    - ``alt_vars`` -- a list of indices of to be anti-commutative variables (odd variables)
    - ``order`` -- ordering to be used for the constructed algebra

    EXAMPLES::

        sage: from sage.rings.polynomial.plural import SCA
        sage: E = SCA(QQ, ['x', 'y', 'z'], [0, 1], order = 'degrevlex')
        sage: E
        Quotient of Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {y*x: -x*y} by the ideal (y^2, x^2)
        sage: E.inject_variables()
        Defining xbar, ybar, zbar
        sage: x,y,z = (xbar,ybar,zbar)
        sage: y*x
        -x*y
        sage: z*x
        x*z
        sage: x^2
        0
        sage: y^2
        0
        sage: z^2
        z^2
        sage: E.one()
        1
    """
    n = len(names)
    alt_start = min(alt_vars)
    alt_end = max(alt_vars)
    assert( alt_start >= 0 )
    assert( (alt_end >= alt_start) and (alt_end < n) )

    relations = {} # {y*x:-x*y}
    from sage.algebras.free_algebra import FreeAlgebra
    A = FreeAlgebra(base_ring, n, names)
    for r in range(0, n-1, 1):
        for c in range(r+1, n, 1):
            if (r in alt_vars) and (c in alt_vars):
                relations[ A.gen(c) * A.gen(r) ] = - A.gen(r) * A.gen(c)

    cdef NCPolynomialRing_plural H = A.g_algebra(relations=relations, order=order)
    I = H.ideal([H.gen(i) * H.gen(i) for i in alt_vars]).twostd()
    return H.quotient(I)

def ExteriorAlgebra(base_ring, names,order='degrevlex'):
    """
    Return the exterior algebra on some generators

    This is also known as a Grassmann algebra. This is a finite
    dimensional algebra, where all generators anti-commute.

    See :wikipedia:`Exterior algebra`

    INPUT:

    - ``base_ring`` -- the ground ring
    - ``names`` -- a list of variable names

    EXAMPLES::

        sage: from sage.rings.polynomial.plural import ExteriorAlgebra
        sage: E = ExteriorAlgebra(QQ, ['x', 'y', 'z']) ; E #random
        Quotient of Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field, nc-relations: {z*x: -x*z, z*y: -y*z, y*x: -x*y} by the ideal (z^2, y^2, x^2)
        sage: sorted(E.cover().domain().relations().iteritems(),key=str)
        [(y*x, -x*y), (z*x, -x*z), (z*y, -y*z)]
        sage: sorted(E.cover().kernel().gens(),key=str)
        [x^2, y^2, z^2]
        sage: E.inject_variables()
        Defining xbar, ybar, zbar
        sage: x,y,z = (xbar,ybar,zbar)
        sage: y*x
        -x*y
        sage: all(v^2==0 for v in E.gens())
        True
        sage: E.one()
        1
    """
    n = len(names)
    relations = {} # {y*x:-x*y}
    from sage.algebras.free_algebra import FreeAlgebra
    A = FreeAlgebra(base_ring, n, names)
    for r in range(n-1):
        for c in range(r+1, n):
            relations[ A.gen(c) * A.gen(r) ] = - A.gen(r) * A.gen(c)

    cdef NCPolynomialRing_plural H = A.g_algebra(relations=relations,order=order)
    I = H.ideal([H.gen(i) * H.gen(i) for i in range(n)]).twostd()
    return H.quotient(I)

cdef poly *addwithcarry(poly *tempvector, poly *maxvector, int pos, ring *_ring):
    if p_GetExp(tempvector, pos, _ring) < p_GetExp(maxvector, pos, _ring):
      p_SetExp(tempvector, pos, p_GetExp(tempvector, pos, _ring)+1, _ring)
    else:
      p_SetExp(tempvector, pos, 0, _ring)
      tempvector = addwithcarry(tempvector, maxvector, pos + 1, _ring)
    p_Setm(tempvector, _ring)
    return tempvector
