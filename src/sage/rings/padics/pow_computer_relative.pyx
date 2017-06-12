from __future__ import absolute_import

from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.libs.gmp.mpz cimport mpz_init, mpz_clear, mpz_pow_ui
from sage.libs.flint.padic cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.nmod_vec cimport *
from sage.libs.flint.fmpz_vec cimport *
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_one, fmpz_mul, fmpz_set, fmpz_get_mpz, fmpz_clear, fmpz_pow_ui, fmpz_set_mpz, fmpz_fdiv_q_2exp

from cpython.object cimport Py_EQ, Py_NE
from sage.structure.richcmp cimport richcmp_not_equal
from sage.rings.integer cimport Integer
from sage.rings.all import ZZ
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint


cdef class PowComputer_relative(PowComputer_class):
    """
    A PowComputer for use in `p`-adics implemented by Sage polynomials.

    For a description of inputs see :func:`PowComputer_relative_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative
        sage: 1+1
        3
    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed = None):
        """
        Memory initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative
            sage: 1+1
            3
        """
        self.__allocated = 4

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: A = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod') # indirect doctest
            sage: TestSuite(A).run()
        """
        PowComputer_class.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)
        self.modulus = poly
        self.powhelper_oneunit = poly.parent()(0)
        self.powhelper_teichdiff = poly.parent()(0)
        self.base_ring = poly.base_ring()
        self.poly_ring = poly.parent()

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: A = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod')
            sage: del A
        """
        if self.__allocated >= 4:
            pass

    def __reduce__(self):
        """
        Pickling.

        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: A = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod') # indirect doctest
            sage: A._test_pickling() # indirect doctest

        """
        return PowComputer_relative_maker, (self.prime, self.cache_limit, self.prec_cap, self.ram_prec_cap, self.in_field, self.polynomial(), self._prec_type)

    def _repr_(self):
        """
        String representation of this powcomputer.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: A = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod'); A
            Relative PowComputer for 5
        """
        return "Relative PowComputer for %s" % self.prime

    cdef unsigned long capdiv(self, unsigned long n):
        """
        Returns ceil(n / e).
        """
        if self.e == 1: return n
        if n == 0: return 0
        return (n-1) / self.e + 1

    def polynomial(self, n=None, var='x'):
        """
        Returns ``None``.

        For consistency with subclasses.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative
            sage: A = PowComputer_relative(5, 20, 20, 60, False, None)
            sage: A.polynomial() is None
            True
        """
        return self.modulus

cdef class PowComputer_relative_unram(PowComputer_relative):
    """
    A PowComputer for a relative extension defined by an unramified polynomial.

    For a description of inputs see :func:`PowComputer_relative_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_unram
        sage: R.<a> = ZqFM(25)
        sage: S.<x> = R[]; f = x^2 - 5*x - a
        sage: A = PowComputer_relative_unram(5, 20, 20, 60, False, f); A
        Relative PowComputer for 5 with polynomial x^2 - 5*x - a

    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, _poly, shift_seed=None):
        """
        Memory initialization.

        TESTS::

            sage: 1+1
            3
        """
        self.__allocated = 8

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: 1+1
            3
        """
        if self.__allocated >= 8:
            pass

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_unram
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^2 - 5*x - a
            sage: A = PowComputer_relative_unram(5, 20, 20, 60, False, f)
            sage: TestSuite(A).run()

        """
        PowComputer_relative.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

        # Should these be relative or absolute?
        self.e = 1
        self.f = self.modulus.degree()

cdef class PowComputer_relative_eis(PowComputer_relative):
    """
    A PowComputer for a relative extension defined by an Eisenstein polynomial.

    For a description of inputs see :func:`PowComputer_relative_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_eis
        sage: R.<x> = ZZ[]; f = x^3 - 25*x + 5
        sage: A = PowComputer_flint_eis(5, 20, 20, 60, False, f); A
        FLINT PowComputer for 5 with polynomial x^3 - 25*x + 5
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_eis
            sage: R.<x> = ZZ[]; f = x^3 - 25*x + 5
            sage: A = PowComputer_flint_eis(5, 20, 20, 60, False, f)
            sage: type(A)
            <type 'sage.rings.padics.pow_computer_flint.PowComputer_flint_eis'>

        """
        PowComputer_relative.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

        # Should these be relative or absolute?
        self.e = self.modulus.degree()
        self.f = 1
        K = poly.base_ring().fraction_field()
        Qpmodulus = poly.change_ring(K)
        xep = (self.poly_ring.gen()**self.e - Qpmodulus) / K.uniformizer()
        _, _, self.pxe = Qpmodulus.xgcd(xep)
        self.pxe = self.pxe.change_ring(poly.base_ring())

def PowComputer_relative_maker(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, prec_type):
    """
    Return an appropriate relative PowComputer for the given input.

    INPUT:

    - ``prime`` -- an integral prime

    - ``cache_limit`` -- a non-negative integer, controlling the
      caching.  Powers of ``prime``, reductions of ``poly`` modulo
      different powers of ``prime`` and inverses of the leading
      coefficient modulo different powers of ``prime`` are cached.
      Additional data is cached for ramified extensions.

    - ``prec_cap`` -- the power of `p` modulo which elements of
      largest precision are defined.

    - ``ram_prec_cap`` -- Approximately ``e*prec_cap``, where ``e`` is
      the ramification degree of the extension.  For a ramified
      extension this is what Sage calls the precision cap of the ring.
      In fact, it's possible to have rings with precision cap not a
      multiple of `e`, in which case the actual relationship between
      ``ram_prec_cap`` and ``prec_cap`` is that
      ``prec_cap = ceil(n/e)``

    - ``in_field`` -- (boolean) whether the associated ring is
      actually a field.

    - ``poly`` -- the polynomial defining the extension.

    - `prec_type`` -- one of ``"capped-rel"``, ``"capped-abs"`` or
      ``"fixed-mod"``, the precision type of the ring.

    .. NOTE::

        Because of the way templates work, this function imports the
        class of its return value from the appropriate element files.
        This means that the returned PowComputer will have the
        appropriate compile-time-type for Cython.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
        sage: R.<x> = ZZ[]
        sage: A = PowComputer_flint_maker(3, 20, 20, 20, False, x^3 + 2*x + 1, 'capped-rel'); type(A)
        <type 'sage.rings.padics.qadic_flint_CR.PowComputer_'>
        sage: A = PowComputer_flint_maker(3, 20, 20, 20, False, x^3 + 2*x + 1, 'capped-abs'); type(A)
        <type 'sage.rings.padics.qadic_flint_CA.PowComputer_'>
        sage: A = PowComputer_flint_maker(3, 20, 20, 20, False, x^3 + 2*x + 1, 'fixed-mod'); type(A)
        <type 'sage.rings.padics.qadic_flint_FM.PowComputer_'>

    """
    PC = PowComputer_relative_eis(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)
    PC._prec_type = 'fixed-mod'
    return PC
    #if prec_type == 'capped-rel':
    #    pass
    #elif prec_type == 'capped-abs':
    #    pass
    #elif prec_type == 'fixed-mod':
    #    pass
    #elif prec_type == 'floating-point':
    #    pass
    #else:
    #    raise ValueError("unknown prec_type `%s`" % prec_type)
    #return PowComputer_(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)
