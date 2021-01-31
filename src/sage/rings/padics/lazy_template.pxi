r"""
Lazy template for `p`-adic rings and fields

In order to use this template you need to write a linkage file and
gluing file.

The linkage file implements a common API that is then used in the
class FMElement defined here.  
See sage/libs/linkages/padics/lazy/API.pxi for the functions needed.

The gluing file does the following:
- ctypedef's ``cdigit``, ``cdigit_ptr``, ``celement`` and ``randgen``
  to be the appropriate types
- includes the linkage file
- includes this template
- defines classes inheriting for classes defined in this file,
  and implements any desired extra methods
See padic_lazy_element.pxd and padic_lazy_element.pyx for an example.

AUTHOR:

- Xavier Caruso (2021-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import re
from sage.misc import persist

from libc.stdlib cimport malloc, free
from sage.libs.gmp.mpz cimport mpz_sizeinbase, mpz_tstbit

from sage.structure.element import coerce_binop
from sage.structure.element cimport have_same_parent
from sage.structure.coerce cimport coercion_model
from sage.misc.prandom import randint

from sage.rings.all import ZZ
from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.padics.padic_lazy_errors cimport *
from sage.rings.padics.padic_lazy_errors import raise_error

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1


cdef cdigit tmp_digit
cdef celement tmp_poly
digit_init(tmp_digit)
element_init(tmp_poly)


cdef class LazyElement(pAdicGenericElement):
    r"""
    Template class for lazy `p`-adic elements.

    EXAMPLES:

        sage: from sage.rings.padics.padic_lazy_element import LazyElement
        sage: R = ZpL(5)
        sage: a = R(1)
        sage: isinstance(a, LazyElement)
        True
    """
    def __init__(self, parent):
        r"""
        Initialize this element.

        TESTS::

            sage: R = ZpL(5)
            sage: a = R(1/2)  # indirect doctest
            sage: TestSuite(a).run()
        """
        pAdicGenericElement.__init__(self, parent)
        self.prime_pow = <PowComputer_class>self._parent.prime_pow
        self._valuation = 0
        self._precrel = 0
        self._precbound = maxordp

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle
        this element.

        TESTS::

            sage: R = ZpL(5)
            sage: a = R(0)
            sage: loads(dumps(a)) == a
            True
        """
        raise NotImplementedError("must be implemented in subclasses")

    cpdef bint _is_base_elt(self, p) except -1:
        r"""
        Return ``True`` if this element is an element of Zp or Qp (rather than
        an extension).

        INPUT:

        - ``p`` -- a prime, which is compared with the parent of this element.

        EXAMPLES::

            sage: a = ZpL(5)(3)
            sage: a._is_base_elt(5)
            True
            sage: a._is_base_elt(17)
            False
        """
        return p == self._parent.prime()

    cdef cdigit_ptr _getdigit_relative(self, long i):
        r"""
        Return a pointer on the `i`-th significant digit of this number.

        .. NOTE:

            This function do not check that the requested digit
            has been already computed.

        INPUT:

        - ``i`` -- a positive integer

        """
        pass

    cdef cdigit_ptr _getdigit_absolute(self, long i):
        r"""
        Return a pointer on the digit in position `i` of
        this number.

        .. NOTE:

            This function do not check that the requested digit
            has been already computed.

        INPUT:

        - ``i`` -- an integer

        """
        pass

    cdef void _getslice_relative(self, celement slice, long start, long length):
        r"""
        Select a slice of the digits of this number.

        INPUT:

        - ``slice`` -- a ``celement`` to store the slice

        - ``start`` -- a positive integer, the starting position of the slice
          in relative precision

        - ``length`` -- a positive integer, the length of the slice


        .. NOTE::

            This methods only sets up a pointer to the requested slice
            (the slice is not copied). Hence any future modification
            of the slice ``slice`` will affect this number.
        """
        pass

    cdef int _next_c(self):
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code which is a superposition of the following:

        - ``0`` -- no error
        - ``ERROR_ABANDON    = 1`` -- computation has been abandonned
        - ``ERROR_NOTDEFINED = 2`` -- a number is not defined
        - ``ERROR_PRECISION  = 4`` -- not enough precision
        - ``ERROR_OVERFLOW   = 8`` -- overflow
        - ``ERROR_NOTSQUARE  = 16`` -- not a square
        - ``ERROR_INTEGRAL   = 32`` -- not in the integer ring
        - ``ERROR_DIVISION   = 64`` -- division by zero (or something indistinguishable from zero)
        - ``ERROR_CIRCULAR   = 128`` -- circular definition
        """
        raise NotImplementedError("must be implemented in subclasses")

    cdef int _jump_c(self, long prec):
        r"""
        Compute the digits of this number until the absolute precision
        ``prec``.

        INPUT:

        - ``prec`` -- an integer

        OUTPUT:

        An error code (see :meth:`_next_c` for details).

        """
        cdef int error
        cdef long pr = min(prec, self._precbound)
        while self._precrel + self._valuation < pr:
            error = self._next_c()
            if error:
                return error
        if prec > self._precbound:
            return ERROR_PRECISION
        return 0

    cdef int _jump_relative_c(self, long prec, long halt):
        r"""
        Compute the digits of this number until the relative precision
        ``prec``.

        INPUT:

        - ``prec`` -- an integer

        - ``halt`` -- an integer; the absolute precision after which the 
          computation is abandonned if the first significant digit has not 
          been found yet

        OUTPUT:

        An error code (see :meth:`_next_c` for details).

        """
        if self._valuation >= maxordp:
            return 0
        cdef int error = 0
        if self._valuation <= -maxordp:
            error = self._next_c()
        cdef long valhalt = min(halt, self._precbound)
        while not error and self._valuation < valhalt and self._precrel == 0:
            error = self._next_c()
        if self._valuation >= self._precbound:
            error |= ERROR_PRECISION
        elif self._valuation >= halt:
            error |= ERROR_ABANDON
        if not error:
            error = self._jump_c(self._valuation + prec)
        return error

    cdef int _init_jump(self) except -1:
        r"""
        Start the computation of the digits of this number.
        This method should be called at initialization.

        OUTPUT:

        An error code (see :meth:`_next_c` for details).

        """
        cdef int error = 0
        if self._precbound < maxordp:
            error = self._jump_c(self._precbound)
            raise_error(error, True)
        return error

    def digit(self, i):
        r"""
        Return the coefficient of `p^i` in the `p`-adic expansion
        of this number.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: R = ZpL(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.digit(0)
            0
            sage: a.digit(1)
            4
            sage: a.digit(7)
            3
            sage: a.digit(100)
            1

        As a shortcut, one can use the bracket operator::

            sage: a[100]
            ...
            1

        If the digit is not known, an error is raised::

            sage: b = a.add_bigoh(10)
            sage: b.digit(20)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision

        TESTS::

            sage: a.digit(-1)
            0
        """
        cdef int error = self._jump_c(i+1)
        if error:
            raise_error(error)
        cdef cdigit_ptr coeff = self._getdigit_absolute(i)
        return digit_get_sage(coeff)

    def expansion(self, n=None, lift_mode='simple'):
        r"""
        Return an iterator over the list of coefficients in a `p`-adic 
        expansion of this element, that is the list of `a_i` so that 
        this element can be expressed as

        .. MATH::

            \pi^v \cdot \sum_{i=0}^\infty a_i \pi^i

        where `v` is the valuation of this element when the parent is
        a field, and `v = 0` otherwise and the `a_i` are between `0`
        and `p - 1`.

        INPUT:

        - ``n`` -- an integer or ``None`` (default ``None``); if
          given, returns the corresponding entries in the expansion.

        - ``lift_mode`` -- only ``'simple'`` is accepted

        OUTPUT:

        - If ``n`` is ``None``, an iterable giving a `p`-adic expansion
          of this element.

        - If ``n`` is an integer, the coefficient of `p^n` in the
          `p`-adic expansion of this element.

        EXAMPLES::

            sage: R = ZpL(7, print_mode="digits")
            sage: a = R(1/2021); a
            ...23615224635636163463

        Without any argument, this method returns an iterator over the
        digits of this number::

            sage: E = a.expansion()
            sage: E
            7-adic expansion of ...23615224635636163463
            sage: next(E)
            3
            sage: next(E)
            6
            sage: next(E)
            4

        On the contrary, passing in an integer returns the digit at the
        given position::

            sage: a.expansion(5)
            1
        """
        if lift_mode != 'simple':
            raise NotImplementedError
        if n is None:
            if self.prime_pow.in_field:
                start = self.valuation()
            else:
                start = 0
            return ExpansionIter(self, start, self._precbound)
        elif isinstance(n, slice):
            if n.step is not None:
                raise NotImplementedError("step is not allowed")
            return self.slice(n.start, n.stop, True)
        else:
            n = Integer(n)
            if n >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            return self.digit(n)

    def slice(self, start=None, stop=None, bound=False):
        r"""
        Return a slice of this number.

        INPUT:

        - ``start`` -- an integer or ``None`` (default: ``None``),
          the first position of the slice

        - ``stop`` -- an integer or ``None`` (default: ``None``),
          the first position not included in the slice

        - ``bound`` -- a boolean (default: ``False``); whether the
          precision on the output should be bounded or unbounded

        EXAMPLES::

            sage: K = QpL(7, prec=10)
            sage: a = K(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...

            sage: s = a.slice(3, 6)
            sage: s
            3*7^3 + 6*7^4 + 7^5 + ...

        In the above example, the precision on `b` remains unbounded::

            sage: s.precision_absolute()
            +Infinity

        Passing in ``bound=True`` changes this behaviour::

            sage: a.slice(3, 6, bound=True)
            3*7^3 + 6*7^4 + 7^5 + O(7^6)

        When ``start`` is omitted, the slice starts at the beginning of the number:

            sage: a.slice(stop=6)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + ...

            sage: b = K(20/21); b
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + 2*7^5 + 2*7^6 + 2*7^7 + 2*7^8 + ...
            sage: b.slice(stop=6)
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + 2*7^5 + ...

        As a shortcut, one can use the bracket operator.
        However, in this case, the precision is bounded::

            sage: a[3:6]
            3*7^3 + 6*7^4 + 7^5 + O(7^6)
            sage: b[:6]
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + 2*7^5 + O(7^6)

        TESTS::

            sage: x = a[:3]
            sage: x.slice(stop=3)
            3 + 6*7 + 4*7^2 + ...
            sage: x.slice(stop=4)
            3 + 6*7 + 4*7^2 + O(7^3)
        """
        if start is None:
            start = -maxordp
        elif start >= maxordp:
            raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
        if stop is None:
            stop = maxordp
        elif stop >= maxordp:
            raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
        cdef LazyElement x = element_class_slice(self._parent, self, start, stop, 0)
        if bound and x._precbound > stop:
            x = element_class_bound(self._parent, x, stop)
        return x

    def _repr_(self):
        r"""
        Return a string representation of this element.

        EXAMPLES:

        For unbounded elements, the number of printed terms is given by
        the default precision of the ring::

            sage: R = ZpL(5, 10)
            sage: a = R(1/2021)
            sage: a   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + 5^5 + 4*5^6 + 4*5^7 + 2*5^8 + 3*5^9 + ...

            sage: S = ZpL(5, 5)
            sage: b = S(1/2021)
            sage: b   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + ...

        On the contrary, bounded elements are printed until their own
        precision, regardless the default precision of the parent::

            sage: a[:15]   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + 5^5 + 4*5^6 + 4*5^7 + 2*5^8 + 3*5^9 + 2*5^10 + 5^11 + 5^12 + 5^13 + O(5^15)

            sage: b[:15]   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + 5^5 + 4*5^6 + 4*5^7 + 2*5^8 + 3*5^9 + 2*5^10 + 5^11 + 5^12 + 5^13 + O(5^15)

        """
        # This code should be integrated to the p-adic printer
        if self._valuation <= -maxordp:
            return "valuation not known"
        if self._is_exact_zero():
            return "0"
        if self._precbound >= maxordp:
            unbounded = True
            try:
                if self.prime_pow.in_field:
                    x = self.at_precision_relative(permissive=False)
                else:
                    x = self.at_precision_absolute(permissive=False)
            except (ValueError, PrecisionError, RecursionError):
                unbounded = False
                x = element_class_bound(self._parent, self, self._valuation + self._precrel)
        else:
            unbounded = False
            x = self
        s = pAdicGenericElement._repr_(x)
        mode = self._parent._printer.dict()['mode']
        if unbounded:
            s = re.sub(r'O\(.*\)$', '...', s)
        elif mode == "digits":
            s = re.sub(r'^\.\.\.\??', '...?', s)
        elif mode == "bars":
            s = re.sub(r'^\.\.\.\|?', '...?|', s)
        if s == "" or s == "...":
            s = "0 + ..."
        return s

    cdef bint _is_equal(self, LazyElement right, long prec, bint permissive) except -1:
        r"""
        C function for checking equality at a given precision.

        INPUT:

        - ``right`` -- the second element involved in the comparison

        - ``prec`` -- an integer, the precision at which the equality is checked

        - ``permissive`` -- a boolean; if ``True``, be silent if the precision
          on one input are less than ``prec``; otherwise, raise an error

        """
        cdef int error
        cdef long i
        if self._valuation <= -maxordp:
            error = self._next_c()
            raise_error(error)
        if right._valuation <= -maxordp:
            error = right._next_c()
            raise_error(error)
        for i in range(min(self._valuation, right._valuation), prec):
            error = self._jump_c(i+1)
            raise_error(error, permissive)
            if error:   # not enough precision
                return True
            error = right._jump_c(i+1)
            raise_error(error, permissive)
            if error:   # not enough precision
                return True
            if not digit_equal(self._getdigit_absolute(i), right._getdigit_absolute(i)):
                return False
        return True

    @coerce_binop
    def is_equal_at_precision(self, right, prec):
        r"""
        Compare this element with ``right`` at precision ``prec``.

        INPUT:

        - ``right`` -- a lazy `p`-adic number

        - ``prec`` -- an integer

        EXAMPLES::

            sage: R = ZpL(7, prec=10)
            sage: a = R(1/2); a
            4 + 3*7 + 3*7^2 + 3*7^3 + 3*7^4 + 3*7^5 + 3*7^6 + 3*7^7 + 3*7^8 + 3*7^9 + ...
            sage: b = R(99/2); b
            4 + 3*7 + 4*7^2 + 3*7^3 + 3*7^4 + 3*7^5 + 3*7^6 + 3*7^7 + 3*7^8 + 3*7^9 + ...

            sage: a.is_equal_at_precision(b, 1)
            True
            sage: a.is_equal_at_precision(b, 2)
            True
            sage: a.is_equal_at_precision(b, 3)
            False
        """
        return self._is_equal(right, min(prec, maxordp), False)

    @coerce_binop
    def is_equal_to(self, LazyElement right, prec=None):
        r"""
        Compare this element with ``right``.

        INPUT:

        - ``right`` -- a lazy `p`-adic number

        - ``prec`` -- an integer or ``None`` (default: ``None``); if 
          given, compare the two elements at this precision; otherwise
          use the default precision of the parent and raise an error
          if equality could not be decided

        EXAMPLES::

            sage: R = ZpL(7)

            sage: a = R(1/2)
            sage: b = R(1/3)
            sage: c = R(1/6)

            sage: a.is_equal_to(b)
            False

        When equality indeed holds, it is not possible to conclude by
        comparing more and more accurate approximations.
        In this case, the computation is abandonned::

            sage: a.is_equal_to(b + c)
            Traceback (most recent call last):
            ...
            PrecisionError: unable to decide equality; try to bound precision

        Bounding precision gives a 

            sage: a.is_equal_to(b + c, prec=20)
            True
            sage: a.is_equal_to(b + c, prec=100)
            True

        In some cases, equality tests may fail when the first different
        digit is beyond the default precision of the parent::

            sage: a.is_equal_to(a + 7^50)
            Traceback (most recent call last):
            ...
            PrecisionError: unable to decide equality; try to bound precision

        """
        cdef long default_prec
        if self is right:
            return True
        if self._valuation >= maxordp and right._valuation >= maxordp:
            return True
        if prec is None:
            default_prec = self._parent.default_prec()
            if self.prime_pow.in_field:
                self._jump_relative_c(default_prec, self._valuation + default_prec)
                right._jump_relative_c(default_prec, right._valuation + default_prec)
            else:
                self._jump_c(default_prec)
                right._jump_c(default_prec)
            prec = min(self._valuation + self._precrel, right._valuation + right._precrel)
            if not self._is_equal(right, prec, True):
                return False
            prec = min(self._precbound, right._precbound)
        else:
            prec = Integer(prec)
        if prec < maxordp:
            return self._is_equal(right, prec, True)
        else:
            raise PrecisionError("unable to decide equality; try to bound precision")

    def __eq__(self, other):
        if not have_same_parent(self, other):
            try:
                a, b = coercion_model.canonical_coercion(self, other)
            except TypeError:
                return False
            return a == b
        return self.is_equal_to(other)

    cpdef bint _is_exact_zero(self) except -1:
        r"""
        Return ``True`` if this element is an exact zero.

        EXAMPLES::

            sage: R = ZpL(5)
            sage: a = R(0); a
            0
            sage: a._is_exact_zero()
            True

            sage: b = a.add_bigoh(20)
            sage: b
            O(5^20)
            sage: b._is_exact_zero()
            False

        """
        return self._valuation >= maxordp

    cpdef bint _is_inexact_zero(self) except -1:
        r"""
        Return ``True`` if, at the current stage of computations, this
        number cannot be distinguished from zero.

        EXAMPLES::

            sage: R = ZpL(5, print_mode="digits")
            sage: a = R(20/21)

        Computations have not started yet, 

            sage: a._is_inexact_zero()
            True

        When we print `a`, the first digits are computed, after what
        we know that `a` is not zero::

            sage: a
            ...34010434010434010440
            sage: a._is_inexact_zero()
            False
        """
        return self._precrel == 0

    def is_exact(self):
        r"""
        Return ``True`` if this element is exact, that is if its
        precision is unbounded.

        EXAMPLES::

            sage: R = ZpL(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.is_exact()
            True

            sage: b = a.add_bigoh(10)
            sage: b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + O(5^10)
            sage: b.is_exact()
            False
        """
        return self._precbound >= maxordp

    def precision_absolute(self):
        r"""
        Return the absolute precision of this element.

        This is the power of `p` modulo which this element is known.
        For unbounded elements, this methods return `+\infty`.

        EXAMPLES::

            sage: R = ZpL(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.precision_absolute()
            +Infinity

            sage: b = a.add_bigoh(10)
            sage: b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + O(5^10)
            sage: b.precision_absolute()
            10

        TESTS::

            sage: s = R.selfref()
            sage: (1/s).precision_absolute()
            Traceback (most recent call last):
            ...
            PrecisionError: no lower bound on the valuation is known
        """
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if self._precbound >= maxordp:
            return Infinity
        return Integer(self._precrel + self._valuation)

    def precision_relative(self):
        """
        Returns the relative precision of this element.

        This is the power of `p` modulo which the unit part of this
        element is known.
        For unbounded nonzero elements, this methods return `+\infty`.

        EXAMPLES::

            sage: R = ZpL(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.precision_relative()
            +Infinity

            sage: b = a.add_bigoh(10)
            sage: b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + O(5^10)
            sage: b.precision_relative()
            9

        The relative precision of (exact and inexact) `0` is `0`::

            sage: x = R(0); x
            0
            sage: x.precision_relative()
            0

            sage: y = R(0, 10); y
            O(5^10)
            sage: y.precision_relative()
            0

        TESTS::

            sage: s = R.selfref()
            sage: (1/s).precision_relative()
            Traceback (most recent call last):
            ...
            PrecisionError: no lower bound on the valuation is known
        """
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if self._precbound >= maxordp and self._valuation < maxordp:
            return Infinity
        return Integer(self._precrel)

    def at_precision_absolute(self, prec=None, permissive=None):
        r"""
        Return this element bounded at the given precision.

        INPUT:

        - ``prec`` -- an integer or ``None`` (default: ``None``);
          if ``None``, use the default precision of the parent

        - ``permissive`` -- a boolean (default: ``False`` if ``prec``
          is given, ``True`` otherwise); if ``False``, raise an error
          if the precision of this element is not sufficient

        EXAMPLES::

            sage: R = ZpL(7, prec=10)
            sage: a = R(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...
            sage: a.at_precision_absolute(5)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)
            sage: a.at_precision_absolute(10)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + O(7^10)

            sage: b = a.add_bigoh(5)
            sage: b.at_precision_absolute(10)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision

            sage: b.at_precision_absolute(10, permissive=True)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)

        Bounding at a negative precision is not permitted over `\ZZ_p`::

            sage: a.at_precision_absolute(-1)
            Traceback (most recent call last):
            ...
            ValueError: precision must be nonnegative

        but, of course, it is over `\QQ_p`::

            sage: K = R.fraction_field()
            sage: K(a).at_precision_absolute(-1)
            O(7^-1)

        .. SEEALSO::

            :meth:`at_precision_relative`, :meth:`add_bigoh`
        """
        if prec is Infinity:
            if permissive is False and self._precbound < maxordp:
                raise_error(ERROR_PRECISION)
            return self
        if prec is None:
            if permissive is None:
                permissive = True
            prec = self._parent.default_prec()
        else:
            prec = Integer(prec)
            if not self.prime_pow.in_field and prec < 0:
                raise ValueError("precision must be nonnegative")
            if prec > maxordp:
                raise OverflowError("beyond maximal precision (which is %s)" % maxordp)
        error = self._jump_c(prec)
        if permissive is None:
            permissive = False
        raise_error(error, permissive)
        return element_class_bound((<LazyElement>self)._parent, self, prec)

    def add_bigoh(self, absprec):
        r"""
        Return a new element with absolute precision decreased to ``absprec``.

        INPUT:

        - ``absprec`` -- an integer or infinity

        EXAMPLES::

            sage: R = ZpL(7, prec=10)
            sage: a = R(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...
            sage: a.add_bigoh(5)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)

        When ``absprec`` is negative, we return an element in the fraction
        field::

            sage: b = a.add_bigoh(-1)
            sage: b
            O(7^-1)
            sage: b.parent()
            7-adic Field with lazy precision

        .. SEEALSO::

            :meth:`at_precision_absolute`, :meth:`at_precision_relative`
        """
        if absprec < 0 and (not self.prime_pow.in_field):
            self = element_class_bound(self._parent.fraction_field(), self)
        return self.at_precision_absolute(absprec, True)

    def at_precision_relative(self, prec=None, halt=True, permissive=None):
        r"""
        Return this element bounded at the given precision.

        INPUT:

        - ``prec`` -- an integer or ``None`` (default: ``None``);
          if ``None``, use the default precision of the parent

        - ``halt`` -- an integer or a boolean (default: ``True``);
          the absolute precision after which the computation is abandonned
          if the first significant digit has not been found yet;
          if ``True``, the default precision of the parent is used;
          if ``False``, the computation is never abandonned

        - ``permissive`` -- a boolean (default: ``False`` if ``prec``
          is given, ``True`` otherwise); if ``False``, raise an error
          if the precision of this element is not sufficient

        EXAMPLES::

            sage: R = ZpL(5, prec=10)
            sage: a = R(20/21); a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.at_precision_relative(5)
            4*5 + 4*5^2 + 5^4 + O(5^6)

        We illustrate the behaviour of the parameter ``halt``.
        We create a very small number whose first significant is bar beyond
        the default precision::

            sage: b = R(5^20)
            sage: b
            0 + ...

        Without any help, Sage does not run the computation far enough to determine
        the valuation and an error is raised::

            sage: b.at_precision_relative(5)
            Traceback (most recent call last):
            ...
            PrecisionError: computation has been abandonned; try to increase precision

        By setting the argument ``halt``, one can force the computation to continue
        until a prescribed limit::

            sage: b.at_precision_relative(5, halt=20)   # not enough to find the valuation
            Traceback (most recent call last):
            ...
            PrecisionError: computation has been abandonned; try to increase precision

            sage: b.at_precision_relative(5, halt=21)   # now, we're okay
            5^20 + O(5^25)

        .. NOTE:

            It is also possible to pass in ``halt=False`` but it is not recommended
            because the computation can hang forever if this element were `0`.

        .. SEEALSO::

            :meth:`at_precision_absolute`, :meth:`add_bigoh`

        TESTS::

            sage: a.at_precision_relative(-1)
            Traceback (most recent call last):
            ...
            ValueError: precision must be nonnegative
        """
        default_prec = self._parent.default_prec()
        if prec is None:
            if permissive is None:
                permissive = True
            prec = default_prec
        else:
            prec = Integer(prec)
            if prec < 0:
                raise ValueError("precision must be nonnegative")
            if prec > maxordp:
                raise OverflowError("beyond maximal precision (which is %s)" % maxordp)
        if halt is True:
            halt = self._valuation + default_prec
        elif halt is False:
            halt = maxordp
        else:
            halt = min(maxordp, halt)
        error = self._jump_relative_c(prec, halt)
        if permissive is None:
            permissive = False
        raise_error(error, permissive)
        return element_class_bound((<LazyElement>self)._parent, self, self._valuation +  prec)

    def lift_to_precision(self, absprec=None):
        """
        Return another element of the same parent, lifting this element
        and having absolute precision at least ``absprec``.

        INPUT:

        - ``absprec`` -- an integer or ``None`` (default: ``None``), the
          absolute precision of the result. If ``None``, the default
          precision of the parent is used.

        EXAMPLES::

            sage: R = ZpL(5, prec=10)
            sage: a = R(20/21, 5); a
            4*5 + 4*5^2 + 5^4 + O(5^5)

            sage: a.lift_to_precision(20)
            4*5 + 4*5^2 + 5^4 + O(5^20)

        When the precision is omitted, the default precision of the parent
        is used::

            sage: a.lift_to_precision()
            4*5 + 4*5^2 + 5^4 + O(5^10)

        When the parent is a field, the behaviour is slightly different since
        the default precision of the parent becomes the relative precision
        of the lifted element::

            sage: K = R.fraction_field()
            sage: K(a).lift_to_precision()
            4*5 + 4*5^2 + 5^4 + O(5^11)

        Note that the precision never decreases::

            sage: a.lift_to_precision(2)
            4*5 + 4*5^2 + 5^4 + O(5^5)

        In particular, unbounded element are not affected by this method::

            sage: b = R(20/21); b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: b.lift_to_precision()
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: b.lift_to_precision(2)
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
        """
        if self._precbound >= maxordp:
            return self
        cdef long prec
        cdef long default_prec = self._parent.default_prec()
        if absprec is None:
            if self.prime_pow.in_field:
                self._jump_relative_c(1, self._precbound)
                if self._precrel == 0:
                    return self._parent.zero()
                prec = self._valuation + default_prec
            else:
                prec = default_prec
        else:
            if absprec > maxordp:
                raise OverflowError("beyond the maximal precision (which is %s)" % maxordp)
            prec = absprec
        if prec <= self._precbound:
            return self
        cdef LazyElement ans = element_class_slice(self._parent, self, -maxordp, self._precbound, 0)
        ans._precbound = prec
        ans._init_jump()
        return ans

    cdef long valuation_c(self):
        r"""
        Return the best lower bound we have on the valuation of
        this element at the current stage of the computation.
        """
        return self._valuation

    def valuation(self, halt=True):
        r"""
        Return the valuation of this element.

        INPUT:

        - ``halt`` -- an integer or a boolean (default: ``True``);
          the absolute precision after which the computation is abandonned
          if the first significant digit has not been found yet;
          if ``True``, the default precision of the parent is used;
          if ``False``, the computation is never abandonned

        EXAMPLES::

            sage: R = ZpL(5, 10)
            sage: a = R(2001); a
            1 + 5^3 + 3*5^4 + ...
            sage: a.valuation()
            0

            sage: b = a - 1/a; b
            2*5^3 + 5^4 + 5^5 + 4*5^6 + 3*5^7 + 4*5^8 + 3*5^9 + 3*5^10 + 3*5^11 + 5^12 + ...
            sage: b.valuation()
            3

        The valuation of an exact zero is `+\infty`::

            sage: R(0).valuation()
            +Infinity

        The valuation of an inexact zero is its absolute precision::

            sage: R(0, 20).valuation()
            20

        We illustrate the behaviour of the parameter ``halt``.
        We create a very small number whose first significant is bar beyond
        the default precision::

            sage: z = R(5^20)
            sage: z
            0 + ...

        Without any help, Sage does not run the computation far enough to determine
        the valuation and an error is raised::

            sage: z.valuation()
            Traceback (most recent call last):
            ...
            PrecisionError: cannot determine the valuation; try to increase the halting precision

        By setting the argument ``halt``, one can force the computation to continue
        until a prescribed limit::

            sage: z.valuation(halt=20)   # not enough to find the valuation
            Traceback (most recent call last):
            ...
            PrecisionError: cannot determine the valuation; try to increase the halting precision

            sage: z.valuation(halt=21)   # now, we're okay
            20

        .. NOTE:

            It is also possible to pass in ``halt=False`` but it is not recommended
            because the computation can hang forever if this element were `0`.

        TESTS::

            sage: x = R.selfref()
            sage: (~x).valuation()
            Traceback (most recent call last):
            ...
            PrecisionError: no lower bound on the valuation is known
        """
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if halt is True:
            halt = self._valuation + self._parent.default_prec()
        elif halt is False:
            halt = maxordp
        else:
            halt = min(maxordp, halt)
        self._jump_relative_c(1, halt)
        if self._precbound >= maxordp and self._precrel == 0:
            raise PrecisionError("cannot determine the valuation; try to increase the halting precision")
        return Integer(self._valuation)

    def unit_part(self, halt=True):
        r"""
        Return the unit part of this element.

        INPUT:

        - ``halt`` -- an integer or a boolean (default: ``True``);
          the absolute precision after which the computation is abandonned
          if the first significant digit has not been found yet;
          if ``True``, the default precision of the parent is used;
          if ``False``, the computation is never abandonned

        EXAMPLES::

            sage: R = ZpL(5, 10)
            sage: a = R(20/21); a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.unit_part()
            4 + 4*5 + 5^3 + 4*5^5 + 3*5^6 + 4*5^7 + 5^9 + ...

            sage: b = 1/a; b
            4*5^-1 + 4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + ...
            sage: b.unit_part()
            4 + 4*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + ...

        The unit part of `0` is not defined::

            sage: R(0).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

            sage: R(0, 20).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

        See :meth:`valuation` for more details on the paramter ``halt``.

        """
        val = self.valuation(halt)
        if self._valuation >= self._precbound:
            raise ValueError("unit part of 0 not defined")
        return self >> val

    def val_unit(self, halt=True):
        r"""
        Return the valuation and the unit part of this element.

        INPUT:

        - ``halt`` -- an integer or a boolean (default: ``True``);
          the absolute precision after which the computation is abandonned
          if the first significant digit has not been found yet;
          if ``True``, the default precision of the parent is used;
          if ``False``, the computation is never abandonned

        EXAMPLES::

            sage: R = ZpL(5, 10)
            sage: a = R(20/21); a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.val_unit()
            (1, 4 + 4*5 + 5^3 + 4*5^5 + 3*5^6 + 4*5^7 + 5^9 + ...)

            sage: b = 1/a; b
            4*5^-1 + 4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + ...
            sage: b.val_unit()
            (-1, 4 + 4*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + ...)

        If this element is indistinguishable from zero, an error is raised
        since the unit part of `0` is not defined::

            sage: R(0).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

            sage: R(0, 20).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

        See :meth:`valuation` for more details on the paramter ``halt``.

        """
        val = self.valuation(halt)
        if self._valuation >= self._precbound:
            raise ValueError("unit part of 0 not defined")
        return val, self >> val

    def residue(self, absprec=1, field=True, check_prec=True):
        r"""
        Return the image of this element in the quotient
        `\ZZ/p^\mathrm{absprec}\ZZ`.

        INPUT:

        - ``absprec`` -- a non-negative integer (default: ``1``)

        - ``field`` -- boolean (default ``True``); when ``absprec`` is ``1``,
          whether to return an element of GF(p) or Zmod(p).

        - ``check_prec`` -- boolean (default ``True``); whether to raise an error 
          if this element has insufficient precision to determine the reduction.

        EXAMPLES::

            sage: R = ZpL(7, 10)
            sage: a = R(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...
            sage: a.residue()
            3
            sage: a.residue(2)
            45

        If this element has negative valuation, an error is raised::

            sage: K = R.fraction_field()
            sage: b = K(20/21)
            sage: b.residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have non-negative valuation in order to compute residue
        """
        if absprec >= maxordp:
            raise OverflowError
        if absprec < 0:
            raise ValueError("cannot reduce modulo a negative power of p")
        error = self._jump_c(absprec)
        if error:
            raise_error(error, not check_prec)
        if self._valuation < 0:
            raise ValueError("element must have non-negative valuation in order to compute residue")
        cdef celement digits
        self._getslice_relative(digits, 0, min(self._precrel, absprec - self._valuation))
        cdef Integer ans = element_get_sage(digits, self.prime_pow)
        if field and absprec == 1:
            return self._parent.residue_class_field()(ans)
        else:
            return self._parent.residue_ring(absprec)(ans)

    def lift(self, absprec=None):
        r"""
        Return a rational number which is congruent to this element modulo
        `p^\mathrm{prec}`.

        INPUT:

        - ``absprec`` -- an integer or ``None`` (default: ``None``); if ``None``,
          the absolute precision of this element is used

        EXAMPLES::

            sage: R = ZpL(7, 10)
            sage: a = R(1/2021, 5); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)
            sage: a.lift()
            15676
            sage: a.lift(2)
            45

        Here is another example with an element of negative valuation::

            sage: K = R.fraction_field()
            sage: b = K(20/21, 5); b
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + O(7^5)
            sage: b.lift()
            39223/7

        For unbounded elements, we must specify a precision::

            sage: c = R(1/2021)
            sage: c.lift()
            Traceback (most recent call last):
            ...
            ValueError: you must specify a precision for unbounded elements

            sage: c.lift(5)
            15676
        """
        if absprec is None:
            if self._precbound < maxordp:
                absprec = self._precbound
            else:
                raise ValueError("you must specify a precision for unbounded elements")
        else:
            absprec = Integer(absprec)
        cdef int error = self._jump_c(absprec)
        if error:
            raise_error(error)
        if absprec < self._valuation:
            return Integer(0)
        cdef celement digits
        self._getslice_relative(digits, 0, absprec - self._valuation)
        ans = element_get_sage(digits, self.prime_pow)
        if self._valuation:
            ans *= self._parent.prime() ** self._valuation
        return ans

    def __rshift__(self, s):
        cdef long start
        cdef long shift = long(s)
        if shift:
            if (<LazyElement>self)._parent.is_field():
                start = -maxordp
            else:
                start = shift
            return element_class_slice((<pAdicLazyElement>self)._parent, self, start, maxordp, shift)
        else:
            return self

    def __lshift__(self, s):
        return self.__rshift__(-s)

    cpdef _add_(self, other):
        if isinstance(self, LazyElement_zero):
            return other
        if isinstance(other, LazyElement_zero):
            return self
        return element_class_add(self._parent, self, <LazyElement>other)

    cpdef _sub_(self, other):
        if self is other:
            return self._parent.zero()
        if isinstance(self, LazyElement_zero):
            return other
        if isinstance(other, LazyElement_zero):
            return self
        return element_class_sub(self._parent, self, <LazyElement>other)

    cpdef _neg_(self):
        if isinstance(self, LazyElement_zero):
            return self
        return element_class_sub(self._parent, self._parent.zero(), self)

    cpdef _mul_(self, other):
        if isinstance(self, LazyElement_zero) or isinstance(other, LazyElement_one):
            return self
        if isinstance(self, LazyElement_one) or isinstance(other, LazyElement_zero):
            return other
        return element_class_mul(self._parent, self, <LazyElement>other)

    cpdef _div_(self, other):
        if isinstance(other, LazyElement_zero):
            raise ZeroDivisionError("cannot divide by zero")
        if (<LazyElement>other)._precbound < maxordp and (<LazyElement>other)._precrel == 0:
            raise ZeroDivisionError("cannot divide by something indistinguishable from zero")
        if isinstance(other, LazyElement_one):
            if self.prime_pow.in_field:
                return self
            else:
                return element_class_bound(self._parent.fraction_field(), self)
        return element_class_div(self._parent.fraction_field(), self, <LazyElement>other)

    def __invert__(self):
        if isinstance(self, LazyElement_zero):
            raise ZeroDivisionError("cannot divide by zero")
        if self._precbound < maxordp and self._precrel == 0:
            raise ZeroDivisionError("cannot divide by something indistinguishable from zero")
        if isinstance(self, LazyElement_one):
            return self
        return element_class_div(self._parent.fraction_field(), self._parent.one(), self)

    def inverse_of_unit(self):
        if isinstance(self, LazyElement_one):
            return self
        if self.prime_pow.in_field:
            if isinstance(self, LazyElement_zero):
                raise ZeroDivisionError("cannot divide by zero")
        else:
            if self.valuation() > 0:
                raise ZeroDivisionError("not a unit")
        return element_class_div(self._parent, self._parent.one(), self)        

    def sqrt(self):
        return element_class_sqrt(self._parent, self)

    def _test_pickling(self, **options):
        tester = self._tester(**options)
        from sage.misc.all import loads, dumps
        if self._precbound >= maxordp:
            tester.assertEqual(loads(dumps(self)), self.at_precision_relative())
        else:
            tester.assertEqual(loads(dumps(self)), self)

cdef class LazyElement_abandon(LazyElement):
    def __init__(self):
        self._valuation = -maxordp

    cdef int _next_c(self):
        return ERROR_ABANDON

cdef lazyelement_abandon = LazyElement_abandon()


cdef class LazyElement_init(LazyElement):
    def __cinit__(self):
        element_init(self._digits)

    def __dealloc__(self):
        element_clear(self._digits)

    cdef cdigit_ptr _getdigit_relative(self, long i):
        return element_get_digit(self._digits, i)

    cdef cdigit_ptr _getdigit_absolute(self, long i):
        return element_get_digit(self._digits, i - self._valuation)

    cdef void _getslice_relative(self, celement slice, long start, long length):
        element_get_slice(slice, self._digits, start, length)



# Assignment
############

# Zero

cdef class LazyElement_zero(LazyElement):
    def __init__(self, parent):
        LazyElement.__init__(self, parent)
        self._valuation = maxordp

    def __reduce__(self):
        return self.__class__, (self._parent,)

    cdef cdigit_ptr _getdigit_relative(self, long i):
        return digit_zero

    cdef cdigit_ptr _getdigit_absolute(self, long i):
        return digit_zero

    cdef void _getslice_relative(self, celement slice, long start, long length):
        element_init(slice)

    cdef int _jump_c(self, long prec):
        return 0

    cdef int _next_c(self):
        return 0


# One

cdef class LazyElement_one(LazyElement_init):
    def __init__(self, parent):
        LazyElement.__init__(self, parent)
        element_set_digit_ui(self._digits, 1, 0)
        self._precrel = maxordp

    def __reduce__(self):
        return self.__class__, (self._parent,)

    cdef int _jump_c(self, long prec):
        return 0

    cdef int _next_c(self):
        return 0


# Bound

cdef class LazyElement_bound(LazyElement):
    def __init__(self, parent, LazyElement x, precbound=None):
        LazyElement.__init__(self, parent)
        self._x = x
        if precbound is None:
            self._precbound = x._precbound
        else:
            self._precbound = min(x._precbound, precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._x, self._precbound)

    cdef cdigit_ptr _getdigit_relative(self, long i):
        return self._x._getdigit_relative(i)

    cdef cdigit_ptr _getdigit_absolute(self, long i):
        return self._x._getdigit_absolute(i)

    cdef void _getslice_relative(self, celement slice, long start, long length):
        self._x._getslice_relative(slice, start, length)

    cdef int _jump_c(self, long prec):
        cdef LazyElement x = self._x
        cdef int error
        if prec > self._precbound:
            error = ERROR_PRECISION | x._jump_c(self._precbound)
        else:
            error = x._jump_c(prec)
        self._precbound = min(self._precbound, x._precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        return error

    cdef int _next_c(self):
        cdef LazyElement x = self._x
        cdef int error = x._next_c()
        self._precbound = min(self._precbound, x._precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        return error


# Value

cdef class LazyElement_value(LazyElement_init):
    def __init__(self, parent, value, long shift=0, precbound=None):
        LazyElement.__init__(self, parent)
        element_set_digit_sage(self._digits, value, 0)
        self._value = value
        self._shift = shift
        self._valuation = -shift
        if precbound is not None and precbound is not Infinity:
            self._precbound = min(maxordp, precbound)
        self._valuebound = maxordp
        self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._value, self._shift, self._precbound)

    cdef int _jump_c(self, long prec):
        if self._valuebound >= maxordp:
            return LazyElement._jump_c(self, prec)
        if (self._precbound is not None) and (prec > self._precbound):
            self._precrel = self._precbound - self._valuation
            return ERROR_PRECISION
        cdef long precrel = min(prec, maxordp) - self._valuation
        if precrel > self._precrel:
            self._precrel = precrel
        if prec >= maxordp:
            return ERROR_OVERFLOW
        return 0

    cdef int _next_c(self):
        # The algorithm is not optimal (quadratic in the precision),
        # but it sounds okay
        element_reduce_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        if digit_is_zero(self._getdigit_relative(self._precrel)):
            self._valuebound = self._valuation + self._precrel
            if self._precrel == 0:
                self._valuation = self._precbound
            elif self._precbound >= maxordp:
                self._precrel = maxordp
            else:
                self._precrel = self._precbound - self._valuation
        return 0


# Random

cdef class LazyElement_random(LazyElement_init):
    def __init__(self, parent, valuation, precbound, seed=None):
        LazyElement.__init__(self, parent)
        if seed is None:
            self._seed = randint(0, 2*maxordp)
        else:
            self._seed = seed
        digit_random_init(self._generator, self._seed)
        if valuation is None:
            self._initialvaluation = ZZ.random_element()
        else:
            self._initialvaluation = valuation
        self._valuation = self._initialvaluation
        if precbound is not None:
            self._precbound = min(maxordp, precbound)
        self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._initialvaluation, self._precbound, self._seed)

    cdef int _next_c(self):
        cdef cdigit r
        digit_random(r, self.prime_pow, self._generator)
        if self._precrel == 0 and digit_is_zero(r):
            self._valuation += 1
            digit_clear(r)
        else:
            element_set_digit(self._digits, r, self._precrel)
            self._precrel += 1
        return 0


# Operations
############

# Slice and shift

cdef class LazyElement_slice(LazyElement):
    def __init__(self, parent, LazyElement x, long start, long stop, long shift):
        # self[0] = x[shift]
        LazyElement.__init__(self, parent)
        self._x = x
        self._start = start
        self._shift = shift
        self._valuation = max(x._valuation, start) - shift
        self._precrel = min(x._precrel + x._valuation, stop) - self._valuation - shift
        if x._precbound < stop:
            self._precbound = min(maxordp, x._precbound - shift)
        self._stop = stop
        if self._precrel < 0:
            self._precrel = 0
        while self._precrel > 0 and digit_is_zero(self._getdigit_relative(0)):
            self._precrel -= 1
            self._valuation += 1
        self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._x, self._start, self._stop, self._shift)

    cdef cdigit_ptr _getdigit_relative(self, long i):
        return self._getdigit_absolute(i + self._valuation)

    cdef cdigit_ptr _getdigit_absolute(self, long i):
        cdef long j = i + self._shift
        if j < self._start or j >= self._stop:
            return digit_zero
        else:
            return self._x._getdigit_absolute(j)

    cdef void _getslice_relative(self, celement slice, long start, long length):
        cdef LazyElement x = self._x
        cdef long s = start + self._valuation + self._shift
        cdef long start_absolute = max(self._start, s)
        cdef long stop_absolute = min(self._stop, s + length)
        x._getslice_relative(slice, start_absolute - x._valuation, stop_absolute - start_absolute)

    cdef int _jump_c(self, long prec):
        cdef LazyElement x = self._x
        cdef int error = 0
        cdef long pr
        if prec <= self._valuation + self._precrel:
            return 0
        if prec > self._precbound:
            prec = self._precbound
            error = ERROR_PRECISION
        cdef int errorx = x._jump_c(min(prec + self._shift, self._stop))
        pr = max(self._valuation, x._valuation + x._precrel - self._shift)
        if self._precrel == 0:
            while self._valuation < pr and digit_is_zero(self._getdigit_relative(0)):
                self._valuation += 1
        self._precrel = prec - self._valuation
        if errorx:
            return errorx
        else:
            return error

    cdef int _next_c(self):
        cdef int error
        cdef long n = self._precrel + self._valuation
        n += self._shift
        if n <= self._stop:
            error = self._x._jump_c(n+1)
            if error:
                return error
        if self._precrel == 0 and (n > self._stop or digit_is_zero(self._getdigit_relative(0))):
            self._valuation += 1
        else:
            self._precrel += 1
        return 0


# Addition

cdef class LazyElement_add(LazyElement_init):
    def __init__(self, parent, LazyElement x, LazyElement y):
        LazyElement.__init__(self, parent)
        self._x = x
        self._y = y
        self._valuation = min(x._valuation, y._valuation)
        self._precbound = min(x._precbound, y._precbound)
        self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._x, self._y)

    cdef int _jump_c(self, long prec):
        # We reimplement _jump_c for better performances
        cdef long n = self._valuation + self._precrel
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef int error = x._jump_c(prec) | y._jump_c(prec)
        prec = min(prec, x._valuation + x._precrel, y._valuation + y._precrel)
        while n < prec:
            element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
            element_iadd_digit(self._digits, y._getdigit_absolute(n), self._precrel)
            element_reducesmall_digit(self._digits, self._precrel, self.prime_pow)
            if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
                self._valuation += 1
                element_shift_right(self._digits)
            else:
                self._precrel += 1
            n += 1
        return error

    cdef int _next_c(self):
        cdef long n = self._valuation + self._precrel
        if n >= self._valuebound:
            if self._precrel == 0:
                self._valuation = self._precbound
            elif self._precbound < maxordp:
                self._precrel = self._precbound - self._valuation
            else:
                self._precrel = maxordp
            return 0
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef int error = x._jump_c(n+1) | y._jump_c(n+1)
        if error:
            return error
        element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
        element_iadd_digit(self._digits, y._getdigit_absolute(n), self._precrel)
        element_reducesmall_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        return 0


# Subtraction

cdef class LazyElement_sub(LazyElement_init):
    def __init__(self, parent, LazyElement x, LazyElement y):
        LazyElement.__init__(self, parent)
        self._x = x
        self._y = y
        self._valuation = min(x._valuation, y._valuation)
        self._precbound = min(x._precbound, y._precbound)
        self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._x, self._y)

    cdef int _jump_c(self, long prec):
        # We reimplement _jump_c for better performances
        cdef long n = self._valuation + self._precrel
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef int error = x._jump_c(prec) | y._jump_c(prec)
        prec = min(prec, x._valuation + x._precrel, y._valuation + y._precrel)
        while n < prec:
            element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
            element_isub_digit(self._digits, y._getdigit_absolute(n), self._precrel)
            element_reduceneg_digit(self._digits, self._precrel, self.prime_pow)
            if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
                self._valuation += 1
                element_shift_right(self._digits)
            else:
                self._precrel += 1
            n += 1
        return error

    cdef int _next_c(self):
        cdef long n = self._valuation + self._precrel
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef int error = x._jump_c(n+1) | y._jump_c(n+1)
        if error:
            return error
        element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
        element_isub_digit(self._digits, y._getdigit_absolute(n), self._precrel)
        element_reduceneg_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        return 0



# Multiplication

cdef class LazyElement_mul(LazyElement_init):
    def __cinit__(self):
        digit_init(self._lastdigit_x)
        digit_init(self._lastdigit_y)

    def __dealloc__(self):
        digit_clear(self._lastdigit_x)
        digit_clear(self._lastdigit_y)

    def __init__(self, parent, LazyElement x, LazyElement y):
        LazyElement.__init__(self, parent)
        self._x = x
        self._y = y
        self._valuation = min(maxordp, x._valuation + y._valuation)
        if x._precbound < maxordp:
            y._jump_relative_c(1, y._valuation + self._parent.default_prec())
            self._precbound = min(self._precbound, y._valuation + x._precbound)
        if y._precbound < maxordp:
            x._jump_relative_c(1, x._valuation + self._parent.default_prec())
            self._precbound = min(self._precbound, x._valuation + y._precbound)
        self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._x, self._y)

    cdef int _next_c(self):
        global tmp_digit, tmp_poly
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef long n = self._valuation + self._precrel

        cdef int errorx = x._jump_c(n - y._valuation + 1)
        cdef int errory = y._jump_c(n - x._valuation + 1)
        cdef int error = errorx | errory
        if self._precrel == 0:
            self._valuation = x._valuation + y._valuation
            if self._valuation > n:
                return 0
            if self._valuation < n or x._precrel == 0 or y._precrel == 0:
                return error | ERROR_PRECISION
        elif error:
            return error

        n = self._precrel
        digit_set(self._lastdigit_x, x._getdigit_relative(n))
        digit_set(self._lastdigit_y, y._getdigit_relative(n))
        digit_mul(tmp_digit, x._getdigit_relative(0), self._lastdigit_y)
        element_iadd_digit(self._digits, tmp_digit, n)
        if n:
            digit_mul(tmp_digit, self._lastdigit_x, y._getdigit_relative(0))
            element_iadd_digit(self._digits, tmp_digit, n)

        cdef long m = n + 2
        cdef long len = 1
        cdef celement slicex, slicey
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
            x._getslice_relative(slicex, len - 1, len)
            y._getslice_relative(slicey, (m-1)*len - 1, len)
            element_mul(tmp_poly, slicex, slicey)
            element_iadd_slice(self._digits, tmp_poly, n)
            if m > 2:
                x._getslice_relative(slicex, (m-1)*len - 1, len)
                y._getslice_relative(slicey, len - 1, len)
                element_mul(tmp_poly, slicex, slicey)
                element_iadd_slice(self._digits, tmp_poly, n)

        element_reduce_digit(self._digits, n, self.prime_pow)
        self._precrel += 1
        return 0

    cdef int _update_last_digit(self):
        if self._precrel == 0:
            return ERROR_UNEXPECTED
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef long n = self._precrel - 1
        cdef long m = n + 2
        cdef long len = 2
        cdef celement slice
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
        len -= 1

        digit_sub(tmp_digit, x._getdigit_relative(n), self._lastdigit_x)
        y._getslice_relative(slice, 0, len)
        element_scalarmul(tmp_poly, slice, tmp_digit)
        element_iadd_slice(self._digits, tmp_poly, n)
        if m == 2:
            len -= 1
        digit_sub(tmp_digit, y._getdigit_relative(n), self._lastdigit_y)
        x._getslice_relative(slice, 0, len)
        element_scalarmul(tmp_poly, slice, tmp_digit)
        element_iadd_slice(self._digits, tmp_poly, n)
        if m == 2:
            digit_mul(tmp_digit, tmp_digit, self._lastdigit_x)
            element_iadd_digit(self._digits, tmp_digit, 2*len)
        element_reduce_digit(self._digits, n, self.prime_pow)

        digit_set(self._lastdigit_x, x._getdigit_relative(n))
        digit_set(self._lastdigit_y, y._getdigit_relative(n))
        return 0


cdef class LazyElement_muldigit(LazyElement_init):
    # compute x*y where x is assumed to be in (0,p)
    def __init__(self, parent, LazyElement_div x, LazyElement y):
        LazyElement.__init__(self, parent)
        self._x = <cdigit_ptr>x._inverse
        self._y = y
        self._valuation = y._valuation
        self._init_jump()

    cdef int _next_c(self):
        cdef long n = self._valuation + self._precrel
        cdef int error = self._y._jump_c(n+1)
        if error:
            return error
        digit_mul(tmp_digit, self._x, self._y._getdigit_absolute(n))
        element_iadd_digit(self._digits, tmp_digit, self._precrel)
        element_reduce_digit(self._digits, self._precrel, self.prime_pow)
        self._precrel += 1
        return 0

    cdef void _erase_first_digit(self):
        # This is a awful hack,
        # but it's very useful for division
        if self._precrel:
            self._valuation += 1
            self._precrel -= 1
            element_shift_right(self._digits)


# Division

cdef class LazyElement_div(LazyElement_init):
    def __cinit__(self):
        digit_init(self._inverse)

    def __dealloc__(self):
        digit_clear(self._inverse)

    def __init__(self, parent, LazyElement num, LazyElement denom, precbound=None):
        LazyElement.__init__(self, parent)
        self._num = num
        self._denom = denom
        if denom._valuation <= -maxordp:
            self._maxprec = maxordp + 1
        else:
            self._maxprec = denom._valuation + max(1, self._parent.default_prec())
        if precbound is not None:
            self._precbound = min(maxordp, precbound)
        cdef int error = self._bootstrap_c()
        if num._precbound < maxordp:
            self._precbound = min(self._precbound, num._precbound - denom._valuation)
        if denom._precbound < maxordp:
            self._precbound = min(self._precbound, denom._precbound + num._valuation - 2*denom._valuation)
        if error:
            self._valuation = -maxordp
        else:
            self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._num, self._denom, self._precbound)

    cdef int _bootstrap_c(self):
        cdef int error
        cdef LazyElement num = self._num
        cdef LazyElement denom = self._denom

        while denom._valuation < self._maxprec and denom._precrel == 0:
            error = denom._next_c()
            if error:
                if error & ERROR_PRECISION:
                    error |= ERROR_DIVISION
                return error
            if self._maxprec > maxordp and denom._valuation > -maxordp:
                self._maxprec = denom._valuation + max(1, self._parent.default_prec())
        if denom._precrel == 0:
            return ERROR_ABANDON

        self._valuation = num._valuation - denom._valuation
        digit_inv(self._inverse, denom._getdigit_relative(0), self.prime_pow)
        self._definition = lazyelement_abandon
        cdef parent = self._parent
        cdef LazyElement a = element_class_muldigit(parent, self, num)
        cdef LazyElement b = element_class_muldigit(parent, self, denom)
        cdef LazyElement c = element_class_slice(parent, b, denom._valuation + 1, maxordp, 0)
        cdef LazyElement d = element_class_mul(parent, c, self)
        self._definition = element_class_sub(parent, a, d)
        return 0

    cdef int _next_c(self):
        cdef LazyElement definition = self._definition
        cdef long val
        if definition is None:
            error = self._bootstrap_c()
            if error:
                return error
        else:
            val = self._valuation + self._denom._valuation
            error = definition._jump_c(val + self._precrel + 1)
            if error:
                return error
            if definition._valuation > val:
                self._valuation = min(self._precbound, definition._valuation - self._denom._valuation)
                self._precbound = min(self._precbound, definition._precbound - self._denom._valuation)
            else:
                digit = definition._getdigit_relative(self._precrel)
                element_set_digit(self._digits, digit, self._precrel)
                self._precrel += 1
        return 0


# Square root

cdef class LazyElement_sqrt(LazyElement_init):
    def __init__(self, parent, LazyElement x):
        LazyElement.__init__(self, parent)
        self._x = x
        if x._valuation <= -maxordp:
            self._valuation = -maxordp
            self._maxprec = maxordp + 1
        else:
            self._valuation = x._valuation >> 1
            self._maxprec = x._valuation + 2*self._parent.default_prec()
        cdef int error = self._bootstrap_c()
        if error & ERROR_NOTSQUARE:
            raise ValueError("not a square")
        if not error:
            if x._precbound < maxordp:
                self._precbound = x._precbound - x._valuation / 2
                if self._parent.prime() == 2:
                    self._precbound -= 1
                self._precbound = min(maxordp, self._precbound)
            self._init_jump()

    def __reduce__(self):
        return self.__class__, (self._parent, self._x)

    cdef int _bootstrap_c(self):
        cdef LazyElement x = self._x
        while x._valuation < self._maxprec and x._precrel == 0:
            error = x._next_c()
            if error:
                return error
            if self._maxprec > maxordp and x._valuation > -maxordp:
                self._maxprec = x._valuation + 2*self._parent.default_prec()
        if x._valuation > -maxordp:
            if x._valuation & 1:
                self._valuation = (x._valuation + 1) >> 1
            else:
                self._valuation = x._valuation >> 1
        if x._precrel == 0:
            return ERROR_ABANDON
        if x._valuation & 1 != 0:
            return ERROR_NOTSQUARE

        cdef parent = self._parent
        cdef long val = self._valuation
        cdef cdigit digit
        cdef Integer zd, p = self.prime_pow.prime
        cdef LazyElement u, y, c, d

        if p == 2:
            element_set_digit_ui(self._digits, 1, 0)
            self._precrel = 1
            if x._precrel == 1:
                error = x._next_c()
                if error:
                    return error
            if not digit_equal_ui(x._getdigit_relative(1), 0):
                return ERROR_NOTSQUARE
            if x._precrel == 2:
                error = x._next_c()
                if error:
                    return error
            if not digit_equal_ui(x._getdigit_relative(2), 0):
                return ERROR_NOTSQUARE
            zd = Integer(1)
            self._definition = lazyelement_abandon
            u = element_class_slice(parent, self, val + 2, maxordp, val)
            y = element_class_slice(parent, x, -maxordp, maxordp, val + 1)
            c = element_class_value(parent, zd, shift=-val+1) 
            d = element_class_slice(parent, u*u, -maxordp, maxordp, -val + 1)
            self._definition = y + c - d
        else:
            digit_init(digit)
            if digit_sqrt(digit, x._getdigit_relative(0), self.prime_pow):
                digit_clear(digit)
                return ERROR_NOTSQUARE
            element_set_digit(self._digits, digit, 0)
            self._precrel = 1
            zd = digit_get_sage(digit)
            self._definition = lazyelement_abandon
            u = element_class_slice(parent, self, val + 1, maxordp, val)
            y = element_class_slice(parent, x, -maxordp, maxordp, 2*val)
            c = element_class_value(parent, zd*zd)
            d = element_class_value(parent, 2*zd, shift=val)
            self._definition = (y + c - u*u) / d
        return 0

    cdef int _next_c(self):
        cdef LazyElement x = self._x
        cdef LazyElement definition = self._definition
        cdef long n = self._valuation + self._precrel
        cdef int error
        if definition is None:
            error = self._bootstrap_c()
            if error:
                return error
        else:
            error = definition._jump_c(n+1)
            if error:
                return error
            element_set_digit(self._digits, definition._getdigit_relative(self._precrel), self._precrel)
            self._precrel += 1
        return 0

    def next(self):
        return self._next_c()


# Teichmüller lifts

cdef class LazyElement_teichmuller(LazyElement_init):
    def __init__(self, parent, Integer xbar):
        LazyElement.__init__(self, parent)
        cdef cdigit digit
        digit_init(digit)
        digit_set_sage(digit, xbar)
        digit_mod(digit, digit, self.prime_pow)
        if digit_equal_ui(digit, 0):
            digit_clear(digit)
            self._trivial = True
            self._valuation = maxordp
        else:
            element_set_digit(self._digits, digit, 0)
            self._trivial = digit_equal_ui(digit, 1)
        self._precrel += 1

        cdef LazyElement xn
        cdef Integer p
        cdef int size, i        
        if not self._trivial:
            xn = self
            p = self.prime_pow.prime
            size = mpz_sizeinbase(p.value, 2)
            i = size - 2
            self._xns = [ ]
            while i >= 0:
                xn = element_class_mul(parent, xn, xn)
                self._xns.append(xn)
                if mpz_tstbit(p.value, i):
                    xn = element_class_mul(parent, xn, self)
                    self._xns.append(xn)
                i -= 1
            self._xp = xn
        self._ready = True

    def __reduce__(self):
        xbar = digit_get_sage(element_get_digit(self._digits, 0))
        return self.__class__, (self._parent, xbar)

    cdef int _jump_c(self, long prec):
        if not self._ready:
            return ERROR_ABANDON
        if self._trivial:
            if self._valuation == 0 and self._precrel < prec:
                self._precrel = prec
            return 0
        return LazyElement._jump_c(<LazyElement>self, prec)

    cdef int _next_c(self):
        if self._trivial:
            if self._valuation:
                self._precrel += 1
            return 0
        cdef int error
        cdef LazyElement xp = self._xp
        cdef LazyElement_mul xn
        self._precrel += 1
        xp._jump_c(self._precrel)
        element_set_digit(self._digits, xp._getdigit_relative(self._precrel - 1), self._precrel - 1)
        for xn in self._xns:
            error = xn._update_last_digit()
            if error:
                return error | ERROR_UNEXPECTED
        return 0


# Self-referent definitions
###########################

cdef class LazyElement_selfref(LazyElement_init):
    def __init__(self, parent, long valuation, digits=None):
        LazyElement.__init__(self, parent)
        if valuation >= maxordp:
            raise OverflowError("valuation is too large (maximum is %s)" % maxordp)
        self._valuation = valuation
        self._definition = None
        self._next = maxordp
        cdef cdigit digit
        if digits is not None:
            digits = [ Integer(d) for d in digits ]
            for d in digits:
                digit_set_sage(digit, d)
                element_iadd_digit(self._digits, digit, self._precrel)
                element_reduce_digit(self._digits, self._precrel, self.prime_pow)
                if self._precrel == 0 and digit_is_zero(element_get_digit(self._digits, 0)):
                    self._valuation += 1
                    element_shift_right(self._digits)
                else:
                    self._precrel += 1
        self._initialvaluation = self._valuation
        self._initialprecrel = self._precrel

    def __reduce__(self):
        digits = [ ]
        for i in range(self._initialprecrel):
            digits.append(digit_get_sage(element_get_digit(self._digits, i)))
        definition = None
        if id(self) not in persist.already_pickled:
            persist.already_pickled[id(self)] = True
            definition = self._definition
        return unpickle_selfref, (id(self), self.__class__, self._parent, self._initialvaluation, digits, definition)

    cpdef set(self, LazyElement definition):
        if self._definition is not None:
            raise ValueError("this self-referent number is already defined")
        self._definition = definition
        self._precbound = max(self._valuation + self._precrel, definition._precbound)
        self._init_jump()

    cdef int _next_c(self):
        cdef LazyElement definition = self._definition
        cdef cdigit_ptr digit
        cdef long n = self._valuation + self._precrel
        cdef long diffval
        cdef int error

        if definition is None:
            return ERROR_NOTDEFINED
        if n >= self._next:
            return ERROR_CIRCULAR

        cdef long svenext = self._next
        self._next = n
        error = definition._jump_c(n+1)
        if not error:
            digit = definition._getdigit_absolute(n)
            if self._precrel == 0 and digit_is_zero(digit):
                self._valuation += 1
            else:
                element_set_digit(self._digits, digit, self._precrel)
                self._precrel += 1
        self._next = svenext
        return error

def unpickle_selfref(uid, cls, parent, valuation, digits, definition):
    if uid in persist.already_unpickled:
        elt = persist.already_unpickled[uid]
    else:
        elt = cls(parent, valuation, digits)
        persist.already_unpickled[uid] = elt
    if definition is not None:
        elt.set(definition)
    return elt


# Expansion
###########

cdef class ExpansionIter(object):
    cdef LazyElement elt
    cdef long start
    cdef long stop
    cdef long current

    def __init__(self, LazyElement elt, start, stop):
        self.elt = elt
        self.start = start
        self.stop = stop
        self.current = self.start

    def __repr__(self):
        return "%s-adic expansion of %s" % (self.elt._parent.prime(), self.elt)

    def __len__(self):
        if self.stop >= maxordp:
            raise NotImplementedError("infinite sequence")
        return Integer(self.stop - self.start)

    def __iter__(self):
        return self

    def __next__(self):
        if self.current < self.stop:
            digit = self.elt.digit(self.current)
            self.current += 1
            return digit
        else:
            raise StopIteration
