r"""
Tests for the Sage <-> PARI interface

Deprecation checks::

    sage: pari.poltchebi(10)
    doctest:...: DeprecationWarning: poltchebi is deprecated. Please use polchebyshev instead.
    See http://trac.sagemath.org/18203 for details.
    512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
    sage: pari.nth_prime(10)
    doctest:...: DeprecationWarning: nth_prime is deprecated. Please use prime instead.
    See http://trac.sagemath.org/20216 for details.
    29
    sage: pari.prime_list(10)
    doctest:...: DeprecationWarning: prime_list is deprecated. Please use primes instead.
    See http://trac.sagemath.org/20216 for details.
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    sage: pari.primes_up_to_n(20)
    doctest:...: DeprecationWarning: pari.primes_up_to_n(n) is deprecated, use pari.primes(end=n) instead
    See http://trac.sagemath.org/20216 for details.
    [2, 3, 5, 7, 11, 13, 17, 19]
"""
