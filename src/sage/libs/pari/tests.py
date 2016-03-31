r"""
Tests for the Sage <-> PARI interface

Deprecation checks::

    sage: pari.poltchebi(10)
    doctest:...: DeprecationWarning: poltchebi is deprecated. Please use polchebyshev instead.
    See http://trac.sagemath.org/18203 for details.
    512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
    sage: pari("x^3 + 1").polsturm(-1, 1)
    doctest:...: DeprecationWarning: argument 2 of the PARI/GP function polsturm is undocumented and deprecated
    1
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
    sage: pari.polcyclo_eval(8, 2)
    doctest:...: DeprecationWarning: polcyclo_eval is deprecated. Please use polcyclo instead.
    See http://trac.sagemath.org/20217 for details.
    17

A long list of doctests which used to be part of manually written code
which is now automatically generated:

Create a gp file::

    sage: import tempfile
    sage: gpfile = tempfile.NamedTemporaryFile(mode="w")
    sage: gpfile.file.write("mysquare(n) = {\n")
    sage: gpfile.file.write("    n^2;\n")
    sage: gpfile.file.write("}\n")
    sage: gpfile.file.write("polcyclo(5)\n")
    sage: gpfile.file.flush()

Read it in Sage, we get the result of the last line::

    sage: pari.read(gpfile.name)
    x^4 + x^3 + x^2 + x + 1

Call the function defined in the gp file::

    sage: pari('mysquare(12)')
    144

Constants::

    sage: pari.euler()
    0.577215664901533
    sage: pari.euler(precision=100).python()
    0.577215664901532860606512090082...
    sage: pari.pi()
    3.14159265358979
    sage: pari.pi(precision=100).python()
    3.1415926535897932384626433832...

Polynomial functions::

    sage: pari.pollegendre(7)
    429/16*x^7 - 693/16*x^5 + 315/16*x^3 - 35/16*x
    sage: pari.pollegendre(7, 'z')
    429/16*z^7 - 693/16*z^5 + 315/16*z^3 - 35/16*z
    sage: pari.pollegendre(0)
    1

    sage: pari.polcyclo(8)
    x^4 + 1
    sage: pari.polcyclo(7, 'z')
    z^6 + z^5 + z^4 + z^3 + z^2 + z + 1
    sage: pari.polcyclo(1)
    x - 1

Random seed::

    sage: a = pari.getrand()
    sage: a.type()
    't_INT'
"""
