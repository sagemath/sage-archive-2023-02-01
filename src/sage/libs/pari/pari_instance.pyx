"""
PARI C-library interface

AUTHORS:

- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta

- William Stein (2006-03-06): added newtonpoly

- Justin Walker: contributed some of the function definitions

- Gonzalo Tornaria: improvements to conversions; much better error
  handling.

- Robert Bradshaw, Jeroen Demeyer, William Stein (2010-08-15):
  Upgrade to PARI 2.4.3 (#9343)

- Jeroen Demeyer (2011-11-12): rewrite various conversion routines
  (#11611, #11854, #11952)

- Peter Bruin (2013-11-17): split off this file from gen.pyx (#15185)


EXAMPLES::

    sage: pari('5! + 10/x')
    (120*x + 10)/x
    sage: pari('intnum(x=0,13,sin(x)+sin(x^2) + x)')
    85.1885681951527
    sage: f = pari('x^3-1')
    sage: v = f.factor(); v
    [x - 1, 1; x^2 + x + 1, 1]
    sage: v[0]   # indexing is 0-based unlike in GP.
    [x - 1, x^2 + x + 1]~
    sage: v[1]
    [1, 1]~

Arithmetic obeys the usual coercion rules::

    sage: type(pari(1) + 1)
    <type 'sage.libs.pari.gen.gen'>
    sage: type(1 + pari(1))
    <type 'sage.libs.pari.gen.gen'>

GUIDE TO REAL PRECISION AND THE PARI LIBRARY

The default real precision in communicating with the Pari library
is the same as the default Sage real precision, which is 53 bits.
Inexact Pari objects are therefore printed by default to 15 decimal
digits (even if they are actually more precise).

Default precision example (53 bits, 15 significant decimals)::

    sage: a = pari(1.23); a
    1.23000000000000
    sage: a.sin()
    0.942488801931698

Example with custom precision of 200 bits (60 significant
decimals)::

    sage: R = RealField(200)
    sage: a = pari(R(1.23)); a   # only 15 significant digits printed
    1.23000000000000
    sage: R(a)         # but the number is known to precision of 200 bits
    1.2300000000000000000000000000000000000000000000000000000000
    sage: a.sin()      # only 15 significant digits printed
    0.942488801931698
    sage: R(a.sin())   # but the number is known to precision of 200 bits
    0.94248880193169751002382356538924454146128740562765030213504

It is possible to change the number of printed decimals::

    sage: R = RealField(200)    # 200 bits of precision in computations
    sage: old_prec = pari.set_real_precision(60)  # 60 decimals printed
    sage: a = pari(R(1.23)); a
    1.23000000000000000000000000000000000000000000000000000000000
    sage: a.sin()
    0.942488801931697510023823565389244541461287405627650302135038
    sage: pari.set_real_precision(old_prec)  # restore the default printing behavior
    60

Unless otherwise indicated in the docstring, most Pari functions
that return inexact objects use the precision of their arguments to
decide the precision of the computation. However, if some of these
arguments happen to be exact numbers (integers, rationals, etc.),
an optional parameter indicates the precision (in bits) to which
these arguments should be converted before the computation. If this
precision parameter is missing, the default precision of 53 bits is
used. The following first converts 2 into a real with 53-bit
precision::

    sage: R = RealField()
    sage: R(pari(2).sin())
    0.909297426825682

We can ask for a better precision using the optional parameter::

    sage: R = RealField(150)
    sage: R(pari(2).sin(precision=150))
    0.90929742682568169539601986591174484270225497

Warning regarding conversions Sage - Pari - Sage: Some care must be
taken when juggling inexact types back and forth between Sage and
Pari. In theory, calling p=pari(s) creates a Pari object p with the
same precision as s; in practice, the Pari library's precision is
word-based, so it will go up to the next word. For example, a
default 53-bit Sage real s will be bumped up to 64 bits by adding
bogus 11 bits. The function p.python() returns a Sage object with
exactly the same precision as the Pari object p. So
pari(s).python() is definitely not equal to s, since it has 64 bits
of precision, including the bogus 11 bits. The correct way of
avoiding this is to coerce pari(s).python() back into a domain with
the right precision. This has to be done by the user (or by Sage
functions that use Pari library functions in gen.pyx). For
instance, if we want to use the Pari library to compute sqrt(pi)
with a precision of 100 bits::

    sage: R = RealField(100)
    sage: s = R(pi); s
    3.1415926535897932384626433833
    sage: p = pari(s).sqrt()
    sage: x = p.python(); x  # wow, more digits than I expected!
    1.7724538509055160272981674833410973484
    sage: x.prec()           # has precision 'improved' from 100 to 128?
    128
    sage: x == RealField(128)(pi).sqrt()  # sadly, no!
    False
    sage: R(x)               # x should be brought back to precision 100
    1.7724538509055160272981674833
    sage: R(x) == s.sqrt()
    True

Elliptic curves and precision: If you are working with elliptic
curves and want to compute with a precision other than the default
53 bits, you should use the precision parameter of ellinit()::

    sage: R = RealField(150)
    sage: e = pari([0,0,0,-82,0]).ellinit(precision=150)
    sage: eta1 = e.elleta()[0]
    sage: R(eta1)
    3.6054636014326520859158205642077267748102690

Number fields and precision: TODO

TESTS:

Check that output from PARI's print command is actually seen by
Sage (ticket #9636)::

    sage: pari('print("test")')
    test

"""
