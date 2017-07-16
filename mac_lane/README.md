Mac Lane's Algorithms in Sage
=============================
This package implements most of Mac Lane's algorithms [1,2] and related
structures to represent discrete valuations and discrete pseudo-valuations on
rings in Sage.

The package should run on an unmodified Sage 8.0 and can be imported with
```
sage: from mac_lane import *
```

To run the included tests, execute `sage -tp --optional=sage,standalone mac_lane/`.

High-Level Interface
====================
Valuations can be defined conveniently on some Sage rings such as p-adic rings and function fields.

p-adic valuations
-----------------
Valuations on number fields can be easily specified if they uniquely extend the valuation of a rational prime:
```
sage: v = pAdicValuation(QQ, 2)
sage: v(1024)
10
```

They are normalized such that the rational prime has valuation 1.
```
sage: K.<a> = NumberField(x^2 + x + 1)
sage: v = pAdicValuation(K, 2)
sage: v(1024)
10
```

If there are multiple valuations over a prime, they can be obtained by extending a valuation from a smaller ring.
```
sage: K.<a> = NumberField(x^2 + x + 1)
sage: v = pAdicValuation(QQ, 7)
sage: v.extensions(K)
[[ 7-adic valuation, v(x + 3) = 1 ]-adic valuation,
[ 7-adic valuation, v(x + 5) = 1 ]-adic valuation]
sage: w,ww = _
sage: w(a + 3), ww(a + 3)
(1, 0)
sage: w(a + 5), ww(a + 5)
(0, 1)
```

Function Field valuations
-------------------------
Similarly, valuations can be defined on function fields:
```
sage: K.<x> = FunctionField(QQ)
sage: v = FunctionFieldValuation(K, x)
sage: v(1/x)
-1

sage: v = FunctionFieldValuation(K, 1/x)
sage: v(1/x)
1
```

On extensions of function fields, valuations can be specified explicitly by
providing a prime on the underlying rational function field when the extension
is unique::
```
sage: K.<x> = FunctionField(QQ)
sage: R.<y> = K[]
sage: L.<y> = K.extension(y^2 - x)
sage: v = FunctionFieldValuation(L, x)
sage: v(x)
1
```

Valuations can also be extended from smaller function fields.
```
sage: K.<x> = FunctionField(QQ)
sage: v = FunctionFieldValuation(K, x - 4)
sage: R.<y> = K[]
sage: L.<y> = K.extension(y^2 - x)
sage: v.extensions(L)
[ (x - 4)-adic valuation, v(y - 2) = 1 ]-adic valuation,
 [ (x - 4)-adic valuation, v(y + 2) = 1 ]-adic valuation]
```

Low-Level Interface
===================

Mac Lane valuations
-------------------
Internally, all the above is backed by the algorithms described in [1,2]. Let
us consider the extensions of `FunctionFieldValuation(K, x - 4)` to the field
`L` above to outline how this works internally.

First, the valuation on K is induced by a valuation on ℚ[x]. To construct this
valuation, we start from the trivial valuation on ℚ and consider its induced
Gauss valuation on ℚ[x], i.e., the valuation that assigns to a polynomial the
minimum of the coefficient valuations.
```
sage: R.<x> = QQ[]
sage: v = TrivialValuation(QQ)
sage: v = GaussValuation(R, v)
```
The Gauss valuation can be augmented by specifying that x - 4 has valuation 1.
```
sage: v = v.augmentation(x - 4, 1); v
[ Gauss valuation induced by Trivial valuation on Rational Field, v(x - 4) = 1 ]
```

This valuation then extends uniquely to the fraction field.
```
sage: K.<x> = FunctionField(QQ)
sage: v = v.extension(K); v
(x - 4)-adic valuation
```

Over the function field we repeat the above process, i.e., we define the Gauss
valuation induced by it and augment it to approximate an extension to L.
```
sage: R.<y> = K[]
sage: w = GaussValuation(R, v)
sage: w = w.augmentation(y - 2, 1); w
[ Gauss valuation induced by (x - 4)-adic valuation, v(y - 2) = 1 ]
sage: L.<y> = K.extension(y^2 - x)
sage: ww = w.extension(L); ww
[ (x - 4)-adic valuation, v(y - 2) = 1 ]-adic valuation
```

Limit valuations
----------------
In the previous example the final valuation `ww` is not merely given by evaluating `w` on the ring `K[y]`.
```
sage: ww(y^2 - x)
+Infinity
sage: y = R.gen()
sage: w(y^2 - x)
1
```

Instead `ww` is given by a limit, i.e., an infinite sequence of augmentations of valuations.
```
sage: ww._base_valuation
[ Gauss valuation induced by (x - 4)-adic valuation, v(y - 2) = 1 , … ]
```

The terms of this infinite sequence are computed on demand.
```
sage: ww._base_valuation._approximation
[ Gauss valuation induced by (x - 4)-adic valuation, v(y - 2) = 1 ]
sage: ww(y - 1/4*x - 1)
2
sage: ww._base_valuation._approximation
[ Gauss valuation induced by (x - 4)-adic valuation, v(y - 1/4*x - 1) = 2 ]
```

Non-classical valuations
------------------------
Using the low-level interface we are not limited to classical valuations on
function fields that correspond to points on the corresponding curves. Instead
we can start with a non-trivial valuation on the field of constants.
```
sage: v = pAdicValuation(QQ, 2)
sage: R.<x> = QQ[]
sage: w = GaussValuation(R, v) # v is not trivial
sage: K.<x> = FunctionField(QQ)
sage: w = w.extension(K)
sage: w.residue_field()
Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)
```

Mac Lane Approximants
=====================
The main tool underlying this package is an algorithm by Mac Lane to compute,
starting from a Gauss valuation on a polynomial ring and a monic squarefree
polynomial G, approximations to the limit valuation which send G to infinity.
```
sage: v = pAdicValuation(QQ, 2)
sage: R.<x> = QQ[]
sage: f = x^5 + 3*x^4 + 5*x^3 + 8*x^2 + 6*x + 12
sage: v.mac_lane_approximants(f)
[[ Gauss valuation induced by 2-adic valuation, v(x^2 + x + 1) = 3 ],
 [ Gauss valuation induced by 2-adic valuation, v(x) = 1/2 ],
 [ Gauss valuation induced by 2-adic valuation, v(x) = 1 ]]
```

Note that from these approximants one can already see the residual degrees and
ramification indices of the corresponding extensions. The approximants can be
pushed to arbitrary precision.
```
sage: sage: v.mac_lane_approximants(f, required_precision=10)
[[ Gauss valuation induced by 2-adic valuation, v(x^2 + 193*x + 13/21) = 10 ],
 [ Gauss valuation induced by 2-adic valuation, v(x + 86) = 10 ],
 [ Gauss valuation induced by 2-adic valuation, v(x) = 1/2, v(x^2 + 36/11*x + 2/17) = 11 ]]
```

Note that in the limit they are factors of `f`.

References
==========

[1] Mac Lane, S. (1936). A construction for prime ideals as absolute values of
an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

[2] MacLane, S. (1936). A construction for absolute values in polynomial rings.
Transactions of the American Mathematical Society, 40(3), 363-395.

[3] Rüth, J. (2014). Models of Curves and Valuations (PhD thesis) Chapter 4
"Mac Lane Valuations".
