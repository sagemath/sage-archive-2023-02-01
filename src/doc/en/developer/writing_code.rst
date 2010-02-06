Writing Code for Sage
=====================

If there is something you would like to implement and make available
in Sage, you have a wide range of options:

#. Implement it as Sage scripts.

#. Implement it as Python scripts that use the Sage library.

#. Implement it in C/C++ and make the result accessible to Sage using
   Cython.

#. Implement it using Cython.

#. Implement it using one or more of the following: Flint, FpLLL, GAP,
   GSL, IML, LinBox, M4RI, Matplotlib, Maxima, MWRank, ECLib,
   NetworkX, NTL, Numpy, PARI/GP, PolyBoRi, R, Scipy, Singular, Sympy
   or any of the other libraries included with Sage [1]_.

#. Or any combination of the above.

If you have Magma, Maple or Mathematica and do not mind restricting
who can use your code, you could also implement parts of your
program in one of these systems and make it available in Sage.

Flint, FpLLL, GAP, GSL, IML, LinBox, M4RI, Matplotlib, Maxima, MWRank,
ECLib, NetworkX, NTL, Numpy, PARI/GP, PolyBoRi, R, Scipy, Singular,
Sympy are all included with all distributions of Sage. GAP, Singular,
and PARI are very mature, and each implements a great amount of
functionality, though in different domains. GAP addresses group theory
well, Singular attacks polynomial computation, and PARI contains
sophisticated, optimized number theory algorithms. Notably absent from
this triad is a good system for exact linear algebra (something Magma
does extremely well), but this gap is being filled by code written for
Sage or covered by specialized C/C++ libraries like LinBox, IML and
M4RI.

.. Say something about GSL, Matplotlib, Maxima, MWRANK, NetworkX, NTL, Numpy

Sage is not just about gathering together functionality. It is about
providing a clear, systematic and consistent way to access a large
number of algorithms, in a coherent framework that makes sense
mathematically. In the design of Sage, the semantics of objects, the
definitions, etc., are informed by how the corresponding objects are
used in everyday mathematics.

This document was authored by William Stein, David Joyner, John
Palmieri and others with the editorial help of Iftikhar Burhanuddin
and Martin Albrecht.

.. [1]
   See http://www.sagemath.org/links-components.html for a full list
   of packages shipped with every copy of Sage

Contents:

.. toctree::

   conventions
   coding_in_python
   coding_in_other
   doctesting
   sage_manuals
