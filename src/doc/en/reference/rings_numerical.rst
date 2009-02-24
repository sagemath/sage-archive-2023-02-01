Fixed and Arbitrary Precision Numerical Fields
==============================================

Sage supports two optimized fixed precision fields for the numerical
computation, the real double (RealDoubleField) and complex double
fields (ComplexDoubleField).
Sage also supports arbitrary precision real (RealField) and
complex fields (ComplexField).  Sage also supports optimized
interval arithmetic (RealIntervalField).


Real and complex double elements are optimized implementations that
use the GNU Scientific Library for arithmetic and some special
functions.  Arbitrary precision real and complex numbers are
implemented using the MPFR library, which builds on GMP.  (Note that
Sage doesn't currently use the MPC library.)  The interval arithmetic
field is implemented using the MPFI library.

In many cases the PARI C-library is used to compute special functions
when implementations aren't otherwise available.

.. toctree::
   :maxdepth: 2

   sage/rings/real_double
   sage/rings/complex_double
   sage/rings/real_mpfr
   sage/rings/complex_field
   sage/rings/complex_number
   sage/rings/real_mpfi