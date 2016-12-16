********
Appendix
********

.. _section-precedence:

Arithmetical binary operator precedence
=======================================

What is ``3^2*4 + 2%5``? The value (38) is determined by this
"operator precedence table". The table below is based on the table
in § 5.14 of the *Python Language Reference Manual* by G. Rossum
and F. Drake. the operations are listed here in increasing order of
precedence.


==========================  =================
Operators                   Description
==========================  =================
or                          boolean or
and                         boolean and
not                         boolean not
in, not in                  membership
is, is not                  identity test
>, <=, >, >=, ==, !=        comparison
+, -                        addition, subtraction
\*, /, %                    multiplication, division, remainder
\*\*, ^                     exponentiation
==========================  =================

Therefore, to compute ``3^2*4 + 2%5``, Sage brackets the
computation this way: ``((3^2)*4) + (2%5)``. Thus, first compute
``3^2``, which is ``9``, then compute both ``(3^2)*4`` and ``2%5``,
and finally add these.
