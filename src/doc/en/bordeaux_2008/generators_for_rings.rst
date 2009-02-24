Generators for Rings of Modular Forms
=====================================

Computing Generators
--------------------

For any congruence subgroup :math:`\Gamma`, the direct sum

.. math::

   M(\Gamma) =  \bigoplus_{k\geq 0} M_k(\Gamma)

is a ring, since the product of modular forms
:math:`f\in M_k(\Gamma)` and :math:`g \in M_{k'}(\Gamma)` is
an element :math:`fg \in M_{k+k'}(\Gamma)`. Sage can compute
likely generators for rings of modular forms, but currently doesn't
prove any of these results.

We verify the statement proved in Serre's "A Course in Arithmetic"
that :math:`E_4` and :math:`E_6` generate the space of level
one modular forms.

.. skip

::

    sage: from sage.modular.modform.find_generators import modform_generators
    sage: modform_generators(1)
    [(4, 1 + 240*q + 2160*q^2 + 6720*q^3 + O(q^4)),
     (6, 1 - 504*q - 16632*q^2 - 122976*q^3 + O(q^4))]

Have you ever wondered which forms generate the ring
:math:`M(\Gamma_0(2))`? it turns out a form of weight 2 and two
forms of weight 4 together generate.

.. skip

::

    sage: modform_generators(2)
    [(2, 1 + 24*q + 24*q^2 + ... + 288*q^11 + O(q^12)),
     (4, 1 + 240*q^2 + .. + 30240*q^10 + O(q^12)),
     (4, q + 8*q^2 + .. + 1332*q^11 + O(q^12))]

Here's generators for :math:`M(\Gamma_0(3))`. Notice that
elements of weight :math:`6` are now required, in addition to
weights :math:`2` and :math:`4`.

.. skip

::

    sage: modform_generators(3)
    [(2, 1 + 12*q + 36*q^2 + .. + 168*q^13 + O(q^14)),
     (4, 1 + 240*q^3 + 2160*q^6 + 6720*q^9 + 17520*q^12 + O(q^14)),
     (4, q + 9*q^2 + 27*q^3 + 73*q^4 + .. + O(q^14)),
     (6, q - 6*q^2 + 9*q^3 + 4*q^4 + .. + O(q^14)),
     (6, 1 - 504*q^3 - 16632*q^6 .. + O(q^14)),
     (6, q + 33*q^2 + 243*q^3 + .. + O(q^14))]
