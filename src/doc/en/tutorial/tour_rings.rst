***********
Basic Rings
***********

We illustrate some basic rings in Sage. For example, the field of
rational numbers may be referred to using either ``RationalField()``
or ``QQ``:

::

    sage: RationalField()
    Rational Field
    sage: QQ
    Rational Field
    sage: 1/2 in QQ
    True

The decimal number ``1.2`` is considered to be in ``QQ``, since there is
a coercion map from the reals to the rationals:

::

    sage: 1.2 in QQ
    True

However, there is no coercion map from the finite field with 3
elements to the rationals:

::

    sage: c = GF(3)(1)   # c is the element 1 of the field GF(3)
    sage: c in QQ
    False

Also, of course, the symbolic constant :math:`\pi` is not in the rationals:

::

    sage: pi in QQ
    False

The symbol ``I`` represents the square root of :math:`-1`; ``i`` is a synonym for
``I``. Of course, this is not in the rationals:

::

    sage: i  # square root of -1
    I
    sage: i in QQ
    False

By the way, some other pre-defined Sage rings are the integers
``ZZ``, the real numbers ``RR``, and the complex numbers ``CC``. We
discuss polynomial rings in :ref:`section-poly`.

Now we illustrate some arithmetic.

::

    sage: a, b = 4/3, 2/3
    sage: a + b
    2
    sage: 2*b == a
    True
    sage: parent(2/3)
    Rational Field
    sage: parent(4/2)
    Rational Field
    sage: 2/3 + 0.1       # automatic coercion before addition
    0.766666666666667
    sage: 0.1 + 2/3       # coercion rules are symmetric in SAGE
    0.766666666666667

There is one subtlety in defining complex numbers: as mentioned
above, the symbol ``i`` represents a square root of :math:`-1`, but it is a
*formal* square root of :math:`-1`.  Calling
``CC(i)`` returns the complex square root of :math:`-1`.

.. link

::

    sage: i = CC(i)       # floating point complex number
    sage: z = a + b*i
    sage: z
    1.33333333333333 + 0.666666666666667*I
    sage: z.imag()        # imaginary part
    0.666666666666667
    sage: z.real() == a   # automatic coercion before comparison
    True
    sage: QQ(11.1)
    111/10