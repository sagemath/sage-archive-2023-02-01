One initial typo causes many failures::

    sage: a = binomiak(10,5)
    sage: a
    252
    sage: a == factorial(10)/factorial(5)^2
    True
    sage: a//12
    21

But this is a new failure::

    sage: binomial(10,5)
    255
