r"""
Yamanouchi Words

A right (respectively left) Yamanouchi word on a completely ordered
alphabet, for instance [1,2,...,n], is a word math such that any
right (respectively left) factor of math contains more entries math
than math. For example, the word [2, 3, 2, 2, 1, 3, 1, 2, 1, 1] is
a right Yamanouchi one.

The evaluation of a word math encodes the number of occurrences of
each letter of math. In the case of Yamanouchi words, the
evaluation is a partition. For example, the word [2, 3, 2, 2, 1, 3,
1, 2, 1, 1] has evaluation [4, 4, 2].

Yamanouchi words can be useful in the computation of
Littlewood-Richardson coefficients `c_{\lambda, \mu}^\nu`.
According to the Littlewood-Richardson
rule, `c_{\lambda, \mu}^\nu` is the number of skew tableaux
of shape `\nu / \lambda` and evaluation `\mu`,
whose row readings are Yamanouchi words.
"""

