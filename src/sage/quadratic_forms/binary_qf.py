"""
Binary Quadratic Forms with integer coefficients.

The form $a x^2 + b x y + d y^2$ is stored as a triple (a, b, c) of integers.

EXAMPLES:
    sage: Q = BinaryQF([1,2,3])
    sage: Q
    x^2 + 2*x*y  + 3*y^2
    sage: Q.discriminant()
    -8
    sage: Q.reduced_form()
    x^2 + 2*y^2
    sage: Q(1, 1)
    6

TESTS:
    sage: Q == loads(dumps(Q))
    True

AUTHORS:
    -- Jon Hanke (2006-08-08):
        * Appended to add the reduced_representatives, dyadic_trace,
          is_reduced, and + on 8-3-2006 for Coding Sprint #2.
        * Added Documentation and complex_point() method on 8-8-2006.
    -- Nick Alexander: add doctests and clean code for Doc Days 2
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein and Jon Hanke
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.pari.all import pari
from sage.rings.all import (is_fundamental_discriminant, ZZ, divisors)
from sage.structure.sage_object import SageObject
class BinaryQF(SageObject):
    """
    BinaryQF([a,b,c])

    INPUT:
        v -- a list of 3 entries:  [a,b,c]

    OUTPUT:
        the binary quadratic form a*x^2 + b*x*y + c*y^2.

    EXAMPLES:
        sage: b = BinaryQF([1,2,3])
        sage: b.discriminant()
        -8
    """
    # Initializes the form with a 3-element list
    def __init__(self, abc_triple):
        """
        Creates the binary quadratic form $ax^2 + bxy + cy^2$ from the
        triple [a,b,c] over IntegerRing().

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q
            x^2 + 2*x*y  + 3*y^2
            sage: Q = BinaryQF([1,2])
            Traceback (most recent call last):
            ...
            ValueError: Binary quadratic form must be given by a list of three coefficients
        """
        if len(abc_triple) != 3:
            # Check we have three coefficients
            raise ValueError, "Binary quadratic form must be given by a list of three coefficients"
        self.a = ZZ(abc_triple[0])
        self.b = ZZ(abc_triple[1])
        self.c = ZZ(abc_triple[2])

    def __getitem__(self, n):
        """
        Return the n-th component of this quadratic form.

        If this form is $a x^2 + b x y + c y^2$, the 0-th component is $a$,
        the 1-st component is $b$, and $2$-nd component is $c$.

        Indexing is like lists -- negative indices and slices are allowed.

        EXAMPLES:
            sage: Q = BinaryQF([2,3,4])
            sage: Q[0]
            2
            sage: Q[2]
            4
            sage: Q[:2]
            (2, 3)
            sage: tuple(Q)
            (2, 3, 4)
            sage: list(Q)
            [2, 3, 4]
        """
        return (self.a, self.b, self.c)[n]

    def __call__(self, *args):
        r"""
        Evaluate this quadratic form at a point.

        INPUT:
            args -- x and y values, often as a pair x, y or a list [x, y]

        EXAMPLES:
            sage: Q = BinaryQF([2,3,4])
            sage: Q(1, 2)
            24
        """
        return self.polynomial()(*args)

    def __cmp__(self, right):
        """
        Returns True if self and right are identical: the same coefficients.

        EXAMPLES:
            sage: P = BinaryQF([2,2,3])
            sage: Q = BinaryQF([2,2,3])
            sage: R = BinaryQF([1,2,3])
            sage: P == Q # indirect doctest
            True
            sage: P == R # indirect doctest
            False

        TESTS:
            sage: P == P
            True
            sage: Q == P
            True
            sage: R == P
            False
            sage: P == 2
            False
        """
        if not isinstance(right, BinaryQF):
            return cmp(type(self), type(right))
        return cmp((self.a,self.b,self.c), (right.a,right.b,right.c))

    def __add__(self, Q):
        """
        Returns the component-wise sum of two forms.

        That is, given $a_1 x^2 + b_1 x y + c_1 y^2$ and $a_2 x^2 + b_2 x y +
        c_2 y^2$, returns the form
        $$(a_1 + a_2) x^2 + (b_1 + b_2) x y + (c_1 + c_2) y^2 .$$

        EXAMPLES:
            sage: P = BinaryQF([2,2,3]); P
            2*x^2 + 2*x*y + 3*y^2
            sage: Q = BinaryQF([-1,2,2]); Q
            -x^2 + 2*x*y + 2*y^2
            sage: P + Q
            x^2 + 4*x*y + 5*y^2
            sage: P + Q == BinaryQF([1,4,5]) # indirect doctest
            True

        TESTS:
            sage: Q + P == BinaryQF([1,4,5]) # indirect doctest
            True
        """
        return BinaryQF([self.a + Q.a, self.b + Q.b, self.c + Q.c])

    def __sub__(self, Q):
        """
        Returns the component-wise difference of two forms.

        That is, given $a_1 x^2 + b_1 x y + c_1 y^2$ and $a_2 x^2 + b_2 x y +
        c_2 y^2$, returns the form
        $$(a_1 - a_2) x^2 + (b_1 - b_2) x y + (c_1 - c_2) y^2 .$$

        EXAMPLES:
            sage: P = BinaryQF([2,2,3]); P
            2*x^2 + 2*x*y + 3*y^2
            sage: Q = BinaryQF([-1,2,2]); Q
            -x^2 + 2*x*y + 2*y^2
            sage: P - Q
            3*x^2 + y^2
            sage: P - Q == BinaryQF([3,0,1]) # indirect doctest
            True

        TESTS:
            sage: Q - P == BinaryQF([3,0,1]) # indirect doctest
            False
            sage: Q - P != BinaryQF([3,0,1]) # indirect doctest
            True
        """
        return BinaryQF([self.a - Q.a, self.b - Q.b, self.c - Q.c])

    def _repr_(self):
        """
        Display the quadratic form.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3]); Q # indirect doctest
            x^2 + 2*x*y + 3*y^2

            sage: Q = BinaryQF([-1,2,3]); Q
            -x^2 + 2*x*y + 3*y^2

            sage: Q = BinaryQF([0,0,0]); Q
            0
        """
        return repr(self.polynomial())

    def polynomial(self):
        """
        Returns the binary quadratic form as a homogeneous 2-variable
        polynomial.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.polynomial()
            x^2 + 2*x*y + 3*y^2

            sage: Q = BinaryQF([-1,-2,3])
            sage: Q.polynomial()
            -x^2 - 2*x*y + 3*y^2

            sage: Q = BinaryQF([0,0,0])
            sage: Q.polynomial()
            0
        """
        M = ZZ['x,y']
        (x,y) = M.gens()
        return self.a * x**2  + self.b* x*y + self.c * y**2

    def discriminant(self):
        """
        Returns the discriminant $b^2 - 4ac$ of the binary
        form $ax^2 + bxy + cy^2$.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.discriminant()
            -8
        """
        return self.b**2 - 4 * self.a * self.c

    def has_fundamental_discriminant(self):
        """
        Checks if the discriminant D of this form is a fundamantal
        discriminant (i.e. D is the smallest element of its
        squareclass with D = 0 or 1 mod 4).

        EXAMPLES:
            sage: Q = BinaryQF([1,0,1])
            sage: Q.discriminant()
            -4
            sage: Q.has_fundamental_discriminant()
            True

            sage: Q = BinaryQF([2,0,2])
            sage: Q.discriminant()
            -16
            sage: Q.has_fundamental_discriminant()
            False
        """
        return is_fundamental_discriminant(self.discriminant())

    def is_weakly_reduced(self):
        """
        Checks if the form $ax^2 + bxy + cy^2$  satisfies
        $|b| \leq a \leq c$, i.e., is weakly reduced.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.is_weakly_reduced()
            False

            sage: Q = BinaryQF([2,1,3])
            sage: Q.is_weakly_reduced()
            True

            sage: Q = BinaryQF([1,-1,1])
            sage: Q.is_weakly_reduced()
            True
        """
        if self.discriminant() >= 0:
            raise NotImplementedError, "only implemented for nagative discriminants"
        return (abs(self.b) <= self.a) and (self.a <= self.c)

    def reduced_form(self):
        """
        EXAMPLES:
            sage: a = BinaryQF([33,11,5])
            sage: a.is_reduced()
            False
            sage: b = a.reduced_form(); b
            5*x^2 - x*y + 27*y^2
            sage: b.is_reduced()
            True

            sage: a = BinaryQF([15,0,15])
            sage: a.is_reduced()
            True
            sage: b = a.reduced_form(); b
            15*x^2 + 15*y^2
            sage: b.is_reduced()
            True
        """
        if self.discriminant() >= 0 or self.a < 0:
            raise NotImplementedError, "only implemented for positive definite forms"
        if not self.is_reduced():
            v = list(pari('Vec(qfbred(Qfb(%s,%s,%s)))'%(self.a,self.b,self.c)))
            return BinaryQF(v)
        else:
            return self

    def reduce(self):
        """
        EXAMPLES:
            sage: a = BinaryQF([37,17,2])
            sage: a.is_reduced()
            False
            sage: a.reduce(); a
            x^2 + x*y + 2*y^2
            sage: a.is_reduced()
            True
        """
        if not self.is_reduced():
            red = self.reduced_form()
            self.a, self.b, self.c = red.a, red.b, red.c

    def is_reduced(self):
        """
        Checks if the quadratic form is reduced, i.e., if the form
        $ax^2 + bxy + cy^2$ satisfies $|b|\leq a \leq c$, and
        that $b\geq 0$ if either $a = b$ or $a = c$.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.is_reduced()
            False

            sage: Q = BinaryQF([2,1,3])
            sage: Q.is_reduced()
            True

            sage: Q = BinaryQF([1,-1,1])
            sage: Q.is_reduced()
            False

            sage: Q = BinaryQF([1,1,1])
            sage: Q.is_reduced()
            True
        """
        return (-self.a < self.b <= self.a < self.c) or \
               (ZZ(0) <= self.b <= self.a == self.c)

    def complex_point(self):
        """
        Returns the point in the complex upper half-plane associated
        to this (positive definite) quadratic form).

        For positive definite forms with negative discriminants, this is a
        root $\tau$ of $a x^2 + b x + c$ with the imaginary part of $\tau$
        greater than 0.

        EXAMPLES:
            sage: Q = BinaryQF([1,0,1])
            sage: Q.complex_point()
            1.00000000000000*I
        """
        if self.discriminant() >= 0:
            raise NotImplementedError, "only implemented for negative discriminant"
        R = ZZ['x']
        x = R.gen()
        Q1 = R(self.polynomial()(x,1))
        return [z  for z in Q1.complex_roots()  if z.imag() > 0][0]


def BinaryQF_reduced_representatives(D):
    """
    Returns a list of inequivalent reduced representatives for the equivalence
    classes of positive definite binary forms of discriminant D.

    NOTE: The order of the representatives is unspecified but deterministic.

    WARNING: The representatives are not necessarily primitive, unless the
    discriminant is fundamental!

    EXAMPLES:
        sage: BinaryQF_reduced_representatives(-4)
        [x^2 + y^2]

        sage: BinaryQF_reduced_representatives(-163)
        [x^2 + x*y + 41*y^2]

        sage: BinaryQF_reduced_representatives(-12)
        [x^2 + 3*y^2, 2*x^2 + 2*x*y + 2*y^2]

        sage: BinaryQF_reduced_representatives(-16)
        [x^2 + 4*y^2, 2*x^2 + 2*y^2]

        sage: BinaryQF_reduced_representatives(-63)
        [x^2 + x*y + 16*y^2, 2*x^2 - x*y + 8*y^2, 2*x^2 + x*y + 8*y^2, 3*x^2 + 3*x*y + 6*y^2, 4*x^2 + x*y + 4*y^2]

        The number of inequivalent reduced binary forms with a fixed negative
        fundamental discriminant D is the class number of the quadratic field
        $Q(\sqrt{D})$.

        sage: len(BinaryQF_reduced_representatives(-13*4))
        2
        sage: QuadraticField(-13*4, 'a').class_number()
        2
        sage: p=next_prime(2^20); p
        1048583
        sage: len(BinaryQF_reduced_representatives(-p))
        689
        sage: QuadraticField(-p, 'a').class_number()
        689
    """
    D = ZZ(D)
    if not ( D < 0 and (D % 4 in [0,1])):
        raise ValueError, "discriminant is not valid and positive definite"

    form_list = []

    from sage.misc.all import xsrange
    for a in xsrange(1,1+(-3*D).isqrt()):
        a4 = 4*a
        #Only iterate over b of the same parity as D
        for b in xsrange(D%2,a+1,2):
            t = b*b-D
            if t % a4 == 0:
                c = t // a4
                if a < c:
                    form_list.append(BinaryQF([a,b,c]))
                    if b>0 and a>b:
                        form_list.append(BinaryQF([a,-b,c]))
                else:
                    if (a == c):
                        form_list.append(BinaryQF([a,b,c]))

    form_list.sort()
    return form_list
