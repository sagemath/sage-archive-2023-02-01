
####################################################################################
## Code to make binary quadratic forms (for use in computing Siegel modular forms)
## Written by Jonathan Hanke on 8-1-2006 for Coding Sprint #1.
##     Appended to add the reduced_representatives, dyadic_trace, is_reduced, and + on 8-3-2006 for Coding Sprint #2.
##     Added Documentation and complex_point() method on 8-8-2006.
##
## Sent to: davidgru@maths.usyd.edu.au and nathan@cedar.math.ucla.edu and wstein@gmail.com
####################################################################################

from sage.rings.all import (is_fundamental_discriminant, ZZ, divisors)


class BinaryQF:

    ## Initializes the form with a 3-element list
    def __init__(self, abc_triple):
        """
        Creates a binary quadratic form ax^2 + bxy + cy^2 from the triple [a,b,c] over IntegerRing().

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q
            x^2  + 2xy  + 3y^2
        """
        assert len(abc_triple) == 3  ## Check we have three coefficients
        assert (abc_triple[0] in ZZ) and (abc_triple[1] in ZZ) and (abc_triple[2] in ZZ)  ## Check we have integer coefficients
        self.a = ZZ(abc_triple[0])
        self.b = ZZ(abc_triple[1])
        self.c = ZZ(abc_triple[2])


    ## Returns the sum of the two forms (componentwise)
    def __add__(self, Q):
        return BinaryQF([self.a + Q.a, self.b + Q.b, self.c + Q.c])


    ## Displays the output
    def __repr__(self):

        ## Deal with the zero form
        if self.a==0 and self.b==0 and self.c==0:
            return "0"

        ## Account for leading coeff
        lc_flag = True

        ## Print the first coefficient
        out_str = ""
        if abs(self.a) > 1:
            out_str +=  self.a
        if self.a != 0:
            out_str += "x^2 "
            lc_flag = False

        ## Print the second coefficient
        if self.b < 0:
            out_str += " - "
        elif self.b > 0 and lc_flag == False:
            out_str += " + "
        if abs(self.b) > 1:
            out_str += str(abs(self.b))
        if self.b != 0:
            out_str += "xy "
            lc_flag = False

        ## Print the third coefficient
        if self.c < 0:
            out_str += " - "
        elif self.c > 0 and lc_flag == False:
            out_str += " + "
        if abs(self.c) > 1:
            out_str += str(abs(self.c))
        if self.c != 0:
            out_str += "y^2 "
            lc_flag = False

        return out_str



    ## Returns the associated homogeneous quadratic polynomial
    def polynomial(self):
        """
        Returns the binary quadratic form as a 2-variable polynomial.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.polynomial()
            3*y^2 + 2*x*y + x^2
        """
        M = ZZ['x,y']
        (x,y) = M.gens()
        return self.a * x**2  + self.b* x*y + self.c * y**2

    ## Gives the dyadic trace
    def dyadic_trace(self):
        """
        Returns the "dyadic trace" of a*c - b^2 of the binary form ax^2 + bxy + cy^2

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.dyadic_trace()
            2
        """
        return self.a + self.c - self.b    ## Hopefully this is correct?!? =)


    ## Gives the discriminant
    def discriminant(self):
        """
        Returns the discriminant b^2 - 4ac of the binary form ax^2 + bxy + cy^2

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.discriminant()
            -8
        """
        return self.b**2 - 4 * self.a * self.c

    ## Checks if the discriminant D is fundamental
    def has_fundamental_discriminant(self):
        """
        Checks if the discriminant D of this form is a fundamantal discriminant
        (i.e. D is the smallest element of its squareclass with  D = 0 or 1 mod 4).

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


    ## Checks if the quadratic form is weakly reduced
    def is_weakly_reduced(self):
        """
        Checks if the form ax^2 + bxy + cy^2  satisfies |b| <= a <= c.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.is_weakly_reduced()
            True

            sage: Q = BinaryQF([2,1,3])
            sage: Q.is_weakly_reduced()
            False

            sage: Q = BinaryQF([1,-1,1])
            sage: Q.is_weakly_reduced()
            True
        """
        assert self.discriminant() < 0  ## Only consider positive definite forms for now...
        return (self.a <= abs(self.b)) and (abs(self.b) <= self.c)


    ## Checks if the quadratic form is reduced
    def is_reduced(self):
        """
        Checks if the form ax^2 + bxy + cy^2  satisfies |b| <= a <= c,
        and that b>=0  if  either a = b or a = c.

        EXAMPLES:
            sage: Q = BinaryQF([1,2,3])
            sage: Q.is_reduced()
            True

            sage: Q = BinaryQF([2,1,3])
            sage: Q.is_reduced()
            False

            sage: Q = BinaryQF([1,-1,1])
            sage: Q.is_reduced()
            False

            sage: Q = BinaryQF([1,1,1])
            sage: Q.is_reduced()
            True
        """
        assert self.discriminant() < 0  ## Only consider positive definite forms for now...

        ## Check that the form is weakly reduced
        if self.is_weakly_reduced() == False:
            return False

        ## Check that b >= 0 when a = |b| or c = |b|, since it's weakly reduced
        if self.b >= 0:
            return True
        else:
            return (self.a != abs(self.b)) and (self.c != abs(self.b))


    ## Returns the associated point in the complex upper-halfplane.
    def complex_point(self):
        """
        Returns the point in the complex upper half-plane associated to this (positive definite) quadratic form).

        EXAMPLES:
            sage: Q = BinaryQF([1,0,1])
            sage: Q.complex_point()
             1.0000000000000000*I
        """
        assert self.discriminant() < 0
        R = ZZ['x']
        x = R.gen()
        Q1 = R(self.polynomial()(x,1))
        return [z  for z in Q1.complex_roots()  if z.imag() > 0][0]



## Find (inequivalent) reduced representatives for all classes of discriminant D
def BinaryQF_reduced_representatives(D):
    """
    Returns reduced representatives for the equivalence classes of
    positive definite bianry forms of discriminant D.

    Note: These representatives are not necessarily primitive, unless
    the discriminant is fundamental!

    EXAMPLES:
        sage: BinaryQF_reduced_representatives(-4)
        [x^2  + y^2 ]

        sage: BinaryQF_reduced_representatives(-163)
        [x^2  + xy  + 41y^2 ]

        sage: BinaryQF_reduced_representatives(-12)
        [x^2  + 3y^2 , 2x^2  + 2xy  + 2y^2 ]

        sage: BinaryQF_reduced_representatives(-16)
        [x^2  + 4y^2 , 2x^2  + 2y^2 ]
    """
    assert D < 0  and  D in ZZ  and ((D % ZZ(4) == 0) or (D % ZZ(4) == 1))  ## Check that our discriminant is allowed (and is for positive definite forms)

    ## Find the range of allowed b's
    bmax = (-D / ZZ(3)).sqrt().ceil()
    b_range = range(-bmax, bmax+1)

    ## Find the set of all (possibly mutually equivalent) quadratic forms
    form_list = []
    for b in b_range:
        b_plus = abs(b)
        tmp_num = b**2 - D
        if (b**2 - D) % ZZ(4) == 0:
            tmp_num_4 = tmp_num / ZZ(4)
            b_divs__a = [a  for a in divisors(tmp_num_4) if a <= tmp_num_4.sqrt() and b_plus <= a]   ## Look for |b| <= a <= c
            for a in b_divs__a:
                c = tmp_num_4 / a
                if (a != b_plus and c != b_plus) or (b >= 0):         ## Require b>0 if a = c
                    form_list.append(BinaryQF([a,b,c]))

    ## Return the list
    return form_list
