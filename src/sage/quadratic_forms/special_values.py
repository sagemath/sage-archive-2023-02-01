##########################################################
## Routines for computing special values of L-functions ##
##########################################################

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.arith import kronecker_symbol, bernoulli, factorial, fundamental_discriminant
from sage.rings.all import RealField
from sage.combinat.combinat import bernoulli_polynomial
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.symbolic.constants import pi, I
from sage.rings.real_mpfr import is_RealField
from sage.misc.functional import denominator, numerator
from sage.rings.infinity import infinity


## ---------------- The Gamma Function  ------------------

def gamma__exact(n):
    """
    Evaluates the exact value of the gamma function at an integer or
    half-integer argument.

    EXAMPLES:
        sage: gamma__exact(4)
        6
        sage: gamma__exact(3)
        2
        sage: gamma__exact(2)
        1
        sage: gamma__exact(1)
        1

        sage: gamma__exact(1/2)
        sqrt(pi)
        sage: gamma__exact(3/2)
        1/2*sqrt(pi)
        sage: gamma__exact(5/2)
        3/4*sqrt(pi)
        sage: gamma__exact(7/2)
        15/8*sqrt(pi)

        sage: gamma__exact(-1/2)
        -2*sqrt(pi)
        sage: gamma__exact(-3/2)
        4/3*sqrt(pi)
        sage: gamma__exact(-5/2)
        -8/15*sqrt(pi)
        sage: gamma__exact(-7/2)
        16/105*sqrt(pi)

    """
    from sage.all import sqrt
    ## SANITY CHECK
    if (not n in QQ) or (denominator(n) > 2):
        raise TypeError, "Oops!  You much give an integer or half-integer argument."

    if (denominator(n) == 1):
        if n <= 0:
            return infinity
        if n > 0:
            return factorial(n-1)
    else:
        ans = QQ(1)
        while (n != QQ(1)/2):
            if (n<0):
                ans *= QQ(1)/n
                n = n + 1
            elif (n>0):
                n = n - 1
                ans *= n

        ans *= sqrt(pi)
        return ans



## ------------- The Riemann Zeta Function  --------------

def zeta__exact(n):
    """
    Returns the exact value fo the Riemann Zeta function

    References:
        Iwasawa's "Lectures on p-adic L-functions", p13
            Special value of zeta(2k)
        Ireland and Rosen's "A Classical Introduction to Modern Number Theory"
        Washington's "Cyclotomic Fields"

    EXAMPLES:
        sage: ## Testing the accuracy of the negative special values
        sage: RR = RealField(100)
        sage: for i in range(1,10):
        ...       print "zeta(" + str(1-2*i) + "): ", RR(zeta__exact(1-2*i)) - zeta(RR(1-2*i))
        zeta(-1):  0.00000000000000000000000000000
        zeta(-3):  0.00000000000000000000000000000
        zeta(-5):  0.00000000000000000000000000000
        zeta(-7):  0.00000000000000000000000000000
        zeta(-9):  0.00000000000000000000000000000
        zeta(-11):  0.00000000000000000000000000000
        zeta(-13):  0.00000000000000000000000000000
        zeta(-15):  0.00000000000000000000000000000
        zeta(-17):  0.00000000000000000000000000000

        sage: RR = RealField(100)
        sage: for i in range(1,10):
        ...       print "zeta(" + str(1-2*i) + "): ", RR(zeta__exact(1-2*i)) - zeta(RR(1-2*i))
        zeta(-1):  0.00000000000000000000000000000
        zeta(-3):  0.00000000000000000000000000000
        zeta(-5):  0.00000000000000000000000000000
        zeta(-7):  0.00000000000000000000000000000
        zeta(-9):  0.00000000000000000000000000000
        zeta(-11):  0.00000000000000000000000000000
        zeta(-13):  0.00000000000000000000000000000
        zeta(-15):  0.00000000000000000000000000000
        zeta(-17):  0.00000000000000000000000000000

    """
    if n<0:
        k = 1-n
        return -bernoulli(k)/k
    elif n>1:
        if (n % 2 == 0):
            return ZZ(-1)**(n/2 + 1) * ZZ(2)**(n-1) * pi**n * bernoulli(n) / factorial(n)
        else:
            raise TypeError, "n must be a critical value! (I.e. even > 0 or odd < 0.)"
    elif n==1:
        return infinity
    elif n==0:
        return -1/2

## ---------- Dirichlet L-functions with quadratic characters ----------

def QuadraticBernoulliNumber(k, d):
    """
    Compute k-th Bernoulli number for the primitive
    quadratic character associated to chi(x) = (d/x).

    References:  Iwasawa's "Lectures on p-adic L-functions", pp7-16.

    EXAMPLES:
        sage: ## Makes a set of odd fund discriminants < -3
        sage: Fund_odd_test_set = [D  for D in range(-163, -3, 4)  if is_fundamental_discriminant(D)]

        sage: ## In general, we have B_{1, \chi_d} = -2h/w  for odd fund disc < 0
        sage: for D in Fund_odd_test_set:
        ...      if len(BinaryQF_reduced_representatives(D)) != -QuadraticBernoulliNumber(1, D):
        ...          print "Oops!  There is an error at D = ", D



    """
    ## Ensure the character is primitive
    d1 = fundamental_discriminant(d)
    f = abs(d1)

    ## Make the (usual) k-th Bernoulli polynomial
    x =  PolynomialRing(QQ, 'x').gen()
    bp = bernoulli_polynomial(x, k)

    ## Make the k-th quadratic Bernoulli number
    total = sum([kronecker_symbol(d1, i) * bp(i/f)  for i in range(f)])
    total *= (f ** (k-1))

    return total


def quadratic_L_function__exact(n, d):
    """
    Returns the exact value of a quadratic twist of the Riemann Zeta function by chi_d(x) = (d/x).

    References:
        Iwasawa's "Lectures on p-adic L-functions", p16-17
            Special values of L(1-n, chi) and L(n, chi)
        Ireland and Rosen's "A Classical Introduction to Modern Number Theory"
        Washington's "Cyclotomic Fields"

    EXAMPLES:
        sage: bool(quadratic_L_function__exact(1, -4) == pi/4)
        True

    """
    from sage.all import SR, sqrt
    if n<=0:
        k = 1-n
        return -QuadraticBernoulliNumber(k,d)/k
    elif n>=1:
        ## Compute the kind of critical values (p10)
        if kronecker_symbol(fundamental_discriminant(d), -1) == 1:
            delta = 0
        else:
            delta = 1

        ## Compute the positive special values (p17)
        if ((n - delta) % 2 == 0):
            f = abs(fundamental_discriminant(d))
            if delta == 0:
                GS = sqrt(f)
            else:
                GS = I * sqrt(f)
            ans = SR(ZZ(-1)**(1+(n-delta)/2))
            ans *= (2*pi/f)**n
            ans *= GS     ## Evaluate the Gauss sum here! =0
            ans *= 1/(2 * I**delta)
            ans *= QuadraticBernoulliNumber(n,d)/factorial(n)
            return ans
        else:
            if delta == 0:
                raise TypeError, "n must be a critical value!\n" + "(I.e. even > 0 or odd < 0.)"
            if delta == 1:
                raise TypeError, "n must be a critical value!\n" + "(I.e. odd > 0 or even <= 0.)"



def quadratic_L_function__numerical(n, d, num_terms=1000):
    """
    Evaluate the Dirichlet L-function (for quadratic character) numerically
    (in a very naive way).

    EXAMPLES:
        sage:  ## Test several values for a given character
        sage: RR = RealField(100)
        sage: for i in range(5):
        ...       print "L(" + str(1+2*i) + ", (-4/.)): ", RR(quadratic_L_function__exact(1+2*i, -4)) - quadratic_L_function__numerical(RR(1+2*i),-4, 10000)
        L(1, (-4/.)):  0.000049999999500000024999996962707
        L(3, (-4/.)):  4.9999997000000374953360932931e-13
        L(5, (-4/.)):  4.9999992362712983899678889809e-21
        L(7, (-4/.)):  2.9187853493177436800233956915e-29
        L(9, (-4/.)):  0.00000000000000000000000000000

        sage: ## Testing the accuracy of the negative special values
        sage: ## ---- THIS FAILS SINCE THE DIRICHLET SERIES DOESN'T CONVERGE HERE! ----

        sage: ## Test several characters agree with the exact value, to a given accuracy.
        sage: for d in range(-20,0):
        ...       if abs(RR(quadratic_L_function__numerical(1, d, 10000) - quadratic_L_function__exact(1, d))) > 0.001:
        ...           print "Oops!  We have a problem at d = ", d, "    exact = ", RR(quadratic_L_function__exact(1, d)), "    numerical = ", RR(quadratic_L_function__numerical(1, d))
        ...

    """
    ## Set the correct precision if it's given (for n).
    if is_RealField(n.parent()):
        R = n.parent()
    else:
        R = RealField()

    d1 = fundamental_discriminant(d)
    ans = R(0)
    for i in range(1,num_terms):
        ans += R(kronecker_symbol(d1,i) / R(i)**n)
    return ans
