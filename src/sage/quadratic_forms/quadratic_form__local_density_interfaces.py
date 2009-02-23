
## // This is needed in the filter for primitivity...
## #include "../max-min.h"

from copy import deepcopy

from sage.rings.arith import valuation
from sage.rings.rational_field import QQ, RationalField


#  ////////////////////////////////
#  // Private Front-end Routines //
#/////////////////////////////////////////////////////////////////////////////////////////////

def local_good_density(self, p, m):
    """
    Finds the Good-type local density of Q representing m at p.
    (Front end routine for its congruence counterpart.)
    """
    #print "Doing Good Density with p = " + str(p) +  " and m = " + str(m)
    return self.local_good_density_congruence(p, m, [], [])


def local_zero_density(self, p, m):
    """
    Finds the Zero-type local density of Q representing m at p.
    (Front end routine for its congruence counterpart.)
    """
    #print "Doing Good Density with p = " + str(p) +  " and m = " + str(m)
    return self.local_zero_density_congruence(p, m, [], [])


def local_bad_density(self, p, m):
    """
    Finds the Bad-type local density of Q representing m at p.
    (Front end routine for its congruence counterpart.)
    """
    #print "Doing Bad Density with p = " + str(p) +  " and m = " + str(m)
    return self.local_bad_density_congruence(p, m, [], [])


def local_badI_density(self, p, m):
    """
    Finds the Bad-type I local density of Q representing m at p.
    (Front end routine for its congruence counterpart.)
    """
    #print "Doing Bad I Density with p = " + str(p) +  " and m = " + str(m)
    return self.local_badI_density_congruence(p, m, [], [])


def local_badII_density(self, p, m):
    """
    Finds the Bad-type II local density of Q representing m at p.
    (Front end routine for its congruence counterpart.)
    """
    #print "Doing Bad II Density with p = " + str(p) +  " and m = " + str(m)
    return self.local_badII_density_congruence(p, m, [], [])


## ---------------  These are the important ones, which we'll filter for primitive forms!!!  ------------------


### TODO:  THESE TWO ROUTINES HAVE **A LOT** OF CODE IN COMMON, BUT NEITHER PUTS THE FORM IN LOCAL NORMAL FORM FIRST!


def local_density(self, p, m):
    """
    Gives the local density -- should be called by the user. =)

    NOTE: This screens for imprimitive forms, but *doesn't* put the
    quadratic form in local normal form, which is a *requirement* of
    the routines performing the computations!

    mpq_class Matrix_mpz::local_density(const mpz_class & p, const mpz_class & m) const {

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])   ## NOTE: This is already in local normal form for *all* primes p!
        sage: Q.local_density(p=2, m=1)
        1
        sage: Q.local_density(p=3, m=1)
        8/9
        sage: Q.local_density(p=5, m=1)
        24/25
        sage: Q.local_density(p=7, m=1)
        48/49
        sage: Q.local_density(p=11, m=1)
        120/121

    """
    n = self.dim()
    if (n == 0):
        raise TypeError, "Oops!  We currently don't handle 0-dim'l forms. =("


    ## Set the modulus to check for imprimitive forms at p
    no_val_flag = True

    ## Check for imprimitive forms at p -- ASSUMES THE FORM IS UPPER TRIANGULAR!
    ## (NOTE: We could do better if we know it's normalized!)
    for i in range(n):
        for j in range(i, n):
            if (self[i,j] != 0):
                if (no_val_flag == True):
                    no_val_flag = False
                    p_valuation = valuation(self[i,j], p)
                else:
                    p_valuation = min(p_valuation, valuation(self[i,j], p))


    ## DIAGNOSTIC
    #cout << " Using the matrix: \n " << (*this) << endl;
    #cout << "Valuation(m,p) = " << Valuation(m,p) << endl;
    #cout << "p_valuation = " << p_valuation << endl;


    ## If m is less p-divisible than the matrix, return zero
    if ((m != 0) and (valuation(m,p) < p_valuation)):   ## Note: The (m != 0) condition protects taking the valuation of zero.
        return QQ(0)

    ## If the form is imprimitive, divide it (and m) by p-powers to get a primitive form
    else:
        if (p_valuation > 0):

            ## Make a new (primitive) matrix
            Q1 = deepcopy(self)
            #Q1 = QuadraticForm(self.base_ring(), self.dim())

            ## DIAGNOSTIC
            #print " p = " << p
            #print " p_valuation = " << p_valuation
            #print " p_mod = " << (p ** p_valuation)

            p_mod = p ** p_valuation      ## This should give a power...
            for i in range(n):
                for j in range(i, n):
                    Q1[i,j] = self[i,j] / p_mod

            ## Make a new number mm
            mm = m / p_mod

            ## Then return the densities for the reduced problem
            return Q1.local_good_density(p, mm) + Q1.local_zero_density(p, mm) + Q1.local_bad_density(p, mm)


        ## Otherwise, proceed as usual... =)
        else:
            return self.local_good_density(p, m) + self.local_zero_density(p, m) + self.local_bad_density(p, m)




def local_primitive_density(self, p, m):
    """
    Gives the local primitive density -- should be called by the user. =)

    NOTE: This screens for imprimitive forms, but *doesn't* put the
    quadratic form in local normal form, which is a *requirement* of
    the routines performing the computations!

    mpq_class Matrix_mpz::local_density(const mpz_class & p, const mpz_class & m) const {
    """

    n = self.dim()
    if (n == 0):
        raise TypeError, "Oops!  We currently don't handle 0-dim'l forms. =("


    ## Set the modulus to check for imprimitive forms at p
    no_val_flag = True

    ## Check for imprimitive forms at p -- ASSUMES THE FORM IS UPPER TRIANGULAR!
    ## (NOTE: We could do better if we know it's normalized!)
    for i in range(n):
        for j in range(i, n):
            if (self[i,j] != 0):
                if (no_val_flag == True):
                    no_val_flag = False
                    p_valuation = valuation(self[i,j], p)
                else:
                    p_valuation = min(p_valuation, valuation(self[i,j], p))


    ## DIAGNOSTIC
    #cout << " Using the matrix: \n " << (*this) << endl;
    #cout << "Valuation(m,p) = " << Valuation(m,p) << endl;
    #cout << "p_valuation = " << p_valuation << endl;


    ## If m is less p-divisible than the matrix, return zero
    if ((m != 0) and (valuation(m,p) < p_valuation)):   ## Note: The (m != 0) condition protects taking the valuation of zero.
        return QQ(0)

    ## If the form is imprimitive, divide it (and m) by p-powers to get a primitive form
    else:
        if (p_valuation > 0):

            ## Make a new (primitive) matrix
            Q1 = deepcopy(self)
            #Q1 = QuadraticForm(self.base_ring(), self.dim())

            ## DIAGNOSTIC
            #print " p = " << p
            #print " p_valuation = " << p_valuation
            #print " p_mod = " << (p ** p_valuation)

            p_mod = p ** p_valuation      ## This should give a power...
            for i in range(n):
                for j in range(i, n):
                    Q1[i,j] = self[i,j] / p_mod

            ## Make a new number mm
            mm = m / p_mod

            ## Then return the densities for the reduced problem
            return Q1.local_good_density(p, mm) + Q1.local_bad_density(p, mm)


        ## Otherwise, proceed as usual... =)
        else:
            return self.local_good_density(p, m) + self.local_bad_density(p, m)

