
include "../ext/cdefs.pxi"
include "../ext/gmp.pxi"

from sage.rings.arith import valuation
from sage.rings.integer_mod import IntegerMod, Mod
from sage.rings.integer_mod_ring import IntegerModRing

from sage.rings.integer_ring import ZZ

from sage.rings.integer_mod cimport IntegerMod_gmp





cdef CountAllLocalTypesNaive_cdef(Q, p, k, m, zvec, nzvec):
    """
    ///////////////////////////////////////////////////////////////////
    /// Naively counts the number of solutions of Q(x) = m mod p^k   //
    /// of type solntype, satisfying the mod p congruence conditions //
    /// at the indices of the vectors "zero" and "nonzero"           //
    ///////////////////////////////////////////////////////////////////

    valarray <mpz_class> Matrix_mpz::CountAllLocalTypesNaive(const mpz_class & p, unsigned long k, const mpz_class & m,
                                             const valarray<size_t> & zero, const valarray<size_t> & nonzero) const
    """

    ## DIAGNOSTIC
    #print "   --> CountAllLocalTypesNaive is using the form Q \n" + str(Q)
    #print "       p = " + str(p) + "  and   m = " + str(m)

    #cdef mpz_t* v
    cdef long n, i
    cdef long a, b    ## Used to quickly evaluate Q(v)
    cdef long ptr     ## Used to increment the vector
    cdef long solntype    ## Used to store the kind of solution we find




    n = Q.dim()
    R = p ** k

    ## Cython Variables
    cdef IntegerMod_gmp zero, one
    zero = IntegerMod_gmp(IntegerModRing(R), 0)
    one = IntegerMod_gmp(IntegerModRing(R), 1)


    Q1 = Q.base_change_to(IntegerModRing(R))



    ###########################################
    #v = <mpz_t*> sage_malloc(sizeof(mpz_t)*n)
    #for i from 0 <= i < n:
    #    mpz_init(v[i])
    ###########################################

    ## Initialize the counting vector
    count_vector = [0  for i in range(6)]

    ## Initialize v = (0, ... , 0)
    v = [Mod(0, R)  for i in range(n)]
    #v = []
    #for i in 1 <= i < n:
    #    v.append(zero)


    ## Some declarations to speed up the loop
    R_n = R ** n
    m1 = Mod(m, R)

    ## Count the local solutions
    for i from 0 <= i < R_n:

        ## Perform a carry (when value = R-1) until we can increment freely
        ptr = len(v)
        while ((ptr > 0) and (v[ptr-1] == R-1)):
            v[ptr-1] += 1
            ptr += -1

        ## Only increment if we're not already at the zero vector =)
        if (ptr > 0):
            v[ptr-1] += 1


        ## Evaluate Q(v) quickly
        tmp_val = Mod(0, R)
        for a from 0 <= a < n:
            for b from a <= b < n:
                tmp_val += Q1[a,b] * v[a] * v[b]

        ## Sort the solution by it's type
        #if (Q1(v) == m1):
        if (tmp_val == m1):
            solntype = local_solution_type_cdef(Q1, p, v, zvec, nzvec)
            if (solntype != 0):
                count_vector[solntype] += 1



    ## Generate the Bad-type and Total counts
    count_vector[3] = count_vector[4] + count_vector[5]
    count_vector[0] = count_vector[1] + count_vector[2] + count_vector[3]


    ## DIAGNOSTIC
    #cout << "R = " << R << "\n";
    #cout << "n = " << n << "\n";
    #cout << "R_n = " << R_n << "\n";
    #
    #for(i=1; i<=25; i++) {
    #    cout << "v = " << v << "\n";
    #    Increment(v,R);
    #}
    #
    #cout << "Q1 = " << Q1 << "\n";
    #cout << "v = " << v << "\n";
    #cout << "Q1 * v = " << Q1 * v<< "\n";

    ###########################
    #for i from 0 <= i < n:
    #    mpz_clear(v[i])
    #sage_free(v)
    ###########################

    return count_vector



def CountAllLocalTypesNaive(Q, p, k, m, zvec, nzvec):
    """
    Python wrapper function for this (now) cython function.
    """
    return CountAllLocalTypesNaive_cdef(Q, p, k, m, zvec, nzvec)







cdef local_solution_type_cdef(Q, p, w, zvec, nzvec):
    """
    ////////////////////////////////////////////////////////////////////////////////////
    /// Private routine to check if a given solution vector w (of Q(w) = m mod p^k)   //
    /// is of a certain local type and satisfies certain congruence conditions mod p. //
    ///   (Personal Note: For p=2, we should still use p=2 and not p=8.)              //
    ////////////////////////////////////////////////////////////////////////////////////

    size_t Matrix_mpz::local_solution_type(const mpz_class & p, const valarray<mpz_class> & w,
               const valarray<size_t> & zero, const valarray<size_t> & nonzero) const

    """

    ## Note: Here p is assumed to be a prime >= 2, though the routine still works if not...

    ## ToDo?: Add a check that Q is square and has the same size as w.

    cdef long i
    cdef long n

    n = Q.dim()

    zero_flag = False        ## Tests the zero mod p congruence conditions
    nonzero_flag = False     ## Tests the nonzero congruence conditions


    ## Check if the solution satisfies the zvec "zero" congruence conditions
    ## (either zvec is empty or its components index the zero vector mod p)
    if (len(zvec) == 0):
        zero_flag = True
    else:
        i = 0
        while ( (i < len(zvec)) and ((w[zvec[i]] % p) == 0) ):
            i += 1
        if (i == len(zvec)):
            zero_flag = True


    ## DIAGNOSTIC
    #print "IsLocalSolutionType: Finished the Zero congruence condition test \n"

    if (zero_flag == False):
        return 0

    ## DIAGNOSTIC
    #print "IsLocalSolutionType: Passed the Zero congruence condition test \n"


    ## Check if the solution satisfies the nzvec "nonzero" congruence conditions
    ## (either nzvec is empty or its components index a non-zero vector mod p)
    if (len(nzvec) == 0):
        nonzero_flag = True
    else:
        i = 0
        while ((nonzero_flag == False) and (i < len(nzvec))):
            if ((w[nzvec[i]] % p) != 0):
                nonzero_flag = True
            i += 1

    if (nonzero_flag == False):
        return <long> 0


    ## Check if the solution has the appropriate (local) type


    ## 1: Check Good-type
    for i from 0 <= i < n:
        if (((w[i] % p) != 0)  and ((Q[i,i] % p) != 0)):
            return <long> 1
    if (p == 2):
        for i from 0 <= i < (n - 1):
            if (((Q[i,i+1] % p) != 0) and (((w[i] % p) != 0) or ((w[i+1] % p) != 0))):
                return <long> 1


    ## 2: Check Zero-type
    Zero_flag = True
    for i from 0 <= i < n:
        if ((w[i] % p) != 0):
            Zero_flag = False
    if (Zero_flag == True):
        return <long> 2


    ## Check if wS1 is zero or not
    wS1_nonzero_flag = False
    for i from 0 <= i < n:

        ## Compute the valuation of each index, allowing for off-diagonal terms
        if (Q[i,i] == 0):
            if (i == 0):
                val = valuation(Q[i,i+1], p)    ## Look at the term to the right
            elif (i == n - 1):
                val = valuation(Q[i-1,i], p)    ## Look at the term above
            else:
                val = valuation(Q[i,i+1] + Q[i-1,i], p)    ## Finds the valuation of the off-diagonal term since only one isn't zero
        else:
            val = valuation(Q[i,i], p)

        ## Test each index
        if ((val == 1) and ((w[i] % p) != 0)):
            wS1_nonzero_flag = True


    ## 4: Check Bad-type I
    if (wS1_nonzero_flag == True):
        #print " Bad I Soln :  " + str(w)
        return <long> 4


    ##    cout << " Bad II Soln :  " << w << "  wS1_nonzero_flag = " << wS1_nonzero_flag << endl;

    ## 5: Check Bad-type II
    if (wS1_nonzero_flag == False):
        #print " Bad II Soln :  " + str(w)
        return <long> 5


    ## Error if we get here! =o
    print "   Solution vector is " + str(w)
    print "   and Q is \n" + str(Q) + "\n"
    raise RuntimeError, "Error in IsLocalSolutionType: Should not execute this line... =( \n"

