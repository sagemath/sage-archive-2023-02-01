

from sage.quadratic_forms.count_local_2 import CountAllLocalTypesNaive

from sage.rings.arith import valuation





def VecIncrement__deprecated(self, v, R):
    """
    Performs an in-place imcrement operation on the vector v, whose
    entries are satisfy 0 <= v[i] <= R-1.  No values are returned.

void Matrix_mpz::Increment(valarray<mpz_class> & v, const mpz_class & R) const
00008 {
00009     size_t i;
00010     i = v.size();
00011
00012     // Do the carry operations
00013     while ((i > 0) && (v[i-1] == R-1))   // Assumes that all components satisfy 0 <= v[i] <= R-1
00014       {
00015         v[i-1] = 0;
00016         i--;
00017       }
00018
00019     // Only increment if we're not already at the zero vector =)
00020     if (i > 0)
00021       v[i-1]++;
00022 }
    """
    i = len(v)

    ## Do the carry operations
    while ((i > 0) and (v[i-1] == R-1)):   ## Assumes that all components satisfy 0 <= v[i] <= R-1
        v[i-1] = 0
        i += -1

    ## Only increment if we're not already at the zero vector =)
    if (i > 0):
        v[i-1] += 1




def local_solution_type__deprecated(self, p, w, zvec, nzvec):
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
        return 0


    ## Check if the solution has the appropriate (local) type


    ## 1: Check Good-type
    for i in range(len(w)):
        if (((w[i] % p) != 0)  and ((self[i,i] % p) != 0)):
            return 1
    if (p == 2):
        for i in range(len(w) - 1):
            if (((self[i,i+1] % p) != 0) and (((w[i] % p) != 0) or ((w[i+1] % p) != 0))):
                return 1


    ## 2: Check Zero-type
    Zero_flag = True
    for i in range(len(w)):
        if ((w[i] % p) != 0):
            Zero_flag = False
    if (Zero_flag == True):
        return 2


    ## Check if wS1 is zero or not
    wS1_nonzero_flag = False
    for i in range(self.dim()):

        ## Compute the valuation of each index, allowing for off-diagonal terms
        if (self[i,i] == 0):
            if (i == 0):
                val = valuation(self[i,i+1], p)    ## Look at the term to the right
            elif (i == self.dim() - 1):
                val = valuation(self[i-1,i], p)    ## Look at the term above
            else:
                val = valuation(self[i,i+1] + self[i-1,i], p)    ## Finds the valuation of the off-diagonal term since only one isn't zero
        else:
            val = valuation(self[i,i], p)

        ## Test each index
        if ((val == 1) and ((w[i] % p) != 0)):
            wS1_nonzero_flag = True


    ## 4: Check Bad-type I
    if (wS1_nonzero_flag == True):
        #print " Bad I Soln :  " + str(w)
        return 4


    ##    cout << " Bad II Soln :  " << w << "  wS1_nonzero_flag = " << wS1_nonzero_flag << endl;

    ## 5: Check Bad-type II
    if (wS1_nonzero_flag == False):
        #print " Bad II Soln :  " + str(w)
        return 5


    ## Error if we get here! =o
    print "   Solution vector is " + str(w)
    print "   and Q is \n" + str(self) + "\n"
    raise RuntimeError, "Error in IsLocalSolutionType: Should not execute this line... =( \n"






def CountAllLocalTypesNaive__deprecated(self, p, k, m, zvec, nzvec):
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


    if (True or (p == 2)):          ##  <----------------- THIS IS WIERD...

        n = self.dim()
        R = p ** k

        ## Initialize the counting vector
        count_vector = [0  for i in range(6)]

        ## Initialize v = (0, ... , 0)
        v = [0  for i in range(n)]


        ## Some declarations to speed up the loop
        R_n = R ** n
        m1 = m % R

        ## Count the local solutions
        #for(size_t i=1; i<=(R_n).get_ui(); i++):
        for i in range(R_n):
#            self.VecIncrement(v, R)                          ## Increments the vector v in-place.
            VecIncrement_cdef(v, R)                          ## Increments the vector v in-place.
            if (self(v) % R  == m1):                        ## Evaluates the quadratic form (mod R) at the vector v -- can be sped up!
                solntype = self.local_solution_type(p, v, zvec, nzvec)
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


        return count_vector

    else:
        raise RuntimeError, "Error in count_local_typeNaive: Matrix \n" + str(self) + "\n is not symmetric!"





##///////////////////////////////////////////
##/// Front-ends for our counting routines //
##///////////////////////////////////////////

def count_local_type(self, p, k, m, solntype, zvec, nzvec):
    """
    mpz_class Matrix_mpz::count_local_type(const mpz_class & p, long k, const mpz_class & m, size_t solntype,
                         const valarray<size_t> & zero, const valarray<size_t> & nonzero) const

    // Ideally this would use the count_local_typeWithSymmetry routine, but this is fine for now. =)
    """
    #return self.CountAllLocalTypesNaive(p, k, m, zvec, nzvec)[0]
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[0]


def count_local_good_type(self, p, k, m, zvec, nzvec):
    """
    mpz_class Matrix_mpz::count_local_good_type(const mpz_class & p, long k, const mpz_class & m,
                             const valarray<size_t> & zero, const valarray<size_t> & nonzero) const

    """
    #return self.CountAllLocalTypesNaive(p, k, m, zvec, nzvec)[1]
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[1]


def count_local_zero_type(self, p, k, m, zvec, nzvec):
    """
    mpz_class Matrix_mpz::count_local_zero_type(const mpz_class & p, long k, const mpz_class & m,
                             const valarray<size_t> & zero, const valarray<size_t> & nonzero) const

    """
    #return self.CountAllLocalTypesNaive(p, k, m, zvec, nzvec)[2]
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[2]


def count_local_bad_type(self, p, k, m, zvec, nzvec):
    """
    mpz_class Matrix_mpz::count_local_bad_type(const mpz_class & p, long k, const mpz_class & m,
                             const valarray<size_t> & zero, const valarray<size_t> & nonzero) const

    """
    #return self.CountAllLocalTypesNaive(p, k, m, zvec, nzvec)[3]
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[3]


def count_local_bad_typeI(self, p, k, m, zvec, nzvec):
    """
    mpz_class Matrix_mpz::count_local_bad_typeI(const mpz_class & p, long k, const mpz_class & m,
                             const valarray<size_t> & zero, const valarray<size_t> & nonzero) const

    """
    #return self.CountAllLocalTypesNaive(p, k, m, zvec, nzvec)[4]
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[4]


def count_local_bad_typeII(self, p, k, m, zvec, nzvec):
    """
    mpz_class Matrix_mpz::count_local_bad_typeII(const mpz_class & p, long k, const mpz_class & m,
                             const valarray<size_t> & zero, const valarray<size_t> & nonzero) const

    """
    #return self.CountAllLocalTypesNaive(p, k, m, zvec, nzvec)[5]
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[5]
