
## Include headers from the front-end (repeated execution) routines.
#include <assert.h>



#####################################################
## Historical Note:
## ----------------
## (Much) Older versions of the NTL code can be found at
## Old-Laptop-Odyssues_Extended-Home-directory_8-28-2003/home/jonhanke/TEMP/C++/Modular_Project/
##
##  Updated for GMP libraries on 1/21/04
#####################################################

from copy import deepcopy

from sage.sets.set import Set
from sage.rings.rational_field import QQ
from sage.rings.arith import valuation, kronecker_symbol
from sage.rings.integer_ring import ZZ
from sage.misc.misc import prod, verbose

##############################################################################
### Creates a new index vector which points to the extracted rows/columns
### in the extracted matrix where the ...?
###
### This is used internally in our reduction procedure which allows us
### to keep track of a list of (row/column) indices through the
### process of extracting a submatrix corresponding to some other set
### of (row/column) indices.
###
### (Note: This is probably not very efficient,
### and could be improved by using vectors,
### but we're using valarrays for now...)
#############################################################################


def reindex_vector_from_extraction(self, Original, Extracted):
    """
    This takes a vector of indices, and applies a  truncated permutation to them,
    defined by Extracted.

    TO DO: Please revisit this routine, and eliminate it!

    valarray<size_t> Matrix_mpz::reindex_vector_from_extraction(const valarray<size_t> & Original, const valarray<size_t> & Extracted) const
    """
    ## Set some default values (all to be overwritten)
    Reindex = [-1  for i in range(len(Original))]

    ## Replace all Original indices entries with the position of the matching Extraced index
    ind = 0
    for i in range(len(Original)):
        for i in range(len(Extracted)):
            if (Original[i] == Extracted[j]):
                Reindex[ind] = j
                ind += 1

    ## Copy these to a new vecarray of the appropriate length -- Since Reindex may be too big
    Final = [Reindex[i]  for i in range(ind)]

    return Final




def count_modp__by_gauss_sum(self, n, p, m, Qdet):
    """
    Returns the number of solutions of Q(x) = m over the finite field F_p,
    where p is a prime number > 2 and Q has determinant Qdet.

    These are defined in Table 1 on p363 of Hanke's "Local Densities..." paper.
    """
    neg1 = -1

    ## DIAGNOSTIC
    verbose(" n = " + str(n))
    verbose(" neg1 = " + str(neg1))
    verbose(" p**(n-1) = " + str(p**(n-1)))

    if ((p != 2) and (n >= 1)):   ## TO DO: Check that p is an odd prime and n >= 1  --  To Do: Need to check if p is a prime...

        if (m % p == 0):
            if (n % 2 != 0):
                count = (p**(n-1))
            else:
                count = (p**(n-1)) + (p-1) * (p**((n-2)/2)) * kronecker_symbol(((neg1**(n/2)) * Qdet) % p, p)
        else:
            if (n % 2 != 0):
                count = (p**(n-1)) + (p**((n-1)/2)) * kronecker_symbol(((neg1**((n-1)/2)) * Qdet * m) % p, p)
            else:
                count = (p**(n-1)) - (p**((n-2)/2)) * kronecker_symbol(((neg1**(n/2)) * Qdet) % p, p)

    else:
        raise RuntimeError, "\n Error in count_modp__by_gauss_sum: Either p is not prime or n<1. \n"

    return count





def local_good_density_congruence_odd(self, p, m, Zvec, NZvec):
    """
    Finds the Good-type local density of Q representing m at p.
    (Assuming that p > 2 and Q is given in local diagonal form.)

    mpq_class Matrix_mpz::local_good_density_congruence_odd(const mpz_class & p, const mpz_class & m,
                                            const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const
    """

    ## To Do: Check to see if Q_int is diagonal (only ok for p>2)


    n = self.dim()


    ## Assuming Q is diagonal, find the indices of the p-unit (diagonal) entries
    Qtrimvec = [i  for i in range(n)  if (self.__getitem__((i,i)) % p) != 0]



    ## DEBUGGING:  Check the matrix is primitive --
    ##   (WARNING: We may be passed imprimitive form from the BI-reduction...)  <==  This should be fixed now. =)
    if (Qtrimvec == []):
        raise RuntimeError, "Uh oh.  There should be some p-unit diagonal elements at this stage...  but Qtrimvec is empty. =("



    ## DIAGNOSTIC
    #print " Stage 1"
    #print " Q = " + str(Q)
    #print " n = " + str(n) + "  and  Qtrimvec = " + str(Qtrimvec)


    Qtrim = self.extract_variables(Qtrimvec)


    ## print "We're here now..." << endl;


    ## Construct the big and small matrices:
    ## -------------------------------------

    ## Construct new congruence condition indices for the trimmed matrix
    trimZvec = self.reindex_vector_from_extraction(Zvec, Qtrimvec)


    ## DIAGNOSTIC
    #print "\n Zvec = " + str(Zvec)
    #print " NZvec = " + str(NZvec)
    #print " Qtrim = " + str(Qtrim)
    #print " Qtrimvec = " + str(Qtrimvec)
    #print " trimZvec = " + str(trimZvec)
    #####  print " Vector indexing Qtrim = " << range(len(Qtrimvec))
    #####  print " and its complement by trimZvec  = " << VectorComplement(range(len(Qtrimvec)), trimZvec)



    ## print "Getting closer..."


    ## Make the big trim vector (all trim indices not somewhere in Zvec)
    ##   and the big free vector (all non-trim indices not somewhere in Zvec)
    ##   [To Do: This could probably be faster if we assumed Zvec was ordered.]
    bigtrimvec = [ i  for i in range(len(Qtrimvec))  if i not in trimZvec]
    bigfreelength = (n - len(Qtrimvec)) - (len(Zvec) - len(trimZvec))




    ## Make the small vector from the big one (all indices not somewhere in Zvec or NZvec)
    new_vec = list(Set(Zvec + NZvec))
    trim_new_vec = self.reindex_vector_from_extraction(new_vec, Qtrimvec)


    smalltrimvec = [ i  for i in range(len(Qtrimvec))  if i not in new_vec]
    smallfreelength = (n - len(Qtrimvec)) - (len(new_vec) - len(trim_new_vec))



    ## Make the big and small matrices
    bigmatrix = Qtrim.extract_variables(bigtrimvec)
    smallmatrix = Qtrim.extract_variables(smalltrimvec)


    ## DIAGNOSTIC
    #print "\n Q is : " + str(self)
    #print " m is : " + str(m)
    #print " p is : " + str(p)
    #print " Qtrim is : " + str(Qtrim)
    #print " bigtrimvec is : " + str(bigtrimvec)
    #print " bigfreelength is : " + str(bigfreelength)
    #print " smalltrimvec is : " + str(smalltrimvec)
    #print " smallfreelength is : " + str(smallfreelength)
    #print " bigmatrix is : \n" + str(bigmatrix)
    #print " smallmatrix is : \n" + str(smallmatrix)


    ## print "And closer still..."


    ## Number of representations


    ## Adjust the local solutions to count only the good-type ones
    if (bigtrimvec == []):    ## Check if we have the empty matrix...
        big_factor = 0

    else:
        big_dim = bigmatrix.dim()
        big_det = prod([bigmatrix[i,i]  for i in range(bigmatrix.dim())])   ## This gives the product of the diagonal entries, not the Gram matrix det

        if (m % p != 0):
            big_factor = (p**bigfreelength) \
                * self.count_modp__by_gauss_sum(big_dim, p, m, big_det)
        else:
            big_factor = (p**bigfreelength) \
                * (self.count_modp__by_gauss_sum(big_dim, p, m, big_det) - 1)


    ## Similarly for the smallmatrix if it exists
    if (smalltrimvec == []):    ## Check if we have the empty matrix...
        small_factor = 0

    else:
        if (len(NZvec) > 0):
            small_dim = smallmatrix.dim()
            small_det = prod([smallmatrix[i,i]  for i in range(smallmatrix.dim())])   ## This gives the product of the diagonal entries, not the Gram matrix det

            if (m % p != 0):
                small_factor = (p**smallfreelength) \
                    * self.count_modp__by_gauss_sum(small_dim, p, m, small_det)
            else:
                small_factor = (p**smallfreelength) \
                    * (self.count_modp__by_gauss_sum(small_dim, p, m, small_det) - 1)

        else:
            small_factor = 0

    total = big_factor - small_factor
    good_density = QQ(total / p**(n-1))


    ## DIAGNOSTIC
    #print "\n big factor = " + str(big_factor)
    #print " small factor = " + str(small_factor)
    #print " big_det = " + str(prod([bigmatrix[i,i]  for i in range(bigmatrix.dim())]))
    #print " smal_det = " + str(prod([smallmatrix[i,i]  for i in range(smallmatrix.dim())]))
    #print " total = " + str(total)
    #print " denominator = " + str(p**(n-1))
    #print " Good Density = " + str(good_density)
    #print " count_modp__by_gauss_sum output = " + str(self.count_modp__by_gauss_sum(big_dim, p, m, prod([bigmatrix[i,i]  for i in range(bigmatrix.dim())])))

    return good_density





def local_good_density_congruence_even(self, p, m, Zvec, NZvec):
    """
    Finds the Good-type local density of Q representing m at p.
    (Assuming that p = 2 and Q is given in local normal form.)

    mpq_class Matrix_mpz::local_good_density_congruence_even(const mpz_class & p, const mpz_class & m,
                                             const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const
    """
    #print " Break 0"
    #print "\n Q is : " + str(Q)

    n = self.dim()


    ## Assuming Q is diagonal, trim it to ensure it's non-degenerate mod 8      <--- ??? Do we really assume that Q is diagonal?!?!?  =(
    Qtrimvec = []

    #print " Break 0.1"

    ## Find the indices of the non-zero blocks mod 8
    for i in range(n):

        ## DIAGNOSTIC
        verbose(" i = " + str(i))
        verbose(" n = " + str(n))
        verbose(" Qtrimvec = " + str(Qtrimvec))

        nz_flag = False

        if  ((self[i,i] % 8) != 0):
            nz_flag = True
        else:
            #print " here 1"
            if ((i == 0) and ((self[i,i+1] % 8) != 0)):
                nz_flag = True
            else:
                #print " here 2" << endl;
                if ((i == n-1) and ((self[i-1,i] % 8) != 0)):
                    nz_flag = True
                else:
                    #print " here 3" << endl;
                    if ( (i > 0)  and  (i < n-1)  and  (((self[i,i+1] % 8) != 0) or ((self[i-1,i] % 8) != 0)) ):
                        nz_flag = True


        if (nz_flag == True):
            Qtrimvec += [i]


    #print " Break 1"


    ## DEBUGGING: Quick tests to make sure the form isn't zero mod 8
    if (Qtrimvec == []):
        raise RuntimeError, "Oh no!  The quadratic form should always be non-degenerate mod 8... =("


    Qtrim = self.extract_variables(Qtrimvec)


    ## DEBUGGING: Quick tests to make sure the extracted form isn't zero mod 8
    assert(Qtrim.dim() > 0);

    #print " Break 2"


    ## Construct the big and small matrices:
    ## -------------------------------------

    ## Construct new congruence condition indices for the trimmed matrix
    trimZvec = self.reindex_vector_from_extraction(Zvec, Qtrimvec)
    trimNZvec = self.reindex_vector_from_extraction(NZvec, Qtrimvec)

    ## Make the trimmed congruence vector
    new_vec = list(Set(Zvec + NZvec))
    trim_vec = self.reindex_vector_from_extraction(new_vec, Qtrimvec)

    ## DIAGNOSTIC
    verbose("")
    verbose("Zvec = " + str(Zvec))
    verbose("NZvec = " + str(NZvec))
    verbose("new_vec = " + str(new_vec))
    verbose("")


    ## DIAGNOSTIC
    verbose("\n Q is : " + str(self))
    verbose(" m is : " + str(m))
    verbose(" Qtrim is : " + str(Qtrim))
    verbose(" Qtrimvec is : " + str(Qtrimvec))
    verbose(" trimZvec is : " + str(trimZvec))
    verbose(" trimNZvec is : " + str(trimNZvec))


    ## DEBUGGING: Check that partlyfreenum is in range...
    if not (len(new_vec) >= len(trim_vec)):
        raise RuntimeError, "Oh no!  Some variable called 'partlyfreenum' is out of range. =("


    ## Compute the number of different free components
    partlyfreenum = len(new_vec) - len(trim_vec)
    veryfreenum = (n - len(Qtrimvec)) - partlyfreenum


    ## In the free part, each component with a congruence condition contrubite a factor of 4,
    ## while components with no congruence conditions contribute a factor of 8.

    total = (4 ** partlyfreenum) * (8 ** veryfreenum) \
        * Qtrim.count_local_good_type(2, 3, m, trimZvec, trimNZvec)
    good_density = ZZ(total) / ZZ(8**(n-1))


    ## DIAGNOSTIC
    verbose(" partlyfreenum = " + str(partlyfreenum))
    verbose(" veryfreenum = " + str(veryfreenum))
    verbose(" Qtrim.count_local_good_type(2, 3, m, trimZvec, trimNZvec) = " + str(Qtrim.count_local_good_type(2, 3, m, trimZvec, trimNZvec)))
    verbose(" total = " + str(total))
    verbose("    total has type ", type(total))
    verbose(" denominator = " + str(8**(n-1)))
    verbose(" Good Density = " + str(good_density))

    return good_density




def local_good_density_congruence_even__experimental(self, p, m, Zvec, NZvec):
    """
    Finds the Good-type local density of Q representing m at p.
    (Assuming that p = 2 and Q is given in local normal form.)

    mpq_class Matrix_mpz::local_good_density_congruence_even(const mpz_class & p, const mpz_class & m,
                                             const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const
    """
    #print " Break 0"
    #print "\n Q is : " + str(Q)

    n = self.dim()


    ## Assuming Q is diagonal, trim it to ensure it's non-degenerate mod 8      <--- ??? Do we really assume that Q is diagonal?!?!?  =(
    Qtrimvec = []

    #print " Break 0.1"

    ## Find the indices of the non-zero blocks mod 8
    for i in range(n):

        ## DIAGNOSTIC
        verbose(" i = " + str(i))
        verbose(" n = " + str(n))
        verbose(" Qtrimvec = " + str(Qtrimvec))

        nz_flag = False

        if  ((self[i,i] % 8) != 0):
            nz_flag = True
        else:
            #print " here 1"
            if ((i == 0) and ((self[i,i+1] % 8) != 0)):
                nz_flag = True
            else:
                #print " here 2" << endl;
                if ((i == n-1) and ((self[i-1,i] % 8) != 0)):
                    nz_flag = True
                else:
                    #print " here 3" << endl;
                    if ( (i > 0)  and  (i < n-1)  and  (((self[i,i+1] % 8) != 0) or ((self[i-1,i] % 8) != 0)) ):
                        nz_flag = True


        if (nz_flag == True):
            Qtrimvec += [i]


    #print " Break 1"


    ## DEBUGGING: Quick tests to make sure the form isn't zero mod 8
    if (Qtrimvec == []):
        raise RuntimeError, "Oh no!  The quadratic form should always be non-degenerate mod 8... =("


    Qtrim = self.extract_variables(Qtrimvec)


    ## DEBUGGING: Quick tests to make sure the extracted form isn't zero mod 8
    assert(Qtrim.dim() > 0);

    #print " Break 2"


    ## Construct the big and small matrices:
    ## -------------------------------------

    ## Construct new congruence condition indices for the trimmed matrix
    trimZvec = self.reindex_vector_from_extraction(Zvec, Qtrimvec)
    trimNZvec = self.reindex_vector_from_extraction(NZvec, Qtrimvec)

    ## Make the trimmed congruence vector
    new_vec = list(Set(Zvec + NZvec))
    trim_vec = self.reindex_vector_from_extraction(new_vec, Qtrimvec)

    ## DIAGNOSTIC
    verbose("")
    verbose("Zvec = " + str(Zvec))
    verbose("NZvec = " + str(NZvec))
    verbose("new_vec = " + str(new_vec))
    verbose("")


    ## DIAGNOSTIC
    verbose("\n Q is : " + str(self))
    verbose(" m is : " + str(m))
    verbose(" Qtrim is : " + str(Qtrim))
    verbose(" Qtrimvec is : " + str(Qtrimvec))
    verbose(" trimZvec is : " + str(trimZvec))
    verbose(" trimNZvec is : " + str(trimNZvec))


    ## DEBUGGING: Check that partlyfreenum is in range...
    if not (len(new_vec) >= len(trim_vec)):
        raise RuntimeError, "Oh no!  Some variable called 'partlyfreenum' is out of range. =("


    ## Compute the number of different free components
    partlyfreenum = len(new_vec) - len(trim_vec)
    veryfreenum = (n - len(Qtrimvec)) - partlyfreenum


    ## In the free part, each component with a congruence condition contrubite a factor of 4,
    ## while components with no congruence conditions contribute a factor of 8.

    total = (4 ** partlyfreenum) * (8 ** veryfreenum) \
        * Qtrim.count_local_good_type(2, 3, m, trimZvec, trimNZvec)
    good_density = ZZ(total) / ZZ(8**(n-1))


    ## DIAGNOSTIC
    verbose(" partlyfreenum = " + str(partlyfreenum))
    verbose(" veryfreenum = " + str(veryfreenum))
    verbose(" Qtrim.count_local_good_type(2, 3, m, trimZvec, trimNZvec) = " + str(Qtrim.count_local_good_type(2, 3, m, trimZvec, trimNZvec)))
    verbose(" total = " + str(total))
    verbose("    total has type ", type(total))
    verbose(" denominator = " + str(8**(n-1)))
    verbose(" Good Density = " + str(good_density))

    return good_density



    ############################################
     ## Note: Assumes all forms are primitive       <== What does this mean?!?
    ############################################



def local_good_density_congruence(self, p, m, Zvec, NZvec):
    """
    Finds the Good-type local density of Q representing m at p.
    (Front end routine for parity specific routines for p.)

    mpq_class Matrix_mpz::local_good_density_congruence(const mpz_class & p, const mpz_class & m,
                                        const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const
    """
    ## DIAGNOSTIC
    verbose(" In local_good_density_congruence with ")
    verbose(" Q is: \n" + str(self))
    verbose(" p = " + str(p))
    verbose(" m = " + str(m))
    verbose(" Zvec = " + str(Zvec))
    verbose(" NZvec = " + str(NZvec))


    ## Check that Q is in local normal form -- should replace this with a diagonalization check?
    ##   (it often may not be since the reduction procedure
    ##   often mixes up the order of the valuations...)
    #
    #if (self != self.local_normal_form(p))
    #    print "Warning in local_good_density_congruence: Q is not in local normal form! \n";



    ## Check that the congruence conditions don't overlap
    if (len(Set(Zvec).intersection(Set(NZvec))) != 0):
        return QQ(0)

    #print "We're here!"


    ## Decide which routine to use to compute the Good-type density
    if (p > 2):
        return self.local_good_density_congruence_odd(p, m, Zvec, NZvec)

    if (p == 2):
        #print "\n Using the (p=2) Local_Good_Density_Even routine! \n"
        return self.local_good_density_congruence_even(p, m, Zvec, NZvec)

    raise RuntimeError, "\n Error in Local_Good_Density: The 'prime' p = " + str(p) + " is < 2. \n"







def local_zero_density_congruence(self, p, m, Zvec, NZvec):
    """
    Finds the Zero-type local density of Q representing m at p,
    allowing certain congruence conditions mod p.

    mpq_class Matrix_mpz::local_zero_density_congruence(const mpz_class & p, const mpz_class & m,
                                        const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const
    """
    ## DIAGNOSTIC
    #print " In local_zero_density_congruence with "
    #print " Q is: \n" + str(self)
    #print " p = " + str(p)
    #print " m = " + str(m)
    #print " Zvec = " + str(Zvec)
    #print " NZvec = " + str(NZvec)

    p2 = p * p

    if ((m % (p2) != 0) or (len(NZvec) > 0)):
        return 0
    else:
        ## Need 2 steps since we can't take negative powers... =|
        EmptyVec = []
        if (self.dim() > 2):
            return QQ(1 / p**(self.dim() - 2)) \
                * self.local_density_congruence(p, m / (p2), EmptyVec, EmptyVec);
        else:
            return QQ(p**(2 - self.dim())) \
                * self.local_density_congruence(p, m / (p2), EmptyVec, EmptyVec);






def local_badI_density_congruence(self, p, m, Zvec, NZvec):
    """
    Finds the Bad-type I local density of Q representing m at p.
    (Assuming that p > 2 and Q is given in local diagonal form.)

    mpq_class Matrix_mpz::local_badI_density_congruence(const mpz_class & p, const mpz_class & m,
                                        const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const
    """
    ## DIAGNOSTIC
    #print " In local_badI_density_congruence with "
    #print " Q is: \n" + str(self)
    #print " p = " + str(p)
    #print " m = " + str(m)
    #print " Zvec = " + str(Zvec)
    #print " NZvec = " + str(NZvec)


    n = self.dim()
    if n == 0:
        raise RuntimeException, "Oops!  We don't currently allow 0-dim'l forms... =( "


    ## Define the indexing sets S_i
    S0 = []
    S1_empty_flag = True    ## This is used to check if we should be computing BI solutions at all!
                            ## (We should really to this earlier, but S1 must be non-zero to proceed.)


    ## Find the valuation of each variable (which will be the same over 2x2 blocks),
    ## remembering those of valuation 0 and if an entry of valuation 1 exists.
    for i in range(n):

        ## Compute the valuation of each index, allowing for off-diagonal terms
        if (self[i,i] == 0):
            if (i == 0):
                val = valuation(self[i,i+1], p)    ## Look at the term to the right
            else:
                if (i == n-1):
                    val = valuation(self[i-1,i], p)    ## Look at the term above
                else:
                    val = valuation(self[i,i+1] + self[i-1,i], p)    ## Finds the valuation of the off-diagonal term since only one isn't zero
        else:
            val = valuation(self[i,i], p)

        if (val == 0):
            S0 += [i]
        elif (val == 1):
            S1_empty_flag = False    ## Need to have a non-empty S1 set to proceed with Bad-type I reduction...

    ## Check that S1 is non-empty to proceed, otherwise return no solutions.
    if (S1_empty_flag == True):
        return 0


    ## Check that the form is primitive...
    if (S0 == []):
        print " Using Q = " + str(self)
        print " and p = " + str(p)
        raise RuntimeError, "Oops! The form is not primitive!"


    S0New = S0    ## REDUNDANT ALERT!!



    ## Note: The following lines assume that NZvec is an ordered vector.  Should check this here...

    ## DIAGNOSTIC
    #print " m = " + str(m) + "   p = " + str(p)
    #print " S0 = " + str(S0) + "   NZvec = " str(NZvec) "   IsDisjoint = " + str(IsDisjointOrdered(S0, NZvec))
    #print " len(S0) = " + str(len(S0))


    ## Make the form Qnew for the reduction procedure
    if ((m % p == 0) and (len(Set(S0New).intersection(Set(NZvec))) == 0) and (len(S0New) != n)):
        #print "Test 1.1 "
        Qnew = deepcopy(self)
        #Qnew = QuadraticForm(self.base_ring(), self.dim(), self.entries)


        ## Run over all indices with i, while keeping track of the next index in S0 as S0[j].
        j = 0       ## The index for following the increasing set (of indices) S0
        for i in range(n):
            ## Note: Short circuit && is necessary here:
            if ((j < len(S0New)) and (i == S0New[j])):                 ## i is in S0
                j += 1
                Qnew[i,i] = p * Qnew[i,i]
                if ((p == 2) and (i < n)):
                    Qnew[i+1,i] = p * Qnew[i+1,i]
                    Qnew[i,i+1] = p * Qnew[i,i+1]
            else:                             ## i not in S0;  Since S0 is increasing, we know i < S0(j)
                Qnew[i,i] = Qnew[i,i] / p

                #print "  dividing in row " + str(i)

                if ((p == 2) and (i < n-1)):
                    Qnew[i+1,i] = Qnew[i+1,i] / p
                    Qnew[i,i+1] = Qnew[i,i+1] / p

        ## DIAGNOSTIC
        #print "\n\n Check of Bad-type I reduction: \n";
        #print " Q is " + str(self)
        #print " Qnew is " + str(Qnew)
        #print " p = " + str(p)
        #print " m / p = " + str(m/p)
        #print " VectorComplement(Zvec, S0) is " + str(VectorComplement(Zvec, S0))       <------  FIX THIS!
        #print " NZvec " << + str(NZvec)


        ## Do the reduction
        ## (Need 2 steps since we can't take negative powers... =| )
        VC = list(Set([i  for i in Zvec  if i not in S0New]))
        if (len(S0New) > 1):
            return QQ(1 / p**(len(S0New) - 1)) * Qnew.local_good_density_congruence(p, m / p, VC, NZvec)
        else:
            return QQ(p**(1 - len(S0New))) * Qnew.local_good_density_congruence(p, m / p, VC, NZvec)

    else:
        return 0



def local_badII_density_congruence(self, p, m, Zvec, NZvec):
    """
    Finds the Bad-type II local density of Q representing m at p.
    (Assuming that p > 2 and Q is given in local diagonal form.)

    mpq_class Matrix_mpz::local_badII_density_congruence(const mpz_class & p, const mpz_class & m,
                                        const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const
    """
    ## DIAGNOSTIC
    #print " In local_badII_density_congruence with "
    #print " Q is: \n" + str(self)
    #print " p = " + str(p)
    #print " m = " + str(m)
    #print " Zvec = " + str(Zvec)
    #print " NZvec = " + str(NZvec)


    n = self.dim()

    ## Define the indexing sets S_i
    S0 = []
    S1 = []
    S2 = []

    for i in range(n):

        ## Compute the valuation of each index, allowing for off-diagonal terms
        if (self[i,i] == 0):
            if (i == 0):
                val = valuation(self[i,i+1], p)    ## Look at the term to the right
            elif (i == n-1):
                val = valuation(self[i-1,i], p)    ## Look at the term above
            else:
                val = valuation(self[i,i+1] + self[i-1,i], p)    ## Finds the valuation of the off-diagonal term since only one isn't zero
        else:
            val = valuation(self[i,i], p)

        ## Sort the indices into disjoint sets by their valuation
        if (val == 0):
            S0 += [i]
        elif (val == 1):
            S1 += [i]
        elif (val >= 2):
            S2 += [i]



    ## Check that the form is primitive...
    if (S0 == []):
        print " Using Q = " + str(self)
        print " and p = " + str(p)
        raise RuntimeError, "Oops! The form is not primitive!"


    ## DIAGNOSTIC
    #print "\n Entering BII routine "
    #print " S0 is " + str(S0)
    #print " S1 is " + str(S1)
    #print " S2 is " + str(S2)



    ## Note: The following lines assume that NZvec is an ordered vector.  Should check this here...

    ## DIAGNOSTIC
    #print " m = " + str(m) + "   p = " + str(p)
    #print " S0 = " + str(S0) + "   NZvec = " str(NZvec) "   IsDisjoint = " + str(len(Set(S0 + S1).intersection(Set(NZvec))) == 0)
    #print " len(S0) = " + str(len(S0))


    p2 = p * p

    ## Make the form Qnew for the reduction procedure
    if ((m % (p2) == 0) and (S2 != []) and (len(Set(S0 + S1).intersection(Set(NZvec))) == 0)):   ##  <=====  CHECK THIS!!! ****  Should this be (S0 U S1) instead of S0 ???
        #print "Test 1.1 "

        Qnew = deepcopy(self)
        #Qnew = QuadraticForm(self.base_ring(), self.dim(), self.entries)

        ## Run over all indices with i, while keeping track of the next index in S2 as S2[j].
        j = 0       ## The index for following the increasing set (of indices) S2
        for i in range(n):
            ## Note: Short circuit && is necessary here:
            if ((j < len(S2)) and (i == S2[j])):                  ## i is in S2
                j += 1
                Qnew[i,i] = Qnew[i,i] / p2
                if ((p == 2) and (i < n-1)):
                    Qnew[i+1,i] = Qnew[i+1,i] / p2
                    Qnew[i,i+1] = Qnew[i,i+1] / p2

        ## DIAGNOSTIC
        #print "\n\n Check of Bad-type II reduction: \n";
        #print " Q is " + str(self)
        #print " Qnew is " + str(Qnew)


        ## Do the reduction
        ## (Need 2 steps for each case since we can't take negative powers... =| )
        ## ------------------------------------------------------------------------
        ##new_Zvec = VectorComplement(Zvec, list(Set(S0 + S1)))
        new_Zvec = list(Set([i  for i in Zvec  if i not in S0 + S1]))


        ## DIAGNOSTIC
        #print "  m = " << m << ",  m/p2 = " + str(m / p2)
        #print "  new_Zvec = " + str(new_Zvec)
        #print "  NZvec = " + str(NZvec)
        #print "  local_density_congruence(Qnew, p, m / p2, new_Zvec, NZvec) = " \
        #   + str(local_density_congruence(Qnew, p, m / p2, new_Zvec, NZvec))
        #print "  local_density_congruence(Qnew, p, m / p2, List(Set(S2 + new_Zvec)), NZvec)) = " \
        #   + str(local_density_congruence(Qnew, p, m / p2, List(Set(S2 + new_Zvec)), NZvec))
        #print "  count_local_type(Qnew, 2, 2, m=1, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Qnew, 2, 2, 1, 0, Zvec, NZvec))
        #print "  count_local_type(Qnew, 2, 3, m=1, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Qnew, 2, 3, 1, 0, Zvec, NZvec))
        #print "  count_local_type(Qnew, 2, 4, m=1, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Qnew, 2, 4, 1, 0, Zvec, NZvec))
        #print "  count_local_type(Qnew, 2, 5, m=1, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Qnew, 2, 5, 1, 0, Zvec, NZvec))
        #print "  count_local_type(Q, 2, 2, m=2, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Q, 2, 2, 2, 0, Zvec, NZvec))
        #print "  count_local_type(Q, 2, 3, m=2, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Q, 2, 3, 2, 0, Zvec, NZvec))
        #print "  count_local_type(Q, 2, 4, m=2, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Q, 2, 4, 2, 0, Zvec, NZvec))
        #print "  count_local_type(Q, 2, 5, m=2, 0, trimZvec, trimNZvec) = " \
        #   + str(count_local_type(Q, 2, 5, 2, 0, Zvec, NZvec))


        if (n > len(S2) + 2):
            return QQ(1 / p**(n - len(S2) - 2)) \
                * (Qnew.local_density_congruence(p, m / p2, new_Zvec, NZvec) \
                - Qnew.local_density_congruence(p, m / p2, List(Set(S2 + new_Zvec)), NZvec))
        else:
            return QQ(p**(len(S2) + 2 - n)) \
                * (Qnew.local_density_congruence(p, m / p2, new_Zvec, NZvec) \
                - Qnew.local_density_congruence(p, m / p2, List(Set(S2 + new_Zvec)), NZvec))

    else:
        return 0





def local_bad_density_congruence(self, p, m, Zvec, NZvec):
    """
    Finds the Bad-type local density of Q representing
    m at p, allowing certain congruence conditions mod p.
    """
    return self.local_badI_density_congruence(p, m, Zvec, NZvec) + self.local_badII_density_congruence(p, m, Zvec, NZvec)




#########################################################
## local_density and local_density_congruence routines ##
#########################################################

def local_density_congruence(self, p, m, Zvec, NZvec):
    """
    Finds the local density of Q representing m at p,
    allowing certain congruence conditions mod p.

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_density_congruence(p=2, m=1, Zvec=[], NZvec=[])
        1
        sage: Q.local_density_congruence(p=3, m=1, Zvec=[], NZvec=[])
        8/9
        sage: Q.local_density_congruence(p=5, m=1, Zvec=[], NZvec=[])
        24/25
        sage: Q.local_density_congruence(p=7, m=1, Zvec=[], NZvec=[])
        48/49
        sage: Q.local_density_congruence(p=11, m=1, Zvec=[], NZvec=[])
        120/121

    """
    return self.local_good_density_congruence(p, m, Zvec, NZvec) \
                + self.local_zero_density_congruence(p, m, Zvec, NZvec) \
                + self.local_bad_density_congruence(p, m, Zvec, NZvec)



def local_primitive_density_congruence(self, p, m, Zvec, NZvec):
    """
    Finds the primitive local density of Q representing
    m at p, allowing certain congruence conditions mod p.

    Note: The following routine is not used internally, but is included for consistency.
    """
    return self.local_good_density_congruence(p, m, Zvec, NZvec) \
                + self.local_bad_density_congruence(p, m, Zvec, NZvec)


