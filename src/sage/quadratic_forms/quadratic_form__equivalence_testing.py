"""
Equivalence Testing
"""
from sage.misc.mrange import mrange
from sage.rings.arith import hilbert_symbol, prime_divisors, is_prime, valuation, GCD, legendre_symbol
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from quadratic_form import is_QuadraticForm


from sage.env import SAGE_LOCAL

import tempfile, os, copy

from random import random


################################################################################
## Routines to test if two quadratic forms over ZZ are globally equivalent.   ##
## (For now, we require both forms to be positive definite.)                  ##
################################################################################

def is_globally_equivalent__souvigner(self, other, return_transformation=False):
    """
    Uses the Souvigner code to compute the number of automorphisms.

    INPUT:
        a QuadraticForm

    OUTPUT:
        boolean, and optionally a matrix

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q1 = QuadraticForm(ZZ, 3, [8, 6, 5, 3, 4, 2])
        sage: M = Q.is_globally_equivalent__souvigner(Q1, True) ; M   # optional -- souvigner
        [ 0  0 -1]
        [ 1  0  0]
        [-1  1  1]
        sage: Q1(M) == Q                                              # optional -- souvigner
        True

    """
    ## Write an input text file
    F_filename = '/tmp/tmp_isom_input' + str(random()) + ".txt"
    F = open(F_filename, 'w')
    #F = tempfile.NamedTemporaryFile(prefix='tmp_isom_input', suffix=".txt")  ## This failed because it may have hyphens, which are interpreted badly by the Souvigner code.
    F.write("\n #1 \n")

    ## Write the first form
    n = self.dim()
    F.write(str(n) + "x0 \n")      ## Use the lower-triangular form
    for i in range(n):
        for j in range(i+1):
            if i == j:
                F.write(str(2 * self[i,j]) + " ")
            else:
                F.write(str(self[i,j]) + " ")
        F.write("\n")

    ## Write the second form
    F.write("\n")
    n = self.dim()
    F.write(str(n) + "x0 \n")      ## Use the lower-triangular form
    for i in range(n):
        for j in range(i+1):
            if i == j:
                F.write(str(2 * other[i,j]) + " ")
            else:
                F.write(str(other[i,j]) + " ")
        F.write("\n")
    F.flush()
    #print "Input filename = ", F.name
    #os.system("less " + F.name)

    ## Call the Souvigner automorphism code
    souvigner_isom_path = os.path.join(SAGE_LOCAL,'bin','Souvigner_ISOM')
    G1 = tempfile.NamedTemporaryFile(prefix='tmp_isom_ouput', suffix=".txt")
    #print "Output filename = ", G1.name
    #print  "Executing the shell command:   " + souvigner_isom_path + " '" +  F.name + "' > '" + G1.name + "'"
    os.system(souvigner_isom_path + " '" +  F.name + "' > '" + G1.name +"'")


    ## Read the output
    G2 = open(G1.name, 'r')
    line = G2.readline()
    if line.startswith("Error:"):
        raise RuntimeError("There is a problem using the souvigner code...  " + line)
    elif line.find("not isomorphic") != -1:     ## Checking if this text appears, if so then they're not isomorphic!
        return False
    else:
        ## Decide whether to read the transformation matrix, and return true
        if not return_transformation:
            F.close()
            G1.close()
            G2.close()
            os.system("rm -f " + F_filename)
            return True
        else:
            ## Try to read the isomorphism matrix
            M = Matrix(ZZ, n, n)
            for i in range(n):
                new_row_text = G2.readline().split()
                #print new_row_text
                for j in range(n):
                    M[i,j] = new_row_text[j]

            ## Remove temporary files and return the value
            F.close()
            G1.close()
            G2.close()
            os.system("rm -f " + F_filename)
            if return_transformation:
                return M.transpose()
            else:
                return True
            #return True, M

    ## Raise and error if we're here:
    raise RuntimeError("Oops! There is a problem...")



def is_globally_equivalent_to(self, other, return_matrix=False, check_theta_to_precision='sturm', check_local_equivalence=True):
    """
    Determines if the current quadratic form is equivalent to the
    given form over ZZ.  If return_matrix is True, then we also return
    the transformation matrix M so that self(M) == other.

    INPUT:
        a QuadraticForm

    OUTPUT:
        boolean, and optionally a matrix

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: M = Matrix(ZZ, 4, 4, [1,2,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1])
        sage: Q1 = Q(M)
        sage: Q.(Q1)                                                    # optional -- souvigner
        True
        sage: MM = Q.is_globally_equivalent_to(Q1, return_matrix=True)  # optional -- souvigner
        sage: Q(MM) == Q1                                               # optional -- souvigner
        True

    ::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q2 = QuadraticForm(ZZ, 3, [2, 1, 2, 2, 1, 3])
        sage: Q3 = QuadraticForm(ZZ, 3, [8, 6, 5, 3, 4, 2])
        sage: Q1.is_globally_equivalent_to(Q2)                          # optional -- souvigner
        False
        sage: Q1.is_globally_equivalent_to(Q3)                          # optional -- souvigner
        True
        sage: M = Q1.is_globally_equivalent_to(Q3, True) ; M            # optional -- souvigner
        [-1 -1  0]
        [ 1  1  1]
        [-1  0  0]
        sage: Q1(M) == Q3                                               # optional -- souvigner
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1, -1])
        sage: Q.is_globally_equivalent_to(Q)
        Traceback (most recent call last):
        ...
        ValueError: not a definite form in QuadraticForm.is_globally_equivalent_to()

    """
    ## only for definite forms
    if not self.is_definite():
        raise ValueError("not a definite form in QuadraticForm.is_globally_equivalent_to()")

    ## Check that other is a QuadraticForm
    #if not isinstance(other, QuadraticForm):
    if not is_QuadraticForm(other):
        raise TypeError("Oops!  You must compare two quadratic forms, but the argument is not a quadratic form. =(")


    ## Now use the Souvigner code by default! =)
    return other.is_globally_equivalent__souvigner(self, return_matrix)    ## Note: We switch this because the Souvigner code has the opposite mapping convention to us.  (It takes the second argument to the first!)



    ## ----------------------------------  Unused Code below  ---------------------------------------------------------

    ## Check if the forms are locally equivalent
    if (check_local_equivalence == True):
        if not self.is_locally_equivalent_to(other):
            return False

    ## Check that the forms have the same theta function up to the desired precision (this can be set so that it determines the cusp form)
    if check_theta_to_precision is not None:
        if self.theta_series(check_theta_to_precision, var_str='', safe_flag=False) != other.theta_series(check_theta_to_precision, var_str='', safe_flag=False):
            return False


    ## Make all possible matrices which give an isomorphism -- can we do this more intelligently?
    ## ------------------------------------------------------------------------------------------

    ## Find a basis of short vectors for one form, and try to match them with vectors of that length in the other one.
    basis_for_self, self_lengths = self.basis_of_short_vectors(show_lengths=True)
    max_len = max(self_lengths)
    short_vectors_of_other = other.short_vector_list_up_to_length(max_len + 1)

    ## Make the matrix A:e_i |--> v_i to our new basis.
    A = Matrix(basis_for_self).transpose()
    Q2 = A.transpose() * self.matrix() * A       ## This is the matrix of 'self' in the new basis
    Q3 = other.matrix()

    ## Determine all automorphisms
    n = self.dim()
    Auto_list = []

    ## DIAGNOSTIC
    #print "n = " + str(n)
    #print "pivot_lengths = " + str(pivot_lengths)
    #print "vector_list_by_length = " + str(vector_list_by_length)
    #print "length of vector_list_by_length = " + str(len(vector_list_by_length))

    for index_vec in mrange([len(short_vectors_of_other[self_lengths[i]])  for i in range(n)]):
        M = Matrix([short_vectors_of_other[self_lengths[i]][index_vec[i]]   for i in range(n)]).transpose()
        if M.transpose() * Q3 * M == Q2:
            if return_matrix:
                return A * M.inverse()
            else:
                return True

    ## If we got here, then there is no isomorphism
    return False


def is_locally_equivalent_to(self, other, check_primes_only=False, force_jordan_equivalence_test=False):
    """
    Determines if the current quadratic form (defined over ZZ) is
    locally equivalent to the given form over the real numbers and the
    `p`-adic integers for every prime p.

    This works by comparing the local Jordan decompositions at every
    prime, and the dimension and signature at the real place.

    INPUT:
        a QuadraticForm

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q2 = QuadraticForm(ZZ, 3, [2, 1, 2, 2, 1, 3])
        sage: Q1.is_globally_equivalent_to(Q2)                           # optional -- souvigner
        False
        sage: Q1.is_locally_equivalent_to(Q2)
        True

    """
    ## TO IMPLEMENT:
    if self.det() == 0:
        raise NotImplementedError("OOps!  We need to think about whether this still works for degenerate forms...  especially check the signature.")

    ## Check that both forms have the same dimension and base ring
    if (self.dim() != other.dim()) or (self.base_ring() != other.base_ring()):
        return False

    ## Check that the determinant and level agree
    if (self.det() != other.det()) or (self.level() != other.level()):
        return False

    ## -----------------------------------------------------

    ## Test equivalence over the real numbers
    if self.signature() != other.signature():
        return False

    ## Test equivalence over Z_p for all primes
    if (self.base_ring() == ZZ) and (force_jordan_equivalence_test == False):

        ## Test equivalence with Conway-Sloane genus symbols (default over ZZ)
        if self.CS_genus_symbol_list() != other.CS_genus_symbol_list():
            return False
    else:
        ## Test equivalence via the O'Meara criterion.
        for p in prime_divisors(ZZ(2) * self.det()):
            #print "checking the prime p = ", p
            if not self.has_equivalent_Jordan_decomposition_at_prime(other, p):
                return False

    ## All tests have passed!
    return True




def has_equivalent_Jordan_decomposition_at_prime(self, other, p):
    """
    Determines if the given quadratic form has a Jordan decomposition
    equivalent to that of self.

    INPUT:
        a QuadraticForm

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 0, -1, 1, 0, 3])
        sage: Q2 = QuadraticForm(ZZ, 3, [1, 0, 0, 2, -2, 6])
        sage: Q3 = QuadraticForm(ZZ, 3, [1, 0, 0, 1, 0, 11])
        sage: [Q1.level(), Q2.level(), Q3.level()]
        [44, 44, 44]
        sage: Q1.has_equivalent_Jordan_decomposition_at_prime(Q2,2)
        False
        sage: Q1.has_equivalent_Jordan_decomposition_at_prime(Q2,11)
        False
        sage: Q1.has_equivalent_Jordan_decomposition_at_prime(Q3,2)
        False
        sage: Q1.has_equivalent_Jordan_decomposition_at_prime(Q3,11)
        True
        sage: Q2.has_equivalent_Jordan_decomposition_at_prime(Q3,2)
        True
        sage: Q2.has_equivalent_Jordan_decomposition_at_prime(Q3,11)
        False

    """
    ## Sanity Checks
    #if not isinstance(other, QuadraticForm):
    if not isinstance(other, type(self)):
        raise TypeError("Oops!  The first argument must be of type QuadraticForm.")
    if not is_prime(p):
        raise TypeError("Oops!  The second argument must be a prime number.")

    ## Get the relevant local normal forms quickly
    self_jordan = self.jordan_blocks_by_scale_and_unimodular(p, safe_flag= False)
    other_jordan = other.jordan_blocks_by_scale_and_unimodular(p, safe_flag=False)

    ## DIAGNOSTIC
    #print "self_jordan = ", self_jordan
    #print "other_jordan = ", other_jordan


    ## Check for the same number of Jordan components
    if len(self_jordan) != len(other_jordan):
        return False


    ## Deal with odd primes:  Check that the Jordan component scales, dimensions, and discriminants are the same
    if p != 2:
        for i in range(len(self_jordan)):
            if (self_jordan[i][0] != other_jordan[i][0]) \
               or (self_jordan[i][1].dim() != other_jordan[i][1].dim()) \
               or (legendre_symbol(self_jordan[i][1].det() * other_jordan[i][1].det(), p) != 1):
                return False

        ## All tests passed for an odd prime.
        return True


    ## For p = 2:  Check that all Jordan Invariants are the same.
    elif p == 2:

        ## Useful definition
        t = len(self_jordan)          ## Define t = Number of Jordan components


        ## Check that all Jordan Invariants are the same (scale, dim, and norm)
        for i in range(t):
            if (self_jordan[i][0] != other_jordan[i][0]) \
               or (self_jordan[i][1].dim() != other_jordan[i][1].dim()) \
               or (valuation(GCD(self_jordan[i][1].coefficients()), p) != valuation(GCD(other_jordan[i][1].coefficients()), p)):
                return False

        ## DIAGNOSTIC
        #print "Passed the Jordan invariant test."


        ## Use O'Meara's isometry test 93:29 on p277.
        ## ------------------------------------------

        ## List of norms, scales, and dimensions for each i
        scale_list = [ZZ(2)**self_jordan[i][0]  for i in range(t)]
        norm_list = [ZZ(2)**(self_jordan[i][0] + valuation(GCD(self_jordan[i][1].coefficients()), 2))  for i in range(t)]
        dim_list = [(self_jordan[i][1].dim())  for i in range(t)]

        ## List of Hessian determinants and Hasse invariants for each Jordan (sub)chain
        ## (Note: This is not the same as O'Meara's Gram determinants, but ratios are the same!)  -- NOT SO GOOD...
        ## But it matters in condition (ii), so we multiply all by 2 (instead of dividing by 2 since only square-factors matter, and it's easier.)
        j = 0
        self_chain_det_list = [ self_jordan[j][1].Gram_det() * (scale_list[j]**dim_list[j])]
        other_chain_det_list = [ other_jordan[j][1].Gram_det() * (scale_list[j]**dim_list[j])]
        self_hasse_chain_list = [ self_jordan[j][1].scale_by_factor(ZZ(2)**self_jordan[j][0]).hasse_invariant__OMeara(2) ]
        other_hasse_chain_list = [ other_jordan[j][1].scale_by_factor(ZZ(2)**other_jordan[j][0]).hasse_invariant__OMeara(2) ]

        for j in range(1, t):
            self_chain_det_list.append(self_chain_det_list[j-1] * self_jordan[j][1].Gram_det() * (scale_list[j]**dim_list[j]))
            other_chain_det_list.append(other_chain_det_list[j-1] * other_jordan[j][1].Gram_det() * (scale_list[j]**dim_list[j]))
            self_hasse_chain_list.append(self_hasse_chain_list[j-1] \
                                         * hilbert_symbol(self_chain_det_list[j-1], self_jordan[j][1].Gram_det(), 2) \
                                         * self_jordan[j][1].hasse_invariant__OMeara(2))
            other_hasse_chain_list.append(other_hasse_chain_list[j-1] \
                                          * hilbert_symbol(other_chain_det_list[j-1], other_jordan[j][1].Gram_det(), 2) \
                                          * other_jordan[j][1].hasse_invariant__OMeara(2))


        ## SANITY CHECK -- check that the scale powers are strictly increasing
        for i in range(1, len(scale_list)):
            if scale_list[i-1] >= scale_list[i]:
                   raise RuntimeError("Oops!  There is something wrong with the Jordan Decomposition -- the given scales are not strictly increasing!")


        ## DIAGNOSTIC
        #print "scale_list = ", scale_list
        #print "norm_list = ", norm_list
        #print "dim_list = ", dim_list
        #print
        #print "self_chain_det_list = ", self_chain_det_list
        #print "other_chain_det_list = ", other_chain_det_list
        #print "self_hasse_chain_list = ", self_hasse_chain_list
        #print "other_hasse_chain_det_list = ", other_hasse_chain_list


        ## Test O'Meara's two conditions
        for i in range(t-1):

            ## Condition (i): Check that their (unit) ratio is a square (but it suffices to check at most mod 8).
            modulus = norm_list[i] * norm_list[i+1] / (scale_list[i] ** 2)
            if modulus > 8:
                   modulus = 8
            if (modulus > 1) and (((self_chain_det_list[i] / other_chain_det_list[i]) % modulus) != 1):
                #print "Failed when i =", i, " in condition 1."
                return False

            ## Check O'Meara's condition (ii) when appropriate
            if norm_list[i+1] % (4 * norm_list[i]) == 0:
                if self_hasse_chain_list[i] * hilbert_symbol(norm_list[i] * other_chain_det_list[i], -self_chain_det_list[i], 2) \
                       != other_hasse_chain_list[i] * hilbert_symbol(norm_list[i], -other_chain_det_list[i], 2):      ## Nipp conditions
                    #print "Failed when i =", i, " in condition 2."
                    return False


        ## All tests passed for the prime 2.
        return True

    else:
        raise TypeError("Oops!  This should not have happened.")



def isometry(self, other):
    """
    Given two equivalent quadratic forms, computes an isometry from one to the other.
    Note: Currently only works for positive definite quadratic forms. This is because isometry() relies on 
    the method short_vector_list_up_to_length() which is only defined for positive definite forms.

    INPUT:
        self: a positive definite QuadraticForm
        other: a positive definite QuadraticForm
        NOTE: This algorithm may not terminate if self is not equivalent to other!

    OUTPUT:
        a matrix T representing the isometry transformation, such that if QM is the gram matrix of Q and
        FM is the gram matrix of F, then QM == T.transpose() * FM * T yields True.

    EXAMPLES::
    sage: Q = DiagonalQuadraticForm(QQ, [1, 1, 2])
    sage: F = DiagonalQuadraticForm(QQ, [2, 2, 2])
    sage: T = Q.isometry(F)
    sage: T
    [ 1/2  1/2    0]
    [-1/2  1/2    0]
    [   0    0    1]
    
    sage: Q.Gram_matrix() == T.transpose() * F.Gram_matrix() * T
    True
    
    sage: T = F.isometry(Q)
    sage: T
    [ 1 -1  0]
    [ 1  1  0]
    [ 0  0  1]
    
    sage: F.Gram_matrix() == T.transpose() * Q.Gram_matrix() * T
    True
    
    sage: A = DiagonalQuadraticForm(QQ, [1,5])
    sage: B = QuadraticForm(QQ, 2, [1, 12, 81])
    sage: T = A.isometry(B)
    sage: T
    [  1  -2]
    [  0 1/3]
    
    sage: A.Gram_matrix() == T.transpose() * B.Gram_matrix() * T
    True
    
    sage: C = DiagonalQuadraticForm(QQ, [1, 5, 9])
    sage: D = DiagonalQuadraticForm(QQ, [6, 30, 1])
    sage: T = D.isometry(C)
    sage: T
    [  1  -5   0]
    [  1   1   0]
    [  0   0 1/3]
    
    sage: D.Gram_matrix() == T.transpose() * C.Gram_matrix() * T
    True
    
    sage: E = DiagonalQuadraticForm(QQ, [1, 1])
    sage: F = QuadraticForm(QQ, 2, [17, 94, 130])
    sage: T = E.isometry(F)
    sage: T
    [189/17 -43/17]
    [    -4      1]

    sage: E.Gram_matrix() == T.transpose() * F.Gram_matrix() * T
    True
    """
    from sage.matrix.constructor import Matrix
    ##Define a method that determines whether or not a matrix is diagonal.
    def is_diagonal(matrix):
        dim = matrix.dimensions()[0]
        for i in range(dim):
            for j in range(dim):
                if (i != j) and (matrix[i][j] != 0):
                    return False
        return True

    Q = copy.deepcopy(self)
    F = copy.deepcopy(other)

    if not Q.is_positive_definite() or not F.is_positive_definite():
        raise TypeError("Both Quadratic Forms must be positive definite.")

    if not is_QuadraticForm(other):
        raise TypeError("First argument must be a Quadratic Form.")

    n = Q.dim()

    ##If either form is not diagonal, diagonalize it.
    q_diagonal_transform = f_diagonal_transform = Matrix.identity(n)
    if not is_diagonal(Q.Gram_matrix()):
        Q, q_diagonal_transform = Q.rational_diagonal_form(True)
    if not is_diagonal(F.Gram_matrix()):
        F, f_diagonal_transform = F.rational_diagonal_form(True)
        
    #This is the method that does all the work to compute the isometry.
    transform = diagonal_isometry(Q,F)
    
    return f_diagonal_transform * transform * q_diagonal_transform.inverse()


##Helper method for isometry.
def diagonal_isometry(V, W):
##Computes an isometry between diagonal equivalent forms V and W.
    from sage.functions.other import sqrt
    from quadratic_form import DiagonalQuadraticForm
    from sage.matrix.constructor import matrix, column_matrix, Matrix
    from sage.modules.free_module_element import vector
    #----HELPER METHODS-----#

    def diagonal_qf_pull_off_n_terms(gram_matrix, n, ring):
        ##Returns the diagonal quadratic form corresponding to the given gram matrix
        ##with the first n terms removed.
        diagonal = gram_matrix.diagonal();
        return DiagonalQuadraticForm(ring, diagonal[n:])

    def compute_gram_matrix_from_basis(Q, basis):
        ##Computes the gram matrix by treating the columns of basis as basis vectors.
        n = Q.dim()
        rows = [];
        for i in range(n):
            rows.append([Q.bilinear_map(basis.column(i), basis.column(j)) for j in range(n)])
        return Matrix(rows)

    def modify_basis(basis, v, i):
        ##Modifies basis to include v (a linear combination of the basis vectors) at column i in the basis.
        b = copy.deepcopy(basis)
        column = vector(QQ, b.dimensions()[0])
        for j in range(len(v)):
            column += b.column(i + j) * v[j]
        b.set_column(i, column)
        return b

    def gram_schmidt(matrix, fixed_vector_index, gram_matrix):
        ##Orthogonalizes the columns of matrix against a fixed column vector.
        n = matrix.dimensions()[0]
        vectors = [0] * n
        fixed_vec = matrix.column(fixed_vector_index)
        vectors[fixed_vector_index] = fixed_vec
    
        for i in range(fixed_vector_index):
            vectors[i] = matrix.column(i)
        for i in range(fixed_vector_index + 1, n):
            vectors[i] = matrix.column(i) - (gram_matrix[i, fixed_vector_index] / gram_matrix[fixed_vector_index, fixed_vector_index]) * fixed_vec

        return column_matrix(vectors)
    
    def vectors_of_common_length_dev(Q, F, qb, fb, iteration):
        ##Returns a pair of vectors v1 and v2, such that Q(v1) = F(v2), 
        ##and qb and fb are bases when their columns at index iteration are replaced with
        ##v1 and v2 respectively. (i.e. qb and fb contain no all-zero rows)
        def zero_row(matrix, col, i):
            ##Sets the ith column of matrix to col and returns true if matrix has an all zero row.
            m = modify_basis(matrix, col, i)
            rows, cols = m.dimensions()
            z = [0] * cols
            for i in range(rows):
                if m.row(i).list() == z:
                    return True
            return False

        i = max(Q.Gram_matrix()[0][0], F.Gram_matrix()[0][0])
        n = 100
        #NOTE: If the given forms are not equivalent, this will be an infinite loop.
        while True:
            q_vecs, f_vecs = Q.short_vector_list_up_to_length(n, real_entry_flag=True), F.short_vector_list_up_to_length(n, real_entry_flag=True)
            while i < n:
                ##Find a value that is represented by both forms
                if q_vecs[i] and f_vecs[i]:
                    ##Find a pair of vectors such that if each becomes a column of qb and fb respectively,
                    ##neither matrix will have an all-zero row.
                    for j in range(len(q_vecs[i])):
                        for k in range(len(f_vecs[i])):
                            if not zero_row(qb, q_vecs[i][j], iteration) and not zero_row(fb, f_vecs[i][k], iteration):
                                return (q_vecs[i][j], f_vecs[i][k])
                i += 1
            n += 100

    #--- END OF HELPER METHODS ---#

    Q, F = copy.deepcopy(V), copy.deepcopy(W)
    n = Q.dim()
    
    q_basis, f_basis = Matrix.identity(n), Matrix.identity(n)
    for i in range(n-1):
        ##Find vectors v and w such that Q(v) = F(w)
        v, w = vectors_of_common_length_dev(Q, F, q_basis, f_basis, i)
        #print("Find vectors {0} and {1} such that Q(v) = F(w)".format(v, w))        

        ##Modify the bases to include v and w.
        q_basis = modify_basis(q_basis, v, i)
        f_basis = modify_basis(f_basis, w, i)
        #print("Modified bases:\nQb =\n{0}, \nFb =\n{1}".format(q_basis, f_basis))
        
        ##Compute the gram matrices with respect to the modified bases.
        QM = compute_gram_matrix_from_basis(V, q_basis)
        FM = compute_gram_matrix_from_basis(W, f_basis)
        #print("Gram matrices with respect to modified bases:\nQM =\n{0}, \nFM =\n{1}".format(QM,FM))
        
        ##Ensure that the bases are orthogonal, so the gram matrices will be diagonal.
        q_basis = gram_schmidt(q_basis, i, QM)
        f_basis = gram_schmidt(f_basis, i, FM)
        #print("Orthogonal bases:\nQb =\n{0}, \nFb =\n{1}".format(q_basis, f_basis))

        ##Compute the gram matrices with respect to the orthogonal bases.
        QM = compute_gram_matrix_from_basis(V, q_basis)
        FM = compute_gram_matrix_from_basis(W, f_basis)
        #print("Gram matrices with respect to orthogonal bases:\nQM =\n{0}, \nFM =\n{1}".format(QM,FM))        

        ##Pull off the first term and continue
        Q = diagonal_qf_pull_off_n_terms(QM, i+1, Q.base_ring())
        F = diagonal_qf_pull_off_n_terms(FM, i+1, F.base_ring())
        #print("Smaller qf's:\nQ =\n{0}, \nF = \n{1}".format(Q, F))
        
    ##Compute the final term as a special case:
    QM, FM = Q.Gram_matrix(), F.Gram_matrix()
    if QM[0][0] != 0:
        q_basis.set_col_to_multiple_of_col(n-1, n-1, sqrt(FM[0][0] / QM[0][0]))
    elif FM[0][0] != 0:
        f_basis.set_col_to_multiple_of_col(n-1, n-1, sqrt(QM[0][0] / FM[0][0]))

    return f_basis * q_basis.inverse()
