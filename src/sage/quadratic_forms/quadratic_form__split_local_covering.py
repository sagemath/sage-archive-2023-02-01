
#########################################################################
## Routines that look for a split local covering for a given quadratic ##
## form in 4 variables.                                                ##
#########################################################################

from copy import deepcopy

from sage.quadratic_forms.extras import extend_to_primitive
from sage.quadratic_forms.quadratic_form import QuadraticForm__constructor, is_QuadraticForm

from sage.rings.real_double import RDF
from sage.matrix.matrix_space import MatrixSpace
#from sage.matrix.matrix import Matrix
from sage.matrix.constructor import matrix
from sage.calculus.calculus import sqrt, floor, ceil
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import round
from sage.rings.arith import GCD


def vectors_by_length(self, bound):
    """
    Returns a list of short vectors together with their values.

    This is a maive algorithm which uses the Cholesky decomposition,
    but does not use the LLL-reduction algorithm.

    OUTPUT:
        A list L of length (bound + 1) whose entry L[i] is a list of
        all vectors of length i.

    Reference: This is a slightly modified version of Cohn's Algorithm
    2.7.5 in "A Course in Compuational Number Theory", with the
    increment step moved around and slightly reindexed to allow clean
    looping.

    Note: We could speed this up for very skew matrices by using LLL
    first, and then changing coordinates back, but for our purposes
    the simpler method is efficient enough. =)

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q.vectors_by_length(5)
        [[[0, 0]],
         [[0, -1], [-1, 0]],
         [[-1, -1], [1, -1]],
         [],
         [[0, -2], [-2, 0]],
         [[-1, -2], [1, -2], [-2, -1], [2, -1]]]

        sage: Q1 = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q1.vectors_by_length(5)
        [[[0, 0, 0, 0]],
         [[-1, 0, 0, 0]],
         [],
         [[0, -1, 0, 0]],
         [[-1, -1, 0, 0], [1, -1, 0, 0], [-2, 0, 0, 0]],
         [[0, 0, -1, 0]]]

    """
    Theta_Precision = bound               ## Unsigned long
    n = self.dim()

    ## Make the vector of vectors which have a given value
    ## (So theta_vec[i] will have all vectors v with Q(v) = i.)
    empty_vec_list = [[] for i in range(Theta_Precision + 1)]
    theta_vec = [[] for i in range(Theta_Precision + 1)]

    ## Initialize Q with zeros and Copy the Cholesky array into Q
    Q = self.cholesky_decomposition()


    ## 1. Initialize
    T = n * [RDF(0)]    ## Note: We index the entries as 0 --> n-1
    U = n * [RDF(0)]
    i = n-1
    T[i] = RDF(Theta_Precision)
    U[i] = RDF(0)

    L = n * [0]
    x = n * [0]
    Z = RDF(0)

    ## 2. Compute bounds
    Z = sqrt(T[i] / Q[i][i])
    L[i] = ZZ(floor(Z - U[i]))
    x[i] = ZZ(ceil(-Z - U[i]) - 0)

    done_flag = False
    Q_val_double = RDF(0)
    Q_val = 0 ## WARNING: Still need a good way of checking overflow for this value...

    ## Big loop which runs through all vectors
    while not done_flag:

        ## 3b. Main loop -- try to generate a complete vector x (when i=0)
        while (i > 0):
            #print " i = ", i
            #print " T[i] = ", T[i]
            #print " Q[i][i] = ", Q[i][i]
            #print " x[i] = ", x[i]
            #print " U[i] = ", U[i]
            #print " x[i] + U[i] = ", (x[i] + U[i])
            #print " T[i-1] = ", T[i-1]

            T[i-1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i])

            #print " T[i-1] = ",  T[i-1]
            #print " x = ", x
            #print

            i = i - 1
            U[i] = 0
            for j in range(i+1, n):
                U[i] = U[i] + Q[i][j] * x[j]

            ## Now go back and compute the bounds...
            ## 2. Compute bounds
            Z = sqrt(T[i] / Q[i][i])
            L[i] = ZZ(floor(Z - U[i]))
            x[i] = ZZ(ceil(-Z - U[i]) - 0)


        ## 4. Solution found (This happens when i = 0)
        #print "-- Solution found! --"
        #print " x = ", x
        #print " Q_val = Q(x) = ", Q_val
        Q_val_double = Theta_Precision - T[0] + Q[0][0] * (x[0] + U[0]) * (x[0] + U[0])
        Q_val = ZZ(floor(round(Q_val_double)))

        ## SANITY CHECK: Roundoff Error is < 0.001
        if abs(Q_val_double -  Q_val) > 0.001:
            print " x = ", x
            print " Float = ", Q_val_double, "   Long = ", Q_val
            raise RuntimeError, "The roundoff error is bigger than 0.001, so we should use more precision somehwere..."

        #print " Float = ", Q_val_double, "   Long = ", Q_val, "  XX "
        #print " The float value is ", Q_val_double
        #print " The associated long value is ", Q_val

        if (Q_val <= Theta_Precision):
            #print " Have vector ",  x, " with value ", Q_val
            theta_vec[Q_val].append(deepcopy(x))


        ## 5. Check if x = 0, for exit condition. =)
        j = 0
        done_flag = True
        while (j < n):
            if (x[j] != 0):
                done_flag = False
            j += 1


        ## 3a. Increment (and carry if we go out of bounds)
        x[i] += 1
        while (x[i] > L[i]) and (i < n-1):
            i += 1
            x[i] += 1


    #print " Leaving ThetaVectors()"
    return theta_vec




def complementary_subform_to_vector(self, v):
    """
    Finds the (n-1)-dim'l quadratic form orthogonal to the vector v.

    Note: This is usually not a direct summand!

    Technical Notes: There is a minor difference in the cancellation
    code here (form the C++ version) since the notation Q[i,j] indexes
    coefficients of the quadratic polynomial here, not the symmetric
    matrix.  Also, it produces a better splitting now, for the full
    lattice (as opposed to a sublattice in the C++ code) since we
    now extend v to a unimodular matrix.

    EXAMPLES:
    sage: Q1 = DiagonalQuadraticForm(ZZ, [1,3,5,7])
    sage: Q1.complementary_subform_to_vector([1,0,0,0])
    Quadratic form in 3 variables over Integer Ring with coefficients:
    [ 3 0 0 ]
    [ * 5 0 ]
    [ * * 7 ]

    sage: Q1.complementary_subform_to_vector([1,1,0,0])
    Quadratic form in 3 variables over Integer Ring with coefficients:
    [ 12 0 0 ]
    [ * 5 0 ]
    [ * * 7 ]

    sage: Q1.complementary_subform_to_vector([1,1,1,1])
    Quadratic form in 3 variables over Integer Ring with coefficients:
    [ 624 -480 -672 ]
    [ * 880 -1120 ]
    [ * * 1008 ]

    """
    n = self.dim()

    ## Copy the quadratic form
    Q = deepcopy(self)

    ## Find the first non-zero component of v, and call it nz  (Note: 0 <= nz < n)
    nz = 0
    while (nz < n) and (v[nz] == 0):
        nz += 1

    ## Abort if v is the zero vector
    if nz == n:
        raise TypeError, "Oops, v cannot be the zero vector! =("

    ## Make the change of basis matrix
    new_basis = extend_to_primitive(matrix(ZZ,n,1,v))

    ## Change Q (to Q1) to have v as its nz-th basis vector
    Q1 = Q(new_basis)

    ## Pick out the value Q(v) of the vector
    d = Q1[0, 0]

    #print Q1

    ## For each row/column, perform elementary operations to cancel them out.
    for i in range(1,n):

        ## Check if the (i,0)-entry is divisible by d,
        ## and stretch its row/column if not.
        if Q1[i,0] % d != 0:
            Q1 = Q1.multiply_variable(d / GCD(d, Q1[i, 0]/2), i)

        #print "After scaling at i =", i
        #print Q1

        ## Now perform the (symmetric) elementary operations to cancel out the (i,0) entries/
        Q1 = Q1.add_symmetric(-(Q1[i,0]/2) / (GCD(d, Q1[i,0]/2)), i, 0)

        #print "After cancelling at i =", i
        #print Q1

    ## Check that we're done!
    done_flag = True
    for i in range(1, n):
        if Q1[0,i] != 0:
            done_flag = False

    if done_flag == False:
        raise RuntimeError, "There is a problem cancelling out the matrix entries! =O"


    ## Return the complementary matrix
    return Q1.extract_variables(range(1,n))



def split_local_cover(self):
    """
    Tries to find subform of the given (positive definite quaternary)
    quadratic form Q of the form

        d*x^2 + T(y,z,w)

    where d > 0 is as small as possible.

    This is done by exhaustive search on small vectors, and then
    comparing the local conditions of its sum with it's complementary
    lattice and the original quadratic form Q.

    EXAMPLES:
        sage: Q1 = DiagonalQuadraticForm(ZZ, [7,5,3])
        sage: Q1.split_local_cover()
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 0 0 ]
        [ * 7 0 ]
        [ * * 5 ]

    """
    ## 0. If a split local cover already exists, then return it.
    if hasattr(self, "__split_local_cover"):
        if is_QuadraticForm(self.__split_local_cover):  ## Here the computation has been done.
            return self.__split_local_cover
        elif isinstance(self.__split_local_cover, Integer):    ## Here it indexes the values already tried!
            current_length = self.__split_local_cover + 1
            Length_Max = current_length + 5
    else:
        current_length = 1
        Length_Max = 6

    ## 1. Find a range of new vectors
    all_vectors = self.vectors_by_length(Length_Max)
    current_vectors = all_vectors[current_length]

    ## Loop until we find a split local cover...
    while True:

        ## 2. Check if any of the primitive ones produce a split local cover
        for v in current_vectors:
            #print "current length = ", current_length
            #print "v = ", v
            Q = QuadraticForm__constructor(ZZ, 1, [current_length]) + self.complementary_subform_to_vector(v)
            #print Q
            if Q.local_representation_conditions() == self.local_representation_conditions():
                self.__split_local_cover = Q
                return Q

        ## 3. Save what we have checked and get more vectors.
        self.__split_local_cover = current_length
        current_length += 1
        if current_length >= len(all_vectors):
            Length_Max += 5
            all_vectors = self.vectors_by_length(Length_Max)
        current_vectors = all_vectors[current_length]

