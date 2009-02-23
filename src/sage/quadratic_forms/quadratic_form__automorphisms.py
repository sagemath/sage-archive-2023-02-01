
from sage.interfaces.gp import gp
from sage.matrix.constructor import Matrix
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.misc.mrange import mrange

from sage.modules.free_module_element import vector
from sage.rings.arith import GCD
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import SAGE_ROOT

import tempfile, os
from random import random


def basis_of_short_vectors(self, show_lengths=False, safe_flag=True):
    """
    Return a basis for ZZ^n made of vectors with minimal lengths Q(v).

    The safe_flag allows us to select whether we want a copy of the
    output, or the original output.  By default safe_flag = True, so
    we return a copy of the cached information.  If this is set to
    False, then the routine is much faster but the return values are
    vulnerable to being corrupted by the user.

    OUTPUT:
        a list of vectors, and optionally a list of values for each vector.

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.basis_of_short_vectors()
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        sage: Q.basis_of_short_vectors(True)
        ([(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)], [1, 3, 5, 7])
    """
    ## Try to use the cached results
    try:
        ## Return the appropriate result
        if show_lengths:
            if safe_flag:
                return deep_copy(self.__basis_of_short_vectors), deepcopy(self.__basis_of_short_vectors_lengths)
            else:
                return self.__basis_of_short_vectors, self.__basis_of_short_vectors_lengths
        else:
            if safe_flag:
                return deepcopy(self.__basis_of_short_vectors)
            else:
                return deepcopy(self.__basis_of_short_vectors)
    except:
        pass


    ## Set an upper bound for the number of vectors to consider
    Max_number_of_vectors = 10000

    ## Generate a PARI matrix string for the associated Hessian matrix
    M_str = str(gp(self.matrix()))


    ## Run through all possible minimal lengths to find a spanning set of vectors
    n = self.dim()
    #MS = MatrixSpace(QQ, n)
    M1 = Matrix([[0]])
    vec_len = 0
    while M1.rank() < n:

        ## DIAGONSTIC
        #print
        #print "Starting with vec_len = ", vec_len
        #print "M_str = ", M_str

        vec_len += 1
        gp_mat = gp.qfminim(M_str, vec_len, Max_number_of_vectors)[3].mattranspose()
        number_of_vecs = ZZ(gp_mat.matsize()[1])
        vector_list = []
        for i in range(number_of_vecs):
            #print "List at", i, ":", list(gp_mat[i+1,])
            new_vec = vector([ZZ(x)  for x in list(gp_mat[i+1,])])
            vector_list.append(new_vec)


        ## DIAGNOSTIC
        #print "number_of_vecs = ", number_of_vecs
        #print "vector_list = ", vector_list


        ## Make a matrix from the short vectors
        if len(vector_list) > 0:
            M1 = Matrix(vector_list)


        ## DIAGNOSTIC
        #print "matrix of vectors = \n", M1
        #print "rank of the matrix = ", M1.rank()



    #print " vec_len = ", vec_len
    #print M1


    ## Organize these vectors by length (and also introduce their negatives)
    max_len = vec_len/2
    vector_list_by_length = [[]  for _ in range(max_len + 1)]
    for v in vector_list:
        l = self(v)
        vector_list_by_length[l].append(v)
        vector_list_by_length[l].append(vector([-x  for x in v]))


    ## Make a matrix from the column vectors (in order of ascending length).
    sorted_list = []
    for i in range(len(vector_list_by_length)):
        for v in vector_list_by_length[i]:
            sorted_list.append(v)
    sorted_matrix = Matrix(sorted_list).transpose()


    ## Determine a basis of vectors of minimal length
    pivots = sorted_matrix.pivots()
    basis = [sorted_matrix.column(i) for i in pivots]
    pivot_lengths = [self(v)  for v in basis]


    ## DIAGNOSTIC
    #print "basis = ", basis
    #print "pivot_lengths = ", pivot_lengths


    ## Cache the result
    self.__basis_of_short_vectors = basis
    self.__basis_of_short_vectors_lenghts = pivot_lengths


    ## Return the appropriate result
    if show_lengths:
        return basis, pivot_lengths
    else:
        return basis




def short_vector_list_up_to_length(self, len_bound, up_to_sign_flag=False):
    """
    Return a list of lists of short vectors v, sorted by length, with
    Q(v) < len_bound.  The list in output[i] indexes all vectors of
    length i.  If the up_to_sign_flag is set to True, then only one of
    the vectors of the pair [v, -v] is listed.

    Note:  This processes the PARI/GP output to always give elements of type ZZ.

    OUTPUT:
        a list of lists of vectors.

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.short_vector_list_up_to_length(3)
        [[(0, 0, 0, 0)], [(1, 0, 0, 0), [-1, 0, 0, 0]], []]
        sage: Q.short_vector_list_up_to_length(4)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0), [-1, 0, 0, 0]],
         [],
         [(0, 1, 0, 0), [0, -1, 0, 0]]]
        sage: Q.short_vector_list_up_to_length(5)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0), [-1, 0, 0, 0]],
         [],
         [(0, 1, 0, 0), [0, -1, 0, 0]],
         [(1, 1, 0, 0),
         [-1, -1, 0, 0],
         (-1, 1, 0, 0),
         [1, -1, 0, 0],
         (2, 0, 0, 0),
         [-2, 0, 0, 0]]]
        sage: Q.short_vector_list_up_to_length(5, True)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0)],
         [],
         [(0, 1, 0, 0)],
         [(1, 1, 0, 0), (-1, 1, 0, 0), (2, 0, 0, 0)]]

    """
    ## Set an upper bound for the number of vectors to consider
    Max_number_of_vectors = 10000

    ## Generate a PARI matrix string for the associated Hessian matrix
    M_str = str(gp(self.matrix()))

    ## Generate the short vectors
    gp_mat = gp.qfminim(M_str, 2*len_bound-2, Max_number_of_vectors)[3].mattranspose()
    number_of_rows = gp_mat.matsize()[1]
    gp_mat_vector_list = [vector([ZZ(x)  for x in gp_mat[i+1,]])  for i in range(number_of_rows)]

    ## Sort the vectors into lists by their length
    vec_list = [[]  for i in range(len_bound)]
    for v in gp_mat_vector_list:
        vec_list[self(v)].append(v)
        if up_to_sign_flag == False:
            vec_list[self(v)].append([-x  for x in v])

    ## Add the zero vector by hand
    zero_vec = vector([ZZ(0)  for _ in range(self.dim())])
    vec_list[0].append(zero_vec)

    ## Return the sorted list
    return vec_list



def short_primitive_vector_list_up_to_length(self, len_bound, up_to_sign_flag=False):
    """
    Return a list of lists of short primitive vectors v, sorted by length, with
    Q(v) < len_bound.  The list in output[i] indexes all vectors of
    length i.  If the up_to_sign_flag is set to True, then only one of
    the vectors of the pair [v, -v] is listed.

    Note:  This processes the PARI/GP output to always give elements of type ZZ.

    OUTPUT:
        a list of lists of vectors.

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.short_vector_list_up_to_length(5, True)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0)],
         [],
         [(0, 1, 0, 0)],
         [(1, 1, 0, 0), (-1, 1, 0, 0), (2, 0, 0, 0)]]
        sage: Q.short_primitive_vector_list_up_to_length(5, True)
        [[], [(1, 0, 0, 0)], [], [(0, 1, 0, 0)], [(1, 1, 0, 0), (-1, 1, 0, 0)]]


    """
    ## Get a list of short vectors
    full_vec_list = self.short_vector_list_up_to_length(len_bound, up_to_sign_flag)

    ## Make a new list of the primitive vectors
    prim_vec_list = [[v  for v in L  if GCD(list(v)) == 1]   for L in full_vec_list]                 ## TO DO:  Arrange for GCD to take a vector argument!

    ## Return the list of primitive vectors
    return prim_vec_list




## ----------------------------------------------------------------------------------------------------


def automorphisms(self):
    """
    Return a list of the automorphisms of the quadratic form.

    sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
    sage: Q.number_of_automorphisms()
    48
    sage: 2^3 * factorial(3)
    48

    sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
    sage: Q.number_of_automorphisms()
    384
    sage: 2^4 * factorial(4)
    384

    sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
    sage: Q.number_of_automorphisms()
    16

    #sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1])
    #sage: time Q.number_of_automorphisms()
    #CPU times: user 119.44 s, sys: 0.87 s, total: 120.31 s
    #Wall time: 121.80
    #3840
    #sage: 2^5 * factorial(5)
    #3840
    """
    ## Check for a cached value
    try:
        return self.__automorphisms
    except:
        pass


    ## Find a basis of short vectors, and their lengths
    basis, pivot_lengths = self.basis_of_short_vectors(show_lengths=True)

    ## List the relevant vectors by length
    max_len = max(pivot_lengths)
    vector_list_by_length = self.short_primitive_vector_list_up_to_length(max_len + 1)


    ## Make the matrix A:e_i |--> v_i to our new basis.
    A = Matrix(basis).transpose()
    #A1 = A.inverse() * A.det()
    #Q1 = A1.transpose() * self.matrix() * A1       ## This is the matrix of Q
    #Q = self.matrix() * A.det()**2
    Q2 = A.transpose() * self.matrix() * A       ## This is the matrix of Q in the new basis
    Q3 = self.matrix()


    ## Determine all automorphisms
    n = self.dim()
    Auto_list = []
    #ct = 0

    ## DIAGNOSTIC
    #print "n = " + str(n)
    #print "pivot_lengths = " + str(pivot_lengths)
    #print "vector_list_by_length = " + str(vector_list_by_length)
    #print "length of vector_list_by_length = " + str(len(vector_list_by_length))

    for index_vec in mrange([len(vector_list_by_length[pivot_lengths[i]])  for i in range(n)]):
        M = Matrix([vector_list_by_length[pivot_lengths[i]][index_vec[i]]   for i in range(n)]).transpose()
        #Q1 = self.matrix()
        #if self(M) == self:
        #ct += 1
        #print "ct = ", ct, "   M = "
        #print M
        #print
        if M.transpose() * Q3 * M == Q2:       ## THIS DOES THE SAME THING! =(
            Auto_list.append(M)


    ## Cache the answer and return the list
    self.__automorphisms = Auto_list
    self.__number_of_automorphisms = len(Auto_list)
    return Auto_list




def number_of_automorphisms(self, recompute=False):
    """
    Return a list of the number of automorphisms (of det 1 and -1) of
    the quadratic form.

    If recompute is True, then we will recompute the cached result.

    OUTPUT:
        an integer >= 2.

    EXAMPLES:
        sage: Q = QuadraticForm(ZZ, 3, [1, 0, 0, 1, 0, 1], unsafe_initialization=True, number_of_automorphisms=-1)
        sage: Q.list_external_initializations()
        ['number_of_automorphisms']
        sage: Q.number_of_automorphisms()
        -1
        sage: Q.number_of_automorphisms(recompute=True)
        48
        sage: Q.list_external_initializations()
        []
    """
    ## Try to use the cached version if we can
    if not recompute:
        try:
            #print "Using the cached number of automorphisms."
            #print "We stored", self.__number_of_automorphisms
            return self.__number_of_automorphisms
        except:
            pass

    ## Otherwise cache and return the result
    #print "Recomputing the number of automorphisms based on the list of automorphisms."
    #self.__number_of_automorphisms = len(self.automorphisms())                                     ## This is now deprecated.
    self.__number_of_automorphisms = self.number_of_automorphisms__souvigner()
    try:
        self._external_initialization_list.remove('number_of_automorphisms')
    except:
        pass  ## Do nothing if the removal fails, since it might not be in the list (causing an error)!
    return self.__number_of_automorphisms



def number_of_automorphisms__souvigner(self):
    """
    Uses the Souvigner code to compute the number of automorphisms.
    """
    ## Write an input text file
    F_filename = '/tmp/tmp_isom_input' + str(random()) + ".txt"
    F = open(F_filename, 'w')
    #F = tempfile.NamedTemporaryFile(prefix='tmp_auto_input', suffix=".txt")   ## This fails because the Souvigner code doesn't like random filenames (hyphens are bad...)!
    F.write("#1 \n")
    n = self.dim()
    F.write(str(n) + "x0 \n")      ## Use the lower-triangular form
    for i in range(n):
        for j in range(i+1):
            if i == j:
                F.write(str(2 * self[i,j]) + " ")
            else:
                F.write(str(self[i,j]) + " ")
        F.write("\n")
    F.flush()
    #print "Input filename = ", F.name
    #os.system("less " + F.name)

    ## Call the Souvigner automorphism code
    souvigner_auto_path = SAGE_ROOT + "/local/bin/Souvigner_AUTO"                 ## FIX THIS LATER!!!
    G1 = tempfile.NamedTemporaryFile(prefix='tmp_auto_ouput', suffix=".txt")
    #print "Output filename = ", G1.name
    os.system(souvigner_auto_path + " '" +  F.name + "' > '" + G1.name +"'")


    ## Read the output
    G2 = open(G1.name, 'r')
    for line in G2:
        if line.startswith("|Aut| = "):
            num_of_autos = sage_eval(line.lstrip("|Aut| = "))
            F.close()
            G1.close()
            G2.close()
            os.system("rm -f " + F_filename)
            #os.system("rm -f " + G1.name)
            return num_of_autos

    ## Raise and error if we're here:
    raise RuntimeError, "Oops! There is a problem..."



def set_number_of_automorphisms(self, num_autos):
    """
    Set the number of automorphisms to be the value given.  No error
    checking is performed, to this may lead to erroneous results.

    The fact that this result was set externally is recorded in the
    internal list of external initializations, accessible by the
    method list_external_initializations().

    Return a list of the number of
    automorphisms (of det 1 and -1) of the quadratic form.

    OUTPUT:
        None

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1, 1, 1])
        sage: Q.list_external_initializations()
        []
        sage: Q.set_number_of_automorphisms(-3)
        sage: Q.number_of_automorphisms()
        -3
        sage: Q.list_external_initializations()
        ['number_of_automorphisms']

    """
    self.__number_of_automorphisms = num_autos
    text = 'number_of_automorphisms'
    if not text in self._external_initialization_list:
        self._external_initialization_list.append(text)




def Nipp_automorphism_testing(self):
    """
    Testing the automorphism routine agains Nipp's Tables

        --- MOVE THIS ELSEWHERE!!! ---

    """
    for i in range(20):
        Q = QuadraticForm(ZZ, 4, Nipp[i][2])
        my_num = Q.number_of_automorphisms()
        nipp_num = Nipp.number_of_automorphisms(i)
        print "    i = " + str(i) + "  my_num = " + str(my_num) + "  nipp_num = " + str(nipp_num)


