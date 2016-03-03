"""
Equivalence Testing

AUTHORS:

- Anna Haensch (2014-12-01): added test for rational isometry
"""

from sage.arith.all import hilbert_symbol, prime_divisors, is_prime, valuation, GCD, legendre_symbol
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from quadratic_form import is_QuadraticForm


################################################################################
## Routines to test if two quadratic forms over ZZ are globally equivalent.   ##
## (For now, we require both forms to be positive definite.)                  ##
################################################################################

def is_globally_equivalent_to(self, other, return_matrix=False, check_theta_to_precision=None, check_local_equivalence=None):
    """
    Determines if the current quadratic form is equivalent to the
    given form over ZZ.  If ``return_matrix`` is True, then we return
    the transformation matrix `M` so that ``self(M) == other``.

    INPUT:

    - ``self``, ``other`` -- positive definite integral quadratic forms

    - ``return_matrix`` -- (boolean, default ``False``) return
      the transformation matrix instead of a boolean

    OUTPUT:

    - if ``return_matrix`` is ``False``: a boolean

    - if ``return_matrix`` is ``True``: either ``False`` or the
      transformation matrix

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: M = Matrix(ZZ, 4, 4, [1,2,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1])
        sage: Q1 = Q(M)
        sage: Q.is_globally_equivalent_to(Q1)
        True
        sage: MM = Q.is_globally_equivalent_to(Q1, return_matrix=True)
        sage: Q(MM) == Q1
        True

    ::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q2 = QuadraticForm(ZZ, 3, [2, 1, 2, 2, 1, 3])
        sage: Q3 = QuadraticForm(ZZ, 3, [8, 6, 5, 3, 4, 2])
        sage: Q1.is_globally_equivalent_to(Q2)
        False
        sage: Q1.is_globally_equivalent_to(Q2, return_matrix=True)
        False
        sage: Q1.is_globally_equivalent_to(Q3)
        True
        sage: M = Q1.is_globally_equivalent_to(Q3, True); M
        [-1 -1  0]
        [ 1  1  1]
        [-1  0  0]
        sage: Q1(M) == Q3
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1, -1])
        sage: Q.is_globally_equivalent_to(Q)
        Traceback (most recent call last):
        ...
        ValueError: not a definite form in QuadraticForm.is_globally_equivalent_to()

    ALGORITHM: this uses the PARI function ``qfisom()``, implementing
    an algorithm by Plesken and Souvignier.
    """
    if check_theta_to_precision is not None:
        from sage.misc.superseded import deprecation
        deprecation(19111, "The check_theta_to_precision argument is deprecated and ignored")
    if check_local_equivalence is not None:
        from sage.misc.superseded import deprecation
        deprecation(19111, "The check_local_equivalence argument is deprecated and ignored")

    ## Check that other is a QuadraticForm
    if not is_QuadraticForm(other):
        raise TypeError("you must compare two quadratic forms, but the argument is not a quadratic form")

    ## only for definite forms
    if not self.is_definite() or not other.is_definite():
        raise ValueError("not a definite form in QuadraticForm.is_globally_equivalent_to()")

    mat = other._pari_().qfisom(self)
    if not mat:
        return False

    if return_matrix:
        return mat.sage()
    else:
        return True


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
        sage: Q1.is_globally_equivalent_to(Q2)
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

def is_rationally_isometric(self, other):
    """
    Determines if two regular quadratic forms over a number field are isometric.

    INPUT:

    a quadratic form

    OUTPUT:

    boolean

    EXAMPLES::

        sage: V=DiagonalQuadraticForm(QQ,[1,1,2])
        sage: W=DiagonalQuadraticForm(QQ,[2,2,2])
        sage: V.is_rationally_isometric(W)
        True

    ::

        sage: K.<a>=NumberField(x^2-3)
        sage: V=QuadraticForm(K,4,[1,0,0,0,2*a,0,0,a,0,2]);V
        Quadratic form in 4 variables over Number Field in a with defining polynomial x^2 - 3 with coefficients:
        [ 1 0 0 0 ]
        [ * 2*a 0 0 ]
        [ * * a 0 ]
        [ * * * 2 ]
        sage: W=QuadraticForm(K,4,[1,2*a,4,6,3,10,2,1,2,5]);W
        Quadratic form in 4 variables over Number Field in a with defining polynomial x^2 - 3 with coefficients:
        [ 1 2*a 4 6 ]
        [ * 3 10 2 ]
        [ * * 1 2 ]
        [ * * * 5 ]
        sage: V.is_rationally_isometric(W)
        False

    ::

        sage: K.<a>=NumberField(x^4+2*x+6)
        sage: V=DiagonalQuadraticForm(K,[a,2,3,2,1]);V
        Quadratic form in 5 variables over Number Field in a with defining polynomial x^4 + 2*x + 6 with coefficients:
        [ a 0 0 0 0 ]
        [ * 2 0 0 0 ]
        [ * * 3 0 0 ]
        [ * * * 2 0 ]
        [ * * * * 1 ]
        sage: W=DiagonalQuadraticForm(K,[a,a,a,2,1]);W
        Quadratic form in 5 variables over Number Field in a with defining polynomial x^4 + 2*x + 6 with   coefficients:
        [ a 0 0 0 0 ]
        [ * a 0 0 0 ]
        [ * * a 0 0 ]
        [ * * * 2 0 ]
        [ * * * * 1 ]
        sage: V.is_rationally_isometric(W)
        False

    ::

        sage: K.<a>=NumberField(x^2-3)
        sage: V=DiagonalQuadraticForm(K,[-1,a,-2*a])
        sage: W=DiagonalQuadraticForm(K,[-1,-a,2*a])
        sage: V.is_rationally_isometric(W)
        True

    TESTS::

        sage: K.<a>=QuadraticField(3)
        sage: V=DiagonalQuadraticForm(K,[1,2])
        sage: W=DiagonalQuadraticForm(K,[1,0])
        sage: V.is_rationally_isometric(W)
        Traceback (most recent call last):
        ...
        NotImplementedError: This only tests regular forms

    Forms must have the same base ring otherwise a `TypeError` is raised::

        sage: K1.<a> = QuadraticField(5)
        sage: K2.<b> = QuadraticField(7)
        sage: V = DiagonalQuadraticForm(K1,[1,a])
        sage: W = DiagonalQuadraticForm(K2,[1,b])
        sage: V.is_rationally_isometric(W)
        Traceback (most recent call last):
        ...
        TypeError: forms must have the same base ring.

    Forms which have different dimension are not isometric::

        sage: W=DiagonalQuadraticForm(QQ,[1,2])
        sage: V=DiagonalQuadraticForm(QQ,[1,1,1])
        sage: V.is_rationally_isometric(W)
        False

    Forms whose determinants do not differ by a square in the base field are not isometric::

        sage: K.<a>=NumberField(x^2-3)
        sage: V=DiagonalQuadraticForm(K,[-1,a,-2*a])
        sage: W=DiagonalQuadraticForm(K,[-1,a,2*a])
        sage: V.is_rationally_isometric(W)
        False

    ::

        sage: K.<a> = NumberField(x^5 - x + 2, 'a')
        sage: Q = QuadraticForm(K,3,[a,1,0,-a**2,-a**3,-1])
        sage: m = Q.matrix()
        sage: for _ in range(5):
        ....:     t = random_matrix(ZZ,3,algorithm='unimodular')
        ....:     m2 = t*m*t.transpose()
        ....:     Q2 = QuadraticForm(K, 3, [m2[i,j] / (2 if i==j else 1)
        ....:                               for i in range(3) for j in range(i,3)])
        ....:     print Q.is_rationally_isometric(Q2)
        True
        True
        True
        True
        True

    """

    if self.Gram_det() == 0 or other.Gram_det() == 0:
        raise NotImplementedError("This only tests regular forms")

    if self.base_ring() != other.base_ring():
        raise TypeError("forms must have the same base ring.")

    if self.dim() != other.dim():
        return False

    if not ((self.Gram_det()*other.Gram_det()).is_square()):
        return False

    L1=self.Gram_det().support()
    L2=other.Gram_det().support()

    for p in set().union(L1,L2):
        if self.hasse_invariant(p) != other.hasse_invariant(p):
            return False

    if self.base_ring() == QQ:
        if self.signature() != other.signature():
            return False
    else:

        M = self.rational_diagonal_form().Gram_matrix_rational()
        N = other.rational_diagonal_form().Gram_matrix_rational()
        K = self.base_ring()

        Mentries = M.diagonal()
        Nentries = N.diagonal()

        for emb in K.real_embeddings():

            Mpos=0
            for x in Mentries:
                Mpos+= emb(x) >= 0

            Npos=0
            for x in Nentries:
                Npos+= emb(x) >= 0

            if Npos != Mpos:
                return False

    return True



def isometry(self, other):
    r"""
    Given two rationally equivalent quadratic forms, computes a 
    transition matrix mapping from one to the other.

    INPUT:

    - ``self`` -- a quadratic form
    - ``other`` -- a quadratic form

    OUTPUT:
    
    - A matrix ``T`` representing the isometry transformation, such that if 
      ``QM`` is the gram matrix of ``self`` and ``FM`` is the gram matrix of 
      ``other``, then ``QM == T.transpose() * FM * T`` yields ``True``.

    EXAMPLES::
    
        sage: Q = DiagonalQuadraticForm(QQ, [1, 1, 2])
        sage: F = DiagonalQuadraticForm(QQ, [2, 2, 2])
        sage: T = Q.isometry(F); T
        [-1/2 -1/2    0]
        [   0    0   -1]
        [-1/2  1/2    0]
        sage: Q.Gram_matrix() == T.transpose() * F.Gram_matrix() * T
        True
    
        sage: T = F.isometry(Q); T
        [-1/3 -4/3 -1/3]
        [  -1    0    1]
        [ 2/3 -1/3  2/3]
        sage: F.Gram_matrix() == T.T * Q.Gram_matrix() * T
        True

    ::

        sage: L = QuadraticForm(QQ, 3, [2, 2, 0, 2, 2, 5])
        sage: M = QuadraticForm(QQ, 3, [2, 2, 0, 3, 2, 3])
        sage: L.isometry(M)
        Traceback (most recent call last):
        ...
        ArithmeticError: Quadratic form in 3 variables over Rational Field with coefficients: 
        [ 2 2 0 ]
        [ * 2 2 ]
        [ * * 5 ] is not rationally isometric to Quadratic form in 3 variables over Rational Field with coefficients: 
        [ 2 2 0 ]
        [ * 3 2 ]
        [ * * 3 ]

    ::
    
        sage: A = DiagonalQuadraticForm(QQ, [1, 5])
        sage: B = QuadraticForm(QQ, 2, [1, 12, 81])
        sage: T = A.isometry(B); T
        [  1  -2]
        [  0 1/3]
        sage: A.Gram_matrix() == T.T * B.Gram_matrix() * T
        True

    ::
    
        sage: C = DiagonalQuadraticForm(QQ, [1, 5, 9])
        sage: D = DiagonalQuadraticForm(QQ, [6, 30, 1])
        sage: T = C.isometry(D); T
        [ 7/18  5/18     0]
        [-1/18  7/18     0]
        [    0     0     3]
        sage: C.Gram_matrix() == T.T * D.Gram_matrix() * T
        True

    ::
    
        sage: E = DiagonalQuadraticForm(QQ, [1, 1])
        sage: F = QuadraticForm(QQ, 2, [17, 94, 130])
        sage: T = F.isometry(E); T
        [     -4 -189/17]
        [     -1  -43/17]
        sage: F.Gram_matrix() == T.T * E.Gram_matrix() * T
        True
    """
    if not is_QuadraticForm(other):
        raise TypeError("First argument must be a Quadratic Form.")

    if not self.is_rationally_isometric(other):
        raise ArithmeticError("{0} is not rationally isometric to {1}".format(self, other))

    Q, F = self, other
    n = Q.dim()

    # If either form is not diagonal, diagonalize it.
    Q, q_diagonal_transform = self.rational_diagonal_form(True)
    F, f_diagonal_transform = other.rational_diagonal_form(True)
        
    # Call the method that does all the work to compute the isometry.
    transform = _diagonal_isometry(Q,F)
    
    return f_diagonal_transform * transform * q_diagonal_transform.inverse()



def _diagonal_isometry(V, W):
    import copy
    from sage.functions.other import sqrt 
    from quadratic_form import DiagonalQuadraticForm
    from sage.matrix.constructor import Matrix
    from sage.modules.free_module_element import vector
    
    Q, F = copy.deepcopy(V), copy.deepcopy(W)
    #QM, FM = Q.Gram_matrix(), F.Gram_matrix()
    n = Q.dim()

    change_of_basis_matrix = Matrix.identity(QQ, n)
    #column_set? = [false for i in range(n)]

    for i in range(n):
        print "i: {0}".format(i)
        print "F = {0}".format(F)
        print "Q = {0}".format(Q)

        if Q.Gram_matrix()[0][0] != F.Gram_matrix()[0][0]:

            w = F.solve(Q.Gram_matrix()[0][0])
            print "w: {0}".format(w)
            change_of_basis_matrix.set_column(i, w)
            print "matrix (before ortho):\n{0}".format(change_of_basis_matrix)

            FM = _compute_gram_matrix_from_basis(W, change_of_basis_matrix)
            print "gram:\n{0}".format(FM)

            change_of_basis_matrix = _gram_schmidt(change_of_basis_matrix, i, W.bilinear_map)
            print "matrix (after ortho):\n{0}".format(change_of_basis_matrix)
            FM = _compute_gram_matrix_from_basis(W, change_of_basis_matrix)
            print "gram:\n{0}".format(FM)
        
        F = DiagonalQuadraticForm(F.base_ring(), FM.diagonal())
        F = F.extract_variables(range(i+1, F.dim()))
        #Q = DiagonalQuadraticForm(Q.base_ring(), Q.Gram_matrix().diagonal())
        Q = Q.extract_variables(range(1, Q.dim()))

    return change_of_basis_matrix



def _compute_gram_matrix_from_basis(Q, basis):
    r"""
    Computes the gram matrix of a quadratic form with respect to the given basis.

    INPUT:

    - ``Q`` -- a quadratic form of rank n
    - ``basis`` -- an nxn matrix whose columns are treated as basis vectors

    OUTPUT:

    - A matrix representing the gram matrix of ``Q`` with respect to ``basis``
    
    EXAMPLES::

        sage: from sage.quadratic_forms.quadratic_form__equivalence_testing import _compute_gram_matrix_from_basis
        
        sage: Q = QuadraticForm(QQ, 3, [1, 2, 2, 1, 2, 0])
        sage: std_basis = matrix.identity(3)
        sage: QM = _compute_gram_matrix_from_basis(Q, std_basis)
        sage: Q.Gram_matrix() == QM
        True

    ::

        sage: Q = DiagonalQuadraticForm(QQ, [1, 2, 2])
        sage: basis = matrix([[1, 1, 0], [2, 1, 0], [0, 0, 1]])
        sage: QM = _compute_gram_matrix_from_basis(Q, basis); QM
        [9 5 0]
        [5 3 0]
        [0 0 2]
    """
    from sage.matrix.constructor import matrix, Matrix

    n = Q.dim()
    rows = [];
    for i in range(n):
        rows.append([Q.bilinear_map(basis.column(i), basis.column(j)) for j in range(n)])
    return Matrix(rows)



def _gram_schmidt(m, fixed_vector_index, inner_product):
    r"""
    Orthogonalizes a set of vectors, starting at a fixed vector, with respect to a given 
    inner product.

    INPUT:

    - ``m`` -- a square matrix whose columns represent vectors
    - ``fixed_vector_index`` -- any vectors preceding the vector (i.e. to its left) 
        at this index are not changed.
    - ``inner_product`` - a function that takes two vector arguments and returns a scalar,
        representing an inner product.

    OUTPUT:

    - A matrix consisting of orthogonal columns with respect to the given inner product
    
    EXAMPLES:: 

        sage: from sage.quadratic_forms.quadratic_form__equivalence_testing import _gram_schmidt
        sage: Q = QuadraticForm(QQ, 3, [1, 4, 6, 1, 2, 4]); Q
        Quadratic form in 3 variables over Rational Field with coefficients: 
        [ 1 4 6 ]
        [ * 1 2 ]
        [ * * 4 ]
        sage: QM = Q.Gram_matrix(); QM
        [1 2 3]
        [2 1 1]
        [3 1 4]
        sage: om = _gram_schmidt(QM, 0, Q.bilinear_map); om
        [     1 106/79 -50/51]
        [     2 -25/79 -70/51]
        [     3 -77/79  70/51]
        sage: v0 = om.column(0); v1 = om.column(1); v2 = om.column(2)
        sage: Q.bilinear_map(v0, v1) == 0
        True
        sage: Q.bilinear_map(v0, v2) == 0
        True
        sage: Q.bilinear_map(v1, v2) == 0
        True
    """
    from sage.matrix.constructor import column_matrix
      
    n = m.dimensions()[0]
    vectors = [0] * n
    
    for i in range(n):
        vectors[i] = m.column(i)
    for i in range(fixed_vector_index, n):
        for j in range(i+1, n):
            vectors[j] = vectors[j] - (inner_product(vectors[j], vectors[i]) / inner_product(vectors[i], vectors[i])) * vectors[i]

    return column_matrix(vectors)
            
