
from copy import deepcopy

from sage.rings.real_mpfr import RealField
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.calculus.calculus import sqrt, floor, ceil

from sage.interfaces.gp import gp

from sage.modular.dims import sturm_bound


def theta_series(self, Max='sturm', var_str='q', safe_flag=True):
    """
    Compute the theta series as a power series in the variable given
    in var_str (which defaults to 'q'), up to the specified precision
    O(q^max).

    If no precision is specified, then it defaults to computing the
    precision specified by sturm_bound() + 1, which suffices to
    uniquely determine the cuspidal part of the theta series.

    This uses the PARI/GP function qfrep, wrapped by the
    theta_by_pari() method.  This caches teh result for future
    computations.

    WARNING: This may not be the correct default bound for
    odd-dimensional quadratic forms!!  CHECK THIS!!!

    The safe_flag allows us to select whether we want a copy of the
    output, or the original output.  It is only meningful when a
    vector is returned, otherwise a copy is automatically made in
    creating the power series.  By default safe_flag = True, so we
    return a copy of the cached information.  If this is set to False,
    then the routine is much faster but the return values are
    vulnerable to being corrupted by the user.

    """
    ## Sanity Check: Max is an integer or an allowed string:
    try:
        M = ZZ(Max)
    except:
        M = -1

    if (Max not in ['sturm']) and (not M >= 0):
        print Max
        raise TypeError, "Oops!  Max is not an integer >= 0 or an allowed string."

    if Max == 'sturm':
        return self.theta_by_pari(sturm_bound(self.level(), self.dim() / ZZ(2)) + 1, var_str, safe_flag)
    else:
        return self.theta_by_pari(M, var_str, safe_flag)



## -------------  Compute the theta function by using the PARI/GP routine qfrep  ------------

def theta_by_pari(self, Max, var_str='q', safe_flag=True):
    """
    Use PARI/GP to compute the theta function as a power series (or
    vector) up to the precision O(q^Max).  This also caches the result
    for future computations.

    If var_str = '', then we return a vector v where v[i] counts the
    number of vectors of length i.

    The safe_flag allows us to select whether we want a copy of the
    output, or the original output.  It is only meningful when a
    vector is returned, otherwise a copy is automatically made in
    creating the power series.  By default safe_flag = True, so we
    return a copy of the cached information.  If this is set to False,
    then the routine is much faster but the return values are
    vulnerable to being corrupted by the user.


    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Prec = 100
        sage: compute = Q.theta_by_pari(Prec, '')
        sage: exact = [1] + [8 * sum([d  for d in divisors(i)  if d % 4 != 0])  for i in range(1, Prec)]
        sage: compute == exact
        True

    INPUT:
        Max -- an integer >=0
        var_str -- a string

    OUTPUT:
        a power series or a vector
    """
    ## Try to use the cached result if it's enough precision
    try:
        if len(self.__theta_vec) >= Max:
            theta_vec = self.__theta_vec[:Max]
        else:
            raise RuntimeError, ""
    except:
        ## Generate a PARI matrix string for the associated Hessian matrix
        M_str = str(gp(self.matrix()))

        ## Compute the theta function as a vector
        gp_vec = list(gp.qfrep(M_str, 2*Max-2))

        ## Create the "halved" representation/theta vector
        theta_vec = [ZZ(1)] + [ZZ(2) * ZZ(gp_vec[x])  for x in range(1, len(gp_vec), 2)]

        ## Cache the theta vector
        self.__theta_vec = theta_vec


    ## Return the answer
    if var_str == '':
        if safe_flag:
            return deepcopy(theta_vec)         ## We must make a copy here to insure the integrity of the cached version!
        else:
            return theta_vec
    else:
        return PowerSeriesRing(ZZ, var_str)(theta_vec)



## -------------  Compute the theta function by using an explicit Cholesky decomposition ------------

##########################################################################
## Routines to compute the Fourier expansion of the theta function of Q ##
## (to a given precision) via a Cholesky decomposition over RR.         ##
##                                                                      ##
## The Cholesky code was taken from:                                    ##
## ~/Documents/290_Project/C/Ver13.2__3-5-2007/Matrix_mpz/Matrix_mpz.cc ##
##########################################################################

def cholesky_decomposition(self, bit_prec = 53):
    """
    Give the Cholesky decomposition of Q as a real matrix of precision
        bit_prec.

    RESTRICTIONS: Q must be given as a QuadraticForm defnined over ZZ,
        QQ, or some RealField().  If it is over some real field, then
        an error is raised if the precision given is not less than the
        defined precision of the real field defining the quadratic form!

    REFERENCE: From Cohen's "A Course in Computational Algebraic
        Number Theory" book, p 103.

    INPUT:
        bit_prec -- a natural number.

    OUTPUT:
        an upper triangular real matrix of precision bit_prec.


    TO DO: If we only care about working over the real double field
        (RDF), then we can use the cholesky() method present for
        square matrices over that.


    ##/////////////////////////////////////////////////////////////////////////////////////////////////
    ##/// Finds the Cholesky decomposition of a quadratic form -- as an upper-triangular matrix!
    ##/// (It's assumed to be global, hence twice the form it refers to.)  <-- Python revision asks:  Is this true?!? =|
    ##/////////////////////////////////////////////////////////////////////////////////////////////////


    EXAMPLES:


    """

    ## Check that the precision passed is allowed.
    if isinstance(self.base_ring(), RealField) and (self.base_ring().prec() < bit_prec):
        raise RuntimeError, "Oops! The precision requested is greater than that of the given quadratic form!"

    ## 1. Initialization
    n = self.dim()
    R = RealField(bit_prec)
    MS = MatrixSpace(R, n, n)
    Q = MS(R(0.5)) * MS(self.matrix())               ## Initialize the real symmetric matrix A with the matrix for Q(x) = x^t * A * x

    ## DIAGNOSTIC
    #print "After 1:  Q is \n" + str(Q)

    ## 2. Loop on i
    for i in range(n):
        for j in range(i+1, n):
            Q[j,i] = Q[i,j]             ## Is this line redudnant?
            Q[i,j] = Q[i,j] / Q[i,i]

        ## 3. Main Loop
        for k in range(i+1, n):
            for l in range(k, n):
                Q[k,l] = Q[k,l] - Q[k,i] * Q[i,l]

    ## 4. Zero out the strictly lower-triangular entries
    for i in range(n):
        for j in range(i-1):
            Q[i,j] = 0

    return Q


def theta_by_cholesky(self, q_prec):
    """
    Uses the real Choelesky decomposition to compute (the q-expansion
    of) the theta function of the quadratic form as a power series in
    q with terms correct up to the power q ^ q_prec.  (So its error is
    O(q ^ (q_prec + 1)).)

    REFERENCE: From Cohen's "A Course in Computational Algebraic
        Number Theory" book, p 102.


    EXAMPLES:
        ## Check the sum of 4 squares form against Jacobi's formula
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Theta = Q.theta_by_cholesky(10)
        sage: Theta
        1 + 8*q + 24*q^2 + 32*q^3 + 24*q^4 + 48*q^5 + 96*q^6 + 64*q^7 + 24*q^8 + 104*q^9 + 144*q^10
        sage: Expected =  [1] + [8*sum([d for d in divisors(n) if d%4 != 0])  for n in range(1,11)]
        sage: Expected
        [1, 8, 24, 32, 24, 48, 96, 64, 24, 104, 144]
        sage: Theta.list() == Expected
        True

        ## Check the form x^2 + 3y^2 + 5z^2 + 7w^2 represents everything except 2 and 22.
        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Theta = Q.theta_by_cholesky(50)
        sage: Theta_list = Theta.list()
        sage: [m  for m in range(len(Theta_list))  if Theta_list[m] == 0]
        [2, 22]

    """
    ## RAISE AN ERROR -- This routine is deprecated!
    #raise NotImplementedError, "This routine is deprecated.  Try theta_series(), which uses theta_by_pari()."


    n = self.dim()
    theta = [0 for i in range(q_prec+1)]
    PS = PowerSeriesRing(ZZ, 'q')

    bit_prec = 53                                       ## TO DO: Set this precision to reflect the appropriate roundoff
    Cholesky = self.cholesky_decomposition(bit_prec)     ## error estimate, to be confident through our desired q-precision.
    Q = Cholesky      ##  <----  REDUNDANT!!!
    R = RealField(bit_prec)
    half = R(0.5)



    ## 1. Initialize
    i = n - 1
    T = [R(0)  for j in range(n)]
    U = [R(0)  for j in range(n)]
    T[i] = R(q_prec)
    U[i] = 0
    L = [0 for j in range (n)]
    x = [0 for j in range (n)]


    ## 2. Compute bounds
    #Z = sqrt(T[i] / Q[i,i])      ## IMPORTANT NOTE: sqrt preserves the precision of the real number it's given... which is not so good... =|
    #L[i] = floor(Z - U[i])       ## Note: This is a SAGE Integer
    #x[i] = ceil(-Z - U[i]) - 1   ## Note: This is a SAGE Integer too


    done_flag = False
    from_step4_flag = False
    from_step3_flag = True        ## We start by pretending this, since then we get to run through 2 and 3a once. =)

    #double Q_val_double;
    #unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...



    ## Big loop which runs through all vectors
    while (done_flag == False):

        ## Loop through until we get to i=1 (so we defined a vector x)
        while from_step3_flag or from_step4_flag:              ## IMPORTANT WARNING:  This replaces a do...while loop, so it may have to be adjusted!

            ## Go to directly to step 3 if we're coming from step 4, otherwise perform step 2.
            if from_step4_flag:
                from_step4_flag = False
            else:
                ## 2. Compute bounds
                from_step3_flag = False
                Z = sqrt(T[i] / Q[i,i])
                L[i] = floor(Z - U[i])
                x[i] = ceil(-Z - U[i]) - 1



            ## 3a. Main loop

            ## DIAGNOSTIC
            #print
            #print "  L = ", L
            #print "  x = ", x

            x[i] += 1
            while (x[i] > L[i]):

                ## DIAGNOSTIC
                #print "  x = ", x

                i += 1
                x[i] += 1


            ## 3b. Main loop
            if (i > 0):
                from_step3_flag = True

                ## DIAGNOSTIC
                #print " i = " + str(i)
                #print " T[i] = " + str(T[i])
                #print " Q[i,i] = " + str(Q[i,i])
                #print " x[i] = " + str(x[i])
                #print " U[i] = " + str(U[i])
                #print " x[i] + U[i] = " + str(x[i] + U[i])
                #print " T[i-1] = " + str(T[i-1])

                T[i-1] = T[i] - Q[i,i] * (x[i] + U[i]) * (x[i] + U[i])

                # DIAGNOSTIC
                #print " T[i-1] = " + str(T[i-1])
                #print

                i += - 1
                U[i] = 0
                for j in range(i+1, n):
                    U[i] += Q[i,j] * x[j]



        ## 4. Solution found (This happens when i=0)
        from_step4_flag = True
        Q_val_double = q_prec - T[0] + Q[0,0] * (x[0] + U[0]) * (x[0] + U[0])
        Q_val = floor(Q_val_double + half)        ## Note: This rounds the value up, since the "round" function returns a float, but floor returns integer.



        ## DIAGNOSTIC
        #print " Q_val_double = ",  Q_val_double
        #print " Q_val = ",  Q_val
        #raise RuntimeError


        ## OPTIONAL SAFETY CHECK:
        eps = 0.000000001
        if (abs(Q_val_double - Q_val) > eps):
            raise RuntimeError, "Oh No!  We have a problem with the floating point precision... \n" \
                + " Q_val_double = " + str(Q_val_double) + "\n" \
                + " Q_val = " + str(Q_val) + "\n" \
                + " x = " + str(x) + "\n"


        ## DIAGNOSTIC
        #print " The float value is " + str(Q_val_double)
        #print " The associated long value is " + str(Q_val)
        #print

        if (Q_val <= q_prec):
            theta[Q_val] += 2

        ## 5. Check if x = 0, for exit condition. =)
        done_flag = True
        for j in range(n):
            if (x[j] != 0):
                done_flag = False


    ## Set the value: theta[0] = 1
    theta[0] = 1

    ## DIAGNOSTIC
    #print "Leaving ComputeTheta \n"


    ## Return the series, truncated to the desired q-precision
    return PS(theta)

