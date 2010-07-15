r"""
This module contains lattice related functions relevant in cryptography.

Feel free to add more functionality.

AUTHORS:
    Richard Lindner <rlindner@cdc.informatik.tu-darmstadt.de>
    Michael Schneider <mischnei@cdc.informatik.tu-darmstadt.de>
"""

def gen_lattice(type='modular', n=4, m=8, q=11, seed=None, \
                quotient=None, dual=False, ntl=False):
    r"""
    This function generates different types of integral lattice
    bases of row vectors relevant in cryptography.

    Randomness can be set either with ``seed'', or by using
    set_random_seed(...) in SAGE.

    INPUT:

    * ``type`` - one of the following strings
        * ``'modular'`` (default). A class of lattices for which asymptotic
        worst-case to average-case connections hold. For more refer to [A96].
        * ``'random'`` - Special case of modular (n=1). A dense class of
        lattice used for testing basis reduction algorithms proposed by
        Goldstein and Mayer [GM02].
        * ``'ideal'`` - Special case of modular. Allows for a more compact
        representation proposed by [LM06].
        * ``'cyclotomic'`` - Special case of ideal. Allows for efficient
        processing proposed by [LM06].
    * ``n`` - Determinant size, det(L) = q^n. For ideal lattices this is also
    the degree of the quotient polynomial.
    * ``m`` - Lattice dimension, L \subseteq Z^m.
    * ``q`` - Coefficent size, q*Z^m \subseteq L.
    * ``seed`` - Randomness seed.
    * ``quotient`` - For the type ideal, this determines the quotient
    polynomial. Ignored for all other types.
    * ``dual`` - Set this flag if you want a basis for q*dual(L), for example
    for Regev's LWE bases [R05].
    * ``ntl`` - Set this flag if you want the lattice basis in NTL readable
    format.

    OUTPUT: ``B`` a unique size-reduced triangular (lower=primal, upper=dual)
    basis of row vectors for the lattice in question.

    EXAMPLES:

    * Modular basis ::

        sage: sage.crypto.gen_lattice(m=10, seed=42)
        [11  0  0  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0  0  0]
        [ 2  4  3  5  1  0  0  0  0  0]
        [ 1 -5 -4  2  0  1  0  0  0  0]
        [-4  3 -1  1  0  0  1  0  0  0]
        [-2 -3 -4 -1  0  0  0  1  0  0]
        [-5 -5  3  3  0  0  0  0  1  0]
        [-4 -3  2 -5  0  0  0  0  0  1]

    * Random basis ::

        sage: sage.crypto.gen_lattice(type='random', n=1, m=10, q=11^4, \
        seed=42)
        [14641     0     0     0     0     0     0     0     0     0]
        [  431     1     0     0     0     0     0     0     0     0]
        [-4792     0     1     0     0     0     0     0     0     0]
        [ 1015     0     0     1     0     0     0     0     0     0]
        [-3086     0     0     0     1     0     0     0     0     0]
        [-5378     0     0     0     0     1     0     0     0     0]
        [ 4769     0     0     0     0     0     1     0     0     0]
        [-1159     0     0     0     0     0     0     1     0     0]
        [ 3082     0     0     0     0     0     0     0     1     0]
        [-4580     0     0     0     0     0     0     0     0     1]

    * Ideal bases with quotient x^n-1, m=2*n are NTRU bases ::

        sage: sage.crypto.gen_lattice(type='ideal', seed=42, quotient=x^4-1)
        [11  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0]
        [ 4 -2 -3 -3  1  0  0  0]
        [-3  4 -2 -3  0  1  0  0]
        [-3 -3  4 -2  0  0  1  0]
        [-2 -3 -3  4  0  0  0  1]

    * Cyclotomic bases with n=2^k are SWIFFT bases ::

        sage: sage.crypto.gen_lattice(type='cyclotomic', seed=42)
        [11  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0]
        [ 4 -2 -3 -3  1  0  0  0]
        [ 3  4 -2 -3  0  1  0  0]
        [ 3  3  4 -2  0  0  1  0]
        [ 2  3  3  4  0  0  0  1]

    * Dual modular bases are related to Regev's famous public-key encryption
    [R05] ::

        sage: sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)
        [ 1  0  0  0 -2 -1  4  2  5  4]
        [ 0  1  0  0 -4  5 -3  3  5  3]
        [ 0  0  1  0 -3  4  1  4 -3 -2]
        [ 0  0  0  1 -5 -2 -1  1 -3  5]
        [ 0  0  0  0 11  0  0  0  0  0]
        [ 0  0  0  0  0 11  0  0  0  0]
        [ 0  0  0  0  0  0 11  0  0  0]
        [ 0  0  0  0  0  0  0 11  0  0]
        [ 0  0  0  0  0  0  0  0 11  0]
        [ 0  0  0  0  0  0  0  0  0 11]

    * Relation of primal and dual bases ::

        sage: B_primal=sage.crypto.gen_lattice(m=10, q=11, seed=42)
        sage: B_dual=sage.crypto.gen_lattice(m=10, q=11, seed=42, dual=True)
        sage: transpose(B_primal)*B_dual == 11*identity_matrix(10)
        True

    REFERENCES:

    [A96] Mikl{\'o}s Ajtai.
    Generating hard instances of lattice problems (extended abstract).
    STOC, pp. 99--108, ACM, 1996.

    [GM02] Daniel Goldstein and Andrew Mayer.
    On the equidistribution of Hecke points.
    Forum Mathematicum, 15:2, pp. 165--189, De Gruyter, 2003.

    [LM06] Vadim Lyubashevsky and Daniele Micciancio.
    Generalized compact knapsacks are collision resistant.
    ICALP, pp. 144--155, Springer, 2006.

    [R05] Oded Regev.
    On lattices, learning with errors, random linear codes, and cryptography.
    STOC, pp. 84--93, ACM, 2005.
    """
    from sage.rings.finite_rings.integer_mod_ring \
        import IntegerModRing
    from sage.matrix.constructor import Matrix, \
        identity_matrix, block_matrix
    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import IntegerRing
    if seed != None:
        from sage.misc.randstate import set_random_seed
        set_random_seed(seed)

    if type == 'random':
        if n != 1: raise ValueError('random bases require n = 1')

    ZZ = IntegerRing()
    ZZ_q = IntegerModRing(q)
    A = identity_matrix(ZZ_q, n)

    if type == 'random' or type == 'modular':
        R = MatrixSpace(ZZ_q, m-n, n)
        A = A.stack(R.random_element())

    elif type == 'ideal':
        if quotient == None: raise \
            ValueError('ideal bases require a quotient polynomial')
        x = quotient.default_variable()
        if n != quotient.degree(x): raise \
            ValueError('ideal bases require n  = quotient.degree()')
        R = ZZ_q[x].quotient(quotient, x)
        for i in range(m//n):
            A = A.stack(R.random_element().matrix())

    elif type == 'cyclotomic':
        from sage.rings.arith import euler_phi
        from sage.misc.functional import cyclotomic_polynomial

        # we assume that n+1 <= min( euler_phi^{-1}(n) ) <= 2*n
        found = False
        for k in range(2*n,n,-1):
            if euler_phi(k) == n:
                found = True
                break
        if not found: raise \
            ValueError('cyclotomic bases require that n is an image of' + \
                       'Euler\'s totient function')

        R = ZZ_q['x'].quotient(cyclotomic_polynomial(k, 'x'), 'x')
        for i in range(m//n):
            A = A.stack(R.random_element().matrix())

    # switch from representatives 0,...,(q-1) to (1-q)/2,....,(q-1)/2
    def minrep(a):
        if abs(a-q) < abs(a): return a-q
        else: return a
    A_prime = A[n:m].lift().apply_map(minrep)

    if not dual:
        B = block_matrix([ZZ(q), ZZ.zero() , A_prime, ZZ.one() ], 2 , 2 , \
                         False)
    else:
        B = block_matrix([ZZ.one() , -A_prime.transpose(), ZZ.zero() , \
                         ZZ(q)], 2 , 2 , False)

    if not ntl:
        return B
    else:
        from sage.libs.ntl.ntl_ZZ import ntl_ZZ
        from sage.libs.ntl.ntl_mat_ZZ import ntl_mat_ZZ
        return ntl_mat_ZZ(m, m, map(ntl_ZZ, B.list()))
