#*****************************************************************************
#  Copyright (C) 2018 Edgar Costa <edgarcosta@math.dartmouth.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.functions.log import log
from sage.functions.other import ceil


def power_sum(coefficients, length = None):
    """
    Returns the power sum coefficients of a polynomial.
    """
    # f = \prod (x - xi)
    # returns [pk] where pk = \sum x_i ^k

    # e_k = sum x_{i_1} ... x_{i_k}
    # where 1 <= i_1 < i_2 < ... < i_k <= n

    if length == None:
        length = len(coefficients);
    if coefficients == []:
        return [None] + [0]*(length - 1);
    e = coefficients[:];
    e.reverse()
    assert e[0] == 1;
    for i in range(1, len(e), 2):
        e[i] *= -1;

    e.extend([0]*(length - len(e)));

    p = [0] * length;
    p[0] = None;
    p[1] = e[1];
    for k in range(2, length):
        p[k] = (-1)**(k-1) * k * e[k];
        p[k] += sum([(-1)**(k - 1 + i) * e[k - i] * p[i] for i in range(1,k)]);

    assert len(p) == length

    return p;




def charpoly_frobenius(frob_matrix, charpoly_prec, p, weight, a = 1, known_factor = [1]):
    """
    Returns the characteristic polynomial of the given Frobenius matrix.

    INPUT:

    - ``frob_matrix`` -- a matrix representing the p-power Frobenius lift to Z_q up to some precision

    - ``charpoly_prec`` -- a vector ai, such that, frob_matrix.change_ring(ZZ).charpoly()[i] will be correct mod p^ai, this can be easily deduced from the hodge numbers and knowing the q-adic precision of frob_matrix

    - ``p`` -- prime p

    - ``weight`` -- weight of the motive

    - ``a`` -- q = q^a

    - ``known_factor`` -- If there is a known factor, this are its coefficients

    OUTPUT: a list of integers corresponding to the characteristic polynomial of the Frobenius action

    Examples::

        sage: from sage.schemes.cycliccovers.charpoly_frobenius import charpoly_frobenius
        sage: M = Matrix([[O(17), 8 + O(17)], [O(17), 15 + O(17)]])
        sage: charpoly_frobenius(M, [2, 1, 1], 17, 1, 1)
        [17, 2, 1]

        sage: R = Zq(17 ** 2 , names=('a',));
        sage: M = Matrix(R, [[8*17 + 16*17**2 + O(17**3), 8 + 11*17 + O(17**2)], [7*17**2 + O(17**3), 15 + 8*17 + O(17**2)]])
        sage: charpoly_frobenius(M, [3, 2, 2], 17, 1, 2)
        [289, 30, 1]

        sage: M = Matrix([[8*31 + 8*31**2 + O(31**3), O(31**3), O(31**3), O(31**3)], [O(31**3), 23*31 + 22*31**2 + O(31**3), O(31**3), O(31**3)], [O(31**3), O(31**3), 27 + 7*31 + O(31**3), O(31**3)], [O(31**3), O(31**3), O(31**3), 4 + 23*31 + O(31**3)]])
        sage: charpoly_frobenius(M, [4, 3, 2, 2, 2], 31, 1, 1)
        [961, 0, 46, 0, 1]

        sage: M = Matrix([(4*43^2 + O(43^3), 17*43 + 11*43^2 + O(43^3), O(43^3), O(43^3), 17 + 37*43 + O(43^3), O(43^3)),\
        ....:  (30*43 + 23*43^2 + O(43^3), 5*43 + O(43^3), O(43^3), O(43^3), 3 + 38*43 + O(43^3), O(43^3)),\
        ....:  (O(43^3), O(43^3), 9*43 + 32*43^2 + O(43^3), 13 + 25*43 + O(43^3), O(43^3), 17 + 18*43 + O(43^3)),\
        ....:  (O(43^3), O(43^3), 22*43 + 25*43^2 + O(43^3), 11 + 24*43 + O(43^3), O(43^3), 36 + 5*43 + O(43^3)),\
        ....:  (42*43 + 15*43^2 + O(43^3), 22*43 + 8*43^2 + O(43^3), O(43^3), O(43^3), 29 + 4*43 + O(43^3), O(43^3)),\
        ....:  (O(43^3), O(43^3), 6*43 + 19*43^2 + O(43^3), 8 + 24*43 + O(43^3), O(43^3), 31 + 42*43 + O(43^3))])
        sage: charpoly_frobenius(M, [5, 4, 3, 2, 2, 2, 2], 43, 1, 1)
            [79507, 27735, 6579, 1258, 153, 15, 1]

        sage: M = Matrix([(1 + O(4999), O(4999), 0, 0), \
        ....:  (O(4999), 4860 + O(4999), 0, 0), \
        ....:  (0, 0, O(4999), O(4999)), \
        ....:  (0, 0, O(4999), 1 + O(4999))])
        sage: charpoly_frobenius(M, [2, 1, 1], 4999, 1, 1, [1, -2 ,1 ])
        [4999, 139, 1]

    """
    assert known_factor[-1] == 1;
    if a > 1:
        F = frob_matrix
        sigmaF = F;
        for _ in range(a - 1):
            sigmaF = matrix(F.base_ring(), [[elt.frobenius() for elt in row] for row in sigmaF.rows()])
            F = F * sigmaF
        F = F.change_ring(ZZ)
    else:
        F = frob_matrix.change_ring(ZZ)
    cp = F.charpoly().list();
    assert len(charpoly_prec) == len(cp) - (len(known_factor) - 1)
    assert cp[-1] == 1;

    # reduce cp mod prec
    degree = len(charpoly_prec) - 1;
    halfdegree = ceil(degree/2) + 1;
    mod = [0] * (degree + 1)
    for i in range(len(charpoly_prec)):
        mod[-i] = p**charpoly_prec[ -i]
        cp[-i] = cp[-i] % mod[-i]

    # figure out the sign
    if weight % 2 == 1:
        # for odd weight the sign is always 1
        # it's the charpoly of a USp matrix
        # and charpoly of a symplectic matrix is reciprocal
        sign = 1;
    else:
        # For the moment I will not worry about this case
        if known_factor != [1]:
            raise NotImplementedError();
        for i in range(degree/2):
            p_power = p ** min( charpoly_prec[i], charpoly_prec[degree - i] + (a*(degree - 2*i)*weight/2));
            # Note: degree*weight = 0 mod 2
            if cp[i] % p_power != 0 and cp[degree-i] % p_power != 0:
                if 0 == (cp[i] + cp[degree - i] * p**(a*(degree-2*i)*weight/2)) %  p_power:
                    sign = -1;
                else:
                    sign = 1;
                assert 0 == (-sign*cp[i] + cp[degree - i] * p**(a*(degree-2*i)*weight/2)) %  p_power;
                break;
    cp[0] = sign * p**(a*degree*weight/2)

    #calculate the i-th power sum of the roots and correct cp allong the way
    e = cp[-halfdegree:];
    e.reverse();
    for k in range(halfdegree):
        if k % 2 != 0:
            e[k] = -e[k] % mod[degree - k];
        #e[k] = cp[degree - k] if (k%2 ==0) else -cp[degree - k]
        if k > 0:
            # verify if p^charpoly_prec[degree - k] > 2*degree/k * q^(w*k/2)
            assert log(k)/log(p) + charpoly_prec[degree - k] > log(2*degree)/log(p) + a*0.5*weight*k, "log(k)/log(p) + charpoly_prec[degree - k] <= log(2*degree)/log(p) + a*0.5*weight*k, k = %d" % k

    fix_e = known_factor[:];
    fix_e.reverse();
    if len(fix_e) < halfdegree:
        fix_e.extend([0]*(halfdegree - len(fix_e)));
    for i in range(halfdegree):
        if i%2 != 0:
            fix_e[i] *= -1;

    #e[k] = \sum x_{i_1} x_{i_2} ... x_{i_k} # where x_* are eigenvalues
    # and i_1 < i_2 ... < i_k

    #s[k] = \sum x_i ^k for k>0
    s = [None]*(halfdegree)
    res = [None]*len(charpoly_prec);
    res[0] = sign * p**(a*degree*weight/2);
    res[-1] = 1;
    e[1] -= fix_e[1];
    e[1] = e[1] % mod[degree - 1];
    for k in range(1, halfdegree):
        # assume that s[i] and e[i] are correct for i < k
        # e[k] correct modulo mod[degree - k]
        # S = sum (-1)^i e[k-i] * s[i]
        # s[k] = (-1)^(k-1) (k*e[k] + S) ==> (-1)^(k-1) s[k] - S = k*e[k]
        S = sum( (-1)**i * e[k-i] * s[i] for i in range(1,k))
        s[k] = (-1)**(k - 1) * (S + k*e[k]);
        #hence s[k] is correct modulo k*mod[degree - k]
        localmod = k*mod[ degree -k]
        #s[k] +=   (-1)**k * fix_power_sum[k];
        s[k] = s[k] % localmod

        # |x_i| = p^(w*0.5)
        # => s[k] <= degree*p^(a*w*k*0.5)
        # recall, 2*degree*p^(a*w*k*0.5) /k < mod[degree - k]
        if s[k] > degree*p**(a*weight*k*0.5):
            s[k] = -(-s[k] % localmod);

        #now correct e[k] with:
        # (-1)^(k-1) s[k] - S = k*e[k]
        e[k] = (-S + (-1)**(k - 1) * s[k])//k;
        assert (-S + (-1)**(k - 1) * s[k])%k == 0
        res[degree - k] = e[k] if k % 2 == 0 else -e[k]
        res[k] = sign*res[degree - k]*p**(a*(degree-2*k)*weight/2)
        # fix e[k + 1]
        if k + 1 < halfdegree:
            e[k + 1] -= sum([fix_e[k + 1 - i] * e[i] for i in range(k + 1)]);
            e[k + 1] = e[k + 1] % mod[degree-(k+1)];
    return res
