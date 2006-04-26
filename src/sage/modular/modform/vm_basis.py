from eis_series import eisenstein_series
from sage.rings.all import QQ, Integer
from sage.modular.dims import dimension_cusp_forms_gamma0
from sage.matrix.all import MatrixSpace

def victor_miller_basis(k, prec=10, cusp_only=False, var='q'):
    r"""
    Compute and return the Victor-Miller basis for
    modular forms of weight k and level 1 to precision
    O(q^prec).  if \code{cusp_only} is True, return
    only a basis for the cuspidal subspace.

    INPUT:
        k -- an integer
        prec -- (default: 10) a positive integer
        cusp_only -- bool (default: False)
        var -- string (default: 'q'

    OUTPUT:
        list -- entries a power series in var defined over Q.

    EXAMPLES:
        sage: victor_miller_basis(1, 6)
        []
        sage: victor_miller_basis(0, 6)
        [1 + O(q^6)]
        sage: victor_miller_basis(2, 6)
        []
        sage: victor_miller_basis(4, 6)
        [1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)]

        sage: victor_miller_basis(6, 6, var='w')
        [1 - 504*w - 16632*w^2 - 122976*w^3 - 532728*w^4 - 1575504*w^5 + O(w^6)]

        sage: victor_miller_basis(6, 6)
        [1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
        sage: victor_miller_basis(12, 6)

        [1 + 196560*q^2 + 16773120*q^3 + 398034000*q^4 + 4629381120*q^5 + O(q^6),
         q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)]
        sage: victor_miller_basis(12, 6, cusp_only=True)
        [1 - 1008*q + 220752*q^2 + 16519104*q^3 + 399517776*q^4 + 4624512480*q^5 + O(q^6)]
        sage: victor_miller_basis(24, 6, cusp_only=True)

        [1 - 622944*q^2 + 82317312*q^3 + 38334552480*q^4 + 6618389299200*q^5 + O(q^6),
         q - 1032*q^2 + 245196*q^3 + 10965568*q^4 + 60177390*q^5 + O(q^6)]
        sage: victor_miller_basis(24, 6)

        [1 + 52416000*q^3 + 39007332000*q^4 + 6609020221440*q^5 + O(q^6),
         q + 195660*q^3 + 12080128*q^4 + 44656110*q^5 + O(q^6),
         q^2 - 48*q^3 + 1080*q^4 - 15040*q^5 + O(q^6)]
        sage: victor_miller_basis(32, 6)

        [1 + 2611200*q^3 + 19524758400*q^4 + 19715347537920*q^5 + O(q^6),
         q + 50220*q^3 + 87866368*q^4 + 18647219790*q^5 + O(q^6),
         q^2 + 432*q^3 + 39960*q^4 - 1418560*q^5 + O(q^6)]

    """
    k = Integer(k)

    R = QQ[[var]]
    if k == 0:
        return [R(1, prec)]
    elif k%2 == 1 or k < 4:
        return []

    kk = k % 12
    if kk == 2:
        kk += 12
    b = None
    for a in range(15):
        c = kk - 4*a
        if c % 6 == 0:
            b = c // 6
            break
    assert not (b is None), "bug in VM basis"

    F4 = 240*eisenstein_series(4, prec)
    F6 = -504*eisenstein_series(6, prec)
    if var != 'q':
        F4 = R(F4)
        F6 = R(F6)
    Delta = (F4**3 - F6**2)/R(1728,prec)
    d = dimension_cusp_forms_gamma0(1, k)
    g = F6**(2*d + b) * F4**a
    m = Delta / (F6*F6)
    G = []
    for j in range(d):
        G.append(g)
        if j < d-1:
            g *= m

    if not cusp_only:
        G.insert(0, R(eisenstein_series(k, prec)))

    M = MatrixSpace(QQ, len(G), prec)
    # we have to slice since precision in products can increase.
    e = [list(g)[:int(prec)] for g in G]
    A = M(sum(e, []))
    # this is still provably correct -- the guess is still proven right.
    # it's just that naive guess based on coefficients is way too big.
    E = A.echelon_form(height_guess=10**(k))
    return [R(list(v), prec) for v in E.rows()]


def delta_qexp(prec=10, var='q'):
    """
    Return the q-expansion of Delta.
    """
    F4 = 240*eisenstein_series(4, prec)
    F6 = -504*eisenstein_series(6, prec)
    R = QQ[[var]]
    if var != 'q':
        F4 = R(F4)
        F6 = R(F6)
    return (F4**3 - F6**2)/R(1728,prec)
