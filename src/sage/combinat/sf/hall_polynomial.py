
from skew_tableau import SkewTableau

def q_power(q, n):
    return 1-q**n


def q_product(q, n):
    if n > 1:
        return prod([1-q^i for i in range(1,n+1)])
    else:
        return 1

def q_binomial(q, n, r):
    if (r >= 0 and r <= n):
        return q_product(q, n+r)/( q_product(q,r) * q_product(q, n) )
    else:
        return 0

def phi(a, b, c, N, q):
    res = 0
    for r in range(min(a,b)):
        res += (-1)**r * q**(N*r+r*(r+1)/2) * q_binomial(q,a,r)*q_binomial(q,c-r,b-r)
    return res

def n_function(part):
    res = 0
    for i in range(len(part)):
        res += i*part[i]
    return res

def hall_polynomial_s(skt, q):
    skt = SkewTableau(skt)
    chain = skt.to_chain()

    r = len(chain) - 1
    l = len(chain[-1])

    tmp = []
    for i in range(r+1):
        tmp.append( )


    tmp[r+1] = tmp[r]



def weight_yam(w, l):
    return len( filter( lambda x: l in x, w ) )

def is_new_yamanouchi_word(w):
    w = [i for i in reversed(w)]
    for i in range(1, len(w)+1):
        wtemp = w[:i]
        for l in range(1, len(w)+1):
            if weight_yam(temp, l) < weight_yam(wtemp, l+1):
                return False

    return True



def hall_polynomial(mu, nu, lam, q):
    res = 0
    mu = Partition(mu).conjugate()
    nu = Partition(nu).conjugate()
    lam = Partition(lam).conjugate()

    for skt in SemistandardSkewTableau([lam, mu], nu):
        if is_new_yamanouchi_word(skt.to_word()):
            res += hall_polynomial_s(skt, q)

    return res
