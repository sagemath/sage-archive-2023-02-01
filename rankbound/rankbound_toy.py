import numpy as np
from scipy.special import erfcx,spence
from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF
from sage.functions.log import log, exp
from sage.symbolic.constants import euler_gamma

npi = np.float(np.pi)
eg = np.float(euler_gamma)

def c_n(E,n,N=None):
    if N is None:
        N = E.conductor()

    n = ZZ(n)
    if n==1:
        return 0
    if not n.is_prime_power():
        return 0

    n_float = RDF(n)
    logn = log(n_float)

    if n.is_prime():
        return -E.ap(n)*logn/n_float
    else:
        p,e = n.factor()[0]
        ap = E.ap(p)
        logp = log(RDF(p))
        if p.divides(N):
            return - ap**e*logp/n_float
        c = np.complex(ap,np.sqrt(4*p-ap**2))/2
        aq = np.round(2*np.real(c**e))
        return RDF(-aq*logp/n_float)

def rankbound(E,Delta=1,N=None):
    if N is None:
        N = E.conductor()

    t = RDF(Delta*2*npi)
    expt = RDF(exp(t))

    C0 = -eg+np.log(np.sqrt(N)/(2*npi))
    u = RDF(t*C0)

    w = RDF(npi**2/6-spence(1-RDF(1)/expt))

    y = RDF(0)
    n = int(0)
    #for n in xrange(np.floor(np.exp(t))):
    while n < expt:
        n += 1
        cn  = c_n(E,n,N)
        if cn!=0:
            logn = log(RDF(n))
            y += cn*(t-logn)

    return 2*(u+w+y)/(t**2)


def rankbound_gaussian(E,t=4,N=None):
    if N is None:
        N = E.conductor()

    q = np.sqrt(2/npi)
    u = q/t*(-eg+np.log(np.sqrt(N)/(2*npi)))
    w = np.float(0)
    for k in range(1,1001):
        w += q/k-t*erfcx(t*k/np.sqrt(2))
    w = w/t

    y = 0.
    for n in range(RDF(20)**t):
        #cn = c_n(E,n)
        cn = c_n(E,n,N)
        if cn != 0:
            logn = np.log(n)
            y += cn*np.exp(-(logn/t)**2/2)
    y = y*q/t

    #print(u,w,y)

    return RDF(u+w+y+0.1)
