
########################################################################
## Routines from Gonzalo Tornaria (7/9/07)
## for computing with ternary quadratic forms.
#######################################################################


#from sage.rings.rational_field import QQ

from sage.rings.integer_ring import ZZ
from sage.misc.functional import is_odd
from sage.rings.power_series_ring import PowerSeriesRing

from sage.libs.pari.all import pari
from sage.misc.misc import prod
from sage.rings.arith import factor


## TO DO -- Add second argument
#  def __call__(self,v,w=None):
#    if w==None:
#        return half(v * self._matrix_() * v)
#    else:
#      return v * self._matrix_() * w



def disc_Tornaria(self):
    """
    Returns the discriminant of the quadratc form,
    defined as

        det(B)      for even dimension
        det(B)/2    for odd dimension

    where 2Q(x) = x^t * B * x.
    """
    if is_odd(self.dim()):
      return  self.base_ring()(self.det() / 2)      ## This is not so good for characteristic 2.
    else:
      return self.det()


#def content(self):
#    """
#    Returns the GCD of the coefficients fo the quadratic form.
#
#    Warning: Only works over Euclidean domains... probably just ZZ. =|
#    """
#    return gcd(self.__coeffs)


#def is_primitive(self):
#    """
#    Checks if the form is a multiple of another form... only over ZZ for now.
#    """
#    return self.content() == 1



#def primitive(self):
#    """
#    Returns a primitive quadratic forms in the similarity class of the given form.
#
#    This only works when we have GCDs... so over ZZ.
#    """
#    c=self.content()
#    new_coeffs = [self.base_ring()(a/c)  for a in self.__coeffs]
#    return QuadraticForm(self.base_ring(), self.dim(), new_coeffs)



#def adjoint(self):
#    """
#    This gives the adjoint (integral) quadratic form associated to the
#    given form, essentially defined by taking the adjoint of the matrix.
#    """
#    if is_odd(self.dim()):
#        return QuadraticForm(self.matrix().adjoint()*2)
#    else:
#        return QuadraticForm(self.matrix().adjoint())


## Perhaps this is not needed...
#
#  def antiadjoint(self):
#    """
#
#
#    """
#    n=self.dim
#    try:
#      d=(polygen(self.R)**(n-1)-self.disc()).roots()[0][0]
#      if is_odd(n):
#        return self.adjoint() / d**(n-2) / 4
#      else:
#        return self.adjoint() / d**(n-2)
#    except IndexError:
#      raise ValueError, "not an adjoint"


## See above...
#
#  def is_adjoint(self):
#    try:
#      self.antiadjoint()
#    except ValueError:
#      return False
#    return True


def reciprocal(self):
    return self.adjoint().primitive() * self.content()


def omega(self):
    """
    This is the content of the adjont of the primitive associated quadratic form.

    Ref: See Dickson's "Studies in Number Theory".
    """
    return self.primitive().adjoint().content()

def delta(self):
    """
    This is the omega of the adjoint form,
    which is the same as the omega of the reciprocal form.
    """
    return self.adjoint().omega()


def level__Tornaria(self):
    """
    Hopefully this agrees with the usual level...
    """
    return self.base_ring()(self.disc()/self.omega()/self.content()**self.dim())


def discrec(self):
    return self.reciprocal().disc()




### Rational equivalence

def hasse_conductor(self):
    """
    This is the product of all primes where the Hasse invariant is -1.

    Note: For ternary forms, this is the discriminant of the associated
    quaternion algebra associated to the quadratic space
    """
    # * self.hasse(-1) ???
    return prod(filter(lambda(p):self.hasse_invariant(p)==-1, \
             map(lambda(x):x[0],factor(2*self.level()))))


### Genus theory

def basiclemma(self,M):
    a=self(self.basiclemmavec(M))
    assert gcd(a,M) == 1
    return a

def basiclemmavec(self,M):
    """
    Finds a vector where the value of the quadratic form is coprime to M.
    """
    V=FreeModule(self.base_ring(),self.dim())
    mat = self.matrix()
    vec = []
    mod = []
    M0 = abs(M)
    if M0 == 1:
        return V(0)

    for i in range(self.dim()):
        M1 = prime_to_m_part(M0, self[i,i])
        if M1 != 1:
            vec.append(V.i)
            mod.append(M1)
        M0 = M0/M1
        if M0 == 1:
            return __crt_list(vec,mod)

    for i in range(self.dim()):
        for j in range(i):
            M1 = prime_to_m_part(M0, self[i,j])
            if M1 != 1:
                vec.append(V.i + V.j)
                mod.append(M1)
            M0 = M0/M1
            if M0 == 1:
                return __crt_list(vec,mod)

    raise ValueError, "not primitive form"

def __crt_list(ls, ms):
    """
    Chinese remainder theorem for lists; find an element l
    such that l = ls[i] mod ms[i] for all i.
    """
    return sum(map(prod,zip,(ls,crt_basis(ms))))


### FIXME: get the rules for validity of characters straight...
### p=2 might be bad!!!
def xi(self,p):
    """
    Return the value of the genus characters Xi_p... which may be missing one character.
    We allow -1 as a prime.

    Reference: Dickson's "Studies in the Theory of Numbers"
    """
    if self.dim() == 2 and self.disc() % p:
        raise ValueError, "not a valid character"
    if self.dim() >= 3 and self.omega() % p:
        raise ValueError, "not a valid character"
    if (p == -1) or (p == 2):
        return kronecker(p, self.basiclemma(2))
    return kronecker(self.basiclemma(p), p)


def xi_rec(self,p):
    """
    Returns Xi(p) for the reciprocal form.
    """
    return self.reciprocal().xi(p)


def lll(self):
    """
    Returns an LLL-reduced form of Q (using Pari).
    """
    return self(self.matrix().lllgram())


def representation_number_list(self, B):
    """
    Returns the vector of representation numbers < B.
    """
    return pari(1).concat(self._pari_().qfrep(B-1, 1) * 2)


def ThetaByPari(self, B):
    """
    Returns the theta function up to O(q^B).
    """
    PSR = PowerSeriesRing(ZZ,'q')
    return PSR(self.representation_number_list(B), B)



def representation_vector_list(self, B):
    """
    Find all vectors v where Q(v) <= B.
    """
    n,m,vs = self._pari_().qfminim(2*(B-1), 10**8)
    if n != 2 * len(vs):
        raise RuntimeError("insufficient number of vectors")
    ms = [[] for _ in xrange(m/2+1)]
    for v in vs:
        ms[int(self(v))].append(v)
    return ms



### zeros

def is_singular_vector(self, v, p=0):
    """
    Determines if the vector v is on the conic Q(x) = 0 (mod p),
    and that this point is non-singular point of the conic.
    """
    if not self.is_zero(v, p):
        return False
    vm = vector(self.base_ring(), v) * self.matrix()
    if p != 0:
        vm = vm % p

    return (vm == 0)



