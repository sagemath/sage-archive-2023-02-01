import sage.rings.polynomial_element

#import sage.rings.padics.padic_ring
#import sage.rings.padics.padic_ring_element
#import sage.rings.padics.padic_field
#import sage.rings.padics.padic_field_element

#class pAdicPolynomial(sage.rings.polynomial_element.Polynomial):
#    def _mul_(self):
#        """
#        Override the default Karatsuba because it leaks precision
#        """
#        if right == 0 or self == 0:
#           return self.polynomial(0)
#        return self._mul_generic(right)
#
#    def factor(self):
#        """
#        Factors self.
#
#        Algorithm -- Round4
#        """

#def _round4_(f):
#Can assume f is square-free, monic


#def _hensel_lift(f, g, h):


#def _hensel_internal(f, g, h, s, t, p, n): #f = gh (mod p^n), sg + ht = 1 (


# root finding over p-adic fields
def newton_lift(f,a):
        """
        If v(f(a))>v(f'(a)) then lift a to a root of f
        """
        f_   = f.derivative()
        vfa  = f(a).valuation()
        vf_a = f_(a).valuation()
        if vfa <=  vf_a:
            raise TypeError, "the approximation of the root is not sufficient for newton lifting"

        for i in range(0, ceil(log(a.parent().precision_cap(),2))):
          a = a - f(a)//f_(a)
        return a

def content(f):
        """
        The content of the polynomial f.
        The content is the greatest common divisor of the coefficients of f.
        """
        v = Min([a.valuation() for a in f.coeffs()])
        return f.base_ring().uniformizer()**v


def roots_squarefree(f):
        """
        find the roots of a squarefree polynomial
        """
        Lx   = f.parent()
        x    = Lx.gen()
        L    = Lx.base_ring()
        pi   = L.uniformizer()
        rcfx = Lx.change_ring(L.residue_class_field())

        #C   = [[f//f.content(),0,0]]
        C   = [[f//content(f),0,0]]
        G   = []

        while not C == []:
            for c in C:
                C.remove(c)
                psi = c[0]
                delta = c[1]
                s = c[2]
                R = rcfx(psi).roots()
                for b in R:
                    beta = L(b[0])
                    npsi = psi(pi*x+beta)
                    npsi = npsi//content(npsi)
                    if Min([a.precision_absolute() for a in psi.coeffs()])<1:
                        raise TypeError, "insufficient precision in padic root finding"
                    #print "v(f(a))",f(delta+pi**s*beta).valuation()
                    if rcfx(npsi).degree() == 1:
                        G.append(newton_lift(f,delta+pi**s*beta))
                    else:
                        C.append([npsi,delta+pi**s*beta,s+1])
        return G

def roots_padic(f):
    if f == 0:
      raise TypeError, "hmm, how many roots does the zero polynomial have ?"
    if f.degree == 0:
      return []
    n = f.degree()
    f_ = f.derivative()
    g = gcd(f,f_)
    r = roots_squarefree(f//g)
    rm = [[a,1] for a in r]
    nr = len(r)
    #print "////////////////////////////////"
    if g.degree() != 0:
      while nr<n and f_.degree()<>0:
        for i in range(0,len(r)-1):
           if f_(r[i]) == 0:
              rm[i][1]+=1
              nr+=1
           f_=f_.derivative()
    return [(a[0],a[1]) for a in rm]



#    gr = roots_padic(g)
#      for a in gr:
#        if a[0] in r:
#          rm[r.index(a[0])][1]+=a[1]
#    return [(a[0],a[1]) for a in rm]
