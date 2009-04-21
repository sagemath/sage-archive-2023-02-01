from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense, Polynomial
from sage.rings.infinity import infinity
from sage.libs.all import pari_gen
from sage.structure.factorization import Factorization
import sage.rings.padics.misc

class Polynomial_padic_flat(Polynomial_generic_dense):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False, absprec=None):
        """
        Initialization function for the class  Polynomial_padic_flat.
        """

        if x is None:
            Polynomial_generic_dense.__init__(self, parent, x = None, is_gen = is_gen)
            return
        R = parent.base_ring()
        if sage.rings.fraction_field_element.is_FractionFieldElement(x):
            if x.denominator() != 1:
                raise TypeError, "denominator must be 1"
            else:
                x = x.numerator()
        if isinstance(x, Polynomial):
            if x.base_ring() is R:
                x = list(x.list())
            else:
                x = [R(a) for a in x.list()]
        elif isinstance(x, list):
            if check:
                x = [R(a) for a in x]
        elif isinstance(x, dict):
            if check:
                m = infinity
                zero = R(0)
                n = max(x.keys())
                v = [zero for _ in xrange(n+1)]
                for i, z in x.iteritems():
                    v[i] = R(z)
                    m = min(m, v[i].precision_absolute())
                x = v
            else:
                m = sage.rings.padics.misc.min(a.precision_absolute() for a in x.values())
            if not absprec is None:
                m = min(m, absprec)
            Polynomial_generic_dense.__init__(self, parent, x, absprec = m)
            return
        elif isinstance(x, pari_gen):
            x = [R(w) for w in x.Vecrev()]
        else:
            x = [R(x)]
        if absprec is None:
            absprec = infinity
        m = min([a.precision_absolute() for a in x] + [absprec])
        Polynomial_generic_dense.__init__(self, parent, x, absprec = m)

    def _mul_(self, right):
        """
        Returns the product of this Polynomial_padic_flat by right.
        """
        return self._mul_generic(right)

    def _repr(self, name=None):
        r"""
        EXAMPLES:
            sage: R.<w> = PolynomialRing(Zp(5, prec=5, type = 'capped-abs', print_mode = 'val-unit'))
            sage: f = 24 + R(4/3)*w + w^4
            sage: f._repr()
            '(1 + O(5^5))*w^4 + (1043 + O(5^5))*w + (24 + O(5^5))'
            sage: f._repr(name='z')
            '(1 + O(5^5))*z^4 + (1043 + O(5^5))*z + (24 + O(5^5))'

        AUTHOR:
            -- David Roe (2007-03-03), based on Polynomial_generic_dense._repr()
        """
        s = " "
        n = m = self.degree()
        if name is None:
            name = self.parent().variable_name()
        #atomic_repr = self.parent().base_ring().is_atomic_repr()
        coeffs = self.list()
        for x in reversed(coeffs):
            if x != 0:
                if n != m:
                    s += " + "
                #x = repr(x)
                x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
            n -= 1
        if s==" ":
            return "0"
        return s[1:]

    def factor(self, absprec = None):
        r"""
        Returns the factorization of this Polynomial_padic_flat.

        EXAMPLES:
            sage: R.<w> = PolynomialRing(Zp(5, prec=5, type = 'capped-abs', print_mode = 'val-unit'))
            sage: f = w^5-1
            sage: f.factor()
            ((1 + O(5^5))*w + (624 + O(5^5))) * ((1 + O(5^5))*w^4 + (2501 + O(5^5))*w^3 + (1876 + O(5^5))*w^2 + (1251 + O(5^5))*w + (626 + O(5^5)))

        See \#4038:
            sage: E = EllipticCurve('37a1')
            sage: K =Qp(7,10)
            sage: EK = E.base_extend(K)
            sage: E = EllipticCurve('37a1')
            sage: K = Qp(7,10)
            sage: EK = E.base_extend(K)
            sage: g = EK.division_polynomial_0(3)
            sage: g.factor()
            (3 + O(7^10)) * ((1 + O(7^10))*x + (1 + 2*7 + 4*7^2 + 2*7^3 + 5*7^4 + 7^5 + 5*7^6 + 3*7^7 + 5*7^8 + 3*7^9 + O(7^10))) * ((1 + O(7^10))*x^3 + (6 + 4*7 + 2*7^2 + 4*7^3 + 7^4 + 5*7^5 + 7^6 + 3*7^7 + 7^8 + 3*7^9 + O(7^10))*x^2 + (6 + 3*7 + 5*7^2 + 2*7^4 + 7^5 + 7^6 + 2*7^8 + 3*7^9 + O(7^10))*x + (2 + 5*7 + 4*7^2 + 2*7^3 + 6*7^4 + 3*7^5 + 7^6 + 4*7^7 + O(7^10)))
        """
        if self == 0:
            raise ValueError, "Factorization of 0 not defined"
        if absprec is None:
            absprec = min([x.precision_absolute() for x in self.list()])
        else:
            absprec = Integer(absprec)
        if absprec <= 0:
            raise ValueError, "absprec must be positive"
        G = self._pari_().factorpadic(self.base_ring().prime(), absprec)
        pols = G[0]
        exps = G[1]
        F = []
        R = self.parent()
        for i in xrange(len(pols)):
            f = R(pols[i], absprec = absprec)
            e = int(exps[i])
            F.append((f,e))

        #if R.base_ring().is_field():
        #    # When the base ring is a field we normalize
        #    # the irreducible factors so they have leading
        #    # coefficient 1.
        #    for i in range(len(F)):
        #        cur = F[i][0].leading_coefficient()
        #        if cur != 1:
        #            F[i] = (F[i][0].monic(), F[i][1])
        #    return Factorization(F, self.leading_coefficient())
        #else:
        #    # When the base ring is not a field, we normalize
        #    # the irreducible factors so that the leading term
        #    # is a power of p.  We also ensure that the gcd of
        #    # the coefficients of each term is 1.
        c = self.leading_coefficient().valuation()
        u = self.base_ring()(1)
        for i in range(len(F)):
            upart = F[i][0].leading_coefficient().unit_part().lift_to_precision(absprec)
            lval = F[i][0].leading_coefficient().valuation()
            if upart != 1:
                F[i] = (F[i][0] / upart, F[i][1])
                u *= upart ** F[i][1]
            c -= lval * F[i][1]
        if c != 0:
            F.append((self.parent()(self.base_ring().prime_pow(c)), 1))
            u = u.add_big_oh(absprec - c)
        return Factorization(F, u)
