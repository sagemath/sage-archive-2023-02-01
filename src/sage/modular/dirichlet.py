r"""
Dirichlet characters

A \class{DirichletCharacter} is the extension of a homomorphism
$$
  (\Z/N\Z)^* \to R^*,
$$
for some ring $R$, to the map $\Z/N\Z \to R$ obtained by sending
those $x\in\Z/N\Z$ with $\gcd(N,x)>1$ to $0$.

EXAMPLES:
    sage: G, x = DirichletGroup(35).objgens()
    sage: e = x[0]*x[1]^2; e
    [zeta_12^3, zeta_12^2 - 1]
    sage: e.order()
    12

AUTHORS:
    -- William Stein (2006-01-07): added more examples
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
#
# CHANGE LOG:
#
# 2005-09-02 (was): Fixed bug in comparison of Dirichlet characters.
# It was checking that their values were the same, but not checking
# that they had the same level!
#
#*****************************************************************************

import random
import weakref

import sage.rings.arith as arith
import sage.misc.misc as misc
import sage.rings.all as rings
import sage.structure.gens as gens
import sage.rings.number_field.number_field as number_field
import sage.rings.coerce as coerce
from sage.structure.element import MultiplicativeGroupElement

def TrivialCharacter(N, base_ring=rings.RationalField()):
    return DirichletGroup(N, base_ring)(1)

def is_DirichletCharacter(x):
    return isinstance(x, DirichletCharacter)

class DirichletCharacter(MultiplicativeGroupElement):
    """
    A Dirichlet character
    """
    def __init__(self, parent, values_on_gens):
        r"""
        Create with \code{DirichletCharacter(parent, values_on_gens)}

        INPUT:
            parent -- DirichletGroup, a group of Dirichlet characters
            values_on_gens -- list of ring elements, the values of the
                              Dirichlet character on the chosen generators
                              of $(\Z/N\Z)^*$.
        OUTPUT:
            DirichletCharacter -- a Dirichlet character

        EXAMPLES:
            sage: G, e = DirichletGroup(13).objgen()
            sage: G
            Group of Dirichlet characters of modulus 13 over Cyclotomic Field of order 12 and degree 4
            sage: e
            [zeta_12]
            sage: loads(e.dumps()) == e
            True

            sage: G, x = DirichletGroup(35).objgens()
            sage: e = x[0]*x[1]; e
            [zeta_12^3, zeta_12^2]
            sage: e.order()
            12
            sage: loads(e.dumps()) == e
            True
        """
        MultiplicativeGroupElement.__init__(self, parent)
        if len(values_on_gens) != len(parent.unit_gens()):
            raise ValueError, \
                  "wrong number of values(=%s) on unit gens (want %s)"%( \
                   values_on_gens,len(parent.unit_gens()))
        R = parent.base_ring()
        self.__values_on_gens = [R(x) for x in values_on_gens]
        self.__modulus = parent.modulus()

    def __eval_at_minus_one(self):
        """
        Efficiently evalute the character at -1 using knowledge of its
        order.   This is potentially much more efficient than computing
        the value of -1 directly using dlog and a large power of the
        image root of unity.

        We use the following.
        Proposition: Suppose eps is a character mod p^n, where p is a prime.
        Then eps(-1) = -1 if and only if
               p = 2 and the factor of eps at 4 is nontrivial
           or
               p > 2 and 2 does not divide euler_phi(p^n)/order(eps).
        """
        try:
            return self.__value_at_minus_one
        except AttributeError:
            D = self.decomposition()
            val = self.base_ring()(1)
            for e in D:
                if e.modulus() % 2 == 0:
                    val *= e.__values_on_gens[0]
                elif (arith.euler_phi(e.parent().modulus()) / e.order()) % 2 != 0:
                    val *= -1
            self.__value_at_minus_one = val
        return self.__value_at_minus_one

    def __call__(self, m):
        m = int(m%self.__modulus)
        try:
            return self.__values[m]
        except AttributeError:
            pass
        if m == self.__modulus - 1:
            return self.__eval_at_minus_one()
        self.values()  # compute all values
        return self.__values[m]

    def change_base(self, R):
        """
        Returns the base extension of self to the ring R.
        """
        G = self.parent().change_base(R)
        return G(self)

    def __cmp__(self, other):
        if not isinstance(other, DirichletCharacter) or other.parent() != self.parent():
            return coerce.cmp(self, other)
        return misc.generic_cmp(self.__values_on_gens, other.__values_on_gens)

    def __hash__(self):
        return hash(str(self.__values_on_gens))

    def __invert__(self):
        values_on_gens = [1/x for x in self.__values_on_gens]
        return DirichletCharacter(self.parent(), values_on_gens)

    def _mul_(self,  other):
        values_on_gens = [self.__values_on_gens[i]*other.__values_on_gens[i]
                          for i in range(len(self.__values_on_gens))]
        return DirichletCharacter(self.parent(), values_on_gens)

    def __pow__(self, exp):
        vals = [z**exp for z in self.__values_on_gens]
        return DirichletCharacter(self.parent(), vals)

    def __repr__(self):
        return str(self.__values_on_gens)

    def base_ring(self):
        """
        Returns the base ring of this Dirichlet character.

        EXAMPLES:
            sage: G = DirichletGroup(11)
            sage: G.gen(0).base_ring()
            Cyclotomic Field of order 10 and degree 4
            sage: G = DirichletGroup(11, RationalField())
            sage: G.gen(0).base_ring()
            Rational Field
        """
        return self.parent().base_ring()

    def bernoulli(self, k):
        r"""
        Returns the generalized Bernoulli number $B_{k,eps}$.

        Let eps be this character (not necessarily primitive), and
        let $k \geq 0$ be an integer weight.  This function computes
        the (generalized) Bernoulli number $B_{k,eps}$, e.g., as defined
        on page 44 of Diamond-Im:
        $$
          \sum_{a=1}^{N} \eps(a) t*e^{at} / (e^{Nt}-1)
                 = sum_{k=0}^{\infty} B_{k,eps}/{k!} t^k.
        $$
        where $N$ is the modulus of $\eps$.

        EXAMPLES:
            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.bernoulli(5)
            7430/13*zeta_12^3 - 34750/13*zeta_12^2 - 11380/13*zeta_12 + 9110/13
        """
        try:
            self.__bernoulli
        except AttributeError:
            self.__bernoulli = {}
        if self.__bernoulli.has_key(k):
            return self.__bernoulli[k]
        N = self.modulus()
        K = self.base_ring()

        if N != 1 and self(-1) != K((-1)**k):
            return K(0)

        # TODO: The following is very slow, since poly rings over a
        # very slow field are very slow... (this could change as SAGE
        # evolves).
        if False:
            R = rings.PowerSeriesRing(K, variable="t")
            t = R.gen()
            prec = k+2   # todo: fix this
            F = sum([(self(a) * t * (a*t).exp(prec)) / ((N*t).exp(prec) - 1) \
                     for a in range(1,N)])
            self.__bernoulli[k] = F[k]*arith.factorial(k)
            return self.__bernoulli[k]


        # This is better since it computes the same thing, but requires
        # no arith in a poly ring over a number field.
        prec = k+2
        R = rings.PowerSeriesRing(rings.RationalField())
        t = R.gen()
        # g(t) = t/(e^{Nt}-1)
        g = t/((N*t).exp(prec) - 1)

        # h(n) = g(t)*e^{nt}
        h = [0] + [g * ((n*t).exp(prec)) for n in range(1,N+1)]

        ber = sum([self(a)*h[a][k] for a in range(1,N+1)]) * arith.factorial(k)

        self.__bernoulli[k] = ber
        return self.__bernoulli[k]

    def change_base_ring(self, R):
        """
        Tries to compute the Dirichlet character over R that has the
        same values as this character, and if successful returns that
        character.
        """
        if self.base_ring() == R:
            return self
        return DirichletGroup(self.modulus(),R)(self)

    def conductor(self):
        """
        Computes and returns the conductor of this character.
        """
        try:
            return self.__conductor
        except AttributeError:
            pass
        if self.modulus() == 1 or self.is_trivial():
            self.__conductor = 1
            return self.__conductor
        F = arith.factor(self.modulus())
        if len(F) > 1:
            self.__conductor = misc.mul([d.conductor() for d in self.decomposition()])
            return self.__conductor
        p = F[0][0]
        # When p is odd, and x =/= 1, the conductor is the smallest p**r such that
        #   Order(x) divides EulerPhi(p**r) = p**(r-1)*(p-1).
        # For a given r, whether or not the above divisiblity holds
        # depends only on the factor of p**(r-1) on the right hand side.
        # Since p-1 is coprime to p, this smallest r such that the
        # divisibility holds equals Valuation(Order(x),p)+1.
        self.__conductor = p**(arith.valuation(self.order(),p) + 1)
        if p == 2 and F[0][1] > 2 and self.__values_on_gens[1].multiplicative_order() != 1:
            self.__conductor *= 2;
        return self.__conductor

    def decomposition(self):
        """
        Return the decomposition of self as a product of Dirichlet characters
        of prime power modulus, where the prime powers exactly divide the
        modulus of this character.
        """
        try:
            return self.__decomp
        except AttributeError:
            pass
        D = self.parent().decomposition()
        vals = [[z] for z in self.__values_on_gens]
        R = self.base_ring()
        if self.modulus()%8 == 0:   # 2 factors at 2.
            vals[0].append(vals[1][0])
            del vals[1]
        self.__decomp = [D[i](vals[i]) for i in range(len(D))]
        return self.__decomp

    def extend(self, M):
        """
        Returns the extension of this character to a Dirichlet character
        modulo the multiple M of the modulus.
        """
        if M % self.modulus() != 0:
            raise ArithmeticError, "M(=%s) must be a multiple of the modulus(=%s)"%(M,self.modulus())
        H = DirichletGroup(M, self.base_ring())
        return H(self)

    def gauss_sum(self, a=1):
        r"""
        Return a Gauss sum associated to this Dirichlet character.

        The Gauss sum associated to $\chi$ is
         $$
              g_a(\chi) = \sum_{r \in \Z/m\Z} \chi(r)\,\zeta^{ar},
         $$
        where $m$ is the modulus of $\chi$ and $\zeta$ is a primitive
        $m$th root of unity, i.e., $\zeta$ is
        \code{self.parent().zeta()}.

        FACTS: If the modulus is a prime $p$ and the character is
        nontrivial, then the Gauss sum has absolute value $\sqrt{p}$.

        CACHING: Computed Gauss sums are \emph{not} cached with this
        character.

        EXAMPLES:
            sage: G = DirichletGroup(3)
            sage: e = G([-1])
            sage: e.gauss_sum(1)
            2*zeta_6 - 1
            sage: e.gauss_sum(2)
            -2*zeta_6 + 1
            sage: norm(e.gauss_sum())
            3

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.gauss_sum()
            -zeta_156^46 + zeta_156^45 + zeta_156^42 + zeta_156^41 + 2*zeta_156^40 + zeta_156^37 - zeta_156^36 - zeta_156^34 - zeta_156^33 - zeta_156^31 + 2*zeta_156^30 + zeta_156^28 - zeta_156^24 - zeta_156^22 + zeta_156^21 + zeta_156^20 - zeta_156^19 + zeta_156^18 - zeta_156^16 - zeta_156^15 - 2*zeta_156^14 - zeta_156^10 + zeta_156^8 + zeta_156^7 + zeta_156^6 + zeta_156^5 - zeta_156^4 - zeta_156^2 - 1
            sage: factor(norm(e.gauss_sum()))
            13^24
        """
        G = self.parent()
        K = G.base_ring()
        if not rings.is_CyclotomicField(K) or rings.is_RationalField(K):
            raise NotImplementedError, "Gauss sums only currently implemented when the base ring is a cyclotomic field or QQ."
        g = 0
        m = G.modulus()
        L = rings.CyclotomicField(arith.lcm(m,G.zeta_order()))
        zeta = L.gen(0)
        n = zeta.multiplicative_order()
        zeta = zeta ** (n // m)
        if a != 1:
            zeta = zeta**a
        z = 1
        for c in self.values()[1:]:
            z *= zeta
            g += L(c)*z
        return g

    def gauss_sum_numerical(self, prec=53, a=1):
        r"""
        Return a Gauss sum associated to this Dirichlet character as
        an approximate complex number with prec bits of precision.

        INPUT:
            prec -- integer (deafault: 53), *bits* of precision
            a -- integer, as for gauss_sum.

        The Gauss sum associated to $\chi$ is
         $$
              g_a(\chi) = \sum_{r \in \Z/m\Z} \chi(r)\,\zeta^{ar},
         $$
        where $m$ is the modulus of $\chi$ and $\zeta$ is a primitive
        $m$th root of unity, i.e., $\zeta$ is
        \code{self.parent().zeta()}.

        EXAMPLES:
            sage: G = DirichletGroup(3)
            sage: e = G.0
            sage: e.gauss_sum_numerical()
            0.00000000000000055511151231257827 + 1.7320508075688772*I
            sage: abs(e.gauss_sum_numerical())
            1.7320508075688772
            sage: sqrt(3)
            1.7320508075688772
            sage: e.gauss_sum_numerical(a=2)
            -0.0000000000000011102230246251565 - 1.7320508075688772*I
            sage: e.gauss_sum_numerical(a=2, prec=100)
            0.0000000000000000000000000000047331654313260708324703713916967 - 1.7320508075688772935274463415062*I
            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.gauss_sum_numerical()
            -3.0749720589952387 + 1.8826966926190174*I
            sage: abs(e.gauss_sum_numerical())
            3.6055512754639896
            sage: sqrt(13)
            3.6055512754639891
        """
        G = self.parent()
        K = G.base_ring()
        if not rings.is_CyclotomicField(K) or rings.is_RationalField(K):
            raise NotImplementedError, "Gauss sums only currently implemented when the base ring is a cyclotomic field or QQ."
        phi = K.complex_embedding(prec)
        CC = phi.codomain()

        g = 0
        m = G.modulus()
        zeta = CC.zeta(m)
        if a != 1:
            zeta = zeta**a
        z = 1
        for c in self.values()[1:]:
            z *= zeta
            g += phi(c)*z
        return g




    def is_even(self):
        r"""
        Return \code{True} if and only if $\eps(-1) = 1$.
        """
        return not self.is_odd()

    def is_odd(self):
        r"""
        Return \code{True} if and only if $\eps(-1) \neq 1$.
        """
        try:
            return self.__is_odd
        except AttributeError:
            pass
        self.__is_odd = (self(-1) != self.base_ring()(1))
        return self.__is_odd

    def is_primitive(self):
        try:
            return self.__is_primitive
        except AttributeError:
            pass
        self.__is_primitive = (self.conductor() == self.modulus())
        return self.__is_primitive

    def is_trivial(self):
        r"""
        Returns \code{True} if this is the trivial character, i.e., has
        order 1.
        """
        try:
            self.__is_trivial
        except AttributeError:
            pass
        self.__is_trivial = True
        R = self.base_ring()
        for x in self.__values_on_gens:
            if x != R(1):
                self.__is_trivial = False
                return False
        return self.__is_trivial

    def maximize_base_ring(self):
        r"""
        Let
        $$
               \eps : (\Z/N\Z)^* \to \Q(\zeta_n)
        $$
        be a Dirichlet character.  This function returns an equal
        Dirichlet character
        $$
               \chi : (\Z/N\Z)^* \to \Q(\zeta_m)
        $$
        where $m$ is the least common multiple of $n$ and the
        exponent of $(\Z/N\Z)^*$.
        """
        g = rings.IntegerModRing(self.modulus()).unit_group_exponent()
        if g == 1:
            g = 2
        z = self.base_ring().zeta()
        n = z.multiplicative_order()
        m = arith.LCM(g,n)
        if n == m:
            return self
        K = rings.CyclotomicField(m)
        return self.change_base_ring(K)

    def minimize_base_ring(self):
        r"""
        Return a Dirichlet character that equals this one, but over as
        small a subfield (or subring) of the base ring as possible.

        \note{This function is currently only implemented when the
        base ring is a number field.}
        """

        if self.is_trivial(): return self
        if isinstance(self.base_ring(),rings.RationalField): return self

        if isinstance(self.base_ring(),number_field.NumberField_generic):
            if self.order() <= 2:
                return self.change_base(rings.RationalField())
            if arith.euler_phi(self.order()) == self.base_ring().degree():
                return self
            K = rings.CyclotomicField(self.order())
            return self.change_base(K)

        raise NotImplementedError, "minimize_base_ring is currently " + \
              "only implemented when the base ring is a number field."

    def modulus(self):
        """
        The modulus of this character.
        """
        return self.__modulus

    def multiplicative_order(self):
        """
        The order of this character.
        """
        try:
            return self.__order
        except AttributeError:
            pass
        self.__order = int(arith.LCM([x.multiplicative_order() \
                                      for x in self.__values_on_gens]))
        return self.__order

    def restrict(self, M):
        """
        Returns the restriction of this character to a Dirichlet
        character modulo the divisor M of the modulus, which must also
        be a multiple of the conductor of this character.
        """
        M = int(M)
        if self.modulus()%M != 0:
            raise ValueError, "M(=%s) must divide the modulus(=%s)"%(M,self.modulus())
        if M%self.conductor() != 0:
            raise ValueError, "conductor(=%s) must divide M(=%s)"%(self.conductor(),M)
        H = DirichletGroup(M, self.base_ring())
        return H(self)

    def values(self):
        """
        Returns a list of the values of this character on each integer
        between 0 and the modulus.
        """
        try:
            return self.__values
        except AttributeError:
            pass
        # Build cache of all values of the Dirichlet character.
        # I'm going to do it this way, since in my app the modulus
        # is *always* small and we want to evaluate the character
        # *a* *lot*.
        R = self.parent().base_ring()
        zero = R(0)
        mod = self.__modulus
        x = [zero for _ in range(mod)]
        if self.is_trivial():  # easy special case
            for n in range(mod):
                # todo: optimization idea -- factor phi(n), then
                # "cross out" all multiples of each prime factor of phi(n).
                if arith.GCD(n,mod) == 1:
                    x[n] = 1
            self.__values = x
            return
        #end
        # we lift the gens to the int type here, since it is a lot faster
        gens = [z.lift() for z in self.parent().unit_gens()]
        exponents = [0 for _ in range(len(gens))]
        Z = self.parent().integers_mod()
        n = 1
        last = [g.multiplicative_order()-1 for g in self.parent().unit_gens()]
        stop = list(last)
        stop[len(stop)-1] += 1
        value = R(1)
        val_on_gen = self.__values_on_gens
        only_int_vals = True
        tmp = []
        for z in val_on_gen:
            if z == R(1):
                tmp.append(1)
            elif z == R(-1):
                tmp.append(-1)
            else:
                only_int_vals = False
                break
        #end
        if only_int_vals:
            val_on_gen = tmp
            value = 1
        while exponents != stop:
            # record character value on n
            x[int(n)] = value
            # iterate:
            #   increase the exponent vector by 1,
            #   increase n accordingly, and increase value
            exponents[0] += 1   # inc exponent
            value *= val_on_gen[0]  # inc value
            n *= gens[0]
            n %= mod
            i = 0
            while i < len(exponents)-1 and exponents[i] > last[i]:
                exponents[i] = 0
                # now increment position i+1:
                exponents[i+1] += 1
                value *= val_on_gen[i+1]
                n *= gens[i+1]
                n %= mod
                i += 1
            #end
        #end
        if only_int_vals:
            x = [R(z) for z in x]
        self.__values = x
        return self.__values

    def values_on_gens(self):
        """
        Returns a list of the values of this character on each of the
        minimal generators of $(\Z/N\Z)^*$, where $N$ is the modulus.
        """
        return self.__values_on_gens

# This is so there is only one DirichletGroup with given parameters
# at a time.
_objsDirichletGroup = {}
def DirichletGroup(modulus, base_ring=None, zeta=None, zeta_order=None):
        key = (base_ring, modulus, zeta, zeta_order)
        if _objsDirichletGroup.has_key(key):
            x = _objsDirichletGroup[key]()
            if x != None: return x
        R = DirichletGroup_class(modulus, base_ring, zeta, zeta_order)
        _objsDirichletGroup[key] = weakref.ref(R)
        return R

def is_DirichletGroup(x):
    return isinstance(x, DirichletGroup_class)

class DirichletGroup_class(gens.MultiplicativeAbelianGenerators):
    """
    Group of Dirichlet characters modulo $N$ over a given base ring $R$.
    """
    def __init__(self, modulus, base_ring=None, zeta=None, zeta_order=None):
        r"""
        One can create the group of Dirichlet characters modulo~$N$
        with values in the subgroup $\langle \zeta_n\rangle$ of the
        multiplicative group of the \code{base_ring}.  If the
        base_ring is omitted then we use $\Q(\zeta_n)$, where $n$ is
        the exponent of $(\Z/N\Z)^*$.  If $\zeta$ is omitted then we
        compute and use a maximal-order zeta in base_ring, if
        possible.

        INPUT:
            modulus -- int
            base_ring -- Ring (optional), where characters take their values
                         (should be an integral domain).
            zeta -- Element (optional), element of base_ring; zeta is a root of unity
            zeta_order -- int (optional), the order of zeta

        OUTPUT:
            DirichletGroup -- a group of Dirichlet characters.

        NOTES:
            Uniqueness -- If a group is created with the same parameters
            as another DirichletGroup still in memory, then the same group
            is returned instead of a new group defined by the same
            parameters.

        EXAMPLES:
            The default base ring is a cyclotomic field of order the exponent
            of $(\Z/N\Z)^*$.

                sage: DirichletGroup(20)
                Group of Dirichlet characters of modulus 20 over Cyclotomic Field of order 4 and degree 2

            We create the group of Dirichlet character mod 20 with values
            in the rational numbers:
                sage: G = DirichletGroup(20, Q); G
                Group of Dirichlet characters of modulus 20 over Rational Field
                sage: G.order()
                4
                sage: G.base_ring()
                Rational Field

            The elements of G print as lists giving the values of the
            character on the generators of $(Z/NZ)^*$:
                sage: list(G)
                [[1, 1], [-1, 1], [1, -1], [-1, -1]]

            Next we construct the group of Dirichlet character mod 20, but
            with values in Q(zeta_n):
                sage: G = DirichletGroup(20)
                sage: G.list()
                [[1, 1], [-1, 1], [1, zeta_4], [-1, zeta_4], [1, -1], [-1, -1], [1, -zeta_4], [-1, -zeta_4]]

            We next compute several invariants of G:
                sage: G.gens()
                ([-1, 1], [1, zeta_4])
                sage: G.unit_gens()
                [11, 17]
                sage: G.zeta()
                zeta_4
                sage: G.zeta_order()
                4

            In this example we create a Dirichlet character with values in a
            number field.  We have to give zeta, but not its order.
                sage: R = PolynomialRing(Q); x = R.gen()
                sage: K = NumberField(x^4 + 1); a = K.gen(0)
                sage: G = DirichletGroup(5, K, a); G
                Group of Dirichlet characters of modulus 5 over Number Field in a with defining polynomial x^4 + 1
                sage: G.list()
                [[1], [a^2], [-1], [-a^2]]


            sage: G, e = DirichletGroup(13).objgens()
            sage: loads(G.dumps()) == G
            True

            sage: G = DirichletGroup(19, GF(5))
            sage: loads(G.dumps()) == G
            True
        """
        if base_ring == None:
            if not (zeta==None and zeta_order==None):
                raise ValueError, "zeta and zeta_order must be None if base_ring not specified."
            e = rings.IntegerModRing(modulus).unit_group_exponent()
            base_ring = rings.CyclotomicField(e)

        if zeta == None:
            e = rings.IntegerModRing(modulus).unit_group_exponent()
            try:
                zeta = base_ring.zeta(e)
            except ValueError:
                zeta = base_ring.zeta()
            zeta_order = zeta.multiplicative_order()

        elif zeta_order == None:
            zeta_order = zeta.multiplicative_order()

        #self._base_ring = base_ring
        self._zeta = zeta
        self._zeta_order = zeta_order
        self._modulus = modulus
        self._integers = rings.IntegerModRing(modulus)

    def change_base(self, R):
        """
        Returns the Dirichlet group over R with the same modulus as self.
        """
        return DirichletGroup(self.modulus(), R)

    def __call__(self, x):
        if isinstance(x, (int,rings.Integer)) and x==1:
            R = self.base_ring()
            x = [R(1) for _ in self.unit_gens()]
        if isinstance(x, list):  # list of values on each unit generator
            return DirichletCharacter(self, x)
        elif isinstance(x, DirichletCharacter):  # coercion
            if x.parent() == self:
                return x
            if self.modulus() % x.conductor() != 0:
                raise TypeError, "conductor (=%s) must divide modulus (=%s)"%(\
                    x.conductor(), self.modulus())
            a = []
            R = self.base_ring()
            for u in self.unit_gens():
                v = u.lift()
                # have to do this, since e.g., unit gens mod 11 are not units mod 22.
                while arith.GCD(x.modulus(),int(v)) != 1:
                    v += self.modulus()
                a.append(R(x(v)))
            return self(a)
        raise TypeError, "No coercion of %s into %s defined."%(x, self)

    def __cmp__(self, other):
        if not is_DirichletGroup(other):
            return -1
        if other.modulus() < self.modulus():
            return -1
        if other.modulus() > self.modulus():
            return 1
        if other.base_ring() < self.base_ring():
            return -1
        if other.base_ring() > self.base_ring():
            return 1
        if other._zeta != self._zeta:
            return -1
        return 0

    def __len__(self):
        return self.order()

    def __repr__(self):
        return "Group of Dirichlet characters of modulus %s over %s"%\
               (self.modulus(),self.base_ring())

    def primitive_character(self):
        """
        Returns the primitive character associated to self.
        """
        return self.restrict(self.conductor())

    def base_ring(self):
        """
        Returns the base ring of self.
        """
        #return self._base_ring
        return self._zeta.parent()

    def decomposition(self):
        """
        Returns the Dirichlet groups of prime power modulus
        corresponding to primes dividing modulus.
        """
        try:
            return self._decomp
        except AttributeError:
            pass
        R = self.base_ring()
        self._decomp = [DirichletGroup(p**r,R) for p, r \
                           in arith.factor(self.modulus())]
        return self._decomp

    def exponent(self):
        return self._zeta_order

    def galois_orbits(self):
        """
        Return a list of the Galois orbits of Dirichlet characters
        in self.

        The Galois group is the absolute Galois group of the prime
        subfield of Frac(R).
        """
        try:
            return self._galois_orbits
        except AttributeError:
            pass
        L = list(self)
        G = []
        n = self.zeta_order()
        R = self.base_ring()
        p = R.characteristic()
        if p == 0:
            Auts = [e for e in xrange(2,n) if arith.GCD(e,n) == 1]
        else:
            # The nontrivial automorphisms in characteristic p are
            # k-th powering for
            #         k = p, p^2, ..., p^(r-1),
            # where p^r = 1 (mod n), so r is the mult order of p modulo n.
            r = rings.IntegerModRing(n)(p).multiplicative_order()
            Auts = [p**m for m in xrange(1,r)]
        while len(L) > 0:
            eps = L[0]
            orbit = [eps]
            del L[0]
            for s in Auts:
                chi = eps**s
                try:
                    L.remove(chi)
                    orbit.append(chi)
                except ValueError:
                    pass
            G.append(orbit)
        self._galois_orbits = G
        return G

    def gen(self, n=0):
        n = int(n)
        g = self.gens()
        if n<0 or n>=len(g):
            raise IndexError, "n(=%s) must be between 0 and %s"%(n,len(g)-1)
        return g[n]

    def gens(self):
        """
        Returns generators of self.

        EXAMPLES:
            sage: G = DirichletGroup(20)
            sage: G.gens()
            ([-1, 1], [1, zeta_4])
        """
        try:
            return self._gens
        except AttributeError:
            pass
        self._gens = []
        ug = self.unit_gens()
        R = self.base_ring()
        one = [R(1) for i in range(len(ug))]
        zeta = self.zeta()
        ord = self.zeta_order()
        for i in range(len(ug)):
            vals = list(one)
            vals[i] = zeta**(ord//arith.GCD(ord,ug[i].multiplicative_order()))
            self._gens.append(self(vals))
        self._gens = tuple(self._gens)
        return self._gens

    def integers_mod(self):
        r"""
        Returns the group of integers $\Z/N\Z$ where $N$ is the
        modulus of self.
        """
        return self._integers

    def modulus(self):
        """
        Returns the modulus of self.
        """
        return self._modulus

    def ngens(self):
        """
        Returns the number of generators of self.
        """
        return len(self.gens())

    def order(self):
        try:
            return self._order
        except AttributeError:
            pass
        self._order = rings.Integer(1)
        for g in self.gens():
            self._order *= int(g.order())
        return self._order

    def random_element(self):
        e = self(1)
        for i in range(self.ngens()):
            g = self.gen(i)
            n = random.randrange(g.order())
            e *= g**n
        return e

    def unit_gens(self):
        r"""
        Returns the minimal generators for the units of $(\Z/N\Z)^*$,
        where $N$ is the modulus of self.
        """
        return self._integers.unit_gens()

    def zeta(self):
        """
        Returns the chosen root zeta of unity in the base ring $R$.
        """
        return self._zeta

    def zeta_order(self):
        """
        Returns the order of the chosen root zeta of unity in the base
        ring $R$.
        """
        return self._zeta_order




