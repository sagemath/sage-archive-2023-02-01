

class ModularFormElement:
    """
    An element of a space of modular forms.
    """
    def __init__(self, parent, vector, check=True):
        """
        INPUT:
            parent -- ModularForms (an ambient space of modular forms)
            vector -- a vector on the basis for parent
        OUTPUT:
            ModularFormElement -- a modular form

        EXAMPLES:
            sage: M = ModularForms(Gamma0(11),2)
            sage: f = M.0
            sage: f.parent()
            Space of modular forms on Congruence Subgroup Gamma0(11) of weight 2 and dimension 2 over Rational Field
        """

        if not isinstance(parent, ModularForms):
            raise TypeError, "First argument must be an ambient space of modular forms."
        self.__parent = parent
        if not isinstance(vector, free_module_element.FreeModuleElement):
            raise TypeError, "Second argument must be a vector."
        self.__vector = vector
        if not vector in parent.vector_space():
            raise ArithmeticError, "Vector must be an element of the " +\
                      " vector space of the ambient space."

    def __ensure_is_compatible(self, other):
        if not isinstance(other, ModularFormElement):
            raise TypeError, "Second argument must be a modular form."
        if self.ambient_space() != other.ambient_space():
            raise ArithmeticError, "Modular forms must be in the same ambient space."

    def __add__(self, other):
        self.__ensure_is_compatible(other)
        return ModularFormElement(self.ambient_space(), self.vector() + other.vector())

    def __radd__(self, left):
        return self+left

    def __cmp__(self, other):
        self.__ensure_is_compatible(other)
        return self.__vector == other.__vector

    def coefficients(self, X):
        """
        The coefficients a_n of self, for integers n>=0 in the list X.

        This function caches the results of the compute function.
        """
        try:
            self.__coefficients
        except AttributeError:
            self.__coefficients = {}
        Y = [n for n in X   if    not (n in self.__coefficients.keys())]
        v = self.compute(Y)
        for i in range(len(v)):
            self.__coefficients[X[i]] = v[i]
        return v

    def compute(self, X):
        """
        Compute the coefficients a_n of self, for integers n>=0 in the list X.
        The results need not be cached; use the coefficients method instead
        for cached results.
        """
        return self.__parent.compute_coefficients(self.__vector, X)

    def __getitem__(self, n):
        return self.coefficients([n])[0]

    def __getslice__(self, i, j):
        return self.coefficients(range(i,j))

    def __neg__(self):
        return (-1)*self

    def __pos__(self):
        return self

    def __repr__(self):
        return str(self.qexp())

    def _latex_(self):
        return latex(self.qexp())

    def __sub__(self, right):
        return self + (-right)

    def __rsub__(self, left):
        raise TypeError

    def base_field(self):
        return self.parent().base_field()

    def character(self):
        chi = self.parent().character()
        if chi == None:
            raise NotImplementedError, "Determination of character in this " + \
                  "case not implemented yet."
        return chi

    def is_zero(self):
        return self.vector().is_zero()

    def level(self):
        return self.parent().level()

    def parent(self):
        return self.__parent

    def prec(self):
        try:
            self.__qexp
        except AttributeError:
            return self.parent().prec()
        return self.__qexp.prec()

    def qexp(self, prec=None):
        r"""
        The $q$-expansion of the modular form to precision $O(q^\text{prec})$.
        This function takes one argument, which is the integer prec.
        """
        if prec == None: prec = self.parent().prec()

        R = rings.PowerSeriesRing(self.base_field(), name='q')
        q = R.gen(0)

        try:
            self.__qexp
        except AttributeError:
            self.__qexp = R(0, prec=0)

        pr = self.__qexp.prec()
        if prec > pr:
            v = self.compute(range(pr,prec))
            self.__qexp = R(list(self.__qexp) + v, prec=prec)

        return self.__qexp.O(prec)

    def weight(self):
        return self.parent().weight()

    def valuation(self):
        if self.__valuation != None:
            return self.__valuation
        v = self.qexp().valuation()
        if v != rings.infinity:
            self.__valuation = v
            return v
        v = self.qexp(self.parent().sturm_bound()).valuation()
        self.__valuation = v
        return v


class EisensteinSeries(ModularFormElement):
    def __init__(self, parent, vector, t, chi, psi):
        N = parent.level()
        K = parent.base_field()
        if chi.parent().modulus() != N or psi.parent().modulus() != N:
            raise ArithmeticError, "Incompatible moduli"
        if chi.parent().base_ring() != K or psi.parent().base_ring() != K:
            raise ArithmeticError, "Incompatible base rings"
        t = int(t)
        #if not isinstance(t, int): raise TypeError, "weight must be an int"
        if parent.weight() == 2 and chi.is_trivial() \
               and psi.is_trivial() and t==1:
            raise ArithmeticError, "If chi and psi are trivial and k=2, then t must be >1."
        ModularFormElement.__init__(self, parent, vector)
        self.__chi = chi
        self.__psi = psi
        self.__t   = t

    def compute(self, X):
        """
        Compute the coefficients of $q^n$ of the power series of self,
        for $n$ in the list $X$.  The results are not cached.  (Use
        coefficients for cached results).
        """
        if self.weight() == 2 and (self.__chi.is_trivial() and self.__psi.is_trivial()):
            return self.__compute_weight2_trivial_character(X)
        else: # general case
            return self.__compute_general_case(X)

    def __compute_weight2_trivial_character(self, X):
        """
        Compute E_2 - t*E_2(q^t)
        """
        F = self.base_field()
        v = []
        t = self.__t
        for n in X:
            if n <= 0:
                v.append(F(t-1)/F(24))
            else:
                an = arith.sigma(n,1)
                if n%t==0:
                    an -= t*arith.sigma(n/t,1)
                v.append(an)
        return v


    def __compute_general_case(self, X):
        """
        General case (except weight 2, trivial character, where this is wrong!)
        $\chi$ is a primitive character of conductor $L$
        $\psi$ is a primitive character of conductor $M$
        We have $MLt \mid N$, and
        $$
          E_k(chi,psi,t) =
           c_0 + sum_{m \geq 1}[sum_{n|m}psi(n)n^{k-1}chi(m/n)] q^{mt},
        $$
        with $c_0=0$ if $L>1$,
         and
        $c_0=L(1-k,psi)/2$ if $L=1$ (that second $L$ is an $L$-function $L$).
        """
        c0, chi, psi, K, n, t, L, M = self.__defining_parameters()
        zero = K(0)
        k = self.weight()
        v = [c0]
        for i in X:
            if i==0: continue
            if i%t != 0:
                v.append(zero)
            else:
                m = i/t
                v.append(sum([psi(n)*chi(m/n)*n**(k-1) for n in arith.divisors(m)]))
        return v

    def __defining_parameters(self):
        try:
            return self.__defining_params
        except AttributeError:
            k = self.weight()
            chi = self.__chi.minimize_base_ring()
            psi = self.__psi.minimize_base_ring()
            n = arith.LCM(chi.base_ring().zeta().multiplicative_order(),\
                          psi.base_ring().zeta().multiplicative_order())
            K = rings.CyclotomicField(n)
            chi = chi.change_ring(K)
            psi = psi.change_ring(K)
            t = self.__t
            L = chi.conductor()
            M = psi.conductor()
            if L == 1:
                c0 = K(-psi.bernoulli(k))/K(2*k)
            else:
                c0 = K(0)
            self.__defining_params = (c0, chi, psi, K, n, t, L, M)
        return self.__defining_params

    def chi(self):
        return self.__chi

    def psi(self):
        return self.__psi

    def t(self):
        return self.__t

    def L(self):
        return self.__chi.conductor()

    def M(self):
        return self.__psi.conductor()

    def character(self):
        try:
            return self.__character
        except AttributeError:
            self.__character = self.__chi * (~self.__psi)
        return self.__character

    def level_at_which_new(self):
        if self.__chi.is_trivial() and self.__psi.is_trivial() and self.weight() == 2:
            return arith.factor(self.__t)[0][0]
        return self.L()*self.M()


def __find_eisen_chars(character, k):
    N = character.modulus()
    if character.is_trivial():
        V = [(character, character, t) for t in arith.divisors(N) if t>1]
        if k != 2:
            V.insert(0,(character, character, 1))
        return V

    eps = character
    if eps(-1) != (-1)**k:
        return []
    eps = eps.maximize_base_ring()
    G = eps.parent()

    # Find all pairs chi, psi such that:
    #
    #  (1) cond(chi)*cond(psi) divides the level, and
    #
    #  (2) chi == eps*psi, where eps is the nebentypus character of self.
    #
    # See [Miyake, Modular Forms] Lemma 7.1.1.

    K = G.base_ring()
    C = {}

    t0 = misc.cputime()

    for e in G:
        m = rings.Integer(e.conductor())
        if C.has_key(m):
            C[m].append(e)
        else:
            C[m] = [e]

    misc.verbose("Enumeration with conductors.",t0)

    params = []
    for L in arith.divisors(N):
        misc.verbose("divisor %s"%L)
        if not C.has_key(L):
            continue
        GL = C[L]
        for R in arith.divisors(N/L):
            if not C.has_key(R):
                continue
            GR = C[R]
            for chi in GL:
                for psi in GR:
                    if chi == eps*psi:
                        for t in arith.divisors(N/(R*L)):
                            params.append( (chi,psi,t) )
    return params


def __find_eisen_chars_gamma1(N, k):
    pairs = []
    s = (-1)**k
    G = dirichlet.DirichletGroup(N)
    E = list(G)
    parity = [c(-1) for c in E]
    for i in range(len(E)):
        for j in range(i,len(E)):
            if parity[i]*parity[j] == s and N % (E[i].conductor()*E[j].conductor()) == 0:
                pairs.append((E[i],E[j]))
                if i!=j: pairs.append((E[j],E[i]))
        #endfors
    #end if

    triples = []
    D = arith.divisors(N)
    for chi, psi in pairs:
        c_chi = chi.conductor()
        c_psi = psi.conductor()
        D = arith.divisors(N/(c_chi * c_psi))
        if (k==2 and chi.is_trivial() and psi.is_trivial()):
            D.remove(1)
        for t in D:
            triples.append((chi, psi, t))
    return triples



def compute_eisenstein_params(character, k):
    """
    Compute and return a list of all parameters $(\chi,\psi,t)$ that
    define the Eisenstein series with given character and weight $k$.

    Only the parity of $k$ is relevant.

    If character is an integer $N$, then the parameters for
    $\Gamma_1(N)$ are computed instead.  Then the condition is that
    $\chi(-1)*\psi(-1) =(-1)^k$.
    """

    if isinstance(character, (int,long)):
        N = character
        character = None
    else:
        N = character.modulus()

    if character != None:
        return __find_eisen_chars(character, k)
    else:
        return __find_eisen_chars_gamma1(N, k)
