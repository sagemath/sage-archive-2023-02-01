"""
Generic spaces of modular forms
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

class ModularFormsSpace(hecke.HeckeModule_generic):
    """
    A generic space of modular forms.
    """
    def __init__(self, group, weight, character, base_field):
        if not isinstance(group, congroup.CongruenceSubgroup):
            raise TypeError, "group (=%s) must be a congruence subroup"%group
        weight = int(weight)
        #if not isinstance(weight, int):
        #    raise TypeError, "weight must be an int"
        if not ((character is None) or isinstance(character, dirichlet.DirichletCharacter)):
            raise TypeError, "character must be a Dirichlet character"
        if not isinstance(base_field, rings.Ring):
            raise TypeError, "base_field must be a ring"
        if not base_field.is_field():
            raise ArithmeticError, "base_field must be a field."
        self.__weight, self.__group, self.__character, self.__base_field = \
                      weight, group, character, base_field
        hecke.HeckeModule_generic.__init__(self, base_field, group.level())

    def ambient_space(self):
        raise NotImplementedError   # do not implement this

    def change_ring(self):
        raise NotImplementedError

    def weight(self): return self.__weight

    def group(self): return self.__group

    def character(self): return self.__character

    def base_field(self): return self.__base_field

    def has_character(self): return self.character() != None

    def is_ambient(self): raise NotImplementedError

    def __add__(self, right):
        if self.ambient_space() != right.ambient_space():
            raise ArithmeticError, ("Intersection of %s and %s not defined because they " + \
                                    "do not lie in a common ambient space.")%\
                                   (self, right)
        if self.is_ambient(): return self
        if right.is_ambient(): return right
        V = self.vector_space() + right.vector_space()
        return ModularFormsSubmodule(self.ambient_space(), V)


    def __and__(self, right):
        return self.intersect(right)

    def __call__(self, x, check=True):
        if is_instance(x, ModularForm):
            if x.parent() == self:
                return x
            if not check:
                f = x.copy()
                f.set_parent(self)
                return f
            assert NotImplementedError
        return ModularFormElement(self, x, check=check)

    def __cmp__(self, x):
        if not isinstance(x, ModularFormsSpace):
            return -1
        if self.is_ambient() or x.is_ambient():
            if not (self.is_ambient() and x.is_ambient()): return -1
            if (self.__group, self.__weight, self.__character, self.__base_field) == \
               (x.__group, x.__weight, x.__character, x.__base_field):
                return 0
            else:
                return -1
        if self.vector_space() != x.vector_space():
            return -1
        return 0

    def __contains__(self, x):
        """
        True if x is an element or submodule of self.
        """
        if self.is_ambient() and x.is_ambient():
            return self.key() == x.key()
        raise NotImplementedError

    def __create_newspace(self, basis, level, t, is_cuspidal):
        V = self.vector_space().submodule(basis, check=False)
        S = ModularForms(self.ambient_space(), V)
        S.__newspace_params = {'level':level, 't':t}
        S.__is_cuspidal = is_cuspidal
        S.__is_eisenstein = not is_cuspidal
        return S

    def __newspace_bases(self):
        if hasattr(self, "__newspace_bases_list"):
            return self.__newspace_bases_list
        assert self.is_ambient()
        V = self.vector_space()
        eps, k, N = self.__character, self.__weight, self.__level
        # First the cuspidal new spaces.
        m = eps.conductor()
        levels = [M for M in arith.divisors(N) if M%m==0]
        levels.reverse()
        B = []; i = 0
        for M in levels:
            n = dims.dimension_new_cusp_forms(eps.restrict(M), k)
            for t in arith.divisors(N/M):
                basis = [V.gen(i+j) for j in range(n)]
                print M, basis
                if len(basis) > 0:
                    B.append((M, t, True, basis))
                i += n
        # Now the Eisenstein series
        #x = [0 for _ in range(len(levels))]
        x = {}
        for E in self.eisenstein_series():  # count number of e.s. of each level
            Mt = (E.new_level(), E.t())
            if not x.has_key(Mt):
                x[Mt] = 1
            else:
                x[Mt] += 1
        k = x.keys()
        k.sort()
        k.reverse()
        for M, t in k:
            n = x[(M,t)]
            B.append((M, t, False, [V.gen(i+j) for j in range(n)]))
            i += n
        self.__newspace_bases_list = B
        return self.__newspace_bases_list

    def __submodule_from_subset_of_basis(self, x):
        V = self.vector_space()
        return V.submodule([V.gen(i) for i in x], check=False)

    def basis(self):
        if hasattr(self, "__basis"): return self.__basis
        self.__basis = tuple([ModularFormElement(self, x, check=False) for \
                        x in self.vector_space().basis()])
        return self.__basis

    def gen(self, n):
        return self.basis()[n]

    def gens(self):
        return self.basis()

    def sturm_bound(self, M=None):
        r"""
        For a space M of modular forms, this function returns an integer B
        such that two modular forms in either self or M are equal if and only
        if their q-expansions are equal to precision B.  If M is none, then
        M is set equal to self.

        NOTES:
        Reference for the Sturm bound that we use in the definition of
        of this function:

         J. Sturm, On the congruence of modular forms,
              Number theory (New York, 1984--1985), Springer,
              Berlin, 1987, pp.~275--280.

        Useful Remark:

            Kevin Buzzard pointed out to me (William Stein) in Fall
            2002 that the above bound is fine for Gamma1 with
            character, as one sees by taking a power of $f$.  More
            precisely, if $f\con 0\pmod{p}$ for first $s$
            coefficients, then $f^r = 0 \pmod{p}$ for first $s r$
            coefficents.  Since the weight of $f^r$ is $r
            \text{weight}(f)$, it follows that if $s \geq $ the sturm
            bound for $\Gamma_0$ at weight(f), then $f^r$ has
            valuation large enough to be forced to be $0$ at $r\cdot$
            weight(f) by Sturm bound (which is valid if we choose $r$
            right).  Thus $f \con 0 \pmod{p}$.  Conclusion: For
            $\Gamma_1$ with fixed character, the Sturm bound is
            \emph{exactly} the same as for $\Gamma_0$.  A key point is
            that we are finding $\Z[\eps]$ generators for the Hecke
            algebra here, not $\Z$-generators.  So if one wants
            generators for the Hecke algebra over $\Z$, this bound is
            wrong.

            This bound works over any base, even a finite field.
            There might be much better bounds over $\Q$, or for
            comparing two eigenforms.
        """
        if M != None:
            raise NotImplementedError
        if self.__sturm_bound == None:
            # the +1 below is because O(q^prec) has precision prec.
            self.__sturm_bound = int(\
                math.ceil(self.weight()*dims.idxG0(self.level())/12.0) + 1)
        return self.__sturm_bound

    def character(self):
        return self.__character

    def cuspidal_submodule(self):
        if self.__is_cuspidal == True:
            return self
        if self.__cuspidal_submodule != None:
            return self.__cuspidal_submodule
        if self.is_ambient():
            # By definition the cuspidal submodule of the ambient space
            # is spanned by the first n standard basis vectors, where
            # n is the dimension of the cuspidal submodule.
            n = self.__ambient_cusp_dimension()
            W = self.__submodule_from_subset_of_basis(range(n))
            S = ModularForms(self, W)
            S.__is_cuspidal = True
            S.__is_eisenstein = (n==0)
            self.__cuspidal_submodule = S
            return S
        C = self.ambient_space().cuspidal_submodule()
        S = self.intersect(C)
        if S.dimension() < self.dimension():
            self.__is_cuspidal = False
            self.__cuspidal_submodule = S
        else:
            assert S.dimension() == self.dimension()
            self.__is_cuspidal = True
        S.__is_eisenstein = (S.dimension()==0)

    def decomposition(self):
        """

        This function returns a list of submodules $V(f_i,t)$
        corresponding to newforms $f_i$ of some level dividing the
        level of self, such that the direct sum of the submodules
        equals self, if possible.  The space $V(f_i,t)$ is the image
        under $g(q)$ maps to $g(q^t)$ of the intersection with
        $R[[q]]$ of the space spanned by the conjugates of $f_i$,
        where $R$ is the base ring of self.

        """
        raise NotImplementedError

    def newspaces(self):
        r"""
        This function returns a list of submodules $S(M,t)$ and
        $E(M,t)$, corresponding to levels $M$ dividing $N$ and integers $t$
        dividing $N/M$, such that self is the direct sum of these
        spaces, if possible.  Here $S(M,t)$ is by definition
        the image under $f(q) \mapsto f(q^t)$ of the new submodule of
        cusp forms of level $M$, and similarly $E(M,t)$ is the image of
        Eisenstein series.

        Notes: (1) the submodules $S(M,t)$ need not be stable under
        Hecke operators of index dividing $N/M$.  (2) Since self can
        be an arbitrary submodule, there's no guarantee any $S(M,t)$ or
        $E(M,t)$ is in self, so the return list could be empty.
        """
        V = self.embedded_submodule()
        return [self.__create_newspace(basis=B,level=M,t=t,is_cuspidal=is_cuspidal) \
                for M, t, is_cuspidal, B in self.ambient_space().__newspace_bases() \
                if V.contains_each(B)]


    def eisenstein_submodule(self):
        if self.__is_eisenstein == True:
            return self
        if self.__eisenstein_submodule != None:
            return self.__eisenstein_submodule
        if self.is_ambient():
            # By definition the eisenstein submodule of the ambient space
            # is spanned by the n+1 through n+d standard basis vectors, where
            # n is the dimension of the cuspidal submodule and d
            # is the dimension of the eisenstein submodule (i.e., the
            # number of eisenstein series).
            n = self.__ambient_cusp_dimension()
            d = self.__ambient_eis_dimension()
            W = self.__submodule_from_subset_of_basis(range(n,n+d))
            E = ModularForms(self, W)
            E.__is_eisenstein = True
            E.__is_cuspidal = (d==0)
            self.__eisenstein_submodule = E
            return E
        A = self.ambient_space().eisenstein_submodule()
        E = self.intersect(A)
        if E.dimension() < self.dimension():
            self.__is_eisenstein = False
            self.__eisenstein_submodule = E
        else:
            assert E.dimension() == self.dimension()
            self.__is_eisenstein = True
        E.__is_cuspidal = (E.dimension()==0)

    def embedded_submodule(self):
        if self.is_ambient():
            return self.vector_space()
        return self.__embedded_submodule

    def hecke_matrix(self, n):
        raise NotImplementedError

    def intersect(self, right):
        if self.ambient_space() != right.ambient_space():
            raise ArithmeticError, "Intersection of %s and %s not defined."%\
                                   (self, right)
        V = self.embedded_submodule().intersect(right.embedded_submodule())
        return ModularForms(self.ambient_space(),V)

    def is_ambient(self):
        return self.__ambient == None

    def key(self):
        if self.is_ambient():
            return self.__key
        return self.__ambient

    def level(self):
        return self.group().level()

    def modular_symbols(self):
        raise NotImplementedError

