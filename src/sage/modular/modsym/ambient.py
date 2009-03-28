"""
Ambient spaces of modular symbols.

EXAMPLES:

We compute a space of modular symbols modulo 2. The dimension is
different than that of the corresponding space in characteristic
0::

    sage: M = ModularSymbols(11,4,base_ring=GF(2)); M
    Modular Symbols space of dimension 7 for Gamma_0(11) of weight 4
    with sign 0 over Finite Field of size 2
    sage: M.basis()
    ([X*Y,(1,0)], [X*Y,(1,8)], [X*Y,(1,9)], [X^2,(0,1)], [X^2,(1,8)], [X^2,(1,9)], [X^2,(1,10)])
    sage: M0 =  ModularSymbols(11,4,base_ring=QQ); M0
    Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4
    with sign 0 over Rational Field
    sage: M0.basis()
    ([X^2,(0,1)], [X^2,(1,6)], [X^2,(1,7)], [X^2,(1,8)], [X^2,(1,9)], [X^2,(1,10)])

The charpoly of the Hecke operator `T_2` has an extra
factor `x`.

::

    sage: M.T(2).matrix().fcp('x')
    (x + 1)^2 * x^5
    sage: M0.T(2).matrix().fcp('x')
    (x - 9)^2 * (x^2 - 2*x - 2)^2
"""

################################################################################
#       SAGE: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
################################################################################

# SAGE packages
from   sage.misc.search import search
import sage.misc.latex as latex
import sage.misc.misc as misc

import sage.matrix.matrix_space as matrix_space
import sage.modules.free_module_element as free_module_element
import sage.modules.free_module as free_module
import sage.misc.misc as misc
import sage.modular.congroup as congroup
import sage.modular.dims as dims
import sage.modular.dirichlet as dirichlet
import sage.modular.hecke.all as hecke
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.all as rings
import sage.rings.arith as arith
import sage.structure.formal_sum as formal_sum
import sage.categories.all as cat
from sage.modular.cusps import Cusp
from sage.rings.arith import binomial
import sage.structure.all

import boundary
import element
import heilbronn
import manin_symbols
import modular_symbols
import modsym
import p1list
import relation_matrix
import space
import subspace

QQ = rings.Rational
ZZ = rings.Integers


class ModularSymbolsAmbient(space.ModularSymbolsSpace, hecke.AmbientHeckeModule):
    """
    An ambient space of modular symbols for a congruence subgroup of
    SL_2(Z).

    This class is an abstract base class, so only derived classes
    should be instantiated. INPUT:


    -  ``weight`` - an integer = 2

    -  ``group`` - a congruence subgroup.

    -  ``sign`` - an integer, either -1, 0, or 1

    -  ``base_ring`` - a commutative ring
    """
    def __init__(self, group, weight, sign, base_ring,
                 character = None):
        """
        Initialize a space of modular symbols.
        """
        weight = int(weight)
        if weight <= 1:
            raise ValueError, "Weight (=%s) Modular symbols of weight <= 1 not defined."%weight
        if not congroup.is_CongruenceSubgroup(group):
            raise TypeError, "group must be a congruence subgroup"

        sign = int(sign)
        if not isinstance(base_ring, rings.Ring) and base_ring.is_field():
            raise TypeError, "base_ring must be a commutative ring"

        if character == None and congroup.is_Gamma0(group):
            character = dirichlet.TrivialCharacter(group.level(), base_ring)

        space.ModularSymbolsSpace.__init__(self, group, weight,
                                           character, sign, base_ring)

        try:
            formula = self._dimension_formula()
        except NotImplementedError:
            formula = None
        rank = self.rank()
        if formula != None:
            assert rank == formula, \
                   "Computed dimension (=%s) of ambient space \"%s\" doesn't match dimension formula (=%s)!\n"%(d, self, formula) + \
                   "ModularSymbolsAmbient: group = %s, weight = %s, sign = %s, base_ring = %s, character = %s"%(
                         group, weight, sign, base_ring, character)

        hecke.AmbientHeckeModule.__init__(self, base_ring, rank, group.level(), weight)

    def __cmp__(self, other):
        if not isinstance(other, space.ModularSymbolsSpace):
            return cmp(type(self), type(other))
        if isinstance(other, ModularSymbolsAmbient):
            return misc.cmp_props(self, other, ['group', 'weight', 'sign', 'base_ring', 'character'])
        c = cmp(self, other.ambient_hecke_module())
        if c: return c
        if self.free_module() == other.free_module():
            return 0
        return -1

    def new_submodule(self, p=None):
        """
        Returns the new or p-new submodule of self.

        INPUT:


        -  ``p`` - (default: None); if not None, return only
           the p-new submodule.


        OUTPUT: the new or p-new submodule of self

        EXAMPLES::

            sage: ModularSymbols(100).new_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 31 for Gamma_0(100) of weight 2 with sign 0 over Rational Field
            sage: ModularSymbols(389).new_submodule()
            Modular Symbols space of dimension 65 for Gamma_0(389) of weight 2 with sign 0 over Rational Field
        """
        # Check for special cases where the answer is easy.
        # If not in one of those cases, use the generic code.
        if self.level().is_prime() and self.weight() == 2:
            return self
        return hecke.AmbientHeckeModule.new_submodule(self, p=p)

    def manin_symbols(self):
        raise NotImplementedError

    def manin_generators(self):
        """
        Return list of all Manin symbols for this space. These are the
        generators in the presentation of this space by Manin symbols.

        EXAMPLES::

            sage: M = ModularSymbols(2,2)
            sage: M.manin_generators()
            [(0,1), (1,0), (1,1)]

        ::

            sage: M = ModularSymbols(1,6)
            sage: M.manin_generators()
            [[Y^4,(0,0)], [X*Y^3,(0,0)], [X^2*Y^2,(0,0)], [X^3*Y,(0,0)], [X^4,(0,0)]]
        """
        return self._manin_generators

    def manin_basis(self):
        r"""
        Return a list of indices into the list of Manin generators (see
        ``self.manin_generators()``) such that those symbols
        form a basis for the quotient of the `\mathbb{Q}`-vector
        space spanned by Manin symbols modulo the relations.

        EXAMPLES::

            sage: M = ModularSymbols(2,2)
            sage: M.manin_basis()
            [1]
            sage: [M.manin_generators()[i] for i in M.manin_basis()]
            [(1,0)]
            sage: M = ModularSymbols(6,2)
            sage: M.manin_basis()
            [1, 10, 11]
            sage: [M.manin_generators()[i] for i in M.manin_basis()]
            [(1,0), (3,1), (3,2)]
        """
        try:
            return self._manin_basis
        except AttributeError:
            self.compute_presentation()
        return self._manin_basis

    def p1list(self):
        try:
            return self.__p1list
        except AttributeError:
            self.__p1list = p1list.P1List(self.level())
        return self.__p1list

    def relation_matrix(self):
        raise NotImplementedError

    def compute_presentation(self):
        B, basis, mod = relation_matrix.compute_presentation(
                self.manin_symbols(), self.sign(),
                self.base_ring(), self.weight())
        self._manin_generators = self.manin_symbols().manin_symbol_list()
        self._manin_basis = basis
        self._manin_gens_to_basis = B
        self._mod2term = mod

    def manin_gens_to_basis(self):
        try:
            return self._manin_gens_to_basis
        except AttributeError:
            self.compute_presentation()
            return self._manin_gens_to_basis


    #####################################################################
    # Coercion
    #####################################################################
    def __call__(self, x, computed_with_hecke=False):
        r"""
        Coerce x into this modular symbols space (self). The result is
        either an element of self or a subspace of self.

        The allowed inputs for x are as follows:


        -  ``Vector`` - a vector of the same degree. This
           defines the corresponding linear combination of the basis of self.

        -  ``ManinSymbol`` - a Manin symbol of the same weight
           as the space

        -  ``ModularSymbolsElement`` - a modular symbol whose
           ambient parent is this space of modular symbols. (TODO: make more
           sophisticated)

        -  0 - the integer 0; results in the 0 modular symbol.

        -  3-tuple - Given a 3-tuple (i,u,v), returns the modular symbol
           element defined by the Manin symbol
           `[X^{i}\cdot Y^{k-2-i}, (u,v)]`, where k is the weight.
           Note that we must have `0\leq i \leq k-2`.

        -  2-tuple - Given a 2-tuple (u,v), returns the element defined by
           the Manin symbol `[X^0 \cdot Y^{2-k}, (u,v)]`.

        -  2-elements list - Given a list ``[alpha, beta]``,
           where `\alpha` and `\beta` are (coercible to)
           cusps, return the modular symbol `\{\alpha, \beta\}`. When
           the the weight `k > 2` return
           `Y^{k-2} \{\alpha, \beta\}`.

        -  3-element list - Given a list ``[i, alpha, beta]``,
           where `i` is an integer, and `\alpha`,
           `\beta` are (coercible to) cusps, return the modular symbol
           `X^i Y^{k-2-i} \{\alpha, \beta\}`.

           If our list is ``[f, alpha, beta]``, where `f`
           is a homogeneous polynomial in two variables of degree k-2 with
           integer coefficients, and alpha and beta are cusps, return the
           corresponding sum of modular symbols as an element of self. So if
           `f = \sum_{i=0}^{k-2} a_i X^i Y^{k-2-i}`, return
           `\sum_{i=0}^{k-2} a_i * [ i, alpha, beta ]`.
        """
        if isinstance(x, free_module_element.FreeModuleElement):
            if x.degree() != self.dimension():
                raise TypeError, "Incompatible degrees: x has degree %s but modular symbols space has dimension %s"%(
                    x.degree(), self.dimension())
            #if x.parent().base_ring() != self.base_ring():
            #    raise TypeError, "Vector x is over %s, but modular symbols space is over %s."%(
            #        x.parent().base_ring(), self.base_ring())
            return element.ModularSymbolsElement(self, x)

        elif isinstance(x, (manin_symbols.ManinSymbol, element.ModularSymbolsElement)):
            return self.element(x)

        elif isinstance(x, modular_symbols.ModularSymbol):
            return self(x.manin_symbol_rep())

        elif isinstance(x, (int, rings.Integer)) and x==0:
            return element.ModularSymbolsElement(self, self.free_module()(0))

        elif isinstance(x, tuple):
            return self.manin_symbol(x)

        elif isinstance(x, formal_sum.FormalSum):
            return sum([c*self(y) for c, y in x], self(0))

        elif isinstance(x, list):
            if len(x) == 3 and rings.is_MPolynomial(x[0]):
                return self.modular_symbol_sum(x)
            else:
                return self.modular_symbol(x)

        raise TypeError, "No coercion of %s into %s defined."%(x, self)

    def _action_on_modular_symbols(self, g):
        """
        Compute the matrix of the action of the 2x2 integer matrix
        g=[a,b,c,d] (which must be specified as an integer list) on self
        with respect to the standard basis.

        Use _matrix_of_operator_on_modular_symbols for more general
        operators.
        """
        if not isinstance(g, list):
            raise TypeError, "g must be a list"
        if not len(g) == 4:
            raise TypeError, "g must be a list of length 4"
        return self._matrix_of_operator_on_modular_symbols(self, [g])

    def manin_symbol(self, x, check=True):
        if check:
            if len(x) == 2:
                x = (0,x[0],x[1])
            if len(x) == 3:
                # Manin symbol of the form (i, u, v), which corresponds to [X^i*Y^(k-2-i), (u,v)].
                if x[0] < 0 or x[0] > self.weight()-2:
                    raise ValueError, "The first entry of the tuple (=%s) must be an integer between 0 and k-2 (=%s)."%(
                        x, self.weight()-2)
            else:
                raise ValueError, "x (=%s) must be of length 2 or 3"%x
        # end check

        N = self.level()
        x = (x[0], x[1]%N, x[2]%N)
        try:
            return self.__manin_symbol[x]
        except AttributeError:
            self.__manin_symbol = {}
        except KeyError:
            pass
        y = manin_symbols.ManinSymbol(self.manin_symbols(), x)
        z = self(y)
        self.__manin_symbol[x] = z
        return z

    def _modular_symbol_0_to_alpha(self, alpha, i=0):
        if alpha.is_infinity():
            return self.manin_symbol((i,0,1), check=False)
        v, c = arith.continued_fraction_list(alpha._rational_(), partial_convergents=True)
        a = self(0)
        zero = rings.ZZ(0)
        one = rings.ZZ(1)
        two = rings.ZZ(2)
        if self.weight() > two:
            R = rings.ZZ['X']
            X = R.gen(0)
            ## need to add first two terms, which aren't necessarily
            ## zero in this case. we do the first here, and the
            ## second in the k=0 case below, so as to avoid code
            ## duplication
            a += self.manin_symbol((i,0,1), check=False)
            for k in range(0,len(c)):
                ## matrix entries associated to this partial sum
                if k == 0:
                    x = c[0][0]
                    y = -1
                    z = 1
                    w = 0
                else:
                    x = c[k][0]
                    y = c[k-1][0]
                    z = c[k][1]
                    w = c[k-1][1]
                    if k%2 == 0:
                        y = -y
                        w = -w

                ## two options here: write out the polynomial directly,
                ## and deal with all the separate cases, or create two
                ## polynomials and then exponentiate and multiply them.
                ## given how fast ntl/flint/etc are, the second may
                ## be faster.

                ## method 1: write out solution. this is currently
                ## incorrect, because it ends up doing 0^0 in the sum,
                ## so i'll fix it and do timings soon.
##                for s in range(0,self.weight()-two+1):
##                    coeff = sum([ binomial(i,t)*binomial(self.weight()-two-i,s-t)*
##                                  x**t * y**(i-t) * z**(s-t) *
##                                  w**(self.weight()-two-i-s+t) for t in range(0,s) ])
##                    m = coeff * self.manin_symbol((s, y, w), check=False)
##                    a += m

                ## method 2
                p1 = x*X+y
                p2 = z*X+w
                if i == 0:
                    p1 = R(one)
                if (self.weight()-2-i == 0):
                    p2 = R(one)
                poly = (p1**i) * (p2**(self.weight()-2-i))
                for s in range(0,self.weight()-1): ## k-2+1 = k-1
                    a += poly[s] * self.manin_symbol((s,z,w), check=False)
        else:
            for k in range(1,len(c)):
                u = c[k][1]
                v = c[k-1][1]
                if k % 2 == 0:
                    v = -v
                x = self.manin_symbol((i, u, v), check=False)
                a += x
        return a

    def modular_symbol(self, x, check=True):
        r"""
        Create a modular symbol in this space.

        INPUT:


        -  ``x`` - a list of either 2 or 3 entries

        -  ``2 entries: [`\alpha`, `\beta`]`` -
           creates the modular symbol alpha, beta, or, if the weight is 2 the
           symbol `Y^(k-2){\alpha,\beta}`.

        -  ``3 entries: [i, `\alpha`,
           `\beta`]`` - create the modular symbol
           `X^i Y^{k-2-i}{\alpha,\beta}`.


        EXAMPLES::

            sage: set_modsym_print_mode('modular')
            sage: M = ModularSymbols(11)
            sage: M.modular_symbol([2/11, oo])
            -{-1/9,0}
            sage: M.1
            {-1/8,0}
            sage: M.modular_symbol([-1/8, 0])
            {-1/8,0}
            sage: M.modular_symbol([0, -1/8, 0])
            {-1/8,0}
            sage: M.modular_symbol([10, -1/8, 0])
            Traceback (most recent call last):
            ...
            ValueError: The first entry of the tuple (=[10, -1/8, 0]) must be an integer between 0 and k-2 (=0).

        ::

            sage: N = ModularSymbols(6,4)
            sage: set_modsym_print_mode('manin')
            sage: N([1,Cusp(-1/4),Cusp(0)])
            17/2*[X^2,(2,3)] - 9/2*[X^2,(2,5)] + 15/2*[X^2,(3,1)] - 15/2*[X^2,(3,2)]
            sage: N([1,Cusp(-1/2),Cusp(0)])
            1/2*[X^2,(2,3)] + 3/2*[X^2,(2,5)] + 3/2*[X^2,(3,1)] - 3/2*[X^2,(3,2)]

        Use check=False for efficiency if the input x is a list of length 3
        whose first entry is an Integer, and whose second and third entries
        are cusps::

            sage: M.modular_symbol([0, Cusp(2/11), Cusp(oo)], check=False)
            -(1,9)

        ::

            sage: set_modsym_print_mode()   # return to default.
        """

        if check:
            if len(x) == 2:
                x = [0,x[0],x[1]]
            elif len(x) == 3:
                if x[0] < 0 or x[0] > self.weight()-2:
                    raise ValueError, "The first entry of the tuple (=%s) must be an integer between 0 and k-2 (=%s)."%(
                        x, self.weight()-2)
            else:
                raise ValueError, "x (=%s) must be of length 2 or 3"%x
            i = rings.Integer(x[0])
            alpha = Cusp(x[1])
            beta = Cusp(x[2])
        else:
            i = x[0]
            alpha = x[1]
            beta = x[2]

        # Compute {0,beta} - {0,alpha}
        a = self._modular_symbol_0_to_alpha(alpha, i)
        b = self._modular_symbol_0_to_alpha(beta, i)
        return b - a

    def modular_symbol_sum(self, x, check=True):
        """
        Take a 3-element list consisting of a homogeneous polynomial in
        degree self.weight()-2 over ZZ, and return the corresponding sum of
        modular symbols as an element of self. That is, if x is ``[f,
        alpha, beta]``, where
        `f = \sum_{i=0}^{k-2} a_i X^i Y^{k-2-i}`, return
        `\sum_{i=0}^{k-2} a_i * [ i, alpha, beta ]`.

        EXAMPLES:
        """
        if check:
            if len(x) != 3:
                raise ValueError, "%s must have length 3"%x
            f = x[0]
            R = self.base_ring()['X','Y']
            X = R.gen(0)
            try:
                f = R(f)
            except TypeError:
                raise ValueError, \
                      "f must be coercible to a polynomial over %s"%self.base_ring()
            if (not f.is_homogeneous()) or (f.degree() != self.weight()-2):
                raise ValueError, "f must be a homogeneous polynomial of degree k-2"
            alpha = Cusp(x[1])
            beta = Cusp(x[2])
        else:
            f = x[0]
            R = self.base_ring()
            X = R.gen(0)
            alpha = x[1]
            beta = x[2]

        s = self(0)

        for term in f.monomials():
            deg = term.degree(X)
            a = self._modular_symbol_0_to_alpha(alpha, deg)
            b = self._modular_symbol_0_to_alpha(beta, deg)
            s += f.monomial_coefficient(term) * (b-a)

        return s


    def _compute_dual_hecke_matrix(self, n):
        return self.hecke_matrix(n).transpose()

    def _compute_hecke_matrix_prime(self, p, rows=None):
        """
        Compute and return the matrix of the p-th Hecke operator.

        EXAMPLES: We first compute some examples for Gamma0(N)::

            sage: m = ModularSymbols(2, weight=4)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^2 - 9*x + 8

        ::

            sage: m = ModularSymbols(1,weight=12)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^3 - 2001*x^2 - 97776*x - 1180224
            sage: m._compute_hecke_matrix_prime(13).charpoly('x')
            x^3 - 1792159238562*x^2 - 2070797989680255444*x - 598189440899986203208472

        ::

            sage: m = ModularSymbols(1,weight=12, sign=-1)
            sage: m._compute_hecke_matrix_prime(5)
            [4830]
            sage: m._compute_hecke_matrix_prime(23)
            [18643272]

        ::

            sage: m = ModularSymbols(3,4)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^2 - 18*x + 81

        ::

            sage: m = ModularSymbols(6,4)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^6 - 14*x^5 + 29*x^4 + 172*x^3 - 124*x^2 - 320*x + 256
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^6 - 50*x^5 + 511*x^4 + 3012*x^3 - 801*x^2 - 9234*x + 6561

        ::

            sage: m = ModularSymbols(15,4, sign=-1)
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^4 - 2*x^3 + 18*x^2 + 18*x - 243

        ::

            sage: m = ModularSymbols(6,4)
            sage: m._compute_hecke_matrix_prime(7).charpoly('x')
            x^6 - 1344*x^5 + 666240*x^4 - 140462080*x^3 + 8974602240*x^2 + 406424518656*x + 3584872677376

        ::

            sage: m = ModularSymbols(4,4)
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^3 - 84*x^2 + 2352*x - 21952

        We now compute some examples for modular symbols on Gamma1(N)::

            sage: m = ModularSymbols(Gamma1(13),2, sign=-1)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^2 + 3*x + 3

        The following is an example with odd weight::

            sage: m = ModularSymbols(Gamma1(5),3)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^4 - 10*x^3 + 50*x^2 - 170*x + 289

        This example has composite conductor and weight2 dividing the
        conductor and nontrivial sign::

            sage: m = ModularSymbols(Gamma1(9),3, sign=1)
            sage: m._compute_hecke_matrix_prime(3).charpoly('x')
            x^6 + 3*x^4 - 19*x^3 + 24*x^2 - 9*x
        """
        # note -- p doesn't have to be prime despite the function name
        p = int(rings.Integer(p))   # go through Integer so p = 2.5 gives an error.
        # NOTE -- it is actually NOT necessary that p be prime.
        if isinstance(rows, list):
            rows = tuple(rows)
        try:
            return self._hecke_matrices[(p,rows)]
        except AttributeError:
            self._hecke_matrices = {}
        except KeyError:
            pass
        tm = misc.verbose("Computing Hecke operator T_%s"%p)

        if arith.is_prime(p):
            H = heilbronn.HeilbronnCremona(p)
        else:
            H = heilbronn.HeilbronnMerel(p)

        B = self.manin_basis()
        if not rows is None:
            B = [B[i] for i in rows]
        cols = []
        mod2term = self._mod2term
        R = self.manin_gens_to_basis()
        K = self.base_ring()
        W = R.new_matrix(nrows=len(B), ncols = R.nrows())
        syms = self.manin_symbols()
        n = len(syms)
        j = 0
        for i in B:
            for h in H:
                entries = syms.apply(i,h)
                for k, x in entries:
                    f, s = mod2term[k]
                    if s != 0:
                        W[j,f] = W[j,f] + s*K(x)
            j += 1
        tm = misc.verbose("start matrix multiply",tm)
        if hasattr(W, '_matrix_times_matrix_dense'):
            Tp = W._matrix_times_matrix_dense(R)
            misc.verbose("done matrix multiply and computing Hecke operator",tm)
        else:
            Tp = W * R
            tm = misc.verbose("done matrix multiply",tm)
            Tp = Tp.dense_matrix()
            misc.verbose("done making Hecke operator matrix dense",tm)
        self._hecke_matrices[(p,rows)] = Tp
        return Tp


    def __heilbronn_operator(self, M, H, t=1):
        """
        Returns the matrix function from self to M defined by the pair (H,
        t), where H is a list of matrices and t is an integer.

        INPUT:


        -  ``self`` - ModularSymbols , domain (an ambient space
           of modular symbols),

        -  ``M`` - ModularSymbols, codomain (a space of modular
           symbols),

        -  ``H`` - list, a list of matrices in M_2(Z),

        -  ``t`` - int, an integer.


        OUTPUT:


        -  ``free module morphism`` - A function from self to M
           defined by t and the matrices in H.
        """

        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        hom = hecke.HeckeModuleHomspace(self, M)
        if self.dimension() == 0 or M.dimension() == 0:
            A = MS(0)
            phi = hom(A, "Heilbronn operator(%s,%s)"%(H,t))
            return phi

        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        k = self.weight()
        for n in B:
            z = M(0)
            i, u, v = syms[n]
            # We apply each Heilbronn matrix to the
            #    Manin symbol [X^i*Y^(k-2-i), (u,v)]
            for h in H:
                # Apply h to the polynomial part
                (a,b,c,d) = tuple(h)
                # P gives the ordered coefficients of (a*X+b*Y)^i*(c*X+d*Y)^(j-i)
                P = manin_symbols.apply_to_monomial(i, k-2, a,b,c,d)
                # Apply h to the (u,v) part of the Manin symbol
                (uu,vv) = (u*a+v*c, u*b+v*d)

                # For the generalized Heilbronn operator, we through away any
                # symbols for which the (u,v) part of the symbol doesn't have
                # both entries divisible by t.
                if t != 1:
                    if uu%t != 0 or vv%t != 0:
                        continue
                    uu = uu//t
                    vv = vv//t

                # Now coerce each Manin symbol
                #
                #         P[m]*[X^m*Y^(k-2-m), (uu,vv)],    for m=0,...,len(P)
                #
                # into the image space M and add that to z.
                # Note that we coerce in Manin symbols as tuples.
                for m in range(len(P)):
                    x = M((m,uu,vv))
                    z += x*P[m]

            rows.append(z.element())

        A = MS(rows)
        return hom(A, "Heilbronn operator(%s,%s)"%(H,t))

    def _repr_(self):
        return "Modular Symbols space of dimension %s and weight %s for %s with sign %s and character %s over %s"%(
                self.dimension(), self.weight(), self.group(), self.sign(), self.character(), self.base_ring())

    def _latex_(self):
        return "\\mathrm{ModSym}_{%s}(%s,%s;%s)"%(self.weight(),
                                                     latex.latex(self.group()),
                                                     latex.latex(self.character()),
                                                     latex.latex(self.base_ring()))

    def _matrix_of_operator_on_modular_symbols(self, codomain, R):
        """
        INPUT:


        -  ``self`` - this space of modular symbols

        -  ``codomain`` - space of modular symbols

        -  ``R`` - list of lists [a,b,c,d] of length 4, which
           we view as elements of GL_2(Q).


        OUTPUT: A matrix, which represents the operator

        .. math::

                            x \mapsto \sum_{g in R} g.x


        where g.x is the formal linear fractional transformation on modular
        symbols.
        """
        rows = []
        for b in self.basis():
            v = formal_sum.FormalSum(0, check=False)
            for c, x in b.modular_symbol_rep():
                for g in R:
                    y = x.apply(g)
                    v += y*c
            w = codomain(v).element()
            rows.append(w)
        M = matrix_space.MatrixSpace(self.base_ring(), len(rows), codomain.degree(), sparse=False)
        return M(rows)

    def _compute_atkin_lehner_matrix(self, d):
        r"""
        INPUT:


        -  ``d`` - integer that divides level


        OUTPUT: matrix

        EXAMPLES: An example at level 29::

            sage: M = ModularSymbols((DirichletGroup(29,QQ).0), 2,1); M
            Modular Symbols space of dimension 4 and level 29, weight 2, character [-1], sign 1, over Rational Field
            sage: w = M._compute_atkin_lehner_matrix(29)
            sage: w^2 == 1
            True
            sage: w.fcp()
            (x - 1)^2 * (x + 1)^2

        This doesn't work since the character has order 2::

            sage: M = ModularSymbols((DirichletGroup(13).0), 2,1); M
            Modular Symbols space of dimension 0 and level 13, weight 2, character [zeta12], sign 1, over Cyclotomic Field of order 12 and degree 4
            sage: M._compute_atkin_lehner_matrix(13)
            Traceback (most recent call last):
            ...
            ValueError: Atkin-Lehner only leaves space invariant when character is trivial or quadratic.  In general it sends M_k(chi) to M_k(1/chi)

        Note that Atkin-Lehner does make sense on `\Gamma_1(13)`,
        but doesn't commute with the Hecke operators::

            sage: M = ModularSymbols(Gamma1(13),2)
            sage: w = M.atkin_lehner_operator(13).matrix()
            sage: t = M.T(2).matrix()
            sage: t*w == w*t
            False
            sage: w^2 == 1
            True
        """
        chi = self.character()
        if chi is not None and chi.order() > 2:
            raise ValueError, "Atkin-Lehner only leaves space invariant when character is trivial or quadratic.  In general it sends M_k(chi) to M_k(1/chi)"

        N = self.level()
        k = self.weight()
        R = self.base_ring()
        if N%d != 0:
            raise ValueError, "d must divide N"

        g, x, y = arith.xgcd(d, -N//d)
        g = [d*x, y, N, d]
        A = self._action_on_modular_symbols(g)
        scale = R(d)**(1 - k//2)
        Wmat = scale * A
        return Wmat

    def boundary_map(self):
        """
        The boundary map to the corresponding space of boundary modular
        symbols.
        """
        try:
            return self.__boundary_map
        except AttributeError:
            # compute boundary map
            B = self.boundary_space()
            I = [B(b) for b in self.basis()]
            W = matrix_space.MatrixSpace(self.base_ring(), len(I), B.rank(), sparse=True)

            # Note -- the underlying elements have degree the number of distinct
            # cusps known when the element was computed.  This isn't constant,
            # so we pad the elements.
            E = [x.element() for x in I]
            zero = self.base_ring()(0)
            n = int(B.dimension())
            E = sum([ list(x) + [zero]*(n - len(x)) for x in E ], [])

            A = W( E )
            H = cat.Hom(self, B)
            self.__boundary_map = H(A, "boundary map")
            return self.__boundary_map

    def cusps(self):
        try:
            return self.__cusps
        except AttributeError:
            f = self.boundary_map()
            B = f.codomain()
            C = B._known_cusps()
            self.__cusps = C
            return C

    def boundary_space(self):
        raise NotImplementedError

    def cuspidal_submodule(self):
        """
        The cuspidal submodule of self.

        EXAMPLES::

            sage: M = ModularSymbols(12,2,0,GF(5)) ; M
            Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Finite Field of size 5
            sage: M.cuspidal_submodule()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Finite Field of size 5
            sage: ModularSymbols(1,24,-1).cuspidal_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 24 with sign -1 over Rational Field

        The cuspidal submodule of the cuspidal submodule is itself::

            sage: M = ModularSymbols(389)
            sage: S = M.cuspidal_submodule()
            sage: S.cuspidal_submodule() is S
            True
        """
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            try:
                if self.__is_cuspidal:
                    return self
            except AttributeError:
                pass
            S = self.boundary_map().kernel()
            S._set_is_cuspidal(True)
            S._is_full_hecke_module = True
            ## We know the cuspidal subspace is stable, so
            ## if it's one-dimensional, it must be simple
            if S.dimension() == 1:
                S._is_simple = True
            if self.base_ring().characteristic() == 0:
                d = self._cuspidal_submodule_dimension_formula()
                if not d is None:
                    assert d == S.dimension(), "According to dimension formulas the cuspidal subspace of \"%s\" has dimension %s; however, computing it using modular symbols we obtained %s, so there is a bug (please report!)."%(self, d, S.dimension())
            self.__cuspidal_submodule = S
        return self.__cuspidal_submodule

    def _degeneracy_raising_matrix(self, level):
        raise NotImplementedError

    def _degeneracy_lowering_matrix(self, level, t):
        # Use Proposition 2.6.15 in Merel's 1585 paper (or Prop 15 in
        # electronic version of that paper).
        H = heilbronn.HeilbronnMerel(t)
        M = self.hecke_module_of_level(level)
        return self.__heilbronn_operator(M,H,t).matrix()

    def rank(self):
        """
        Returns the rank of self.

        INPUT:


        -  ``ModularSymbols self`` - arbitrary space of modular
           symbols


        OUTPUT:


        -  ``int`` - the rank


        EXAMPLES::

            sage: M = ModularSymbols(389)
            sage: M.rank()
            65

        ::

            sage: ModularSymbols(11,sign=0).rank()
            3
            sage: ModularSymbols(100,sign=0).rank()
            31
            sage: ModularSymbols(22,sign=1).rank()
            5
            sage: ModularSymbols(1,12).rank()
            3
            sage: ModularSymbols(3,4).rank()
            2
            sage: ModularSymbols(8,6,sign=-1).rank()
            3
        """
        try:
            return self.__rank
        except AttributeError:
            self.__rank = len(self.manin_basis())
        return self.__rank

    def eisenstein_submodule(self):
        """
        Return the Eisenstein submodule of this space of modular symbols.
        """
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = self.cuspidal_submodule().complement()
            return self.__eisenstein_submodule

    def element(self, x):
        """
        Creates and returns an element of self from a modular symbol, if
        possible.

        INPUT:


        -  ``x`` - an object of one of the following types:
           ModularSymbol, ManinSymbol.


        OUTPUT:

        ModularSymbol - a modular symbol with parent self.
        """
        if isinstance(x, manin_symbols.ManinSymbol):
            if not x.parent().weight() == self.weight():
                raise ArithmeticError, "incompatible weights: Manin symbol has weight %s, but modular symbols space has weight %s"%(
                    x.parent().weight(), self.weight())
            t = self.manin_symbols().index(x.tuple())
            if isinstance(t, tuple):
                i, scalar = t
                v = self.manin_gens_to_basis().row(i) * scalar
            else:
                v = self.manin_gens_to_basis().row(t)
            return element.ModularSymbolsElement(self, v)

        elif isinstance(x, element.ModularSymbolsElement):
            M = x.parent()
            if M.ambient_hecke_module() != self:
                # TODO -- sometimes do something more sophisticated here.
                raise TypeError, "Modular symbol (%s) does not lie in this space."%x
            return self(x.element())

        else:
            raise ValueError, "Cannot create element of %s from %s."%(x,self)

    def dual_star_involution_matrix(self):
        """
        Return the matrix of the dual star involution, which is induced by
        complex conjugation on the linear dual of modular symbols.
        """
        try:
            return self.__dual_star_involution_matrix
        except AttributeError:
            pass
        self.__dual_star_involution_matrix = self.star_involution().matrix().transpose()
        return self.__dual_star_involution_matrix

    def factorization(self):
        r"""
        Returns a list of pairs `(S,e)` where `S` is spaces
        of modular symbols and self is isomorphic to the direct sum of the
        `S^e` as a module over the *anemic* Hecke algebra adjoin
        the star involution. The cuspidal `S` are all simple, but
        the Eisenstein factors need not be simple.

        EXAMPLES::

            sage: ModularSymbols(Gamma0(22), 2).factorization()
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field)^2 *
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field)^2 *
            (Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field)

        ::

            sage: ModularSymbols(1,6,0,GF(2)).factorization()
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 6 with sign 0 over Finite Field of size 2) *
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 6 with sign 0 over Finite Field of size 2)

        ::

            sage: ModularSymbols(18,2).factorization()
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 for Gamma_0(18) of weight 2 with sign 0 over Rational Field) *
            (Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 7 for Gamma_0(18) of weight 2 with sign 0 over Rational Field)

        ::

            sage: M = ModularSymbols(DirichletGroup(38,CyclotomicField(3)).1^2,  2, +1); M
            Modular Symbols space of dimension 7 and level 38, weight 2, character [1, zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2
            sage: M.factorization()                    # long time (about 8 seconds)
            (Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 7 and level 38, weight 2, character [1, zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2) *
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 and level 38, weight 2, character [1, zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2) *
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 and level 38, weight 2, character [1, zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2) *
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 and level 38, weight 2, character [1, zeta3], sign 1, over Cyclotomic Field of order 3 and degree 2)
        """

##         EXAMPLES:
##             sage: M = ModularSymbols(Gamma0(22), 2); M
##             Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field
##             sage: M.factorization():
##             ...    print b.dimension(), b.level(), e
##             1 11 2
##             1 11 2
##             1 11 2
##             1 22 1

##         An example with sign 1:
##             sage: M = ModularSymbols(Gamma0(22), 2, sign=1); M
##             Modular Symbols space of dimension 5 for Gamma_0(22) of weight 2 with sign 1 over Rational Field
##             sage: for b, e in M.factorization():
##             ...    print b.dimension(), b.level(), e
##             1 11 2
##             1 11 2
##             1 22 1

##         An example for Gamma1:
##             sage: M = ModularSymbols(Gamma1(26), 2, sign=1); M
##             Modular Symbols space of dimension 33 for Gamma_1(26) of weight 2 with sign 1 and over Rational Field
##             sage: for b, e in M.factorization():
##             ...    print b.dimension(), b.level(), e
##             1 13 2
##             1 13 2
##             1 13 2
##             2 13 2
##             2 13 2
##             2 13 2
##             2 13 2
##             2 13 2
##             1 26 1
##             1 26 1
##             1 26 1
##             2 26 1
##             2 26 1

##         An example with level divisible by a square:
##             sage: M = ModularSymbols(Gamma0(2*9),2); M
##             ???
##             sage: for b, e in M.factorization():
##             ...    print b.dimension(), b.level(), e
##             ???
        try:
            return self._factorization
        except AttributeError:
            pass

        try:
            if self._is_simple:
                return [self]
        except AttributeError:
            pass

        D = []

        # Treat the cuspidal and eisenstein parts separately.  The
        # cuspidal part is very straightforward because of
        # Atkin-Lehner-Li theory.  The eisenstein part is trickier,
        # because of E2 and that the new and old Eisenstein subspaces
        # can intersect (e.g., they do for M_2(Gamma_0(6))), even
        # in a way that involves forms other than E_2 (i.e., twists
        # of E2).

        # 1. Cuspidal part -- compute the factors and their multiplicities
        #                     using Atkin-Lehner-Li.

        # 2. Eisenstein part -- just call normal decomposition.

        # In the special case of weight 2 we have to do a bunch of
        # annoying extra work below to deal with the Eisenstein series E_2.
        k = self.weight()

        ## If the characteristic of the base ring is 2,
        ## the star involution is the identity, so we
        ## want to avoid adding each cuspidal submodule
        ## twice.
        if self.base_ring().characteristic() == 2:
            skip_minus = True
        else:
            skip_minus = False

        # The cuspidal part
        # We only run through spaces of level a multiple of the conductor of the character, which
        # we compute below, or set to 1 in case of Gamma_H or Gamma_1
        chi = self.character()
        cond = 1 if chi is None   else   chi.conductor()
        # Now actually run through the divisor levels, taking only the ones with that are
        # a multiple of the conductor.
        for d in reversed(arith.divisors(self.level())):
            if d%cond != 0: continue
            n = arith.number_of_divisors(self.level() // d)
            M = self.modular_symbols_of_level(d)
            N = M.new_submodule().cuspidal_submodule().decomposition()
            for A in N:
                if self.sign() == 0:
                    V = A.plus_submodule()
                    V._is_simple = True
                    D.append((V,n))
                    if skip_minus:
                        continue
                    V = A.minus_submodule()
                    V._is_simple = True
                    D.append((V,n))
                else:
                    A._is_simple = True
                    D.append((A,n))
        # The eisenstein part
        for E in self.eisenstein_submodule().decomposition(anemic=True):
            D.append((E,1))

        r = self.dimension()
        s = sum([A.rank()*mult for A, mult in D])
        D = sage.structure.all.Factorization(D, cr=True, sort=False)
        D.sort(_cmp = cmp)
        assert r == s, "bug in factorization --  self has dimension %s, but sum of dimensions of factors is %s\n%s"%(r, s, D)
        self._factorization = D
        return self._factorization

    def factor(self):
        """
        Synonym for self.factorization().
        """
        return self.factorization()

    def is_cuspidal(self):
        try:
            return self.__is_cuspidal
        except AttributeError:
            S = self.ambient_hecke_module().cuspidal_submodule()
            self.__is_cuspidal = (S.dimension() == self.dimension())
        return self.__is_cuspidal

    def is_eisenstein(self):
        try:
            return self.__is_eisenstein
        except AttributeError:
            S = self.ambient_hecke_module().eisenstein_submodule()
            self.__is_eisenstein = self.is_subspace(S)
        return self.__is_eisenstein

    def manin_symbols_basis(self):
        """
        A list of Manin symbols that form a basis for the ambient space
        self. INPUT:


        -  ``ModularSymbols self`` - an ambient space of
           modular symbols


        OUTPUT:


        -  ``list`` - a list of 2-tuples (if the weight is 2)
           or 3-tuples, which represent the Manin symbols basis for self.


        EXAMPLES::

            sage: m = ModularSymbols(23)
            sage: m.manin_symbols_basis()
            [(1,0), (1,17), (1,19), (1,20), (1,21)]
            sage: m = ModularSymbols(6, weight=4, sign=-1)
            sage: m.manin_symbols_basis()
            [[X^2,(2,1)]]
        """
        s = self.manin_symbols()
        return [s.manin_symbol(i) for i in self.manin_basis()]


    def modular_symbols_of_sign(self, sign):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, level, etc.) as this space except with given
        sign.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(11),2,sign=0)
            sage: M
            Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: M.modular_symbols_of_sign(-1)
            Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field
            sage: M = ModularSymbols(Gamma1(11),2,sign=0)
            sage: M.modular_symbols_of_sign(-1)
            Modular Symbols space of dimension 1 for Gamma_1(11) of weight 2 with sign -1 and over Rational Field
        """
        if sign == self.sign():
            return self
        return modsym.ModularSymbols(self.group(), self.weight(), sign=sign, base_ring=self.base_ring())


    def modular_symbols_of_weight(self, k):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with weight
        k.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(6),2,sign=0)
            sage: M.modular_symbols_of_weight(3)
            Modular Symbols space of dimension 4 for Gamma_1(6) of weight 3 with sign 0 and over Rational Field
        """
        if k == self.weight():
            return self
        return modsym.ModularSymbols(self.group(), weight=k, sign=self.sign(), base_ring=self.base_ring())

    def _compute_sign_submodule(self, sign, compute_dual=True):
        """
        Return the subspace of self that is fixed under the star
        involution.

        INPUT:


        -  ``sign`` - int (either -1 or +1)

        -  ``compute_dual`` - bool (default: True) also
           compute dual subspace. This are useful for many algorithms.


        OUTPUT: subspace of modular symbols

        EXAMPLES::

            sage: ModularSymbols(1,12,0,GF(5)).minus_submodule() ## indirect doctest
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Finite Field of size 5
        """
        S = self.star_involution().matrix() - self.base_ring()(sign)
        V = S.kernel()
        if compute_dual:
            Vdual = S.transpose().kernel()
            M = self.submodule(V, Vdual, check=False)
        else:
            M = self.submodule(V, check=False)
        M._set_sign(sign)
        return M

    def star_involution(self):
        """
        Return the star involution on self, which is induced by complex
        conjugation on modular symbols.
        """
        try:
            return self.__star_involution
        except AttributeError:
            pass
        S = self.__heilbronn_operator(self, [[-1,0, 0,1]], 1)
        S.name("Star involution on %s"%self)
        self.__star_involution = S
        return self.__star_involution

    def submodule(self, M, dual_free_module=None, check=True):
        """
        Return the submdoule of self with given generators or free module
        M.

        INPUT:


        -  ``M`` - a submodule of the ambient free module or
           generators for a submodule

        -  ``dual_free_module (optional)`` - this may be
           useful to speed up certain calculations; it is the corresponding
           submodule of the ambient dual module

        -  ``check`` - bool (optional: default - True); if
           True, actually check that M is a module, which means it is
           invariant under all Hecke operators.


        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: M.submodule([M.0])
            Traceback (most recent call last):
            ...
            ValueError: The submodule must be invariant under all Hecke operators.
            sage: M.eisenstein_submodule().basis()
            ((1,0) - 1/5*(1,9),)
            sage: M.basis()
            ((1,0), (1,8), (1,9))
            sage: M.submodule([M.0 - 1/5*M.2])
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field

        .. note::

           It would make more sense to only check that M is invariant
           under the Hecke operators with index coprime to the level.
           Unfortunately, I do not know a reasonable algorithm for
           determining whether a module is invariant under just the
           anemic Hecke algebra, since I do not know an analogue of
           the Sturm bound for the anemic Hecke algebra. - William
           Stein, 2007-07-27
        """
        if check:
            if not free_module.is_FreeModule(M):
                V = self.free_module()
                if isinstance(M, (list,tuple)):
                    M = V.span([V(x.element()) for x in M])
                else:
                    M = V.span(M)
        return subspace.ModularSymbolsSubspace(self, M, dual_free_module=dual_free_module, check=check)

    def twisted_winding_element(self, i, eps):
        r"""
        Given a space of modular symbols, an integer 0 = i = k-2, and a
        Dirichlet character eps, return the so-called 'twisted winding
        element':

        .. math::

                       \sum_{a \in (\mathbb{Z}/m\mathbb{Z})^\times} \varepsilon(a) * [ i, 0, a/m ].
        """

        if not dirichlet.is_DirichletCharacter(eps):
            raise TypeError, "eps must be a Dirichlet character."
        if (i < 0) or (i > self.weight()-2):
            raise ValueError, "i must be between 0 and k-2."

        m = eps.modulus()
        s = self(0)

        for a in ([ x for x in range(1,m) if rings.gcd(x,m) == 1 ]):
            s += eps(a) * self.modular_symbol([i, Cusp(0), Cusp(a/m)])

        return s

    ######################################################################
    # Z-module of integral modular symbols.
    #######################################################################
    def integral_structure(self, algorithm='default'):
        r"""
        Return the `\mathbb{Z}`-structure of this modular symbols
        spaces generated by all integral modular symbols.

        INPUT:


        -  ``algorithm`` - string (default: 'default' - choose
           heuristically)

        -  ``'pari'`` - use pari for the HNF computation

        -  ``'padic'`` - use p-adic algorithm (only good for
           dense case)


        ALGORITHM: It suffices to consider lattice generated by the free
        generating symbols `X^iY^{k-2-i}.(u,v)` after quotienting
        out by the `S` (and `I`) relations, since the
        quotient by these relations is the same over any ring.

        EXAMPLES: In weight 2 the rational basis is often integral.

        ::

            sage: M = ModularSymbols(11,2)
            sage: M.integral_structure()
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]

        This is rarely the case in higher weight::

            sage: M = ModularSymbols(6,4)
            sage: M.integral_structure()
            Free module of degree 6 and rank 6 over Integer Ring
            Echelon basis matrix:
            [  1   0   0   0   0   0]
            [  0   1   0   0   0   0]
            [  0   0 1/2 1/2 1/2 1/2]
            [  0   0   0   1   0   0]
            [  0   0   0   0   1   0]
            [  0   0   0   0   0   1]

        Here is an example involving `\Gamma_1(N)`.

        ::

            sage: M = ModularSymbols(Gamma1(5),6)
            sage: M.integral_structure()
            Free module of degree 10 and rank 10 over Integer Ring
            Echelon basis matrix:
            [     1      0      0      0      0      0      0      0      0      0]
            [     0      1      0      0      0      0      0      0      0      0]
            [     0      0  1/102      0  5/204  1/136  23/24   3/17 43/136 69/136]
            [     0      0      0   1/48      0   1/48  23/24    1/6    1/8  17/24]
            [     0      0      0      0   1/24      0  23/24    1/3    1/6    1/2]
            [     0      0      0      0      0   1/24  23/24    1/3  11/24   5/24]
            [     0      0      0      0      0      0      1      0      0      0]
            [     0      0      0      0      0      0      0    1/2      0    1/2]
            [     0      0      0      0      0      0      0      0    1/2    1/2]
            [     0      0      0      0      0      0      0      0      0      1]
        """
        if not self.base_ring() == rational_field.RationalField():
            raise NotImplementedError

        try:
            return self.__integral_structure
        except AttributeError:
            pass

        # The attribute _mod2term is set by self.compute_presentation().
        # It is a list of pairs (n, c), such that the ith element of the list
        # is equivalent to c times the n-th basis Manin symbol.
        G = set([i for i, _ in self._mod2term])

        # Now G is a set of integer i such that these integers gives
        # indices of Manin symbols that together generate the integral
        # structure.  We next obtain the corresponding list of elements
        # by passing to the quotient by the remaining relations
        # via the _manin_gens_to_basis attribute.

        # Next we take each element of X, which gives a linear combination
        # of the basis of the underlying vector space of self, and compute
        # the Z-module they span.

        G = list(G)
        G.sort()
        B = self._manin_gens_to_basis.matrix_from_rows(list(G)).dense_matrix()
        B, d = B._clear_denom()
        if algorithm == 'default':
            # pari is much better in the weight 2 case when the input
            # matrix is extremely sparse; the p-adic algorithm is
            # terrible in the sparse case.
            if self.weight() == 2:
                algorithm = 'pari'
            else:
                algorithm = 'padic'
        if algorithm == 'pari':
            B = B.echelon_form(algorithm='pari', include_zero_rows=False)
        elif algorithm == 'padic':
            B = B.echelon_form(algorithm='padic', include_zero_rows=False)
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm
        W = B.row_module()
        if d != 1:
            W = W.scale(1/d)
        self.__integral_structure = W
        assert W.rank() == self.rank(), "there is a bug in computing integral structure on modular symbols"
        return self.__integral_structure

    ######################################################################
    # Eigenvalues
    #######################################################################
    def compact_newform_eigenvalues(self, v, names='alpha'):
        r"""
        Return compact systems of eigenvalues for each Galois conjugacy
        class of cuspidal newforms in this ambient space.

        INPUT:


        -  ``v`` - list of positive integers


        OUTPUT:


        -  ``list`` - of pairs (E, x), where E\*x is a vector
           with entries the eigenvalues `a_n` for
           `n \in v`.


        EXAMPLES::

            sage: M = ModularSymbols(43,2,1)
            sage: X = M.compact_newform_eigenvalues(prime_range(10))
            sage: X[0][0] * X[0][1]
            (-2, -2, -4, 0)
            sage: X[1][0] * X[1][1]
            (alpha1, -alpha1, -alpha1 + 2, alpha1 - 2)

        ::

            sage: M = ModularSymbols(DirichletGroup(24,QQ).1,2,sign=1)
            sage: M.compact_newform_eigenvalues(prime_range(10),'a')
            [([-1/2 -1/2]
            [ 1/2 -1/2]
            [  -1    1]
            [  -2    0], (1, -2*a0 - 1))]
            sage: a = M.compact_newform_eigenvalues([1..10],'a')[0]
            sage: a[0]*a[1]
            (1, a0, a0 + 1, -2*a0 - 2, -2*a0 - 2, -a0 - 2, -2, 2*a0 + 4, -1, 2*a0 + 4)
            sage: M = ModularSymbols(DirichletGroup(13).0^2,2,sign=1)
            sage: M.compact_newform_eigenvalues(prime_range(10),'a')
            [([  -zeta6 - 1]
            [ 2*zeta6 - 2]
            [-2*zeta6 + 1]
            [           0], (1))]
            sage: a = M.compact_newform_eigenvalues([1..10],'a')[0]
            sage: a[0]*a[1]
            (1, -zeta6 - 1, 2*zeta6 - 2, zeta6, -2*zeta6 + 1, -2*zeta6 + 4, 0, 2*zeta6 - 1, -zeta6, 3*zeta6 - 3)
        """
        if self.sign() == 0:
            raise ValueError, "sign must be nonzero"
        v = list(v)

        # Get decomposition of this space
        D = self.cuspidal_submodule().new_subspace().decomposition()
        for A in D:
            # since sign is zero and we're on the new cuspidal subspace
            # each factor is definitely simple.
            A._is_simple = True
        B = [A.dual_free_module().basis_matrix().transpose() for A in D]

        # Normalize the names strings.
        names = ['%s%s'%(names,i) for i in range(len(B))]

        # Find an integer i such that the i-th columns of the basis for the
        # dual modules corresponding to the factors in D are all nonzero.
        nz = None
        for i in range(self.dimension()):
            # Decide if this i works, i.e., ith row of every element of B is nonzero.
            bad = False
            for C in B:
                if C.row(i) == 0:
                    # i is bad.
                    bad = True
                    continue
            if bad: continue
            # It turns out that i is not bad.
            nz = i
            break

        if nz is not None:
            R = self.hecke_images(nz, v)
            return [(R*m, D[i].dual_eigenvector(names=names[i], lift=False, nz=nz)) for i, m in enumerate(B)]
        else:
            # No single i works, so we do something less uniform.
            ans = []
            cache = {}
            for i in range(len(D)):
                nz = D[i]._eigen_nonzero()
                if cache.has_key(nz):
                     R = cache[nz]
                else:
                     R = self.hecke_images(nz, v)
                     cache[nz] = R
                ans.append((R*B[i], D[i].dual_eigenvector(names=names[i], lift=False, nz=nz)))
            return ans


class ModularSymbolsAmbient_wtk_g0(ModularSymbolsAmbient):
    r"""
    Modular symbols for `\Gamma_0(N)` of integer weight
    `k > 2` over the field `F`.
    """
    def __init__(self, N, k, sign, F):
        r"""
        Initialize a space of modular symbols of weight `k` for
        `\Gamma_0(N)`, over `\mathbb{Q}`.

        For weight `2`, it is faster to use
        ``ModularSymbols_wt2_g0``.

        INPUT:


        -  ``N`` - int, the level

        -  ``k`` - integer weight = 2.

        -  ``sign`` - int, either -1, 0, or 1

        -  ``F`` - field


        EXAMPLES::

            sage: ModularSymbols(1,12)
            Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
            sage: ModularSymbols(1,12, sign=1).dimension()
            2
            sage: ModularSymbols(15,4, sign=-1).dimension()
            4
            sage: ModularSymbols(6,6).dimension()
            10
            sage: ModularSymbols(36,4).dimension()
            36
        """
        N = int(N)
        k = int(k)
        sign = int(sign)
        if not sign in [-1,0,1]:
            raise TypeError, "sign must be an int in [-1,0,1]"

        ModularSymbolsAmbient.__init__(self, weight=k, group=congroup.Gamma0(N),
                                       sign=sign, base_ring=F)


    def _dimension_formula(self):
        if self.base_ring().characteristic() == 0:
            N, k, sign = self.level(), self.weight(), self.sign()
            if sign != 0: return None
            if k%2 == 1:
                return 0
            elif k > 2:
                return 2*dims.dimension_cusp_forms_gamma0(N,k) + dims.c0(N)
            else:
                return 2*dims.dimension_cusp_forms_gamma0(N,k) + dims.c0(N)-1
        else:
            raise NotImplementedError

    def _repr_(self):
        return ("Modular Symbols space of dimension %s for Gamma_0(%s) of weight %s with sign %s " + \
                "over %s")%(self.dimension(), self.level(),self.weight(), self.sign(),
                            self.base_ring())

    def _cuspidal_submodule_dimension_formula(self):
        if self.base_ring().characteristic() == 0:
            N, k, sign = self.level(), self.weight(), self.sign()
            if sign == 0:
                m = 2
            else:
                m = 1
            return m * dims.dimension_cusp_forms_gamma0(N, k)
        else:
            raise NotImplementedError


    def _degeneracy_raising_matrix(self, level):
        level = int(level)
        N = self.level()
        M = self.hecke_module_of_level(level)

        # 1. Find coset representatives H for Gamma_0(M.level()) \ Gamma_0(self.level())
        H = congroup.degeneracy_coset_representatives_gamma0(level, N, 1)
        # 2. The map is
        #        [P,pi(g)] |--> sum_{h in H} [P, pi(h*g)]
        #
        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        if self.dimension() == 0 or M.dimension() == 0:
            return MS(0)
        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        k = self.weight()
        G = matrix_space.MatrixSpace(integer_ring.IntegerRing(),2)
        H = [G(h) for h in H]
        for n in B:
            z = M(0)
            s = syms.manin_symbol(n)
            g = G(list(s.lift_to_sl2z(N)))
            i = s.i
            # We apply each matrix in H according to the above formula
            for h in H:
                hg = h*g
                z += M((i, hg[1,0], hg[1,1]))
            rows.append(z.element())

        A = MS(rows)
        return A


    def _cuspidal_new_submodule_dimension_formula(self):
        if self.base_ring().characteristic() == 0:
            N, k, sign = self.level(), self.weight(), self.sign()
            if sign == 0:
                m = 2
            else:
                m = 1
            return m * dims.dimension_new_cusp_forms_gamma0(N, k)
        else:
            raise NotImplementedError

    def boundary_space(self):
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_g0(
            self.level(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def manin_symbols(self):
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = manin_symbols.ManinSymbolList_gamma0(
                level=self.level(), weight=self.weight())
        return self.__manin_symbols

    def modular_symbols_of_level(self, N):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with the level
        N.

        For example, if self is the space of modular symbols of weight 2
        for Gamma_0(22), and level is 11, then this function returns
        modular symbols of weight 2 for Gamma_0(11).

        EXAMPLES::

            sage: M = ModularSymbols(11)
            sage: M.modular_symbols_of_level(22)
            Modular Symbols space of dimension 7 for Gamma_0(22) of weight 2 with sign 0 over Rational Field
            sage: M = ModularSymbols(Gamma1(6))
            sage: M.modular_symbols_of_level(12)
            Modular Symbols space of dimension 9 for Gamma_1(12) of weight 2 with sign 0 and over Rational Field
        """
        return modsym.ModularSymbols(congroup.Gamma0(rings.Integer(N)),
                                     self.weight(), sign=self.sign(),
                                     base_ring=self.base_ring())


    def _hecke_images(self, i, v):
        """
        Return matrix whose rows are the images of the `i`-th
        standard basis vector under the Hecke operators `T_p` for
        all integers in `v`.

        INPUT:


        -  ``i`` - nonnegative integer

        -  ``v`` - a list of positive integer


        OUTPUT:


        -  ``matrix`` - whose rows are the Hecke images


        EXAMPLES::

            sage: M = ModularSymbols(11,4,1)
            sage: M._hecke_images(0,[1,2,3,4])
            [ 1  0  0  0]
            [ 9  0  1 -1]
            [28  2 -1 -1]
            [73  2  5 -7]
            sage: M.T(1)(M.0).element()
            (1, 0, 0, 0)
            sage: M.T(2)(M.0).element()
            (9, 0, 1, -1)
            sage: M.T(3)(M.0).element()
            (28, 2, -1, -1)
            sage: M.T(4)(M.0).element()
            (73, 2, 5, -7)
            sage: M = ModularSymbols(12,4)
            sage: M._hecke_images(0,[1,2,3,4])
            [  1   0   0   0   0   0   0   0   0   0   0   0]
            [  8   1  -1  -2   2   2  -3   1  -2   3  -1   0]
            [ 27   4  -4  -8   8  10 -14   4  -9  14  -5   0]
            [ 64  10 -10 -20  20  26 -36  10 -24  38 -14   0]
            sage: M.T(1)(M.0).element()
            (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: M.T(2)(M.0).element()
            (8, 1, -1, -2, 2, 2, -3, 1, -2, 3, -1, 0)
            sage: M.T(3)(M.0).element()
            (27, 4, -4, -8, 8, 10, -14, 4, -9, 14, -5, 0)
            sage: M.T(4)(M.0).element()
            (64, 10, -10, -20, 20, 26, -36, 10, -24, 38, -14, 0)
        """
        # Find basis vector for ambient space such that it is not in
        # the kernel of the dual space corresponding to self.
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        return heilbronn.hecke_images_gamma0_weight_k(c.u,c.v, c.i, N, self.weight(),
                                                      v, self.manin_gens_to_basis())

class ModularSymbolsAmbient_wt2_g0(ModularSymbolsAmbient_wtk_g0):
    """
    Modular symbols for Gamma_0(N) of integer weight 2 over the field
    F.
    """
    def __init__(self, N, sign, F):
        """
        Initialize a space of modular symbols. INPUT:


        -  ``N`` - int, the level

        -  ``sign`` - int, either -1, 0, or 1


        OUTPUT: The space of modular symbols of weight 2, trivial
        character, level N and given sign.

        EXAMPLES::

            sage: M = ModularSymbols(Gamma0(12),2)
        """
        ModularSymbolsAmbient_wtk_g0.__init__(self,
                                                N=N, k=2, sign=sign, F=F)

    def _dimension_formula(self):
        if self.base_ring().characteristic() == 0:
            N, sign = self.level(), self.sign()
            if sign != 0: return None
            return 2*dims.dimension_cusp_forms_gamma0(N,2) + dims.c0(N) - 1
        else:
            raise NotImplementedError

    def _cuspidal_submodule_dimension_formula(self):
        if self.base_ring().characteristic() == 0:
            if self.sign() == 0:
                m = 2
            else:
                m = 1
            return m * dims.dimension_cusp_forms_gamma0(self.level(), 2)
        else:
            raise NotImplementedError

    def _cuspidal_new_submodule_dimension_formula(self):
        if self.base_ring().characteristic() == 0:
            if self.sign() == 0:
                m = 2
            else:
                m = 1
            return m * dims.dimension_new_cusp_forms_gamma0(self.level(), 2)
        else:
            raise NotImplementedError


    def _compute_hecke_matrix_prime(self, p, rows=None):
        """
        Compute and return the matrix of the p-th Hecke operator.
        EXAMPLES::

            sage: m = ModularSymbols(37,2)
            sage: m._compute_hecke_matrix_prime(2).charpoly('x')
            x^5 + x^4 - 8*x^3 - 12*x^2
        """
        # note -- p doesn't have to be prime.
        if isinstance(rows, list):
            rows = tuple(rows)
        try:
            return self._hecke_matrices[(p,rows)]
        except AttributeError:
            self._hecke_matrices = {}
        except KeyError:
            pass
        tm = misc.verbose("Computing Hecke operator T_%s"%p)

        H = heilbronn.HeilbronnCremona(p)
        ##H = heilbronn.HeilbronnMerel(p)
        B = self.manin_basis()
        if not rows is None:
            B = [B[i] for i in rows]

        cols = []
        N = self.level()
        P1 = self.p1list()
        mod2term = self._mod2term
        R = self.manin_gens_to_basis()
        W = R.new_matrix(nrows=len(B), ncols = R.nrows())  # the 0 with given number of rows and cols.
        j = 0
        tm = misc.verbose("Matrix non-reduced", tm)
        for i in B:
            # The following step is where most of the time is spent.
            c,d = P1[i]
            v = H.apply(c,d, N)

            # v is now a list of pairs ((c,d),m), where m is the
            # number of times that (c,d) appears in the image of x
            # under the matrices in H.  Also, the pairs (c,d) are
            # normalized.
            # Let ind(c,d) denote the index of the normalized pair
            # (c,d) in the fixed ordered list of elements of
            # P1(Z/NZ).  Then the list of pairs (ind(c,d), m)
            # obtained from the above list defines a sparse vector
            # s, and the image of x under T_p is the product
            # of s with the matrix R defined above.
            for z, m in v:
                k = P1.index_of_normalized_pair(z[0],z[1])
                if k != -1:
                    f, s = mod2term[k]
                    if s != 0:
                        W[j,f] = W[j,f] + s*m
            j += 1
        tm = misc.verbose("done making non-reduced matrix",tm)
        misc.verbose("start matrix-matrix (%s x %s) times (%s x %s) multiply to get Tp"%(W.nrows(), W.ncols(),
                                                                                         R.nrows(), R.ncols()))
        if hasattr(W, '_matrix_times_matrix_dense'):
            Tp = W._matrix_times_matrix_dense(R)
            misc.verbose("done matrix multiply and computing Hecke operator",tm)
        else:
            Tp = W * R
            tm = misc.verbose("done multiplying",tm)
            Tp = Tp.dense_matrix()
            misc.verbose("done making hecke operator dense",tm)
        if rows is None:
            self._hecke_matrices[(p,rows)] = Tp
        return Tp

    def boundary_space(self):
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_g0(
            self.level(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def _hecke_image_of_ith_basis_vector(self, n, i):
        """
        Return `T_n(e_i)`, where `e_i` is the
        `i`th basis vector of this ambient space.

        INPUT:


        -  ``n`` - an integer which should be prime.


        OUTPUT:


        -  ``modular symbol`` - element of this ambient space


        EXAMPLES::

            sage: M = ModularSymbols(43,2,1)
            sage: M._hecke_image_of_ith_basis_vector(2, 0)
            3*(1,0) - 2*(1,33)
            sage: M.hecke_operator(2)(M.0)
            3*(1,0) - 2*(1,33)
            sage: M._hecke_image_of_ith_basis_vector(6, 1)
            -2*(1,33)
            sage: M.hecke_operator(6)(M.1)
            -2*(1,33)
        """
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        I = heilbronn.hecke_images_gamma0_weight2(c.u,c.v,N,[n], self.manin_gens_to_basis())
        return self(I[0])

    def _hecke_images(self, i, v):
        """
        Return images of the `i`-th standard basis vector under the
        Hecke operators `T_p` for all integers in `v`.

        INPUT:


        -  ``i`` - nonnegative integer

        -  ``v`` - a list of positive integer


        OUTPUT:


        -  ``matrix`` - whose rows are the Hecke images


        EXAMPLES::

            sage: M = ModularSymbols(46,2,-1)
            sage: M._hecke_images(1,[3,4,5,6])
            [ 0  1 -2  2  0]
            [ 2 -3  0  0  1]
            [ 0 -2  2 -2  0]
            [-5  3 -1  1  0]
            sage: v = M.basis()[1]
            sage: M.T(3)(v).element()
            (0, 1, -2, 2, 0)
            sage: M.T(4)(v).element()
            (2, -3, 0, 0, 1)
            sage: M.T(5)(v).element()
            (0, -2, 2, -2, 0)
            sage: M.T(6)(v).element()
            (-5, 3, -1, 1, 0)
        """
        # Find basis vector for ambient space such that it is not in
        # the kernel of the dual space corresponding to self.
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        return heilbronn.hecke_images_gamma0_weight2(c.u,c.v,N, v, self.manin_gens_to_basis())


class ModularSymbolsAmbient_wtk_g1(ModularSymbolsAmbient):
    def __init__(self, level, weight, sign, F):
        """
        Initialize a space of modular symbols for Gamma1(N).

        INPUT:


        -  ``level`` - int, the level

        -  ``weight`` - int, the weight = 2

        -  ``sign`` - int, either -1, 0, or 1

        -  ``F`` - field


        EXAMPLES::

            sage: ModularSymbols(Gamma1(17),2)
            Modular Symbols space of dimension 25 for Gamma_1(17) of weight 2 with sign 0 and over Rational Field
            sage: [ModularSymbols(Gamma1(7),k).dimension() for k in [2,3,4,5]]
            [5, 8, 12, 16]

        ::

            sage: M = ModularSymbols(Gamma1(7),3)
        """
        ModularSymbolsAmbient.__init__(self,
                weight=weight,
                group=congroup.Gamma1(level),
                sign=sign,
                base_ring=F)


    def _dimension_formula(self):
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError
        level, weight, sign = self.level(), self.weight(), self.sign()
        if sign != 0: return None
        d = 2*dims.dimension_cusp_forms_gamma1(level,weight) + dims.c1(level)
        if level == 1 and weight%2 == 1:
            return 0
        if weight == 2:
            return d - 1
        if weight % 2 == 0:
            return d

        # TODO: I don't know a formula for dim ModSym_k(Gamma_1(N)) for odd k!!!

        return None

    def _repr_(self):
        return ("Modular Symbols space of dimension %s for Gamma_1(%s) of weight %s with sign %s " + \
                "and over %s")%(self.dimension(), self.level(),self.weight(),
                                self.sign(), self.base_ring())

    def _cuspidal_submodule_dimension_formula(self):
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * dims.dimension_cusp_forms_gamma1(self.level(), self.weight())

    def _cuspidal_new_submodule_dimension_formula(self):
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * dims.dimension_new_cusp_forms_gamma1(self.level(), self.weight())


    def _compute_hecke_matrix_prime_power(self, n, p, r):
        return self._compute_hecke_matrix_prime(n)

##     def _xxx_degeneracy_raising_matrix(self, M):
##         R = congroup.degeneracy_coset_representatives_gamma1(M.level(), self.level(), 1)
##         return self._matrix_of_operator_on_modular_symbols(M, R)

    def _degeneracy_raising_matrix(self, level):
        level = int(level)
        N = self.level()
        M = self.hecke_module_of_level(level)

        # 1. Find coset representatives H for Gamma_0(M.level()) \ Gamma_0(self.level())
        H = congroup.degeneracy_coset_representatives_gamma1(M.level(), N, 1)
        # 2. The map is
        #        [P,pi(g)] |--> sum_{h in H} [P, pi(h*g)]
        #
        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        if self.dimension() == 0 or M.dimension() == 0:
            return MS(0)
        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        k = self.weight()
        G = matrix_space.MatrixSpace(integer_ring.IntegerRing(),2)
        H = [G(h) for h in H]
        for n in B:
            z = M(0)
            s = syms.manin_symbol(n)
            g = G(list(s.lift_to_sl2z(N)))
            i = s.i
            # We apply each matrix in H according to the above formula
            for h in H:
                hg = h*g
                z += M((i, hg[1,0], hg[1,1]))
            rows.append(z.element())

        A = MS(rows)
        return A

    def boundary_space(self):
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_g1(
            self.level(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def manin_symbols(self):
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = manin_symbols.ManinSymbolList_gamma1(
                level=self.level(), weight=self.weight())
        return self.__manin_symbols


    def modular_symbols_of_level(self, N):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with the level
        N.
        """
        return modsym.ModularSymbols(congroup.Gamma1(N), self.weight(),self.sign(), self.base_ring())

class ModularSymbolsAmbient_wtk_gamma_h(ModularSymbolsAmbient):
    def __init__(self, group, weight, sign, F):
        r"""
        Initialize a space of modular symbols for `\Gamma_H(N)`.

        INPUT:


        -  ``group`` - a congruence subgroup
           `\Gamma_H(N)`.

        -  ``weight`` - int, the weight = 2

        -  ``sign`` - int, either -1, 0, or 1

        -  ``F`` - field


        EXAMPLES::

            sage: ModularSymbols(GammaH(15,[4]),2)
            Modular Symbols space of dimension 9 for Congruence Subgroup Gamma_H(15) with H generated by [4] of weight 2 with sign 0 and over Rational Field
        """
        ModularSymbolsAmbient.__init__(self,
                weight=weight, group=group,
                sign=sign, base_ring=F)

    def _dimension_formula(self):
        return None

    def _repr_(self):
        return ("Modular Symbols space of dimension %s for %s of weight %s with sign %s " + \
                "and over %s")%(self.dimension(), self.group(),self.weight(),
                                self.sign(), self.base_ring())

    def _cuspidal_submodule_dimension_formula(self):
        return None

    def _cuspidal_new_submodule_dimension_formula(self):
        return None

    def _compute_hecke_matrix_prime_power(self, n, p, r):
        raise NotImplementedError

    def _degeneracy_raising_matrix(self, level):
        raise NotImplementedError

    def boundary_space(self):
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_gamma_h(
            self.group(), self.weight(), self.sign(), self.base_ring())
        return self.__boundary_space

    def manin_symbols(self):
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = manin_symbols.ManinSymbolList_gamma_h(
                group=self.group(), weight=self.weight())
        return self.__manin_symbols

    def modular_symbols_of_level(self, N):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with the level
        N.
        """
        raise NotImplementedError


class ModularSymbolsAmbient_wtk_eps(ModularSymbolsAmbient):
    def __init__(self, eps, weight, sign=0):
        """
        Space of modular symbols with given weight, character, and sign.

        INPUT:


        -  ``eps`` - dirichlet.DirichletCharacter, the
           "Nebentypus" character.

        -  ``weight`` - int, the weight = 2

        -  ``sign`` - int, either -1, 0, or 1


        EXAMPLES::

            sage: eps = DirichletGroup(4).gen(0)
            sage: eps.order()
            2
            sage: ModularSymbols(eps, 2)
            Modular Symbols space of dimension 0 and level 4, weight 2, character [-1], sign 0, over Rational Field
            sage: ModularSymbols(eps, 3)
            Modular Symbols space of dimension 2 and level 4, weight 3, character [-1], sign 0, over Rational Field

        We next create a space with character of order bigger than 2.

        ::

            sage: eps = DirichletGroup(5).gen(0)
            sage: eps     # has order 4
            [zeta4]
            sage: ModularSymbols(eps, 2).dimension()
            0
            sage: ModularSymbols(eps, 3).dimension()
            2

        Here is another example::

            sage: G, e = DirichletGroup(5).objgen()
            sage: M = ModularSymbols(e,3)
            sage: loads(M.dumps()) == M
            True
        """
        level = eps.modulus()
        ModularSymbolsAmbient.__init__(self,
                weight = weight,
                group = congroup.Gamma1(level),
                sign = sign,
                base_ring = eps.base_ring(),
                character = eps)

    def _repr_(self):
        return ("Modular Symbols space of dimension %s and level %s, weight %s, character %s, sign %s, " + \
                "over %s")%(self.dimension(), self.level(), self.weight(),
                    self.character(), self.sign(), self.base_ring())


    def _cuspidal_submodule_dimension_formula(self):
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * dims.dimension_cusp_forms_eps(self.character(), self.weight())

    def _cuspidal_new_submodule_dimension_formula(self):
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError
        if self.sign() == 0:
            m = 2
        else:
            m = 1
        return m * dims.dimension_new_cusp_forms(self.character(), self.weight())

    def _matrix_of_operator_on_modular_symbols(self, codomain, R, character_twist=False):
        """
        INPUT:


        -  ``self`` - this space of modular symbols

        -  ``codomain`` - space of modular symbols

        -  ``R`` - list of lists [a,b,c,d] of length 4, which
           we view as elements of GL_2(Q).


        OUTPUT: a matrix, which represents the operator

        .. math::

                            x \mapsto \sum_{g in R} g.x


        where g.x is the formal linear fractional transformation on modular
        symbols.
        """
        eps = self.character()
        rows = []
        for b in self.basis():
            v = formal_sum.FormalSum(0, check=False)
            for c, x in b.modular_symbol_rep():
                for g in R:
                    y = x.apply(g)
                    if character_twist:
                        v += y*c*eps(g[0])
                    else:
                        v += y*c
            w = codomain(v).element()
            rows.append(w)
        M = matrix_space.MatrixSpace(self.base_ring(), len(rows), codomain.degree(), sparse=False)
        return M(rows)

##     def _xxx_degeneracy_raising_matrix(self, M):
##         R = congroup.degeneracy_coset_representatives_gamma0(M.level(), self.level(), 1)
##         return self._matrix_of_operator_on_modular_symbols(M, R, character_twist = True)

    def _degeneracy_raising_matrix(self, level):
        level = int(level)
        N = self.level()
        M = self.hecke_module_of_level(level)

        # 1. Find coset representatives H for Gamma_0(M.level()) \ Gamma_0(self.level())
        H = congroup.degeneracy_coset_representatives_gamma0(M.level(), N, 1)
        # 2. The map is
        #        [P,pi(g)] |--> sum_{h in H} [P, pi(h*g)]
        #
        MS = matrix_space.MatrixSpace(self.base_ring(), self.dimension(), M.dimension())
        if self.dimension() == 0 or M.dimension() == 0:
            return MS(0)
        rows = []
        B = self.manin_basis()
        syms = self.manin_symbols()
        k = self.weight()
        G = matrix_space.MatrixSpace(integer_ring.IntegerRing(),2)
        H = [G(h) for h in H]
        eps = self.character()  # note: in my thesis I twisted by eps^(-1), which is definitely a mistake
                                # since twisting by eps gives the right answer and by eps^(-1) does not.
        for n in B:
            z = M(0)
            s = syms.manin_symbol(n)
            g = G(list(s.lift_to_sl2z(N)))
            i = s.i
            # We apply each matrix in H according to the above formula
            for h in H:
                hg = h*g
                z += eps(h[0,0])*M((i, hg[1,0], hg[1,1]))
            rows.append(z.element())
        A = MS(rows)
        return A

    def _dimension_formula(self):
        return None

    def boundary_space(self):
        try:
            return self.__boundary_space
        except AttributeError:
            pass
        self.__boundary_space = boundary.BoundarySpace_wtk_eps(
            self.character(), self.weight(), self.sign())
        return self.__boundary_space

    def manin_symbols(self):
        try:
            return self.__manin_symbols
        except AttributeError:
            self.__manin_symbols = manin_symbols.ManinSymbolList_character(
                character=self.character(), weight=self.weight())
        return self.__manin_symbols

    def modular_symbols_of_level(self, N):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with the level
        N.

        For example, if self is the space of modular symbols of weight 2
        for Gamma_0(22), and level is 11, then this function returns
        modular symbols of weight 2 for Gamma_0(11).
        """
        if self.level() % N == 0:
            eps = self.character().restrict(N)
        elif N % self.level() == 0:
            eps = self.character().extend(N)
        else:
            raise ValueError, "The level N (=%s) must be a divisor or multiple of the modulus of the character (=%s)"%(N, self.level())
        return modsym.ModularSymbols(eps, self.weight(), self.sign(), self.base_ring())

    def modular_symbols_of_sign(self, sign):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with given
        sign.
        """
        return modsym.ModularSymbols(self.character(), self.weight(), sign, self.base_ring())

    def modular_symbols_of_weight(self, k):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, sign, etc.) as this space except with weight
        k.
        """
        return modsym.ModularSymbols(self.character(), k, self.sign(), self.base_ring())

    def _hecke_images(self, i, v):
        """
        Return images of the `i`-th standard basis vector under the
        Hecke operators `T_p` for all integers in `v`.

        INPUT:


        -  ``i`` - nonnegative integer

        -  ``v`` - a list of positive integer


        OUTPUT:


        -  ``matrix`` - whose rows are the Hecke images
        """
        if self.weight() != 2:
            raise NotImplementedError, "hecke images only implemented when the weight is 2"
        chi = self.character()
        # Find basis vector for ambient space such that it is not in
        # the kernel of the dual space corresponding to self.
        c = self.manin_generators()[self.manin_basis()[i]]
        N = self.level()
        if chi.order() > 2:
            return heilbronn.hecke_images_nonquad_character_weight2(c.u,c.v,N,
                                 v, chi, self.manin_gens_to_basis())
        else:
            return heilbronn.hecke_images_quad_character_weight2(c.u,c.v,N,
                                 v, chi, self.manin_gens_to_basis())
        raise NotImplementedError
