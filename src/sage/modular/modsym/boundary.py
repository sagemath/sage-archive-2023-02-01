"""
Space of boundary modular symbols.

Used mainly for computing the cuspidal subspace of modular symbols.
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
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
#*****************************************************************************

__doc_exclude = ['repr_lincomb', 'QQ']

# Python imports
import operator
import weakref

# SAGE imports
import sage.modules.free_module as free_module
from   sage.misc.misc import repr_lincomb

import sage.modular.congroup as congroup
import sage.modular.dirichlet as dirichlet
import sage.rings.all as rings
import sage.rings.arith as arith
import sage.rings.coerce as coerce
import sage.structure.gens as gens

import sage.modular.cusps as cusps

import ambient
import sage.modular.modsym.manin_symbols as manin_symbols

import sage.modular.hecke.all as hecke

import sage.modules.module_element as module_element

from sage.modules.all import is_FreeModuleElement

QQ = rings.RationalField()

import element

class BoundarySpaceElement(hecke.HeckeModuleElement):
    """
    A boundary symbol.
    """
    def __init__(self, parent, x):
        """
        Create a boundary symbol.

        INPUT:
            parent -- BoundarySpace; a space of boundary modular symbols
            x -- a dict with integer keys and values in the base
                 field of parent.
        """
        self.__parent = parent
        self.__x = x
        hecke.HeckeModuleElement.__init__(self, parent, self.element())

    def _repr_(self):
        """
        Returns a string representation for self for printing purposes.
        """
        g = self.__parent._known_gens_repr
        z = [0 for _ in xrange(len(g))]
        for i, c in self.__x.items():
            z[i] = c
        return repr_lincomb(g, z)

    def __add__(self, other):
        if not isinstance(other, BoundarySpaceElement):
            return coerce.bin_op(self, other, operator.add)
        z = dict(other.__x)
        for i, c in self.__x.items():
            if z.has_key(i):
                z[i] += c
            else:
                z[i] = c
        return BoundarySpaceElement(self.__parent, z)


    def __sub__(self, other):
        if not isinstance(other, BoundarySpaceElement):
            return coerce.bin_op(self, other, "-")
        z = dict(self.__x)
        for i, c in other.__x.items():
            if z.has_key(i):
                z[i] -= c
            else:
                z[i] = -c
        return BoundarySpaceElement(self.__parent, z)

    def __mul__(self, other):
        x = {}
        for i, c in self.__x.items():
            x[i] = c*other
        return BoundarySpaceElement(self.__parent, x)

    def __neg__(self):
        return self*(-1)

    def parent(self):
        return self.__parent

    def element(self):
        return self.__parent.free_module()(self.__x)


class BoundarySpace(hecke.HeckeModule_generic):
    """
    Space of boundary symbols for a congruence subgroup of SL_2(Z).

    This class is an abstract base class, so only derived classes should
    be instantiated.

    INPUT:
        weight -- int, the weight
        group -- congroup.CongruenceGroup, a congruence subgroup.
        sign -- int, either -1, 0, or 1
        base_ring -- rings.Ring (defaults to the rational numbers)
    """
    def __init__(self,
                 group = congroup.Gamma0(1),
                 weight = 2,
                 sign = 0,
                 base_ring = QQ,
                 character = None):
        """
        Initialize a space of boundary symbols.
        """
        weight = int(weight)
        if weight <= 1:
            raise ArithmeticError, "weight must be at least 2"
        if not isinstance(group, congroup.CongruenceSubgroup):
            raise TypeError, "group must be a congruence subgroup"
        sign = int(sign)
        if not isinstance(base_ring, rings.Ring) and rings.is_CommutativeRing(base_ring):
            raise TypeError, "base_ring must be a commutative ring"
        if character == None and isinstance(group, congroup.Gamma0):
            character = dirichlet.TrivialCharacter(group.level(), base_ring)
        (self.__group, self.__weight, self.__character,
          self.__sign, self.__base_ring) = (group, weight,
                                             character, sign, base_ring)
        self._known_gens = []
        self._known_gens_repr = []
        self._is_zero = []
        hecke.HeckeModule_generic.__init__(self, base_ring, group.level())

    def is_ambient(self):
        return True

    def group(self):
        """
        Return the congruence subgroup associated to this space of
        boundary modular symbols.
        """
        return self.__group

    def weight(self):
        """
        Return the weight of this space of boundary modular symbols.
        """
        return self.__weight

    def character(self):
        """
        Return the Dirichlet character assocaited to this space of
        boundary modular symbols.
        """
        return self.__character

    def sign(self):
        """
        Return the sign of the complex conjugation involution on this
        space of boundary modular symbols.
        """
        return self.__sign

    def gen(self, i=0):
        """
        Return the i-th generator of this space.
        """
        return BoundarySpaceElement(self, {i:1})

    def __len__(self):
        return len(self._known_gens)

    def free_module(self):
        return free_module.FreeModule(self.__base_ring, len(self._known_gens), sparse=True)

    def rank(self):
        """
        The rank of the space generated by boundary symbols that have
        been found so far in the course of computing the boundary map.

        WARNING: This number may change as more elements are coerced
        into this space!!  (This is an implementation detail that will
        likely change.)
        """
        return len(self._known_gens)

    #####################################################################
    # Coercion
    #####################################################################

    def _coerce_in_manin_symbol(self, x):
        i = x.i
        alpha, beta = x.endpoints(self.level())
        if self.weight() == 2:
            return self(alpha) - self(beta)
        if i == 0:
            return self(alpha)
        elif i == self.weight() - 2:
            return -self(beta)
        else:
            return self(0)

    def __call__(self, x):
        """
        Coerce x into a boundary symbol space.

        If x is a modular symbol (with the same group, weight,
        character, sign, and base field), this returns the image of
        that modular symbol under the boundary map.
        """
        if isinstance(x, int) and x == 0:
            return BoundarySpaceElement(self, {})

        elif isinstance(x, cusps.Cusp):
            return self._coerce_cusp(x)

        elif manin_symbols.is_ManinSymbol(x):
            return self._coerce_in_manin_symbol(x)

        elif element.is_ModularSymbolsElement(x):
            M = x.parent()
            if not isinstance(M, ambient.ModularSymbolsAmbient):
                raise TypeError, "x (=%s) must be an element of a space of modular symbols of type ModularSymbolsAmbient"%x
            if M.level() != self.level():
                raise TypeError, "x (=%s) must have level %s but has level %s"%(
                    x, self.level(), M.level())
            S = x.manin_symbol_rep()
            if len(S) == 0:
                return self(0)
            return sum([c*self._coerce_in_manin_symbol(v) for c, v in S])

        elif is_FreeModuleElement(x):
            y = dict([(i,x[i]) for i in xrange(len(x))])
            return BoundarySpaceElement(self, y)

        raise TypeError, "Coercion of %s (of type %s, parent %s) into %s not (yet) defined."%(x,type(x), x.parent(), self)

    def _repr_(self):
        return ("Space of Boundary Modular Symbols of weight %s for" + \
                " %s with sign %s and character %s over %s")%(
                 self.weight(), self.group(), self.sign(),
                 self.character(), self.base_ring())

    def _cusp_index(self, cusp):
        g = self._known_gens
        N = self.level()
        for i in xrange(len(g)):
            if self._is_equiv(cusp, g[i]):
                return i
        return -1

    def element(self, x):
        """
        Creates and returns an element of self from a modular or
        sage symbol, if possible.
        INPUT:
           x -- an object of one of the following types:
                ModularSymbol, ManinSymbol.  We ignore the group,
                weight, and character of the input modular or sage symbol.
        OUTPUT:
           ModularSymbol -- a modular symbol with parent self.
        """
        raise NotImplementedError

    def filename(self):
        """
        Returns the filename of self that should be used to store self
        in the database.
        INPUT:
           self -- space of modular symbols
        OUTPUT:
           str -- a string
        """
        return "boundsym-%s_%s_%s_%s_%s"%\
               (self.__group.name(),self.__weight,self.__character,\
                self.__sign,self.__base_ring.name())


class BoundarySpace_wtk_g0(BoundarySpace):
    """
    Boundary symbols for Gamma_0(N) of integer weight k > 2 over the field F.
    """
    def __init__(self, level, weight, sign, F):
        """
        Initialize a space of boundary symbols of weight k for
        Gamma_0(N), over F.

        For weight 2, it is faster to use BoundarySpace_wt2_g0.

        INPUT:
            level -- int, the level
            weight -- integer weight >= 2.
            sign -- int, either -1, 0, or 1
            F -- field
        """
        level = int(level)
        sign = int(sign)
        weight = int(weight)
        if not sign in [-1,0,1]:
            raise ArithmeticError, "sign must be an int in [-1,0,1]"
        if level <= 0:
            raise ArithmeticError, "level must be positive"
        BoundarySpace.__init__(self,
                                 weight = weight,
                                 group  = congroup.Gamma0(level),
                                 sign   = sign,
                                 base_ring = F)

    def _repr_(self):
        return ("Space of Boundary Modular Symbols for %s of weight %s " + \
                "and over %s")%(self.group(), self.weight(), self.base_ring())

    def _coerce_cusp(self, c):
        """
        Coerce cusp into a boundary symbol space.
        """
        if self.weight()%2 != 0:
            return self(0)
        N = self.level()
        sign = self.sign()
        i = self._cusp_index(c)
        if i != -1:
            if i in self._is_zero:
                return self(0)
            return BoundarySpaceElement(self, {i:1})

        if sign != 0:
            i2 = self._cusp_index(-c)
            if i2 != -1:
                if i2 in self._is_zero:
                    return self(0)
                return BoundarySpaceElement(self, {i2:sign})

        # found a new cusp class
        g = self._known_gens
        g.append(c)
        self._known_gens_repr.append("[%s]"%c)

        if sign == -1:
            # new cusp and nonzero sign, so if the sign is -1
            # its possible that the cusp class is killed by
            # the sign relations
            if self._is_equiv(c, -c):
                self._is_zero.append(len(g)-1)
                return self(0)

        return BoundarySpaceElement(self, {(len(g)-1):1})

    def _is_equiv(self, c1, c2):
        return c1.is_gamma0_equiv(c2, self.level())


class BoundarySpace_wtk_g1(BoundarySpace):
    def __init__(self, level, weight, sign, F):
        """
        Initialize a space of boundary modular symbols for Gamma1(N).

        INPUT:
            level -- int, the level
            weight -- int, the weight >= 2
            sign -- int, either -1, 0, or 1
            F -- base ring

        EXAMPLES:
            sage: from sage.modular.modsym.boundary import BoundarySpace_wtk_g1
            sage: BoundarySpace_wtk_g1(17, 2, 0, QQ)
            Boundary Modular Symbols space for Gamma_1(17) of weight 2 over Rational Field
        """
        level = int(level)
        sign = int(sign)
        if not sign in [-1,0,1]:
            raise ArithmeticError, "sign must be an int in [-1,0,1]"
        if level <= 0:
            raise ArithmeticError, "level must be positive"

        BoundarySpace.__init__(self,
                weight = weight,
                group  = congroup.Gamma1(level),
                sign   = sign,
                base_ring = F)

    def _repr_(self):
        return ("Boundary Modular Symbols space for Gamma_1(%s) of weight %s " + \
                "over %s")%(self.level(),self.weight(), self.base_ring())


    def _is_equiv(self, c1, c2):
        return c1.is_gamma1_equiv(c2, self.level())

    def _cusp_index(self, cusp):
        g = self._known_gens
        N = self.level()
        for i in xrange(len(g)):
            t, eps = self._is_equiv(cusp, g[i])
            if t:
                return i, eps
        return -1, 0

    def _coerce_cusp(self, c):
        """
        Coerce symbol into a boundary symbol space.
        """
        N    = self.level()
        k    = self.weight()
        sign = self.sign()
        i, eps = self._cusp_index(c)
        if i != -1:
            if i in self._is_zero:
                return self(0)
            return BoundarySpaceElement(self, {i : eps**k})

        if sign != 0:
            i2, eps = self._cusp_index(-c)
            if i2 != -1:
                if i2 in self._is_zero:
                    return self(0)
                return BoundarySpaceElement(self, {i2:sign*(eps**k)})

        # found a new cusp class
        g = self._known_gens
        g.append(c)
        self._known_gens_repr.append("[%s]"%c)

        ################################################################
        #
        # The set of boundary modular symbols for Gamma_1(N) is the
        # free abelian group on the set of pairs [P, [(u,v)]], where
        # the [(u,v)] are pairs with gcd(u,v) = 1 modulo the relations:
        #
        #        [(-u, -v)] = (-1)^k [(u,v)]
        #        [gamma(u,v)] = [(u,v)]  all gamma in Gamma_1(N).
        #
        # It's possible for the first two relations to kill a class,
        # i.e., for a pair [(u,v)] to be 0.  For example, when N=4,
        # u=1, v=2 and k=3 then (-1,-2) is equiv mod Gamma_1(4) to
        # (1,2) since v=-v (mod 4) and u=-u (mod 2).  But since k is
        # odd, [(-1,-2)] is also equivalent to -[(1,2)].  Thus this
        # symbol is equivalent to its negative, hence 0 (note: this
        # wouldn't be the case in char 2).  (See also prop 2.30 of
        # Stein Phd. thesis).
        #
        # When the sign is nonzero, we have the additional relations
        #
        #        [(-u,v)] = sign*[(u,v)]
        #
        ################################################################

        # Does cusp class vanish because of - relations (see above comment)?
        if k % 2 != 0:
            (u, v) = (c.numerator(), c.denominator())
            if (2*v) % N == 0:
                if (2*u) % v.gcd(N) == 0:
                    self._is_zero.append(len(g)-1)
                    return self(0)

        if sign == -1:
            # new cusp and nonzero sign, so its possible that the cusp
            # class is killed by the sign relations.
            t, eps = self._is_equiv(c, -c)
            if t:
                if sign != eps:
                    self._is_zero.append(len(g)-1)
                    return self(0)

        return BoundarySpaceElement(self, {(len(g)-1):1})


class BoundarySpace_wtk_eps(BoundarySpace):
    def __init__(self, eps, weight, sign=0):
        """
        Space of boundary modular symbols with given weight, character, and sign.

        INPUT:
            eps -- dirichlet.DirichletCharacter, the "Nebentypus" character.
            weight -- int, the weight >= 2
            sign -- int, either -1, 0, or 1
        EXAMPLES:

        """
        level = eps.modulus()
        sign = int(sign)
        self.__eps = eps
        if not sign in [-1,0,1]:
            raise ArithmeticError, "sign must be an int in [-1,0,1]"
        if level <= 0:
            raise ArithmeticError, "level must be positive"
        BoundarySpace.__init__(self,
                weight = weight,
                group = congroup.Gamma1(level),
                sign = sign,
                base_ring = eps.base_ring(),
                character = eps)

    def _repr_(self):
        return ("Boundary Modular Symbols space of level %s, weight %s, character %s " + \
                "and dimension %s over %s")%(self.level(), self.weight(),
                    self.character(), self.rank(), self.base_ring())


    def _is_equiv(self, c1, c2):
        return c1.is_gamma0_equiv(c2, self.level(), transformation=True)

    def _cusp_index(self, cusp):
        g = self._known_gens
        N = self.level()
        for i in xrange(len(g)):
            t, s = self._is_equiv(cusp, g[i])
            if t:
                return i, self.__eps(s)
        return -1, 0

    def _coerce_cusp(self, c):
        """
        Coerce symbol into a boundary symbol space.
        """
        N    = self.level()
        k    = self.weight()
        sign = self.sign()
        i, eps = self._cusp_index(c)
        if i != -1:
            if i in self._is_zero:
                return self(0)
            return BoundarySpaceElement(self, {i : eps})

        if sign != 0:
            i2, eps = self._cusp_index(-c)
            if i2 != -1:
                if i2 in self._is_zero:
                    return self(0)
                return BoundarySpaceElement(self, {i2:sign*eps})

        # found a new cusp class
        g = self._known_gens
        g.append(c)
        self._known_gens_repr.append("[%s]"%c)

        #############################################################
        # Does cusp class vanish because of the character relations
        # (see Prop 2.30 of Stein, Ph.D. thesis)?
        #
        #      TODO?: This is a very dumb way to check for solutions
        #      to an equation (seep Prop 2.30 for which equation);
        #      however, computing the cusp equivalence for the
        #      boundary map takes much less time than computing the
        #      kernel of the boundary map, so it's not worth
        #      optimizing this now.
        #############################################################

        (u, v) = (c.numerator(), c.denominator())
        gcd = arith.gcd
        d = gcd(v,N)
        x = N//d

        for j in range(d):
            alpha = 1 - j*x
            if gcd(alpha, N) == 1:
                if (v*(1-alpha))%N == 0 and (u*(1-alpha))%d == 0:
                    if self.__eps(alpha) != 1:
                        self._is_zero.append(len(g)-1)
                        return self(0)

        if sign != 0:
            # new cusp and nonzero sign, so its possible that the cusp
            # class is killed by the sign relations.
            t, s = self._is_equiv(c, -c)
            if t:
                if sign != self.__eps(s):
                    self._is_zero.append(len(g)-1)
                    return self(0)

        return BoundarySpaceElement(self, {(len(g)-1):1})


