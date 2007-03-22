r"""
Homology of modular abelian varieties.

\sage can compute with homology groups associated to modular abelian
varieties with coefficients in any commutative ring.  Supported
operations include computing matrices and characteristic polynomials
of Hecke operators, rank, and rational decomposition as a direct sum
of factors (obtained by cutting out kernels of Hecke operators).

AUTHOR:
    -- William Stein (2007-03)


EXAMPLES:
    sage: J = J0(43)
    sage: H = J.integral_homology()
    sage: H
    Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
    sage: H.hecke_matrix(19)
    [ 0  0 -2  0  2  0]
    [ 2 -4 -2  0  2  0]
    [ 0  0 -2 -2  0  0]
    [ 2  0 -2 -4  2 -2]
    [ 0  2  0 -2 -2  0]
    [ 0  2  0 -2  0  0]
    sage: H.base_ring()
    Integer Ring
    sage: d = H.decomposition(); d
    [
    Submodule of rank 2 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43),
    Submodule of rank 4 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
    ]
    sage: a = d[0]
    sage: a.hecke_matrix(5)
    [-4  0]
    [ 0 -4]
    sage: a.T(7)
    Hecke operator T_7 on Submodule of rank 2 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################


from sage.modular.hecke.all import HeckeModule_free_module
from sage.rings.all import Integer, ZZ, QQ


class Homology(HeckeModule_free_module):
    def hecke_polynomial(self, n, var):
        return self.hecke_matrix(n).charpoly(var)


class Homology_abvar(Homology):
    def __init__(self, abvar, base):
        HeckeModule_free_module.__init__(
            self, base, abvar.level(), weight=2)
        self.__abvar = abvar

    def _repr_(self):
        return "Homology of %s"%self.__abvar

    def abelian_variety(self):
        """
        Return the abelian variety that this is the homology of.

        EXAMPLES:
            sage: H = J0(48).homology()
            sage: H.abelian_variety()
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(48)
        """
        return self.__abvar

    def ambient_hecke_module(self):
        """
        Return the ambient Hecke module that this homology is contained in.

        EXAMPLES:
            sage: H = J0(48).homology(); H
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(48)
            sage: H.ambient_hecke_module()
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(48)
        """
        return self

    def free_module(self):
        """
        Return the underlying free module of the homology group.

        EXAMPLES:
            sage: H = J0(48).homology()
            sage: H.free_module()
            Ambient free module of rank 6 over the principal ideal domain Integer Ring
        """
        return self.base_ring()**self.rank()

    def hecke_bound(self):
        r"""
        Return bound on the number of Hecke operators needed to generate
        the Hecke algebra as a $\ZZ$-module acting on this space.

        EXAMPLES:
            sage: J0(48).homology().hecke_bound()
            16
            sage: J1(15).homology().hecke_bound()
            4
        """
        return self.__abvar.modular_symbols(sign=1).hecke_bound()

    def hecke_matrix(self, n):
        """
        Return the matrix of the n-th Hecke operator acting on
        this homology group.
        """
        raise NotImplementedError

    def rank(self):
        """
        Return the rank as a module or vector space of this homology group.

        EXAMPLES:
            sage: H = J0(5077).homology(); H
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(5077)
            sage: H.rank()
            844
        """
        return self.__abvar.dimension() * 2

    def submodule(self, U):
        r"""
        Return the submodule of this homology group given by U, which should
        be a submodule of the free module associated to this homology group.

        NOTE: We do not check that U is invariant under all Hecke operators.

        EXAMPLES:
            sage: H = J0(23).homology(); H
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(23)
            sage: F = H.free_module()
            sage: U = F.span([[1,2,3,4]])
            sage: M = H.submodule(U); M
            Submodule of rank 1 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(23)

        Note that the submodule command doesn't actually check that the object
        defined is a homology group or is invariant under the Hecke operators.
        For example, the fairly random $M$ that we just defined is not invariant
        under the Hecke operators, so it is not a Hecke submodule -- it is only
        a $\ZZ$-submodule.
            sage: M.hecke_matrix(3)
            Traceback (most recent call last):
            ...
            ArithmeticError: subspace is not invariant under matrix
        """
        return Homology_submodule(self, U)


class IntegralHomology(Homology_abvar):
    def __init__(self, abvar):
        Homology_abvar.__init__(self, abvar, ZZ)

    def _repr_(self):
        return "Integral Homology of %s"%self.abelian_variety()

    def hecke_matrix(self, n):
        """
        Return the matrix of the n-th Hecke operator acting on
        this homology group.

        EXAMPLES:
            sage: J0(48).homology().hecke_bound()
            16
            sage: t = J1(13).homology().hecke_matrix(3); t
            [ 0  0  2 -2]
            [-2 -2  0  2]
            [-2 -2  0  0]
            [ 0 -2  2 -2]
            sage: t.base_ring()
            Integer Ring
        """
        n = Integer(n)
        return self.abelian_variety()._integral_hecke_matrix(n)

    def hecke_polynomial(self, n, var):
        """
        Return the n-th Hecke polynomial on this rational homology group.

        EXAMPLES:
            sage: f = J0(43).integral_homology().hecke_polynomial(2)
            sage: f.base_ring()
            Integer ring
            sage: factor(f)
            (x + 2)^2 * (x^2 - 2)^2
        """
        n = Integer(n)
        M = self.abelian_variety().modular_symbols(sign=1)
        f = (M.hecke_polynomial(n, var)**2).change_ring(ZZ)
        return f

class RationalHomology(Homology_abvar):
    def __init__(self, abvar):
        Homology_abvar.__init__(self, abvar, QQ)

    def _repr_(self):
        return "Rational Homology of %s"%self.abelian_variety()

    def hecke_matrix(self, n):
        """
        Return the matrix of the n-th Hecke operator acting on
        this homology group.

        EXAMPLES:
            sage: t = J1(13).homology(QQ).hecke_matrix(3); t
            [ 0  0  2 -2]
            [-2 -2  0  2]
            [-2 -2  0  0]
            [ 0 -2  2 -2]
            sage: t.base_ring()
            Rational Field
            sage: t = J1(13).homology(GF(3)).hecke_matrix(3); t
            [0 0 2 1]
            [1 1 0 2]
            [1 1 0 0]
            [0 1 2 1]
            sage: t.base_ring()
            Finite Field of size 3
        """
        n = Integer(n)
        return self.abelian_variety()._rational_hecke_matrix(n)

    def hecke_polynomial(self, n, var='x'):
        """
        Return the n-th Hecke polynomial on this rational homology group.

        EXAMPLES:
            sage: f = J0(43).rational_homology().hecke_polynomial(2)
            sage: f.base_ring()
            Rational field
            sage: factor(f)
            (x + 2)^2 * (x^2 - 2)^2
        """
        n = Integer(n)
        M = self.abelian_variety().modular_symbols(sign=1)
        f = M.hecke_polynomial(n, var)**2
        return f


class Homology_over_base(Homology_abvar):
    def __init__(self, abvar, base_ring):
        Homology_abvar.__init__(self, abvar, base_ring)

    def _repr_(self):
        return "Homology with coefficients in %s of %s"%(self.base_ring(), self.abelian_variety())

    def hecke_matrix(self, n):
        """
        Return the matrix of the n-th Hecke operator acting on
        this homology group.

        EXAMPLES:
            sage: t = J1(13).homology(GF(3)).hecke_matrix(3); t
            [0 0 2 1]
            [1 1 0 2]
            [1 1 0 0]
            [0 1 2 1]
            sage: t.base_ring()
            Finite Field of size 3
        """
        n = Integer(n)
        return self.abelian_variety()._integral_hecke_matrix(n).change_ring(self.base_ring())

class Homology_submodule(Homology):
    def __init__(self, ambient, submodule):
        self.__ambient = ambient
        self.__submodule = submodule
        HeckeModule_free_module.__init__(
            self, ambient.base_ring(), ambient.level(), weight=2)

    def _repr_(self):
        return "Submodule of rank %s of %s"%(self.rank(), self.__ambient)

    def ambient_hecke_module(self):
        """
        Return the ambient Hecke module that this homology is contained in.

        EXAMPLES:
            sage: H = J0(48).homology(); H
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(48)
            sage: d = H.decomposition(); d
            [
            Submodule of rank 2 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(48),
            Submodule of rank 4 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(48)
            ]
            sage: d[0].ambient_hecke_module()
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(48)
        """
        return self.__ambient

    def free_module(self):
        """
        Return the underlying free module of the homology group.

        EXAMPLES:
            sage: H = J0(48).homology()
            sage: K = H.decomposition()[1]; K
            Submodule of rank 4 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(48)
            sage: K.free_module()
            Free module of degree 6 and rank 4 over Integer Ring
            Echelon basis matrix:
            [ 1  0  0  0  0  0]
            [ 0  1  0  0  1 -1]
            [ 0  0  1  0 -1  1]
            [ 0  0  0  1  0 -1]
        """
        return self.__submodule

    def hecke_bound(self):
        """
        Return a bound on the number of Hecke operators needed to
        generate the Hecke algebra acting on this homology group.

        EXAMPLES:
            sage: d = J0(43).homology().decomposition(2); d
            [
            Submodule of rank 2 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43),
            Submodule of rank 4 of Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            ]

        Because this factor has dimension 2 it corresponds to an elliptic curve,
        so we have a Hecke bound of 1.
            sage: d[0].hecke_bound()
            1
            sage: d[1].hecke_bound()
            8
        """
        if self.rank() <= 2:
            return ZZ(1)
        else:
            return self.__ambient.hecke_bound()

    def hecke_matrix(self, n):
        """
        Return the matrix of the n-th Hecke operator acting on
        this homology group.

        EXAMPLES:
            sage: d = J0(125).homology(GF(17)).decomposition(2); d
            [
            Submodule of rank 4 of Homology with coefficients in Finite Field of size 17 of Jacobian of the modular curve associated to the congruence subgroup Gamma0(125),
            Submodule of rank 4 of Homology with coefficients in Finite Field of size 17 of Jacobian of the modular curve associated to the congruence subgroup Gamma0(125),
            Submodule of rank 8 of Homology with coefficients in Finite Field of size 17 of Jacobian of the modular curve associated to the congruence subgroup Gamma0(125)
            ]
            sage: t = d[0].hecke_matrix(17); t
            [16 15 15  0]
            [ 0  5  0  2]
            [ 2  0  5 15]
            [ 0 15  0 16]
            sage: t.base_ring()
            Finite Field of size 17
            sage: t.fcp()
            (x^2 + 13*x + 16)^2
        """
        n = Integer(n)
        try:
            return self.__hecke_matrix[n]
        except AttributeError:
            self.__hecke_matrix = {}
        except KeyError:
            pass
        t = self.__ambient.hecke_matrix(n)
        s = t.restrict(self.__submodule)
        self.__hecke_matrix[n] = s
        return s

    def rank(self):
        """
        Return the rank of this homology group.

        EXAMPLES:
            sage: d = J0(43).homology().decomposition(2)
            sage: [H.rank() for H in d]
            [2, 4]
        """
        return self.__submodule.rank()


