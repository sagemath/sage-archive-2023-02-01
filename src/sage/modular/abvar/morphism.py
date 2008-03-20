r"""
Morphisms between modular abelian varieties, including Hecke operators
acting on modular abelian varieties.

\sage can compute with Hecke operators on modular abelian varieties.
A Hecke operator is defined by given a modular abelian variety and an
index.  Given a Hecke operator, \sage can compute the characteristic
polynomial, and the action of the Hecke operator on various homology
groups.

AUTHOR:
    -- William Stein (2007-03)
    -- Craig Citro (2008-03)

EXAMPLES:
    sage: A = J0(54)
    sage: t5 = A.hecke_operator(5); t5
    Hecke operator T_5 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(54)
    sage: t5.charpoly().factor()
    (x - 3)^2 * (x + 3)^2 * x^4
    sage: B = A.new_subvariety(); B
    Abelian variety factor of dimension 2 of J0(54)
    sage: t5 = B.hecke_operator(5); t5
    Hecke operator T_5 on Abelian variety factor of dimension 2 of J0(54)
    sage: t5.charpoly().factor()
    (x - 3)^2 * (x + 3)^2
    sage: t5.action_on_homology().matrix()
    [ 0  3  3 -3]
    [-3  3  3  0]
    [ 3  3  0 -3]
    [-3  6  3 -3]
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.rings.all import ZZ, QQ
import abvar as abelian_variety
import sage.modules.matrix_morphism
import sage.structure.element

from finite_subgroup import FiniteSubgroupElement

class Morphism_abstract(sage.modules.matrix_morphism.MatrixMorphism_abstract):
    """
    A morphism between modular abelian varieties.
    EXAMPLES:
        sage: t = J0(11).hecke_operator(2)
        sage: from sage.modular.abvar.morphism import Morphism
        sage: isinstance(t, Morphism)
        True
    """
    def complementary_isogeny(self):
        """
        Returns the complementary isogeny of self.

        EXAMPLES:
            sage: J = J0(43)
            sage: A = J[1]
            sage: T5 = A.hecke_operator(5)
            sage: T5.is_isogeny()
            True
            sage: T5.complementary_isogeny()
            Morphism defined by the matrix
            [ 1 -1  0  0]
            [-1  3  0  0]
            [-1  0  3  1]
            [ 0 -1  1  1]
            sage: T5.complementary_isogeny() * T5
            Morphism defined by the matrix
            [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2]
        """
        if not self.is_isogeny():
            raise ValueError, "self is not an isogeny"
        M = self.matrix()
        try:
            iM, denom = M._invert_iml()
        except AttributeError:
            iM = M.matrix_over_field().invert()
            iM, denom = iM._clear_denom()
        return Morphism(self.parent().reversed(), iM)

    def is_isogeny(self):
        """
        Return True if this morphism is an isogeny of abelian varieties.

        EXAMPLES:
            sage: J = J0(39)
            sage: Id = J.hecke_operator(1)
            sage: Id.is_isogeny()
            True
            sage: J.hecke_operator(19).is_isogeny()
            False
        """
        M = self.matrix()
        return M.nrows() == M.ncols() == M.rank()

    def cokernel(self):
        try:
            return self.__cokernel
        except AttributeError:
            I = self.image()
            C = self.codomain().quotient(I)
            self.__cokernel = C
            return C


    def kernel(self):
        """
        Return the kernel of this morphism.

        OUTPUT:
            G -- a finite group
            A -- an abelian variety (identity component of the kernel)

        EXAMPLES:
        We compute the kernel of a projection map.  Notice that the
        kernel has a nontrivial abelian variety part.
            sage: A, B, C = J0(33)
            sage: pi = J0(33).projection(B)
            sage: pi.kernel()
            (Finite subgroup with invariants [20] over QQbar of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 2 of J0(33))

        We compute the kernels of some Hecke operators:
            sage: t2 = J0(33).hecke_operator(2)
            sage: t2
            Hecke operator T_2 on Abelian variety J0(33) of dimension 3
            sage: t2.kernel()
            (Finite subgroup with invariants [2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 0 of J0(33))
            sage: t3 = J0(33).hecke_operator(3)
            sage: t3.kernel()
            (Finite subgroup with invariants [3, 3] over QQ of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 0 of J0(33))
        """
        A = self.matrix()
        L = A.image().change_ring(ZZ)
        # Saturate the image of the matrix corresponding to self.
        Lsat = L.saturation()
        # Now find a matrix whose rows map exactly onto the
        # saturation of L.
        X = A.solve_left(Lsat.basis_matrix())
        D = self.domain()
        V = (A.kernel().basis_matrix() * D.vector_space().basis_matrix()).row_module()
        Lambda = V.intersection(D._ambient_lattice())
        from abvar import ModularAbelianVariety
        abvar = ModularAbelianVariety(D.groups(), Lambda, D.base_ring())

        if Lambda.rank() == 0:
            field_of_definition = QQ
        else:
            field_of_definition = None
        K = D.finite_subgroup(X.rows(), field_of_definition=field_of_definition)

        return K, abvar

    def factor_out_component_group(self):
        r"""
        View self as a morphism $f:A \to B$.  Then $\ker(f)$ is an
        extension of an abelian variety $C$ by a finite component
        group $G$.  This function constructs a morphism $g$ with
        domain $A$ and codomain Q isogenous to $C$ such that $\ker(g)$
        is equal to $C$.

        OUTPUT:
            a morphism

        EXAMPLES:
            sage: A,B,C = J0(33)
            sage: pi = J0(33).projection(A)
            sage: pi.kernel()
            (Finite subgroup with invariants [5] over QQbar of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 2 of J0(33))
            sage: psi = pi.factor_out_component_group()
            sage: psi.kernel()
            (Finite subgroup with invariants [] over QQbar of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 2 of J0(33))


        ALGORITHM: We compute a subgroup $G$ of $B$ so that the
        composition $h: A\to B \to B/G$ has kernel that contains
        $A[n]$ and component group isomorphic to $(\ZZ/n\ZZ)^{2d}$,
        where $d$ is the dimension of $A$.  Then $h$ factors through
        multiplication by $n$, so there is a morphism $g: A\to B/G$
        such that $g \circ [n] = h$.  Then $g$ is the desired
        morphism.  We give more details below about how to transform
        this into linear algebra.
        """
        try:
            return self.__factor_out
        except AttributeError:
            A = self.matrix()
            L = A.image()
            # Saturate the image of the matrix corresponding to self.
            Lsat = L.saturation()
            if L == Lsat: # easy case
                self.__factor_out = self
                return self
            # Now find a matrix whose rows map exactly onto the
            # saturation of L.
            X = A.solve_left(Lsat.basis_matrix())

            # Find an integer n such that 1/n times the lattice Lambda of A
            # contains the row span of X.
            n = X.denominator()

            # Let M be the lattice of the codomain B of self.
            # Now 1/n * Lambda contains Lambda and maps
            # via the matrix of self to a lattice L' that
            # contains Lsat.   Consider the lattice
            #                R = M + L'.
            # This is a lattice that contains the lattice M of B.
            # Also 1/n*Lambda maps exactly to L' in R.
            # We have
            #     R/L' = (M+L')/L' = M/(L'/\M) = M/Lsat
            # which is torsion free!

            Q      = self.codomain()
            M      = Q.lattice()
            one_over_n = ZZ(1)/n
            Lprime = (one_over_n * self.matrix() * M.basis_matrix()).row_module(ZZ)

            # This R is a lattice in the ambient space for B.
            R = Lprime + M

            from abvar import ModularAbelianVariety
            C = ModularAbelianVariety(Q.groups(), R, Q.base_field())

            # We have to change the basis of the representation of A
            # to the basis for R instead of the basis for M.  Each
            # row of A is written in terms of M, but needs to be
            # in terms of R's basis, which contains M with finite index.
            change_basis_from_M_to_R = R.basis_matrix().solve_left(M.basis_matrix())
            matrix = one_over_n * A * change_basis_from_M_to_R

            # Finally
            g = Morphism(self.domain().Hom(C), matrix)
            self.__factor_out = g
            return g

    def image(self):
        """
        Return the image of this morphism.

        OUTPUT:
            an abelian variety

        EXAMPLES:
        We compute the image of projection onto a factor of $J_0(33)$:
            sage: A,B,C = J0(33)
            sage: A
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: f = J0(33).projection(A)
            sage: f.image()
            Abelian subvariety of dimension 1 of J0(33)
            sage: f.image() == A
            True

        We compute the image of a Hecke operator:
            sage: t2 = J0(33).hecke_operator(2); t2.fcp()
            (x - 1)^2 * (x + 2)^4
            sage: phi = t2 + 2
            sage: phi.image()
            Abelian subvariety of dimension 1 of J0(33)

        The sum of the image and the kernel is the whole space:
            sage: phi.kernel()[1] + phi.image() == J0(33)
            True
        """
        return self(self.domain())

    def __call__(self, X):
        """
        INPUT:
            X -- abelian variety, finite group, or torsion element

        OUTPUT:
            abelian variety, finite group, torsion element

        EXAMPLES:
        We apply morphisms to elements:
            sage: t2 = J0(33).hecke_operator(2)
            sage: G  = J0(33).n_torsion_subgroup(2); G
            Finite subgroup with invariants [2, 2, 2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2(G.0)
            [(-1/2, 0, 1/2, -1/2, 1/2, -1/2)]
            sage: t2(G.0) in G
            True
            sage: t2(G.1)
            [(0, -1, 1/2, 0, 1/2, -1/2)]
            sage: t2(G.2)
            [(0, 0, 0, 0, 0, 0)]
            sage: K = t2.kernel()[0]; K
            Finite subgroup with invariants [2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2(K.0)
            [(0, 0, 0, 0, 0, 0)]

        We apply morphisms to subgroups:
            sage: t2 = J0(33).hecke_operator(2)
            sage: G  = J0(33).n_torsion_subgroup(2); G
            Finite subgroup with invariants [2, 2, 2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2(G)
            Finite subgroup with invariants [2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2.fcp()
            (x - 1)^2 * (x + 2)^4

        We apply morphisms to abelian subvarieties:
            sage: E11a0, E11a1, B = J0(33)
            sage: t2 = J0(33).hecke_operator(2)
            sage: t3 = J0(33).hecke_operator(3)
            sage: E11a0
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: t3(E11a0)
            Abelian subvariety of dimension 1 of J0(33)
            sage: t3(E11a0).decomposition()
            [
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)
            ]
            sage: t3(E11a0) == E11a1
            True
            sage: t2(E11a0) == E11a0
            True
            sage: t3(E11a0) == E11a0
            False
            sage: t3(E11a0 + E11a1) == E11a0 + E11a1
            True
        """
        from abvar import is_ModularAbelianVariety
        from finite_subgroup import FiniteSubgroup
        if isinstance(X, FiniteSubgroupElement):
            return self._image_of_element(X)
        elif is_ModularAbelianVariety(X):
            return self._image_of_abvar(X)
        elif isinstance(X, FiniteSubgroup):
            return self._image_of_finite_subgroup(X)
        else:
            raise TypeError, "X must be an abelian variety or finite subgroup"

    def _image_of_element(self, x):
        """
        Return the image of the torsion point x under this morphism.

        The parent of the image element is always the group of all
        torsion elements of the abelian variety.

        INPUT:
            x -- a torsion point on an abelian variety

        OUTPUT:
            a torsion point

        EXAMPLES:
        """
        v = x.element() * self.matrix()
        T = self.codomain().qbar_torsion_subgroup()
        return T(v)

    def _image_of_finite_subgroup(self, G):
        """
        Return the image of the finite group $G$ under this morphism.

        INPUT:
            G -- a finite subgroup of the domain of this morphism

        OUTPUT:
            a finite subgroup of the codomain
        """
        H = [self._image_of_element(x) for x in G.gens()]
        return self.codomain().finite_subgroup(H, field_of_definition = G.field_of_definition())

    def _image_of_abvar(self, A):
        D = self.domain()
        C = self.codomain()
        if A is D:
            B = self.matrix()
        else:
            if not A.is_subvariety(D):
                raise ValueError, "A must be an abelian subvariety of self."
            # Write the vector space corresponding to A in terms of self's
            # vector space, then take the image under self.
            B = D.vector_space().coordinate_module(A.vector_space()).basis_matrix() * self.matrix()

        V = (B * C.vector_space().basis_matrix()).row_module(QQ)

        lattice = V.intersection(C.lattice())
        base_field = C.base_field()
        return abelian_variety.ModularAbelianVariety(C.groups(), lattice, base_field)


class Morphism(Morphism_abstract, sage.modules.matrix_morphism.MatrixMorphism):
    pass

class HeckeOperator(Morphism):
    """
    A Hecke operator acting on a modular abelian variety.
    """
    def __init__(self, abvar, n):
        """
        Create the Hecke operator of index $n$ acting on the abelian
        variety abvar.

        INPUT:
            abvar -- a modular abelian variety
            n -- a positive integer

        EXAMPLES:
            sage: J = J0(37)
            sage: T2 = J.hecke_operator(2); T2
            Hecke operator T_2 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
        """
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n must be positive"
        if not abelian_variety.is_ModularAbelianVariety(abvar):
            raise TypeError, "abvar must be a modular abelian variety"
        self.__abvar = abvar
        self.__n = n
        sage.modules.matrix_morphism.MatrixMorphism_abstract.__init__(self, abvar._Hom_(abvar))

    def _repr_(self):
        """
        String representation of this Hecke operator.

        EXAMPLES:
            sage: J = J0(37)
            sage: J.hecke_operator(2)._repr_()
            'Hecke operator T_2 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)'
        """
        return "Hecke operator T_%s on %s"%(self.__n, self.__abvar)

    def index(self):
        """
        Return the index of this Hecke operator. (For example, if this
        is the operator $T_n$, then the index is the integer $n$.)

        OUTPUT:
            n -- a (Sage) Integer

        EXAMPLES:
            sage: J = J0(15)
            sage: t = J.hecke_operator(53)
            sage: t
            Hecke operator T_53 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(15)
            sage: t.index()
            53
            sage: t = J.hecke_operator(54)
            sage: t
            Hecke operator T_54 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(15)
            sage: t.index()
            54

        TODO
        This is an EXTREMELY long doctest in the current model, but
        instant in the previous model.
            J = J1(12345)
            t = J.hecke_operator(997)
            t
            Hecke operator T_997 on Jacobian of the modular curve associated to the congruence subgroup Gamma1(12345)
            t.index()
            997
            type(t.index())
            <type 'sage.rings.integer.Integer'>
        """
        return self.__n

    def n(self):
        r"""
        Alias for \code{self.index()}.

        EXAMPLES:
            sage: J = J0(17)
            sage: J.hecke_operator(5).n()
            5
        """
        return self.index()

    def characteristic_polynomial(self, var='x'):
        """
        Return the characteristic polynomial of this Hecke operator in
        the given variable.

        INPUT:
            var -- a string (default: 'x')

        OUTPUT:
            a polynomial in var over the rational numbers.

        EXAMPLES:
            sage: A = J0(43)[1]; A
            Abelian variety factor of dimension 2 of J0(43)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Abelian variety factor of dimension 2 of J0(43)
            sage: f = t2.characteristic_polynomial(); f
            x^4 - 4*x^2 + 4
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: f.factor()
            (x^2 - 2)^2
            sage: t2.characteristic_polynomial('y')
            y^4 - 4*y^2 + 4
        """
        return self.__abvar.rational_homology().hecke_polynomial(self.__n, var).change_ring(ZZ)

    def charpoly(self, var='x'):
        r"""
        Synonym for \code{self.characteristic_polynomial(var)}.

        INPUT:
            var -- string (default: 'x')

        EXAMPLES:
            sage: A = J1(13)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Jacobian of the modular curve associated to the congruence subgroup Gamma1(13)
            sage: f = t2.charpoly(); f
            x^4 + 6*x^3 + 15*x^2 + 18*x + 9
            sage: f.factor()
            (x^2 + 3*x + 3)^2
            sage: t2.charpoly('y')
            y^4 + 6*y^3 + 15*y^2 + 18*y + 9
        """
        return self.characteristic_polynomial(var)

    def action_on_homology(self, R=ZZ):
        """
        Return the action of this Hecke operator on the homology
        $H_1(A; R)$ of this abelian variety with coefficients in $R$.

        EXAMPLES:
            sage: A = J0(43)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            sage: h2 = t2.action_on_homology(); h2
            Hecke operator T_2 on Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            sage: h2.matrix()
            [-2  1  0  0  0  0]
            [-1  1  1  0 -1  0]
            [-1  0 -1  2 -1  1]
            [-1  0  1  1 -1  1]
            [ 0 -2  0  2 -2  1]
            [ 0 -1  0  1  0 -1]
            sage: h2 = t2.action_on_homology(GF(2)); h2
            Hecke operator T_2 on Homology with coefficients in Finite Field of size 2 of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            sage: h2.matrix()
            [0 1 0 0 0 0]
            [1 1 1 0 1 0]
            [1 0 1 0 1 1]
            [1 0 1 1 1 1]
            [0 0 0 0 0 1]
            [0 1 0 1 0 1]
        """
        return self.__abvar.homology(R).hecke_operator(self.index())

    def matrix(self):
        try:
            return self._matrix
        except AttributeError:
            pass
        self._matrix = self.action_on_homology().matrix()
        return self._matrix
