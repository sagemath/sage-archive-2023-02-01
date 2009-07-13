"""
Morphisms of free modules.

AUTHOR:
    - William Stein

TESTS::

    sage: V = ZZ^2; f = V.hom([V.1,-2*V.0])
    sage: loads(dumps(f))
    Free module morphism defined by the matrix
    [ 0  1]
    [-2  0]
    Domain: Ambient free module of rank 2 over the principal ideal domain ...
    Codomain: Ambient free module of rank 2 over the principal ideal domain ...
    sage: loads(dumps(f)) == f
    True
"""

####################################################################################
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
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
####################################################################################

# A matrix morphism is a morphism that is defined by multiplication by a
# matrix.  Elements of domain must either have a method "vector()" that
# returns a vector that the defining matrix can hit from the left, or
# be coercible into vector space of appropriate dimension.

import sage.matrix.all as matrix
import sage.misc.misc as misc
import sage.modules.free_module as free_module
import matrix_morphism

import free_module_homspace

def is_FreeModuleMorphism(x):
    """
    EXAMPLES::

        sage: V = ZZ^2; f = V.hom([V.1,-2*V.0])
        sage: sage.modules.free_module_morphism.is_FreeModuleMorphism(f)
        True
        sage: sage.modules.free_module_morphism.is_FreeModuleMorphism(0)
        False
    """
    return isinstance(x, FreeModuleMorphism)

class FreeModuleMorphism(matrix_morphism.MatrixMorphism):
    def __init__(self, parent, A):
        """
        INPUT:

            -  ``parent`` - a homspace in a (sub) category of free modules

            -  ``A`` - matrix

        EXAMPLES::

            sage: V = QQ^3; W = span([[1,2,3],[-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ,3,[1..9]))
            sage: type(phi)
            <class 'sage.modules.free_module_morphism.FreeModuleMorphism'>
        """
        if not free_module_homspace.is_FreeModuleHomspace(parent):
            raise TypeError, "parent (=%s) must be a free module hom space"%parent
        if isinstance(A, matrix_morphism.MatrixMorphism):
            A = A.matrix()
        A = parent._matrix_space()(A)
        matrix_morphism.MatrixMorphism.__init__(self, parent, A)

    def __call__(self, x):
        """
        Evaluate this matrix morphism at x, which is either an element
        that can be coerced into the domain or a submodule of the domain.

        EXAMPLES::

            sage: V = QQ^3; W = span([[1,2,3],[-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ,3,[1..9]))

        We compute the image of some elements::

            sage: phi(V.0)
            (1, 2, 3)
            sage: phi(V.1)
            (4, 5, 6)
            sage: phi(V.0  - 1/4*V.1)
            (0, 3/4, 3/2)

        We compute the image of a *subspace*::

            sage: V = QQ^3; W = span([[1,2,3],[-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ,3,[1..9]))
            sage: phi.rank()
            2
            sage: phi(V)
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]

        We restrict phi to W and compute the image of an element::

            sage: psi = phi.restrict_domain(W)
            sage: psi(W.0) == phi(W.0)
            True
            sage: psi(W.1) == phi(W.1)
            True

        We compute the image of a submodule of a ZZ-module embedded in
        a rational vector space::

            sage: V = QQ^3; W = V.span_of_basis([[2,2,3],[-1,2,5/3]], ZZ)
            sage: phi = W.hom([W.0, W.0-W.1]); phi
            Free module morphism defined by the matrix
            [ 1  0]
            [ 1 -1]...
            sage: phi(span([2*W.1],ZZ))
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [  6   0 8/3]
            sage: phi(2*W.1)
            (6, 0, 8/3)
        """
        if free_module.is_FreeModule(x):
            V = self.domain().submodule(x)
            return self.restrict_domain(V).image()
        return matrix_morphism.MatrixMorphism.__call__(self, x)

    def _repr_(self):
        """
        Return string representation of this morphism of free modules.

        EXAMPLES::

            sage: V = QQ^3; W = span([[1,2,3],[-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ,3,[1..9]))
            sage: phi._repr_()
            'Free module morphism defined by the matrix\n[1 2 3]\n[4 5 6]\n[7 8 9]\nDomain: Vector space of dimension 3 over Rational Field\nCodomain: Vector space of dimension 3 over Rational Field'
        """
        if max(self.matrix().nrows(),self.matrix().ncols()) > 5:
            mat = "(not printing %s x %s matrix)"%(self.matrix().nrows(), self.matrix().ncols())
        else:
            mat = str(self.matrix())
        return "Free module morphism defined by the matrix\n%s\nDomain: %s\nCodomain: %s"%(\
            mat, misc.strunc(self.domain()), misc.strunc(self.codomain()))

    def change_ring(self, R):
        """
        Change the ring over which this morphism is defined.  This changes the ring of the
        domain, codomain, and underlying matrix.

        EXAMPLES::

            sage: V0 = span([[0,0,1],[0,2,0]],ZZ); V1 = span([[1/2,0],[0,2]],ZZ); W = span([[1,0],[0,6]],ZZ)
            sage: h = V0.hom([-3*V1.0-3*V1.1, -3*V1.0-3*V1.1])
            sage: h.base_ring()
            Integer Ring
            sage: h
            Free module morphism defined by the matrix
            [-3 -3]
            [-3 -3]...
            sage: h.change_ring(QQ).base_ring()
            Rational Field
            sage: f = h.change_ring(QQ); f
            Free module morphism defined by the matrix
            [-3 -3]
            [-3 -3]
            Domain: Vector space of degree 3 and dimension 2 over Rational Field
            Basis ...
            Codomain: Vector space of degree 2 and dimension 2 over Rational Field
            Basis ...
            sage: f = h.change_ring(GF(7)); f
            Free module morphism defined by the matrix
            [4 4]
            [4 4]
            Domain: Vector space of degree 3 and dimension 2 over Finite Field of ...
            Codomain: Vector space of degree 2 and dimension 2 over Finite Field of ...
        """
        D = self.domain().change_ring(R)
        C = self.codomain().change_ring(R)
        A = self.matrix().change_ring(R)
        return D.hom(A, C)

    def inverse_image(self, V):
        """
        Given a submodule V of the codomain of self, return the
        inverse image of V under self, i.e., the biggest submodule of
        the domain of self that maps into V.

        EXAMPLES:

        We test computing inverse images over a field::

            sage: V = QQ^3; W = span([[1,2,3],[-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ,3,[1..9]))
            sage: phi.rank()
            2
            sage: I = phi.inverse_image(W); I
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    0]
            [   0    1 -1/2]
            sage: phi(I.0) in W
            True
            sage: phi(I.1) in W
            True
            sage: W = phi.image()
            sage: phi.inverse_image(W) == V
            True

        We test computing inverse images between two spaces embedded in different
        ambient spaces.::

            sage: V0 = span([[0,0,1],[0,2,0]],ZZ); V1 = span([[1/2,0],[0,2]],ZZ); W = span([[1,0],[0,6]],ZZ)
            sage: h = V0.hom([-3*V1.0-3*V1.1, -3*V1.0-3*V1.1])
            sage: h.inverse_image(W)
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [0 2 1]
            [0 0 2]
            sage: h(h.inverse_image(W)).is_submodule(W)
            True
            sage: h(h.inverse_image(W)).index_in(W)
            +Infinity
            sage: h(h.inverse_image(W))
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 3 12]


        We test computing inverse images over the integers::

            sage: V = QQ^3; W = V.span_of_basis([[2,2,3],[-1,2,5/3]], ZZ)
            sage: phi = W.hom([W.0, W.0-W.1])
            sage: Z = W.span([2*W.1]); Z
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [    2    -4 -10/3]
            sage: Y = phi.inverse_image(Z); Y
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [  6   0 8/3]
            sage: phi(Y) == Z
            True
        """
        if self.rank() == 0:
            # Special case -- if this is the 0 map, then the only possibility
            # for the inverse image is that it is the whole domain.
            return self.domain()

        R = self.base_ring()
        A = self.matrix()

        # Replace the module V that we are going to pullback by a
        # submodule that is contained in the image of self, since our
        # plan is to lift all generators of V.
        V = self.image().intersection(V)
        # Write V in terms of the basis for the codomain.
        V = self.codomain().coordinate_module(V)
        B = V.basis_matrix()

        # Compute the kernel, which is contained in the inverse image.
        K = self.kernel()

        if R.is_field():
            # By solving, find lifts of each of the basis elements of V.
            # Each row of C gives a linear combination of the basis for the domain
            # that maps to one of the basis elements V.
            C = A.solve_left(B)

        else:
            if not hasattr(A, 'hermite_form'):
                raise NotImplementedError, "base ring (%s) must have hermite_form algorithm in order to compute inverse image"%R

            # 1. Compute H such that U*A = H = hnf(A) without zero
            # rows. What this "does" is find a basis for the image of
            # A and explicitly represents each element in this basis
            # as the image of some element of the domain (the rows of
            # U give these elements of the domain).
            H, U = A.hermite_form(transformation=True,include_zero_rows=False)

            # 2. Next we find the unique solution to the equation
            #    Y*H = B.  This writes each basis element of V in
            # terms of our image basis found in the previous step.
            Y = H.solve_left(B)

            # 3. Multiply Y by U then takes those same linear combinations
            # from step 2 above and lifts them to coefficients that define
            # linear combinations of the basis for the domain.
            C = Y*U

        # Finally take the linear combinations of the basis for the
        # domain defined by C. Together with the kernel K, this spans
        # the inverse image of V.
        if self.domain().is_ambient():
            L = C.row_module(R)
        else:
            L = (C*self.domain().basis_matrix()).row_module(R)

        return K + L

    def lift(self, x):
        r"""
        Given an element of the image, return an element of the codomain that maps onto it.

        EXAMPLE::

            sage: X = QQ**2
            sage: V = X.span([[2, 0], [0, 8]], ZZ)
            sage: W = (QQ**1).span([[1/6]], ZZ)
            sage: f = V.hom([W([1/3]), W([1/2])], W)
            sage: f.lift([1/3])
            (8, -16)
            sage: f.lift([1/2])
            (12, -24)
            sage: f.lift([1/6])
            (4, -8)
            sage: f.lift([1/12])
            Traceback (most recent call last):
            ...
            TypeError: element (= [1/12]) is not in free module

        """
        from free_module_element import vector
        x = self.codomain()(x)
        A = self.matrix()
        H, U = A.hermite_form(transformation=True,include_zero_rows=False)
        Y = H.solve_left(vector(self.codomain().coordinates(x)))
        C = Y*U
        t = self.domain().linear_combination_of_basis(C.row(0))
        assert self(t) == x
        return t

