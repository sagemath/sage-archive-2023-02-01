"""
Classical Lie Algebras

These are the Lie algebras corresponding to types `A_n`, `B_n`, `C_n`, `D_n`.

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.misc import repr_lincomb
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import RingElement
from sage.structure.indexed_generators import IndexedGenerators
from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras

from sage.algebras.algebra import Algebra
from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra, LieAlgebraFromAssociative
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.arith import binomial
from sage.sets.family import Family

class ClassicalMatrixLieAlgebra(LieAlgebraFromAssociative):
    """
    A classical Lie algebra represented using matrices.
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type):
        """
        Return the correct parent based on input.

        INPUT:

        - ``R`` -- the base ring
        - ``ct`` -- the Cartan type

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.classical_lie_algebra import ClassicalMatrixLieAlgebra
            sage: ClassicalMatrixLieAlgebra(QQ, ['A', 4])
            Special linear Lie algebra of rank 5 over Rational Field
            sage: ClassicalMatrixLieAlgebra(QQ, CartanType(['B',4]))
            Special orthogonal Lie algebra of rank 9 over Rational Field
            sage: ClassicalMatrixLieAlgebra(QQ, 'C4')
            Symplectic Lie algebra of rank 8 over Rational Field
            sage: ClassicalMatrixLieAlgebra(QQ, cartan_type=['D',4])
            Special orthogonal Lie algebra of rank 8 over Rational Field
        """
        cartan_type = CartanType(cartan_type)

        if not cartan_type.is_finite():
            raise ValueError("only for finite types")

        if cartan_type.type() == 'A':
            return sl(R, cartan_type.rank() + 1, 'matrix')
        if cartan_type.type() == 'B':
            return so(R, 2*cartan_type.rank() + 1, 'matrix')
        if cartan_type.type() == 'C':
            return sp(R, 2*cartan_type.rank(), 'matrix')
        if cartan_type.type() == 'D':
            return so(R, 2*cartan_type.rank(), 'matrix')
        raise NotImplementedError("only implemented for types A, B, C, D")

    def __init__(self, R, ct, e, f, h):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: TestSuite(g).run()
        """
        n = len(e)
        names = ['e%s'%i for i in range(1, n+1)]
        names += ['f%s'%i for i in range(1, n+1)]
        names += ['h%s'%i for i in range(1, n+1)]
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        LieAlgebraFromAssociative.__init__(self, e[0].parent(),
                                           gens=tuple(e + f + h),
                                           names=tuple(names),
                                           index_set=tuple(names),
                                           category=category)
        self._cartan_type = ct

        gens = tuple(self.gens())
        i_set = ct.index_set()
        self._e = Family(dict( (i, gens[c]) for c,i in enumerate(i_set) ))
        self._f = Family(dict( (i, gens[n+c]) for c,i in enumerate(i_set) ))
        self._h = Family(dict( (i, gens[2*n+c]) for c,i in enumerate(i_set) ))

    def e(self, i):
        r"""
        Return the generator `e_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.e(2)
            [0 0 0]
            [0 0 1]
            [0 0 0]
        """
        return self._e[i]

    def f(self, i):
        r"""
        Return the generator `f_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.f(2)
            [0 0 0]
            [0 0 0]
            [0 1 0]
        """
        return self._f[i]

    def h(self, i):
        """
        Return the generator `h_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.h(2)
            [ 0  0  0]
            [ 0  1  0]
            [ 0  0 -1]
        """
        return self._h[i]

    @cached_method
    def index_set(self):
        """
        Return the index_set of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.index_set()
            (1, 2)
        """
        return self._cartan_type.index_set()

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def epsilon(self, i, h):
        r"""
        Return the action of the functional
        `\varepsilon_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``h``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.epsilon(1, g.h(1))
            1
            sage: g.epsilon(2, g.h(1))
            -1
            sage: g.epsilon(3, g.h(1))
            0
        """
        return h[i-1,i-1]

    @abstract_method
    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``h``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.simple_root(1, g.h(1))
            2
            sage: g.simple_root(1, g.h(2))
            -1
        """

    def highest_root_basis_elt(self, pos=True):
        r"""
        Return the basis element corresponding to the highest root `\theta`.
        If ``pos`` is ``True``, then returns `e_{\theta}`, otherwise it
        returns `f_{\theta}`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.highest_root_basis_elt()
            [0 0 1]
            [0 0 0]
            [0 0 0]
        """
        RL = self._cartan_type.root_system().root_lattice()
        coroots = RL.simple_coroots()
        theta = RL.highest_root()
        i,w = theta.to_simple_root(True)
        r = RL.simple_root(i)
        if pos:
            gens = self._e
        else:
            gens = self._f
        cur = gens[i]
        for j in w:
            for k in range(-r.scalar(coroots[j])):
                cur = self.bracket(gens[j], cur)
        return cur

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: tuple(g.basis())
            (
            [0 0 0]  [0 0 0]  [ 0  0  0]  [ 1  0  0]  [0 1 0]  [0 0 0]
            [1 0 0]  [0 0 0]  [ 0  1  0]  [ 0 -1  0]  [0 0 0]  [0 0 1]
            [0 0 0], [0 1 0], [ 0  0 -1], [ 0  0  0], [0 0 0], [0 0 0]
            )
        """
        d = {}
        for i in self.index_set():
            d['e{}'.format(i)] = self._e[i]
            d['f{}'.format(i)] = self._f[i]
            d['h{}'.format(i)] = self._h[i]
        return Family(d)

    # TODO: Uncomment once #16825 is done
    #def affine(self, kac_moody=False):
    #    """
    #    Return the (untwisted) affine (Kac-Moody) Lie algebra of ``self``.
    #    """
    #    from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra
    #    return AffineLieAlgebra(self, kac_moody)

class gl(LieAlgebraFromAssociative):
    r"""
    The Lie algebra `\mathfrak{gl}_n` which consists of all `n \times n`
    matrices.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 4)
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        names = []
        gens = []
        for i in range(n):
            for j in range(n):
                names.append('E_{0}_{1}'.format(i,j))
                mat = MS({(i,j):one})
                mat.set_immutable()
                gens.append(mat)
        self._n = n
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        LieAlgebraFromAssociative.__init__(self, MS, tuple(gens),
                                           names=tuple(names),
                                           index_set=tuple(names),
                                           category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.gl(QQ, 4)
            General linear Lie algebra of rank 4 over Rational Field
        """
        return "General linear Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{gl}_n` is:

        .. MATH::

            \langle x \mid y \rangle = 2n \mathrm{tr}(xy) - 2 \mathrm{tr}(x)
            \mathrm{tr}(y).

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 4)
            sage: g.killing_form()
        """
        return 2 * self._n * (x.value * y.value).trace() \
            - 2 * x.value.trace() * y.value.trace()

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 2)
            sage: g.basis()
            (
            [1 0]  [0 1]  [0 0]  [0 0]
            [0 0], [0 0], [1 0], [0 1]
            )
        """
        return self.gens()

class sl(ClassicalMatrixLieAlgebra):
    r"""
    The Lie algebra `\mathfrak{sl}_n` which consists of `n \times n` matrices
    with trace 0. This is the Lie algebra of type `A_{n-1}`.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    - ``representation`` -- (default: ``'bracket'``) can be one of
      the following:

      * ``'bracket'`` - use brackets and the Chevalley basis
      * ``'matrix'`` - use matrices

    EXAMPLES:

    We first construct `\mathfrak{sl}_2` using the Chevalley basis::

        sage: sl2 = lie_algebras.sl(QQ, 2); sl2
        Lie algebra of ['A', 1] in the Chevalley basis
        sage: E,F,H = sl2.gens()
        sage: E.bracket(F) == H
        True
        sage: H.bracket(E) == 2*E
        True
        sage: H.bracket(F) == -2*F
        True

    We now construct `\mathfrak{sl}_2` as a matrix Lie algebra::

        sage: sl2 = lie_algebras.sl(QQ, 2, representation='matrix')
        sage: E,F,H = sl2.gens()
        sage: E.bracket(F) == H
        True
        sage: H.bracket(E) == 2*E
        True
        sage: H.bracket(F) == -2*F
        True
    """
    # FIXME: Move this and the same methods for all other types into separate
    #   functions for the Lie algebra examples because it breaks pickling
    @staticmethod
    def __classcall_private__(cls, R, n, representation='bracket'):
        """
        Choose the representation.
        """
        if representation == 'bracket':
            return LieAlgebraChevalleyBasis(R, ['A', n-1])
        if representation == 'matrix':
            return super(sl, cls).__classcall__(cls, R, n)
        raise ValueError("invalid representation")

    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        e = [MS({(i,i+1):one}) for i in range(n-1)]
        f = [MS({(i+1,i):one}) for i in range(n-1)]
        h = [MS({(i,i):one, (i+1,i+1):-one}) for i in range(n-1)]
        self._n = n
        ClassicalMatrixLieAlgebra.__init__(self, R, CartanType(['A', n-1]), e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.sl(QQ, 5, representation='matrix')
            Special linear Lie algebra of rank 5 over Rational Field
        """
        return "Special linear Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{sl}_n` is:

        .. MATH::

            \langle x \mid y \rangle = 2n \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
        """
        return 2 * self._n * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]
        """
        i = self.index_set().index(i)
        return h[i,i] - h[i+1,i+1]

class so(ClassicalMatrixLieAlgebra):
    r"""
    The Lie algebra `\mathfrak{so}_n` which consists of orthogonal `n \times n`
    matrices. This is the Lie algebra of type `B_{(n-1)/2}` or `D_{n/2}` if `n`
    is odd or even respectively.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    - ``representation`` -- (default: ``'bracket'``) can be one of
      the following:

      * ``'bracket'`` - use brackets and the Chevalley basis
      * ``'matrix'`` - use matrices
    """
    @staticmethod
    def __classcall_private__(cls, R, n, representation='bracket'):
        """
        Choose the representation.
        """
        if representation == 'bracket':
            if n % 2 == 0:
                return LieAlgebraChevalleyBasis(R, ['D', n//2])
            else:
                return LieAlgebraChevalleyBasis(R, ['B', (n-1)//2])
        if representation == 'matrix':
            return super(so, cls).__classcall__(cls, R, n)
        raise ValueError("invalid representation")

    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: TestSuite(g).run()
            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        self._n = n
        if n % 2 == 0: # Even
            m = n / 2 - 1 # -1 for indexing
            n -= 1
            e = [MS({(m-1,n):one, (m,n-1):-one})]
            f = [MS({(n,m-1):one, (n-1,m):-one})]
            h = [MS({(m-1,m-1):one, (m,m):one, (n-1,n-1):-one, (n,n):-one})]
            m += 1
            ct = CartanType(['D', m])
        else: # Odd
            m = (n-1) / 2 - 1 # -1 for indexing
            n -= 1
            e = [MS({(m,n):2, (n,n-1):-2})]
            f = [MS({(n,m):one, (n-1,n):-one})]
            h = [MS({(m,m):2, (n-1,n-1):-2})]
            m += 1
            ct = CartanType(['B', m])
        e = [MS({(i,i+1):one, (m+i+1,m+i):-one}) for i in range(m-1)] + e
        f = [MS({(i+1,i):one, (m+i,m+i+1):-one}) for i in range(m-1)] + f
        h = [MS({(i,i):one, (i+1,i+1):-one, (m+i,m+i):-one, (m+i+1,m+i+1):one}) for i in range(m-1)] + h
        ClassicalMatrixLieAlgebra.__init__(self, R, ct, e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['B', 4], representation='matrix')
            Special orthogonal Lie algebra of rank 9 over Rational Field
            sage: LieAlgebra(QQ, cartan_type=['D', 4], representation='matrix')
            Special orthogonal Lie algebra of rank 8 over Rational Field
        """
        return "Special orthogonal Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{so}_n` is:

        .. MATH::

            \langle x \mid y \rangle = (n - 2) \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
        """
        return 2 * self._n * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES:

        The even or type `D` case::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1 -1]
            [ 0 -1  2  0]
            [ 0 -1  0  2]

        The odd or type `B` case::

            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -2  2]
        """
        i = self.index_set().index(i)
        if i == len(self.index_set()) - 1:
            if self._n % 2 == 0:
                return h[i-1,i-1] + h[i,i]
            # otherwise we are odd
            return h[i,i]
        return h[i,i] - h[i+1,i+1]

class sp(ClassicalMatrixLieAlgebra):
    r"""
    The Lie algebra `\mathfrak{sp}_{2k}` which consists of `2k \times 2k`
    matrices `X` that satisfy the equation:

    .. MATH::

        X^T M - M X = 0

    where

    .. MATH::

        M = \begin{pmatrix}
        0 & I_k \\
        -I_k & 0
        \end{pmatrix}.

    This is the Lie algebra of type `C_k`.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    - ``representation`` -- (default: ``'bracket'``) can be one of
      the following:

      * ``'bracket'`` - use brackets and the Chevalley basis
      * ``'matrix'`` - use matrices
    """
    @staticmethod
    def __classcall_private__(cls, R, n, representation='bracket'):
        """
        Choose the representation.
        """
        if n % 2 != 0:
            raise ValueError("n must be even")

        if representation == 'bracket':
            return LieAlgebraChevalleyBasis(R, ['C', n//2])
        if representation == 'matrix':
            return super(sp, cls).__classcall__(cls, R, n)
        raise ValueError("invalid representation")

    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        self._n = n
        n = n // 2
        e = [MS({(i,i+1):one, (n+i+1,n+i):-one}) for i in range(n-1)]
        e.append(MS({(n-1,2*n-1):one})) # -1 for indexing
        f = [MS({(i+1,i):one, (n+i,n+i+1):-one}) for i in range(n-1)]
        f.append(MS({(2*n-1,n-1):one})) # -1 for indexing
        h = [MS({(i,i):one, (i+1,i+1):-one, (n+i,n+i):-one, (n+i+1,n+i+1):one}) for i in range(n-1)]
        h.append(MS({(n-1,n-1):one, (2*n-1,2*n-1):-one})) # -1 for indexing
        ClassicalMatrixLieAlgebra.__init__(self, R, CartanType(['C', n]), e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.sp(QQ, 8, representation='matrix')
            Symplectic Lie algebra of rank 8 over Rational Field
        """
        return "Symplectic Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{sp}_n` is:

        .. MATH::

            \langle x \mid y \rangle = (2n + 2) \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
        """
        return (2 * self._n + 2) * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -2]
            [ 0  0 -1  2]
        """
        i = self.index_set().index(i)
        if i == self._n / 2 - 1:
            return 2*h[i,i]
        return h[i,i] - h[i+1,i+1]


#######################################
## Chevalley Basis

class LieAlgebraChevalleyBasis(LieAlgebraWithStructureCoefficients):
    r"""
    A simple finite dimensional Lie algebra in the Chevalley basis.

    Let `L` be a simple complex Lie algebra with roots `\Phi`, then the
    Chevalley basis is given by `e_{\alpha}` for all `\alpha \in \Phi` and
    `h_{\alpha_i} := h_i` where `\alpha_i` is a simple root subject. These
    generators are subject to the relations:

    .. MATH::

        \begin{aligned}
        [h_i, h_j] & = 0
        \\ [h_i, e_{\beta}] & = A_{\alpha_i, \beta} e_{\beta}
        \\ [e_{\beta}, e_{-\beta}] & = \sum_i A_{\beta, \alpha_i} h_i
        \\ [e_{\beta}, e_{\gamma}] & = \begin{cases}
        N_{\beta,\gamma} e_{\beta + \gamma} & \beta + \gamma \in \Phi \\
        0 & \text{otherwise.} \end{cases}
        \end{aligned}

    where `A_{\alpha, \beta} = \frac{2 (\alpha, \beta)}{(\alpha, \alpha)}` and
    `N_{\alpha, \beta}` is the maximum such that
    `\alpha - N_{\alpha, \beta} \beta \in \Phi`.

    For computing the signs of the coefficients, see Section 3 of [CMT]_.

    REFERNCES:

    .. [CMT] A. M. Cohen, S. H. Murray, D. E. Talyor. *Groups of Lie type*.
       http://www.win.tue.nl/~amc/pub/papers/cmt.pdf
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type):
        """
        Normalize ``self`` to ensure a unique represntation.

        TESTS::

            sage: L1 = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L2 = LieAlgebra(QQ, cartan_type=CartanType(['A', 2]))
            sage: L3 = LieAlgebra(QQ, cartan_type=CartanMatrix(['A', 2]))
            sage: L1 is L2 and L2 is L3
            True
        """
        return super(LieAlgebraChevalleyBasis, cls).__classcall__(
            cls, R, CartanType(cartan_type))

    def __init__(self, R, cartan_type):
        r"""
        Initialize ``self``.

        TESTS::

            sage: L = LieAlgebra(QQ, cartan_type=['A',2])
            sage: TestSuite(L).run()
        """
        self._cartan_type = cartan_type
        RL = cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        p_roots = list(RL.positive_roots_by_height())
        n_roots = map(lambda x: -x, p_roots)
        alphacheck = RL.simple_coroots()
        roots = RL.roots()
        num_sroots = len(alpha)

        # Determine the signs for the structure coefficients from the root system
        # We first create the special roots
        sp_sign = {}
        for i,a in enumerate(p_roots):
            for b in p_roots[i+1:]:
                if a + b not in p_roots:
                    continue

                # Compute the sign for the extra special pair
                x, y = (a + b).extraspecial_pair()

                if (x, y) == (a, b): # If it already is an extra special pair
                    sp_sign[(x, y)] = 1
                    sp_sign[(y, x)] = -1
                    continue

                if b - x in roots:
                    t1 = (b-x).norm_squared() / b.norm_squared() * sp_sign[(x, b-x)] * sp_sign[(a, y-a)]
                else:
                    t1 = 0
                if a - x in roots:
                    t2 = (a-x).norm_squared() / a.norm_squared() * sp_sign[(x, a-x)] * sp_sign[(b, y-b)]
                else:
                    t2 = 0

                if t1 - t2 > 0:
                    sp_sign[(a,b)] = 1
                else:
                    sp_sign[(a,b)] = -1
                sp_sign[(b,a)] = -sp_sign[(a,b)]

        # Function to construct the structure coefficients (up to sign)
        def e_coeff(r, s):
            p = 1
            while r - p*s in roots:
                p += 1
            return p

        # Now we can compute all necessary structure coefficients
        s_coeffs = {}
        for i,r in enumerate(p_roots):
            # [e_r, h_i] and [h_i, f_r]
            for ac in alphacheck:
                c = r.scalar(ac)
                if c == 0:
                    continue
                s_coeffs[(r, ac)] = {r: -c}
                s_coeffs[(ac, -r)] = {-r: -c}

            # [e_r, f_r]
            h_sum = {}
            for j, c in r.associated_coroot():
                h_sum[alphacheck[j]] = c
            s_coeffs[(r, -r)] = h_sum

            # [e_r, e_s] and [e_r, f_s] with r != +/-s
            for j, s in enumerate(p_roots[i+1:]):
                j += i+1
                # Since h(s) > h(r), we will always have s - r > 0 (if it is a root)
                # [e_r, f_s]
                if s - r in p_roots:
                    c = e_coeff(r, -s)
                    a,b = s-r, r
                    if p_roots.index(a) + 1 > p_roots.index(b): # Note a != b
                        c = -c * sp_sign[(b, a)]
                    else:
                        c *= sp_sign[(a, b)]
                    s_coeffs[(-r, s)] = {a: -c}
                    s_coeffs[(r, -s)] = {a: c}

                # [e_r, e_s]
                a = r + s
                if a in p_roots:
                    # (r, s) is a special pair
                    c = e_coeff(r, s) * sp_sign[(r, s)]
                    s_coeffs[(r, s)] = {a: c}
                    s_coeffs[(-r, -s)] = {-a: -c}

        # Lastly, make sure a < b for all (a, b) in the coefficients and flip if necessary
        for k in s_coeffs.keys():
            a,b = k[0], k[1]
            if self._basis_cmp(a, b) > 0:
                s_coeffs[(b,a)] = {k:-v for k,v in s_coeffs[k].items()}
                del s_coeffs[k]

        names = ['e{}'.format(i) for i in range(1, num_sroots+1)]
        names += ['f{}'.format(i) for i in range(1, num_sroots+1)]
        names += ['h{}'.format(i) for i in range(1, num_sroots+1)]
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        index_set = p_roots + n_roots + list(alphacheck)
        LieAlgebraWithStructureCoefficients.__init__(self, R, s_coeffs, names, index_set,
                                                     category, prefix='E', bracket='[',
                                                     generator_cmp=self._basis_cmp)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['A', 2])
            Lie algebra of ['A', 2] in the Chevalley basis
        """
        return "Lie algebra of {} in the Chevalley basis".format(self._cartan_type)

    def _repr_generator(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.
        """
        if m in self._cartan_type.root_system().root_lattice().simple_coroots():
            return "h{}".format(m.support()[0])
        return IndexedGenerators._repr_generator(self, m)

    def _basis_cmp(self, x, y):
        """
        Compare two basis element indices. We order the basis elements by
        positive roots, coroots, and negative roots and then according to
        height.

        OUTPUT:

        If ``x == y``, return 0. If ``x < y``, return -1. Else return 1.
        """
        if x == y:
            return 0

        RL = self._cartan_type.root_system().root_lattice()
        p_roots = list(RL.positive_roots_by_height())
        n_roots = map(lambda x: -x, p_roots)
        alphacheck = RL.simple_coroots()

        if x in p_roots:
            if y in p_roots:
                return cmp(p_roots.index(x), p_roots.index(y))
            return -1

        if x in alphacheck:
            if y in p_roots:
                return 1
            if y in alphacheck:
                return cmp(x, y)
            return -1

        # x is in n_roots
        if y not in n_roots:
            return 1
        return cmp(n_roots.index(x), n_roots.index(y))

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.
        """
        return self._cartan_type

    def weyl_group(self):
        """
        Return the Weyl group of ``self``.
        """
        from sage.combinat.root_system.weyl_group import WeylGroup
        return WeylGroup(self._cartan_type)

    # TODO: Uncomment once #16825 is done
    #def affine(self, kac_moody=False):
    #    """
    #    Return the (untwisted) affine Lie algebra of ``self``.
    #    """
    #    from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra
    #    return AffineLieAlgebra(self, kac_moody)

    # Useful in creating the UEA
    @cached_method
    def indices_to_positive_roots_map(self):
        """
        Return the map from indices to positive roots.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.indices_to_positive_roots_map()
            {1: alpha[1], 2: alpha[2], 3: alpha[1] + alpha[2]}
        """
        RL = self._cartan_type.root_system().root_lattice()
        return {i+1: r for i,r in enumerate(RL.positive_roots())}

    def lie_algebra_generators(self, str_keys=False):
        r"""
        Return the Chevalley Lie algebra generators of ``self``.

        INPUT:

        - ``str_keys`` -- (default: ``False``) set to ``True`` to have the
          indices indexed by strings instead of simple (co)roots

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: L.lie_algebra_generators()
            Finite family {-alpha[1]: E[-alpha[1]], alpha[1]: E[alpha[1]], alphacheck[1]: h1}
            sage: L.lie_algebra_generators(True)
            Finite family {'f1': E[-alpha[1]], 'h1': h1, 'e1': E[alpha[1]]}
        """
        index_set = self._cartan_type.index_set()
        RL = self._cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        alphacheck = RL.simple_coroots()
        B = self.basis()
        ret = {}

        if str_keys:
            for i in index_set:
                al = alpha[i]
                ret['e{}'.format(i)] = B[al]
                ret['f{}'.format(i)] = B[-al]
                ret['h{}'.format(i)] = B[alphacheck[i]]
        else:
            for i in index_set:
                al = alpha[i]
                ret[al] = B[al]
                ret[-al] = B[-al]
                ret[alphacheck[i]] = B[alphacheck[i]]

        return Family(ret)

    @cached_method
    def gens(self):
        """
        Return the generators of ``self`` in the order of `e_i`, `f_i`,
        and `h_i`.
        """
        index_set = self._cartan_type.index_set()
        RL = self._cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        alphacheck = RL.simple_coroots()
        B = self.basis()

        ret = []
        for i in index_set:
            ret.append(B[alpha[i]])
        for i in index_set:
            ret.append(B[-alpha[i]])
        for i in index_set:
            ret.append(B[alphacheck[i]])
        return tuple(ret)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.
        """
        return self.gens()[i]

    def highest_root_basis_elt(self, pos=True):
        r"""
        Return the basis element corresponding to the highest root `\theta`.
        If ``pos`` is ``True``, then returns `e_{\theta}`, otherwise it
        returns `f_{\theta}`.
        """
        RL = self._cartan_type.root_system().root_lattice()
        theta = RL.highest_root()
        B = self.basis()
        if pos:
            return B[theta]
        return B[-theta]

