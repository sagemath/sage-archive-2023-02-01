"""
Quasisymmetric functions

REFERENCES:

.. [Ges] I. Gessel, Multipartite P-partitions and inner products of skew Schur
   functions, Contemp. Math. 34 (1984), 289-301.

.. [MR] C. Malvenuto and C. Reutenauer, Duality between quasi-symmetric
   functions and the Solomon descent algebra, J. Algebra 177 (1995), no. 3, 967-982.

AUTHOR:

- Jason Bandlow
- Franco Saliola
- Chris Berg
"""
#*****************************************************************************
#       Copyright (C) 2010 Jason Bandlow <jbandlow@gmail.com>,
#                     2012 Franco Saliola <saliola@gmail.com>,
#                     2012 Chris Berg <chrisjamesberg@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.bindable_class import BindableClass
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.all import CommutativeRings
from sage.categories.rings import Rings
from sage.categories.realizations import Category_realization_of_parent
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.matrix.constructor import matrix
from sage.combinat.subset import Subsets
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.composition import Composition, Compositions
from sage.combinat.composition_tableau import CompositionTableaux
from sage.combinat.partition import Partitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.ncsf_qsym.generic_basis_code import BasesOfQSymOrNCSF
from sage.combinat.ncsf_qsym.combinatorics import *
from sage.combinat.ncsf_qsym.ncsf import *
from sage.sets.set import Set

class QuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    .. rubric:: The Hopf algebra of quasisymmetric functions.

    The ring of quasi-symmetric functions may be realized as a subring
    of polynomials of `n` variables `{\bf k}[x_1, x_2, \ldots, x_n]` where the
    dimension of the subspace of elements of degree `k` is equal to
    the number of compositions of `k` (with less than `n` parts).  The two
    classical bases, the Monomial and Fundamental, are defined by the formulas:

    `M_I = \sum_{1 \leq i_1 < i_2 < \cdots < i_\ell \leq n} x_{i_1}^{I_1}
    x_{i_2}^{I_2} \cdots x_{i_\ell}^{I_\ell}`

    and

    `F_I = \sum_{1 \leq i_1 \leq i_2 \leq \cdots \leq i_\ell \leq n} x_{i_1}^{I_1}
    x_{i_2}^{I_2} \cdots x_{i_\ell}^{I_\ell}`

    where in the sum for the Fundamental basis there is strict inequality from `i_r` to
    `i_{r+1}` if `i_r` is a descent of the composition `I`.

    These bases are related by the formula

    `F_I = \sum_{J \leq I} M_J`

    where the inequality `J \leq I` indicates that `J` is finer than `I`.
    By taking the limit on `n`, the ring of quasi-symmetric polynomials is a Hopf algebra
    that inherits the product and co-product structure from the polynomial ring where it
    lives.

    .. rubric:: The implementation of the quasi-symmetric function Hopf algebra

    We realize the ring of quasi-symmetric functions in Sage as a graded Hopf
    algebra with basis elements indexed by compositions.
    ::

        sage: QSym = QuasiSymmetricFunctions(QQ)
        sage: QSym.category()
        Join of Category of graded hopf algebras over Rational Field and Category of monoids with realizations and Category of coalgebras over Rational Field with realizations

    Currently the quasi-symmetric functions have a minimal implementation and the only
    two bases implemented natively are the Monomial and Fundamental bases.
    ::

        sage: M = QSym.M()
        sage: F = QSym.F()
        sage: M(F[2,1,2])
        M[1, 1, 1, 1, 1] + M[1, 1, 1, 2] + M[2, 1, 1, 1] + M[2, 1, 2]
        sage: F(M[2,1,2])
        F[1, 1, 1, 1, 1] - F[1, 1, 1, 2] - F[2, 1, 1, 1] + F[2, 1, 2]

    The product on this space is commutative and is inherited from the product
    by the realization within the polynomial ring.
    ::

        sage: M[3]*M[1,1] == M[1,1]*M[3]
        True
        sage: M[3]*M[1,1]
        M[1, 1, 3] + M[1, 3, 1] + M[1, 4] + M[3, 1, 1] + M[4, 1]
        sage: F[3]*F[1,1]
        F[1, 1, 3] + F[1, 2, 2] + F[1, 3, 1] + F[1, 4] + F[2, 1, 2] + F[2, 2, 1] + F[2, 3] + F[3, 1, 1] + F[3, 2] + F[4, 1]
        sage: M[3]*F[2]
        M[1, 1, 3] + M[1, 3, 1] + M[1, 4] + M[2, 3] + M[3, 1, 1] + M[3, 2] + M[4, 1] + M[5]
        sage: F[2]*M[3]
        F[1, 1, 1, 2] - F[1, 2, 2] + F[2, 1, 1, 1] - F[2, 1, 2] - F[2, 2, 1] + F[5]

    There is a coproduct on this ring as well, which in the Monomial basis acts by
    cutting the composition into a left half and a right half.  The
    co-product is non-co-commutative.
    ::

        sage: M[1,3,1].coproduct()
        M[] # M[1, 3, 1] + M[1] # M[3, 1] + M[1, 3] # M[1] + M[1, 3, 1] # M[]
        sage: F[1,3,1].coproduct()
        F[] # F[1, 3, 1] + F[1] # F[3, 1] + F[1, 1] # F[2, 1] + F[1, 2] # F[1, 1] + F[1, 3] # F[1] + F[1, 3, 1] # F[]

    .. rubric:: The duality pairing with non-commutative symmetric functions

    These two operations endow the quasi-symmetric functions `QSym` with the
    structure of a Hopf algebra. It is the dual Hopf algebra of the
    non-commutative symmetric functions `NCSF`. Under this duality, the
    Monomial basis of `QSym` is dual to the Complete basis of `NCSF`, and the
    Fundamental basis of `QSym` is dual to the Ribbon basis of `NCSF` (see
    [MR]_).

    ::

        sage: S = M.dual(); S
        Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
        sage: M[1,3,1].duality_pairing( S[1,3,1] )
        1
        sage: M.duality_pairing_matrix( S, degree=4 )
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 1 0 0 0 0]
        [0 0 0 0 1 0 0 0]
        [0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 1]
        sage: F.duality_pairing_matrix( S, degree=4 )
        [1 0 0 0 0 0 0 0]
        [1 1 0 0 0 0 0 0]
        [1 0 1 0 0 0 0 0]
        [1 1 1 1 0 0 0 0]
        [1 0 0 0 1 0 0 0]
        [1 1 0 0 1 1 0 0]
        [1 0 1 0 1 0 1 0]
        [1 1 1 1 1 1 1 1]
        sage: NCSF = M.realization_of().dual()
        sage: R = NCSF.Ribbon()
        sage: F.duality_pairing_matrix( R, degree=4 )
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 1 0 0 0 0]
        [0 0 0 0 1 0 0 0]
        [0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 1]
        sage: M.duality_pairing_matrix( R, degree=4 )
        [ 1  0  0  0  0  0  0  0]
        [-1  1  0  0  0  0  0  0]
        [-1  0  1  0  0  0  0  0]
        [ 1 -1 -1  1  0  0  0  0]
        [-1  0  0  0  1  0  0  0]
        [ 1 -1  0  0 -1  1  0  0]
        [ 1  0 -1  0 -1  0  1  0]
        [-1  1  1 -1  1 -1 -1  1]

    Let `H` and `G` be elements of `QSym` and `h` an element of `NCSF`. Then if
    we represent the duality pairing with the mathematical notation `[ \cdot,
    \cdot ]`,

    `[H G, h] = [H \otimes G, \Delta(h)]~.`

    For example, the coefficient of ``M[2,1,4,1]`` in ``M[1,3]*M[2,1,1]`` may be
    computed with the duality pairing::

        sage: I, J = Composition([1,3]), Composition([2,1,1])
        sage: (M[I]*M[J]).duality_pairing(S[2,1,4,1])
        1

    And the coefficient of ``S[1,3] # S[2,1,1]`` in ``S[2,1,4,1].coproduct()`` is
    equal to this result.
    ::

        sage: S[2,1,4,1].coproduct()
        S[] # S[2, 1, 4, 1] + ... + S[1, 3] # S[2, 1, 1] + ... + S[4, 1] # S[2, 1]

    The duality pairing on the tensor space is another way of getting this
    coefficient, but currently the method ``duality_pairing`` is not defined on
    the tensor squared space. However, we can extend this functionality by
    applying a linear morphism to the terms in the coproduct, as follows.
    ::

        sage: X = S[2,1,4,1].coproduct()
        sage: def linear_morphism(x, y):
        ...     return x.duality_pairing(M[1,3]) * y.duality_pairing(M[2,1,1])
        sage: X.apply_multilinear_morphism(linear_morphism, codomain=ZZ)
        1

    Similarly, if `H` is an element of `QSym` and `g` and `h` are elements of `NCSF`,
    then

    `[ H, g h ] = [ \Delta(H), g \otimes h ]~.`

    For example, the coefficient of ``R[2,3,1]`` in ``R[2,1]*R[2,1]`` is computed with
    the duality pairing by the following command.
    ::

        sage: (R[2,1]*R[2,1]).duality_pairing(F[2,3,1])
        1
        sage: R[2,1]*R[2,1]
        R[2, 1, 2, 1] + R[2, 3, 1]

    This coefficient should then be equal to the coefficient of ``F[2,1] # F[2,1]``
    in ``F[2,3,1].coproduct()``.
    ::

        sage: F[2,3,1].coproduct()
        F[] # F[2, 3, 1] + ... + F[2, 1] # F[2, 1]  + ... + F[2, 3, 1] # F[]

    This can also be computed by the duality pairing on the tensor space,
    as above.
    ::

        sage: X = F[2,3,1].coproduct()
        sage: def linear_morphism(x, y):
        ...     return x.duality_pairing(R[2,1]) * y.duality_pairing(R[2,1])
        sage: X.apply_multilinear_morphism(linear_morphism, codomain=ZZ)
        1

    .. rubric:: The operation dual to multiplication by a non-commutative symmetric function

    Let `g \in NCSF` and consider the linear endomorphism of `NCSF` defined by
    left (respectively, right) multiplication by `g`. Since there is a duality
    between `QSym` and `NCSF`, this linear transformation induces an operator
    `g^\perp` on `QSym` satisfying

    `[ g^\perp(H), h ] = [ H, gh ]~.`

    for any non-commutative symmetric function `h`.

    This is implemented by the method :meth:`~sage.combinat.ncsf_qsym.generic_basis_code.BasesOfQSymOrNCSF.ElementMethods.skew_by`.
    Explicitly, if ``H`` is a quasi-symmetric function and ``g``
    a non-commutative symmetric function, then ``H.skew_by(g)`` and
    ``H.skew_by(g, side='right')`` are expressions that satisfy
    for any non-commutative symmetric function ``h``.

    ::

        H.skew_by(g).duality_pairing(h) == H.duality_pairing(g*h)
        H.skew_by(g, side='right').duality_pairing(h) == H.duality_pairing(h*g)

    For example, ``M[J].skew_by(S[I])`` is `0` unless the composition ``J``
    begins with ``I`` and ``M(J).skew_by(S(I), side='right')`` is `0` unless
    the composition ``J`` ends with ``I``.

    ::

        sage: M[3,2,2].skew_by(S[3])
        M[2, 2]
        sage: M[3,2,2].skew_by(S[2])
        0
        sage: M[3,2,2].coproduct().apply_multilinear_morphism( lambda x,y: x.duality_pairing(S[3])*y )
        M[2, 2]
        sage: M[3,2,2].skew_by(S[3], side='right')
        0
        sage: M[3,2,2].skew_by(S[2], side='right')
        M[3, 2]

    .. rubric:: The counit

    The counit is defined by sending all elements of positive degree to zero::

        sage: M[3].degree(), M[3,1,2].degree(), M.one().degree()
        (3, 6, 0)
        sage: M[3].counit()
        0
        sage: M[3,1,2].counit()
        0
        sage: M.one().counit()
        1
        sage: (M[3] - 2*M[3,1,2] + 7).counit()
        7
        sage: (F[3] - 2*F[3,1,2] + 7).counit()
        7

    .. rubric:: The antipode

    The antipode sends the Fundamental basis element indexed by the
    composition `I` to `-1` to the size of `I` times the Fundamental
    basis element indexed by the conjugate composition to `I`.
    ::

        sage: F[3,2,2].antipode()
        -F[1, 2, 2, 1, 1]
        sage: Composition([3,2,2]).conjugate()
        [1, 2, 2, 1, 1]
        sage: M[3,2,2].antipode()
        -M[2, 2, 3] - M[2, 5] - M[4, 3] - M[7]

    We demonstrate here the defining relation of the antipode::

        sage: X = F[3,2,2].coproduct()
        sage: X.apply_multilinear_morphism(lambda x,y: x*y.antipode())
        0
        sage: X.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
        0

    .. rubric:: The relation with symmetric functions

    The quasi-symmetric functions are a ring which contain the symmetric functions
    as a subring.  The Monomial quasi-symmetric functions are related to the
    monomial symmetric functions by

    `m_\lambda = \sum_{sort(I) = \lambda} M_I~.`

    There are methods to test if an expression in the quasi-symmetric functions is a symmetric
    function and, if it is, send it to an expression in the symmetric functions.
    ::

        sage: f = M[1,1,2] + M[1,2,1]
        sage: f.is_symmetric()
        False
        sage: g = M[3,1] + M[1,3]
        sage: g.is_symmetric()
        True
        sage: g.to_symmetric_function()
        m[3, 1]

    The expansion of the Schur function in terms of the Fundamental quasi-symmetric
    functions is due to [Ges]_.  There is one term in the expansion for each standard
    tableau of shape equal to the partition indexing the Schur function.
    ::

        sage: f = F[3,2] + F[2,2,1] + F[2,3] + F[1,3,1] + F[1,2,2]
        sage: f.is_symmetric()
        True
        sage: f.to_symmetric_function()
        5*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + 2*m[2, 2, 1] + m[3, 1, 1] + m[3, 2]
        sage: s = SymmetricFunctions(QQ).s()
        sage: s(f.to_symmetric_function())
        s[3, 2]

    It is also possible to convert a symmetric function to a quasi-symmetric function.
    ::

        sage: m = SymmetricFunctions(QQ).m()
        sage: M( m[3,1,1] )
        M[1, 1, 3] + M[1, 3, 1] + M[3, 1, 1]
        sage: F( s[2,2,1] )
        F[1, 1, 2, 1] + F[1, 2, 1, 1] + F[1, 2, 2] + F[2, 1, 2] + F[2, 2, 1]

    It is possible to experiment with the quasi-symmetric function expansion of other
    bases, but it is important that the base ring be the same for both algebras.
    ::

        sage: R = QQ['t']
        sage: Qp = SymmetricFunctions(R).hall_littlewood().Qp()
        sage: QSymt = QuasiSymmetricFunctions(R)
        sage: Ft = QSymt.F()
        sage: Ft( Qp[2,2] )
        F[1, 2, 1] + t*F[1, 3] + (t+1)*F[2, 2] + t*F[3, 1] + t^2*F[4]

    ::

        sage: K = QQ['q','t'].fraction_field()
        sage: Ht = SymmetricFunctions(K).macdonald().Ht()
        sage: Fqt = QuasiSymmetricFunctions(Ht.base_ring()).F()
        sage: Fqt(Ht[2,1])
        q*t*F[1, 1, 1] + (q+t)*F[1, 2] + (q+t)*F[2, 1] + F[3]

    The following will raise an error because the base ring of ``F`` is not
    equal to the base ring of ``Ht``.
    ::

        sage: F(Ht[2,1])
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= McdHt[2, 1]) an element of self (=Quasisymmetric functions over the Rational Field in the Fundamental basis)

    .. rubric:: The map to the ring of polynomials

    Since the quasi-symmetric functions are the limit of a subring of polynomials
    as the number of variables increases, there exists a projection
    from the quasi-symmetric functions into the polynomial ring `{\bf k}[x_1, x_2, \ldots, x_n]`.
    Although not precise, we may think of the quasi-symmetric functions indexed by a
    composition as a function in an infinite number of variables and project this into the
    ring of polynomials by setting `x_{n+1} = x_{n+2} = \cdots = 0`.  If the the number
    of variables is smaller than the length of the composition, the result is `0`.
    ::

        sage: M[1,3,1].expand(4)
        x0*x1^3*x2 + x0*x1^3*x3 + x0*x2^3*x3 + x1*x2^3*x3
        sage: F[1,3,1].expand(4)
        x0*x1^3*x2 + x0*x1^3*x3 + x0*x1^2*x2*x3 + x0*x1*x2^2*x3 + x0*x2^3*x3 + x1*x2^3*x3
        sage: M[1,3,1].expand(2)
        0

    TESTS::

        sage: QSym = QuasiSymmetricFunctions(QQ); QSym
        Quasisymmetric functions over the Rational Field
        sage: QSym.base_ring()
        Rational Field
    """
    def __init__(self, R):
        """
        The Hopf algebra of quasi-symmetric functions.
        See ``QuasiSymmetricFunctions`` for full documentation.

        EXAMPLES::

            sage: QuasiSymmetricFunctions(QQ)
            Quasisymmetric functions over the Rational Field
            sage: TestSuite(QuasiSymmetricFunctions(QQ)).run()

        """
        assert R in Rings()
        self._base = R # Won't be needed once CategoryObject won't override base_ring
        category = GradedHopfAlgebras(R)  # TODO: .Commutative()
        Parent.__init__(self, category = category.WithRealizations())

        # Bases
        Monomial    = self.Monomial()
        Fundamental = self.Fundamental()
        dualImmaculate = self.dualImmaculate()
        QS          = self.Quasisymmetric_Schur()

        # Change of bases
        Fundamental.module_morphism(Monomial.sum_of_finer_compositions,
                                    codomain=Monomial, category=category
                                    ).register_as_coercion()
        Monomial   .module_morphism(Fundamental.alternating_sum_of_finer_compositions,
                                    codomain=Fundamental, category=category
                                    ).register_as_coercion()
        #This changes dualImmaculate into Monomial
        dualImmaculate.module_morphism(dualImmaculate._to_Monomial_on_basis,
                                          codomain = Monomial, category = category
                                          ).register_as_coercion()
        #This changes Monomial into dualImmaculate
        Monomial.module_morphism(dualImmaculate._from_Monomial_on_basis,
                                          codomain = dualImmaculate, category = category
                                          ).register_as_coercion()
        #This changes Quasisymmetric Schur into Fundamental
        QS         .module_morphism(QS._to_fundamental_on_basis,
                                    codomain=Fundamental, category=category
                                    ).register_as_coercion()
        #This changes Fundamental into Quasisymmetric Schur
        Fundamental.module_morphism(QS._from_fundamental_on_basis,
                                    codomain=QS, category=category
                                    ).register_as_coercion()

        # Embedding of Sym into QSym in the monomial bases
        Sym = SymmetricFunctions(self.base_ring())
        Sym_m_to_M = Sym.m().module_morphism(Monomial.sum_of_partition_rearrangements,
                                           triangular='upper', inverse_on_support=Monomial._comp_to_par,
                                           codomain=Monomial, category=category)
        Sym_m_to_M.register_as_coercion()
        self.to_symmetric_function = Sym_m_to_M.section()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: M = QuasiSymmetricFunctions(ZZ).M()
            sage: M._repr_()
            'Quasisymmetric functions over the Integer Ring in the Monomial basis'
        """
        return "Quasisymmetric functions over the %s"%self.base_ring()

    def a_realization(self):
        r"""
        Returns the realization of the Monomial basis of the ring of quasi-symmetric functions.

        OUTPUT:

        - The Monomial basis of quasi-symmetric functions.

        EXAMPLES::

            sage: QuasiSymmetricFunctions(QQ).a_realization()
            Quasisymmetric functions over the Rational Field in the Monomial basis
        """
        return self.Monomial()

    _shorthands = tuple(['M', 'F', 'dI', 'QS'])

    def dual(self):
        r"""
        Returns the dual Hopf algebra of the quasi-symmetric functions, which is the
        non-commutative symmetric functions.

        OUTPUT:

        - The non-commutative symmetric functions.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QSym.dual()
            Non-Commutative Symmetric Functions over the Rational Field
        """
        from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
        return NonCommutativeSymmetricFunctions(self.base_ring())

    def from_polynomial(self, f, check=True):
        """
        Returns the quasi-symmetric function in the Monomial basis
        corresponding to the quasi-symmetric polynomial ``f``.

        INPUT:

        - ``f`` -- a polynomial in finitely many variables over the same base
          ring as ``self``. It is assumed that this polynomial is
          quasi-symmetric.
        - ``check`` -- boolean (default: ``True``), checks whether the
          polynomial is indeed quasi-symmetric.

        OUTPUT:

        - quasi-symmetric function in the Monomial basis

        EXAMPLES::

            sage: P = PolynomialRing(QQ, 'x', 3)
            sage: x = P.gens()
            sage: f = x[0] + x[1] + x[2]
            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QSym.from_polynomial(f)
            M[1]

        Beware of setting ``check=False``::

            sage: f = x[0] + 2*x[1] + x[2]
            sage: QSym.from_polynomial(f, check=True)
            Traceback (most recent call last):
            ...
            ValueError: x0 + 2*x1 + x2 is not a quasi-symmetric polynomial
            sage: QSym.from_polynomial(f, check=False)
            M[1]

        To expand the quasi-symmetric function in a basis other than the
        Monomial basis, the following shorthands are provided::

            sage: M = QSym.Monomial()
            sage: f = x[0]**2+x[1]**2+x[2]**2
            sage: g = M.from_polynomial(f); g
            M[2]
            sage: F = QSym.Fundamental()
            sage: F(g)
            -F[1, 1] + F[2]
            sage: F.from_polynomial(f)
            -F[1, 1] + F[2]

        """
        assert self.base_ring() == f.base_ring()
        exponent_coefficient = f.dict()
        z = {}
        for (e, c) in exponent_coefficient.iteritems():
            I = Composition([ei for ei in e if ei > 0])
            if I not in z:
                z[I] = c
        out = self.Monomial()._from_dict(z)
        if check and out.expand(f.parent().ngens(), f.parent().gens()) != f:
            raise ValueError("%s is not a quasi-symmetric polynomial" % f)
        return out

    class Bases(Category_realization_of_parent):
        r"""
        Category of bases of quasi-symmetric functions.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QSym.Bases()
            Category of bases of Quasisymmetric functions over the Rational Field
        """
        def super_categories(self):
            r"""
            Returns the super categories of bases of the Quasi-symmetric functions.

            OUTPUT:

            - a list of categories

            TESTS::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QSym.Bases().super_categories()
                [Category of bases of Non-Commutative Symmetric Functions or Quasisymmetric functions over the Rational Field, Category of commutative rings]

            """
            return [BasesOfQSymOrNCSF(self.base()), CommutativeRings()]

        class ParentMethods:
            def from_polynomial(self, f, check=True):
                """
                The quasi-symmetric function expanded in this basis
                corresponding to the quasi-symmetric polynomial ``f``.

                This is a default implementation that computes
                the expansion in the Monomial basis and converts
                to this basis.

                INPUT:

                - ``f`` -- a polynomial in finitely many variables over the same base
                  ring as ``self``. It is assumed that this polynomial is
                  quasi-symmetric.
                - ``check`` -- boolean (default: ``True``), checks whether the
                  polynomial is indeed quasi-symmetric.

                OUTPUT:

                - quasi-symmetric function

                EXAMPLES::

                    sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                    sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                    sage: P = PolynomialRing(QQ, 'x', 3)
                    sage: x = P.gens()
                    sage: f = x[0] + x[1] + x[2]
                    sage: M.from_polynomial(f)
                    M[1]
                    sage: F.from_polynomial(f)
                    F[1]
                    sage: f = x[0]**2+x[1]**2+x[2]**2
                    sage: M.from_polynomial(f)
                    M[2]
                    sage: F.from_polynomial(f)
                    -F[1, 1] + F[2]

                If the polynomial is not quasi-symmetric, an error
                is raised::

                    sage: f = x[0]^2+x[1]
                    sage: M.from_polynomial(f)
                    Traceback (most recent call last):
                    ...
                    ValueError: x0^2 + x1 is not a quasi-symmetric polynomial
                    sage: F.from_polynomial(f)
                    Traceback (most recent call last):
                    ...
                    ValueError: x0^2 + x1 is not a quasi-symmetric polynomial


                TESTS:

                We convert some quasi-symmetric functions to quasi-symmetric
                polynomials and back::

                    sage: f = (M[1,2] + M[1,1]).expand(3); f
                    x0*x1^2 + x0*x2^2 + x1*x2^2 + x0*x1 + x0*x2 + x1*x2
                    sage: M.from_polynomial(f)
                    M[1, 1] + M[1, 2]
                    sage: f = (2*M[2,1]+M[1,1]+3*M[3]).expand(3)
                    sage: M.from_polynomial(f)
                    M[1, 1] + 2*M[2, 1] + 3*M[3]
                    sage: f = (F[1,2] + F[1,1]).expand(3); f
                    x0*x1^2 + x0*x1*x2 + x0*x2^2 + x1*x2^2 + x0*x1 + x0*x2 + x1*x2
                    sage: F.from_polynomial(f)
                    F[1, 1] + F[1, 2]
                    sage: f = (2*F[2,1]+F[1,1]+3*F[3]).expand(3)
                    sage: F.from_polynomial(f)
                    F[1, 1] + 2*F[2, 1] + 3*F[3]
                """
                g = self.realization_of().from_polynomial(f, check=check)
                return self(g)

        class ElementMethods:
            def expand(self,n,alphabet='x'):
                r"""
                Expands the quasi-symmetric function into ``n`` variables in an alphabet,
                which by default is 'x'.

                INPUT:

                - ``n`` -- An integer `>0`; the number of variables in the expansion.
                - ``alphabet`` -- (default:'x'); the alphabet in which ``self`` is to be
                  expanded.

                OUTPUT:

                - An expansion of ``self`` into the ``n`` variables specified by ``alphabet``.

                EXAMPLES::

                    sage: F=QuasiSymmetricFunctions(QQ).Fundamental()
                    sage: F[3].expand(3)
                    x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3
                    sage: F[2,1].expand(3)
                    x0^2*x1 + x0^2*x2 + x0*x1*x2 + x1^2*x2

                One can use a different set of variable by adding an optional
                argument alphabet=...
                ::

                    sage: F=QuasiSymmetricFunctions(QQ).Fundamental()
                    sage: F[3].expand(2,alphabet='y')
                    y0^3 + y0^2*y1 + y0*y1^2 + y1^3

                TESTS::

                    sage: (3*F([])).expand(2)
                    3
                    sage: F[4,2].expand(0)
                    0
                    sage: F([]).expand(0)
                    1
                """
                M = self.parent().realization_of().Monomial()
                return M(self).expand(n,alphabet)

            def is_symmetric( self ):
                r"""
                Returns True if ``self`` is an element of the symmetric
                functions.  It tests this by looking at the expansion in
                the Monomial basis and testing if the coefficients are
                the same if the indexing compositions are permutations
                of each other.

                OUTPUT:

                - True if ``self`` is symmetric. False if ``self`` is not symmetric.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: F = QSym.Fundamental()
                    sage: (F[3,2] + F[2,3]).is_symmetric()
                    False
                    sage: (F[1, 1, 1, 2] + F[1, 1, 3] + F[1, 3, 1] + F[2, 1, 1, 1] + F[3, 1, 1]).is_symmetric()
                    True
                    sage: F([]).is_symmetric()
                    True
                """
                M = self.parent().realization_of().Monomial()
                return M(self).is_symmetric()

            def to_symmetric_function(self):
                r"""
                Converts a quasi-symmetric function to a symmetric function.

                OUTPUT:

                - If ``self`` is a symmetric function, then return the expansion
                  in the monomial basis.  Otherwise raise an error.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: F = QSym.Fundamental()
                    sage: (F[3,2] + F[2,3]).to_symmetric_function()
                    Traceback (most recent call last):
                    ...
                    ValueError: F[2, 3] + F[3, 2] is not a symmetric function
                    sage: m = SymmetricFunctions(QQ).m()
                    sage: s = SymmetricFunctions(QQ).s()
                    sage: F(s[3,1,1]).to_symmetric_function()
                    6*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + m[2, 2, 1] + m[3, 1, 1]
                    sage: m(s[3,1,1])
                    6*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + m[2, 2, 1] + m[3, 1, 1]

                """
                if self.is_symmetric():
                    M = self.parent().realization_of().Monomial()
                    return M( self ).to_symmetric_function()
                else:
                    raise ValueError, "%s is not a symmetric function"%self

    class Monomial(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric function in the Monomial basis.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: M = QSym.M()
            sage: F = QSym.F()
            sage: M(F[2,2])
            M[1, 1, 1, 1] + M[1, 1, 2] + M[2, 1, 1] + M[2, 2]
            sage: m = SymmetricFunctions(QQ).m()
            sage: M(m[3,1,1])
            M[1, 1, 3] + M[1, 3, 1] + M[3, 1, 1]
            sage: (1+M[1])^3
            M[] + 3*M[1] + 6*M[1, 1] + 6*M[1, 1, 1] + 3*M[1, 2] + 3*M[2] + 3*M[2, 1] + M[3]
            sage: M[1,2,1].coproduct()
            M[] # M[1, 2, 1] + M[1] # M[2, 1] + M[1, 2] # M[1] + M[1, 2, 1] # M[]

        The following is an alias for this basis::

            sage: QSym.Monomial()
            Quasisymmetric functions over the Rational Field in the Monomial basis

        TESTS::

            sage: M(F([]))
            M[]
            sage: M(F(0))
            0
            sage: M(m([]))
            M[]
        """
        def __init__(self, QSym):
            """
            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).Monomial(); M
                Quasisymmetric functions over the Rational Field in the Monomial basis
                sage: TestSuite(M).run()
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='M', bracket=False,
                                             category=QSym.Bases())

        def dual(self):
            r"""
            Returns the dual basis to the Monomial basis. This is the complete basis of the
            non-commutative symmetric functions.

            OUTPUT:

            - The complete basis of the non-commutative symmetric functions.

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M.dual()
                Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
            """
            return self.realization_of().dual().Complete()

        def product_on_basis(self, I, J):
            """
            The product on Monomial basis elements. The product of the basis elements
            indexed by two compositions `I` and `J` is the sum of the basis elements
            indexed by compositions in the shuffle product of `I` and `J`.

            INPUT:

            - ``I``, ``J`` -- compositions

            OUTPUT:

            - The product of the Monomial quasi-symmetric functions indexed by ``I`` and
              ``J``, expressed in the Monomial basis.

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: c1 = Composition([2])
                sage: c2 = Composition([1,3])
                sage: M.product_on_basis(c1, c2)
                M[1, 2, 3] + M[1, 3, 2] + M[1, 5] + M[2, 1, 3] + M[3, 3]
                sage: M.product_on_basis(c1, Composition([]))
                M[2]
            """
            return self.sum_of_monomials(I.shuffle_product(J, overlap=True))

        def antipode_on_basis(self,compo):
            r"""
            Returns the result of the antipode applied to a quasi-symmetric Monomial basis
            element.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The result of the antipode applied to the composition ``compo``, expressed
              in the Monomial basis.

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M.antipode_on_basis(Composition([2,1]))
                M[1, 2] + M[3]
                sage: M.antipode_on_basis(Composition([]))
                M[]
            """
            return (-1)**(len(compo))*self.sum_of_fatter_compositions(compo.reversed())

        def coproduct_on_basis(self, compo):
            r"""
            Returns the coproduct of a Monomial basis element.

            Combinatorial rule: deconcatenation.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The coproduct applied to the Monomial quasi-symmetric function indexed by
              ``compo``, expressed in the Monomial basis.

            EXAMPLES::

                sage: M=QuasiSymmetricFunctions(QQ).Monomial()
                sage: M[4,2,3].coproduct()
                M[] # M[4, 2, 3] + M[4] # M[2, 3] + M[4, 2] # M[3] + M[4, 2, 3] # M[]
                sage: M.coproduct_on_basis(Composition([]))
                M[] # M[]

            """
            return self.tensor_square().sum_of_monomials((Composition(compo[:i]), Composition(compo[i:]))
                                                         for i in range(0,len(compo)+1))

        class Element(CombinatorialFreeModule.Element):

            def expand(self, n, alphabet='x'):
                r"""
                Expands the quasi-symmetric function written in the monomial basis in
                `n` variables.

                INPUT:

                - ``n`` -- an integer
                - ``alphabet`` -- (default:'x') a string

                OUTPUT:

                - The quasi-symmetric function ``self`` expressed in the ``n`` variables
                  described by ``alphabet``.

                .. TODO:: accept an *alphabet* as input

                EXAMPLES::

                    sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                    sage: M[4,2].expand(3)
                    x0^4*x1^2 + x0^4*x2^2 + x1^4*x2^2

                One can use a different set of variable by using the
                optional argument ``alphabet``::

                    sage: M=QuasiSymmetricFunctions(QQ).Monomial()
                    sage: M[2,1,1].expand(4,alphabet='y')
                    y0^2*y1*y2 + y0^2*y1*y3 + y0^2*y2*y3 + y1^2*y2*y3

                TESTS::

                    sage: (3*M([])).expand(2)
                    3
                    sage: M[4,2].expand(0)
                    0
                    sage: M([]).expand(0)
                    1
                """
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                M = self.parent()
                P = PolynomialRing(M.base_ring(), n, alphabet)
                x = P.gens()
                def on_basis(comp, i):
                    if not comp:
                        return P.one()
                    elif len(comp) > i:
                        return P.zero()
                    else:
                        return x[i-1]**comp[-1] * on_basis(comp[:-1], i-1) + \
                                                  on_basis(comp,      i-1)
                return M._apply_module_morphism(self, lambda comp: on_basis(comp,n),
                                                codomain = P)

            def is_symmetric( self ):
                r"""
                Determines if a quasi-symmetric function, written in the Monomial basis,
                is symmetric. It tests this by looking at the expansion in the Monomial
                basis and testing if the coefficients are the same if the indexing
                compositions are permutations of each other.

                OUTPUT:

                - Returns True if ``self`` is an element of the symmetric functions and
                  False otherwise.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.Monomial()
                    sage: (M[3,2] + M[2,3] + M[4,1]).is_symmetric()
                    False
                    sage: (M[3,2] + M[2,3]).is_symmetric()
                    True
                    sage: (M[1,2,1] + M[1,1,2]).is_symmetric()
                    False
                    sage: (M[1,2,1] + M[1,1,2] + M[2,1,1]).is_symmetric()
                    True

                """
                # We must check that every rearrangement of a composition
                # that appears in self appears with the same coefficient.
                # We use a dictionary to keep track of the coefficient
                # and how many rearrangements of the composition we've seen.
                d = {}
                for (I, coeff) in self:
                    partition = I.to_partition()
                    if partition not in d:
                        d[partition] = [coeff, 1]
                    else:
                        if d[partition][0] != coeff:
                            return False
                        else:
                            d[partition][1] += 1
                # make sure we've seen each rearrangement of the composition
                return all(d[partition][1] == Permutations(partition).cardinality()
                            for partition in d)

            def to_symmetric_function( self ):
                r"""
                Takes a quasi-symmetric function, expressed in the monomial basis, and
                returns its symmetric realization, when possible, expressed in the
                monomial basis of symmetric functions.

                OUTPUT:

                - If ``self`` is a symmetric function, then the expansion
                  in the monomial basis of the symmetric functions is returned.
                  Otherwise an error is raised.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.Monomial()
                    sage: (M[3,2] + M[2,3] + M[4,1]).to_symmetric_function()
                    Traceback (most recent call last):
                    ...
                    ValueError: M[2, 3] + M[3, 2] + M[4, 1] is not a symmetric function
                    sage: (M[3,2] + M[2,3] + 2*M[4,1] + 2*M[1,4]).to_symmetric_function()
                    m[3, 2] + 2*m[4, 1]
                    sage: m = SymmetricFunctions(QQ).m()
                    sage: M(m[3,1,1]).to_symmetric_function()
                    m[3, 1, 1]
                    sage: (M(m[2,1])*M(m[2,1])).to_symmetric_function()-m[2,1]*m[2,1]
                    0

                TESTS::

                    sage: (M(0)).to_symmetric_function()
                    0
                    sage: (M([])).to_symmetric_function()
                    m[]
                    sage: (2*M([])).to_symmetric_function()
                    2*m[]

                """
                m = SymmetricFunctions(self.parent().base_ring()).monomial()
                if self.is_symmetric():
                    return m.sum_of_terms([(I, coeff) for (I, coeff) in self
                        if list(I) in Partitions()], distinct=True)
                else:
                    raise ValueError, "%s is not a symmetric function"%self

    M = Monomial

    class Fundamental(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric function in the Fundamental basis.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: F = QSym.F()
            sage: M = QSym.M()
            sage: F(M[2,2])
            F[1, 1, 1, 1] - F[1, 1, 2] - F[2, 1, 1] + F[2, 2]
            sage: s = SymmetricFunctions(QQ).s()
            sage: F(s[3,2])
            F[1, 2, 2] + F[1, 3, 1] + F[2, 2, 1] + F[2, 3] + F[3, 2]
            sage: (1+F[1])^3
            F[] + 3*F[1] + 3*F[1, 1] + F[1, 1, 1] + 2*F[1, 2] + 3*F[2] + 2*F[2, 1] + F[3]
            sage: F[1,2,1].coproduct()
            F[] # F[1, 2, 1] + F[1] # F[2, 1] + F[1, 1] # F[1, 1] + F[1, 2] # F[1] + F[1, 2, 1] # F[]

        The following is an alias for this basis::

            sage: QSym.Fundamental()
            Quasisymmetric functions over the Rational Field in the Fundamental basis

        TESTS::

            sage: F(M([]))
            F[]
            sage: F(M(0))
            0
            sage: F(s([]))
            F[]
            sage: F(s(0))
            0
        """
        def __init__(self, QSym):
            """
            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).Fundamental(); F
                Quasisymmetric functions over the Rational Field in the Fundamental basis
                sage: TestSuite(F).run()
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='F', bracket=False,
                                             category=QSym.Bases())

        def dual(self):
            r"""
            Returns the dual basis to the Fundamental basis. This is the ribbon
            basis of the non-commutative symmetric functions.

            OUTPUT:

            - The ribbon basis of the non-commutative symmetric functions.

            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F.dual()
                Non-Commutative Symmetric Functions over the Rational Field in the Ribbon basis
            """
            return self.realization_of().dual().Ribbon()

        def antipode_on_basis(self, compo):
            r"""
            Returns the antipode to a Fundamental quasi-symmetric basis element.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - Returns the result of the antipode applied to the quasi-symmetric
              Fundamental basis element indexed by ``compo``.

            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F.antipode_on_basis(Composition([2,1]))
                -F[2, 1]
            """
            return (-1)**(compo.size()) * self.monomial(compo.conjugate())

        def coproduct_on_basis(self, compo):
            r"""
            Returns the coproduct to a Fundamental quasi-symmetric basis element.

            Combinatorial rule: quasi deconcatenation.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - Returns the application of the coproduct to the Fundamental quasi-symmetric
              function indexed by the composition ``compo``.

            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: F[4].coproduct()
                F[] # F[4] + F[1] # F[3] + F[2] # F[2] + F[3] # F[1] + F[4] # F[]
                sage: F[2,1,3].coproduct()
                F[] # F[2, 1, 3] + F[1] # F[1, 1, 3] + F[2] # F[1, 3] + F[2, 1] # F[3] + F[2, 1, 1] # F[2] + F[2, 1, 2] # F[1] + F[2, 1, 3] # F[]

            TESTS::

                sage: F.coproduct_on_basis(Composition([2,1,3]))
                F[] # F[2, 1, 3] + F[1] # F[1, 1, 3] + F[2] # F[1, 3] + F[2, 1] # F[3] + F[2, 1, 1] # F[2] + F[2, 1, 2] # F[1] + F[2, 1, 3] # F[]
                sage: F.one().coproduct()          # generic for graded / graded connected
                F[] # F[]
            """
            T = self.tensor_square()
            C = Composition
            return T.sum_of_monomials( ( C(compo[:i]), C(compo[i:]) ) for i in range(len(compo)+1) ) + \
                   T.sum_of_monomials( ( C(compo[:i]+[j]), C([compo[i]-j]+compo[i+1:]) )
                                       for i in range(len(compo))
                                       for j in range(1, compo[i]) )

    F = Fundamental

    class Quasisymmetric_Schur(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric function in the Quasisymmetric
        Schur basis.

        The basis of Quasisymmetric Schur functions is defined in [QSCHUR]_.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QS = QSym.QS()
            sage: F = QSym.F()
            sage: M = QSym.M()
            sage: F(QS[1,2])
            F[1, 2]
            sage: M(QS[1,2])
            M[1, 1, 1] + M[1, 2]
            sage: s = SymmetricFunctions(QQ).s()
            sage: QS(s[2,1,1])
            QS[1, 1, 2] + QS[1, 2, 1] + QS[2, 1, 1]
        """
        def __init__(self, QSym):
            r"""
            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: F = QSym.Fundamental()
                sage: QS = QSym.Quasisymmetric_Schur()
                sage: QS(F(QS.an_element())) == QS.an_element()
                True
                sage: F(QS(F.an_element())) == F.an_element()
                True
                sage: M = QSym.Monomial()
                sage: QS(M(QS.an_element())) == QS.an_element()
                True
                sage: M(QS(M.an_element())) == M.an_element()
                True
                sage: TestSuite(QS).run()
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='QS', bracket=False,
                                             category=QSym.Bases())

        def _realization_name(self):
            r"""
            Return a nicer name for ``self`` than what is inherited
            from :mod:`sage.categories.sets_cat`.

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._realization_name()
                'Quasisymmetric Schur'
            """
            return "Quasisymmetric Schur"

        def _to_fundamental_on_basis(self, comp):
            r"""
            Map the quasi-symmetric Schur function indexed by ``comp`` to
            the Fundamental basis.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT:

            - a quasi-symmetric function in the Fundamental basis

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._to_fundamental_on_basis([1,3,1])
                F[1, 3, 1] + F[2, 2, 1]
            """
            F = self.realization_of().Fundamental()
            return F.sum_of_monomials(T.descent_composition() for T in CompositionTableaux(comp) if T.is_standard())

        ##########################################################################
        # Implementation of the from_fundamental by inverting to_fundamental
        # TODO: discard once we have inverse isomorphisms for graded vector spaces
        # and throw instead, in QuasiSymmetricFunctions.__init__, something like:
        #  Schur_to_F.inverse().register_as_coercion()

        @cached_method
        def _to_fundamental_transition_matrix(self, n):
            r"""
            Return the transition matrix from the basis of Quasisymmetric
            Schur functions to the Fundamental basis.

            INPUT:

            - ``n`` -- an integer

            OUTPUT:

            - a square matrix

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._to_fundamental_transition_matrix(4)
                [1 0 0 0 0 0 0 0]
                [0 1 0 0 0 0 0 0]
                [0 0 1 1 0 0 0 0]
                [0 0 0 1 0 1 0 0]
                [0 0 0 0 1 0 0 0]
                [0 0 0 0 0 1 0 0]
                [0 0 0 0 0 0 1 0]
                [0 0 0 0 0 0 0 1]
            """
            if n == 0:
                return matrix([[]])

            ranks = dict((comp,rank) for (rank,comp) in enumerate(compositions_order(n)))
            d = {}
            for T in CompositionTableaux(n):
                if T.is_standard():
                    I = T.shape_composition()
                    J = T.descent_composition()
                    if (I,J) in d:
                        d[I,J] += 1
                    else:
                        d[I,J] = 1
            m = {}
            for (I,J) in d:
                m[ranks[I], ranks[J]] = d[I,J]
            return matrix(len(ranks), len(ranks), m)

        @cached_method
        def _from_fundamental_transition_matrix(self, n):
            r"""
            Return the transition matrix from the Fundamental basis of
            quasi-symmetric functions to the Quasisymmetric Schur basis.

            INPUT:

            - ``n`` -- an integer

            OUTPUT:

            - a square matrix

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._from_fundamental_transition_matrix(4)
                [ 1  0  0  0  0  0  0  0]
                [ 0  1  0  0  0  0  0  0]
                [ 0  0  1 -1  0  1  0  0]
                [ 0  0  0  1  0 -1  0  0]
                [ 0  0  0  0  1  0  0  0]
                [ 0  0  0  0  0  1  0  0]
                [ 0  0  0  0  0  0  1  0]
                [ 0  0  0  0  0  0  0  1]
            """
            if n == 0:
                return matrix([[]])
            return self._to_fundamental_transition_matrix(n).inverse()

        def _from_fundamental_on_basis(self, comp):
            r"""
            Maps the Fundamental quasi-symmetric function indexed by ``comp`` to the Quasisymmetric
            Schur basis.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT:

            - a quasi-symmetric function in the Quasisymmetric Schur basis

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._from_fundamental_on_basis([1,3,1])
                QS[1, 2, 1, 1] + QS[1, 3, 1] - QS[2, 2, 1]
            """
            comp = Composition(comp)
            if comp == []:
                return self.one()
            comps = compositions_order(comp.size())
            T = self._from_fundamental_transition_matrix(comp.size())
            return self.sum_of_terms( zip(comps, T[comps.index(comp)]) )

    QS = Quasisymmetric_Schur

    class dualImmaculate(CombinatorialFreeModule, BindableClass):
        def __init__(self, QSym):
            r"""
            The dual immaculate basis of the non-commutative symmetric functions. This basis first
            appears in Berg, Bergeron, Saliola, Serrano and Zabrocki's " A lift of the Schur and Hall-Littlewood
            bases to non-commutative symmetric functions".

            EXAMPLES ::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: dI = QSym.dI()
                sage: dI([1,3,2])*dI([1])  # long time (6s on sage.math, 2013)
                dI[1, 1, 3, 2] + dI[2, 3, 2]
                sage: dI([1,3])*dI([1,1])
                dI[1, 1, 1, 3] + dI[1, 1, 4] + dI[1, 2, 3] - dI[1, 3, 2] - dI[1, 4, 1] - dI[1, 5] + dI[2, 3, 1] + dI[2, 4]
                sage: dI([3,1])*dI([2,1])  # long time (7s on sage.math, 2013)
                dI[1, 1, 5] - dI[1, 4, 1, 1] - dI[1, 4, 2] - 2*dI[1, 5, 1] - dI[1, 6] - dI[2, 4, 1] - dI[2, 5] - dI[3, 1, 3] + dI[3, 2, 1, 1] + dI[3, 2, 2] + dI[3, 3, 1] + dI[4, 1, 1, 1] + 2*dI[4, 2, 1] + dI[4, 3] + dI[5, 1, 1] + dI[5, 2]
                sage: F = QSym.F()
                sage: dI(F[1,3,1])
                -dI[1, 1, 1, 2] + dI[1, 1, 2, 1] - dI[1, 2, 2] + dI[1, 3, 1]
                sage: F(dI(F([2,1,3])))
                F[2, 1, 3]
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='dI', bracket=False,
                                             category=QSym.Bases())

        def _to_Monomial_on_basis(self, J):
            r"""
            Given a dual immaculate function, this method returns the expansion of the function in the quasi-symmetric monomial basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT:

            - A quasi-symmetric function in the monomial basis.

            EXAMPLES::

            sage: dI = QuasiSymmetricFunctions(QQ).dI()
            sage: dI._to_Monomial_on_basis(Composition([1,3]))
            M[1, 1, 1, 1] + M[1, 1, 2] + M[1, 2, 1] + M[1, 3]
            sage: dI._to_Monomial_on_basis(Composition([]))
            M[]
            sage: dI._to_Monomial_on_basis(Composition([2,1,2]))
            4*M[1, 1, 1, 1, 1] + 3*M[1, 1, 1, 2] + 2*M[1, 1, 2, 1] + M[1, 1, 3] + M[1, 2, 1, 1] + M[1, 2, 2] + M[2, 1, 1, 1] + M[2, 1, 2]
            """
            M = self.realization_of().Monomial()
            if J == []:
                return M([])
            C = Compositions()
            C_size = Compositions(J.size())
            return M.sum_of_terms((C(I), number_of_fCT(C(I), J)) for I in C_size)

        @cached_method
        def _matrix_monomial_to_dual_immaculate(self, n):
            r"""
            This function caches the change of basis matrix from the quasisymmetric monomial basis to the dual immaculate basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT:

            - A list. Each entry in the list is a row in the change of basis matrix.

            EXAMPLES::

                sage: dI = QuasiSymmetricFunctions(QQ).dI()
                sage: dI._matrix_monomial_to_dual_immaculate(3)
                [[1, -1, -1, 1], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 0, 1]]
                sage: dI._matrix_monomial_to_dual_immaculate(0)
                [[1]]
            """
            N = NonCommutativeSymmetricFunctions(self.base_ring())
            I = N.I()
            S = N.S()
            mat = []
            C = Compositions()
            for alp in Compositions(n):
                row = []
                expansion = S(I(C(alp)))
                for bet in Compositions(n):
                    row.append(expansion.coefficient(C(bet)))
                mat.append(row)
            return mat

        def _from_Monomial_on_basis(self, J):
            r"""
            Given a quasi-symmetric monomial function, this method returns the expansion into the dual immaculate basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT:

            - A quasi-symmetric function in the dual immaculate basis.

            EXAMPLES:

                sage: dI = QuasiSymmetricFunctions(QQ).dI()
                sage: dI._from_Monomial_on_basis(Composition([]))
                dI[]
                sage: dI._from_Monomial_on_basis(Composition([2,1]))
                -dI[1, 1, 1] - dI[1, 2] + dI[2, 1]
                sage: dI._from_Monomial_on_basis(Composition([3,1,2]))
                -dI[1, 1, 1, 1, 1, 1] + dI[1, 1, 1, 1, 2] + dI[1, 1, 1, 3] - dI[1, 1, 4] - dI[1, 2, 1, 1, 1] + dI[1, 2, 3] + dI[2, 1, 1, 1, 1] - dI[2, 1, 1, 2] + dI[2, 2, 1, 1] - dI[2, 2, 2] - dI[3, 1, 1, 1] + dI[3, 1, 2]
            """
            n = J.size()
            C = Compositions()
            C_n = Compositions(n)
            mat = self._matrix_monomial_to_dual_immaculate(n)
            column = C_n.list().index(J)
            return self.sum_of_terms( (C(I), mat[C_n.list().index(I)][column])
                                            for I in C_n)

    dI = dualImmaculate
