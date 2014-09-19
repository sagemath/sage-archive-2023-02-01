"""
Lie Algebras Given By Structure Coefficients

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

from copy import copy
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc import repr_lincomb
from sage.structure.indexed_generators import IndexedGenerators
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper

from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.finite_dimensional_lie_algebras_with_basis import FiniteDimensionalLieAlgebrasWithBasis

from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import (LieGenerator,
    LieBracket, LieAlgebraElement)
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
from sage.algebras.lie_algebras.ideal import LieAlgebraIdeal
from sage.algebras.lie_algebras.quotient import QuotientLieAlgebra
from sage.rings.all import ZZ
from sage.rings.ring import Ring
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.infinity import infinity
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.modules.free_module import FreeModule, span
from sage.sets.family import Family, AbstractFamily

# TODO: Much of this could be moved to (FiniteDimensional)LieAlgebrasWithBasis
class LieAlgebraWithStructureCoefficients(FinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    A Lie algebra with a set of specified structure coefficients.

    The structure coefficients are specified as a dictionary whose keys are
    pairs of generators and values are dictionaries of generators mapped
    to coefficients.

    EXAMPLES:

    We create the Lie algebra of `\QQ^3` under the Lie bracket defined
    by `\times` (cross-product)::

        sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
        sage: (x,y,z) = L.gens()
        sage: L.bracket(x, y)
        z
        sage: L.bracket(y, x)
        -z
    """
    @staticmethod
    def __classcall_private__(cls, R, s_coeff, names=None, index_set=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L2 = LieAlgebra(QQ, 'x,y', {('y','x'):{'x':-1}})
            sage: L is L2
            True
        """
        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(s_coeff)
        if len(s_coeff) == 0:
            return AbelianLieAlgebra(R, names, index_set, **kwds)

        if names is None:
            if index_set is None:
                raise ValueError("either the names or the index set must be specified")
            if len(index_set) <= 1:
                return AbelianLieAlgebra(R, names, index_set, **kwds)
        elif len(names) <= 1:
            return AbelianLieAlgebra(R, names, index_set, **kwds)

        return super(LieAlgebraWithStructureCoefficients, cls).__classcall__(
            cls, R, s_coeff, tuple(names), index_set, **kwds)

    @staticmethod
    def _standardize_s_coeff(s_coeff):
        """
        Helper function to standardize ``s_coeff`` into the appropriate tuple
        of tuples. Strips items with coefficients of 0 and duplicate entries.
        This does not check the Jacobi relation (nor antisymmetry if the
        cardinality is infinite).
        """
        # Try to handle infinite basis (once/if supported)
        if isinstance(s_coeff, AbstractFamily) and s_coeff.cardinality() == infinity:
            return s_coeff

        sc = {}
        # Make sure the first gen is smaller than the second in each key
        for k in s_coeff.keys():
            v = s_coeff[k]
            if isinstance(v, dict):
                v = v.items()

            if k[0] > k[1]:
                key = LieBracket(k[1], k[0])
                vals = tuple((g, -val) for g, val in v if val != 0)
            else:
                key = LieBracket(*k)
                vals = tuple((g, val) for g, val in v if val != 0)

            if key in sc.keys() and sorted(sc[key]) != sorted(vals):
                raise ValueError("non-equal brackets")

            if len(vals) > 0:
                sc[key] = vals
        return Family(sc)

    def __init__(self, R, s_coeff, names, index_set, category=None, prefix='',
                 bracket=False, latex_bracket=False, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: TestSuite(L).run()
        """
        cat = LieAlgebras(R).WithBasis().or_subcategory(category)
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, index_set, cat)
        IndexedGenerators.__init__(self, self._indices, prefix=prefix, bracket=bracket,
                                   latex_bracket=latex_bracket, **kwds)

        # Transform the values in the structure coefficients to elements
        self._s_coeff = Family({k: self._from_dict(dict(s_coeff[k])) for k in s_coeff.keys()})

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.basis()
            Finite family {'y': y, 'x': x}
        """
        return Family({i: self.monomial(i) for i in self._indices})

    def structure_coefficients(self):
        """
        Return the dictonary of structure coefficients of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.structure_coefficients()
            Finite family {[x, y]: x}
        """
        return self._s_coeff

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.dimension()
            2
        """
        return self.basis().cardinality()

    # Move to FiniteDimensionalLieAlgebrasWithBasis category?
    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket of ``[self, y]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
            sage: L.bracket(x, y) # indirect doctest
            z
            sage: L.bracket(y, x)
            -z
            sage: L.bracket(x + y - z, x - y + z)
            -2*y - 2*z
        """
        comp = self._print_options['generator_cmp']
        ordered = True
        if comp(x, y) > 0: # x > y
            x,y = y,x
            ordered = False
        b = LieBracket(x, y)
        try:
            val = self._s_coeff[b]
        except KeyError:
            return self.zero()
        if ordered:
            return val
        return -val

    def free_module(self, sparse=True):
        """
        Return ``self`` as a free module.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
            sage: L.free_module()
            Sparse vector space of dimension 3 over Rational Field
        """
        return FreeModule(self.base_ring(), self.dimension(), sparse=sparse)

    def subalgebra(self, gens, names=None, index_set=None, category=None):
        """
        Return the subalgebra of ``self`` generated by ``gens``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
            sage: L.subalgebra([x+y])
            Subalgebra generated of Lie algebra on 3 generators (x, y, z) over Rational Field with basis:
            (x + y,)
            sage: L.subalgebra([x+y, y+z])
            Lie algebra on 3 generators (x, y, z) over Rational Field

        TESTS::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
            sage: L.subalgebra([x+y, x+y])
            Subalgebra generated of Lie algebra on 3 generators (x, y, z) over Rational Field with basis:
            (x + y,)
        """
        B, s_coeff = LieSubalgebraWithStructureCoefficients._compute_basis_structure_coeff(self, gens)
        if len(B) == self.dimension():
            return self
        return LieSubalgebraWithStructureCoefficients(self, B, s_coeff, names, index_set, category)

    # TODO:
    # - Make a centralizer of a set version
    # - Move to FiniteDimensionalLieAlgebrasWithBasis once implemented
    # - Return a subalgebra
    #@cached_method
    #def center(self):
    #    """
    #    Return a list of elements which correspond to a basis for the center
    #    of ``self``.
    #    """
    #    R = self.base_ring()
    #    B = self.basis()
    #    K = list(B.keys())
    #    k = len(K)
    #    d = {}
    #    for a,i in enumerate(K):
    #        Bi = B[i]
    #        for b,j in enumerate(K):
    #            Bj = B[j]
    #            for m,c in Bi.bracket(Bj):
    #                d[(a, K.index(m)+k*b)] = c
    #    m = Matrix(R, d, nrows=k, ncols=k*k, sparse=True)
    #    from_vector = lambda x: self.sum_of_terms( ((K[i], c) for i,c in x.iteritems()),
    #                                               distinct=True)
    #    return tuple(map( from_vector, m.kernel().basis() ))
        # Dense version
        # R = self.base_ring()
        # B = self.basis()
        # K = list(B.keys())
        # eqns = [[] for dummy in range(k)]
        # for a,i in enumerate(K):
        #     for b,j in enumerate(K):
        #         v = B[i]*B[j] - B[j]*B[i]
        #         eqns[a].extend([v[k] for k in K])
        # m = Matrix(R, eqns)
        # from_vector = lambda x: self.sum_of_terms(((K[i], c) for i,c in x.iteritems()),
        #                                           distinct=True)
        # return tuple(map( from_vector, m.kernel().basis() ))

    class Element(LieAlgebraElement):
        """
        An element of a Lie algebra given by structure coefficients.
        """
        def _repr_(self):
            r"""
            Return a string representation of ``self``.
            """
            prefix = self.parent()._print_options['prefix']
            if self.parent()._names is not None and not prefix:
                return LieAlgebraElement._repr_(self)
            return repr_lincomb(self._sorted_items_for_printing(),
                                scalar_mult=self.parent()._print_options['scalar_mult'],
                                repr_monomial = self.parent()._repr_generator,
                                strip_one = True)

        def _latex_(self):
            r"""
            Return a `\LaTeX` representation of ``self``.
            """
            prefix = self.parent()._print_options['prefix']
            if self.parent()._names is not None and not prefix:
                return LieAlgebraElement._latex_(self)
            return repr_lincomb(self._sorted_items_for_printing(),
                                scalar_mult       = self.parent()._print_options['scalar_mult'],
                                latex_scalar_mult = self.parent()._print_options['latex_scalar_mult'],
                                repr_monomial = self.parent()._latex_generator,
                                is_latex=True, strip_one = True)

        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.
            """
            P = self.parent()
            return P.sum(cx * cy * P.bracket_on_basis(mx, my)
                         for mx,cx in self for my,cy in y)

        def to_vector(self):
            """
            Return ``self`` as a vector.
            """
            V = self.parent().free_module()
            return V([self[k] for k in self.parent()._ordered_indices])

class LieSubalgebraWithStructureCoefficients(LieSubalgebra):
    """
    A Lie subalgebra of a Lie algebra given by structure coefficients.
    """
    @staticmethod
    def _compute_basis_structure_coeff(A, gens):
        """
        Compute an echonalized basis of ``self`` and the corresponding
        structure coefficients.

        INPUT:

        - ``A`` -- the ambient Lie algebra
        - ``gens`` -- a set of generators as elements in ``A``
        """
        I = A._ordered_indices
        R = A.base_ring()
        M = A.free_module()
        zero = A.zero()

        # Setup the current basis matrix
        cur = matrix(R, map(lambda x: x.to_vector(), gens), sparse=M.is_sparse())
        cur.echelonize()
        # Remove all zero rows from cur
        cur = cur.delete_rows([i for i,row in enumerate(cur.rows()) if row.is_zero()])

        added = True
        while added:
            added = False
            cur.echelonize()
            basis = tuple(A._from_dict({I[i]: c for i,c in row.iteritems()}) for row in cur)
            s_coeff = {}

            for i in range(len(basis)):
                for j in range(i+1, len(basis)):
                    b = basis[i].bracket(basis[j])
                    if b == zero:
                        # There is a dependence and no resulting structure
                        #   coefficients, so we don't have to do anything more
                        continue

                    bv = b.to_vector()
                    dep = M.linear_dependence(cur.rows() + [bv])
                    if dep:
                        # The dependency is the first in the sequence, there's only
                        #   one because the original set is linearly independent
                        dep = dep[0]
                        # Rewrite as a linear combination of the current basis
                        s = -dep[-1] # Guaranteed to be non-zero
                        d = {k: R(c*~s) for k,c in enumerate(dep[:-1]) if c != 0}
                        s_coeff[LieBracket(i,j)] = d
                    else:
                        added = True
                        cur = cur.stack(bv)
        return (basis, Family(s_coeff))

    def __init__(self, ambient, basis, s_coeff, names=None, index_set=None, category=None):
        r"""
        Initialize ``self``.
        """
        self._s_coeff = s_coeff # TODO: convert the values to elements of ``self``
        # TODO: make sure the category is in subobjects?
        LieSubalgebra.__init__(self, ambient, basis, names, index_set, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Subalgebra generated of {} with basis:\n{}".format(self._ambient, self.gens())

    def basis(self):
        """
        Return a basis of ``self``.
        """
        I = self._indices
        B = self._gens
        return Family({ I[i]: b for i,b in enumerate(B) })

    def subalgebra(self, gens, names=None, index_set=None, category=None):
        """
        Return the subalgebra of ``self`` generated by ``gens``.
        """
        gens = map(lambda x: x.value, gens)
        B, s_coeff = LieSubalgebraWithStructureCoefficients._compute_basis_structure_coeff(self._ambient, gens)
        assert len(B) != self._ambient.dimension()
        return LieSubalgebraWithStructureCoefficients(self._ambient, B, s_coeff,
                                                      names, index_set, category)

class LieAlgebraIdealWithStructureCoefficients(LieAlgebraIdeal,
                        LieSubalgebraWithStructureCoefficients):
    """
    A Lie algebra ideal of a Lie algebra given by structure coefficients.
    """
    def reduce(self, y):
        """
        Return ``y`` modulo ``self``.
        """
        M = self._ambient.free_module()
        I = M.subspace(map(lambda x: x.to_vector(), self._gens))
        R = I.complement()
        ld = M.linear_dependence(R.basis() + I.basis() + [y.to_vector()])[0]
        normalizer = -ld[-1]
        coeffs = map(lambda x: x / normalizer, ld[:len(R)])
        vec = M(R.linear_combination_of_basis(coeffs))
        return self._ambient(vec)

# This should either not inherit from QuotientLieAlgebra or refactor out common code
# This should inherit from LieAlgebraWithStructureCoefficients
class QuotientLieAlgebraWithStructureCoefficients(QuotientLieAlgebra):
    """
    A quotient Lie algebra of a Lie algebra given by structure coefficients.

    INPUT:

    - ``I`` -- an ideal of a Lie algebra given by structure coefficients
    - ``names`` -- the names of the generators
    - ``index_set`` -- the indexing set for the generators
    - ``category`` -- the category of the quotient Lie algebra
    """
    def __init__(self, lie, I, names=None, index_set=None, category=None):
        """
        Initialize ``self``.
        """
        R = lie.base_ring()
        #cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        cat = FiniteDimensionalLieAlgebrasWithBasis(R)#.Quotients()
        QuotientLieAlgebra.__init__(self, lie, I, names, index_set, cat)

        # Construct the structure coefficients
        M = lie.free_module()
        IM = M.subspace(map(lambda x: x.to_vector(), I._gens))
        RM = IM.complement()
        B = RM.basis() + IM.basis()
        dim = len(RM)
        self._s_coeff = {}
        for i,kl,zl in enumerate(gens):
            for kr,zr in gens[i+1:]:
                b = LieBracket(zl, zr)
                ld = M.linear_dependence(B + [y.to_vector()])[0]
                normalizer = -ld[-1]
                self._s_coeffs[b] = {indices[j]: val / normalizer
                                      for j,val in enumerate(ld[:dim])}

    Element = LieAlgebraElement

class AbelianLieAlgebra(LieAlgebraWithStructureCoefficients):
    r"""
    An abelian Lie algebra.

    A Lie algebra `\mathfrak{g}` is abelian if `[x, y] = 0` for all
    `x, y \in \mathfrak{g}`.

    EXAMPLES::

        sage: L.<x, y> = LieAlgebra(QQ, abelian=True)
        sage: L.bracket(x, y)
        0
    """
    def __init__(self, R, names=None, index_set=None, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: TestSuite(L).run()
        """
        LieAlgebraWithStructureCoefficients.__init__(self, R, Family({}), names, index_set, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, 3, 'x', abelian=True)
            Abelian Lie algebra on 3 generators (x0, x1, x2) over Rational Field
        """
        gens = self.lie_algebra_generators()
        if gens.cardinality() == 1:
            return "Abelian Lie algebra on generator {} over {}".format(tuple(gens)[0], self.base_ring())
        return "Abelian Lie algebra on {} generators {} over {}".format(
            gens.cardinality(), tuple(gens), self.base_ring())

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: L._construct_UEA()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return PolynomialRing(self.base_ring(), self.variable_names())

    def is_abelian(self):
        """
        Return ``True`` since this is an abelian Lie algebra.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: L.is_abelian()
            True
        """
        return True

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket of basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.bracket_on_basis(x.leading_support(), y.leading_support())
            0
        """
        return self.zero()

    class Element(LieAlgebraElement):
        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::

                sage: L.<x, y> = LieAlgebra(QQ, abelian=True)
                sage: L.bracket(x, y)
                0
            """
            return self.parent().zero()
