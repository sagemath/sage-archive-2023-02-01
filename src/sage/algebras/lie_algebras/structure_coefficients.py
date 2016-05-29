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

#from copy import copy
#from sage.misc.cachefunc import cached_method
#from sage.misc.lazy_attribute import lazy_attribute
#from sage.misc.misc import repr_lincomb
from sage.structure.indexed_generators import IndexedGenerators
#from sage.structure.parent import Parent
#from sage.structure.unique_representation import UniqueRepresentation
#from sage.structure.element_wrapper import ElementWrapper

#from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras

#from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra, FinitelyGeneratedLieAlgebra
#from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
#from sage.algebras.lie_algebras.ideal import LieAlgebraIdeal
#from sage.algebras.lie_algebras.quotient import QuotientLieAlgebra
#from sage.rings.all import ZZ
#from sage.rings.ring import Ring
#from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
#from sage.rings.infinity import infinity
#from sage.matrix.matrix_space import MatrixSpace
#from sage.matrix.constructor import matrix
#from sage.modules.free_module_element import vector
from sage.modules.free_module import FreeModule #, span
from sage.sets.family import Family #, AbstractFamily

class LieAlgebraWithStructureCoefficients(FinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    A Lie algebra with a set of specified structure coefficients.

    The structure coefficients are specified as a dictionary `d` whose
    keys are pairs of basis indices, and whose values are
    dictionaries which in turn are indexed by basis indices. The value
    of `d` at a pair `(u, v)` of basis indices is the dictionary whose
    `w`-th entry (for `w` a basis index) is the coefficient of `b_w`
    in the Lie bracket `[b_u, b_v]` (where `b_x` means the basis
    element with index `x`).

    INPUT:

    - ``R`` -- a ring, to be used as the base ring

    - ``s_coeff`` -- a dictionary, indexed by pairs of basis indices
      (see below), and whose values are dictionaries which are
      indexed by (single) basis indices and whose values are elements
      of `R`

    - ``names`` -- list or tuple of strings

    - ``index_set`` -- (default: ``names``) list or tuple of hashable
      and comparable elements

    OUTPUT:

    A Lie algebra over ``R`` which (as an `R`-module) is free with
    a basis indexed by the elements of ``index_set``. The `i`-th
    basis element is displayed using the name ``names[i]``.
    If we let `b_i` denote this `i`-th basis element, then the Lie
    bracket is given by the requirement that the `b_k`-coefficient
    of `[b_i, b_j]` is ``s_coeff[(i, j)][k]`` if
    ``s_coeff[(i, j)]`` exists, otherwise ``-s_coeff[(j, i)][k]``
    if ``s_coeff[(j, i)]`` exists, otherwise `0`.

    EXAMPLES:

    We create the Lie algebra of `\QQ^3` under the Lie bracket defined
    by `\times` (cross-product)::

        sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
        sage: (x,y,z) = L.gens()
        sage: L.bracket(x, y)
        z
        sage: L.bracket(y, x)
        -z

    TESTS:

    We can input structure coefficients that fail the Jacobi
    identity, but the test suite will call us out on it::

        sage: Fake = LieAlgebra(QQ, 'x,y,z', {('x','y'):{'z':3}, ('y','z'):{'z':1}, ('z','x'):{'y':1}})
        sage: TestSuite(Fake).run()
        Failure in _test_jacobi_identity:
        ...
    """
    @staticmethod
    def __classcall_private__(cls, R, s_coeff, names=None, index_set=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L1 = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L2 = LieAlgebra(QQ, 'x,y', {('y','x'):{'x':-1}})
            sage: L1 is L2
            True

        Check that we convert names to the indexing set::

            sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}}, index_set=range(3))
            sage: (x,y,z) = L.gens()
            sage: L[x,y]
            L[2]
        """
        names, index_set = LieAlgebra._standardize_names_index_set(names, index_set)

        # Make sure the structure coefficients are given by the indexing set
        if names is not None and names != tuple(index_set):
            d = {x: index_set[i] for i,x in enumerate(names)}
            get_pairs = lambda X: X.items() if isinstance(X, dict) else X
            s_coeff = {(d[k[0]], d[k[1]]): [(d[x], y) for x,y in get_pairs(s_coeff[k])]
                       for k in s_coeff}

        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(s_coeff)
        if s_coeff.cardinality() == 0:
            return AbelianLieAlgebra(R, names, index_set, **kwds)

        if (names is None and len(index_set) <= 1) or len(names) <= 1:
            return AbelianLieAlgebra(R, names, index_set, **kwds)

        return super(LieAlgebraWithStructureCoefficients, cls).__classcall__(
            cls, R, s_coeff, names, index_set, **kwds)

    @staticmethod
    def _standardize_s_coeff(s_coeff):
        """
        Helper function to standardize ``s_coeff`` into the appropriate form
        (dictionary indexed by pairs, whose values are dictionaries).
        Strips items with coefficients of 0 and duplicate entries.
        This does not check the Jacobi relation (nor antisymmetry if the
        cardinality is infinite).

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
            sage: d = {('y','x'): {'x':-1}}
            sage: LieAlgebraWithStructureCoefficients._standardize_s_coeff(d)
            Finite family {('x', 'y'): (('x', 1),)}
        """
        # Try to handle infinite basis (once/if supported)
        #if isinstance(s_coeff, AbstractFamily) and s_coeff.cardinality() == infinity:
        #    return s_coeff

        sc = {}
        # Make sure the first gen is smaller than the second in each key
        for k in s_coeff.keys():
            v = s_coeff[k]
            if isinstance(v, dict):
                v = v.items()

            if k[0] > k[1]:
                key = (k[1], k[0])
                vals = tuple((g, -val) for g, val in v if val != 0)
            else:
                if not k[0] < k[1]:
                    if k[0] == k[1]:
                        if not all(val == 0 for g, val in v):
                            raise ValueError("elements {} are equal but their bracket is not set to 0".format(k))
                        continue
                    if not k[0] <= k[1]:
                        raise ValueError("elements {} are not comparable".format(k))
                key = tuple(k)
                vals = tuple((g, val) for g, val in v if val != 0)

            if key in sc.keys() and sorted(sc[key]) != sorted(vals):
                raise ValueError("two distinct values given for one and the same bracket")

            if len(vals) > 0:
                sc[key] = vals
        return Family(sc)

    def __init__(self, R, s_coeff, names, index_set, category=None, prefix=None,
                 bracket=None, latex_bracket=None, string_quotes=None, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'): {'x':1}})
            sage: TestSuite(L).run()
        """
        default = (names != tuple(index_set))
        if prefix is None:
            if default:
                prefix = 'L'
            else:
                prefix = ''
        if bracket is None:
            bracket = default
        if latex_bracket is None:
            latex_bracket = default
        if string_quotes is None:
            string_quotes = default

        cat = LieAlgebras(R).WithBasis().FiniteDimensional().or_subcategory(category)
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, index_set, cat)
        IndexedGenerators.__init__(self, self._indices, prefix=prefix,
                                   bracket=bracket, latex_bracket=latex_bracket,
                                   string_quotes=string_quotes, **kwds)

        # Transform the values in the structure coefficients to elements
        self._s_coeff = Family({k: self._from_dict(dict(s_coeff[k]))
                                for k in s_coeff.keys()})

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.basis()
            Finite family {'y': y, 'x': x}
        """
        return Family(self._indices, self.monomial)

    def structure_coefficients(self, include_zeros=False):
        """
        Return the dictonary of structure coefficients of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'):{'x':1}})
            sage: L.structure_coefficients()
            Finite family {('x', 'y'): x}
            sage: S = L.structure_coefficients(True); S
            Finite family {('x', 'y'): x, ('x', 'z'): 0, ('y', 'z'): 0}
            sage: S['x','z'].parent() is L
            True
        """
        if not include_zeros:
            return self._s_coeff
        ret = {}
        zero = self.zero()
        S = dict(self._s_coeff)
        for i,x in enumerate(self._indices):
            for y in self._indices[i+1:]:
                if x < y:
                    b = (x, y)
                else:
                    b = (y, x)
                ret[b] = S.get(b, zero)
        return Family(ret)

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.dimension()
            2
        """
        return self.basis().cardinality()

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket ``[x, y]`` of two basis elements
        (indexed by) ``x`` and ``y`` where ``x < y``.

        (This particular implementation actually does not require
        ``x < y``.)

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
            sage: L.bracket_on_basis('x', 'y')
            z
            sage: L.bracket_on_basis('y', 'x')
            -z
            sage: L.bracket(x + y - z, x - y + z)
            -2*y - 2*z
        """
        ordered = True
        if x > y:
            x,y = y,x
            ordered = False
        b = (x, y)
        try:
            val = self._s_coeff[b]
        except KeyError:
            return self.zero()
        if ordered:
            return val
        return -val

    def module(self, sparse=True):
        """
        Return ``self`` as a free module.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: L.module()
            Sparse vector space of dimension 3 over Rational Field
        """
        return FreeModule(self.base_ring(), self.dimension(), sparse=sparse)

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.three_dimensional(QQ, 4, 1, -1, 2)
            sage: L.some_elements()
            [X, Y, Z, X + Y + Z]
        """
        return list(self.basis()) + [self.sum(self.basis())]

    class Element(LieAlgebraElement):
        """
        An element of a Lie algebra given by structure coefficients.
        """
        def to_vector(self):
            """
            Return ``self`` as a vector.

            EXAMPLES::

                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
                sage: a = x + 3*y - z/2
                sage: a.to_vector()
                (1, 3, -1/2)
            """
            V = self.parent().module()
            return V([self[k] for k in self.parent()._ordered_indices])

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
    @staticmethod
    def __classcall_private__(cls, R, names=None, index_set=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: L1 = LieAlgebra(QQ, 'x,y', {})
            sage: L2.<x, y> = LieAlgebra(QQ, abelian=True)
            sage: L1 is L2
            True
        """
        names, index_set = LieAlgebra._standardize_names_index_set(names, index_set)
        return super(AbelianLieAlgebra, cls).__classcall__(cls, R, names, index_set, **kwds)

    def __init__(self, R, names, index_set, **kwds):
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

