"""
Finitely Generated Matrix Groups

This class is designed for computing with matrix groups defined by a
finite set of generating matrices.

EXAMPLES::

    sage: F = GF(3)
    sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F,2, [1,1,0,1])]
    sage: G = MatrixGroup(gens)
    sage: G.conjugacy_classes_representatives()
    (
    [1 0]  [0 2]  [0 1]  [2 0]  [0 2]  [0 1]  [0 2]
    [0 1], [1 1], [2 1], [0 2], [1 2], [2 2], [1 0]
    )

The finitely generated matrix groups can also be constructed as
subgroups of matrix groups::

    sage: SL2Z = SL(2,ZZ)
    sage: S, T = SL2Z.gens()
    sage: SL2Z.subgroup([T^2])
    Matrix group over Integer Ring with 1 generators (
    [1 2]
    [0 1]
    )

AUTHORS:

- William Stein: initial version

- David Joyner (2006-03-15): degree, base_ring, _contains_, list,
  random, order methods; examples

- William Stein (2006-12): rewrite

- David Joyner (2007-12): Added invariant_generators (with Martin
  Albrecht and Simon King)

- David Joyner (2008-08): Added module_composition_factors (interface
  to GAP's MeatAxe implementation) and as_permutation_group (returns
  isomorphic PermutationGroup).

- Simon King (2010-05): Improve invariant_generators by using GAP
  for the construction of the Reynolds operator in Singular.

- Volker Braun (2013-1) port to new Parent, libGAP.
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################
from __future__ import print_function

from sage.groups.group import Group
from sage.rings.all import ZZ
from sage.rings.integer import is_Integer
from sage.rings.ring import is_Ring
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.interfaces.gap import gap
from sage.matrix.matrix import is_Matrix
from sage.matrix.matrix_space import MatrixSpace, is_MatrixSpace
from sage.matrix.all import matrix
from sage.misc.latex import latex
from sage.structure.sequence import Sequence
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.arith.all import gcd
from sage.rings.fraction_field import FractionField
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.number_field.number_field import CyclotomicField

from sage.groups.matrix_gps.matrix_group import (
    is_MatrixGroup, MatrixGroup_generic, MatrixGroup_gap )
from sage.groups.matrix_gps.group_element import (
    is_MatrixGroupElement, MatrixGroupElement_generic, MatrixGroupElement_gap)




def normalize_square_matrices(matrices):
    """
    Find a common space for all matrices.

    OUTPUT:

    A list of matrices, all elements of the same matrix space.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.finitely_generated import normalize_square_matrices
        sage: m1 = [[1,2],[3,4]]
        sage: m2 = [2, 3, 4, 5]
        sage: m3 = matrix(QQ, [[1/2,1/3],[1/4,1/5]])
        sage: m4 = MatrixGroup(m3).gen(0)
        sage: normalize_square_matrices([m1, m2, m3, m4])
        [
        [1 2]  [2 3]  [1/2 1/3]  [1/2 1/3]
        [3 4], [4 5], [1/4 1/5], [1/4 1/5]
        ]
    """
    deg = []
    gens = []
    for m in matrices:
        if is_MatrixGroupElement(m):
            deg.append(m.parent().degree())
            gens.append(m.matrix())
            continue
        if is_Matrix(m):
            if not m.is_square():
                raise TypeError('matrix must be square')
            deg.append(m.ncols())
            gens.append(m)
            continue
        try:
            m = list(m)
        except TypeError:
            gens.append(m)
            continue
        if isinstance(m[0], (list, tuple)):
            m = [list(_) for _ in m]
            degree = ZZ(len(m))
        else:
            degree, rem = ZZ(len(m)).sqrtrem()
            if rem!=0:
                raise ValueError('list of plain numbers must have square integer length')
        deg.append(degree)
        gens.append(matrix(degree, degree, m))
    deg = set(deg)
    if len(set(deg)) != 1:
        raise ValueError('not all matrices have the same size')
    gens = Sequence(gens, immutable=True)
    MS = gens.universe()
    if not is_MatrixSpace(MS):
        raise TypeError('all generators must be matrices')
    if MS.nrows() != MS.ncols():
        raise ValueError('matrices must be square')
    return gens

def QuaternionMatrixGroupGF3():
    r"""
    The quaternion group as a set of `2\times 2` matrices over `GF(3)`.

    OUTPUT:

    A matrix group consisting of `2\times 2` matrices with
    elements from the finite field of order 3.  The group is
    the quaternion group, the nonabelian group of order 8 that
    is not isomorphic to the group of symmetries of a square
    (the dihedral group `D_4`).

    .. note::
        This group is most easily available via ``groups.matrix.QuaternionGF3()``.

    EXAMPLES:

    The generators are the matrix representations of the
    elements commonly called `I` and `J`, while `K`
    is the product of `I` and `J`. ::

        sage: from sage.groups.matrix_gps.finitely_generated import QuaternionMatrixGroupGF3
        sage: Q = QuaternionMatrixGroupGF3()
        sage: Q.order()
        8
        sage: aye = Q.gens()[0]; aye
        [1 1]
        [1 2]
        sage: jay = Q.gens()[1]; jay
        [2 1]
        [1 1]
        sage: kay = aye*jay; kay
        [0 2]
        [1 0]

    TESTS::

        sage: groups.matrix.QuaternionGF3()
        Matrix group over Finite Field of size 3 with 2 generators (
        [1 1]  [2 1]
        [1 2], [1 1]
        )

        sage: Q = QuaternionMatrixGroupGF3()
        sage: QP = Q.as_permutation_group()
        sage: QP.is_isomorphic(QuaternionGroup())
        True
        sage: H = DihedralGroup(4)
        sage: H.order()
        8
        sage: QP.is_abelian(), H.is_abelian()
        (False, False)
        sage: QP.is_isomorphic(H)
        False
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.matrix.matrix_space import MatrixSpace
    MS = MatrixSpace(FiniteField(3), 2)
    aye = MS([1,1,1,2])
    jay = MS([2,1,1,1])
    return MatrixGroup([aye, jay])

def MatrixGroup(*gens, **kwds):
    r"""
    Return the matrix group with given generators.

    INPUT:

    - ``*gens`` -- matrices, or a single list/tuple/iterable of
      matrices, or a matrix group.

    - ``check`` -- boolean keyword argument (optional, default:
      ``True``). Whether to check that each matrix is invertible.

    EXAMPLES::

        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 5 with 2 generators (
        [1 2]  [1 1]
        [4 1], [0 1]
        )

    In the second example, the generators are a matrix over
    `\ZZ`, a matrix over a finite field, and the integer
    `2`. Sage determines that they both canonically map to
    matrices over the finite field, so creates that matrix group
    there::

        sage: gens = [matrix(2,[1,2, -1, 1]), matrix(GF(7), 2, [1,1, 0,1]), 2]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 7 with 3 generators (
        [1 2]  [1 1]  [2 0]
        [6 1], [0 1], [0 2]
        )

    Each generator must be invertible::

        sage: G = MatrixGroup([matrix(ZZ,2,[1,2,3,4])])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: MatrixGroup([MS.0])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix
        sage: MatrixGroup([MS.0], check=False)  # works formally but is mathematical nonsense
        Matrix group over Finite Field of size 5 with 1 generators (
        [1 0]
        [0 0]
        )

    Some groups are not supported, or do not have much functionality
    implemented::

        sage: G = SL(0, QQ)
        Traceback (most recent call last):
        ...
        ValueError: the degree must be at least 1

        sage: SL2C = SL(2, CC);  SL2C
        Special Linear Group of degree 2 over Complex Field with 53 bits of precision
        sage: SL2C.gens()
        Traceback (most recent call last):
        ...
        AttributeError: 'LinearMatrixGroup_generic_with_category' object has no attribute 'gens'
    """
    if isinstance(gens[-1], dict):   # hack for unpickling
        kwds.update(gens[-1])
        gens = gens[:-1]
    check = kwds.get('check', True)
    if len(gens) == 1:
        if isinstance(gens[0], (list, tuple)):
            gens = list(gens[0])
        else:
            try:
                gens = [g.matrix() for g in gens[0]]
            except AttributeError:
                pass
    if len(gens) == 0:
        raise ValueError('need at least one generator')
    gens = normalize_square_matrices(gens)
    if check and any(not g.is_invertible() for g in gens):
        raise ValueError('each generator must be an invertible matrix')
    MS = gens.universe()
    base_ring = MS.base_ring()
    degree = ZZ(MS.ncols())   # == MS.nrows()
    from sage.libs.gap.libgap import libgap
    try:
        gap_gens = [libgap(matrix_gen) for matrix_gen in gens]
        gap_group = libgap.Group(gap_gens)
        return FinitelyGeneratedMatrixGroup_gap(degree, base_ring, gap_group)
    except (TypeError, ValueError):
        return FinitelyGeneratedMatrixGroup_generic(degree, base_ring, gens)

###################################################################
#
# Matrix group over a generic ring
#
###################################################################

class FinitelyGeneratedMatrixGroup_generic(MatrixGroup_generic):
    """
    TESTS::

        sage: m1 = matrix(SR, [[1,2],[3,4]])
        sage: m2 = matrix(SR, [[1,3],[-1,0]])
        sage: MatrixGroup(m1) == MatrixGroup(m1)
        True
        sage: MatrixGroup(m1) == MatrixGroup(m1.change_ring(QQ))
        False
        sage: MatrixGroup(m1) == MatrixGroup(m2)
        False
        sage: MatrixGroup(m1, m2) == MatrixGroup(m2, m1)
        False

        sage: m1 = matrix(QQ, [[1,2],[3,4]])
        sage: m2 = matrix(QQ, [[1,3],[-1,0]])
        sage: MatrixGroup(m1) == MatrixGroup(m1)
        True
        sage: MatrixGroup(m1) == MatrixGroup(m2)
        False
        sage: MatrixGroup(m1, m2) == MatrixGroup(m2, m1)
        False

        sage: G = GL(2, GF(3))
        sage: H = G.as_matrix_group()
        sage: H == G, G == H
        (True, True)
    """

    def __init__(self, degree, base_ring, generator_matrices, category=None):
        """
        Matrix group generated by a finite number of matrices.

        EXAMPLES::

            sage: m1 = matrix(SR, [[1,2],[3,4]])
            sage: m2 = matrix(SR, [[1,3],[-1,0]])
            sage: G = MatrixGroup(m1, m2)
            sage: TestSuite(G).run()
            sage: type(G)
            <class 'sage.groups.matrix_gps.finitely_generated.FinitelyGeneratedMatrixGroup_generic_with_category'>

            sage: from sage.groups.matrix_gps.finitely_generated import \
            ....:     FinitelyGeneratedMatrixGroup_generic
            sage: G = FinitelyGeneratedMatrixGroup_generic(2, QQ, [matrix(QQ,[[1,2],[3,4]])])
            sage: G.gens()
            (
            [1 2]
            [3 4]
            )
        """
        self._gens_matrix = generator_matrices
        MatrixGroup_generic.__init__(self, degree, base_ring, category=category)

    @cached_method
    def gens(self):
        """
        Return the generators of the matrix group.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]), MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: gens[0] in G
            True
            sage: gens = G.gens()
            sage: gens[0] in G
            True
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators (
            [1 0]  [1 2]
            [0 1], [3 4]
            )
            sage: G.gens()
            (
            [1 0]  [1 2]
            [0 1], [3 4]
            )
        """
        return tuple(self.element_class(self, x, check=False, convert=False)
                     for x in self._gens_matrix)

    def gen(self, i):
        """
        Return the `i`-th generator

        OUTPUT:

        The `i`-th generator of the group.

        EXAMPLES::

            sage: H = GL(2, GF(3))
            sage: h1, h2 = H([[1,0],[2,1]]), H([[1,1],[0,1]])
            sage: G = H.subgroup([h1, h2])
            sage: G.gen(0)
            [1 0]
            [2 1]
            sage: G.gen(0).matrix() == h1.matrix()
            True
        """
        return self.gens()[i]

    def ngens(self):
        """
        Return the number of generators

        OUTPUT:

        An integer. The number of generators.

        EXAMPLES::

            sage: H = GL(2, GF(3))
            sage: h1, h2 = H([[1,0],[2,1]]), H([[1,1],[0,1]])
            sage: G = H.subgroup([h1, h2])
            sage: G.ngens()
            2
        """
        return len(self._gens_matrix)

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: G = MatrixGroup([matrix(CC, [[1,2],[3,4]]),
            ....:                  matrix(CC, [[1,3],[-1,0]])])
            sage: loads(dumps(G)) == G
            True

        Check that :trac:`22128` is fixed::

            sage: R = MatrixSpace(SR, 2)
            sage: G = MatrixGroup([R([[1, 1], [0, 1]])])
            sage: G.register_embedding(R)
            sage: loads(dumps(G))
            Matrix group over Symbolic Ring with 1 generators (
            [1 1]
            [0 1]
            )
        """
        return MatrixGroup, (self._gens_matrix, {'check': False})

    def _test_matrix_generators(self, **options):
        """
        EXAMPLES::

            sage: m1 = matrix(SR, [[1,2],[3,4]])
            sage: m2 = matrix(SR, [[1,3],[-1,0]])
            sage: G = MatrixGroup(m1, m2)
            sage: G._test_matrix_generators()
        """
        tester = self._tester(**options)
        for g,h in zip(self.gens(), MatrixGroup(self.gens()).gens()):
            tester.assertEqual(g.matrix(), h.matrix())

###################################################################
#
# Matrix group over a ring that GAP understands
#
###################################################################

class FinitelyGeneratedMatrixGroup_gap(MatrixGroup_gap):
    """
    Matrix group generated by a finite number of matrices.

    EXAMPLES::

        sage: m1 = matrix(GF(11), [[1,2],[3,4]])
        sage: m2 = matrix(GF(11), [[1,3],[10,0]])
        sage: G = MatrixGroup(m1, m2);  G
        Matrix group over Finite Field of size 11 with 2 generators (
        [1 2]  [ 1  3]
        [3 4], [10  0]
        )
        sage: type(G)
        <class 'sage.groups.matrix_gps.finitely_generated.FinitelyGeneratedMatrixGroup_gap_with_category'>
        sage: TestSuite(G).run()
    """

    def __reduce__(self):
        """
        Implement pickling.

        EXAMPLES::

            sage: m1 = matrix(QQ, [[1,2],[3,4]])
            sage: m2 = matrix(QQ, [[1,3],[-1,0]])
            sage: loads(MatrixGroup(m1, m2).dumps())
            Matrix group over Rational Field with 2 generators (
            [1 2]  [ 1  3]
            [3 4], [-1  0]
            )
        """
        return (MatrixGroup,
                tuple(g.matrix() for g in self.gens()) + ({'check':False},))

    def as_permutation_group(self, algorithm=None):
        r"""
        Return a permutation group representation for the group.

        In most cases occurring in practice, this is a permutation
        group of minimal degree (the degree begin determined from
        orbits under the group action). When these orbits are hard to
        compute, the procedure can be time-consuming and the degree
        may not be minimal.

        INPUT:

        - ``algorithm`` -- ``None`` or ``'smaller'``. In the latter
          case, try harder to find a permutation representation of
          small degree.

        OUTPUT:

        A permutation group isomorphic to ``self``. The
        ``algorithm='smaller'`` option tries to return an isomorphic
        group of low degree, but is not guaranteed to find the
        smallest one.

        EXAMPLES::
        
            sage: MS = MatrixSpace(GF(2), 5, 5)
            sage: A = MS([[0,0,0,0,1],[0,0,0,1,0],[0,0,1,0,0],[0,1,0,0,0],[1,0,0,0,0]])
            sage: G = MatrixGroup([A])
            sage: G.as_permutation_group()
            Permutation Group with generators [(1,2)]
            sage: MS = MatrixSpace( GF(7), 12, 12)
            sage: GG = gap("ImfMatrixGroup( 12, 3 )")
            sage: GG.GeneratorsOfGroup().Length()
            3
            sage: g1 = MS(eval(str(GG.GeneratorsOfGroup()[1]).replace("\n","")))
            sage: g2 = MS(eval(str(GG.GeneratorsOfGroup()[2]).replace("\n","")))
            sage: g3 = MS(eval(str(GG.GeneratorsOfGroup()[3]).replace("\n","")))
            sage: G = MatrixGroup([g1, g2, g3])
            sage: G.cardinality()
            21499084800
            sage: set_random_seed(0); current_randstate().set_seed_gap()
            sage: P = G.as_permutation_group()
            sage: P.cardinality()
            21499084800
            sage: P.degree()  # random output
            144
            sage: set_random_seed(3); current_randstate().set_seed_gap()
            sage: Psmaller = G.as_permutation_group(algorithm="smaller")
            sage: Psmaller.cardinality()
            21499084800
            sage: Psmaller.degree()  # random output
            108

        In this case, the "smaller" option returned an isomorphic group of
        lower degree. The above example used GAP's library of irreducible
        maximal finite ("imf") integer matrix groups to construct the
        MatrixGroup G over GF(7). The section "Irreducible Maximal Finite
        Integral Matrix Groups" in the GAP reference manual has more
        details.

        TESTS::

            sage: A= matrix(QQ, 2, [0, 1, 1, 0])
            sage: B= matrix(QQ, 2, [1, 0, 0, 1])
            sage: a, b= MatrixGroup([A, B]).as_permutation_group().gens()
            sage: a.order(), b.order()
            (2, 1)
        """
        # Note that the output of IsomorphismPermGroup() depends on
        # memory locations and will change if you change the order of
        # doctests and/or architecture
        from sage.groups.perm_gps.permgroup import PermutationGroup
        if not self.is_finite():
            raise NotImplementedError("Group must be finite.")
        n = self.degree()
        MS = MatrixSpace(self.base_ring(), n, n)
        mats = [] # initializing list of mats by which the gens act on self
        for g in self.gens():
            p = MS(g.matrix())
            m = p.rows()
            mats.append(m)
        mats_str = str(gap([[list(r) for r in m] for m in mats]))
        gap.eval("iso:=IsomorphismPermGroup(Group("+mats_str+"))")
        if algorithm == "smaller":
            gap.eval("small:= SmallerDegreePermutationRepresentation( Image( iso ) );")
            C = gap("Image( small )")
        else:
            C = gap("Image( iso )")
        return PermutationGroup(gap_group=C, canonicalize=False)

    def module_composition_factors(self, algorithm=None):
        r"""
        Return a list of triples consisting of [base field, dimension,
        irreducibility], for each of the Meataxe composition factors
        modules. The ``algorithm="verbose"`` option returns more information,
        but in Meataxe notation.

        EXAMPLES::

            sage: F=GF(3);MS=MatrixSpace(F,4,4)
            sage: M=MS(0)
            sage: M[0,1]=1;M[1,2]=1;M[2,3]=1;M[3,0]=1
            sage: G = MatrixGroup([M])
            sage: G.module_composition_factors()
            [(Finite Field of size 3, 1, True),
             (Finite Field of size 3, 1, True),
             (Finite Field of size 3, 2, True)]
            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[0,1],[-1,0]]),MS([[1,1],[2,3]])]
            sage: G = MatrixGroup(gens)
            sage: G.module_composition_factors()
            [(Finite Field of size 7, 2, True)]

        Type ``G.module_composition_factors(algorithm='verbose')`` to get a
        more verbose version.

        For more on MeatAxe notation, see
        http://www.gap-system.org/Manuals/doc/ref/chap69.html
        """
        from sage.misc.sage_eval import sage_eval
        F = self.base_ring()
        if not(F.is_finite()):
            raise NotImplementedError("Base ring must be finite.")
        q = F.cardinality()
        gens = self.gens()
        n = self.degree()
        MS = MatrixSpace(F,n,n)
        mats = [] # initializing list of mats by which the gens act on self
        W = self.matrix_space().row_space()
        for g in gens:
            p = MS(g.matrix())
            m = p.rows()
            mats.append(m)
        mats_str = str(gap([[list(r) for r in m] for m in mats]))
        gap.eval("M:=GModuleByMats("+mats_str+", GF("+str(q)+"))")
        gap.eval("MCFs := MTX.CompositionFactors( M )")
        N = eval(gap.eval("Length(MCFs)"))
        if algorithm == "verbose":
            print(gap.eval('MCFs') + "\n")
        L = []
        for i in range(1,N+1):
            gap.eval("MCF := MCFs[%s]"%i)
            L.append(tuple([sage_eval(gap.eval("MCF.field")),
                            eval(gap.eval("MCF.dimension")),
                            sage_eval(gap.eval("MCF.IsIrreducible")) ]))
        return sorted(L)

    def invariant_generators(self):
        r"""
        Return invariant ring generators.

        Computes generators for the polynomial ring
        `F[x_1,\ldots,x_n]^G`, where `G` in `GL(n,F)` is a finite matrix
        group.

        In the "good characteristic" case the polynomials returned
        form a minimal generating set for the algebra of `G`-invariant
        polynomials.  In the "bad" case, the polynomials returned
        are primary and secondary invariants, forming a not
        necessarily minimal generating set for the algebra of
        `G`-invariant polynomials.

        ALGORITHM:

        Wraps Singular's ``invariant_algebra_reynolds`` and ``invariant_ring``
        in ``finvar.lib``.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[0,1],[-1,0]]),MS([[1,1],[2,3]])]
            sage: G = MatrixGroup(gens)
            sage: G.invariant_generators()
            [x1^7*x2 - x1*x2^7, 
             x1^12 - 2*x1^9*x2^3 - x1^6*x2^6 + 2*x1^3*x2^9 + x2^12, 
             x1^18 + 2*x1^15*x2^3 + 3*x1^12*x2^6 + 3*x1^6*x2^12 - 2*x1^3*x2^15 + x2^18]

            sage: q = 4; a = 2
            sage: MS = MatrixSpace(QQ, 2, 2)
            sage: gen1 = [[1/a,(q-1)/a],[1/a, -1/a]]; gen2 = [[1,0],[0,-1]]; gen3 = [[-1,0],[0,1]]
            sage: G = MatrixGroup([MS(gen1),MS(gen2),MS(gen3)])
            sage: G.cardinality()
            12
            sage: G.invariant_generators()
            [x1^2 + 3*x2^2, x1^6 + 15*x1^4*x2^2 + 15*x1^2*x2^4 + 33*x2^6]

            sage: F = CyclotomicField(8)
            sage: z = F.gen()
            sage: a = z+1/z
            sage: b = z^2
            sage: MS = MatrixSpace(F,2,2)
            sage: g1 = MS([[1/a, 1/a], [1/a, -1/a]])
            sage: g2 = MS([[-b, 0], [0, b]])
            sage: G=MatrixGroup([g1,g2])
            sage: G.invariant_generators()
            [x1^4 + 2*x1^2*x2^2 + x2^4,
             x1^5*x2 - x1*x2^5,
             x1^8 + 28/9*x1^6*x2^2 + 70/9*x1^4*x2^4 + 28/9*x1^2*x2^6 + x2^8]

        AUTHORS:

        - David Joyner, Simon King and Martin Albrecht.

        REFERENCES:

        - Singular reference manual

        - [Stu1993]_

        - S. King, "Minimal Generating Sets of non-modular invariant
          rings of finite groups", :arxiv:`math/0703035`.
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.interfaces.singular import singular
        gens = self.gens()
        singular.LIB("finvar.lib")
        n = self.degree() #len((gens[0].matrix()).rows())
        F = self.base_ring()
        q = F.characteristic()
        ## test if the field is admissible
        if F.gen()==1: # we got the rationals or GF(prime)
            FieldStr = str(F.characteristic())
        elif hasattr(F,'polynomial'): # we got an algebraic extension
            if len(F.gens())>1:
                raise NotImplementedError("can only deal with finite fields and (simple algebraic extensions of) the rationals")
            FieldStr = '(%d,%s)'%(F.characteristic(),str(F.gen()))
        else: # we have a transcendental extension
            FieldStr = '(%d,%s)'%(F.characteristic(),','.join([str(p) for p in F.gens()]))

        ## Setting Singular's variable names
        ## We need to make sure that field generator and variables get different names.
        if str(F.gen())[0]=='x':
            VarStr = 'y'
        else:
            VarStr = 'x'
        VarNames='('+','.join((VarStr+str(i+1) for i in range(n)))+')'
        R=singular.ring(FieldStr,VarNames,'dp')
        if hasattr(F,'polynomial') and F.gen()!=1: # we have to define minpoly
            singular.eval('minpoly = '+str(F.polynomial()).replace('x',str(F.gen())))
        A = [singular.matrix(n,n,str((x.matrix()).list())) for x in gens]
        Lgens = ','.join((x.name() for x in A))
        PR = PolynomialRing(F,n,[VarStr+str(i) for i in range(1,n+1)])

        if q == 0 or (q > 0 and self.cardinality()%q != 0):
            from sage.all import Integer, Matrix
            try:
                elements = [ g.matrix() for g in self.list() ]
            except (TypeError,ValueError):
                elements
            if elements is not None:
                ReyName = 't'+singular._next_var_name()
                singular.eval('matrix %s[%d][%d]'%(ReyName,self.cardinality(),n))
                for i in range(1,self.cardinality()+1):
                    M = Matrix(F, elements[i-1])
                    D = [{} for foobar in range(self.degree())]
                    for x,y in M.dict().items():
                        D[x[0]][x[1]] = y
                    for row in range(self.degree()):
                        for t in D[row].items():
                            singular.eval('%s[%d,%d]=%s[%d,%d]+(%s)*var(%d)'
                                          %(ReyName,i,row+1,ReyName,i,row+1, repr(t[1]),t[0]+1))
                foobar = singular(ReyName)
                IRName = 't'+singular._next_var_name()
                singular.eval('matrix %s = invariant_algebra_reynolds(%s)'%(IRName,ReyName))
            else:
                ReyName = 't'+singular._next_var_name()
                singular.eval('list %s=group_reynolds((%s))'%(ReyName,Lgens))
                IRName = 't'+singular._next_var_name()
                singular.eval('matrix %s = invariant_algebra_reynolds(%s[1])'%(IRName,ReyName))

            OUT = [singular.eval(IRName+'[1,%d]'%(j))
                   for j in range(1,1+singular('ncols('+IRName+')'))]
            return [PR(gen) for gen in OUT]
        if self.cardinality()%q == 0:
            PName = 't'+singular._next_var_name()
            SName = 't'+singular._next_var_name()
            singular.eval('matrix %s,%s=invariant_ring(%s)'%(PName,SName,Lgens))
            OUT = [
                singular.eval(PName+'[1,%d]'%(j))
                for j in range(1,1+singular('ncols('+PName+')'))
                ] + [
                singular.eval(SName+'[1,%d]'%(j)) for j in range(2,1+singular('ncols('+SName+')'))
                ]
            return [PR(gen) for gen in OUT]

    def molien_series(self, xi=None, return_series=True, prec=20, variable='t'):
        r"""
        Compute the Molien series of this finite group with respect to the
        character ``xi``. It can be returned either as a rational function
        in one variable or a power series in one variable. The base field
        must be a finite field, the rationals, or a cyclotomic field.

        Note that the base field characteristic cannot divide the group
        order (i.e., the non-modular case).

        ALGORITHM:

        For a finite group `G` in characteristic zero we construct the Molien series as

        .. MATH::

            \frac{1}{|G|}\sum_{g \in G} \frac{\xi(g)}{\text{det}(I-tg)},

        where `I` is the indentity matrix and `t` an indeterminant.

        For characteristic `p` not dividing the order of `G`, let `k` be the base field
        and `N` the order of `G`. Define `\lambda` as a primitive `N`-th root of unity over `k`
        and `\omega` as a primitive `N`-th root of unity over `\QQ`. For each `g \in G`
        define `k_i(g)` to be the positive integer such that
        `e_i = \lambda^{k_i(g)}` for each eigenvalue `e_i` of `g`. Then the Molien series
        is computed as

        .. MATH::

            \frac{1}{|G|}\sum_{g \in G} \frac{\xi(g)}{\prod_{i=1}^n(1 - t\omega^{k_i(g)})},

        where `t` is an indeterminant. [Dec1998]_

        INPUT:

        - ``xi`` -- (default: trivial character) a linear group character of this group

        - ``return_series`` -- boolean (default: ``True``) if ``True``, then returns
          the Molien series as a power series, ``False`` as a rational function

        - ``prec`` -- integer (default: 20); power series default precision

        - ``variable`` -- string (default: ``'t'``); Variable name for the Molien series

        OUTPUT: single variable rational function or power series with integer cofficients

        EXAMPLES::

            sage: MatrixGroup(matrix(QQ,2,2,[1,1,0,1])).molien_series()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for finite groups
            sage: MatrixGroup(matrix(GF(3),2,2,[1,1,0,1])).molien_series()
            Traceback (most recent call last):
            ...
            NotImplementedError: characteristic cannot divide group order

        Tetrahedral Group::

            sage: K.<i> = CyclotomicField(4)
            sage: Tetra =  MatrixGroup([(-1+i)/2,(-1+i)/2, (1+i)/2,(-1-i)/2], [0,i, -i,0])
            sage: Tetra.molien_series(prec=30)
            1 + t^8 + 2*t^12 + t^16 + 2*t^20 + 3*t^24 + 2*t^28 + O(t^30)
            sage: mol = Tetra.molien_series(return_series=False); mol
            (t^8 - t^4 + 1)/(t^16 - t^12 - t^4 + 1)
            sage: mol.parent()
            Fraction Field of Univariate Polynomial Ring in t over Integer Ring
            sage: xi = Tetra.character(Tetra.character_table()[1])
            sage: Tetra.molien_series(xi, prec=30, variable='u')
            u^6 + u^14 + 2*u^18 + u^22 + 2*u^26 + 3*u^30 + 2*u^34 + O(u^36)
            sage: xi = Tetra.character(Tetra.character_table()[2])
            sage: Tetra.molien_series(xi)
            t^10 + t^14 + t^18 + 2*t^22 + 2*t^26 + O(t^30)

        ::

            sage: S3 = MatrixGroup(SymmetricGroup(3))
            sage: mol = S3.molien_series(prec=10); mol
            1 + t + 2*t^2 + 3*t^3 + 4*t^4 + 5*t^5 + 7*t^6 + 8*t^7 + 10*t^8 + 12*t^9 + O(t^10)
            sage: mol.parent()
            Power Series Ring in t over Integer Ring

        Octahedral Group::

            sage: K.<v> = CyclotomicField(8)
            sage: a = v-v^3 #sqrt(2)
            sage: i = v^2
            sage: Octa = MatrixGroup([(-1+i)/2,(-1+i)/2, (1+i)/2,(-1-i)/2], [(1+i)/a,0, 0,(1-i)/a])
            sage: Octa.molien_series(prec=30)
            1 + t^8 + t^12 + t^16 + t^18 + t^20 + 2*t^24 + t^26 + t^28 + O(t^30)

        Icosahedral Group::

            sage: K.<v> = CyclotomicField(10)
            sage: z5 = v^2
            sage: i = z5^5
            sage: a = 2*z5^3 + 2*z5^2 + 1 #sqrt(5)
            sage: Ico = MatrixGroup([[z5^3,0, 0,z5^2], [0,1, -1,0], [(z5^4-z5)/a, (z5^2-z5^3)/a, (z5^2-z5^3)/a, -(z5^4-z5)/a]])
            sage: Ico.molien_series(prec=40)
            1 + t^12 + t^20 + t^24 + t^30 + t^32 + t^36 + O(t^40)

        ::

            sage: G = MatrixGroup(CyclicPermutationGroup(3))
            sage: xi = G.character(G.character_table()[1])
            sage: G.molien_series(xi, prec=10)
            t + 2*t^2 + 3*t^3 + 5*t^4 + 7*t^5 + 9*t^6 + 12*t^7 + 15*t^8 + 18*t^9 + 22*t^10 + O(t^11)

        ::

            sage: K = GF(5)
            sage: S = MatrixGroup(SymmetricGroup(4))
            sage: G = MatrixGroup([matrix(K,4,4,[K(y) for u in m.list() for y in u])for m in S.gens()])
            sage: G.molien_series(return_series=False)
            1/(t^10 - t^9 - t^8 + 2*t^5 - t^2 - t + 1)

        ::

            sage: i = GF(7)(3)
            sage: G = MatrixGroup([[i^3,0,0,-i^3],[i^2,0,0,-i^2]])
            sage: xi = G.character(G.character_table()[4])
            sage: G.molien_series(xi)
            3*t^5 + 6*t^11 + 9*t^17 + 12*t^23 + O(t^25)
        """
        if not self.is_finite():
            raise NotImplementedError("only implemented for finite groups")
        if xi is None:
            xi = self.trivial_character()
        M = self.matrix_space()
        R = FractionField(self.base_ring())
        N = self.order()
        if R.characteristic() == 0:
            P = PolynomialRing(R, variable)
            t = P.gen()
            #it is possible the character is over a larger cyclotomic field
            K = xi.values()[0].parent()
            if K.degree() != 1:
                if R.degree() != 1:
                    L = K.composite_fields(R)[0]
                else:
                    L = K
            else:
                L = R
            mol = P(0)
            for g in self:
                mol += L(xi(g)) / (M.identity_matrix()-t*g.matrix()).det().change_ring(L)
        elif R.characteristic().divides(N):
            raise NotImplementedError("characteristic cannot divide group order")
        else: #char p>0
            #find primitive Nth roots of unity over base ring and QQ
            F = cyclotomic_polynomial(N).change_ring(R)
            w = F.roots(ring=R.algebraic_closure(), multiplicities=False)[0]
            #don't need to extend further in this case since the order of
            #the roots of unity in the character divide the order of the group
            L = CyclotomicField(N, 'v')
            v = L.gen()
            #construct Molien series
            P = PolynomialRing(L, variable)
            t = P.gen()
            mol = P(0)
            for g in self:
                #construct Phi
                phi = L(xi(g))
                for e in g.matrix().eigenvalues():
                    #find power such that w**n  = e
                    n = 1
                    while w**n != e and n < N+1:
                        n += 1
                    #raise v to that power
                    phi *= (1-t*v**n)
                mol += P(1)/phi
        #We know the coefficients will be integers
        mol = mol.numerator().change_ring(ZZ) / mol.denominator().change_ring(ZZ)
        #divide by group order
        mol /= N
        if return_series:
            PS = PowerSeriesRing(ZZ, variable, default_prec=prec)
            return PS(mol)
        return mol

    def reynolds_operator(self, poly, xi = None):
        r"""
        Compute the Reynolds Operator of this finite group `G`. This is th eprojection from the polynomial ring to the ring of relative invariants. [Stu1993]_.

        ALGORITHM:

        Let `K[x]` be a polynomial ring. Let `K[x]^G= \{f\in K[x]  |  \pi f=f \forall \pi\in G\}` be the ring of invarants of `G`.
            Then the Reynold's operator is a map `R`, from `K[x]` into itself for the nonmodular case,
            defined by `f \rightarrow \frac{1}{|G|} \sum_{ \pi \in G} \pi f`.

        INPUT:

        -`poly`-- a polynomial from K[x]

        - ``xi`` -- (default: trivial character) a linear group character of this group

        OUTPUT:

        - the invariant polynomial under K[x]^G

        AUTHORS:
            Rebecca Lauren Miller and Ben Hutz

        EXAMPLES::

            sage: S3 = MatrixGroup(SymmetricGroup(3))
            sage: R.<x,y,z> = QQ[]
            sage: f = x*y*z^3
            sage: S3.reynolds_operator(f)
            1/3*x^3*y*z + 1/3*x*y^3*z + 1/3*x*y*z^3

        ::

            sage: G = MatrixGroup(CyclicPermutationGroup(4))
            sage: X = G.character_table()
            sage: xi = G.character(X[2])
            sage: W.<v> = CyclotomicField(4)
            sage: J.<x,y,z,w> = W[]
            sage: f = x
            sage: G.reynolds_operator(f,xi)
            1/4*x + (-1/4*v)*y - 1/4*z + (1/4*v)*w

        ::

            sage: G = MatrixGroup(CyclicPermutationGroup(3))
            sage: xi = G.character(G.character_table()[1])
            sage: R.<x,y,z> = QQ[]
            sage: f = x*y
            sage: G.reynolds_operator(f,xi)
            1/3*x*y + (-1/3*zeta3 - 1/3)*x*z + (1/3*zeta3)*y*z

        ::

            sage: i = GF(7)(3)
            sage: G = MatrixGroup([[i^3,0,0,-i^3],[i^2,0,0,-i^2]])
            sage: xi = G.character(G.character_table()[4])
            sage: R.<w,x>=GF(7)[]
            sage: f=w^5*x + x^6
            sage: G.reynolds_operator(f,xi)
            Traceback (most recent call last):
            ...
            NotImplementedError: Not implemented for this base field's charateristic
            sage: G.reynolds_operator(f)
            x^6

        ::

            sage: K.<i> = CyclotomicField(4)
            sage: Tetra =  MatrixGroup([(-1+i)/2,(-1+i)/2, (1+i)/2,(-1-i)/2], [0,i, -i,0])
            sage: xi = Tetra.character(Tetra.character_table()[1])
            sage: R.<x,y>=K[]
            sage: print Tetra.reynolds_operator(x*y^5,xi)
            -1/2*x^5*y + 1/2*x*y^5
            sage: xi = Tetra.character(Tetra.character_table()[2])
            sage: Tetra.reynolds_operator(x*y^9,xi)
            -1/4*x^9*y + (2*izeta3^3 + 3*izeta3^2 + 8*izeta3 + 7/2)*x^7*y^3
            + (-2*izeta3^3 - 3*izeta3^2 - 8*izeta3 - 7/2)*x^3*y^7 + 1/4*x*y^9
        """
        f=poly
        G=self
        F=0
        R = FractionField(f.base_ring())
        if xi is None:
            if R.characteristic()==0 or (not R.characteristic().divides(G.order())):
                F=0
                for g in G:
                    F += f(*g.matrix()*vector(f.parent().gens()))
                F = F/G.order()
                return F
        K = xi.values()[0].parent()
        if R.characteristic()==0:
            if K!=R: # extend base ring of G to include all xi values
                if K.degree() != 1:
                    if R.degree() != 1:
                        L = R.composite_fields(K)[0]
                    else:
                        L = K
                else:
                    L = R
            else:
                L = R
        else:
            raise NotImplementedError("Not implemented for this base field's charateristic")
        f=f.change_ring(L)
        F=L(0)
        for g in G:
            F += L(xi(g))*f(*g.matrix().change_ring(L)*vector(f.parent().gens()))
        F = F/G.order()
        return  F
