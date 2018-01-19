r"""
Orthogonal Linear Groups

The general orthogonal group `GO(n,R)` consists of all `n\times n`
matrices over the ring `R` preserving an `n`-ary positive definite
quadratic form. In cases where there are multiple non-isomorphic
quadratic forms, additional data needs to be specified to
disambiguate. The special orthogonal group is the normal subgroup of
matrices of determinant one.

In characteristics different from 2, a quadratic form is equivalent to
a bilinear symmetric form. Furthermore, over the real numbers a
positive definite quadratic form is equivalent to the diagonal
quadratic form, equivalent to the bilinear symmetric form defined by
the identity matrix. Hence, the orthogonal group `GO(n,\RR)` is the
group of orthogonal matrices in the usual sense.

In the case of a finite field and if the degree `n` is even, then
there are two inequivalent quadratic forms and a third parameter ``e``
must be specified to disambiguate these two possibilities. The index
of `SO(e,d,q)` in `GO(e,d,q)` is `2` if `q` is odd, but `SO(e,d,q) =
GO(e,d,q)` if `q` is even.)

.. warning::

   GAP and Sage use different notations:

   * GAP notation: The optional ``e`` comes first, that is, ``GO([e,]
     d, q)``, ``SO([e,] d, q)``.

   * Sage notation: The optional ``e`` comes last, the standard Python
     convention: ``GO(d, GF(q), e=0)``, ``SO( d, GF(q), e=0)``.

EXAMPLES::

    sage: GO(3,7)
    General Orthogonal Group of degree 3 over Finite Field of size 7

    sage: G = SO( 4, GF(7), 1); G
    Special Orthogonal Group of degree 4 and form parameter 1 over Finite Field of size 7
    sage: G.random_element()   # random
    [4 3 5 2]
    [6 6 4 0]
    [0 4 6 0]
    [4 4 5 1]

TESTS::

    sage: G = GO(3, GF(5))
    sage: latex(G)
    \text{GO}_{3}(\Bold{F}_{5})
    sage: G = SO(3, GF(5))
    sage: latex(G)
    \text{SO}_{3}(\Bold{F}_{5})
    sage: G = SO(4, GF(5), 1)
    sage: latex(G)
    \text{SO}_{4}(\Bold{F}_{5}, +)

AUTHORS:

- David Joyner (2006-03): initial version

- David Joyner (2006-05): added examples, _latex_, __str__, gens,
  as_matrix_group

- William Stein (2006-12-09): rewrite

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.

- Simon Brandhorst (2017-9) added OrthogonalMatrixGroup_with_gap
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2017 Simon Brandhorst <sbrandhorst@web.de>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, NamedMatrixGroup_generic, NamedMatrixGroup_gap )
from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap
from sage.categories.action import Action

def normalize_args_e(degree, ring, e):
    """
    Normalize the arguments that relate the choice of quadratic form
    for special orthogonal groups over finite fields.

    INPUT:

    - ``degree`` -- integer. The degree of the affine group, that is,
      the dimension of the affine space the group is acting on.

    - ``ring`` -- a ring. The base ring of the affine space.

    - ``e`` -- integer, one of `+1`, `0`, `-1`.  Only relevant for
      finite fields and if the degree is even. A parameter that
      distinguishes inequivalent invariant forms.

    OUTPUT:

    The integer ``e`` with values required by GAP.

    TESTS::

        sage: from sage.groups.matrix_gps.orthogonal import normalize_args_e
        sage: normalize_args_e(2, GF(3), +1)
        1
        sage: normalize_args_e(3, GF(3), 0)
        0
        sage: normalize_args_e(3, GF(3), +1)
        0
        sage: normalize_args_e(2, GF(3), 0)
        Traceback (most recent call last):
        ...
        ValueError: must have e=-1 or e=1 for even degree
    """
    if is_FiniteField(ring) and degree%2 == 0:
        if e not in (-1, +1):
            raise ValueError('must have e=-1 or e=1 for even degree')
    else:
        e = 0
    return ZZ(e)





########################################################################
# General Orthogonal Group
########################################################################

def GO(n, R, e=0, var='a'):
    """
    Return the general orthogonal group.

    The general orthogonal group `GO(n,R)` consists of all `n\times n`
    matrices over the ring `R` preserving an `n`-ary positive definite
    quadratic form. In cases where there are multiple non-isomorphic
    quadratic forms, additional data needs to be specified to
    disambiguate.

    In the case of a finite field and if the degree `n` is even, then
    there are two inequivalent quadratic forms and a third parameter
    ``e`` must be specified to disambiguate these two possibilities.

    .. note::

        This group is also available via ``groups.matrix.GO()``.

   INPUT:

    - ``n`` -- integer. The degree.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``e`` -- ``+1`` or ``-1``, and ignored by default. Only relevant
      for finite fields and if the degree is even. A parameter that
      distinguishes inequivalent invariant forms.

    OUTPUT:

    The general orthogonal group of given degree, base ring, and
    choice of invariant form.

    EXAMPLES::

        sage: GO( 3, GF(7))
        General Orthogonal Group of degree 3 over Finite Field of size 7
        sage: GO( 3, GF(7)).order()
        672
        sage: GO( 3, GF(7)).gens()
        (
        [3 0 0]  [0 1 0]
        [0 5 0]  [1 6 6]
        [0 0 1], [0 2 1]
        )

    TESTS::

        sage: groups.matrix.GO(2, 3, e=-1)
        General Orthogonal Group of degree 2 and form parameter -1 over Finite Field of size 3
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    e = normalize_args_e(degree, ring, e)
    if e == 0:
        name = 'General Orthogonal Group of degree {0} over {1}'.format(degree, ring)
        ltx  = r'\text{{GO}}_{{{0}}}({1})'.format(degree, latex(ring))
    else:
        name = 'General Orthogonal Group of degree' + \
            ' {0} and form parameter {1} over {2}'.format(degree, e, ring)
        ltx  = r'\text{{GO}}_{{{0}}}({1}, {2})'.format(degree, latex(ring), '+' if e == 1 else '-')
    if is_FiniteField(ring):
        cmd  = 'GO({0}, {1}, {2})'.format(e, degree, ring.characteristic())
        return OrthogonalMatrixGroup_gap(degree, ring, False, name, ltx, cmd)
    else:
        return OrthogonalMatrixGroup_generic(degree, ring, False, name, ltx)



########################################################################
# Special Orthogonal Group
########################################################################

def SO(n, R, e=None, var='a'):
    """
    Return the special orthogonal group.

    The special orthogonal group `GO(n,R)` consists of all `n\times n`
    matrices with determinant one over the ring `R` preserving an
    `n`-ary positive definite quadratic form. In cases where there are
    multiple non-isomorphic quadratic forms, additional data needs to
    be specified to disambiguate.

    .. note::

        This group is also available via ``groups.matrix.SO()``.

    INPUT:

    - ``n`` -- integer. The degree.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``e`` -- ``+1`` or ``-1``, and ignored by default. Only relevant
      for finite fields and if the degree is even. A parameter that
      distinguishes inequivalent invariant forms.

    OUTPUT:

    The special orthogonal group of given degree, base ring, and choice of
    invariant form.

    EXAMPLES::

        sage: G = SO(3,GF(5))
        sage: G
        Special Orthogonal Group of degree 3 over Finite Field of size 5

        sage: G = SO(3,GF(5))
        sage: G.gens()
        (
        [2 0 0]  [3 2 3]  [1 4 4]
        [0 3 0]  [0 2 0]  [4 0 0]
        [0 0 1], [0 3 1], [2 0 4]
        )
        sage: G = SO(3,GF(5))
        sage: G.as_matrix_group()
        Matrix group over Finite Field of size 5 with 3 generators (
        [2 0 0]  [3 2 3]  [1 4 4]
        [0 3 0]  [0 2 0]  [4 0 0]
        [0 0 1], [0 3 1], [2 0 4]
        )

    TESTS::

        sage: groups.matrix.SO(2, 3, e=1)
        Special Orthogonal Group of degree 2 and form parameter 1 over Finite Field of size 3
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    e = normalize_args_e(degree, ring, e)
    if e == 0:
        name = 'Special Orthogonal Group of degree {0} over {1}'.format(degree, ring)
        ltx  = r'\text{{SO}}_{{{0}}}({1})'.format(degree, latex(ring))
    else:
        name = 'Special Orthogonal Group of degree' + \
            ' {0} and form parameter {1} over {2}'.format(degree, e, ring)
        ltx  = r'\text{{SO}}_{{{0}}}({1}, {2})'.format(degree, latex(ring), '+' if e == 1 else '-')
    if is_FiniteField(ring):
        cmd  = 'SO({0}, {1}, {2})'.format(e, degree, ring.characteristic())
        return OrthogonalMatrixGroup_gap(degree, ring, True, name, ltx, cmd)
    else:
        return OrthogonalMatrixGroup_generic(degree, ring, True, name, ltx)



########################################################################
# Orthogonal Group class
########################################################################

class OrthogonalMatrixGroup_generic(NamedMatrixGroup_generic):

    @cached_method
    def invariant_bilinear_form(self):
        """
        Return the symmetric bilinear form preserved by the orthogonal
        group.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: GO(2,3,+1).invariant_bilinear_form()
            [0 1]
            [1 0]
            sage: GO(2,3,-1).invariant_bilinear_form()
            [2 1]
            [1 1]
        """
        from sage.matrix.constructor import identity_matrix
        m = identity_matrix(self.base_ring(), self.degree())
        m.set_immutable()
        return m

    @cached_method
    def invariant_quadratic_form(self):
        """
        Return the quadratic form preserved by the orthogonal group.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: GO(2,3,+1).invariant_quadratic_form()
            [0 1]
            [0 0]
            sage: GO(2,3,-1).invariant_quadratic_form()
            [1 1]
            [0 2]
        """
        from sage.matrix.constructor import identity_matrix
        m = identity_matrix(self.base_ring(), self.degree())
        m.set_immutable()
        return m

    def _check_matrix(self, x, *args):
        """a
        Check whether the matrix ``x`` is symplectic.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = GO(4, GF(5), +1)
            sage: G._check_matrix(G.an_element().matrix())
        """
        if self._special and x.determinant() != 1:
            raise TypeError('matrix must have determinant one')
        F = self.invariant_bilinear_form()
        if x * F * x.transpose() != F:
            raise TypeError('matrix must be orthogonal with respect to the invariant form')
        # TODO: check that quadratic form is preserved in characteristic two

class OrthogonalMatrixGroup_gap(OrthogonalMatrixGroup_generic, NamedMatrixGroup_gap):

    @cached_method
    def invariant_bilinear_form(self):
        """
        Return the symmetric bilinear form preserved by the orthogonal
        group.

        OUTPUT:

        A matrix `M` such that, for every group element g, the
        identity `g m g^T = m` holds. In characteristic different from
        two, this uniquely determines the orthogonal group.

        EXAMPLES::

            sage: G = GO(4, GF(7), -1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]

            sage: G = GO(4, GF(7), +1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 6 0]
            [0 0 0 2]

            sage: G = GO(4, QQ)
            sage: G.invariant_bilinear_form()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = SO(4, GF(7), -1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]
        """
        m = self.gap().InvariantBilinearForm()['matrix'].matrix()
        m.set_immutable()
        return m

    @cached_method
    def invariant_quadratic_form(self):
        """
        Return the quadratic form preserved by the orthogonal group.

        OUTPUT:

        The matrix `Q` defining "orthogonal" as follows. The matrix
        determines a quadratic form `q` on the natural vector space
        `V`, on which `G` acts, by `q(v) = v Q v^t`. A matrix `M' is
        an element of the orthogonal group if `q(v) = q(v M)` for all
        `v \in V`.

        EXAMPLES::

            sage: G = GO(4, GF(7), -1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = GO(4, GF(7), +1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 3 0]
            [0 0 0 1]

            sage: G = GO(4, QQ)
            sage: G.invariant_quadratic_form()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = SO(4, GF(7), -1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        m = self.gap().InvariantQuadraticForm()['matrix'].matrix()
        m.set_immutable()
        return m

class OrthogonalMatrixGroup_with_gap(FinitelyGeneratedMatrixGroup_gap):
    r"""
    A base class for Orthogonal matrix groups with a gap backend.

    It remembers the bilinear form.
    The difference to `OrthogonalMatrixGroup_gap` is that our groups do not have
    a specific name.

    INPUT:

    - ``degree`` -- integer, the degree (matrix size) of the matrix
    - ``base_ring`` -- ring, the base ring of the matrices
    - ``gens`` -- a list of matrices over the base ring
    - ``invariant_bilinear_form`` -- a symmetric matrix
    - ``category`` -- (default:``None``) a category of groups
    - ``check`` -- bool (default: ``True``) check if the generators
      preserve the bilinear form
    - ``invariant_submodule`` -- a submodule preserved by the group action
      (default: ``None``) registers an action on this submodule.
    - ``invariant_quotient_module`` -- a quotient module preserved by
      the group action (default: ``None``)
      registers an action on this quotient module.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
        sage: bil = Matrix(ZZ,2,[3,2,2,3])
        sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
        sage: O = OrthogonalMatrixGroup_with_gap(2,ZZ,gens,bil)
        sage: O
        Orthogonal group over Integer Ring with 1 generators (
        [ 0 -1]
        [-1  0]
        )
        sage: O.order()
        2

    Infinite groups are O.K. too::

        sage: bil = Matrix(ZZ,4,[0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0])
        sage: f = Matrix(ZZ,4,[0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, -1, 1, 1, 1])
        sage: O = OrthogonalMatrixGroup_with_gap(2,ZZ,[f],bil)
        sage: O.cardinality()
        +Infinity
    """

    def __init__(self, degree, base_ring,
                 gens, invariant_bilinear_form,
                 category=None, check=True,
                 invariant_submodule=None,
                 invariant_quotient_module=None):
        r"""
        Create this orthogonal group from the input.

        TESTS::

            sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: O = OrthogonalMatrixGroup_with_gap(2,ZZ,gens,bil)
            sage: TestSuite(O).run()
        """
        from copy import copy
        G = copy(invariant_bilinear_form)
        G.set_immutable()
        self._invariant_bilinear_form = G
        self._invariant_submodule = invariant_submodule
        self._invariant_quotient_module = invariant_quotient_module
        if check:
            for f in gens:
                self._check_matrix(f)
        if len(gens) == 0:    # handle the trivial group
            gens = [G.parent().identity_matrix()]
        from sage.libs.gap.libgap import libgap
        gap_gens = [libgap(matrix_gen) for matrix_gen in gens]
        gap_group = libgap.Group(gap_gens)
        FinitelyGeneratedMatrixGroup_gap.__init__(self,
                                                  degree,
                                                  base_ring,
                                                  gap_group,
                                                  category=category)

    def _repr_(self):
        r"""
        Return the string representation of this matrix group.

        OUTPUT:

        - a string

        EXAMPLES::

            sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: O = OrthogonalMatrixGroup_with_gap(2,ZZ,gens,bil)
            sage: O
            Orthogonal group over Integer Ring with 1 generators (
            [ 0 -1]
            [-1  0]
            )
        """
        if self.ngens() > 5:
            return 'Orthogonal group over {0} with {1} generators'.format(
                self.base_ring(), self.ngens())
        else:
            from sage.repl.display.util import format_list
            return 'Orthogonal group over {0} with {1} generators {2}'.format(
                self.base_ring(), self.ngens(), format_list(self.gens()))

    def invariant_bilinear_form(self):
        r"""
        Return the symmetric bilinear form preserved by the orthogonal group.

        OUTPUT:

        - the matrix defining the bilinear form

        EXAMPLES::

            sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: O = OrthogonalMatrixGroup_with_gap(2,ZZ,gens,bil)
            sage: O.invariant_bilinear_form()
            [3 2]
            [2 3]
        """
        return self._invariant_bilinear_form

    def _get_action_(self, S, op, self_on_left):
        """
        Provide the coercion system with an action.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: S = ZZ^2
            sage: O = OrthogonalMatrixGroup_with_gap(2, ZZ, gens, bil, invariant_submodule=S)
            sage: O._get_action_(S, operator.mul, False)
            Right action by Orthogonal group over Integer Ring with 1 generators (
            [ 0 -1]
            [-1  0]
            ) on Ambient free module of rank 2 over the principal ideal domain Integer Ring
        """
        import operator
        if (S is self._invariant_submodule and op == operator.mul
                                           and not self_on_left):
            return GroupActionOnSubmodule(self,S)
        if (S is self._invariant_quotient_module and op == operator.mul
                                                 and not self_on_left):
            return GroupActionOnQuotientModule(self,S)
        return None

    def _check_matrix(self, x, *args):
        r"""
        Check whether the matrix ``x`` preserves the bilinear form.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: O = OrthogonalMatrixGroup_with_gap(2,ZZ,gens,bil)
            sage: g = matrix.identity(2)*2
            sage: O(g)
            Traceback (most recent call last):
            ...
            TypeError: matrix must be orthogonal with respect to the invariant form
        """
        F = self.invariant_bilinear_form()
        if x * F * x.transpose() != F:
            raise TypeError('matrix must be orthogonal '
                'with respect to the invariant form')

class GroupActionOnSubmodule(Action):
    r"""
    Matrix group action on a submodule from the right.

    INPUT:

    - ``MatrixGroup`` --  an instance of :class:`OrthogonalMatrixGroup_with_gap`
    - ``submodule`` -- an invariant submodule
    - ``is_left`` -- bool (default: False)

    EXAMPLES::

        sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
        sage: S = span(ZZ,[[0,1]])
        sage: g = Matrix(QQ,2,[1,0,0,-1])
        sage: G = OrthogonalMatrixGroup_with_gap(2, ZZ, [g], invariant_bilinear_form=matrix.identity(2), invariant_submodule=S)
        sage: g = G.an_element()
        sage: x = S.an_element()
        sage: x*g
        (0, -1)
        sage: (x*g).parent()
        Free module of degree 2 and rank 1 over Integer Ring
        Echelon basis matrix:
        [0 1]
    """
    def __init__(self, MatrixGroup,submodule, is_left=False):
        r"""
        Initialize the action

        TESTS::

            sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap, GroupActionOnSubmodule
            sage: S = span(ZZ,[[0,1]])
            sage: g = Matrix(QQ,2,[1,0,0,-1])
            sage: e = Matrix.identity(2)
            sage: G = OrthogonalMatrixGroup_with_gap(2, ZZ, [g], e)
            sage: GroupActionOnSubmodule(G,S)
            Right action by Orthogonal group over Integer Ring with 1 generators (
            [ 1  0]
            [ 0 -1]
            ) on Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1]
        """
        import operator
        Action.__init__(self, MatrixGroup, submodule, is_left, operator.mul)

    def _call_(self, a, g):
        r"""
        This defines the group action.

        INPUT:

        - ``a`` -- an element of the invariant submodule
        - ``g`` -- an element of the acting group

        OUTPUT:

        - an element of the invariant submodule

        EXAMPLES::

            sage: from sage.groups.matrix_gps.orthogonal import GroupActionOnSubmodule
            sage: S = span(QQ,[[0,1]])
            sage: g = Matrix(QQ,2,[1,1,0,1/2])
            sage: G = MatrixGroup([g])
            sage: A = GroupActionOnSubmodule(G,S)
            sage: A
            Right action by Matrix group over Rational Field with 1 generators (
            [  1   1]
            [  0 1/2]
            ) on Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
            sage: s = S.an_element()
            sage: g = G.an_element()
            sage: A(s,g).parent()
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
        """
        return(a.parent()(a*g.matrix()))

class GroupActionOnQuotientModule(Action):
    r"""
    Matrix group action on a quotient module from the right.

    INPUT:

    - ``MatrixGroup`` --  the group acting
          :class:`OrthogonalMatrixGroup_with_gap`
    - ``submodule`` -- an invariant submodule
    - ``is_left`` -- bool (default: False)

    EXAMPLES::

        sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
        sage: S = span(ZZ,[[0,1]])
        sage: Q = S/(6*S)
        sage: g = Matrix(QQ,2,[1,0,0,-1])
        sage: G = OrthogonalMatrixGroup_with_gap(2, ZZ, [g], invariant_bilinear_form=matrix.identity(2), invariant_quotient_module=Q)
        sage: g = G.an_element()
        sage: x = Q.an_element()
        sage: x*g
        (5)
        sage: (x*g).parent()
        Finitely generated module V/W over Integer Ring with invariants (6)
    """
    def __init__(self, MatrixGroup, quotient_module, is_left=False):
        r"""
        Initialize the action

        TESTS::

            sage: from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_with_gap
            sage: S = span(ZZ,[[0,1]])
            sage: Q = S/(6*S)
            sage: g = Matrix(QQ,2,[1,0,0,-1])
            sage: G = OrthogonalMatrixGroup_with_gap(2, ZZ, [g], invariant_bilinear_form=matrix.identity(2), invariant_quotient_module=Q)
            sage: g = G.an_element()
            sage: x = Q.an_element()
            sage: x, x*g
            ((1), (5))
        """
        import operator
        Action.__init__(self, MatrixGroup, quotient_module, is_left, operator.mul)

    def _call_(self, a, g):
        r"""
        This defines the group action.

        INPUT:

        - ``a`` -- an element of the invariant submodule
        - ``g`` -- an element of the acting group

        OUTPUT:

        - an element of the invariant quotient module

        EXAMPLES::

            sage: from sage.groups.matrix_gps.orthogonal import GroupActionOnQuotientModule
            sage: S = span(ZZ,[[0,1]])
            sage: Q = S/(6*S)
            sage: g = Matrix(QQ,2,[1,1,0,7])
            sage: G = MatrixGroup([g])
            sage: A = GroupActionOnQuotientModule(G,Q)
            sage: A
            Right action by Matrix group over Rational Field with 1 generators (
            [1 1]
            [0 7]
            ) on Finitely generated module V/W over Integer Ring with invariants (6)
            sage: q = Q.an_element()
            sage: g = G.an_element()
            sage: A(q,g).parent()
            Finitely generated module V/W over Integer Ring with invariants (6)
        """
        return(a.parent()(a.lift()*g.matrix()))
