r"""
Finitely generated modules over a PID

You can use Sage to compute with finitely generated modules (FGM's)
over a principal ideal domain R presented as a quotient V/W, where V
and W are free.

NOTE: Currently this is only enabled over R=ZZ, since it has not been
tested and debugged over more general PIDs.  All algorithms make sense
whenever there is a Hermite form implementation.  In theory the
obstruction to extending the implementation is only that one has to
decide how elements print.  If you're annoyed that by this, fix things
and post a patch!

We represent M=V/W as a pair (V,W) with W contained in V, and we
internally represent elements of M non-canonically as elements x of
V.  We also fix independent generators g[i] for M in V, and when we
print out elements of V we print their coordinates with respect to the
g[i]; over `\ZZ` this is canonical, since each coefficient is reduce
modulo the additive order of g[i]. To obtain the vector in V
corresponding to x in M, use x.lift().

Morphisms between finitely generated R modules are well supported.
You create a homomorphism by simply giving the images of generators of
M0 in M1.  Given a morphism phi:M0-->M1, you can compute the image of
phi, the kernel of phi, and using y=phi.lift(x) you can lift an
elements x in M1 to an element y in M0, if such a y exists.

TECHNICAL NOTE: For efficiency, we introduce a notion of optimized
representation for quotient modules.  The optimized representation of
M=V/W is the quotient V'/W' where V' has as basis lifts of the
generators g[i] for M.  We internally store a morphism from M0=V0/W0
to M1=V1/W1 by giving a morphism from the optimized representation V0'
of M0 to V1 that sends W0 into W1.



The following TUTORIAL illustrates several of the above points.

First we create free modules V0 and W0 and the quotient module M0.
Notice that everything works fine even though V0 and W0 are not
contained inside `\ZZ^n`, which is extremely convenient. ::

    sage: V0 = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W0 = V0.span([V0.0+2*V0.1, 9*V0.0+2*V0.1, 4*V0.2])
    sage: M0 = V0/W0; M0
    Finitely generated module V/W over Integer Ring with invariants (4, 16)

The invariants are computed using the Smith normal form algorithm, and
determine the structure of this finitely generated module.

You can get the V and W used in constructing the quotient module using
V() and W() methods::

    sage: M0.V()
    Free module of degree 3 and rank 3 over Integer Ring
    Echelon basis matrix:
    [1/2   0   0]
    [  0   2   0]
    [  0   0   1]
    sage: M0.W()
    Free module of degree 3 and rank 3 over Integer Ring
    Echelon basis matrix:
    [1/2   4   0]
    [  0  32   0]
    [  0   0   4]

We note that the optimized representation of M0, mentioned above in
the technical note has a V that need not be equal to V0, in general. ::

    sage: M0.optimized()[0].V()
    Free module of degree 3 and rank 2 over Integer Ring
    User basis matrix:
    [0 0 1]
    [0 2 0]

Create elements of M0 either by coercing in elements of V0, getting
generators, or coercing in a list or tuple or coercing in 0.  Coercing
in a list or tuple takes the corresponding linear combination of the
generators of M0. ::

    sage: M0(V0.0)
    (0, 14)
    sage: M0(V0.0 + W0.0)  # no difference modulo W0
    (0, 14)
    sage: M0([3,20])
    (3, 4)
    sage: 3*M0.0 + 20*M0.1
    (3, 4)

We make an element of M0 by taking a difference of two generators, and
lift it.  We also illustrate making an element from a list, which
coerces to V0, then take the equivalence class modulo W0. ::

    sage: x = M0.0 - M0.1; x
    (1, 15)
    sage: x.lift()
    (0, -2, 1)
    sage: M0(vector([1/2,0,0]))
    (0, 14)
    sage: x.additive_order()
    16


Similarly, we construct V1 and W1, and the quotient M1, in a completely different
2-dimensional ambient space. ::

    sage: V1 = span([[1/2,0],[3/2,2]],ZZ); W1 = V1.span([2*V1.0, 3*V1.1])
    sage: M1 = V1/W1; M1
    Finitely generated module V/W over Integer Ring with invariants (6)

We create the homomorphism from M0 to M1 that sends both generators of
M0 to 3 times the generator of M1.  This is well defined since 3 times
the generator has order 2. ::

    sage: f = M0.hom([3*M1.0, 3*M1.0]); f
    Morphism from module over Integer Ring with invariants (4, 16) to module with invariants (6,) that sends the generators to [(3), (3)]

We evaluate the homomorphism on our element x of the domain, and on the
first generator of the domain.  We also evaluate at an element of V0,
which is coerced into M0. ::

    sage: f(x)
    (0)
    sage: f(M0.0)
    (3)
    sage: f(V0.1)
    (3)

Here we illustrate lifting an element of the image of f, i.e., finding
an element of M0 that maps to a given element of M1::

    sage: y = f.lift(3*M1.0); y
    (0, 13)
    sage: f(y)
    (3)

We compute the kernel of f, i.e., the submodule of elements of M0 that
map to 0.  Note that the kernel is not explicitly represented as a
submodule, but as another quotient V/W where V is contained in V0.
You can explicitly coerce elements of the kernel into M0 though. ::

    sage: K = f.kernel(); K
    Finitely generated module V/W over Integer Ring with invariants (2, 16)

    sage: M0(K.0)
    (2, 0)
    sage: M0(K.1)
    (3, 1)
    sage: f(M0(K.0))
    (0)
    sage: f(M0(K.1))
    (0)

We compute the image of f. ::

    sage: f.image()
    Finitely generated module V/W over Integer Ring with invariants (2)

Notice how the elements of the image are written as (0) and (1),
despite the image being naturally a submodule of M1, which has
elements (0), (1), (2), (3), (4), (5).  However, below we coerce the
element (1) of the image into the codomain, and get (3)::

    sage: list(f.image())
    [(0), (1)]
    sage: list(M1)
    [(0), (1), (2), (3), (4), (5)]
    sage: x = f.image().0; x
    (1)
    sage: M1(x)
    (3)


TESTS::

    sage: from sage.modules.fg_pid.fgp_module import FGP_Module
    sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ)
    sage: W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
    sage: Q = FGP_Module(V, W); Q
    Finitely generated module V/W over Integer Ring with invariants (4, 12)
    sage: Q([1,3])
    (1, 3)
    sage: Q(V([1,3,4]))
    (0, 11)
    sage: Q(W([1,16,0]))
    (0, 0)
    sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],QQ)
    sage: W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1])
    sage: Q = FGP_Module(V, W); Q
    Finitely generated module V/W over Rational Field with invariants (0)
    sage: q = Q.an_element(); q
    (1)
    sage: q*(1/2)
    (1/2)
    sage: (1/2)*q
    (1/2)

AUTHOR:

    - William Stein, 2009
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

from sage.modules.module import Module
from sage.modules.free_module import is_FreeModule
from sage.structure.sequence import Sequence
from fgp_element  import DEBUG, FGP_Element
from fgp_morphism import FGP_Morphism, FGP_Homset
from sage.rings.all import Integer, ZZ, lcm
from sage.misc.cachefunc import cached_method

import sage.misc.weak_dict
from functools import reduce
_fgp_module = sage.misc.weak_dict.WeakValueDictionary()




def FGP_Module(V, W, check=True):
    """
    INPUT:

    - ``V`` -- a free R-module

    - ``W`` -- a free R-submodule of ``V``

    - ``check`` -- bool (default: ``True``); if ``True``, more checks
      on correctness are performed; in particular, we check the data
      types of ``V`` and ``W``, and that ``W`` is a submodule of ``V``
      with the same base ring.

    OUTPUT:

    - the quotient ``V/W`` as a finitely generated R-module


    EXAMPLES::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
        sage: import sage.modules.fg_pid.fgp_module
        sage: Q = sage.modules.fg_pid.fgp_module.FGP_Module(V, W)
        sage: type(Q)
        <class 'sage.modules.fg_pid.fgp_module.FGP_Module_class_with_category'>
        sage: Q is sage.modules.fg_pid.fgp_module.FGP_Module(V, W, check=False)
        True
    """
    key = (V,V.basis_matrix(),W,W.basis_matrix())
    try:
        return _fgp_module[key]
    except KeyError: pass
    M = FGP_Module_class(V, W, check=check)
    _fgp_module[key] = M
    return M


def is_FGP_Module(x):
    """
    Return true of x is an FGP module, i.e., a finitely generated
    module over a PID represented as a quotient of finitely generated
    free modules over a PID.

    EXAMPLES::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
        sage: sage.modules.fg_pid.fgp_module.is_FGP_Module(V)
        False
        sage: sage.modules.fg_pid.fgp_module.is_FGP_Module(Q)
        True
    """
    return isinstance(x, FGP_Module_class)


class FGP_Module_class(Module):
    """
    A finitely generated module over a PID presented as a quotient V/W.

    INPUT:

    - ``V`` -- an R-module

    - ``W`` -- an R-submodule of V

    - ``check`` -- bool (default: True)

    EXAMPLES::

        sage: A = (ZZ^1)/span([[100]], ZZ); A
        Finitely generated module V/W over Integer Ring with invariants (100)
        sage: A.V()
        Ambient free module of rank 1 over the principal ideal domain Integer Ring
        sage: A.W()
        Free module of degree 1 and rank 1 over Integer Ring
        Echelon basis matrix:
        [100]

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
        sage: Q = V/W; Q
        Finitely generated module V/W over Integer Ring with invariants (4, 12)
        sage: type(Q)
        <class 'sage.modules.fg_pid.fgp_module.FGP_Module_class_with_category'>
    """

    # The class to be used for creating elements of this
    # module. Should be overridden in derived classes.
    Element = FGP_Element

    def __init__(self, V, W, check=True):
        """
        INPUT:

        - ``V`` -- an R-module

        - ``W`` -- an R-submodule of V

        - ``check`` -- bool (default: True); if True, more checks on
                   correctness are performed; in particular, we check
                   the data types of V and W, and that W is a
                   submodule of V with the same base ring.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: type(Q)
            <class 'sage.modules.fg_pid.fgp_module.FGP_Module_class_with_category'>
        """
        if check:
            if not is_FreeModule(V):
                raise TypeError("V must be a FreeModule")
            if not is_FreeModule(W):
                raise TypeError("W must be a FreeModule")
            if not W.is_submodule(V):
                raise ValueError("W must be a submodule of V")
            if V.base_ring() != W.base_ring():
                raise ValueError("W and V must have the same base ring")
        self._W = W
        self._V = V
        Module.__init__(self, base=V.base_ring())

    # Note: There once was a
    # def _subquotient_class():
    # method that returned a functionoid to construct new modules, so
    # you would call module._subquotient_class()(V,W,check). This has
    # been replaced with module._module_constructor(V,W,check).

    def _module_constructor(self, V, W, check=True):
        r"""
        Construct a quotient module ``V/W``.

        This should be overridden in derived classes.

        INPUT:

        - ``V`` -- an R-module.

        - ``W`` -- an R-submodule of ``V``.

        - ``check`` -- bool (default: True).

        OUTPUT:

        The quotient ``V/W``.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q._module_constructor(V,W)
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
        """
        return FGP_Module(V,W,check)

    def _coerce_map_from_(self, S):
        """
        Return whether ``S`` canonically coerces to ``self``a.

        INPUT:

        - ``S`` -- anything.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: V = span([[5, -1/2]],ZZ); W = span([[20,-2]],ZZ); Q = V/W; phi=Q.hom([2*Q.0])
            sage: Q._coerce_map_from_(ZZ)
            False
            sage: Q._coerce_map_from_(phi.kernel())
            True
            sage: Q._coerce_map_from_(Q)
            True
            sage: Q._coerce_map_from_(phi.image())
            True
            sage: Q._coerce_map_from_(V/V.zero_submodule())
            True
            sage: Q._coerce_map_from_(V/V)
            False
            sage: Q._coerce_map_from_(ZZ^2)
            False

        Of course, `V` canonically coerces to `Q`, as does twice `V`::

            sage: Q._coerce_map_from_(V)
            True
            sage: Q._coerce_map_from_(V.scale(2))
            True
        """
        if is_FGP_Module(S):
            return S.has_canonical_map_to(self)
        return self._V.has_coerce_map_from(S)

    def _repr_(self):
        """
        Return string representation of this module.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: (V/W)._repr_()
            'Finitely generated module V/W over Integer Ring with invariants (4, 12)'
        """
        I = str(self.invariants()).replace(',)',')')
        return "Finitely generated module V/W over %s with invariants %s"%(self.base_ring(), I)

    def __div__(self, other):
        """
        Return the quotient self/other, where other must be a
        submodule of self.

        EXAMPLES::

            sage: V = span([[5, -1/2]],ZZ); W = span([[20,-2]],ZZ); Q = V/W; phi=Q.hom([2*Q.0])
            sage: Q
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: Q/phi.kernel()
            Finitely generated module V/W over Integer Ring with invariants (2)
            sage: Q/Q
            Finitely generated module V/W over Integer Ring with invariants ()
        """
        if not is_FGP_Module(other):
            if is_FreeModule(other):
                other = other / other.zero_submodule()
            else:
                raise TypeError("other must be an FGP module")
        if not other.is_submodule(self):
            raise ValueError("other must be a submodule of self")
        return self._module_constructor(self._V, other._V+self._W)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q == Q
            True
            sage: loads(dumps(Q)) == Q
            True
            sage: Q == V
            False
            sage: Q == V/V.zero_submodule()
            False
        """
        if not is_FGP_Module(other):
            return False
        return self._V == other._V and self._W == other._W

    def __ne__(self, other):
        """
        True iff self is not equal to other.

        This may not be needed for modules created using the function
        :func:`FGP_Module`, since those have uniqueness built into
        them, but if the modules are created directly using the
        __init__ method for this class, then this may fail; in
        particular, for modules M and N with ``M == N`` returning
        True, it may be the case that ``M != N`` may also return True.
        In particular, for derived classes whose __init__ methods just
        call the __init__ method for this class need this.  See
        http://trac.sagemath.org/sage_trac/ticket/9940 for
        illustrations.

        EXAMPLES:

        Make sure that the problems in
        http://trac.sagemath.org/sage_trac/ticket/9940 are fixed::

            sage: G = AdditiveAbelianGroup([0,0])
            sage: H = AdditiveAbelianGroup([0,0])
            sage: G == H
            True
            sage: G != H # indirect doctest
            False

            sage: N1 = ToricLattice(3)
            sage: sublattice1 = N1.submodule([(1,1,0), (3,2,1)])
            sage: Q1 = N1/sublattice1
            sage: N2 = ToricLattice(3)
            sage: sublattice2 = N2.submodule([(1,1,0), (3,2,1)])
            sage: Q2 = N2/sublattice2
            sage: Q1 == Q2
            True
            sage: Q1 != Q2
            False
        """
        return not self.__eq__(other)

    # __le__ is a synonym for `is_submodule`: see below

    def __lt__(self, other):
        """
        True iff self is a proper submodule of other.

        EXAMPLES::

            sage: V = ZZ^2; W = V.span([[1,2]]); W2 = W.scale(2)
            sage: A = V/W; B = W/W2
            sage: B < A
            False
            sage: A = V/W2; B = W/W2
            sage: B < A
            True
            sage: A < A
            False
        """
        return self.__le__(other) and not self.__eq__(other)

    def __gt__(self, other):
        """
        True iff other is a proper submodule of self.

        EXAMPLES::

            sage: V = ZZ^2; W = V.span([[1,2]]); W2 = W.scale(2)
            sage: A = V/W; B = W/W2
            sage: A > B
            False
            sage: A = V/W2; B = W/W2
            sage: A > B
            True
            sage: A > A
            False
        """
        return self.__ge__(other) and not self.__eq__(other)

    def __ge__(self, other):
        """
        True iff other is a submodule of self.

        EXAMPLES::

            sage: V = ZZ^2; W = V.span([[1,2]]); W2 = W.scale(2)
            sage: A = V/W; B = W/W2
            sage: A >= B
            False
            sage: A = V/W2; B = W/W2
            sage: A >= B
            True
            sage: A >= A
            True
        """
        return other.is_submodule(self)


    def _element_constructor_(self, x, check=True):
        """
        INPUT:

        - ``x`` -- a vector, an fgp module element, or a list or tuple:

               - list or tuple: take the corresponding linear combination of
                 the generators of self.

               - vector: coerce vector into V and define the
                 corresponding element of V/W

               - fgp module element: lift to element of ambient vector
                 space and try to put into V.  If x is in self already,
                 just return x.

        - `check` -- bool (default: True)

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: x = Q(V.0-V.1); x  # indirect doctest
            (0, 3)
            sage: type(x)
            <class 'sage.modules.fg_pid.fgp_element.FGP_Module_class_with_category.element_class'>
            sage: x is Q(x)
            True
            sage: x.parent() is Q
            True
        """
#        print '_element_constructor_', x, check
        if isinstance(x, (list,tuple)):
            try:
                x = self.optimized()[0].V().linear_combination_of_basis(x)
            except ValueError as msg:
                raise TypeError(msg)
        elif isinstance(x, FGP_Element):
            x = x.lift()
        return self.element_class(self, self._V(x))


    def linear_combination_of_smith_form_gens(self, x):
        r"""
        Compute a linear combination of the optimised generators of this module
        as returned by :meth:`.smith_form_gens`.

        EXAMPLE::

            sage: X = ZZ**2 / span([[3,0],[0,2]], ZZ)
            sage: X.linear_combination_of_smith_form_gens([1])
            (1)

        """
        try:
            x = self.optimized()[0].V().linear_combination_of_basis(x)
        except ValueError as msg:
            raise TypeError(msg)
        return self.element_class(self, self._V(x))

    def __contains__(self, x):
        """
        Return true if x is contained in self.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.0 in Q
            True
            sage: 0 in Q
            True
            sage: vector([1,2,3/7]) in Q
            False
            sage: vector([1,2,3]) in Q
            True
            sage: Q.0 - Q.1 in Q
            True
        """
        if hasattr(x, 'parent') and x.parent() == self:
            return True
        try:
            self(x)
            return True
        except TypeError:
            return False

    def submodule(self, x):
        """
        Return the submodule defined by x.

        INPUT:

        - ``x`` -- list, tuple, or FGP module

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.gens()
            ((1, 0), (0, 1))

        We create submodules generated by a list or tuple of elements::

            sage: Q.submodule([Q.0])
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: Q.submodule([Q.1])
            Finitely generated module V/W over Integer Ring with invariants (12)
            sage: Q.submodule((Q.0, Q.0 + 3*Q.1))
            Finitely generated module V/W over Integer Ring with invariants (4, 4)

        A submodule defined by a submodule::

            sage: A = Q.submodule((Q.0, Q.0 + 3*Q.1)); A
            Finitely generated module V/W over Integer Ring with invariants (4, 4)
            sage: Q.submodule(A)
            Finitely generated module V/W over Integer Ring with invariants (4, 4)

        Inclusion is checked::

            sage: A.submodule(Q)
            Traceback (most recent call last):
            ...
            ValueError: x.V() must be contained in self's V.
        """
        if is_FGP_Module(x):
            if not x._W.is_submodule(self._W):
                raise ValueError("x.W() must be contained in self's W.")

            V = x._V
            if not V.is_submodule(self._V):
                raise ValueError("x.V() must be contained in self's V.")

            return x

        if not isinstance(x, (list, tuple)):
            raise TypeError("x must be a list, tuple, or FGP module")

        x = Sequence(x)
        if is_FGP_Module(x.universe()):
            # TODO: possibly inefficient in some cases
            x = [self(v).lift() for v in x]
        V = self._V.submodule(x) + self._W
        return self._module_constructor(V, self._W)

    def has_canonical_map_to(self, A):
        """
        Return True if self has a canonical map to A, relative to the
        given presentation of A. This means that A is a finitely
        generated quotient module, self.V() is a submodule of A.V()
        and self.W() is a submodule of A.W(), i.e., that there is a
        natural map induced by inclusion of the V's. Note that we do
        *not* require that this natural map be injective; for this use
        :meth:`is_submodule`.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: A = Q.submodule((Q.0, Q.0 + 3*Q.1)); A
            Finitely generated module V/W over Integer Ring with invariants (4, 4)
            sage: A.has_canonical_map_to(Q)
            True
            sage: Q.has_canonical_map_to(A)
            False

        """
        if not is_FGP_Module(A):
            return False
        if self.cardinality() == 0 and self.base_ring() == A.base_ring():
            return True
        return self.V().is_submodule(A.V()) and self.W().is_submodule(A.W())

    def is_submodule(self, A):
        """
        Return True if self is a submodule of A.  More precisely, this returns True if
        if ``self.V()`` is a submodule of ``A.V()``, with ``self.W()`` equal to ``A.W()``.

        Compare :meth:`.has_canonical_map_to`.

        EXAMPLES::

            sage: V = ZZ^2; W = V.span([[1,2]]); W2 = W.scale(2)
            sage: A = V/W; B = W/W2
            sage: B.is_submodule(A)
            False
            sage: A = V/W2; B = W/W2
            sage: B.is_submodule(A)
            True

        This example illustrates that this command works in a subtle cases.::

            sage: A = ZZ^1
            sage: Q3 = A / A.span([[3]])
            sage: Q6 = A / A.span([[6]])
            sage: Q6.is_submodule(Q3)
            False
            sage: Q6.has_canonical_map_to(Q3)
            True
            sage: Q = A.span([[2]]) / A.span([[6]])
            sage: Q.is_submodule(Q6)
            True
        """
        if not self.has_canonical_map_to(A):
            return False

        return self.V().is_submodule(A.V()) and self.W() == A.W()

    __le__ = is_submodule

    def V(self):
        """
        If this module was constructed as a quotient V/W, returns V.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.V()
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1/2   0   0]
            [  0   1   0]
            [  0   0   1]

        """
        return self._V

    def cover(self):
        """
        If this module was constructed as V/W, returns the cover module V.  This is the same as self.V().

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.V()
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1/2   0   0]
            [  0   1   0]
            [  0   0   1]
        """
        return self.V()

    def W(self):
        """
        If this module was constructed as a quotient V/W, returns W.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.W()
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1/2   8   0]
            [  0  12   0]
            [  0   0   4]
        """
        return self._W

    def relations(self):
        """
        If this module was constructed as V/W, returns the relations module V. This is the same as self.W().

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.relations()
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1/2   8   0]
            [  0  12   0]
            [  0   0   4]
        """
        return self.W()

    @cached_method
    def _relative_matrix(self):
        """
        V has a fixed choice of basis, and W has a fixed choice of
        basis, and both V and W are free R-modules.  Since W is
        contained in V, we can write each basis element of W as an
        R-linear combination of the basis for V.  This function
        returns that matrix over R, where each row corresponds to a
        basis element of W.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q._relative_matrix()
            [ 1  8  0]
            [ 0 12  0]
            [ 0  0  4]
        """
        V = self._V
        W = self._W
        A = V.coordinate_module(W).basis_matrix().change_ring(self.base_ring())
        return A

    @cached_method
    def _smith_form(self):
        """
        Return matrices S, U, and V such that S = U*R*V, and S is in
        Smith normal form, and R is the relative matrix that defines
        self (see :meth:`._relative_matrix`).

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q._smith_form()
            (
            [ 1  0  0]  [1 0 0]  [ 1  0 -8]
            [ 0  4  0]  [0 0 1]  [ 0  0  1]
            [ 0  0 12], [0 1 0], [ 0  1  0]
            )
        """
        return self._relative_matrix().smith_form()

    def base_ring(self):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.base_ring()
            Integer Ring
        """
        return self._V.base_ring()

    @cached_method
    def invariants(self, include_ones=False):
        """
        Return the diagonal entries of the smith form of the relative
        matrix that defines self (see :meth:`._relative_matrix`)
        padded with zeros, excluding 1's by default.   Thus if v is the
        list of integers returned, then self is abstractly isomorphic to
        the product of cyclic groups `Z/nZ` where `n` is in `v`.

        INPUT:

            - ``include_ones`` -- bool (default: False); if True, also
              include 1's in the output list.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.invariants()
            (4, 12)

        An example with 1 and 0 rows::

            sage: V = ZZ^3; W = V.span([[1,2,0],[0,1,0], [0,2,0]]); Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (0)
            sage: Q.invariants()
            (0,)
            sage: Q.invariants(include_ones=True)
            (1, 1, 0)

        """
        D,_,_ = self._smith_form()

        v = [D[i,i] for i in range(D.nrows())] + [Integer(0)]*(D.ncols()-D.nrows())
        w = tuple([x for x in v if x != 1])
        v = tuple(v)
        self.invariants.set_cache(v, True)
        self.invariants.set_cache(w, False)
        return self.invariants(include_ones)

    def gens(self):
        """
        Returns tuple of elements `g_0,...,g_n` of self such that the module generated by
        the gi is isomorphic to the direct sum of R/ei*R, where ei are the
        invariants of self and R is the base ring.

        Note that these are not generally uniquely determined, and depending on
        how Smith normal form is implemented for the base ring, they may not
        even be deterministic.

        This can safely be overridden in all derived classes.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.gens()
            ((1, 0), (0, 1))
            sage: Q.0
            (1, 0)
        """
        return self.smith_form_gens()

    @cached_method
    def smith_form_gens(self):
        """
        Return a set of generators for self which are in Smith normal form.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.smith_form_gens()
            ((1, 0), (0, 1))
            sage: [x.lift() for x in Q.smith_form_gens()]
            [(0, 0, 1), (0, 1, 0)]
        """
        # Get the rightmost transformation in the Smith form
        _, _, X = self._smith_form()
        # Invert it to get a matrix whose rows (in terms of the basis for V)
        # are the gi (including 1 invariants).
        Y = X**(-1)
        # Get the basis matrix for V
        B = self._V.basis_matrix()
        # Multiply to express the gi in terms of the ambient vector space.
        Z = Y*B
        # Make gens out of the rows of Z that correspond to non-1 invariants.
        v = self.invariants(include_ones=True)
        non1 = [i for i in range(Z.nrows()) if v[i] != 1]
        Z = Z.matrix_from_rows(non1)
        self._gens = tuple([self(z, check=DEBUG) for z in Z.rows()])
        return self._gens

    def coordinate_vector(self, x, reduce=False):
        """
        Return coordinates of x with respect to the optimized
        representation of self.

        INPUT:

        - ``x`` -- element of self

        - ``reduce`` -- (default: False); if True, reduce
          coefficients modulo invariants; this is
          ignored if the base ring isn't ZZ.

        OUTPUT:

        The coordinates as a vector. That is, the same type as
        ``self.V()``, but in general with fewer entries.

        EXAMPLES::

            sage: V = span([[1/4,0,0],[3/4,4,2],[0,0,2]],ZZ); W = V.span([4*V.0+12*V.1])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 0, 0)
            sage: Q.coordinate_vector(-Q.0)
            (-1, 0, 0)
            sage: Q.coordinate_vector(-Q.0, reduce=True)
            (3, 0, 0)

        If x isn't in self, it is coerced in::

            sage: Q.coordinate_vector(V.0)
            (1, 0, -3)
            sage: Q.coordinate_vector(Q(V.0))
            (1, 0, -3)


        TESTS::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.coordinate_vector(Q.0 - Q.1)
            (1, -1)

            sage: O, X = Q.optimized()
            sage: O.V()
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [0 0 1]
            [0 2 0]
            sage: phi = Q.hom([Q.0, 4*Q.1])
            sage: x = Q(V.0); x
            (0, 4)
            sage: Q.coordinate_vector(x, reduce=True)
            (0, 4)
            sage: Q.coordinate_vector(-x, reduce=False)
            (0, -4)
            sage: x == 4*Q.1
            True
            sage: x = Q(V.1); x
            (0, 1)
            sage: Q.coordinate_vector(x)
            (0, 1)
            sage: x == Q.1
            True
            sage: x = Q(V.2); x
            (1, 0)
            sage: Q.coordinate_vector(x)
            (1, 0)
            sage: x == Q.0
            True
        """
        try: T = self.__T
        except AttributeError:
            self.optimized() # computes T as side effect -- see the "optimized" method.
            T = self.__T

        x = self(x)

        c = self._V.coordinate_vector(x.lift())
        b = (c*T).change_ring(self.base_ring())
        if reduce and self.base_ring() == ZZ:

            I = self.invariants()
            return b.parent()([b[i] if I[i] == 0 else b[i] % I[i]   for i in range(len(I))])

        else:
            # Don't know (or not requested) canonical way to reduce
            # each entry yet, or how to compute invariants.
            return b

    def gen(self, i):
        """
        Return the i-th generator of self.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.gen(0)
            (1, 0)
            sage: Q.gen(1)
            (0, 1)
            sage: Q.gen(2)
            Traceback (most recent call last):
            ...
            ValueError: Generator 2 not defined
            sage: Q.gen(-1)
            Traceback (most recent call last):
            ...
            ValueError: Generator -1 not defined
        """
        v = self.gens()
        if i < 0 or i >= len(v):
            raise ValueError("Generator %s not defined"%i)
        return v[i]

    def smith_form_gen(self, i):
        """
        Return the i-th generator of self. A private name (so we can freely
        override gen() in derived classes).

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.smith_form_gen(0)
            (1, 0)
            sage: Q.smith_form_gen(1)
            (0, 1)
        """
        v = self.smith_form_gens()
        if i < 0 or i >= len(v):
            raise ValueError("Smith form generator %s not defined"%i)
        return v[i]

    def optimized(self):
        """
        Return a module isomorphic to this one, but with V replaced by
        a submodule of V such that the generators of self all lift
        trivially to generators of V.  Replace W by the intersection
        of V and W. This has the advantage that V has small dimension
        and any homomorphism from self trivially extends to a
        homomorphism from V.

        OUTPUT:

        - ``Q`` -- an optimized quotient V0/W0 with V0 a submodule of V
          such that phi: V0/W0 --> V/W is an isomorphism

        - ``Z`` -- matrix such that if x is in self.V() and
          c gives the coordinates of x in terms of the
          basis for self.V(), then c*Z is in V0
          and c*Z maps to x via phi above.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: O, X = Q.optimized(); O
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: O.V()
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [0 0 1]
            [0 1 0]
            sage: O.W()
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 0 12  0]
            [ 0  0  4]
            sage: X
            [0 4 0]
            [0 1 0]
            [0 0 1]
            sage: OV = O.V()
            sage: Q(OV([0,-8,0])) == V.0
            True
            sage: Q(OV([0,1,0])) == V.1
            True
            sage: Q(OV([0,0,1])) == V.2
            True
        """
        try:
            if self.__optimized is True:
                return self, None
            return self.__optimized
        except AttributeError: pass
        V = self._V.span_of_basis([x.lift() for x in self.smith_form_gens()])
        M = self._module_constructor(V, self._W.intersection(V))
        # Compute matrix T of linear transformation from self._V to V.
        # This matrix T gives each basis element of self._V in terms
        # of our new optimized V, modulo the W's.
        A = V.basis_matrix().stack(self._W.basis_matrix())
        B, d = A._clear_denom()
        H, U = B.hermite_form(transformation=True)
        Y = H.solve_left(d*self._V.basis_matrix())
        T = Y * U.matrix_from_columns(range(V.rank()))
        self.__T = T

        # Finally we multiply by V.basis_matrix() so X gives vectors
        # in the ambient space instead of coefficients of linear
        # combinations of the basis for V.
        X = T * V.basis_matrix()

        self.__optimized = M, X
        return M, X


    def hom(self, im_gens, codomain=None, check=True):
        """
        Homomorphism defined by giving the images of ``self.gens()`` in some
        fixed fg R-module.

        .. note ::

            We do not assume that the generators given by ``self.gens()`` are
            the same as the Smith form generators, since this may not be true
            for a general derived class.

        INPUTS:

        - ``im_gens`` -- a list of the images of ``self.gens()`` in some
          R-module


        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: phi = Q.hom([3*Q.1, Q.0])
            sage: phi
            Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(0, 3), (1, 0)]
            sage: phi(Q.0)
            (0, 3)
            sage: phi(Q.1)
            (1, 0)
            sage: Q.0 == phi(Q.1)
            True

        This example illustrates creating a morphism to a free module.
        The free module is turned into an FGP module (i.e., quotient
        V/W with W=0), and the morphism is constructed::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (2, 0, 0)
            sage: phi = Q.hom([0,V.0,V.1]); phi
            Morphism from module over Integer Ring with invariants (2, 0, 0) to module with invariants (0, 0, 0) that sends the generators to [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
            sage: phi.domain()
            Finitely generated module V/W over Integer Ring with invariants (2, 0, 0)
            sage: phi.codomain()
            Finitely generated module V/W over Integer Ring with invariants (0, 0, 0)
            sage: phi(Q.0)
            (0, 0, 0)
            sage: phi(Q.1)
            (1, 0, 0)
            sage: phi(Q.2) == V.1
            True

        Constructing two zero maps from the zero module::

            sage: A = (ZZ^2)/(ZZ^2); A
            Finitely generated module V/W over Integer Ring with invariants ()
            sage: A.hom([])
            Morphism from module over Integer Ring with invariants () to module with invariants () that sends the generators to []
            sage: A.hom([]).codomain() is A
            True
            sage: B = (ZZ^3)/(ZZ^3)
            sage: A.hom([],codomain=B)
            Morphism from module over Integer Ring with invariants () to module with invariants () that sends the generators to []
            sage: phi = A.hom([],codomain=B); phi
            Morphism from module over Integer Ring with invariants () to module with invariants () that sends the generators to []
            sage: phi(A(0))
            ()
            sage: phi(A(0)) == B(0)
            True


        A degenerate case::

            sage: A = (ZZ^2)/(ZZ^2)
            sage: phi = A.hom([]); phi
            Morphism from module over Integer Ring with invariants () to module with invariants () that sends the generators to []
            sage: phi(A(0))
            ()

        The code checks that the morphism is valid.  In the example
        below we try to send a generator of order 2 to an element of
        order 14::

            sage: V = span([[1/14,3/14],[0,1/2]],ZZ); W = ZZ^2
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (2, 14)
            sage: Q([1,11]).additive_order()
            14
            sage: f = Q.hom([Q([1,11]), Q([1,3])]); f
            Traceback (most recent call last):
            ...
            ValueError: phi must send optimized submodule of M.W() into N.W()

        """
        if len(im_gens) == 0:
            # 0 map
            N = self if codomain is None else codomain
        else:
            if codomain is None:
                im_gens = Sequence(im_gens)
                N = im_gens.universe()
            else:
                N = codomain
                im_gens = Sequence(im_gens, universe=N)

        if is_FreeModule(N):
            # If im_smith_gens aren't in an R-module, but are in a Free-module,
            # then we quotient out by the 0 submodule and get an R-module.
            N = FGP_Module(N, N.zero_submodule(), check=DEBUG)
            im_gens = Sequence(im_gens, universe=N)

        if len(im_gens) == 0:
            VO = self.optimized()[0].V()
            H = VO.Hom(N.V())
            return FGP_Morphism(self.Hom(N), H(0), check=DEBUG)

        if self.gens() == self.smith_form_gens():
            return self._hom_from_smith(im_gens, check)
        else:
            return self._hom_general(im_gens, check)


    def _hom_general(self, im_gens, check=True):
        """
        Homomorphism defined by giving the images of ``self.gens()`` in some
        fixed fg R-module. We do not assume that the generators given by
        ``self.gens()`` are the same as the Smith form generators, since this
        may not be true for a general derived class.

        INPUT:

        - ``im_gens`` - a Sequence object giving the images of ``self.gens()``,
          whose universe is some fixed fg R-module

        EXAMPLE::

            sage: class SillyModule(sage.modules.fg_pid.fgp_module.FGP_Module_class):
            ...       def gens(self):
            ...           return tuple(flatten([[x,x] for x in self.smith_form_gens()]))
            sage: A = SillyModule(ZZ**1, span([[3]], ZZ))
            sage: A.gen(0)
            (1)
            sage: A.gen(1)
            (1)
            sage: B = ZZ**1 / span([[3]], ZZ)
            sage: A.hom([B.0, 2*B.0], B)
            Traceback (most recent call last):
            ...
            ValueError: Images do not determine a valid homomorphism
            sage: A.hom([B.0, B.0], B)   # indirect doctest
            Morphism from module over Integer Ring with invariants (3,) to module with invariants (3,) that sends the generators to [(1), (1)]

        """
        m = self.ngens()
        A = ZZ**m
        q = A.hom([x.lift() for x in self.gens()], self.V())
        B = q.inverse_image(self.W())
        N = im_gens.universe()
        r = A.hom([x.lift() for x in im_gens], N.V())
        if check:
            if not r(B).is_submodule(N.W()):
                raise ValueError("Images do not determine a valid homomorphism")
        smith_images = Sequence([N(r(q.lift(x.lift()))) for x in self.smith_form_gens()])
        return self._hom_from_smith(smith_images, check=DEBUG)

    def _hom_from_smith(self, im_smith_gens, check=True):
        """
        Homomorphism defined by giving the images of the Smith-form generators
        of self in some fixed fg R-module.

        INPUTS:

        - ``im_gens`` -- a Sequence object giving the images of the Smith-form
          generators of self, whose universe is some fixed fg R-module

        EXAMPLE::

            sage: class SillyModule(sage.modules.fg_pid.fgp_module.FGP_Module_class):
            ...       def gens(self):
            ...           return tuple(flatten([[x,x] for x in self.smith_form_gens()]))
            sage: A = SillyModule(ZZ**1, span([[3]], ZZ))
            sage: A.gen(0)
            (1)
            sage: A.gen(1)
            (1)
            sage: B = ZZ**1 / span([[3]], ZZ)
            sage: A._hom_from_smith(Sequence([B.0]))
            Morphism from module over Integer Ring with invariants (3,) to module with invariants (3,) that sends the generators to [(1), (1)]
        """
        if len(im_smith_gens) != len(self.smith_form_gens()):
            raise ValueError("im_gens must have length the same as self.smith_form_gens()")

        # replace self by representation in which smith-gens g_i are a basis for V.
        M, _ = self.optimized()
        # Define morphism from M to N
        f = M.V().hom([x.lift() for x in im_smith_gens])
        N = im_smith_gens.universe()
        homspace = self.Hom(N)
        phi = FGP_Morphism(homspace, f, check=DEBUG)
        return phi

    def _Hom_(self, N, category=None):
        """
        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([V.0+2*V.1, 9*V.0+2*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.Hom(Q)     # indirect doctest
            Set of Morphisms from Finitely generated module V/W over Integer Ring with invariants (4, 16) to Finitely generated module V/W over Integer Ring with invariants (4, 16) in Category of modules over Integer Ring
            sage: M = V/V.zero_submodule()
            sage: H = M.Hom(Q); H
            Set of Morphisms from Finitely generated module V/W over Integer Ring with invariants (0, 0, 0) to Finitely generated module V/W over Integer Ring with invariants (4, 16) in Category of modules over Integer Ring
            sage: Hom(M,Q) is H
            True
            sage: type(Hom(M,Q))
            <class 'sage.modules.fg_pid.fgp_morphism.FGP_Homset_class_with_category'>
            sage: H.category()
            Category of hom sets in Category of modules over Integer Ring
            sage: H.homset_category()
            Category of modules over Integer Ring

        """
        return FGP_Homset(self, N)

    def random_element(self, *args, **kwds):
        """
        Create a random element of self=V/W, by creating a random element of V and
        reducing it modulo W.

        All arguments are passed onto the random_element method of V.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: Q.random_element()
            (1, 10)
        """
        return self(self._V.random_element(*args, **kwds))

    def cardinality(self):
        """
        Return the cardinality of this module as a set.

        EXAMPLES::

            sage: V = ZZ^2; W = V.span([[1,2],[3,4]]); A = V/W; A
            Finitely generated module V/W over Integer Ring with invariants (2)
            sage: A.cardinality()
            2
            sage: V = ZZ^2; W = V.span([[1,2]]); A = V/W; A
            Finitely generated module V/W over Integer Ring with invariants (0)
            sage: A.cardinality()
            +Infinity
            sage: V = QQ^2; W = V.span([[1,2]]); A = V/W; A
            Vector space quotient V/W of dimension 1 over Rational Field where
            V: Vector space of dimension 2 over Rational Field
            W: Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 2]
            sage: A.cardinality()
            +Infinity
        """
        try:
            return self.__cardinality
        except AttributeError:
            pass
        from sage.rings.all import infinity
        from sage.misc.all import prod
        v = self.invariants()
        self.__cardinality = infinity if 0 in v else prod(v)
        return self.__cardinality

    def __iter__(self):
        """
        Return iterator over all elements of self.

        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([V.0+2*V.1, 4*V.0+2*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (2, 12)
            sage: z = list(V/W)
            sage: print z
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (0, 10), (0, 11), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11)]
            sage: len(z)
            24

        We test that the trivial module is handled correctly (cf. trac #6561)::

            sage: A = (ZZ**1)/(ZZ**1); list(A) == [A(0)]
            True
        """
        if self.base_ring() != ZZ:
            raise NotImplementedError("only implemented over ZZ")
        v = self.invariants()
        if 0 in v:
            raise NotImplementedError("currently self must be finite to iterate over")
        B = self.optimized()[0].V().basis_matrix()
        V = self.base_ring()**B.nrows()
        from sage.misc.mrange import cartesian_product_iterator
        for a in cartesian_product_iterator([range(k) for k in v]):
            b = V(a) * B
            yield self(b)

    def is_finite(self):
        """
        Return True if self is finite and False otherwise.

        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([V.0+2*V.1, 9*V.0+2*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 16)
            sage: Q.is_finite()
            True
            sage: Q = V/V.zero_submodule(); Q
            Finitely generated module V/W over Integer Ring with invariants (0, 0, 0)
            sage: Q.is_finite()
            False
        """
        return 0 not in self.invariants()

    def annihilator(self):
        """
        Return the ideal of the base ring that annihilates self. This
        is precisely the ideal generated by the LCM of the invariants
        of self if self is finite, and is 0 otherwise.

        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([V.0+2*V.1, 9*V.0+2*V.1, 4*V.2])
            sage: Q = V/W; Q.annihilator()
            Principal ideal (16) of Integer Ring
            sage: Q.annihilator().gen()
            16

            sage: Q = V/V.span([V.0]); Q
            Finitely generated module V/W over Integer Ring with invariants (0, 0)
            sage: Q.annihilator()
            Principal ideal (0) of Integer Ring
        """
        if not self.is_finite():
            g = 0
        else:
            g = reduce(lcm, self.invariants())
        return self.base_ring().ideal(g)

    def ngens(self):
        r"""
        Return the number of generators of self.

        (Note for developers: This is just the length of :meth:`.gens`, rather
        than of the minimal set of generators as returned by
        :meth:`.smith_form_gens`; these are the same in the
        :class:`~sage.modules.fg_pid.fgp_module.FGP_Module_class`, but not
        necessarily in derived classes.)

        EXAMPLES::

            sage: A = (ZZ**2) / span([[4,0],[0,3]], ZZ)
            sage: A.ngens()
            1

        This works (but please don't do it in production code!) ::

            sage: A.gens = lambda: [1,2,"Barcelona!"]
            sage: A.ngens()
            3
        """
        return len(self.gens())

    def __hash__(self):
        r"""
        Calculate a hash for self.

        EXAMPLES::

            sage: A = (ZZ**2) / span([[4,0],[0,3]], ZZ)
            sage: hash(A)
            1328975982 # 32-bit
            -7071641102956720018 # 64-bit
        """
        return hash( (self.V(), self.W()) )

##############################################################
# Useful for testing
##############################################################

def random_fgp_module(n, R=ZZ, finite=False):
    """
    Return a random FGP module inside a rank n free module over R.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``R`` -- base ring (default: ZZ)

    - ``finite`` -- bool (default: True); if True, make the random module finite.

    EXAMPLES::

        sage: import sage.modules.fg_pid.fgp_module as fgp
        sage: fgp.random_fgp_module(4)
        Finitely generated module V/W over Integer Ring with invariants (4)
    """
    K = R.fraction_field()
    V = K**n
    i = ZZ.random_element(max(n,1))
    A = V.span([V.random_element() for _ in range(i)], R)
    if not finite:
        i = ZZ.random_element(i+1)
    while True:
        B = A.span([A.random_element() for _ in range(i)], R)
        #Q = A/B
        Q = FGP_Module_class(A,B,check=DEBUG)
        if not finite or Q.is_finite():
            return Q

def random_fgp_morphism_0(*args, **kwds):
    """
    Construct a random fgp module using random_fgp_module,
    then construct a random morphism that sends each generator
    to a random multiple of itself.  Inputs are the same
    as to random_fgp_module.

    EXAMPLES::

        sage: import sage.modules.fg_pid.fgp_module as fgp
        sage: fgp.random_fgp_morphism_0(4)
        Morphism from module over Integer Ring with invariants (4,) to module with invariants (4,) that sends the generators to [(0)]

    """
    A = random_fgp_module(*args, **kwds)
    return A.hom([ZZ.random_element()*g for g in A.smith_form_gens()])

def test_morphism_0(*args, **kwds):
    """
    EXAMPLES::

        sage: import sage.modules.fg_pid.fgp_module as fgp
        sage: s = 0  # we set a seed so results clearly and easily reproducible across runs.
        sage: set_random_seed(s); v = [fgp.test_morphism_0(1) for _ in range(30)]
        sage: set_random_seed(s); v = [fgp.test_morphism_0(2) for _ in range(30)]
        sage: set_random_seed(s); v = [fgp.test_morphism_0(3) for _ in range(10)]
        sage: set_random_seed(s); v = [fgp.test_morphism_0(i) for i in range(1,20)]
        sage: set_random_seed(s); v = [fgp.test_morphism_0(4) for _ in range(50)]    # long time
    """
    phi = random_fgp_morphism_0(*args, **kwds)
    K = phi.kernel()
    I = phi.image()
    from sage.misc.all import prod
    if prod(K.invariants()):
        assert prod(phi.domain().invariants()) % prod(K.invariants()) == 0
    assert I.is_submodule(phi.codomain())
    if len(I.smith_form_gens()) > 0:
        x = phi.lift(I.smith_form_gen(0))
        assert phi(x) == I.smith_form_gen(0)


