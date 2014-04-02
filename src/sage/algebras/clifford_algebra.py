r"""
Clifford Algebras

AUTHORS:

- Travis Scrimshaw (2013-09-06): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from copy import copy

from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.categories.modules_with_basis import ModuleMorphismByLinearity
from sage.rings.all import ZZ
from sage.modules.free_module import FreeModule, FreeModule_generic
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.subset import Subsets
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.algebras.weyl_algebra import repr_from_monomials

class CliffordAlgebraElement(CombinatorialFreeModule.Element):
    """
    An element in a Clifford algebra.

    TESTS::

        sage: Q = QuadraticForm(ZZ, 3, [1, 2, 3, 4, 5, 6])
        sage: Cl.<x,y,z> = CliffordAlgebra(Q)
        sage: elt = ((x^3-z)*x + y)^2
        sage: TestSuite(elt).run()
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: ((x^3-z)*x + y)^2
            -2*x*y*z - x*z + 5*x - 4*y + 2*z + 2
            sage: Cl.zero()
            0
        """
        return repr_from_monomials(self.list(), self.parent()._repr_term)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        TESTS::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: latex( ((x^3-z)*x + y)^2 )
            -2  x y z -  x z + 5  x - 4  y + 2  z + 2
            sage: Cl.<x0,x1,x2> = CliffordAlgebra(Q)
            sage: latex(  (x1 - x2)*x0 + 5*x0*x1*x2 )
            5  x_{0} x_{1} x_{2} -  x_{0} x_{1} +  x_{0} x_{2} - 1
        """
        return repr_from_monomials(self.list(), self.parent()._latex_term, True)

    def _mul_(self, other):
        """
        Return ``self`` multiplied by ``other``.

        INPUT:

        - ``other`` -- element of the same Clifford algebra as ``self``

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: (x^3 - z*y)*x*(y*z + x*y*z)
            x*y*z + y*z - 24*x + 12*y + 2*z - 24
            sage: y*x
            -x*y + 2
            sage: z*x
            -x*z + 3
            sage: z*z
            6
            sage: x*0
            0
            sage: 0*x
            0
        """
        Q = self.parent()._quadratic_form
        zero = self.parent().base_ring().zero()
        d = {}

        for ml,cl in self:
            # Distribute the current term ``cl`` * ``ml`` over ``other``.
            cur = copy(other._monomial_coefficients) # The current distribution of the term
            for i in reversed(ml):
                # Distribute the current factor ``e[i]`` (the ``i``-th
                # element of the standard basis).
                next = {}
                # At the end of the following for-loop, ``next`` will be
                # the dictionary describing the element
                # ``e[i]`` * (the element described by the dictionary ``cur``)
                # (where ``e[i]`` is the ``i``-th standard basis vector).
                for mr,cr in cur.items():
                    # Commute the factor as necessary until we are in order
                    pos = 0
                    for j in mr:
                        if i <= j:
                            break
                        # Add the additional term from the commutation
                        t = list(mr)
                        t.pop(pos)
                        t = tuple(t)
                        next[t] = next.get(t, zero) + cr * Q[i,j]
                        # Note: ``Q[i,j] == Q(e[i]+e[j]) - Q(e[i]) - Q(e[j])`` for
                        # ``i != j``, where ``e[k]`` is the ``k``-th standard
                        # basis vector.
                        cr = -cr
                        if next[t] == zero:
                            del next[t]
                        pos += 1

                    # Check to see if we have a squared term or not
                    t = list(mr)
                    if i in t:
                        t.remove(i)
                        cr *= Q[i,i]
                        # Note: ``Q[i,i] == Q(e[i])`` where ``e[i]`` is the
                        # ``i``-th standard basis vector.
                    else:
                        t.insert(pos, i)
                        # Note that ``t`` is now sorted.
                    t = tuple(t)
                    next[t] = next.get(t, zero) + cr
                    if next[t] == zero:
                        del next[t]
                cur = next

            # Add the distributed terms to the total
            for index,coeff in cur.items():
                d[index] = d.get(index, zero) + cl * coeff
                if d[index] == zero:
                    del d[index]

        return self.__class__(self.parent(), d)

    def list(self):
        """
        Return the list of monomials and their coefficients in ``self``
        (as a list of `2`-tuples, each of which has the form
        ``(monomial, coefficient)``).

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y
            sage: elt.list()
            [((0,), 5), ((1,), 1)]
        """
        return sorted(self._monomial_coefficients.items(), key=lambda (m,c) : (-len(m), m))

    def support(self):
        """
        Return the support of ``self``.

        This is the list of all monomials which appear with nonzero
        coefficient in ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y
            sage: elt.support()
            [(0,), (1,)]
        """
        return sorted(self._monomial_coefficients.keys(), key=lambda x: (-len(x), x))

    def reflection(self):
        """
        Return the image of the reflection automorphism on ``self``.

        The *reflection automorphism* of a Clifford algebra is defined
        as the linear endomorphism of this algebra which maps

        .. MATH::

            x_1 \wedge x_2 \wedge \cdots \wedge x_m \mapsto
            (-1)^m x_1 \wedge x_2 \wedge \cdots \wedge x_m.

        It is an algebra automorphism of the Clifford algebra.

        :meth:`degree_negation` is an alias for :meth:`reflection`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y + x*z
            sage: r = elt.reflection(); r
            x*z - 5*x - y
            sage: r.reflection() == elt
            True

        TESTS:

        We check that the reflection is an involution::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: all(x.reflection().reflection() == x for x in Cl.basis())
            True
        """
        return self.__class__(self.parent(), {m: (-1)**len(m) * c for m,c in self})

    degree_negation = reflection

    def transpose(self):
        r"""
        Return the transpose of ``self``.

        The transpose is an anti-algebra involution of a Clifford algebra
        and is defined (using linearity) by

        .. MATH::

            x_1 \wedge x_2 \wedge \cdots \wedge x_m \mapsto
            x_m \wedge \cdots \wedge x_2 \wedge x_1.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y + x*z
            sage: t = elt.transpose(); t
            -x*z + 5*x + y + 3
            sage: t.transpose() == elt
            True
            sage: Cl.one().transpose()
            1

        TESTS:

        We check that the transpose is an involution::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: all(x.transpose().transpose() == x for x in Cl.basis())
            True

        Zero is sent to zero::

            sage: Cl.zero().transpose() == Cl.zero()
            True
        """
        P = self.parent()
        if not self._monomial_coefficients:
            return P.zero()
        g = P.gens()
        return P.sum(c * P.prod(g[i] for i in reversed(m)) for m,c in self)

    def conjugate(self):
        r"""
        Return the Clifford conjugate of ``self``.

        The Clifford conjugate of an element `x` of a Clifford algebra is
        defined as

        .. MATH::

            \bar{x} := \alpha(x^t) = \alpha(x)^t

        where `\alpha` denotes the :meth:`reflection <reflection>`
        automorphism and `t` the :meth:`transposition <transpose>`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y + x*z
            sage: c = elt.conjugate(); c
            -x*z - 5*x - y + 3
            sage: c.conjugate() == elt
            True

        TESTS:

        We check that the conjugate is an involution::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: all(x.conjugate().conjugate() == x for x in Cl.basis())
            True
        """
        return self.reflection().transpose()

    clifford_conjugate = conjugate

class CliffordAlgebra(CombinatorialFreeModule):
    r"""
    The Clifford algebra of a quadratic form.

    Let `Q : V \to \mathbf{k}` denote a quadratic form on a vector space `V`
    over a field `\mathbf{k}`. The Clifford algebra `Cl(V, Q)` is defined as
    `T(V) / I_Q` where `T(V)` is the tensor algebra of `V` and `I_Q` is the
    two-sided ideal generated by all elements of the form `v \otimes v - Q(v)`
    for all `v \in V`.

    We abuse notation to denote the projection of a pure tensor
    `x_1 \otimes x_2 \otimes \cdots \otimes x_m \in T(V)` onto
    `T(V) / I_Q = Cl(V, Q)` by `x_1 \wedge x_2 \wedge \cdots \wedge x_m`.
    This is motivated by the fact that `Cl(V, Q)` is the exterior algebra
    `\wedge V` when `Q = 0` (one can also think of a Clifford algebra as
    a quantization of the exterior algebra).

    From the definition, a basis of `Cl(V, Q)` is given by monomials of
    the form

    .. MATH::

        \{ e_{i_1} \wedge \cdots \wedge e_{i_k} \mid 1 \leq i_1 < \cdots <
        i_k \leq n \},

    where `n = \dim(V)` and where `\{ e_1, e_2, \cdots, e_n \}` is any
    fixed basis of `V`. Hence

    .. MATH::

        \dim(Cl(V, Q)) = \sum_{k=0}^n \binom{n}{k} = 2^n.

    .. NOTE::

        The algebra `Cl(V, Q)` is a `\ZZ / 2\ZZ`-graded algebra, but not
        (in general) `\ZZ`-graded (in a reasonable way).

    This construction satisfies the following universal property. Let
    `i : V \to Cl(V, Q)` denote the natural inclusion (which is an
    embedding). Then for every associative `\mathbf{k}`-algebra `A`
    and any `\mathbf{k}`-linear map `j : V \to A` satisfying

    .. MATH::

        j(v)^2 = Q(v) \cdot 1_A

    for all `v \in V`, there exists a unique `\mathbf{k}`-algebra
    homomorphism `f : Cl(V, Q) \to A` such that `f \circ i = j`.
    This property determines the Clifford algebra uniquely up to
    canonical isomorphism. The inclusion `i` is commonly used to
    identify `V` with a vector subspace of `Cl(V)`.

    The Clifford algebra also can be considered as a covariant functor
    from the category of vector spaces equipped with quadratic forms
    to the category of algebras. In fact, if `(V, Q)` and `(W, R)`
    are two vector spaces endowed with quadratic forms, and if
    `g : W \to V` is a linear map preserving the quadratic form,
    then we can define an algebra morphism
    `Cl(g) : Cl(W, R) \to Cl(V, Q)` by requiring that it send every
    `w \in W` to `g(w) \in V`. Since the quadratic form `R` on `W`
    is uniquely determined by the quadratic form `Q` on `V` (due to
    the assumption that `g` preserves the quadratic form), this fact
    can be rewritten as follows: If `(V, Q)` is a vector space with a
    quadratic form, and `W` is another vector space, and
    `\phi : W \to V` is any linear map, then we obtain an algebra
    morphism `Cl(\phi) : Cl(W, \phi(Q)) \to Cl(V, Q)` where
    `\phi(Q) = \phi^T \cdot Q \cdot \phi` (we consider `\phi` as a
    matrix) is the quadratic form `Q` pulled back to `W`. In fact, the
    map `\phi` preserves the quadratic form because of

    .. MATH::

        \phi(Q)(x) = x^T \cdot \phi^T \cdot Q \cdot \phi \cdot x
        = (\phi \cdot x)^T \cdot Q \cdot (\phi \cdot x) = Q(\phi(x)).

    Hence we have `\phi(w)^2 = Q(\phi(w)) = \phi(Q)(w)` for all `w \in W`.

    REFERENCES:

    - :wikipedia:`Clifford_algebra`

    INPUT:

    - ``Q`` -- a quadratic form
    - ``names`` -- (default: ``'e'``) the generator names

    EXAMPLES:

    To create a Clifford algebra, all one needs to do is specify a
    quadratic form::

        sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
        sage: Cl = CliffordAlgebra(Q)
        sage: Cl
        The Clifford algebra of the Quadratic form in 3 variables over Integer Ring with coefficients: 
        [ 1 2 3 ]
        [ * 4 5 ]
        [ * * 6 ]

    We can also explicitly name the generators. In this example, the
    Clifford algebra we construct is an exterior algebra (since we
    choose the quadratic form to be zero)::

        sage: Q = QuadraticForm(ZZ, 4, [0]*10)
        sage: Cl.<a,b,c,d> = CliffordAlgebra(Q)
        sage: a*d
        a*d
        sage: d*c*b*a + a + 4*b*c
        a*b*c*d + 4*b*c + a
    """
    @staticmethod
    def __classcall_private__(cls, Q, names=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl1.<e0,e1,e2> = CliffordAlgebra(Q)
            sage: Cl2 = CliffordAlgebra(Q)
            sage: Cl3 = CliffordAlgebra(Q, ['e0','e1','e2'])
            sage: Cl1 is Cl2 and Cl2 is Cl3
            True
        """
        if not isinstance(Q, QuadraticForm):
            raise ValueError("{} is not a quadratic form".format(Q))
        if names is None:
            names = 'e'
        names = tuple(names)
        if len(names) != Q.dim():
            if len(names) == 1:
                names = tuple( '{}{}'.format(names[0], i) for i in range(Q.dim()) )
            else:
                raise ValueError("the number of variables does not match the number of generators")
        return super(CliffordAlgebra, cls).__classcall__(cls, Q, names)

    def __init__(self, Q, names, category=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl = CliffordAlgebra(Q)
            sage: TestSuite(Cl).run()

        TESTS:

        We check that the basis elements are indeed indexed by
        *strictly increasing* tuples::

            sage: Q = QuadraticForm(ZZ, 9)
            sage: Cl = CliffordAlgebra(Q)
            sage: ba = Cl.basis().keys()
            sage: all( tuple(sorted(S)) in ba
            ....:      for S in Subsets(range(9)) )
            True
        """
        self._quadratic_form = Q
        R = Q.base_ring()
        if category is None:
            category = GradedAlgebrasWithBasis(R)
        indices = map( lambda x: tuple(sorted(x)), Subsets(range(Q.dim())) )
        CombinatorialFreeModule.__init__(self, R, indices, category=category)
        self._assign_names(names)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: CliffordAlgebra(Q)
            The Clifford algebra of the Quadratic form in 3 variables over Integer Ring with coefficients: 
            [ 1 2 3 ]
            [ * 4 5 ]
            [ * * 6 ]
        """
        return "The Clifford algebra of the {}".format(self._quadratic_form)

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl._repr_term((0,2))
            'x*z'
            sage: Cl._repr_term(())
            '1'
            sage: Cl._repr_term((1,))
            'y'
        """
        if len(m) == 0:
            return '1'
        term = ''
        for i in m:
            if len(term) != 0:
                term += '*'
            term += self.variable_names()[i]
        return term

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl._latex_term((0,2))
            ' x z'
        """
        if len(m) == 0:
            return '1'
        term = ''
        for i in m:
            term += ' ' + self.latex_variable_names()[i]
        return term

    def _coerce_map_from_(self, V):
        """
        Return if there is a coerce map from ``V`` into ``self``.

        The things which coerce into ``self`` are:

        - Clifford algebras with the same generator names and an equal
          quadratic form over a ring which coerces into the base
          ring of ``self``.
        - The underlying free module of ``self``.
        - The base ring of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Qp = QuadraticForm(QQ, 3, [1,2,3,4,5,6])
            sage: Cl = CliffordAlgebra(Q)
            sage: Clp = CliffordAlgebra(Qp)
            sage: Cl.has_coerce_map_from(Clp)
            False
            sage: Clp.has_coerce_map_from(Cl)
            True

        Check that we preserve the multiplicative structure::

            sage: all(Clp(b)*Clp(b) == Clp(b*b) for b in Cl.basis())
            True

        Check from the underlying free module::

            sage: M = ZZ^3
            sage: Mp = QQ^3
            sage: Cl.has_coerce_map_from(M)
            True
            sage: Cl.has_coerce_map_from(Mp)
            False
            sage: Clp.has_coerce_map_from(M)
            True
            sage: Clp.has_coerce_map_from(Mp)
            True

        Names matter::

            sage: Cln = CliffordAlgebra(Q, names=['x','y','z'])
            sage: Cln.has_coerce_map_from(Cl)
            False
            sage: Cl.has_coerce_map_from(Cln)
            False
        """
        if isinstance(V, CliffordAlgebra):
            Q = self._quadratic_form
            try:
                return (V.variable_names() == self.variable_names() and
                        V._quadratic_form.base_change_to(self.base_ring()) == Q)
            except Exception:
                return False

        if self.free_module().has_coerce_map_from(V):
            return True

        return super(CliffordAlgebra, self)._coerce_map_from_(V)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Qp = QuadraticForm(QQ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Clp = CliffordAlgebra(Qp, names=['x','y','z'])
            sage: M = ZZ^3
            sage: Mp = QQ^3
            sage: Cl(2/3)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= 2/3) an element of self ... 
            sage: Clp(2/3)
            2/3
            sage: Clp(x)
            x
            sage: M = ZZ^3
            sage: Clp( M((1,-3,2)) )
            x - 3*y + 2*z
        """
        # This is the natural lift morphism of the underlying free module
        if x in self.free_module():
            if x.parent().base_ring() == self.base_ring():
                return self.element_class(self, {(i,): c for i,c in x.iteritems()})
            R = self.base_ring()
            return self.element_class(self, {(i,): R(c) for i,c in x.iteritems()})

        if isinstance(x, CliffordAlgebraElement):
            if x.parent() is self:
                return x
            if self.has_coerce_map_from(x.parent()):
                R = self.base_ring()
                return self.element_class(self, {i: R(c) for i,c in x})

        return super(CliffordAlgebra, self)._element_constructor_(x)

    def gen(self, i):
        """
        Return the ``i``-th standard generator of the algebra ``self``.

        This is the ``i``-th basis vector of the vector space on which
        the quadratic form defining ``self`` is defined, regarded as an
        element of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: [Cl.gen(i) for i in range(3)]
            [x, y, z]
        """
        return self._from_dict({(i,): self.base_ring().one()}, remove_zeros=False)

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.algebra_generators()
            (x, y, z)
        """
        return tuple(self.gen(i) for i in range(self.ngens()))

    gens = algebra_generators

    def ngens(self):
        """
        Return the number of algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.ngens()
            3
        """
        return self._quadratic_form.dim()

    @cached_method
    def one_basis(self):
        """
        Return the basis index of the element `1`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.one_basis()
            ()
        """
        return ()

    def is_commutative(self):
        """
        Check if ``self`` is a commutative algebra.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.is_commutative()
            False
        """
        return self.ngens() < 2

    def quadratic_form(self):
        """
        Return the quadratic form of ``self``.

        This is the quadratic form used to define ``self``. The
        quadratic form on ``self`` is yet to be implemented.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.quadratic_form()
            Quadratic form in 3 variables over Integer Ring with coefficients: 
            [ 1 2 3 ]
            [ * 4 5 ]
            [ * * 6 ]
        """
        return self._quadratic_form

    def degree_on_basis(self, m):
        r"""
        Return the degree of the monomial indexed by ``m``.

        This degree is either `0` or `1`, and should be interpreted as a
        residue class modulo `2`, since we consider ``self`` to be
        `\ZZ_2`-graded (not `\ZZ`-graded, although there is a natural
        *filtration* by the length of ``m``). The degree of the monomial
        ``m`` in this `\ZZ_2`-grading is defined to be the length of ``m``
        taken mod `2`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.degree_on_basis((0,))
            1
            sage: Cl.degree_on_basis((0,1))
            0
        """
        return len(m) % ZZ(2)

    @cached_method
    def free_module(self):
        """
        Return the underlying free module `V` of ``self``.

        This is the free module on which the quadratic form that was
        used to construct ``self`` is defined.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.free_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return FreeModule(self.base_ring(), self._quadratic_form.dim())

    def dimension(self):
        """
        Return the rank of ``self`` as a free module.

        Let `V` be a free `R`-module of rank `n`; then, `Cl(V, Q)` is a
        free `R`-module of rank `2^n`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.dimension()
            8
        """
        return ZZ(2)**self._quadratic_form.dim()

    def lift_module_morphism(self, m, names=None):
        r"""
        Lift the matrix ``m`` to an algebra morphism of Clifford algebras.

        Given a linear map `m : W \to V` (here represented by a matrix
        acting on column vectors), this method returns the algebra
        morphism `Cl(m) : Cl(W, m(Q)) \to Cl(V, Q)`, where `Cl(V, Q)`
        is the Clifford algebra ``self`` and where `m(Q)` is the pullback
        of the quadratic form `Q` to `W`. See the documentation
        of :class:`CliffordAlgebra` for how this pullback and the
        morphism `Cl(m)` are defined.

        .. NOTE::

            This is a map into ``self``.

        INPUT:

        - ``m`` -- a matrix
        - ``names`` -- (default: ``'e'``) the names of the generators of the
          Clifford algebra of the domain of (the map represented by) ``m``

        OUTPUT:

        The algebra morphism `Cl(m)` from `Cl(W, m(Q))` to ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: m = matrix([[1,-1,-1],[0,1,-1],[1,1,1]])
            sage: phi = Cl.lift_module_morphism(m, 'abc')
            sage: phi
            Generic morphism:
              From: The Clifford algebra of the Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 10 17 3 ]
            [ * 11 0 ]
            [ * * 5 ]
              To:   The Clifford algebra of the Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 1 2 3 ]
            [ * 4 5 ]
            [ * * 6 ]
            sage: a,b,c = phi.domain().gens()
            sage: phi(a)
            x + z
            sage: phi(b)
            -x + y + z
            sage: phi(c)
            -x - y + z
            sage: phi(a + 3*b)
            -2*x + 3*y + 4*z
            sage: phi(a) + 3*phi(b)
            -2*x + 3*y + 4*z
            sage: phi(a*b)
            x*y + 2*x*z - y*z + 7
            sage: phi(b*a)
            -x*y - 2*x*z + y*z + 10
            sage: phi(a*b + c)
            x*y + 2*x*z - y*z - x - y + z + 7
            sage: phi(a*b) + phi(c)
            x*y + 2*x*z - y*z - x - y + z + 7

        We check that the map is an algebra morphism::

            sage: phi(a)*phi(b)
            x*y + 2*x*z - y*z + 7
            sage: phi(a*b)
            x*y + 2*x*z - y*z + 7
            sage: phi(a*a)
            10
            sage: phi(a)*phi(a)
            10
            sage: phi(b*a)
            -x*y - 2*x*z + y*z + 10
            sage: phi(b) * phi(a)
            -x*y - 2*x*z + y*z + 10
            sage: phi((a + b)*(a + c)) == phi(a + b) * phi(a + c)
            True

        We can also lift arbitrary linear maps::

            sage: m = matrix([[1,1],[0,1],[1,1]])
            sage: phi = Cl.lift_module_morphism(m, 'ab')
            sage: a,b = phi.domain().gens()
            sage: phi(a)
            x + z
            sage: phi(b)
            x + y + z
            sage: phi(a*b)
            x*y - y*z + 15
            sage: phi(a)*phi(b)
            x*y - y*z + 15
            sage: phi(b*a)
            -x*y + y*z + 12
            sage: phi(b)*phi(a)
            -x*y + y*z + 12

            sage: m = matrix([[1,1,1,2], [0,1,1,1], [0,1,1,1]])
            sage: phi = Cl.lift_module_morphism(m, 'abcd')
            sage: a,b,c,d = phi.domain().gens()
            sage: phi(a)
            x
            sage: phi(b)
            x + y + z
            sage: phi(c)
            x + y + z
            sage: phi(d)
            2*x + y + z
            sage: phi(a*b*c + d*a)
            -x*y - x*z + 21*x + 7
            sage: phi(a*b*c*d)
            21*x*y + 21*x*z + 42
        """
        Q = self._quadratic_form(m)
        # If R is a quadratic form and m is a matrix, then R(m) returns
        # the quadratic form m^t R m.

        if Q == self._quadratic_form and names is None:
            Cl = self
        else:
            Cl = CliffordAlgebra(Q, names)

        n = self._quadratic_form.dim()
        f = lambda x: self.prod(self.sum_of_terms( ((j,), m[j,i]) for j in range(n)
                                                   if m[j,i] != 0 )
                                for i in x)
        return Cl.module_morphism(f, codomain=self,
                                  category=GradedAlgebrasWithBasis(self.base_ring()))

    def lift_isometry(self, m, names=None):
        r"""
        Lift an invertible isometry ``m`` of the quadratric form of
        ``self`` to a Clifford algebra morphism.

        Given an invertible linear map `m : V \to W` (here represented by
        a matrix acting on column vectors), this method returns the
        algebra morphism `Cl(m)` from `Cl(V, Q)` to `Cl(W, m^{-1}(Q))`,
        where `Cl(V, Q)` is the Clifford algebra ``self`` and where
        `m^{-1}(Q)` is the pullback of the quadratic form `Q` to `W` along
        the inverse map `m^{-1} : W \to V`. See the documentation of
        :class:`CliffordAlgebra` for how this pullback and the morphism
        `Cl(m)` are defined.

        INPUT:

        - ``m`` -- an isometry of the quadratic form of ``self``
        - ``names`` -- (default: ``'e'``) the names of the generators of
          the Clifford algebra of the codomain of (the map represented by)
          ``m``

        OUTPUT:

        The algebra morphism `Cl(m)` from ``self`` to `Cl(W, m^{-1}(Q))`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: m = matrix([[1,1,2],[0,1,1],[0,0,1]])
            sage: phi = Cl.lift_isometry(m, 'abc')
            sage: phi(x)
            a
            sage: phi(y)
            a + b
            sage: phi(x*y)
            a*b + 1
            sage: phi(x) * phi(y)
            a*b + 1
            sage: phi(z*y)
            a*b - a*c - b*c
            sage: phi(z) * phi(y)
            a*b - a*c - b*c
            sage: phi(x + z) * phi(y + z) == phi((x + z) * (y + z))
            True
        """
        MS = m.parent()
        if not m.is_invertible():
            raise ValueError('{} is not invertible')
        Q = self._quadratic_form(MS(m.inverse()))

        if Q == self._quadratic_form and names is None:
            Cl = self
        else:
            if names is None:
                names = 'e'
            Cl = CliffordAlgebra(Q, names)

        n = Q.dim()
        f = lambda x: Cl.prod(Cl.sum_of_terms(((j,), m[j,i]) for j in range(n) if m[j,i] != 0)
                              for i in x)
        return self.module_morphism(f, codomain=Cl,
                                    category=GradedAlgebrasWithBasis(self.base_ring()))

    Element = CliffordAlgebraElement

class ExteriorAlgebra(CliffordAlgebra):
    r"""
    An exterior algebra of a free module over a commutative ring.

    Let `V` be a module over a commutative ring `R`. The exterior algebra
    (or Grassmann algebra) `\Lambda(V)` of `V` is defined as the quotient
    of the tensor algebra `T(V)` of `V` modulo the two-sided ideal
    generated by all tensors of the form `x \otimes x` with `x \in V`. The
    multiplication on `\Lambda(V)` is denoted by `\wedge` (so
    `v_1 \wedge v_2 \wedge \cdots \wedge v_n` is the projection of
    `v_1 \otimes v_2 \otimes \cdots \otimes v_n` onto `\Lambda(V)`) and
    called the "exterior product" or "wedge product".

    If `V` is a rank-`n` free `R`-module with a basis
    `\{e_1, \ldots, e_n\}`, then `\Lambda(V)` is the `R`-algebra
    noncommutatively generated by the `n` generators `e_1, \ldots, e_n`
    subject to the relations `e_i^2 = 0` for all `i`, and
    `e_i e_j = - e_j e_i` for all `i < j`. As an `R`-module,
    `\Lambda(V)` then has a basis `(\bigwedge_{i \in I} e_i)` with `I`
    ranging over the subsets of `\{1, 2, \ldots, n\}` (where
    `\bigwedge_{i \in I} e_i` is the wedge product of `e_i` for `i`
    running through all elements of `I` from smallest to largest), and
    hence is free of rank `2^n`.

    The exterior algebra of an `R`-module `V` can also be realized
    as the Clifford algebra of `V` for the quadratic form `Q` given by
    `Q(v) = 0` for all vectors `v \in V`.

    The exterior algebra of an `R`-module `V` is a `\ZZ`-graded connected
    Hopf superalgebra. It is commutative in the super sense (i.e., the
    odd elements anticommute and square to `0`).

    .. WARNING::

        We initialize the exterior algebra as an object of the category
        of Hopf algebras, but this is not really correct, since it is a
        Hopf superalgebra with the odd-degree components forming the odd
        part. So use Hopf-algebraic methods with care!

    .. TODO::

        Add a category for Hopf superalgebras. (Once :trac:`10963`
        is finished...)

    INPUT:

    - ``R`` -- the base ring

    REFERENCES:

    - :wikipedia:`Exterior_algebra`
    """
    @staticmethod
    def __classcall_private__(cls, R, names=None, n=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: E1.<e0,e1,e2> = ExteriorAlgebra(QQ)
            sage: E2 = ExteriorAlgebra(QQ, 3)
            sage: E3 = ExteriorAlgebra(QQ, ['e0','e1','e2'])
            sage: E1 is E2 and E2 is E3
            True
        """
        if names is None:
            names = 'e'
        elif names in ZZ:
            n = names
            names = 'e'

        if isinstance(R, FreeModule_generic):
            if n is not None and n != R.dimension():
                raise ValueError("the number of variables does not match the dimension")
            n = R.dimension()
            R = R.base_ring()

        names = tuple(names)
        if n is not None and len(names) != n:
            if len(names) == 1:
                names = tuple( '{}{}'.format(names[0], i) for i in range(n) )
            else:
                raise ValueError("the number of variables does not match the number of generators")
        return super(ExteriorAlgebra, cls).__classcall__(cls, R, names)

    def __init__(self, R, names):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: TestSuite(E).run()
        """
        CliffordAlgebra.__init__(self, QuadraticForm(R, len(names)), names, GradedHopfAlgebrasWithBasis(R))
        # TestSuite will fail if the HopfAlgebra classes will ever have tests for
        # the coproduct being an algebra morphism -- since this is really a
        # Hopf superalgebra, not a Hopf algebra.

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ExteriorAlgebra(QQ, 3)
            The exterior algebra of rank 3 over Rational Field
        """
        return "The exterior algebra of rank {} over {}".format(self.ngens(), self.base_ring())

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by
        ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E._repr_term((0,1,2))
            'x^y^z'
        """
        if len(m) == 0:
            return '1'
        term = ''
        for i in m:
            if len(term) != 0:
                term += '^'
            term += self.variable_names()[i]
        return term

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E._latex_term((0,1,2))
            ' x \\wedge y \\wedge z'
            sage: E.<x0,x1,x2> = ExteriorAlgebra(QQ)
            sage: E._latex_term((0,1,2))
            ' x_{0} \\wedge x_{1} \\wedge x_{2}'
            sage: E._latex_term(())
            '1'
            sage: E._latex_term((0,))
            ' x_{0}'
        """
        if len(m) == 0:
            return '1'
        term = ''
        for i in m:
            if len(term) != 0:
                term += ' \\wedge'
            term += ' ' + self.latex_variable_names()[i]
        return term

    def lift_morphism(self, phi, names=None):
        r"""
        Lift the matrix ``m`` to an algebra morphism of exterior algebras.

        Given a linear map `\phi : V \to W` (here represented by a matrix
        acting on column vectors over the base ring of `V`), this method
        returns the algebra morphism
        `\Lambda(\phi) : \Lambda(V) \to \Lambda(W)`. This morphism is defined
        on generators `v_i \in \Lambda(V)` by `v_i \mapsto \phi(v_i)`.

        .. NOTE::

            This is the map going out of ``self`` as opposed to
            :meth:`~sage.algebras.clifford_algebra.CliffordAlgebraElement.lift_module_morphism()`
            for general Clifford algebras.

        INPUT:

        - ``phi`` -- a linear map from `V \to W`
        - ``names`` -- (default: ``'e'``) the names of the generators of the
          Clifford algebra of the domain of (the map represented by) ``phi``

        OUTPUT:

        The algebra morphism `\Lambda(\phi)` from ``self`` to `\Lambda(W)`.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: phi = matrix([[0,1],[1,1],[1,2]]); phi
            [0 1]
            [1 1]
            [1 2]
            sage: L = E.lift_morphism(phi, ['a','b','c']); L
            Generic morphism:
              From: The exterior algebra of rank 2 over Rational Field
              To:   The exterior algebra of rank 3 over Rational Field
            sage: L(x)
            b + c
            sage: L(y)
            a + b + 2*c
            sage: L(x*y)
            -a^b - a^c + b^c
            sage: L(x)*L(y)
            -a^b - a^c + b^c
            sage: L(x + y)
            a + 2*b + 3*c
            sage: L(x) + L(y)
            a + 2*b + 3*c
            sage: L(1/2*x + 2)
            1/2*b + 1/2*c + 2
            sage: L(E(3))
            3

            sage: psi = matrix([[1, -3/2]]); psi
            [   1 -3/2]
            sage: Lp = E.lift_morphism(psi, ['a']); Lp
            Generic morphism:
              From: The exterior algebra of rank 2 over Rational Field
              To:   The exterior algebra of rank 1 over Rational Field
            sage: Lp(x)
            a
            sage: Lp(y)
            -3/2*a
            sage: Lp(x + 2*y + 3)
            -2*a + 3
        """
        n = phi.nrows()
        R = self.base_ring()
        E = ExteriorAlgebra(R, names, n)
        f = lambda x: self.prod(self.sum_of_terms( ((j,), phi[j,i]) for j in range(n)
                                                   if phi[j,i] != 0)
                                for i in x)
        return self.module_morphism(f, codomain=E, category=GradedAlgebrasWithBasis(R))

    def volume_form(self):
        """
        Return the volume form of ``self``.

        Given the basis `e_1, e_2, \ldots, e_n` of the underlying
        `R`-module, the volume form is defined as `e_1 \wedge e_2
        \wedge \cdots \wedge e_n`.

        This depends on the choice of basis.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.volume_form()
            x^y^z
        """
        return self.element_class(self, {tuple(range(self.ngens())): self.base_ring().one()})

    def boundary(self, s_coeff):
        r"""
        Return the boundary operator `\partial` defined by the structure
        coefficients ``s_coeff`` of a Lie algebra.

        For more on the boundary operator, see
        :class:`ExteriorAlgebraBoundary`.

        INPUT:

        - ``s_coeff`` -- a dictionary whose keys are in `I \times I`, where
          `I` is the index set of the underlying vector space `V`, and whose
          values can be coerced into 1-forms (degree 1 elements) in ``E``

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.boundary({(0,1): z, (1,2): x, (2,0): y})
            Boundary endomorphism of The exterior algebra of rank 3 over Rational Field
        """
        return ExteriorAlgebraBoundary(self, s_coeff)

    def coboundary(self, s_coeff):
        r"""
        Return the coboundary operator `d` defined by the structure
        coefficients ``s_coeff`` of a Lie algebra.

        For more on the coboundary operator, see
        :class:`ExteriorAlgebraCoboundary`.

        INPUT:

        - ``s_coeff`` -- a dictionary whose keys are in `I \times I`, where
          `I` is the index set of the underlying vector space `V`, and whose
          values can be coerced into 1-forms (degree 1 elements) in ``E``

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            Coboundary endomorphism of The exterior algebra of rank 3 over Rational Field
        """
        return ExteriorAlgebraCoboundary(self, s_coeff)

    def degree_on_basis(self, m):
        r"""
        Return the degree of the monomial indexed by ``m``.

        The degree of ``m`` in the `\ZZ`-grading of ``self`` is defined
        to be the length of ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.degree_on_basis((0,))
            1
            sage: E.degree_on_basis((0,1))
            2
        """
        return ZZ(len(m))

    def coproduct_on_basis(self, a):
        r"""
        Return the coproduct on the basis element indexed by ``a``.

        The coproduct is defined by

        .. MATH::

            \Delta(e_{i_1} \wedge \cdots \wedge e_{i_m}) = \sum_{k=0}^m
            \sum_{\sigma \in Sh_{k,m-k}} (-1)^{\sigma}
            (e_{\sigma(i_1)} \wedge \cdots e_{\sigma(i_k)}) \otimes
            (e_{\sigma(i_{k+1})} \wedge \cdots e_{\sigma(i_m)})

        .. WARNING::

            This coproduct is a homomorphism of superalgebras, not a
            homomorphism of algebras!

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.coproduct_on_basis((0,))
            1 # x + x # 1
            sage: E.coproduct_on_basis((0,1))
            1 # x^y + x # y + x^y # 1 - y # x
        """
        from sage.combinat.words.word import Word
        def shuffle(k):
            sh = Word(a[:k]).shuffle(Word(a[k:]))
            for w in sh:
                descents = [i for i in range(len(w)-1) if w[i] > w[i+1]]
                yield ((tuple(w[:k]), tuple(w[k:])), (-1)**len(descents))
        return self.tensor_square().sum_of_terms([term for k in range(len(a)+1)
                                                  for term in shuffle(k)])

    def antipode_on_basis(self, m):
        r"""
        Return the antipode on the basis element indexed by ``m``.

        Given a basis element `\omega`, the antipode is defined by
        `S(\omega) = (-1)^{\deg(\omega)} \omega`.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.antipode_on_basis((1,))
            -y
            sage: E.antipode_on_basis((1,2))
            y^z
        """
        return self.term(m, (-self.base_ring().one())**len(m))

    def counit(self, x):
        """
        Return the counit of ``x``.

        The counit of a form `\omega` is the constant coefficient.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: elt = x*y - 2*x + 3
            sage: E.counit(elt)
            3
        """
        return x.constant_coefficient()

    def interior_product_on_basis(self, a, b):
        """
        Return the interior product of ``a`` on ``b``.

        This depends on the choice of basis of the vector space
        whose exterior algebra is ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.interior_product_on_basis((0,), (0,))
            1
            sage: E.interior_product_on_basis((0,2), (0,))
            z
            sage: E.interior_product_on_basis((1,), (0,2))
            0
            sage: E.interior_product_on_basis((0,2), (1,))
            0
            sage: E.interior_product_on_basis((0,1,2), (0,2))
            -y
        """
        sgn = 1
        t = list(a)
        for i in b:
            if i not in t:
                return self.zero()
            sgn *= (-1)**t.index(i)
            t.remove(i)
        R = self.base_ring()
        return self.term(tuple(t), R(sgn))

    class Element(CliffordAlgebraElement):
        """
        An element of an exterior algebra.
        """
        def _mul_(self, other):
            """
            Return ``self`` multiplied by ``other``.

            INPUT:

            - ``other`` -- element of the same Clifford algebra as ``self``

            EXAMPLES::

                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: x*y
                x^y
                sage: y*x
                -x^y
                sage: z*y*x
                -x^y^z
                sage: (3*x + y)^2
                0
                sage: (x - 3*y + z/3)^2
                0
            """
            zero = self.parent().base_ring().zero()
            d = {}

            for ml,cl in self:
                for mr,cr in other:
                    # Create the next term
                    t = list(mr)
                    for i in reversed(ml):
                        pos = 0
                        for j in t:
                            if i == j:
                                pos = None
                                break
                            if i < j:
                                break
                            pos += 1
                            cr = -cr
                        if pos is None:
                            t = None
                            break
                        t.insert(pos, i)

                    if t is None: # The next term is 0, move along
                        continue

                    t = tuple(t)
                    d[t] = d.get(t, zero) + cl * cr
                    if d[t] == zero:
                        del d[t]

            return self.__class__(self.parent(), d)

        def interior_product(self, x):
            r"""
            Return the internal product or antiderivation of ``self`` with
            respect to ``x``.

            The *interior product* is a map `i_{\alpha} \colon \Lambda^k(V)
            \to \Lambda^{k-1}(V)`, for a fixed `\alpha \in V^*` (thought of
            as an element in `\Lambda^1(V)`, defined by

            - `i_{\alpha}(v) = \alpha(v)` where `v \in V = \Lambda^1(V)`,
            - it is a graded derivation of degree `-1`:

            .. MATH::

                i_{\alpha}(x \wedge y) = (i_{\alpha} x) \wedge y
                + (-1)^{\deg x} x \wedge (i_{\alpha} y)

            The interior product can also be defined by

            .. MATH::

                (i_{\alpha} \omega)(u_1, \ldots, u_k)
                = \omega(\alpha, u_1, \ldots, u_k),

            where `\omega \in \Lambda^k(V)` is thought of as an
            alternating multilinear mapping from
            `V^* \times \cdots \times V^* \to F` with `F` being the base
            field of `V`.

            EXAMPLES::

                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: x.interior_product(x)
                1
                sage: (x + x*y).interior_product(2*y)
                -2*x
                sage: (x*z + x*y*z).interior_product(2*y - x)
                -2*x^z - y^z - z

            REFERENCES:

            - :wikipedia:`Exterior_algebra#Interior_product`
            """
            P = self.parent()
            return P.sum([c * cx * P.interior_product_on_basis(m, mx)
                          for m,c in self for mx,cx in x])

        antiderivation = interior_product

        def hodge_dual(self):
            r"""
            Return the Hodge dual of ``self``.

            The Hodge dual `\ast` is defined on a basis element `\alpha` by
            `i_{\alpha} \sigma` where `\sigma` is the volume form and
            `i_{\alpha}` denotes the antiderivation function with respect
            to `\alpha`.

            .. NOTE::

                The Hodge dual of the Hodge dual is constant on the `k`-th
                graded part of `\Lambda(V)` up to a sign.

            EXAMPLES::

                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: x.hodge_dual()
                y^z
                sage: (x*z).hodge_dual()
                -y
                sage: (x*y*z).hodge_dual()
                1
                sage: [a.hodge_dual().hodge_dual() for a in E.basis()]
                [1, x, y, z, x^y, x^z, y^z, x^y^z]
                sage: (x + x*y).hodge_dual()
                y^z + z
                sage: (x*z + x*y*z).hodge_dual()
                -y + 1
                sage: E = ExteriorAlgebra(QQ, 'wxyz')
                sage: [a.hodge_dual().hodge_dual() for a in E.basis()]
                [1, -w, -x, -y, -z, w^x, w^y, w^z, x^y, x^z, y^z,
                 -w^x^y, -w^x^z, -w^y^z, -x^y^z, w^x^y^z]
            """
            volume_form = self.parent().volume_form()
            return volume_form.interior_product(self)

        def constant_coefficient(self):
            """
            Return the constant coefficient of ``self``.

            .. TODO::

                Define a similar method for general Clifford algebras once
                the morphism to exterior algebras is implemented.

            EXAMPLES::

                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: elt = 5*x + y + x*z + 10
                sage: elt.constant_coefficient()
                10
                sage: x.constant_coefficient()
                0
            """
            return self._monomial_coefficients.get(self.parent().one_basis(), self.base_ring().zero())

        def scalar(self, other):
            r"""
            Return the Clifford scalar product of ``self`` with ``other``.

            The Clifford scalar (inner) product of `x, y \in Cl(V, Q)` is
            defined by `\langle x, y \rangle = \langle x^t y \rangle` where
            `\langle a \rangle` denotes the degree 0 term of `a`.

            .. TODO::

                Define a similar method for general Clifford algebras once
                the morphism to exterior algebras is implemented.

            EXAMPLES::

                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: elt = 5*x + y + x*z
                sage: elt.scalar(z + 2*x)
                0
                sage: elt.transpose() * (z + 2*x)
                -2*x^y + 5*x^z + y^z
            """
            return (self.transpose() * other).constant_coefficient()

#####################################################################
## Differentials

class ExteriorAlgebraDifferential(ModuleMorphismByLinearity, UniqueRepresentation):
    r"""
    A differential of an exterior algebra `\Lambda(L)` defined by the
    structure coefficients of a Lie algebra `L`.
    """
    @staticmethod
    def __classcall__(cls, E, s_coeff):
        """
        Standardizes the structure coeffcients to ensure a unique
        representation.

        EXAMPLES::

            sage: from sage.algebras.clifford_algebra import ExteriorAlgebraDifferential
            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par1 = ExteriorAlgebraDifferential(E, {(0,1): z, (1,2): x, (2,0): y})
            sage: par2 = ExteriorAlgebraDifferential(E, {(0,1): z, (1,2): x, (0,2): -y})
            sage: par3 = ExteriorAlgebraDifferential(E, {(1,0): {2:-1}, (1,2): {0:1}, (2,0):{1:1}})
            sage: par1 is par2 and par2 is par3
            True
        """
        d = {}

        for k,v in dict(s_coeff).items():
            if isinstance(v, dict):
                R = E.base_ring()
                v = E._from_dict({(i,): R(c) for i,c in v.items()})

            # Make sure all elements are in ``E``
            v = E(v)

            # It's okay if v.degree results in an error
            #   (we'd throw a similar error)
            if v.degree() != 1:
                raise ValueError("elements must be degree 1")

            if k[0] < k[1]:
                d[tuple(k)] = v
            else:
                d[(k[1], k[0])] = -v

        from sage.sets.family import Family
        return super(ExteriorAlgebraDifferential, cls).__classcall__(cls, E, Family(d))

    def __init__(self, E, s_coeff):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2):x, (2,0):y})
            sage: TestSuite(par).run() # known bug - morphisms are properly in a category
        """
        self._s_coeff = s_coeff

        # Technically this preserves the grading but with a shift of -1
        cat = AlgebrasWithBasis(E.base_ring())
        ModuleMorphismByLinearity.__init__(self, domain=E, codomain=E, category=cat)

    def homology(self, deg=None, **kwds):
        """
        Return the homology determined by ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: par.homology()
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 0 over Rational Field,
             2: Vector space of dimension 0 over Rational Field,
             3: Vector space of dimension 1 over Rational Field}
            sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            sage: d.homology()
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 0 over Rational Field,
             2: Vector space of dimension 0 over Rational Field,
             3: Vector space of dimension 1 over Rational Field}
        """
        return self.chain_complex().homology(deg, **kwds)

class ExteriorAlgebraBoundary(ExteriorAlgebraDifferential):
    r"""
    The boundary `\partial` of an exterior algebra `\Lambda(L)` defined
    by the structure coefficients of `L`.

    Let `L` be a Lie algebra. We give an exterior algebra `E = \Lambda(L)`
    a chain complex structure by considering a differential
    `\partial : \Lambda^{k+1}(L) \to \Lambda^k(L)` defined by

    .. MATH::

        \partial(x_1 \wedge x_2 \wedge \cdots \wedge x_{k+1})
        = \sum_{i < j} (-1)^{i+j+1}
        [x_i, x_j] \wedge x_1 \wedge \cdots \wedge \hat{x}_i \wedge \cdots
        \wedge \hat{x}_j \wedge \cdots \wedge x_{k+1}

    where `\hat{x}_i` denotes a missing index. The corresponding homology is
    the Lie algebra homology.

    INPUT:

    - ``E`` -- an exterior algebra
    - ``s_coeff`` -- a dictionary whose keys are in `I \times I`, where
      `I` is the index set of the underlying vector space `V`, and whose
      values can be coerced into 1-forms (degree 1 elements) in ``E``

    EXAMPLES:

    We consider the differential given by Lie algebra given by the cross
    product `\times` of `\RR^3`::

        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
        sage: par(x)
        0
        sage: par(x*y)
        z

    We check that `\partial \circ \partial = 0`::

        sage: p2 = par * par
        sage: all(p2(b) == 0 for b in E.basis())
        True

    REFERENCES:

    - :wikipedia:`Exterior_algebra#Lie_algebra_homology`
    """
    def _repr_type(self):
        """
        TESTS::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: par._repr_type()
            'Boundary'
        """
        return "Boundary"

    def _on_basis(self, m):
        """
        Return the differential on the basis element indexed by ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: par._on_basis(())
            0
            sage: par._on_basis((0,))
            0
            sage: par._on_basis((0,1))
            z
            sage: par._on_basis((0,2))
            -y
            sage: par._on_basis((0,1,2))
            0
        """
        E = self.domain()
        sc = self._s_coeff
        keys = sc.keys()
        return E.sum((-1)**(b+1) * sc[(i,j)]
                      * E.monomial(m[:a] + m[a+1:a+b] + m[a+b+1:])
                     for a,i in enumerate(m) for b,j in enumerate(m[a:]) if (i,j) in keys)

    @cached_method
    def chain_complex(self, R=None):
        """
        Return the chain complex over ``R`` determined by ``self``.

        INPUT:

        - ``R`` -- the base ring; the default is the base ring of
          the exterior algebra

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: par = E.boundary({(0,1): z, (1,2): x, (2,0): y})
            sage: C = par.chain_complex(); C
            Chain complex with at most 4 nonzero terms over Rational Field
            sage: ascii_art(C)
                                      [ 0  0  1]       [0]
                                      [ 0 -1  0]       [0]
                        [0 0 0]       [ 1  0  0]       [0]
             0 <-- C_0 <-------- C_1 <----------- C_2 <---- C_3 <-- 0
        """
        from sage.homology.chain_complex import ChainComplex
        from sage.matrix.constructor import Matrix
        E = self.domain()
        n = E.ngens()
        if R is None:
            R = E.base_ring()

        # Group the basis into degrees
        basis_by_deg = {deg: [] for deg in range(n+1)}
        for b in E.basis().keys():
            basis_by_deg[len(b)].append(b)

        # Construct the transition matrices
        data = {}
        prev_basis = basis_by_deg[0]
        for deg in range(1,n+1):
            # Make sure within each basis we're sorted by lex
            basis = sorted(basis_by_deg[deg])
            mat = []
            for b in basis:
                ret = self._on_basis(b)
                mat.append([ret[p] for p in prev_basis])
            data[deg] = Matrix(mat).transpose().change_ring(R)
            prev_basis = basis

        return ChainComplex(data, degree=-1)

class ExteriorAlgebraCoboundary(ExteriorAlgebraDifferential):
    r"""
    The coboundary `d` of an exterior algebra `\Lambda(V)` defined
    by the structure coefficients of `L`.

    Let `L` be a Lie algebra. We endow its exterior algebra
    `E = \Lambda(L)` with a cochain complex structure by considering a
    differential `d : \Lambda^k(L) \to \Lambda^{k+1}(L)` defined by

    .. MATH::

        d x_i = \sum_{j < k} s_{jk}^i c_j c_k.

    The corresponding homology is the Lie algebra cohomology of `L`.

    This can also be thought of as the exterior derivative, in which case
    the resulting cohomology is the de Rham cohomology of a manifold whose
    exterior algebra of differential forms is ``E``.

    INPUT:

    - ``E`` -- an exterior algebra (understood to be the exterior algebra
      of `L`)
    - ``s_coeff`` -- a dictionary whose keys are in `I \times I`, where
      `I` is the index set of the underlying vector space of `L`, and
      whose values can be coerced into 1-forms (degree 1 elements) in
      ``E``

    EXAMPLES:

    We consider the differential coming from the Lie algebra given by the
    cross product `\times` of `\RR^3`::

        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
        sage: d(x)
        y^z
        sage: d(x*y)
        0

    We check that `d \circ d = 0`::

        sage: d2 = d * d
        sage: all(d2(b) == 0 for b in E.basis())
        True

    REFERENCES:

    - :wikipedia:`Exterior_algebra#Differential_geometry`
    """
    def __init__(self, E, s_coeff):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2):x, (2,0):y})
            sage: TestSuite(d).run() # known bug - morphisms are properly in a category
        """
        # Construct the dictionary of costructure coefficients, i.e. given
        # [x_j, x_k] = \sum_i s_{jk}^i x_i, we get x^i |-> \sum_{j<k} s_{jk}^i x^j x^k
        self._cos_coeff = {}
        zero = E.zero()
        B = E.basis()
        for k,v in dict(s_coeff).items():
            k = B[k]
            for m,c in v:
                self._cos_coeff[m] = self._cos_coeff.get(m, zero) + c * k
        ExteriorAlgebraDifferential.__init__(self, E, s_coeff)

    def _repr_type(self):
        """
        TESTS::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            sage: d._repr_type()
            'Coboundary'
        """
        return "Coboundary"

    def _on_basis(self, m):
        """
        Return the differential on the basis element indexed by ``m``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            sage: d._on_basis(())
            0
            sage: d._on_basis((0,))
            y^z
            sage: d._on_basis((1,))
            -x^z
            sage: d._on_basis((2,))
            x^y
            sage: d._on_basis((0,1))
            0
            sage: d._on_basis((0,2))
            0
            sage: d._on_basis((0,1,2))
            0
        """
        E = self.domain()
        cc = self._cos_coeff
        keys = cc.keys()
        return E.sum((-1)**a * E.monomial(m[:a]) * cc[(i,)] * E.monomial(m[a+1:])
                     for a,i in enumerate(m) if (i,) in keys)

    @cached_method
    def chain_complex(self, R=None):
        """
        Return the chain complex over ``R`` determined by ``self``.

        INPUT:

        - ``R`` -- the base ring; the default is the base ring of
          the exterior algebra

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: d = E.coboundary({(0,1): z, (1,2): x, (2,0): y})
            sage: C = d.chain_complex(); C
            Chain complex with at most 4 nonzero terms over Rational Field
            sage: ascii_art(C)
                                      [ 0  0  1]       [0]
                                      [ 0 -1  0]       [0]
                        [0 0 0]       [ 1  0  0]       [0]
             0 <-- C_3 <-------- C_2 <----------- C_1 <---- C_0 <-- 0
        """
        from sage.homology.chain_complex import ChainComplex
        from sage.matrix.constructor import Matrix
        E = self.domain()
        n = E.ngens()
        if R is None:
            R = E.base_ring()

        # Group the basis into degrees
        basis_by_deg = {deg: [] for deg in range(n+1)}
        for b in E.basis().keys():
            basis_by_deg[len(b)].append(b)

        # Construct the transition matrices
        data = {}
        basis = basis_by_deg[0]
        for deg in range(n):
            # Make sure within each basis we're sorted by lex
            next_basis = sorted(basis_by_deg[deg+1])
            mat = []
            for b in basis:
                ret = self._on_basis(b)
                mat.append([ret[p] for p in next_basis])
            data[deg] = Matrix(mat).transpose().change_ring(R)
            basis = next_basis

        return ChainComplex(data, degree=1)

