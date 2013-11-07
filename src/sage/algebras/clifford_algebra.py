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
from copy import copy

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.rings.all import ZZ
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
            -2*x*y*z + 5*x*z + 5*x + 2*y + 2*z - 31
        """
        return repr_from_monomials(self.list(), self.parent()._repr_term)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        TESTS::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: latex( ((x^3-z)*x + y)^2 )
            -2  x y z + 5  x z + 5  x + 2  y + 2  z - 31
            sage: Cl.<x0,x1,x2> = CliffordAlgebra(Q)
            sage: latex(  (x1 - x2)*x0 + 5*x0*x1*x2 )
            5  x_{0} x_{1} x_{2} -  x_{0} x_{1} +  x_{0} x_{2} - 1
        """
        return repr_from_monomials(self.list(), self.parent()._latex_term, True)

    def _mul_(self, other):
        """
        Return ``self`` multiplied by ``other``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: (x^3 - z*y)*x*(y*z + x*y*z)
            4*x*y*z + 4*y*z - 96*x - 24*y - 34*z + 192
        """
        Q = self.parent()._quadratic_form
        d = {}

        # Create the standard basis vectors for simplicity
        e = []
        for i in range(Q.dim()):
            e.append([0]*Q.dim())
            e[-1][i] = 1

        for ml,cl in self._monomial_coefficients.items():
            # Distribute the current term over the other
            cur = copy(other._monomial_coefficients) # The current distribution of the term
            for i in reversed(ml):
                # Distribute the current factor
                next = {}
                for mr,cr in cur.items():
                    # Commute the factor as necessary until we are in order
                    pos = 0
                    for j in mr:
                        if i < j:
                            break
                        # Add the additional term from the commutation
                        t = list(mr)
                        t.pop(pos)
                        t = tuple(t)
                        uv = [0] * Q.dim()
                        uv[i] = uv[j] = 1
                        next[t] = next.get(t, 0) + cr * (Q(uv) - Q(e[i]) - Q(e[j]))
                        cr = -cr
                        if next[t] == 0:
                            del next[t]
                        pos += 1

                    # Check to see if we have a squared term or not
                    t = list(mr)
                    if i in t:
                        t.remove(i)
                        cr *= Q(e[i])
                    else:
                        t.insert(pos, i)
                        t.sort()
                    t = tuple(t)
                    next[t] = next.get(t, 0) + cr
                    if next[t] == 0:
                        del next[t]
                cur = next

            # Add the distributed terms to the total
            for index,coeff in cur.items():
                d[index] = d.get(index, 0) + cl * coeff
                if d[index] == 0:
                    del d[index]
        return self.__class__(self.parent(), d)

    def list(self):
        """
        Return ``self`` as a list.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y
            sage: elt.list()
            [((0,), 5), ((1,), 1)]
        """
        return sorted(self._monomial_coefficients.items(), key=lambda x: (-len(x[0]), x[0]))

    def support(self):
        """
        Return the support of ``self``.

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
        by

        .. MATH::

            x_1 \wedge x_2 \wedge \cdots \wedge x_m \mapsto
            (-1)^m x_1 \wedge x_2 \wedge \cdots \wedge x_m.

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
        return self.__class__(self.parent(), {m: (-1)**len(m)*c for m,c in self._monomial_coefficients.items()})

    def transpose(self):
        r"""
        Return the transpose of ``self``.

        The transpose is an antilinear involution of a Clifford algebra and
        is defined on basis elements:

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
        if len(self._monomial_coefficients) == 0:
            return self.parent().zero()
        g = self.parent().gens()
        return sum(c * self.parent().prod(g[i] for i in reversed(m))
                   for m,c in self._monomial_coefficients.items())

    def conjugate(self):
        r"""
        Return the Clifford conjugate of ``self``.

        The Clifford conjugate of `x` is defined by:

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

    def constant_coefficient(self):
        """
        Return the constant coefficient of ``self``.

        .. TODO::

            This isn't functorial in the underlying quadratic space.
            Just reordering the basis is enough to break
            functoriality, as witnessed by the following::

                sage: Q = QuadraticForm(ZZ, 2, [1,1,2])
                sage: R = QuadraticForm(ZZ, 2, [2,1,1])    # isomorphic to Q by switching variables
                sage: ClQ = CliffordAlgebra(Q, ('x','y'))
                sage: ClR = CliffordAlgebra(R, ('x','y'))    # isomorphic to ClQ by sending x to y, y to x
                sage: ClQ.inject_variables()
                Defining x, y
                sage: (x*y).constant_coefficient()
                0
                sage: ClR.inject_variables()
                Defining x, y
                sage: (y*x).constant_coefficient()    # y*x is the image of x*y under iso induced by switching basis elements
                1

            The implementation is good for exterior algebras, though.
            I'd suggest moving it there. There *is* a reasonable
            way to define a "constant coefficient" if `2` is
            invertible or, more generally, if we have a bilinear form
            extending the quadratic form, but we might just as well
            leave this for a later patch. -- dg

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
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

        The Clifford scalar (inner) product of `x, y \in Cl(V, Q)` is defined
        by `\langle x, y \rangle = \langle x^t y \rangle` where
        `\langle a \rangle` denotes the constant term of `a`.

        .. TODO::

            Pretty sure this is incorrect since `constant_coefficient` is.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1, 2, 3, 4, 5, 6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: elt = 5*x + y + x*z
            sage: elt.scalar(z + 2*x)
            -16
            sage: elt.transpose() * (z + 2*x)
            -2*x*y + 5*x*z + y*z + 12*x - z - 16
        """
        return (self.transpose() * other).constant_coefficient()

class CliffordAlgebra(CombinatorialFreeModule):
    r"""
    The Clifford algebra of a quadratic form.

    REFERENCES:

    - :wikipedia:`Clifford_algebra`

    INPUT:

    - ``Q`` -- a quadratic form
    - ``names`` -- (default: ``'e'``) the generator names

    EXAMPLES:

    To create a Clifford algebra, all one needs to do is specify a quadratic
    form::

        sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
        sage: Cl = CliffordAlgebra(Q)
        sage: Cl
        The Clifford algebra of the Quadratic form in 3 variables over Integer Ring with coefficients: 
        [ 1 2 3 ]
        [ * 4 5 ]
        [ * * 6 ]

    We can also explicitly name the generators. In this example, we also
    construct a Clifford algebra isomorphic to an exterior algebra::

        sage: Q = QuadraticForm(ZZ, 4, [0]*10)
        sage: Cl.<a,b,c,d> = CliffordAlgebra(Q)
        sage: a*d
        a*d
        sage: d*c*b*a + a + 4*b*c
        a*b*c*d + 4*b*c + a
    """
    @staticmethod
    def __classcall_private__(cls, Q, names='e'):
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
        """
        self._quadratic_form = Q
        R = Q.base_ring()
        if category is None:
            category = GradedAlgebrasWithBasis(R)
        indices = map( tuple, Subsets(range(Q.dim())) )
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

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: [Cl.gen(i) for i in range(3)]
            [x, y, z]
        """
        return self._from_dict({(i,): self.base_ring().one()}, remove_zeros=False)

    def gens(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.gens()
            (x, y, z)
        """
        return tuple(self.gen(i) for i in range(self.ngens()))

    algebra_generators = gens

    def ngens(self):
        """
        Return the number of generators of ``self``.

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
        Return the degree of the monomial index by ``m``.

        .. WARNING::

            This is a degree in the sense of super-algebra, i. e., it
            distinguishes "even" from "odd" elements. It is not a
            `\ZZ`-valued degree. In general, Clifford algebras are not
            `\ZZ`-graded.

        .. TODO::

            As far as I understand, this method spills over to the
            exterior algebra, which IS `\ZZ`-graded, and provides
            unexpected results there. Should we do anything about it?
            -- dg

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Cl.<x,y,z> = CliffordAlgebra(Q)
            sage: Cl.degree_on_basis((0,))
            1
            sage: Cl.degree_on_basis((0,1))
            0
        """
        return len(m) % ZZ(2)

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

    Element = CliffordAlgebraElement

class ExteriorAlgebra(CliffordAlgebra):
    r"""
    An exterior algebra of a free module over a commutative ring.

    Let `V` be a module over a commutative ring `R`. The exterior algebra
    (or Grassmann algebra) `\Lambda(V)` is defined as the quotient of
    tensor algebra `T(V)` of `V` modulo the two-sided ideal generated by
    all tensors of the form `x \otimes x` with `x \in V`. The
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
    `I \subseteq \{e_1, \ldots e_n\}` (where `\bigwedge_{i \in I} e_i`
    is the wedge product of all elements of `I` from smallest to
    largest), and hence is free of rank `2^n`.

    The exterior algebra of an `R`-module `V` can also be realized
    as the Clifford algebra of `V` and the quadratic form `Q(v) = 0`
    for all vectors `v \in V`.

    The exterior algebra of an `R`-module `V` is a graded connected
    Hopf superalgebra.

    INPUT:

    - ``R`` -- the base ring

    REFERENCES:

    - :wikipedia:`Exterior_algebra`
    """
    @staticmethod
    def __classcall_private__(cls, R, names='e', n=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: E1.<e0,e1,e2> = ExteriorAlgebra(QQ)
            sage: E2 = ExteriorAlgebra(QQ, 3)
            sage: E3 = ExteriorAlgebra(QQ, ['e0','e1','e2'])
            sage: E1 is E2 and E2 is E3
            True
        """
        if names in ZZ:
            n = names
            names = 'e'
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
        CliffordAlgebra.__init__(self, QuadraticForm(R, len(names)), names, HopfAlgebrasWithBasis(R))
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
        Return a string representation of the basis element indexed by ``m``.

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
        """
        if len(m) == 0:
            return '1'
        term = ''
        for i in m:
            if len(term) != 0:
                term += ' \\wedge'
            term += ' ' + self.latex_variable_names()[i]
        return term

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

    def internal_product_on_basis(self, a, b):
        """
        Return the internal product of ``a`` on ``b``.

        This depends on the choice of basis of the vector space
        whose exterior algebra is ``self``.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: E.internal_product_on_basis((0,), (0,))
            1
            sage: E.internal_product_on_basis((0,2), (0,))
            z
            sage: E.internal_product_on_basis((1,), (0,2))
            0
            sage: E.internal_product_on_basis((0,2), (1,))
            0
            sage: E.internal_product_on_basis((0,1,2), (0,2))
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
            d = {}

            for ml,cl in self._monomial_coefficients.items():
                for mr,cr in other._monomial_coefficients.items():
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
                    d[t] = d.get(t, 0) + cl * cr
                    if d[t] == 0:
                        del d[t]

            return self.__class__(self.parent(), d)

        def internal_product(self, x):
            """
            Return the internal product or antiderivation of ``self`` with
            respect to ``x``.

            EXAMPLES::

                sage: E.<x,y,z> = ExteriorAlgebra(QQ)
                sage: x.internal_product(x)
                1
                sage: (x + x*y).internal_product(2*y)
                -2*x
                sage: (x*z + x*y*z).internal_product(2*y - x)
                -2*x^z - y^z - z
            """
            P = self.parent()
            return P.sum([c * cx * P.internal_product_on_basis(m, mx)
                          for m,c in self._monomial_coefficients.items()
                          for mx,cx in x._monomial_coefficients.items()])

        antiderivation = internal_product

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
            return volume_form.internal_product(self)

