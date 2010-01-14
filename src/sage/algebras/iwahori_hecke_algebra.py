"""
Iwahori Hecke Algebras
"""
#*****************************************************************************
#  Copyright (C) 2010 Daniel Bump <bump at match.stanford.edu>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.all import AlgebrasWithBasis, FiniteDimensionalAlgebrasWithBasis, CoxeterGroups
import sage.combinat.root_system.cartan_type
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.structure.element import is_Element
from sage.rings.all import ZZ
from sage.misc.misc import repr_lincomb
from sage.algebras.algebra_element import AlgebraElement
from sage.combinat.family import Family
import sage.rings.polynomial.laurent_polynomial
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.misc.cachefunc import cached_method

class IwahoriHeckeAlgebraT(CombinatorialFreeModule):
    r"""
    INPUT:

     - ``W`` -- A CoxeterGroup or CartanType
     - ``q1`` -- a parameter.

    OPTIONAL ARGUMENTS:

     - ``q2`` -- another parameter (default -1)
     - ``base_ring`` -- A ring containing q1 and q2 (default q1.parent())
     - ``prefix`` -- a label for the generators (default "T")

    The Iwahori Hecke algebra is defined in:

    Nagayoshi Iwahori, On the structure of a Hecke ring of a Chevalley group
    over a finite field.  J. Fac. Sci. Univ. Tokyo Sect. I 10 1964 215--236
    (1964).

    The Iwahori Hecke algebra is a deformation of the group algebra of
    the Weyl/Coxeter group. Taking the deformation parameter `q=1` as in the
    following example gives a ring isomorphic to that group
    algebra. The parameter `q` is a deformation parameter.

    EXAMPLES::

        sage: H = IwahoriHeckeAlgebraT("A3",1,prefix = "s")
        sage: [s1,s2,s3] = H.algebra_generators()
        sage: s1*s2*s3*s1*s2*s1 == s3*s2*s1*s3*s2*s3
        True
        sage: w0 = H(H.coxeter_group().long_element())
        sage: w0
        s1*s2*s3*s1*s2*s1
        sage: H.an_element()
        3*s1*s2 + 3*s1 + 1

    Iwahori Hecke algebras have proved to be fundamental. See for example:

    Kazhdan and Lusztig, Representations of Coxeter groups and Hecke algebras.
    Invent. Math. 53 (1979), no. 2, 165--184.

    Iwahori-Hecke Algebras: Thomas J. Haines, Robert E. Kottwitz,
    Amritanshu Prasad, http://front.math.ucdavis.edu/0309.5168

    V. Jones, Hecke algebra representations of braid groups and link
    polynomials.  Ann. of Math. (2) 126 (1987), no. 2, 335--388.

    For every simple reflection `s_i` of the Coxeter group, there is a
    corresponding generator `T_i` of the Iwahori Hecke algebra. These
    are subject to the relations

        `(T_i-q_1)*(T_i-q_2) == 0`

    together with the braid relations `T_i T_j T_i ... == T_j T_i T_j ...`,
    where the number of terms on both sides is `k/2` with `k` the order of
    `s_i s_j` in the Coxeter group.

    Weyl group elements form a basis of the Iwahori Hecke algebra `H`
    with the property that if `w1` and `w2` are Coxeter group elements
    such that ``(w1*w2).length() == w1.length() + w2.length()`` then
    ``H(w1*w2) == H(w1)*H(w2)``.

    With the default value `q_2 = -1` and with `q_1 = q` the
    generating relation may be written `T_i^2 = (q-1)*T_i + q*1` as in
    Iwahori's paper.

    EXAMPLES::

        sage: R.<q>=PolynomialRing(ZZ)
        sage: H=IwahoriHeckeAlgebraT("B3",q); H
        The Iwahori Hecke Algebra of Type B3 in q,-1 over Univariate Polynomial Ring in q over Integer Ring and prefix T
        sage: T1,T2,T3 = H.algebra_generators()
        sage: T1*T1
        (q-1)*T1 + q

    It is useful to define ``T1,T2,T3 = H.algebra_generators()`` as above
    so that H can parse its own output::

        sage: H(T1)
        T1

    The Iwahori Hecke algebra associated with an affine Weyl group is
    called an affine Hecke algebra. These may be implemented as follows::

        sage: R.<q>=QQ[]
        sage: H=IwahoriHeckeAlgebraT(['A',2,1],q)
        sage: [T0,T1,T2]=H.algebra_generators()
        sage: T1*T2*T1*T0*T1*T0
        (q-1)*T1*T2*T0*T1*T0 + q*T1*T2*T0*T1
        sage: T1*T2*T1*T0*T1*T1
        q*T1*T2*T1*T0 + (q-1)*T1*T2*T0*T1*T0
        sage: T1*T2*T1*T0*T1*T2
        T1*T2*T0*T1*T0*T2
        sage: (T1*T2*T0*T1*T0*T2).support_of_term() # get the underlying Weyl group element
        [ 2  1 -2]
        [ 3  1 -3]
        [ 2  2 -3]

        sage: R = IwahoriHeckeAlgebraT("A3",0,0,prefix = "s")
        sage: [s1,s2,s3] = R.algebra_generators()
        sage: s1*s1
        0

    TESTS::

        sage: H1 = IwahoriHeckeAlgebraT("A2",1)
        sage: H2 = IwahoriHeckeAlgebraT("A2",1)
        sage: H3 = IwahoriHeckeAlgebraT("A2",-1)
        sage: H1 == H1, H1 == H2, H1 is H2
        (True, True, True)
        sage: H1 == H3
        False

        sage: R.<q>=PolynomialRing(QQ)
        sage: IwahoriHeckeAlgebraT("A3",q).base_ring() == R
        True

        sage: R.<q>=QQ[]; H=IwahoriHeckeAlgebraT("A2",q)
        sage: 1+H(q)
        (q+1)

        sage: R.<q>=QQ[]
        sage: H = IwahoriHeckeAlgebraT("A2",q)
        sage: T1,T2 = H.algebra_generators()
        sage: H(T1+2*T2)
        T1 + 2*T2

    .. rubric:: Design discussion

    This is a preliminary implementation. For work in progress, see:
    http://wiki.sagemath.org/HeckeAlgebras.

     - Should we use q in QQ['q'] as default parameter for q_1?

    """

    @staticmethod
    def __classcall__(cls, W, q1, q2=-1, base_ring=None, prefix="T"):
        """
        TESTS::

            sage: H = IwahoriHeckeAlgebraT("A2", 1)
            sage: H.coxeter_group() == WeylGroup("A2")
            True
            sage: H.cartan_type() == CartanType("A2")
            True
            sage: H.base_ring() == ZZ
            True
            sage: H._q2 == -1
            True
        """
        if W not in CoxeterGroups():
            W = WeylGroup(W)
        if base_ring is None:
            base_ring = q1.parent()
        q2 = base_ring(q2)
        return super(IwahoriHeckeAlgebraT, cls).__classcall__(cls, W, q1=q1, q2=q2, base_ring=base_ring, prefix=prefix)

    def __init__(self, W, q1, q2, base_ring, prefix):
         """
        EXAMPLES ::

            sage: R.<q1,q2>=QQ[]
            sage: H = IwahoriHeckeAlgebraT("A2",q1,q2,base_ring=Frac(R),prefix="t"); H
            The Iwahori Hecke Algebra of Type A2 in q1,q2 over Fraction Field of Multivariate Polynomial Ring in q1, q2 over Rational Field and prefix t
            sage: TestSuite(H).run()

         """
         self._cartan_type = W.cartan_type()
         self._prefix = prefix
         self._index_set = W.index_set()
         self._q1 = base_ring(q1)
         self._q2 = base_ring(q2)

         if W.is_finite():
             category = FiniteDimensionalAlgebrasWithBasis(base_ring)
         else:
             category = AlgebrasWithBasis(base_ring)
         CombinatorialFreeModule.__init__(self, base_ring, W, category = category)

    def _element_constructor_(self, w):
        """
        Construct a basis element from an element of the Weyl group

        EXAMPLES ::

            sage: R.<q>=QQ[]
            sage: H = IwahoriHeckeAlgebraT("A2",q)
            sage: [H(x) for x in H.coxeter_group()]   # indirect doctest
            [1, T1, T1*T2, T1*T2*T1, T2, T2*T1]

        """
        assert w in self.basis().keys()
        return self.monomial(w)

    @cached_method
    def one_basis(self):
        """
        Returns the unit of the underlying Coxeter group, which indexes
        the one of this algebra, as per
        :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: H = IwahoriHeckeAlgebraT("B3", 1)
            sage: H.one_basis()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: H.one_basis() == H.coxeter_group().one()
            True
            sage: H.one()
            1
        """
        return self.coxeter_group().one()

    def _repr_(self):
        """
        EXAMPLES ::

            sage: R.<q1,q2>=QQ[]
            sage: IwahoriHeckeAlgebraT("A2",q1,-q2,base_ring=Frac(R),prefix="Z") # indirect doctest
            The Iwahori Hecke Algebra of Type A2 in q1,-q2 over Fraction Field of Multivariate Polynomial Ring in q1, q2 over Rational Field and prefix Z

        """
        return "The Iwahori Hecke Algebra of Type %s in %s,%s over %s and prefix %s"%(self._cartan_type._repr_(compact=True), self._q1, self._q2, self.base_ring(), self._prefix)

    def cartan_type(self):
        """
        EXAMPLES ::

            sage: IwahoriHeckeAlgebraT("D4",0).cartan_type()
            ['D', 4]

        """
        return self._cartan_type

    def coxeter_group(self):
        """
        EXAMPLES::

            sage: IwahoriHeckeAlgebraT("B2",1).coxeter_group()
            Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)

        """
        return self.basis().keys()

    def index_set(self):
        """
        EXAMPLES::

            sage: IwahoriHeckeAlgebraT("B2",1).index_set()
            [1, 2]
        """
        return self._index_set

    @cached_method
    def algebra_generators(self):
        """
        Returns the generators. They do not have order two but satisfy
        a quadratic relation. They coincide with the simple
        reflections in the Coxeter group when `q_1=1` and `q_2=-1`. In
        this special case, the Iwahori Hecke algebra is identified
        with the group algebra of the Coxeter group.

        EXAMPLES ::

            sage: R.<q>=QQ[]
            sage: H = IwahoriHeckeAlgebraT("A3",q)
            sage: T=H.algebra_generators(); T
            Finite family {1: T1, 2: T2, 3: T3}
            sage: T.list()
            [T1, T2, T3]
            sage: [T[i] for i in [1,2,3]]
            [T1, T2, T3]
            sage: [T1, T2, T3] = H.algebra_generators()
            sage: T1
            T1
            sage: H = IwahoriHeckeAlgebraT(['A',2,1],q)
            sage: T=H.algebra_generators(); T
            Finite family {0: T0, 1: T1, 2: T2}
            sage: T.list()
            [T0, T1, T2]
            sage: [T[i] for i in [0,1,2]]
            [T0, T1, T2]
            sage: [T0, T1, T2] = H.algebra_generators()
            sage: T0
            T0
        """
        return self.coxeter_group().simple_reflections().map(self.monomial)

    def algebra_generator(self, i):
        """
        EXAMPLES ::

            sage: R.<q>=QQ[]
            sage: H = IwahoriHeckeAlgebraT("A3",q)
            sage: [H.algebra_generator(i) for i in H.index_set()]
            [T1, T2, T3]

        """
        return self.algebra_generators()[i]

    def inverse_generator(self, i):
        """
        This method is only available if q1 and q2 are invertible. In
        that case, the algebra generators are also invertible and this
        method returns the inverse of self.algebra_generator(i). The
        base ring should be either a field or a Laurent polynomial ring.

        EXAMPLES::

            sage: P.<q1,q2>=QQ[]
            sage: F=Frac(P)
            sage: H = IwahoriHeckeAlgebraT("A2",q1,q2,base_ring=F)
            sage: H.base_ring()
            Fraction Field of Multivariate Polynomial Ring in q1, q2 over Rational Field
            sage: H.inverse_generator(1)
            ((-1)/(q1*q2))*T1 + ((q1+q2)/(q1*q2))
            sage: H = IwahoriHeckeAlgebraT("A2",q1,-1,base_ring=F)
            sage: H.inverse_generator(2)
            ((-1)/(-q1))*T2 + ((q1-1)/(-q1))
            sage: P1.<r1,r2>=LaurentPolynomialRing(QQ)
            sage: H1 = IwahoriHeckeAlgebraT("B2",r1,r2,base_ring=P1)
            sage: H1.base_ring()
            Multivariate Laurent Polynomial Ring in r1, r2 over Rational Field
            sage: H1.inverse_generator(2)
            (-r1^-1*r2^-1)*T2 + (r2^-1+r1^-1)
            sage: H2 = IwahoriHeckeAlgebraT("C2",r1,-1,base_ring=P1)
            sage: H2.inverse_generator(2)
            (r1^-1)*T2 + (-1+r1^-1)

        """
        try:
            # This currently works better than ~(self._q1) if
            # self.base_ring() is a Laurent polynomial ring since it
            # avoids accidental coercion into a field of fractions.
            i1 = self._q1.__pow__(-1)
            i2 = self._q2.__pow__(-1)
        except:
            raise ValueError, "%s and %s must be invertible."%(self._q1, self._q2)
        return (-i1*i2)*self.algebra_generator(i)+(i1+i2)

    @cached_method
    def inverse_generators(self):
        """
        This method is only available if q1 and q2 are invertible. In
        that case, the algebra generators are also invertible and this
        method returns their inverses.

        EXAMPLES ::
            sage: P.<q> = PolynomialRing(QQ)
            sage: F = Frac(P)
            sage: H = IwahoriHeckeAlgebraT("A2",q,base_ring=F)
            sage: [T1,T2]=H.algebra_generators()
            sage: [U1,U2]=H.inverse_generators()
            sage: U1*T1,T1*U1
            (1, 1)
            sage: P1.<q> = LaurentPolynomialRing(QQ)
            sage: H1 = IwahoriHeckeAlgebraT("A2",q,base_ring=P1,prefix="V")
            sage: [V1,V2]=H1.algebra_generators()
            sage: [W1,W2]=H1.inverse_generators()
            sage: [W1,W2]
            [(q^-1)*V1 + (-1+q^-1), (q^-1)*V2 + (-1+q^-1)]
            sage: V1*W1, W2*V2
            (1, 1)

        """
        return Family(self.index_set(), self.inverse_generator)

    def product_on_basis(self, w1, w2):
        """

        Returns `T_w1 T_w2`, where `w_1` and `w_2` are in the Coxeter group

        EXAMPLES::

            sage: R.<q>=QQ[]; H = IwahoriHeckeAlgebraT("A2",q)
            sage: s1,s2 = H.coxeter_group().simple_reflections()
            sage: [H.product_on_basis(s1,x) for x in [s1,s2]]
            [(q-1)*T1 + q, T1*T2]

        """
        result = self.monomial(w1)
        for i in w2.reduced_word():
            result = self.product_by_generator(result, i)
        return result

    def product_by_generator_on_basis(self, w, i, side = "right"):
        """
        INPUT:
         - ``w`` - an element of the Coxeter group
         - ``i`` - an element of the index set
         - ``side`` - "left" or "right" (default: "right")

        Returns the product `T_w T_i` (resp. `T_i T_w`) if ``side`` is "right" (resp. "left")

        EXAMPLES::

            sage: R.<q>=QQ[]; H = IwahoriHeckeAlgebraT("A2",q)
            sage: s1,s2 = H.coxeter_group().simple_reflections()
            sage: [H.product_by_generator_on_basis(w, 1) for w in [s1,s2,s1*s2]]
            [(q-1)*T1 + q, T2*T1, T1*T2*T1]
            sage: [H.product_by_generator_on_basis(w, 1, side = "left") for w in [s1,s2,s1*s2]]
            [(q-1)*T1 + q, T1*T2, (q-1)*T1*T2 + q*T2]
        """
        wi = w.apply_simple_reflection(i, side = side)
        if w.has_descent(i, side = side):
            return self.term(w ,  self._q1+self._q2) + \
                   self.term(wi, -self._q1*self._q2)
        else:
            return self.monomial(wi)

    def product_by_generator(self, x, i, side = "right"):
        """
        Returns T_i*x, where T_i is the i-th generator. This is coded
        individually for use in x._mul_().

        EXAMPLES ::
            sage: R.<q>=QQ[]; H = IwahoriHeckeAlgebraT("A2",q)
            sage: [T1,T2] = H.algebra_generators()
            sage: [H.product_by_generator(x, 1, side = "left") for x in [T1,T2]]
            [(q-1)*T1 + q, T2*T1]

        """
        return self.sum(self.product_by_generator_on_basis(w, i)._acted_upon_(c) for (w,c) in x)

    def _repr_term(self, t):
        """
        EXAMPLES ::

            sage: R.<q>=QQ[]
            sage: H = IwahoriHeckeAlgebraT("A3",q)
            sage: W=H.coxeter_group()
            sage: H._repr_term(W.from_reduced_word([1,2,3]))
            'T1*T2*T3'

        """
        redword = t.reduced_word()
        if len(redword) == 0:
            return "1"
        else:
            return "*".join("%s%d"%(self._prefix, i) for i in redword)

    class Element(CombinatorialFreeModuleElement):
        """
        A class for elements of an IwahoriHeckeAlgebra

        TESTS::

            sage: R.<q>=QQ[]
            sage: H=IwahoriHeckeAlgebraT("B3",q)
            sage: [T1, T2, T3] = H.algebra_generators()
            sage: T1+2*T2*T3
            T1 + 2*T2*T3

            sage: R.<q1,q2>=QQ[]
            sage: H=IwahoriHeckeAlgebraT("A2",q1,q2,prefix="x")
            sage: sum(H.algebra_generators())^2
            x1*x2 + x2*x1 + (q1+q2)*x1 + (q1+q2)*x2 + (-2*q1*q2)

            sage: H=IwahoriHeckeAlgebraT("A2",q1,q2,prefix="t")
            sage: [t1,t2] = H.algebra_generators()
            sage: (t1-t2)^3
            (q1^2-q1*q2+q2^2)*t1 + (-q1^2+q1*q2-q2^2)*t2

            sage: R.<q>=QQ[]
            sage: H=IwahoriHeckeAlgebraT("G2",q)
            sage: [T1, T2] = H.algebra_generators()
            sage: T1*T2*T1*T2*T1*T2 == T2*T1*T2*T1*T2*T1
            True
            sage: T1*T2*T1 == T2*T1*T2
            False

            sage: H = IwahoriHeckeAlgebraT("A2",1)
            sage: [T1,T2]=H.algebra_generators()
            sage: T1+T2
            T1 + T2

            sage: -(T1+T2)
            -T1 - T2

            sage: 1-T1
            -T1 + 1

            sage: T1.parent()
            The Iwahori Hecke Algebra of Type A2 in 1,-1 over Integer Ring and prefix T
        """

        def inverse(self):
            """
            If the element is a basis element, that is, an element of the
            form self(w) with w in the Weyl group, this method computes
            its inverse. The base ring must be a field or Laurent
            polynomial ring. Other elements of the ring have inverses but
            the inverse method is only implemented for the basis elements.

            EXAMPLES::

                sage: R.<q>=LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebraT("A2",q)
                sage: [T1,T2]=H.algebra_generators()
                sage: x = (T1*T2).inverse(); x
                (q^-2)*T2*T1 + (-q^-1+q^-2)*T1 + (-q^-1+q^-2)*T2 + (1-2*q^-1+q^-2)
                sage: x*T1*T2
                1

            """
            if len(self) != 1:
                raise NotImplementedError, "inverse only implemented for basis elements (monomials in the generators)"%self
            H = self.parent()
            w = self.support_of_term()

            return H.prod(H.inverse_generator(i) for i in reversed(w.reduced_word()))
