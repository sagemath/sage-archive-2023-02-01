r"""
Group homomorphisms for groups with a libGAP backend

EXAMPLES::

    sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
    sage: A = AbelianGroupGap([2, 4])
    sage: F.<a,b> = FreeGroup()
    sage: f = F.hom([g for g in A.gens()])
    sage: K = f.kernel()
    sage: K
    Group(<free, no generators known>)

AUTHORS:

- Simon Brandhorst (2018-02-08): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import HomsetWithBase
from sage.categories.groups import Groups
from sage.categories.morphism import Morphism
from sage.rings.integer_ring import ZZ
from sage.groups.libgap_wrapper import ParentLibGAP
from sage.misc.latex import latex

class GroupMorphism_libgap(Morphism):
    r"""
    Group morphism specified by the images of generators.

    This wraps GAP's ``GroupHomomorphismByImages`` function.
    Checking if the input defines a group homomorphism can be expensive
    if the group is large.

    INPUT:

    - ``homset`` -- the parent
    - ``imgs`` -- a tuple of generators
    - ``check`` -- (default: ``True``) check if the images define
      a group homomorphism

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2, 4])
        sage: A.hom([g^2 for g in A.gens()])
        Group endomorphism of Abelian group with gap, generator orders (2, 4)

    Homomorphisms can be defined between different kinds of libGAP groups::

        sage: G = MatrixGroup([Matrix(ZZ, 2, [0,1,1,0])])
        sage: f = A.hom([G.0, G(1)])
        sage: f
        Group morphism:
        From: Abelian group with gap, generator orders (2, 4)
        To:   Matrix group over Integer Ring with 1 generators (
        [0 1]
        [1 0]
        )
        sage: G.<a,b> = FreeGroup()
        sage: H = G / (G([1]), G([2])^3)
        sage: f = G.hom(H.gens())
        sage: f
        Group morphism:
          From: Free Group on generators {a, b}
          To:   Finitely presented group < a, b | a, b^3 >

    TESTS:

    Old tests inherited from MatrixGroupMorphisms::

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: phi = G.hom(H.gens())
        sage: phi
        Group morphism:
        From: Matrix group over Finite Field of size 5 with 1 generators (
        [1 1]
        [0 1]
        )
        To:   Matrix group over Finite Field of size 5 with 1 generators (
        [1 0]
        [1 1]
        )
        sage: phi(MS([1,1,0,1]))
        [1 0]
        [1 1]
        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])

    Check that :trac:`19406` is fixed::

        sage: G = GL(2, GF(3))
        sage: H = GL(3, GF(2))
        sage: mat1 = H([[-1,0,0],[0,0,-1],[0,-1,0]])
        sage: mat2 = H([[1,1,1],[0,0,-1],[-1,0,0]])
        sage: phi = G.hom([mat1, mat2])
        Traceback (most recent call last):
        ...
        ValueError: images do not define a group homomorphism

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: phi = G.hom(H.gens())
        sage: phi.gap()
        CompositionMapping( [ (6,7,8,10,9)(11,13,14,12,15)(16,19,20,18,17)(21,25,22,24,23) ]
        -> [ [ [ Z(5)^0, 0*Z(5) ], [ Z(5)^0, Z(5)^0 ] ] ], <action isomorphism> )
        sage: type(_)
        <type 'sage.libs.gap.element.GapElement'>

        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])
        sage: phi.kernel()
        Matrix group over Finite Field of size 7 with 1 generators (
        [6 0]
        [0 1]
        )

        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])
        sage: phi
        Group endomorphism of Matrix group over Finite Field of size 7 with 1 generators (
        [3 0]
        [0 1]
        )
        sage: g = G.gens()[0]
        sage: phi(g)
        [2 0]
        [0 1]
        sage: H = MatrixGroup([MS(a.list())])
        sage: H
        Matrix group over Finite Field of size 7 with 1 generators (
        [2 0]
        [0 1]
        )

    The following tests against :trac:`10659`::

        sage: phi(H)   # indirect doctest
        Matrix group over Finite Field of size 7 with 1 generators (
        [4 0]
        [0 1]
        )
        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: g = MS([1,1,0,1])
        sage: G = MatrixGroup([g])
        sage: phi = G.hom(G.gens())
        sage: phi(G.0)
        [1 1]
        [0 1]
        sage: phi(G(g^2))
        [1 2]
        [0 1]

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: gens = [MS([1,2,  -1,1]),MS([1,1,  0,1])]
        sage: G = MatrixGroup(gens)
        sage: phi = G.hom(G.gens())
        sage: phi(G.0)
        [1 2]
        [4 1]
        sage: phi(G.1)
        [1 1]
        [0 1]

    We check that :trac:`19780` is fixed::

        sage: G = groups.matrix.SO(3, 3)
        sage: H = groups.matrix.GL(3, 3)
        sage: phi = G.hom([H(x) for x in G.gens()])
        sage: phi(G.one()).parent()
        General Linear Group of degree 3 over Finite Field of size 3

        sage: MS = MatrixSpace(SR, 2, 2)
        sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
        sage: G.Hom(G)
        Set of Morphisms from Matrix group over Symbolic Ring with 2 generators (
        [1 0]  [1 2]
        [0 1], [3 4]
        ) to Matrix group over Symbolic Ring with 2 generators (
        [1 0]  [1 2]
        [0 1], [3 4]
        ) in Category of groups

        sage: G = MatrixGroup([matrix(GF(5), [[1,3],[0,1]])])
        sage: H = MatrixGroup([matrix(GF(5), [[1,2],[0,1]])])
        sage: G.hom([H.gen(0)])
        Group morphism:
        From: Matrix group over Finite Field of size 5 with 1 generators (
        [1 3]
        [0 1]
        )
        To:   Matrix group over Finite Field of size 5 with 1 generators (
        [1 2]
        [0 1]
        )

        sage: G = GO(2,2,e=1)
        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2])
        sage: G.hom([A.one(),A.gen(0)])
        Group morphism:
        From: General Orthogonal Group of degree 2 and form parameter 1 over Finite Field of size 2
        To:   Abelian group with gap, generator orders (2,)

    Check that :trac:`19407` is fixed::

        sage: G = GL(2, GF(2))
        sage: H = GL(3, ZZ)
        sage: Hom(G, H)
        Set of Morphisms from General Linear Group of degree 2
         over Finite Field of size 2
         to General Linear Group of degree 3 over Integer Ring in Category of groups
    """
    def __init__(self, homset, imgs, check=True):
        r"""
        Constructor method.

        TESTS:

        The following input does not define a valid group homomorphism::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2])
            sage: G = MatrixGroup([Matrix(ZZ, 2, [0,1,-1,0])])
            sage: A.hom([G.gen(0)])
            Traceback (most recent call last):
            ...
            ValueError: images do not define a group homomorphism
        """
        from sage.libs.gap.libgap import libgap
        Morphism.__init__(self, homset)
        dom = homset.domain()
        codom = homset.codomain()
        gens = [x.gap() for x in dom.gens()]
        imgs = [codom(x).gap() for x in imgs]
        if check:
            if not len(gens) == len(imgs):
                raise ValueError("provide an image for each generator")
            self._phi = libgap.GroupHomomorphismByImages(dom.gap(), codom.gap(), gens, imgs)
            # if it is not a group homomorphism, then
            # self._phi is the gap boolean fail
            if self._phi.is_bool():     # check we did not fail
                raise ValueError("images do not define a group homomorphism")
        else:
            ByImagesNC = libgap.function_factory("GroupHomomorphismByImagesNC")
            self._phi = ByImagesNC(dom.gap(), codom.gap(), gens, imgs)

    def __reduce__(self):
        r"""
        Implements pickling.

        We have to work around the fact that GAP does not provide pickling.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,5])
            sage: f = G.hom(G, G.gens())
            sage: f == loads(dumps(f))
            True
        """
        return self.parent(), (tuple(self(g) for g in self.domain().gens()),)

    def _repr_type(self):
        r"""
        Part of the implementation of :meth:`_repr_`.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2, 4])
            sage: A.hom([g^2 for g in A.gens()]) # indirect doctest
            Group endomorphism of Abelian group with gap, generator orders (2, 4)
        """
        return "Group"

    def gap(self):
        r"""
        Return the underlying LibGAP group homomorphism.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: f.gap()
            [ f1, f2 ] -> [ <identity> of ..., f3 ]
        """
        return self._phi

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: latex(f)
            \text{\texttt{Abelian group with gap, generator orders }} \left(2, 4\right)
             \rightarrow \text{\texttt{Abelian group with gap, generator orders }} \left(2, 4\right)
        """
        return r"{} \rightarrow {}".format(latex(self.domain()), latex(self.codomain()))

    def kernel(self):
        r"""
        Return the kernel of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A1 = AbelianGroupGap([6, 6])
            sage: A2 = AbelianGroupGap([3, 3])
            sage: f = A1.hom(A2.gens())
            sage: f.kernel()
            Subgroup of Abelian group with gap, generator orders (6, 6)
             generated by (f1*f2, f3*f4)
            sage: f.kernel().order()
            4
        """
        dom = self.domain()
        return dom._subgroup_constructor(self.gap().Kernel())

    def pushforward(self, J, *args, **kwds):
        r"""
        The image of an element or a subgroup.

        INPUT:

        - ``J`` -- a subgroup or an element of the domain of ``self``

        OUTPUT:

        The image of ``J`` under ``self``.

        .. NOTE::

            ``pushforward`` is the method that is used when a map is called
            on anything that is not an element of its domain. For historical
            reasons, we keep the alias ``image()`` for this method.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2])^3)
            sage: f = G.hom(H.gens())
            sage: S = G.subgroup([a.gap()])
            sage: f.pushforward(S)
            Group([ a ])
            sage: x = f.image(a)
            sage: x
            a
            sage: x.parent()
            Finitely presented group < a, b | a, b^3 >
        """
        dom = self.domain()
        codom = self.codomain()
        phi = self.gap()
        if isinstance(J, dom.Element) and J in dom:
            return self._call_(dom(J))
        if not isinstance(J, ParentLibGAP):
            raise TypeError("J (={}) must be a libgap group".format(J))
        if dom.gap().IsSubgroup(J.gap()).sage():
            return codom._subgroup_constructor(phi.Image(J.gap()))

    image = pushforward

    def _call_(self, g):
        """
        Call syntax for morphisms.

        Some python code for wrapping GAP's ``Images`` function.

        INPUT:

        - ``g`` -- an element of the domain

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2, 4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: a = A.gens()[1]
            sage: f(a)
            f3

        TESTS:

        The following tests we do fall back behind :trac:`10659`::

        sage: O = WeylGroup(['D',6])
        sage: r = prod(O.gens())
        sage: r_ = r^-1
        sage: f = O.hom([r*x*r_ for x in O.gens()])  # long time (19s on sage.math, 2011)
        sage: [f(x) for x in O.gens()]  # long time
        [
        [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
        [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
        [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
        <BLANKLINE>
        [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
        [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
        [1 0 0 0 0 0], [-1  0  0  0  0  0]
        ]
        sage: f(O)  # long time
        Matrix group over Rational Field with 6 generators
        sage: f(O).gens()   # long time
        (
        [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
        [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
        [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
        <BLANKLINE>
        [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
        [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
        [1 0 0 0 0 0], [-1  0  0  0  0  0]
        )
        """
        return self.codomain()( self.gap().Image(g.gap()) )

    def lift(self, h):
        r"""
        Return an element of the domain that maps to ``h``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: a = A.gens()[1]
            sage: f.lift(a^2)
            f2

        If the element is not in the image, we raise an error::

            sage: f.lift(a)
            Traceback (most recent call last):
            ...
            ValueError: f2 is not an element of the image of Group endomorphism
             of Abelian group with gap, generator orders (2, 4)
        """
        if not h in self.codomain():
            raise TypeError("h (={}) must be an element of the codomain".format(h))
        h = self.codomain()(h)
        phi = self.gap()
        if h.gap() not in phi.Image():
            raise ValueError("{} is not an element of the image of {}".format(h, self))
        return self.domain()( phi.PreImagesRepresentative(h.gap()) )

    def preimage(self, S):
        r"""
        Return the preimage of the subgroup ``S``.

        INPUT:

        - ``S`` -- a subgroup of this group

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: B = AbelianGroupGap([4])
            sage: f = A.hom([B.one(), B.gen(0)^2])
            sage: S = B.subgroup([B.one()])
            sage: f.preimage(S) == f.kernel()
            True
        """
        phi = self.gap()
        if not isinstance(S, ParentLibGAP):
            raise TypeError("%s must be a libGAP group of %s"%(S, self))
        if not self.codomain().gap().IsSubgroup(S.gap()).sage():
            raise ValueError("%s must be a subgroup of %s"%(S, self))
        preimage = phi.PreImage(S.gap())
        return self.domain()._subgroup_constructor(preimage)

class GroupHomset_libgap(HomsetWithBase):
    r"""
    Homsets of groups with a libgap backend.

    Do not call this directly instead use :func:`Hom`.

    INPUT:

    - ``G`` -- a libgap group
    - ``H`` -- a libgap group
    - ``category`` -- a category

    OUTPUT:

    The homset of two libgap groups.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2,4])
        sage: H = A.Hom(A)
        sage: H
        Set of Morphisms from Abelian group with gap, generator orders (2, 4)
         to Abelian group with gap, generator orders (2, 4)
         in Category of finite enumerated commutative groups
    """
    def __init__(self, G, H, category=None, check=True):
        r"""
        Return the homset of two libgap groups.

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: H = A.Hom(A)
            sage: TestSuite(H).run()
        """
        if check:
            if not isinstance(G, ParentLibGAP):
                raise TypeError("G (={}) must be a ParentLibGAP group".format(G))
            if not isinstance(H, ParentLibGAP):
                raise TypeError("H (={}) must be a ParentLibGAP group".format(H))
        HomsetWithBase.__init__(self, G, H, category, check=check, base=ZZ)

    Element = GroupMorphism_libgap

    def _element_constructor_(self, x, check=True, **options):
        r"""
        Handle conversions and coercions.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: H = A.Hom(A)
            sage: H([g^2 for g in A.gens()])
            Group endomorphism of Abelian group with gap, generator orders (2, 4)
            sage: H(2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2 to an element of ...
        """
        if isinstance(x, (tuple, list)):
            codomain = self.codomain()
            im_gens = tuple([codomain(g) for g in x])
            return self.element_class(self, im_gens, check=check, **options)
        return super(GroupHomset_libgap, self)._element_constructor_(x, check=check, **options)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        There is only one group homomorphism which is always defined:
        the trivial one.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: G.<a,b> = FreeGroup()
            sage: A.Hom(G).an_element()
            Group morphism:
              From: Abelian group with gap, generator orders (2, 4)
              To:   Free Group on generators {a, b}
        """
        n = len(self.domain().gens())
        one = self.codomain().one()
        return self.element_class(self, (one,) * n, check=False)

