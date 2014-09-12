r"""
Cartesian product of projective spaces

This class builds on the projective space classe and its point and morphism classes.

EXAMPLES: We construct cartesian products projective spaces of various dimensions over the same ring.::

    sage: P1 = ProjectiveSpace(ZZ,1,'x')
    sage: P2 = ProjectiveSpace(ZZ,2,'y')
    sage: ProjectiveSpace_cartesian_product([P1,P2])
    Product of Projective Space of dimension 1 over Integer Ring and
    Projective Space of dimension 2 over Integer Ring

We can also construct the product by specifying the dimensions and the base ring::

    sage: ProjectiveSpace_cartesian_product([1,2,3],QQ,'z')
    Product of Projective Space of dimension 1 over Rational Field and
    Projective Space of dimension 2 over Rational Field and Projective Space
    of dimension 3 over Rational Field

AUTHORS:

- Ben Hutz (2014)

"""

from copy import copy
from sage.misc.latex import latex
from sage.rings.all import (PolynomialRing, ZZ)

from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective_cartesian_product
from sage.schemes.generic.ambient_space import AmbientSpace
from sage.schemes.generic.morphism import SchemeMorphism
from sage.schemes.generic.morphism import SchemeMorphism_point
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.schemes.projective.projective_space import ProjectiveSpace, ProjectiveSpace_ring
from sage.schemes.projective.projective_homset import SchemeHomset_points_projective_cartesian_product


def is_ProjectiveSpace_cartesian_product(x):
    r"""
    Return True if `x` is a product of projective spaces, i.e., an ambient space
    `\mathbb{P}^n_R \times \cdots \times \mathbb{P}^m_R`, where `R` is a ring and `n,\ldots,m\geq 0`
    are integers. 

    EXAMPLES::

        sage: is_ProjectiveSpace_cartesian_product(ProjectiveSpace(5, names='x'))
        False
        sage: is_ProjectiveSpace_cartesian_product(ProjectiveSpace_cartesian_product([1,2,3], ZZ, 'x'))
        True
    """
    return isinstance(x, ProjectiveSpace_cartesian_product_generic)

def ProjectiveSpace_cartesian_product(n, R=None, names='x'):
    r"""
    Returns the cartesian product of projective spaces. Can input either a list of projective spaces
    over the same base ring or the list of dimensions, the base ring, and the variable names.

    INPUT:

    - ``n`` -- a list of integers or a list of projective spaces

    - ``R`` -- a ring

    - ``names`` -- a string

    EXAMPLES::

        sage: P1 = ProjectiveSpace(QQ,2,'x')
        sage: P2 = ProjectiveSpace(QQ,3,'y')
        sage: ProjectiveSpace_cartesian_product([P1,P2])
        Product of Projective Space of dimension 2 over Rational Field and
        Projective Space of dimension 3 over Rational Field

    ::

        sage: ProjectiveSpace_cartesian_product([2,2],GF(7),'y')
        Product of Projective Space of dimension 2 over Finite Field of size 7
        and Projective Space of dimension 2 over Finite Field of size 7

    ::

        sage: P1 = ProjectiveSpace(ZZ,2,'x')
        sage: P2 = ProjectiveSpace(QQ,3,'y')
        sage: ProjectiveSpace_cartesian_product([P1,P2])
        Traceback (most recent call last):
        ...
        AttributeError: Components must be over the same base ring
    """
    if isinstance(R, (list, tuple)):
        n, R = R, n
    if R is None:
        R = ZZ  # default is the integers
    if isinstance(n[0],ProjectiveSpace_ring):
        #this should be a list of projective spaces
        names = []
        N = []
        R = None
        for PS in n:
            if isinstance(PS,ProjectiveSpace_ring)==False:
                raise TypeError("Must be a list of Projective Spaces or (dimensions,base ring,names)")
            if R == None:
                R = PS.base_ring()
            elif R != PS.base_ring():
                raise AttributeError("Components must be over the same base ring")
            N.append(PS.dimension_relative())
            names += PS.variable_names()
        X = ProjectiveSpace_cartesian_product_generic(N, R, names)
        X._components = n
    else:
        if isinstance(R, (list,tuple)):
           n, R = R, n
        if not isinstance(n,(list,tuple)):
            raise ValueError("Need list or tuple of dimensions")
        from sage.categories.commutative_rings import CommutativeRings
        if R not in CommutativeRings():
            raise ValueError("Must be a commutative ring")
        from sage.structure.parent_gens import normalize_names
        names = normalize_names(sum(n) + len(n), names)
        X = ProjectiveSpace_cartesian_product_generic(n, R, names)
    return(X)

class ProjectiveSpace_cartesian_product_generic(AmbientSpace):
    r"""
    Cartesian product of projective spaces `\mathbb{P}^{n_1} \times \cdots \times \mathbb{P}^{n_r}`.

    EXAMPLES::

        sage: P.<x0,x1,x2,x3,x4> = ProjectiveSpace_cartesian_product([1,2],QQ); P
        Product of Projective Space of dimension 1 over Rational Field and
        Projective Space of dimension 2 over Rational Field
        sage: P.coordinate_ring()
        Multivariate Polynomial Ring in x0, x1, x2, x3, x4 over Rational Field
        sage: P[0]
        Projective Space of dimension 1 over Rational Field
        sage: P[1]
        Projective Space of dimension 2 over Rational Field
        sage: Q = P(6,3,2,2,2); Q
        (2 : 1 , 1 : 1 : 1)
        sage: Q[0]
        (2 : 1)
        sage: H = Hom(P,P)
        sage: f = H([x0^2*x3,x2*x1^2,x2^2,2*x3^2,x4^2])
        sage: f(Q)
        (4 : 1 , 1 : 2 : 1)
    """
    def __init__(self, N, R = ZZ, names = None):
        r"""
        The Python constructor

        INPUT:

        - ``N`` - a list or tuple of positive integers

        - ``R`` - a ring

        - ``names`` - a tuple or list of strings. This must either be a single variable name
                    or the complete list of variables.

        Examples:: 

            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([2,2],QQ)
            sage: T
            Product of Projective Space of dimension 2 over Rational Field and
            Projective Space of dimension 2 over Rational Field 
            sage: T.coordinate_ring()
            Multivariate Polynomial Ring in x, y, z, u, v, w over Rational Field
            sage: T[1].coordinate_ring()
            Multivariate Polynomial Ring in u, v, w over Rational Field

        ::

            sage: ProjectiveSpace_cartesian_product([1,1,1],ZZ, ['x','y','z','u','v','w'])
            Product of Projective Space of dimension 1 over Integer Ring and
            Projective Space of dimension 1 over Integer Ring and Projective Space
            of dimension 1 over Integer Ring

        ::

            sage: T = ProjectiveSpace_cartesian_product([1,1],QQ,'z')
            sage: T.coordinate_ring()
            Multivariate Polynomial Ring in z0, z1, z2, z3 over Rational Field
        """
        AmbientSpace.__init__(self, sum(N), R)
        start = 0
        self._components = []
        for i in range(len(N)):
            self._components.append(ProjectiveSpace(N[i],R,names[start:start+N[i]+1]))
            start += N[i]+1
        #Note that the coordinate ring should really be the tensor product of the component
        #coordinate rings. But we just deal with them as multihomogeneous polynomial rings
        self._coordinate_ring = PolynomialRing(R,sum(N)+ len(N),names)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: String.

        EXAMPLES::

            sage: ProjectiveSpace_cartesian_product([1,1,1],ZZ, ['x','y','z','u','v','w'])
            Product of Projective Space of dimension 1 over Integer Ring and
            Projective Space of dimension 1 over Integer Ring and Projective Space
            of dimension 1 over Integer Ring
        """
        S='Product of '
        for i in range(self.ncomponents()-1):
            S += self[i]._repr_() +' and '
        S += self[self.ncomponents()-1]._repr_()
        return(S)

    def _repr_generic_point(self, v = None):
        """
        Return a string representation of the generic point
        on this product space.

        EXAMPLES::
            sage: T = ProjectiveSpace_cartesian_product([1,2,1], QQ, 'x')
            sage: T._repr_generic_point()
            '(x0 : x1 , x2 : x3 : x4 , x5 : x6)'
        """
        if v is None:
            v = list(self.gens())
        else:
            v = list(v)
        start = 0
        N=self.dimension_relative_components()
        splitv = []
        for i in range(len(N)):
            splitv.append(v[start : start+N[i]+1])
            start += N[i]+1
        return '(%s)'%(" , ".join((" : ".join([str(t) for t in P])) for P in splitv))

    def _latex_(self):
        r"""
        Return a LaTeX representation of this product space.

        EXAMPLES::

            sage: latex(ProjectiveSpace_cartesian_product([1,2,3], ZZ, 'x'))
            {\mathbf P}_{\Bold{Z}}^1 \times {\mathbf P}_{\Bold{Z}}^2 \times {\mathbf
            P}_{\Bold{Z}}^3
        """
        return '%s'%" \\times ".join([latex(PS) for PS in self])

    def _latex_generic_point(self, v = None):
        """
        Return a LaTeX representation of the generic point
        on this product space.

        EXAMPLES::
            sage: T = ProjectiveSpace_cartesian_product([1,2,3], ZZ, 'x')
            sage: T._latex_generic_point()
            '\\left(x_{0} : x_{1} , x_{2} : x_{3} : x_{4} , x_{5} : x_{6} : x_{7} : x_{8}\\right)'
        """
        if v is None:
            v = list(self.gens())
        else:
            v = list(v)
        start = 0
        N=self.dimension_relative_components()
        splitv = []
        for i in range(len(N)):
            splitv.append(v[start : start+N[i]+1])
            start += N[i]+1
        return '\\left(%s\\right)'%(" , ".join((" : ".join([latex(t) for t in P])) for P in splitv))

    def __getitem__(self, i):
        r"""
        Return the `i`-th component of the product.

        INPUT:

        - ``i`` - a positive integer

        OUTPUT: A Projective Space.

        Examples:: 

            sage: T.<a,x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([3,2],QQ)
            sage: T[0]
            Projective Space of dimension 3 over Rational Field 
        """
        return(self._components[i])

    def __cmp__(self, right):
        r"""
        Compare two products of projective spaces.

        INPUT:

        - ``right`` - a product of projective spaces

        OUTPUT: Boolean

        Examples:: 

            sage: S.<a,x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([3,2],QQ)
            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([2,2],QQ)
            sage: S == T
            False

        ::

            sage: S.<a,x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([3,2],QQ)
            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([2,2],QQ)
            sage: S != T
            True
        """
        if not isinstance(right, (ProjectiveSpace_cartesian_product_generic)):
            return -1
        else:
            return(cmp(self._components,right._components))

    def dimension_relative(self):
        r"""
        Return the relative dimension of the product of projective spaces.

        OUTPUT: a positive integer.

        Examples:: 

            sage: T.<a,x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([3,2],QQ)
            sage: T.dimension_relative()
            5
        """ 
        return(sum([self[i].dimension_relative() for i in range(self.ncomponents())]))

    def dimension_absolute(self):
        r"""
        Return the absolute dimension of the product of projective spaces.

        OUTPUT: a positive integer.

        Examples::

            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([2,2],GF(17))
            sage: T.dimension_absolute()
            4
            sage: T.dimension()
            4
        """
        base = self.base_scheme()
        if base.is_noetherian():
            return sum([self[i].dimension_relative() + base.dimension() for i in range(self.ncomponents())])
        raise NotImplementedError("Cannot compute the dimension of this scheme.")

    dimension = dimension_absolute

    def dimension_relative_components(self):
        r"""
        Return the relative dimension of the product of projective spaces.

        OUTPUT: a list of positive integers.

        Examples::

            sage: T.<a,x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([3,2],QQ)
            sage: T.dimension_relative_components()
            [3, 2]
        """
        return([self[i].dimension_relative() for i in range(self.ncomponents())])

    def dimension_absolute_components(self):
        r"""
        Return the absolute dimension of the product of projective spaces.

        OUTPUT: a list of positive integers.

        Examples::

            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([2,2],GF(17))
            sage: T.dimension_absolute_components()
            [2, 2]
            sage: T.dimension_components()
            [2, 2]
        """
        base = self.base_scheme()
        if base.is_noetherian():
            return [self[i].dimension_relative() + base.dimension() for i in range(self.ncomponents())]
        raise NotImplementedError("Cannot compute the dimension of this scheme.")

    dimension_components = dimension_absolute_components

    def ncomponents(self):
        r"""
           Returns the number of components of ``self``.

        OUTPUT: an integer.

        EXAMPLES::

            sage: T = ProjectiveSpace_cartesian_product([1,1,1],GF(5),'x')
            sage: T.ncomponents()
            3
        """
        return(len(self._components))

    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P = ProjectiveSpace_cartesian_product([1,1],QQ, 'z')
            sage: point_homset = P._point_homset(Spec(QQ), P)
            sage: P._point(point_homset, [2,2,1,1])
            (1 : 1 , 1 : 1)
        """
        return ProjectiveSpace_cartesian_product_point(*args, **kwds)

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace_cartesian_product([1,1],QQ)
            sage: P._morphism(P.Hom(P), [x,y,z,w])
            Scheme endomorphism of Product of Projective Space of dimension 1 over
            Rational Field and Projective Space of dimension 1 over Rational Field
              Defn: Defined by sending (x : y , z : w) to
                  (x : y , z : w).
        """
        return ProjectiveSpace_cartesian_product_morphism(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace_cartesian_product([1,1],GF(5))
            sage: P._point_homset(Spec(GF(5)), P)
            Set of rational points of Product of Projective Space of dimension 1
            over Finite Field of size 5 and Projective Space of dimension 1 over
            Finite Field of size 5
        """
        return SchemeHomset_points_projective_cartesian_product(*args, **kwds)

    def _validate(self, polynomials):
        r"""
        If ``polynomials`` is a tuple of valid polynomial functions on self,
        return ``polynomials``, otherwise raise a TypeError.

        Since this is a product of projective spaces, the polynomials must be multi-homogeneous.

        INPUT:

        - ``polynomials`` -- tuple of polynomials in the coordinate ring of ``self``

        OUTPUT:

        - tuple of polynomials in the coordinate ring of ``self``

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: T._validate([x^2*u,y^2*w,z^2*u,w^2,u^2])
            [x^2*u, y^2*w, z^2*u, w^2, u^2]

        ::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: T._validate([x^2+w^2,y^2*w,z^2*u,w^2,u^2])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[x^2 + w^2, y^2*w, z^2*u, w^2, u^2]) must be multi-homogeneous

        ::

            sage: R.<t> = PolynomialRing(GF(5))
            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: T._validate([t,t,t,w^2,u^2])
            Traceback (most recent call last):
            ...
            TypeError: polynomials (=[t, t, t, w^2, u^2]) must be elements of Multivariate
            Polynomial Ring in x, y, z, w, u over Rational Field
        """
        if not isinstance(polynomials, (list, tuple)):
            raise TypeError('The argument polynomials=%s must be a list or tuple'%polynomials)
        #check in the coordinate ring
        source_ring = self.coordinate_ring()
        try:
            polynomials = [source_ring(poly) for poly in polynomials]
        except TypeError:
            raise TypeError("polynomials (=%s) must be elements of %s"%(polynomials,source_ring))
        #check multi-homogeneous by looking at the exponents of each monomial
        dims = self.dimension_relative_components()
        for f in polynomials:
            E = f.exponents()
            index = 0
            for i in range(len(dims)):
                d = sum(E[0][index+j] for j in range(0, dims[i]+1))
                if all([d == sum(E[k][index+j] for j in range(0,dims[i]+1)) for k in range(len(E))]) == False:
                        raise TypeError("polys (=%s) must be multi-homogeneous"%polynomials)
                index += dims[i]+1
        return polynomials

    def _check_satisfies_equations(self,v):
        """
        Return True if ``v`` defines a point on the scheme ``self``; raise a
        TypeError otherwise.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: T._check_satisfies_equations([0,1,1,1,1])
            True

        ::

            sage: R.<t> = PolynomialRing(GF(7))
            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],R)
            sage: T._check_satisfies_equations([1+t,1,0,0,1])
            True

        ::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],ZZ)
            sage: T._check_satisfies_equations([1,1,1,0,0])
            Traceback (most recent call last):
            ...
            TypeError: The zero vector is not a point in projective space

        ::

            sage: T.<x,y,z,w> = ProjectiveSpace_cartesian_product([1,1],ZZ)
            sage: T._check_satisfies_equations([1,1,1,0,0])
            Traceback (most recent call last):
            ...
            TypeError: The list v=[1, 1, 1, 0, 0] must have 4 components

        ::

            sage: T.<x,y,z,w> = ProjectiveSpace_cartesian_product([1,1],ZZ)
            sage: T._check_satisfies_equations([1,1/2,1,0])
            Traceback (most recent call last):
            ...
            TypeError: The components of v=[1, 1/2, 1, 0] must be elements of Integer Ring
        """
        if not isinstance(v, (list, tuple)):
            raise TypeError('The argument v=%s must be a list or tuple'%v)
        n = self.ngens()
        if not len(v) == n:
            raise TypeError('The list v=%s must have %s components'%(v, n))
        R = self.base_ring()
        try:
            n = [R(w) for w in v]
        except TypeError:
            raise TypeError('The components of v=%s must be elements of %s'%(v, R))
        #check if any of the component points are 0
        N = self.dimension_relative_components()
        start = 0
        for i in range(len(N)):
            if v[start:start + N[i]+1] == [R(0)]*(N[i]+1):
                raise TypeError('The zero vector is not a point in projective space')
            start += N[i]+1
        return True

    def subscheme(self, X):
        r"""
        Return the closed subscheme defined by ``X``.

        INPUT:

        - ``X`` - a list or tuple of equations

        OUTPUT:

        - :class:`AlgebraicScheme_subscheme_projective_cartesian_product`

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace_cartesian_product([1,1],GF(5))
            sage: X = P.subscheme([x-y,z-w]);X
            Closed subscheme of Product of Projective Space of dimension 1 over
            Finite Field of size 5 and Projective Space of dimension 1 over Finite
            Field of size 5 defined by:
              x - y,
              z - w
            sage: X.defining_polynomials ()
            [x - y, z - w]
            sage: I = X.defining_ideal(); I
            Ideal (x - y, z - w) of Multivariate Polynomial Ring in x, y, z, w over
            Finite Field of size 5
            sage: X.dimension()
            0
            sage: X.base_ring()
            Finite Field of size 5
            sage: X.base_scheme()
            Spectrum of Finite Field of size 5
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of Product of Projective Space of dimension 1 over
              Finite Field of size 5 and Projective Space of dimension 1 over Finite
              Field of size 5 defined by:
              x - y,
              z - w
              To:   Spectrum of Finite Field of size 5
              Defn: Structure map
        """
        return AlgebraicScheme_subscheme_projective_cartesian_product(self, X)

    def change_ring(self, R):
        r"""
        Return a product of projective spaces over a ring `R` and otherwise the same as ``self``.

        INPUT:

        - ``R`` -- commutative ring

        OUTPUT:

        - product of projective spaces over ``R``

        .. NOTE::

            There is no need to have any relation between `R` and the base ring
            of ``self``, if you want to have such a relation, use
            ``self.base_extend(R)`` instead.

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([2,2],QQ)
            sage: T.change_ring(GF(17))
            Product of Projective Space of dimension 2 over Finite Field of size 17
            and Projective Space of dimension 2 over Finite Field of size 17
        """
        newComponents = [P.change_ring(R) for P in self._components]
        return ProjectiveSpace_cartesian_product(newComponents)

    def affine_patch(self, I, return_embedding = False):
        r"""
        Return the `I^{th}` affine patch of this projective space product
        where ``I`` is a multi-index.

        INPUT:

        - ``I`` -- a list or tuple of positive integers

        - ``return_embedding`` -- Boolean, if true the projective embedding is also returned

        OUTPUT:

        - An affine space

        - An embedding into a product of projective spaces (optional)

        EXAMPLES::

            sage: PP = ProjectiveSpace_cartesian_product([2,2,2],ZZ,'x')
            sage: A, phi = PP.affine_patch([0,1,2],True);A
            Affine Space of dimension 6 over Integer Ring
            sage: phi
            Scheme morphism:
              From: Affine Space of dimension 6 over Integer Ring
              To:   Product of Projective Space of dimension 2 over Integer Ring and
            Projective Space of dimension 2 over Integer Ring and Projective Space
            of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x0, x1, x2, x3, x4, x5) to
                    (1 : x0 : x1 , x2 : 1 : x3 , x4 : x5 : 1)
        """
        if not isinstance(I, (list, tuple)):
            raise TypeError('The argument I=%s must be a list or tuple of positive integers'%I)
        PP = self.ambient_space()
        N = PP.dimension_relative_components()
        if len(I) != len(N):
            raise ValueError('The argument I=%s must have %s entries'%(I,len(N)))
        I = tuple([int(i) for i in I])   # implicit type checking
        for i in range(len(I)):
            if I[i] < 0 or I[i] > N[i]:
                raise ValueError("Argument i (= %s) must be between 0 and %s."%(I[i], N[i]))
        try:
            if return_embedding:
                return self.__affine_patches[I]
            else:
                return self.__affine_patches[I][0]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        from sage.schemes.affine.affine_space import AffineSpace
        AA = AffineSpace(PP.base_ring(),sum(N),'x')
        v = list(AA.gens())
        index = 0
        for i in range(len(I)):
            v.insert(index+I[i],1)
            index += N[i]+1
        phi = AA.hom(v,self)
        self.__affine_patches.update({I:(AA,phi)})
        if return_embedding:
            return AA,phi
        else:
            return AA


    def segre_embedding(self, PP = None):
        r"""
        Returns the Segre embedding of ``self`` into the appropriate projective space.

        INPUT:

        -  ``PP`` - (default: None) ambient image projective space;
            this is constructed if it is not given.

        OUTPUT:

        - Hom -- from self to the appropriate subscheme of projective space

        .. TODO::

            Cartesian products with more than two components

        EXAMPLES::

            sage: X.<y0,y1,y2,y3,y4,y5> = ProjectiveSpace_cartesian_product(ZZ,[2,2])
            sage: phi = X.segre_embedding(); phi
             Scheme morphism:
                  From: Product of Projective Space of dimension 2 over Integer Ring and Projective Space of dimension 2 over Integer Ring
                  To:   Closed subscheme of Projective Space of dimension 8 over Integer Ring defined by:
                  -u5*u7 + u4*u8,
                  -u5*u6 + u3*u8,
                  -u4*u6 + u3*u7,
                  -u2*u7 + u1*u8,
                  -u2*u4 + u1*u5,
                  -u2*u6 + u0*u8,
                  -u1*u6 + u0*u7,
                  -u2*u3 + u0*u5,
                  -u1*u3 + u0*u4
                  Defn: Defined by sending (y0 : y1 : y2 , y3 : y4 : y5) to 
                        (y0*y3 : y0*y4 : y0*y5 : y1*y3 : y1*y4 : y1*y5 : y2*y3 : y2*y4 : y2*y5).

            ::

            sage: T = ProjectiveSpace_cartesian_product([1,2],CC,'z')
            sage: T.segre_embedding()
            Scheme morphism:
                  From: Product of Projective Space of dimension 1 over Complex Field with 53 bits of precision and Projective Space of dimension 2 over Complex Field with 53 bits of precision
                  To:   Closed subscheme of Projective Space of dimension 5 over Complex Field with 53 bits of precision defined by:
                  -u2*u4 + u1*u5,
                  -u2*u3 + u0*u5,
                  -u1*u3 + u0*u4
                  Defn: Defined by sending (z0 : z1 , z2 : z3 : z4) to 
                        (z0*z2 : z0*z3 : z0*z4 : z1*z2 : z1*z3 : z1*z4).
        """
        N = self.dimension_relative_components()
        if len(N) > 2:
            raise NotImplementedError("Cannot have more than two components.")
        M = (N[0]+1)*(N[1]+1)-1

        vars = list(self.coordinate_ring().variable_names()) + ['u' + str(i) for i in range(M+1)]
        from sage.rings.all import PolynomialRing
        R = PolynomialRing(self.base_ring(),self.ngens()+M+1,vars,order='lex')

        #set-up the elimination for the segre embedding
        mapping = []
        k = self.ngens()
        for i in range(N[0]+1):
            for j in range(N[0]+1,N[0]+N[1]+2):
                mapping.append(R.gen(k)-R(self.gen(i)*self.gen(j)))
                k+=1

        #change the defining ideal of the subscheme into the variables
        I = R.ideal(list(self.defining_polynomials()) + mapping)
        J=I.groebner_basis()
        s=set(R.gens()[:self.ngens()])
        n=len(J)-1
        L=[]
        while s.isdisjoint(J[n].variables()):
            L.append(J[n])
            n=n-1

        #create new subscheme
        if PP is None:
            from sage.schemes.projective.projective_space import ProjectiveSpace
            PS = ProjectiveSpace(self.base_ring(),M,R.gens()[self.ngens():])
            Y = PS.subscheme(L)
        else:
            if PP.dimension_relative()!= M:
                raise ValueError("Projective Space %s must be dimension %s")%(PP, M)
            S = PP.coordinate_ring()
            psi = R.hom([0]*(N[0]+N[1]+2) + list(S.gens()),S)
            L = [psi(l) for l in L]
            Y = PP.subscheme(L)

        #create embedding for points
        mapping = []
        for i in range(N[0]+1):
            for j in range(N[0]+1,N[0]+N[1]+2):
                mapping.append(self.gen(i)*self.gen(j))
        phi = self.hom(mapping,Y)

        return phi


class ProjectiveSpace_cartesian_product_point(SchemeMorphism_point):
    r"""
    The class of points on products of projective spaces.
    The components are projective space points.

    Examples::

        sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
        sage: T.point([1,2,3,4,5]);
        (1/3 : 2/3 : 1 , 4/5 : 1)
    """
    def __init__(self, parent, polys, check = True):
        r"""
        The Python constructor.

        INPUT:

            - ``parent`` - Homset

            - ``polys`` - anything that defines a point in the class

            - ``check`` - Boolean. Whether or not to perform input checks (Default: True)

        EXAMPLES::

            sage: P1.<x0,x1,x2> = ProjectiveSpace(QQ,2)
            sage: P2 = ProjectiveSpace(QQ,3,'y')
            sage: T = ProjectiveSpace_cartesian_product([P1,P2])
            sage: Q1 = P1(1,1,1)
            sage: Q2 = P2(1,2,3,4)
            sage: T([Q1,Q2])
            (1 : 1 : 1 , 1/4 : 1/2 : 3/4 : 1)

        ::

            sage: T = ProjectiveSpace_cartesian_product([2,2,2],GF(5),'x')
            sage: T.point([1,2,3,4,5,6,7,8,9])
            (2 : 4 : 1 , 4 : 0 : 1 , 3 : 2 : 1)

        ::

            sage: T.<x,y,z,w> = ProjectiveSpace_cartesian_product([1,1],GF(5))
            sage: X = T.subscheme([x-y,z-2*w])
            sage: X([1,1,2,1])
            (1 : 1 , 2 : 1)
        """
        SchemeMorphism.__init__(self, parent)
        if all(isinstance(P, SchemeMorphism_point) for P in polys):
            if check == True:
                Q = []
                self._points=[]
                for i in range(len(polys)):
                    if polys[i].codomain() != parent.codomain().ambient_space()[i]:
                        raise ValueError("Points must be in correct projective spaces")
                    Q += list(polys[i])
                    self._points.append(polys[i])
                parent.codomain()._check_satisfies_equations(Q)
            self._points = polys
        else:
            N=parent.codomain().ambient_space().dimension_relative_components()
            if check == True:
                parent.codomain()._check_satisfies_equations(polys)
                if sum(N) + len(N) != len(polys):
                   raise TypeError("v (=%s) must have %s components"%(v, sum(N) + len(N)))
            self._points = []
            total = 0
            for i in range(len(N)):
                self._points.append(parent.codomain().ambient_space()[i].point([polys[j] for j in range(total,total+N[i]+1)], check))
                total += N[i]+1

    def __getitem__(self, i):
        r"""
        Return the ``n``-th coordinate point.

        INPUT:

            - ``i`` - integer

        OUTPUT: The projective space point that is the ``n``-th coordinate.

        EXAMPLES::

            sage: T = ProjectiveSpace_cartesian_product([2,2,2],GF(5),'x')
            sage: P = T([1,0,1,1,0,0,0,0,1])
            sage: P[1]
            (1 : 0 : 0)
            sage: P[1].codomain()
            Projective Space of dimension 2 over Finite Field of size 5
            sage: P[1][0]
            1
        """
        return(self._points[i])

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: String.

        EXAMPLES::

            sage: T = ProjectiveSpace_cartesian_product([2,2],ZZ,'x')
            sage: P = T([1,2,3,4,5,6])
            sage: P._repr_()
            '(1 : 2 : 3 , 4 : 5 : 6)'
        """
        return('(%s)'%(" , ".join((" : ".join([repr(f) for f in Q])) for Q in self._points)))

    def __eq__(self, other):
        r"""
        Tests the projective equality of two points.

        INPUT:

        - ``right`` - a point on a product of projective spaces

        OUTPUT:

        - Boolean - True if ``self`` and ``right`` define the same point. False otherwise.

        Examples::

            sage: T = ProjectiveSpace_cartesian_product([2,2],ZZ,'x')
            sage: P = T([1,2,3,4,5,6])
            sage: Q = T([2,4,6,4,5,6])
            sage: P == Q
            True
        """
        for i in range(self.codomain().ncomponents()):
           if self[i] != other[i]:
               return False
        return True

    def __ne__(self, other):
        r"""
        Tests the projective inequality of two points.

        INPUT:

        - ``right`` - a point on a product of projective spaces

        OUTPUT:

        - Boolean - False if ``self`` and ``right`` define the same point. True otherwise.

        Examples::

            sage: T = ProjectiveSpace_cartesian_product([1,1,1],ZZ,'x')
            sage: P = T([1,2,3,4,5,6])
            sage: Q = T([2,4,6,4,1,0])
            sage: P != Q
            True
        """
        for i in range(self.codomain().ncomponents()):
           if self[i] != other[i]:
               return True
        return False

    def __cmp__(self, right):
        r"""
        Compare two points in products of projective spaces.

        INPUT:

        - ``other`` -- anything. To compare against the point ``self``.

        OUTPUT: ``+1``, ``0``, or ``-1``.

        EXAMPLES::

            sage: T = ProjectiveSpace_cartesian_product([1,1,1],GF(5),'x')
            sage: P = T([3,2,3,4,1,0])
            sage: Q = T([1,2,3,4,3,1])
            sage: P.__cmp__(Q)
            1

        ::

            sage: T = ProjectiveSpace_cartesian_product([1,1,1],GF(5),'x')
            sage: P = T([1,2,3,4,1,0])
            sage: Q = T([1,2,3,4,3,0])
            sage: P.__cmp__(Q)
            0

        ::

            sage: T = ProjectiveSpace_cartesian_product([1,1,1],GF(5),'x')
            sage: P = T([1,2,3,4,1,0])
            sage: Q = T([1,2,3,4,3,1])
            sage: P.__cmp__(Q)
            -1
        """
        #needed for Digraph
        if not isinstance(right, (ProjectiveSpace_cartesian_product_point)):
            return -1
        else:
            return(cmp(self._points, right._points))

    def __copy__(self):
        r"""
        Returns a copy of the point ``self``.

        OUTPUT:

        - a point in the same space as ``self``.

        EXAMPLES::

            sage: T = ProjectiveSpace_cartesian_product([1,1],QQ,'x')
            sage: P = T([2,1,0,1])
            sage: Q = P.__copy__()
            sage: P is Q
            False
            sage: P == Q
            True
        """
        P = [copy(self[i]) for i in range(self.codomain().ambient_space().ncomponents())]
        return(self.codomain().point(P, False))

    def __iter__(self):
        r"""
        Iterate over the coordinates of the point.

        OUTPUT: An iterator.

        EXAMPLES::

            sage: T = ProjectiveSpace_cartesian_product([1,1],QQ,'x')
            sage: P = T([2,1,0,1])
            sage: iter = P.__iter__()
            sage: iter.next()
            2
            sage: iter.next()
            1
            sage: list(P)
            [2, 1, 0, 1]
        """
        L = []
        for P in self._points:
            L += P._coords
        return iter(L)

    def normalize_coordinates(self):
        r"""
        Removes common factors (componentwise) from the coordinates of ``self`` (including `-1`).

        OUTPUT: None.

        Examples::

            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([2,2],ZZ)
            sage: P = T.point([5,10,15,4,2,6]);
            sage: P.normalize_coordinates()
            sage: P
            (1 : 2 : 3 , 2 : 1 : 3)
        """
        for i in range(self.codomain().ambient_space().ncomponents()):
            self[i].normalize_coordinates()

    def scale_by(self, t):
        r"""
        Scale the coordinates of the point ``self`` by `t`, done componentwise.
        A ``TypeError`` occurs if the point is not in the base ring of the codomain after scaling.

        INPUT:

        - ``t`` -- a ring element

        Examples::

            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([1,1,1],ZZ)
            sage: P = T.point([5,10,15,4,2,6]);
            sage: P.scale_by([2,1,1])
            sage: P
            (10 : 20 , 15 : 4 , 2 : 6)
        """
        if isinstance(t,(tuple,list)) == False:
            raise TypeError("%s must be a list or tuple"%t)
        if len(t) != self.codomain().ambient_space().ncomponents():
            raise TypeError("%s must have same number of components as %r"%(t,self))
        for i in range(self.codomain().ambient_space().ncomponents()):
            self[i].scale_by(t[i])

    def change_ring(self, R, check = True):
        r"""
        Returns a new :class:`ProjectiveSpace_cartesian_product_point` which is ``self`` coerced to ``R``.
        If ``check`` is True, then the initialization checks are performed.

        INPUT:

        - ``R`` -- a ring

        - ``check`` -- Boolean (optional)

        OUTPUT:

        - :class:`ProjectiveSpace_cartesian_product_point`

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProjectiveSpace_cartesian_product([1,1,1],ZZ)
            sage: P = T.point([5,3,15,4,2,6]);
            sage: P.change_ring(GF(3))
            (1 : 0 , 0 : 1 , 1 : 0)
        """
        S = self.codomain().change_ring(R)
        Q = [P.change_ring(R,check) for P in self._points]
        return(S.point(Q, check))

class ProjectiveSpace_cartesian_product_morphism(SchemeMorphism_polynomial):
    r"""
    The class of morphisms on products of projective spaces.
    The components are projective space morphisms.

    Examples::

        sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
        sage: H = T.Hom(T)
        sage: H([x^2,y^2,z^2,w^2,u^2])
        Scheme endomorphism of Product of Projective Space of dimension 2 over
        Rational Field and Projective Space of dimension 1 over Rational Field
          Defn: Defined by sending (x : y : z , w : u) to 
                (x^2 : y^2 : z^2 , w^2 : u^2).
    """

    def __init__(self, parent, polys, check = True):
        r"""
        The Python constructor.

        INPUT:

        - ``parent`` - Homset

        - ``polys`` - anything that defines a point in the class

        - ``check`` - Boolean. Whether or not to perform input checks (Default: True)

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u,y^2*w,z^2*u,w^2,u^2])
            Scheme endomorphism of Product of Projective Space of dimension 2 over
            Rational Field and Projective Space of dimension 1 over Rational Field
              Defn: Defined by sending (x : y : z , w : u) to 
                    (x^2*u : y^2*w : z^2*u , w^2 : u^2).

        ::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u,y^2*w,z^2*u,w^2,u*z])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[x^2*u, y^2*w, z^2*u, w^2, z*u]) must be
            multi-homogeneous of the same degrees (by component)
        """
        if check == True:
            #check multi-homogeneous
            target = parent.codomain().ambient_space()
            if is_ProjectiveSpace_cartesian_product(target):
                dims = target.dimension_relative_components()
            else:
                dims = [target.dimension_relative()]
            index = 0
            #if self is a subscheme, we may need the lift of the polynomials
            try:
                polys[0].exponents()
            except AttributeError:
                polys = [f.lift() for f in polys]
            splitpolys=[]
            for i in range(len(dims)):
                splitpolys.append([polys[index+j] for j in range(0,dims[i]+1)])
                index+=dims[i]+1
            for m in range(len(splitpolys)):
                E = splitpolys[m][0].exponents()
                index = 0
                d = []
                for i in range(len(dims)):
                    d.append(sum(E[0][index+j] for j in range(0,dims[i]+1)))
                    index += dims[i]+1
                for f in splitpolys[m]:
                    E = f.exponents()
                    index = 0
                    for i in range(len(dims)):
                        if all([d[i] == sum(E[k][index+j] for j in range(0,dims[i]+1)) for k in range(len(E))]) == False:
                            raise  TypeError("polys (=%s) must be multi-homogeneous of the same degrees (by component)"%polys)
                        index += dims[i]+1
        SchemeMorphism_polynomial.__init__(self, parent, polys, check)

    def __getitem__(self, i):
        r"""
        Return the ``n``-th coordinate polynomial.

        INPUT:

            - ``i`` - integer

        OUTPUT:

        The (multi)-homomgeneous polynomial that is the ``n``-th coordinate.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u,y^2*w,z^2*u,w^2,u^2])
            sage: F[2]
            z^2*u
        """
        return(self._polys[i])

    def _repr_defn(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: String.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace_cartesian_product([1,1], QQ)
            sage: H = Hom(P,P)
            sage: f = H([x^2,y^2,z,w])
            sage: f._repr_defn()
            'Defined by sending (x : y , z : w) to \n(x^2 : y^2 , z : w).'
        """
        s  = 'Defined by sending '
        s += self.domain().ambient_space()._repr_generic_point()
        s += ' to \n'
        s += self.codomain().ambient_space()._repr_generic_point(self._polys)
        s += '.'
        return s

    def __call__(self, P, check = True):
        r"""
        Make morphisms of products of projective spaces callable.

        INPUT:

        - ``P`` - a point in the domain.

        - ``check`` - Boolean - whether or not to perform the input checks on the image point (Default: True)

        OUTPUT: The image point in the codomain.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProjectiveSpace_cartesian_product([2,1],QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u,y^2*w,z^2*u,w^2,u^2])
            sage: F(T([2,1,3,0,1]))
            (4/9 : 0 : 1 , 0 : 1)
        """
        A = self.codomain()
        f = self.defining_polynomials()
        Q = P[0]._coords+P[1]._coords
        newP = [f(Q) for f in self.defining_polynomials()]
        return(A.point(newP, check))