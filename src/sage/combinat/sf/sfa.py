r"""
Symmetric Functions

AUTHORS:

- Mike Hansen (2007-06-15)

::

    sage: s = SymmetricFunctionAlgebra(QQ, basis='schur')
    sage: e = SymmetricFunctionAlgebra(QQ, basis='elementary')
    sage: f1 = s([2,1]); f1
    s[2, 1]
    sage: f2 = e(f1); f2
    e[2, 1] - e[3]
    sage: f1 == f2
    True
    sage: f1.expand(3, alphabet=['x','y','z'])
    x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2
    sage: f2.expand(3, alphabet=['x','y','z'])
    x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2

::

    sage: m = SFAMonomial(QQ)
    sage: m([3,1])
    m[3, 1]
    sage: m(4)
    4*m[]
    sage: m([4])
    m[4]
    sage: 3*m([3,1])-1/2*m([4])
    3*m[3, 1] - 1/2*m[4]

Code needs to be added to coerce symmetric polynomials into
symmetric functions.

::

    sage: p = SFAPower(QQ)
    sage: m = p(3)
    sage: m
    3*p[]
    sage: m.parent()
    Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
    sage: m + p([3,2])
    3*p[] + p[3, 2]

::

    sage: s = SFASchur(QQ)
    sage: h = SFAHomogeneous(QQ)
    sage: P = SFAPower(QQ)
    sage: e = SFAElementary(QQ)
    sage: m = SFAMonomial(QQ)
    sage: a = s([3,1])
    sage: s(a)
    s[3, 1]
    sage: h(a)
    h[3, 1] - h[4]
    sage: p(a)
    1/8*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] - 1/8*p[2, 2] - 1/4*p[4]
    sage: e(a)
    e[2, 1, 1] - e[2, 2] - e[3, 1] + e[4]
    sage: m(a)
    3*m[1, 1, 1, 1] + 2*m[2, 1, 1] + m[2, 2] + m[3, 1]
    sage: a.expand(4)
    x0^3*x1 + x0^2*x1^2 + x0*x1^3 + x0^3*x2 + 2*x0^2*x1*x2 + 2*x0*x1^2*x2 + x1^3*x2 + x0^2*x2^2 + 2*x0*x1*x2^2 + x1^2*x2^2 + x0*x2^3 + x1*x2^3 + x0^3*x3 + 2*x0^2*x1*x3 + 2*x0*x1^2*x3 + x1^3*x3 + 2*x0^2*x2*x3 + 3*x0*x1*x2*x3 + 2*x1^2*x2*x3 + 2*x0*x2^2*x3 + 2*x1*x2^2*x3 + x2^3*x3 + x0^2*x3^2 + 2*x0*x1*x3^2 + x1^2*x3^2 + 2*x0*x2*x3^2 + 2*x1*x2*x3^2 + x2^2*x3^2 + x0*x3^3 + x1*x3^3 + x2*x3^3

::

    sage: h(m([1]))
    h[1]
    sage: h( m([2]) +m([1,1]) )
    h[2]
    sage: h( m([3]) + m([2,1]) + m([1,1,1]) )
    h[3]
    sage: h( m([4]) + m([3,1]) + m([2,2]) + m([2,1,1]) + m([1,1,1,1]) )
    h[4]
    sage: k = 5
    sage: h( sum([ m(part) for part in Partitions(k)]) )
    h[5]
    sage: k = 10
    sage: h( sum([ m(part) for part in Partitions(k)]) )
    h[10]

::

    sage: P3 = Partitions(3)
    sage: P3.list()
    [[3], [2, 1], [1, 1, 1]]
    sage: m = SFAMonomial(QQ)
    sage: f = sum([m(p) for p in P3])
    sage: m.get_print_style()
    'lex'
    sage: f
    m[1, 1, 1] + m[2, 1] + m[3]
    sage: m.set_print_style('length')
    sage: f
    m[3] + m[2, 1] + m[1, 1, 1]
    sage: m.set_print_style('maximal_part')
    sage: f
    m[1, 1, 1] + m[2, 1] + m[3]
    sage: m.set_print_style('lex')

::

    sage: s = SFASchur(QQ)
    sage: m = SFAMonomial(QQ)
    sage: m([3])*s([2,1])
    2*m[3, 1, 1, 1] + m[3, 2, 1] + 2*m[4, 1, 1] + m[4, 2] + m[5, 1]
    sage: s(m([3])*s([2,1]))
    s[2, 1, 1, 1, 1] - s[2, 2, 2] - s[3, 3] + s[5, 1]
    sage: s(s([2,1])*m([3]))
    s[2, 1, 1, 1, 1] - s[2, 2, 2] - s[3, 3] + s[5, 1]
    sage: e = SFAElementary(QQ)
    sage: e([4])*e([3])*e([1])
    e[4, 3, 1]

::

    sage: s = SFASchur(QQ)
    sage: z = s([2,1]) + s([1,1,1])
    sage: z.coefficient([2,1])
    1
    sage: z.length()
    2
    sage: z.support()
    [[1, 1, 1], [2, 1]]
    sage: z.degree()
    3
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
#*****************************************************************************
from sage.rings.all import Ring, Integer, PolynomialRing, is_Polynomial, is_MPolynomial, ZZ, QQ
from sage.algebras.algebra import Algebra
import sage.combinat.partition
import sage.combinat.skew_partition
import sage.structure.parent_gens
import sage.libs.symmetrica.all as symmetrica
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
from sage.matrix.constructor import matrix
from sage.misc.misc import repr_lincomb, prod, uniq
from sage.algebras.algebra_element import AlgebraElement
import operator
from functools import partial


def SymmetricFunctionAlgebra(R, basis="schur"):
    """
    Return the free algebra over the ring `R` on `n`
    generators with given names.

    INPUT:


    -  ``R`` - ring with identity basis


    OUTPUT: A SymmetricFunctionAlgebra

    EXAMPLES::

        sage: SymmetricFunctionAlgebra(QQ)
        Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis

    ::

        sage: SymmetricFunctionAlgebra(QQ, basis='m')
        Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis

    ::

        sage: SymmetricFunctionAlgebra(QQ, basis='power')
        Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
    """
    if basis == 'schur' or basis == 's':
        return cache_s(R)
    elif basis == "elementary" or  basis ==  'e':
        return cache_e(R)
    elif basis == "homogeneous" or basis ==  'h':
        return cache_h(R)
    elif basis == 'power' or basis ==  'p':
        return cache_p(R)
    elif basis == 'monomial' or basis ==  'm':
        return cache_m(R)
    else:
        raise ValueError, "unknown basis (= %s)"%basis

def SFAPower(R):
    """
    Returns the symmetric function algebra over R with the power-sum
    symmetric functions as the basis.

    EXAMPLES::

        sage: SFAPower(QQ)
        Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='power')


def SFAElementary(R):
    """
    Returns the symmetric function algebra over R with the elementary
    symmetric functions as the basis.

    EXAMPLES::

        sage: SFAElementary(QQ)
        Symmetric Function Algebra over Rational Field, Elementary symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='elementary')


def SFAHomogeneous(R):
    """
    Returns the symmetric function algebra over R with the Homogeneous
    symmetric functions as the basis.

    EXAMPLES::

        sage: SFAHomogeneous(QQ)
        Symmetric Function Algebra over Rational Field, Homogeneous symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='homogeneous')


def SFASchur(R):
    """
    Returns the symmetric function algebra over R with the Schur
    symmetric functions as the basis.

    EXAMPLES::

        sage: SFASchur(QQ)
        Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='schur')


def SFAMonomial(R):
    """
    Returns the symmetric function algebra over R with the monomial
    symmetric functions as the basis.

    EXAMPLES::

        sage: SFAMonomial(QQ)
        Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='monomial')

def is_SymmetricFunctionAlgebra(x):
    """
    Return True if x is a symmetric function algebra; otherwise, return
    False.

    EXAMPLES::

        sage: sage.combinat.sf.sfa.is_SymmetricFunctionAlgebra(5)
        False
        sage: sage.combinat.sf.sfa.is_SymmetricFunctionAlgebra(ZZ)
        False
        sage: sage.combinat.sf.sfa.is_SymmetricFunctionAlgebra(SymmetricFunctionAlgebra(ZZ,'schur'))
        True
    """
    return isinstance(x, SymmetricFunctionAlgebra_generic)



def zee(part):
    """
    Returns the size of the centralizer of permutations of cycle type
    part. Note that this is the inner product between p(part) and
    itself where p is the power-sum symmetric functions.

    INPUT:


    -  ``part`` - an integer partition (for example,
       [2,1,1])


    EXAMPLES::

        sage: from sage.combinat.sf.sfa import zee
        sage: zee([2,1,1])
        4
    """
    if not isinstance(part, sage.combinat.partition.Partition_class):
        part = sage.combinat.partition.Partition_class(part)
    return part.centralizer_size()


def is_SymmetricFunction(x):
    """
    Returns True if x is a symmetric function.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import is_SymmetricFunction
        sage: s = SFASchur(QQ)
        sage: is_SymmetricFunction(2)
        False
        sage: is_SymmetricFunction(s(2))
        True
        sage: is_SymmetricFunction(s([2,1]))
        True
    """
    return isinstance(x, SymmetricFunctionAlgebraElement_generic)

class SymmetricFunctionAlgebra_generic(CombinatorialAlgebra):
    _print_style = 'lex'


    def _change_by_proportionality(self, x, function):
        """
        Return the symmetric function obtained by scaling each basis
        element corresponding to the partition part by f(part).

        INPUT: x: a symmetric function function: a function which takes in
        a partition and returns a scalar

        OUTPUT: a symmetric function in self which is a scaled version of
        x

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([3])+s([2,1])+s([1,1,1]); a
            s[1, 1, 1] + s[2, 1] + s[3]
            sage: f = lambda part: len(part)
            sage: s._change_by_proportionality(a, f)
            3*s[1, 1, 1] + 2*s[2, 1] + s[3]
        """
        BR = self.base_ring()
        z_elt = {}
        for m, c in x._monomial_coefficients.iteritems():
            coeff = function(m)
            z_elt[m] = BR( c*coeff )
        return self._from_dict(z_elt)

    def _change_by_plethysm(self, x, expr, deg_one):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: a = m([2,1])
            sage: a.omega()
            -m[2, 1] - 2*m[3]
            sage: m._change_by_plethysm(-a,-1,[])
            -m[2, 1] - 2*m[3]

        ::

            sage: s = SFASchur(QQ)
            sage: a = s([3])
            sage: s._change_by_plethysm(-a,-1,[])
            s[1, 1, 1]
        """
        #Covert to the power sum
        p = SFAPower(self.base_ring())
        p_x = p(x)
        expr_k = lambda k: expr.subs(**dict([(str(x),x**k) for x in deg_one]))
        f = lambda m,c: (m, c*prod([expr_k(k) for k in m]))
        return self(p_x.map_mc(f))

    # TODO:
    #  - lift to combinatorial_module
    #  - rename to _apply_bimodule_morphism or generalize to true multi_module
    #  - generalization with a "neighbor function" that given x says
    #    for which y one has f(x,y) != 0
    #  - add option orthonormal
    def _apply_multi_module_morphism(self, x, y, f, orthogonal=False):
        """
        INPUT:


        -   x : a element of self

        -```` - y : a element of self

        -   f : a function that takes in two partitions
           (basis elements) and returns an element of the target domain

        -```` - orthogonal: if orthogonal is set to True,
           then f(part1, part2) is assumed to be 0 if part1 != part2.


        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+s([1,1,1])
            sage: b = s([3])+s([2,1])
            sage: f1 = lambda p1, p2: len(p1)*len(p2)
            sage: f2 = lambda p1, p2: len(p1)+len(p2)
            sage: s._apply_multi_module_morphism(a,b,f1,orthogonal=False) #(2+3)*(2+1)
            15
            sage: s._apply_multi_module_morphism(a,b,f1,orthogonal=True)  #(2)*(2)
            4
            sage: s._apply_multi_module_morphism(a,b,f2,orthogonal=False) #2*(2+3+2+1)
            16
            sage: s._apply_multi_module_morphism(a,b,f2,orthogonal=True)  #2+2
            4
        """
        # broken for most coeff ring
        res = 0
        if orthogonal:
            # could check which of x and y has less terms
            # for mx, cx in x:
            for mx, cx in x._monomial_coefficients.iteritems():
                # if not y.has_key(mx):
                if mx not in y._monomial_coefficients:
                    continue
                else:
                    # cy = y[mx]
                    cy = y._monomial_coefficients[mx]
                # might as well call f(mx)
                res += cx*cy*f(mx, mx)
            return res
        else:
            for mx, cx in x._monomial_coefficients.iteritems():
                for my, cy in y._monomial_coefficients.iteritems():
                    res += cx*cy*f(mx,my)
            return res


    def _coerce_impl(self, x):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: m = SFAMonomial(ZZ)
            sage: s._coerce_impl(m([2,1]))
            -2*s[1, 1, 1] + s[2, 1]
        """
        try:
            R = x.parent()
            #Coerce other symmetric functions in
            if is_SymmetricFunctionAlgebra(R):
                #Only perform the coercion if we can go from the base
                #ring of x to the base ring of self
                if self.base_ring().has_coerce_map_from( R.base_ring() ):
                    return self(x)
        except AttributeError:
            pass

        # any ring that coerces to the base ring of this free algebra.
        return self._coerce_try(x, [self.base_ring()])





    def _from_element(self, x):
        """
        Return the element of self with the same 'internal structure' as
        x.

        EXAMPLES::

            sage: e = SFAElementary(QQ)
            sage: s = SFASchur(QQ)
            sage: a = e([2,1]) + e([1,1,1]); a
            e[1, 1, 1] + e[2, 1]
            sage: s._from_element(a)
            s[1, 1, 1] + s[2, 1]
        """
        return self._from_dict(x.monomial_coefficients())

    def _from_cache(self, element, cache_function, cache_dict, **subs_dict):
        """
        Return an element of self from .

        INPUT:


        -  ``element`` - a symmetric function

        -  ``cache_function`` - a function which accepts an
           integer n as its input and creates the cache for that homogenous
           component

        -  ``cache_dict`` - the dictionary storing the cache;
           it is indexed by the positive integers n, and it's values are
           dictionaries indexed by partitions size n. The values of those
           dictionaries are in dictionaries indexed by partitions of size n.

        -  ``subs_dict`` - a dictionary for any substitutions
           to make after the value is extracted from cache_dict.


        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: s = SFASchur(R)
            sage: p21 = Partition([2,1])
            sage: a = s(p21)
            sage: e = SFAElementary(R)
            sage: cache_dict = {}
            sage: cache_dict[3] = {}
            sage: cache_dict[3][p21] = {}
            sage: cache_dict[3][p21][p21] = x^2
            sage: cache_dict[3][p21][Partition([1,1,1])] = 3*x
            sage: cache_function = lambda n: 0 #do nothing
            sage: e._from_cache(a, cache_function, cache_dict)
            3*x*e[1, 1, 1] + x^2*e[2, 1]
            sage: e._from_cache(a, cache_function, cache_dict, x=2)
            6*e[1, 1, 1] + 4*e[2, 1]
        """
        #Convert x to the monomial basis
        BR = self.base_ring()
        zero = BR(0)
        z_elt = {}
        for part,c in element.monomial_coefficients().iteritems():
            if sum(part) not in cache_dict:
                cache_function(sum(part))
            for part2, c2 in cache_dict[sum(part)][part].iteritems():
                z_elt[ part2 ] = z_elt.get(part2, zero) + BR(c*c2.subs(**subs_dict))
        return self._from_dict(z_elt)



    def _invert_morphism(self, n, base_ring, self_to_other_cache, other_to_self_cache,\
                         to_other_function=None, to_self_function=None, \
                         upper_triangular=False, lower_triangular=False, \
                         ones_on_diagonal=False):
        """
        Compute the inverse of a morphism between self and other. In order
        to use this, you must be able compute the morphism in one
        direction. This method assumes that the morphism is indeed
        invertible.

        INPUT:


        -  ``n`` - an integer, the homogeneous component of
           symmetric functions for which we want to a morphism's inverse

        -  ``base_ring`` - the base ring being worked over

        -  ``self_to_other_cache`` - a dictionary which
           stores the transition from self to other

        -  ``other_to_self_cache`` - a dictionary which
           stores the transition from other to self

        -  ``to_other_function`` - a function which takes in
           a partition and returns a function which gives the coefficients of
           self(part) in the other basis

        -  ``to_self_function`` - a function which takes in a
           partition and returns a function which gives the coefficients of
           other(part) in self

        -  ``upper_triangular`` - a boolean, if True, the
           inverse will be computed by back substitution

        -  ``lower_triangular`` - a boolean, if True, the
           inverse will be computed by forward substitution

        -  ``ones_on_diagonal`` - a boolean, if True, the
           entries on the diagonal of the morphism (and inverse) matrix are
           assumed to be one. This is used to remove divisions from the
           forward and back substitute algorithms.


        EXAMPLES: First, we will do an example of inverting the morphism
        which sends a Schur function to its conjugate Schur function. Note
        that this is an involution.

        ::

            sage: s = SFASchur(QQ)
            sage: conj = lambda p1: lambda p2: QQ(1) if p2 == p1.conjugate() else QQ(0)
            sage: c1 = {}
            sage: c2 = {}
            sage: s._invert_morphism(4, QQ, c1, c2, to_other_function = conj)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(c1[4])
            [([1, 1, 1, 1], [([4], 1)]),
             ([2, 1, 1], [([3, 1], 1)]),
             ([2, 2], [([2, 2], 1)]),
             ([3, 1], [([2, 1, 1], 1)]),
             ([4], [([1, 1, 1, 1], 1)])]
            sage: l(c2[4])
            [([1, 1, 1, 1], [([4], 1)]),
             ([2, 1, 1], [([3, 1], 1)]),
             ([2, 2], [([2, 2], 1)]),
             ([3, 1], [([2, 1, 1], 1)]),
             ([4], [([1, 1, 1, 1], 1)])]
            sage: c2 == c1
            True

        We can check that we get the same results if specify
        to_self_function = conj.

        ::

            sage: d1 = {}
            sage: d2 = {}
            sage: s._invert_morphism(4, QQ, d1, d2, to_self_function = conj)
            sage: d1 == c1
            True
            sage: d2 == c2
            True

        Now we do an example of upper triangularity and check that we get
        the same thing whether or not we specify ones_on_diagonal.

        ::

            sage: f = lambda p1: lambda p2: QQ(1) if p2 <= p1 else QQ(0)
            sage: c1 = {}
            sage: c2 = {}
            sage: s._invert_morphism(3, QQ, c1, c2, to_other_function = f, upper_triangular=True)
            sage: l(c1[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], 1), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 1), ([2, 1], 1), ([3], 1)])]
            sage: l(c2[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], -1), ([2, 1], 1)]),
             ([3], [([2, 1], -1), ([3], 1)])]

        ::

            sage: d1 = {}
            sage: d2 = {}
            sage: s._invert_morphism(3, QQ, d1, d2, to_other_function = f,upper_triangular=True, ones_on_diagonal=True)
            sage: c1 == d1
            True
            sage: c2 == d2
            True

        Finally, we do the same thing for lower triangular matrices.

        ::

            sage: f = lambda p1: lambda p2: QQ(1) if p2 >= p1 else QQ(0)
            sage: c1 = {}
            sage: c2 = {}
            sage: s._invert_morphism(3, QQ, c1, c2, to_other_function = f, lower_triangular=True)
            sage: l(c1[3])
            [([1, 1, 1], [([1, 1, 1], 1), ([2, 1], 1), ([3], 1)]),
             ([2, 1], [([2, 1], 1), ([3], 1)]),
             ([3], [([3], 1)])]

        ::

            sage: l(c2[3])
            [([1, 1, 1], [([1, 1, 1], 1), ([2, 1], -1)]),
             ([2, 1], [([2, 1], 1), ([3], -1)]),
             ([3], [([3], 1)])]

        ::

            sage: d1 = {}
            sage: d2 = {}
            sage: s._invert_morphism(3, QQ, d1, d2, to_other_function = f,lower_triangular=True, ones_on_diagonal=True)
            sage: c1 == d1
            True
            sage: c2 == d2
            True
        """
        #Decide whether we know how to go from self to other or
        #from other to self
        if to_other_function is not None:
            known_cache    = self_to_other_cache  #the known direction
            unknown_cache  = other_to_self_cache  #the unknown direction
            known_function = to_other_function
        else:
            unknown_cache  = self_to_other_cache  #the known direction
            known_cache    = other_to_self_cache  #the unknown direction
            known_function = to_self_function

        #Do nothing if we've already computed the inverse
        #for degree n.
        if n in known_cache and n in unknown_cache:
            return

        #Univariate polynomial arithmetic is faster
        #over ZZ.  Since that is all we need to compute
        #the transition matrices between S and P, we
        #should use that.
        #Zt = ZZ['t']
        #t = Zt.gen()
        one  = base_ring(1)
        zero = base_ring(0)

        #Get and store the list of partition we'll need
        pn = sage.combinat.partition.Partitions_n(n).list()
        len_pn = len(pn)

        #Create the initial cache dictionaries
        known_cache_n = {}
        known_matrix_n = matrix(base_ring, len_pn, len_pn)
        unknown_cache_n = {}
        for i in range(len_pn):
            known_cache_part = {}
            f = known_function(pn[i])
            for j in range(len_pn):
                if lower_triangular and j>i:
                    break
                if upper_triangular and i>j:
                    continue
                value = f(pn[j])
                if value != zero:
                    known_cache_part[ pn[ j ] ] = value
                    known_matrix_n[i,j] = value
            known_cache_n[ pn[i] ] = known_cache_part

            unknown_cache_n[ pn[i] ] = {}

        #Compute the inverse of the matrix
        if upper_triangular is not False and lower_triangular is not False:
            raise ValueError, "only one of upper_triangular and lower_triangular can be specified"
        elif upper_triangular is not False:
            #Compute the inverse of by using back
            #substitution.  We solve a len(pn) systems of
            #equations known_matrix_n*x = b_i for x, where e_i
            #is the ith standard basis vector
            inverse = known_matrix_n.parent().zero_matrix()

            delta = lambda i: lambda j: one if i == j else zero

            for column in range(len_pn):
                e = delta(column)
                x = [0]*len_pn
                for i in range(len_pn-1,-1,-1):
                    value = e(i)
                    if not ones_on_diagonal:
                        value /= known_matrix_n[i,i]
                    for j in range(i+1,len_pn):
                        if ones_on_diagonal:
                            value -= known_matrix_n[i,j]*x[j]
                        else:
                            value -= known_matrix_n[i,j]*x[j]/known_matrix_n[i,i]
                    x[i] = value

                for j in range(column+1):
                    if x[j] != zero:
                        inverse[j,column] = x[j]

        elif lower_triangular is not False:
            #Compute the inverse of by using forward
            #substitution.  We solve a len(pn) systems of
            #equations known_matrix_n*x = b_i for x, where e_i
            #is the ith standard basis vector
            inverse = known_matrix_n.parent().zero_matrix()


            delta = lambda i: lambda j: one if i == j else zero

            for column in range(len_pn):
                e = delta(column)
                x = []
                for i in range(len_pn):
                    value = e(i)
                    if not ones_on_diagonal:
                        value /= known_matrix_n[i,i]
                    for j in range(len(x)):
                        if ones_on_diagonal:
                            value -= known_matrix_n[i,j]*x[j]
                        else:
                            value -= known_matrix_n[i,j]*x[j]/known_matrix_n[i,i]
                    x.append(value)
                for j in range(column,len(x)):
                    if x[j] != zero:
                        inverse[j,column] = x[j]

        else:
            inverse = ~known_matrix_n

        for i in range(len_pn):
            for j in range(len_pn):
                if inverse[i,j] != zero:
                    if hasattr(self, '_normalize_coefficients'):
                        unknown_cache_n[ pn[i] ] [ pn[j] ] = self._normalize_coefficients(inverse[i,j])
                    else:
                        unknown_cache_n[ pn[i] ] [ pn[j] ] = inverse[i,j]

        known_cache[ n ]   = known_cache_n
        unknown_cache[ n ] = unknown_cache_n

    def prefix(self):
        """
        Returns the prefix on the elements of self.

        EXAMPLES::

            sage: schur = SFASchur(QQ)
            sage: schur([3,2,1])
            s[3, 2, 1]
            sage: schur.prefix()
            's'
        """
        return self._prefix

    def transition_matrix(self, basis, n):
        """
        Returns the transitions matrix between self and basis for the
        homogenous component of degree n.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: m = SFAMonomial(QQ)
            sage: s.transition_matrix(m,5)
            [1 1 1 1 1 1 1]
            [0 1 1 2 2 3 4]
            [0 0 1 1 2 3 5]
            [0 0 0 1 1 3 6]
            [0 0 0 0 1 2 5]
            [0 0 0 0 0 1 4]
            [0 0 0 0 0 0 1]

        ::

            sage: p = SFAPower(QQ)
            sage: s.transition_matrix(p, 4)
            [ 1/4  1/3  1/8  1/4 1/24]
            [-1/4    0 -1/8  1/4  1/8]
            [   0 -1/3  1/4    0 1/12]
            [ 1/4    0 -1/8 -1/4  1/8]
            [-1/4  1/3  1/8 -1/4 1/24]
            sage: StoP = s.transition_matrix(p,4)
            sage: a = s([3,1])+5*s([1,1,1,1])-s([4])
            sage: a
            5*s[1, 1, 1, 1] + s[3, 1] - s[4]
            sage: mon = a.support()
            sage: coeffs = a.coefficients()
            sage: coeffs
            [5, 1, -1]
            sage: mon
            [[1, 1, 1, 1], [3, 1], [4]]
            sage: cm = matrix([[-1,1,0,0,5]])
            sage: cm * StoP
            [-7/4  4/3  3/8 -5/4 7/24]
            sage: p(a)
            7/24*p[1, 1, 1, 1] - 5/4*p[2, 1, 1] + 3/8*p[2, 2] + 4/3*p[3, 1] - 7/4*p[4]

        ::

            sage: h = SFAHomogeneous(QQ)
            sage: e = SFAElementary(QQ)
            sage: s.transition_matrix(m,7) == h.transition_matrix(s,7).transpose()
            True

        ::

            sage: h.transition_matrix(m, 7) == h.transition_matrix(m, 7).transpose()
            True

        ::

            sage: h.transition_matrix(e, 7) == e.transition_matrix(h, 7)
            True

        ::

            sage: p.transition_matrix(s, 5)
            [ 1 -1  0  1  0 -1  1]
            [ 1  0 -1  0  1  0 -1]
            [ 1 -1  1  0 -1  1 -1]
            [ 1  1 -1  0 -1  1  1]
            [ 1  0  1 -2  1  0  1]
            [ 1  2  1  0 -1 -2 -1]
            [ 1  4  5  6  5  4  1]

        ::

            sage: e.transition_matrix(m,7) == e.transition_matrix(m,7).transpose()
            True
        """
        P = sage.combinat.partition.Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( map( lambda col_part: z.coefficient(col_part), Plist ) )
        return matrix(m)


    def _gram_schmidt(self, n, source, scalar, cache, leading_coeff=None, upper_triangular=True):
        """
        EXAMPLES::

            sage: cache = {}
            sage: from sage.combinat.sf.sfa import zee
            sage: s = SFASchur(QQ)
            sage: m = SFAMonomial(QQ)
            sage: s._gram_schmidt(3, m, zee, cache)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(cache)
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], 2), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 1), ([2, 1], 1), ([3], 1)])]
        """
        BR = self.base_ring(); one = BR(1)
        p = SFAPower(BR)

        #Create a function which converts x and y to the power-sum basis and applies
        #the scalar product.
        pscalar = lambda x,y: p._apply_multi_module_morphism(p(x), p(y), lambda a,b:scalar(a), orthogonal=True)

        if leading_coeff is None:
            leading_coeff = lambda x: one

        #We are going to be doing everything like we are in the upper-triangular case
        #We list the partitions in "decreasing order" and work from the beginning forward.
        #If we are in the lower-triangular case, then we shouldn't reverse the list
        l = sage.combinat.partition.Partitions(n).list()
        if upper_triangular is True:
            l.reverse()

        #precomputed elements
        precomputed_elements = []

        #Handle the initial case
        cache[l[0]] = { l[0]: leading_coeff(l[0]) }
        precomputed_elements.append(leading_coeff( l[0] )*source(l[0]))

        for i in range(1, len(l)):
            start = leading_coeff( l[i] )*source(l[i])
            sub = 0
            for j in range(i):
                sub += pscalar( start, precomputed_elements[j] ) / pscalar(precomputed_elements[j], precomputed_elements[j])*precomputed_elements[j]
            res = start - sub

            if hasattr(self, '_normalize_coefficients'):
                res = res.map_coefficients(self._normalize_coefficients)
            precomputed_elements.append(res)
            cache[l[i]] = {}
            for j in range(i+1):
                cache[l[i]][l[j]] = res.coefficient(l[j])



    def dual_basis(self, scalar=None, scalar_name="", prefix=None):
        """
        Returns the dual basis of self with respect to the scalar product
        scalar. If scalar is None, then the standard (Hall) scalar product
        is used.

        EXAMPLES: The duals of the elementary symmetric functions with
        respect to the Hall scalar product are the forgotten symmetric
        functions.

        ::

            sage: e = SFAElementary(QQ)
            sage: f = e.dual_basis(prefix='f'); f
            Dual basis to Symmetric Function Algebra over Rational Field, Elementary symmetric functions as basis with respect to the Hall scalar product
            sage: f([2,1])^2
            4*f[2, 2, 1, 1] + 6*f[2, 2, 2] + 2*f[3, 2, 1] + 2*f[3, 3] + 2*f[4, 1, 1] + f[4, 2]
            sage: f([2,1]).scalar(e([2,1]))
            1
            sage: f([2,1]).scalar(e([1,1,1]))
            0

        Since the power-sum symmetric functions are orthogonal, their duals
        with respect to the Hall scalar product are scalar multiples of
        themselves.

        ::

            sage: p = SFAPower(QQ)
            sage: q = p.dual_basis(prefix='q'); q
            Dual basis to Symmetric Function Algebra over Rational Field, Power symmetric functions as basis with respect to the Hall scalar product
            sage: q([2,1])^2
            4*q[2, 2, 1, 1]
            sage: p([2,1]).scalar(q([2,1]))
            1
            sage: p([2,1]).scalar(q([1,1,1]))
            0
        """
        import dual
        if scalar is None:
            scalar = zee
            scalar_name = "Hall scalar product"

        return dual.SymmetricFunctionAlgebra_dual(self, scalar, scalar_name, prefix)


    def basis_name(self):
        """
        Returns the name of the basis of self.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s.basis_name()
            'schur'
            sage: p = SFAPower(QQ)
            sage: p.basis_name()
            'power'
            sage: h = SFAHomogeneous(QQ)
            sage: h.basis_name()
            'homogeneous'
            sage: e = SFAElementary(QQ)
            sage: e.basis_name()
            'elementary'
            sage: m = SFAMonomial(QQ)
            sage: m.basis_name()
            'monomial'
        """
        try:
            return self._basis
        except AttributeError:
            return self._prefix.lower()

    def get_print_style(self):
        """
        Returns the value of the current print style for self.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s.get_print_style()
            'lex'
            sage: s.set_print_style('length')
            sage: s.get_print_style()
            'length'
            sage: s.set_print_style('lex')
        """
        return self._print_style

    def set_print_style(self, ps):
        """
        Set the value of the current print style to ps.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s.get_print_style()
            'lex'
            sage: s.set_print_style('length')
            sage: s.get_print_style()
            'length'
            sage: s.set_print_style('lex')
        """
        styles = ['lex', 'length', 'maximal_part']
        if ps not in styles:
            raise ValueError, "the print style must be one of ", styles
        self._print_style = ps



class SymmetricFunctionAlgebraElement_generic(CombinatorialAlgebraElement):
    def __repr__(self):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: f = sum([m(p) for p in Partitions(3)])
            sage: m.get_print_style()
            'lex'
            sage: f
            m[1, 1, 1] + m[2, 1] + m[3]
            sage: m.set_print_style('length')
            sage: f
            m[3] + m[2, 1] + m[1, 1, 1]
            sage: m.set_print_style('maximal_part')
            sage: f
            m[1, 1, 1] + m[2, 1] + m[3]
            sage: m.set_print_style('lex')
        """
        v = self._monomial_coefficients.items()

        ps = self.parent().get_print_style()
        if ps == 'lex':
            v.sort(key=lambda x: x[0])
        if ps == 'length':
            v.sort(key=lambda x: len(x[0]))
        if ps == 'maximal_part':
             v.sort(key=_lmax)

        prefix = self.parent().prefix()
        mons = [ prefix + repr(m) for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def _latex_(self):
        """
        Returns a string representing the LaTeX version of self.

        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: f = sum([m(p) for p in Partitions(3)])
            sage: m.get_print_style()
            'lex'
            sage: latex(f) #indirect doctest
            m_{1,1,1} + m_{2,1} + m_{3}
            sage: m.set_print_style('length')
            sage: latex(f)
            m_{3} + m_{2,1} + m_{1,1,1}
            sage: m.set_print_style('maximal_part')
            sage: latex(f)
            m_{1,1,1} + m_{2,1} + m_{3}
        """
        v = self._monomial_coefficients.items()

        ps = self.parent().get_print_style()
        if ps == 'lex':
            v.sort(key=lambda x: x[0])
        if ps == 'length':
            v.sort(key=lambda x: len(x[0]))
        if ps == 'maximal_part':
            v.sort(key=_lmax)

        prefix = self.parent().prefix()
        mons = [ prefix + '_{' + ",".join(map(str, m)) + '}' for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs, is_latex=True).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def plethysm(self, x, include=None, exclude=None):
        """
        Returns the outer plethysm of self with x.

        By default, the degree one elements are the generators for the
        self's base ring.

        INPUT:


        -  ``x`` - a symmetric function

        -  ``include`` - a list of variables to be treated as
           degree one elements instead of the default degree one elements

        -  ``exclude`` - a list of variables to be excluded
           from the default degree one elements


        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: h = SFAHomogeneous(QQ)
            sage: s ( h([3])( h([2]) ) )
            s[2, 2, 2] + s[4, 2] + s[6]
            sage: p = SFAPower(QQ)
            sage: p([3])( s([2,1]) )
            1/3*p[3, 3, 3] - 1/3*p[9]
            sage: e = SFAElementary(QQ)
            sage: e([3])( e([2]) )
            e[3, 3] + e[4, 1, 1] - 2*e[4, 2] - e[5, 1] + e[6]

        ::

            sage: R.<t> = QQ[]
            sage: s = SFASchur(R);
            sage: a = s([3])
            sage: f = t*s([2])
            sage: a(f)
            t^3*s[2, 2, 2] + t^3*s[4, 2] + t^3*s[6]
            sage: f(a)
            t*s[4, 2] + t*s[6]
        """
        if not is_SymmetricFunction(x):
            raise TypeError, "only know how to compute plethysms between symmetric functions"
        parent = self.parent()
        R = parent.base_ring()
        p = SFAPower(R)
        p_x = p(x)

        #Handle degree one elements
        if include is not None and exclude is not None:
            raise RuntimeError, "include and exclude cannot both be specified"

        try:
            degree_one = [R(g) for g in R.variable_names_recursive()]
        except AttributeError:
            try:
                degree_one = R.gens()
            except NotImplementedError:
                degree_one= []

        if include:
            degree_one = [R(g) for g in include]
        if exclude:
            degree_one = [g for g in degree_one if g not in exclude]


        #Takes in n, and returns a function which takes in a partition and
        #scales all of the parts of that partition by n
        scale_part = lambda n: lambda m: m.__class__([i*n for i in m])

        raise_c = lambda n: lambda c: c.subs(**dict((str(g),g**n) for g in degree_one if g != 1))

        #Takes n an symmetric function f, and an n and returns the
        #symmetric function with all of its basis partitions scaled
        #by n
        pn_pleth = lambda f, n: f.map_support(scale_part(n))

        #Takes in a partition and applies
        f = lambda part: prod( pn_pleth(p_x.map_coefficients(raise_c(i)), i) for i in part )
        return parent(p._apply_module_morphism(p(self),f))

    __call__ = plethysm


    def _inner_plethysm_pk_g(self, k, g, cache):
        r"""
        Returns the inner plethysm between `p_k` and `g`.

        INPUT:


        -  ``k`` - a positive integer

        -  ``g`` - a symmetric function in the power sum basis

        -  ``cache`` - a dictionary whose keys are (k, g) pairs
           and values are the cached output of this function


        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: _inner_plethysm_pk_g = p(0)._inner_plethysm_pk_g
            sage: _inner_plethysm_pk_g(2, p([1,1,1]), {})
            p[1, 1, 1] + 3*p[2, 1]
            sage: _inner_plethysm_pk_g(5, p([2,2,1,1,1]), {})
            p[2, 2, 1, 1, 1]
        """
        try:
            return cache[(k,g)]
        except KeyError:
            pass

        p = SFAPower(QQ)
        res = 0
        degrees = uniq([ sum(m) for m in g.support() ])
        for d in degrees:
            for mu in sage.combinat.partition.Partitions(d):
                mu_k = mu.power(k)
                if mu_k in g:
                    res += g.coefficient(mu_k)*mu_k.centralizer_size()/mu.centralizer_size()*p(mu)

        cache[(k,g)] = res
        return res

    def _inner_plethysm_pnu_g(self, p_x, cache, nu):
        """
        Returns the inner plethysm of p(nu) with another symmetric function
        p_x in the power-sum basis.

        Note that the order of the arguments is somewhat strange in order
        to facilitate partial function application.

        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: s = SFASchur(QQ)
            sage: _inner_plethysm_pnu_g = p(0)._inner_plethysm_pnu_g
            sage: _inner_plethysm_pnu_g( p([1,1,1]), {}, Partition([2,1]))
            6*p[1, 1, 1]
            sage: _inner_plethysm_pnu_g( p([1,1,1]), {}, Partition([]))
            1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]
            sage: s(_)
            s[3]
        """
        #We handle the constant term case separately.  It should be
        #the case that p([]).inner_tensor(s(mu)) = s([ mu.size() ]).
        #Here, we get the degrees of the homogeneous pieces of
        if len(nu) == 0:
            s = SFASchur(self.base_ring())
            p = SFAPower(self.base_ring())
            degrees = [ part.size() for part in p_x.support() ]
            degrees = uniq(degrees)
            return p(sum([s([n]) for n in degrees]))

        #For each k in nu, we compute the inner plethysm of
        #p_k with p_x
        res = [self._inner_plethysm_pk_g(k, p_x, cache) for k in nu]

        #To get the final answer, we compute the inner tensor product
        #of all the symmetric functions in res
        return reduce(lambda x, y: x.itensor(y), res)

    def inner_plethysm(self, x):
        r"""
        Retuns the inenr plethysm of self with x.

        The result of f.inner_plethysm(g) is linear in f and linear in
        'homogeneous pieces' of g. So, to describe this function, we assume
        without loss that f is some Schur function s(la) and g is a
        homogeneous symmetric function of degree n. The function g can be
        thought of as the character of an irreducible representation, rho,
        of the symmetric group `S_n`. Let N be the dimension of
        this representation. If the number of parts of la is greater then
        N, then f.inner_plethysm(g) = 0 by definition. Otherwise, we can
        interpret f as the character of an irreducible `GL_N`
        representation, call it `\sigma`. Now
        `\sigma \circ \rho` is an `S_n` representation
        and, by definition, the character of this representation is
        f.inner_plethysm(g).

        REFERENCES:

        - King, R. 'Branching rules for `GL_m \supset \Sigma_n`
          and the evaluation of inner plethysms.' J. Math. Phys. 15,
          258 (1974)

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: p = SFAPower(QQ)
            sage: h = SFAHomogeneous(QQ)
            sage: s([2,1]).inner_plethysm(s([1,1,1]))
            0
            sage: h([2]).inner_plethysm(s([2,1]))
            h[2, 1]
            sage: s(_)
            s[2, 1] + s[3]

        ::

            sage: f = s([2,1]) + 2*s([3,1])
            sage: f.itensor(f)
            s[1, 1, 1] + s[2, 1] + 4*s[2, 1, 1] + 4*s[2, 2] + s[3] + 4*s[3, 1] + 4*s[4]
            sage: s( h([1,1]).inner_plethysm(f) )
            s[1, 1, 1] + s[2, 1] + 4*s[2, 1, 1] + 4*s[2, 2] + s[3] + 4*s[3, 1] + 4*s[4]

        ::

            sage: s([]).inner_plethysm(s([1,1]) + 2*s([2,1])+s([3]))
            s[2] + s[3]
            sage: [s([]).inner_plethysm(s(p)) for p in Partitions(4)]
            [s[4], s[4], s[4], s[4], s[4]]
        """
        parent = self.parent()
        base_ring = parent.base_ring()
        p = SFAPower(base_ring)

        p_x = p(x)

        cache = {}
        f = partial(self._inner_plethysm_pnu_g, p_x, cache)

        return parent(parent._apply_module_morphism(p(self), f))


    def omega(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.
        The default implementation converts to the Schurs performs the
        automorphism and changes back.

        EXAMPLES::

            sage: J = JackPolynomialsP(QQ,1)
            sage: a = J([2,1]) + J([1,1,1])
            sage: a.omega()
            JackP[2, 1] + JackP[3]
        """
        parent = self.parent()
        s = SFASchur(parent.base_ring())
        return parent(s(self).omega())

    def theta(self,a):
        """
        Returns the image of self under the theta automorphism which sends
        `p[k]` to `a*p[k]`.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s([2,1]).theta(2)
            2*s[1, 1, 1] + 6*s[2, 1] + 2*s[3]
            sage: p = SFAPower(QQ)
            sage: p([2]).theta(2)
            2*p[2]
        """
        p = SFAPower(self.parent().base_ring())
        p_self = p(self)
        res = p_self.map_mc(lambda m,c: (m, c*a**len(m)))
        return self.parent()(res)

    def theta_qt(self,q,t):
        """
        Returns the image of self under the theta automorphism which sends
        `p[k]` to `(1-q^k)/(1-t^k)*p[k]`.

        EXAMPLES::

            sage: QQqt = QQ['q,t'].fraction_field()
            sage: q,t = QQqt.gens()
            sage: p = SFAPower(QQqt)
            sage: p([2]).theta_qt(q,t)
            ((-q^2+1)/(-t^2+1))*p[2]
            sage: p([2,1]).theta_qt(q,t)
            ((q^3-q^2-q+1)/(t^3-t^2-t+1))*p[2, 1]
        """
        BR = self.parent().base_ring()
        p = SFAPower(BR)
        p_self = p(self)
        res = p_self.map_mc(lambda m,c: (m, BR(prod([(1-q**k)/(1-t**k) for k in m])*c)))
        return self.parent()(res)


    def itensor(self, x):
        """
        Returns the inner tensor product of self and x in the basis of
        self.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: b = s([3])
            sage: a.itensor(b)
            s[2, 1]
            sage: c = s([3,2,1])
            sage: c.itensor(c)
            s[1, 1, 1, 1, 1, 1] + 2*s[2, 1, 1, 1, 1] + 3*s[2, 2, 1, 1] + 2*s[2, 2, 2] + 4*s[3, 1, 1, 1] + 5*s[3, 2, 1] + 2*s[3, 3] + 4*s[4, 1, 1] + 3*s[4, 2] + 2*s[5, 1] + s[6]

        TESTS::

            sage: s = SFASchur(QQ)
            sage: a = s([8,8])
            sage: a.itensor(a) #long
            s[4, 4, 4, 4] + s[5, 5, 3, 3] + s[5, 5, 5, 1] + s[6, 4, 4, 2] + s[6, 6, 2, 2] + s[6, 6, 4] + s[7, 3, 3, 3] + s[7, 5, 3, 1] + s[7, 7, 1, 1] + s[8, 4, 2, 2] + s[8, 4, 4] + s[8, 6, 2] + s[8, 8] + s[9, 3, 3, 1] + s[9, 5, 1, 1] + s[10, 2, 2, 2] + s[10, 4, 2] + s[10, 6] + s[11, 3, 1, 1] + s[12, 2, 2] + s[12, 4] + s[13, 1, 1, 1] + s[14, 2] + s[16]
        """
        #Convert both self and x to the p basis
        p = SFAPower(self.parent().base_ring())
        f = lambda part1, part2: zee(part1)*p(part1)
        return self.parent()(p._apply_multi_module_morphism(p(self),p(x),f,orthogonal=True))


    internal_product = itensor
    kronecker_product = itensor
    inner_tensor = itensor

    def scalar_t(self, x, t=None):
        """
        Returns the standard Hall-Littlewood scalar product of self and x.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: sp = a.scalar_t(a); sp
            (-t^2 - 1)/(t^5 - 2*t^4 + t^3 - t^2 + 2*t - 1)
            sage: sp.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        p = SFAPower(self.parent().base_ring())
        if t is None:
            t = QQ['t'].gen()
        f = lambda part1, part2: part1.centralizer_size(t=t, q=ZZ(0))
        return p._apply_multi_module_morphism(p(self), p(x), f, orthogonal=True)


    scalar_hl = scalar_t

    def derivative_with_respect_to_p1(self, n=1):
        """
        Returns the symmetric function obtained by taking the derivative of
        self with respect to the power-sum symmetric function p([1]) when
        the expansion of self in the power-sum basis is considered as a
        polynomial in p([1])'s.

        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: a = p([1,1,1])
            sage: a.derivative_with_respect_to_p1()
            3*p[1, 1]
            sage: a.derivative_with_respect_to_p1(1)
            3*p[1, 1]
            sage: a.derivative_with_respect_to_p1(2)
            6*p[1]
            sage: a.derivative_with_respect_to_p1(3)
            6*p[]

        ::

            sage: s = SFASchur(QQ)
            sage: s([3]).derivative_with_respect_to_p1()
            s[2]
            sage: s([2,1]).derivative_with_respect_to_p1()
            s[1, 1] + s[2]
            sage: s([1,1,1]).derivative_with_respect_to_p1()
            s[1, 1]
        """
        p = SFAPower(self.parent().base_ring())
        res = p(self)
        for i in range(n):
            res = res._derivative_with_respect_to_p1()
        return self.parent()(res)


    def _expand(self, condition, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n
        variables.

        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: a = p([2])+p([3])
            sage: a._expand(lambda part: False, 3)
            x0^3 + x1^3 + x2^3 + x0^2 + x1^2 + x2^2
            sage: a._expand(lambda part: max(part)>2, 3)
            x0^2 + x1^2 + x2^2
        """
        import classical
        parent = self.parent()
        e = eval('symmetrica.compute_' + str(classical.translate[parent.basis_name()]).lower() + '_with_alphabet')
        resPR = PolynomialRing(parent.base_ring(), n, alphabet)
        f = lambda part: resPR(0) if condition(part) else resPR(e(part, n, alphabet))
        return parent._apply_module_morphism(self, f)

    def scalar(self, x):
        """
        Returns standard scalar product between self and s.

        This is the default implementation that converts both self and x
        into Schur functions and performs the scalar product that basis.

        EXAMPLES::

            sage: e = SFAElementary(QQ)
            sage: h = SFAHomogeneous(QQ)
            sage: m = SFAMonomial(QQ)
            sage: p4 = Partitions(4)
            sage: matrix([ [e(a).scalar(h(b)) for a in p4] for b in p4])
            [ 0  0  0  0  1]
            [ 0  0  0  1  4]
            [ 0  0  1  2  6]
            [ 0  1  2  5 12]
            [ 1  4  6 12 24]
            sage: matrix([ [h(a).scalar(e(b)) for a in p4] for b in p4])
            [ 0  0  0  0  1]
            [ 0  0  0  1  4]
            [ 0  0  1  2  6]
            [ 0  1  2  5 12]
            [ 1  4  6 12 24]
            sage: matrix([ [m(a).scalar(e(b)) for a in p4] for b in p4])
            [-1  2  1 -3  1]
            [ 0  1  0 -2  1]
            [ 0  0  1 -2  1]
            [ 0  0  0 -1  1]
            [ 0  0  0  0  1]
            sage: matrix([ [m(a).scalar(h(b)) for a in p4] for b in p4])
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
        """
        sp = self.parent()
        xp = x.parent()
        BR = sp.base_ring()

        s = SFASchur(BR)
        s_self = s(self)
        s_x = s(x)
        return s_self.scalar(s_x)


    def is_schur_positive(self):
        """
        Returns True if and only if self is Schur positive. If s is the
        space of Schur functions over self's base ring, then this is the
        same as self._is_positive(s).

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: a.is_schur_positive()
            True
            sage: a = s([2,1]) - s([3])
            sage: a.is_schur_positive()
            False

        ::

            sage: QQx = QQ['x']
            sage: s = SFASchur(QQx)
            sage: x = QQx.gen()
            sage: a = (1+x)*s([2,1])
            sage: a.is_schur_positive()
            True
            sage: a = (1-x)*s([2,1])
            sage: a.is_schur_positive()
            False
        """
        return self._is_positive( SFASchur(self.parent().base_ring()) )


    def _is_positive(self, s):
        """
        Returns True if and only if self has nonnegative coefficients in
        the basis s.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: a._is_positive(s)
            True
            sage: a = s([2,1]) - s([3])
            sage: a._is_positive(s)
            False
        """
        s_self = s(self)
        return all([ _nonnegative_coefficients(c) for c in s_self.coefficients()])

    def degree(self):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1]) + 3
            sage: z.degree()
            4
        """
        return max( map( sum, self._monomial_coefficients ) + [0] )

    def restrict_degree(self, d, exact=True):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_degree(2)
            0
            sage: z.restrict_degree(1)
            s[1]
            sage: z.restrict_degree(3)
            s[1, 1, 1] + s[2, 1]
            sage: z.restrict_degree(3, exact=False)
            s[1] + s[1, 1, 1] + s[2, 1]
        """
        if exact:
            res = dict( filter( lambda x: sum(x[0]) == d, self._monomial_coefficients.items()) )
        else:
            res = dict( filter( lambda x: sum(x[0]) <= d, self._monomial_coefficients.items()) )
        return self.parent()._from_dict(res)

    def restrict_partition_lengths(self, l, exact=True):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_partition_lengths(2)
            s[2, 1]
            sage: z.restrict_partition_lengths(2, exact=False)
            s[1] + s[2, 1] + s[4]
        """
        if exact:
            res = dict( filter( lambda x: len(x[0]) == l, self._monomial_coefficients.items()) )
        else:
            res = dict( filter( lambda x: len(x[0]) <= l, self._monomial_coefficients.items()) )
        return self.parent()._from_dict(res)

    def restrict_parts(self, n):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_parts(2)
            s[1] + s[1, 1, 1] + s[2, 1]
            sage: z.restrict_parts(1)
            s[1] + s[1, 1, 1]
        """
        res = dict( filter( lambda x: _lmax(x) <= n, self._monomial_coefficients.items()) )
        return self.parent()._from_dict(res)

    def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n
        variables.

        EXAMPLES::

            sage: J =JackPolynomialsJ(QQ, t=2)
            sage: J([2,1]).expand(3)
            4*x0^2*x1 + 4*x0*x1^2 + 4*x0^2*x2 + 6*x0*x1*x2 + 4*x1^2*x2 + 4*x0*x2^2 + 4*x1*x2^2
        """
        s = SFASchur(self.parent().base_ring())
        condition = lambda part: len(part) > n
        return s(self)._expand(condition, n, alphabet)

    def skew_by(self, x):
        """
        Returns the element whose result is the dual to multiplication by x
        applied to self.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s([3,2]).skew_by(s([2]))
            s[2, 1] + s[3]
            sage: s([3,2]).skew_by(s([1,1,1]))
            0
            sage: s([3,2,1]).skew_by(s([2,1]))
            s[1, 1, 1] + 2*s[2, 1] + s[3]

        ::

            sage: p = SFAPower(QQ)
            sage: p([4,3,3,2,2,1]).skew_by(p([2,1]))
            4*p[4, 3, 3, 2]
            sage: zee = sage.combinat.sf.sfa.zee
            sage: zee([4,3,3,2,2,1])/zee([4,3,3,2])
            4
        """
        s = SFASchur(self.parent().base_ring())
        f = lambda part1, part2: s([part1,part2]) if part1.contains(part2) else 0
        return self.parent()(s._apply_multi_module_morphism(s(self),s(x),f))



    def hl_creation_operator(self, nu):
        """
        This is the vertex operator that generalizes Jing's operator It is
        from: Hall-Littlewood Vertex Operators and Kostka Polynomials,
        Shimizono-Zabrocki, Proposition 5 It is a linear operator that
        rases the degree by sum(nu) This creation operator is a t-analogue
        of multiplication by s(nu)

        INPUT:


        -  ``nu`` - a partition


        EXAMPLES::

            sage: s = SFASchur(QQ['t'])
            sage: s([2]).hl_creation_operator([3,2])
            s[3, 2, 2] + t*s[3, 3, 1] + t*s[4, 2, 1] + t^2*s[4, 3] + t^2*s[5, 2]
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLQp(s([2]).hl_creation_operator([2]).hl_creation_operator([3]))
            Qp[3, 2, 2]
            sage: s([2,2]).hl_creation_operator([2,1])
            t*s[2, 2, 2, 1] + t^2*s[3, 2, 1, 1] + t^2*s[3, 2, 2] + t^3*s[3, 3, 1] + t^3*s[4, 2, 1] + t^4*s[4, 3]
            sage: s(1).hl_creation_operator([2,1,1])
            s[2, 1, 1]
            sage: s(0).hl_creation_operator([2,1,1])
            0
            sage: s([3,2]).hl_creation_operator([2,1,1])
            (t^2-t)*s[2, 2, 2, 2, 1] + t^3*s[3, 2, 2, 1, 1] + (t^3-t^2)*s[3, 2, 2, 2] + t^3*s[3, 3, 1, 1, 1] + t^4*s[3, 3, 2, 1] + t^3*s[4, 2, 1, 1, 1] + t^4*s[4, 2, 2, 1] + 2*t^4*s[4, 3, 1, 1] + t^5*s[4, 3, 2] + t^5*s[4, 4, 1] + t^4*s[5, 2, 1, 1] + t^5*s[5, 3, 1]
        """
        s = SFASchur(self.base_ring())

        t = self.base_ring().gen()

        res = self*s(nu)
        for d in range(self.degree()):
            for mu in sage.combinat.partition.Partitions(d+1).filter(lambda p: len(p) <= len(nu)):
                for lam, c in s(mu)*s(nu):
                    if len(lam) <= len(nu):
                        res += c*s(lam)*self.skew_by(s(mu).plethysm((t-1)*s([1])))
        return res


#############
#   Cache   #
#############
from sage.misc.cache import Cache
import schur, monomial, powersum, elementary, homogeneous
cache_s = Cache(schur.SymmetricFunctionAlgebra_schur)
cache_m = Cache(monomial.SymmetricFunctionAlgebra_monomial)
cache_p = Cache(powersum.SymmetricFunctionAlgebra_power)
cache_e = Cache(elementary.SymmetricFunctionAlgebra_elementary)
cache_h = Cache(homogeneous.SymmetricFunctionAlgebra_homogeneous)


###################
def _lmax(x):
    """
    Returns the max of x[0] where x is a list. If x is the empty list,
    the _lmax returns 0.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import _lmax
        sage: _lmax(([3,2,1], 2))
        3
        sage: _lmax(([],4))
        0
    """
    if x[0] == []:
        return 0
    else:
        return max(x[0])

def _nonnegative_coefficients(x):
    """
    Returns True if x has nonnegative coefficients.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import _nonnegative_coefficients
        sage: _nonnegative_coefficients(2)
        True
        sage: _nonnegative_coefficients(-2)
        False
        sage: R.<x> = ZZ[]
        sage: _nonnegative_coefficients(x^2+4)
        True
        sage: _nonnegative_coefficients(x^2-4)
        False
    """
    if is_Polynomial(x) or is_MPolynomial(x):
        return all([ c >= 0 for c in x.coeffs() ])
    else:
        return x >= 0
