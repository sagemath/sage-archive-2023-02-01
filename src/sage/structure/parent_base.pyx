###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "../ext/stdsage.pxi"

# TODO: Unpickled parents with base sometimes have thier base set to None.
# This causes a segfault in the module arithmatic architecture.
#
# sage: H = HomsetWithBase(QQ, RR, base=ZZ); H
# sage: H0 = loads(dumps(H))
# sage: H.base_ring(), H0.base_ring()
# (Integer Ring, None)
#
# Perhaps the code below would help (why was it commented out?).

## def make_parent_with_base_v0(_class, _dict, base, has_coerce_map_from):
##     """
##     This should work for any Python class deriving from this, as long
##     as it doesn't implement some screwy __new__() method.
##     """
##     new_object = _class.__new__(_class)
##     if base is None:
##         (<ParentWithBase>new_object)._base = new_object
##     else:
##         (<ParentWithBase>new_object)._base = base
##     (<ParentWithBase>new_object)._has_coerce_map_from = has_coerce_map_from
##     if not _dict is None:
##         new_object.__dict__ = _dict
##     return new_object

def is_ParentWithBase(x):
    """
    Return True if x is a parent object with base.
    """
    return bool(PY_TYPE_CHECK(x, ParentWithBase))

cdef class ParentWithBase(parent.Parent):
    def __init__(self, base):
        self._base = base
        self._has_coerce_map_from = {}

##     def x__reduce__(self):
##         if HAS_DICTIONARY(self):
##             _dict = self.__dict__
##         else:
##             _dict = None
##         if self._base is self:
##             return (make_parent_with_base_v0, (self.__class__, _dict, None, self._has_coerce_map_from))
##         else:
##             return (make_parent_with_base_v0, (self.__class__, _dict, self._base, self._has_coerce_map_from))

    def base_ring(self):
        return self._base

    # Derived class *must* define base_extend.
    def base_extend(self, X):
        raise TypeError, "BUG: the base_extend method must be defined for '%s' (class '%s')"%(
            self, type(self))

    def base(self):
        return self._base

    # Canonical base extension by X (recursive)
    def base_extend_canonical(self, X):
        """
        Returns canonical base extension.
        NOTE: this function should not be extended.
        AUTHOR: Gonzalo Tornaria (2007-06-20)
        """
        return  self.base_extend_canonical_c(X)

    # Canonical base extension by X (recursive)
    cdef base_extend_canonical_c(self, ParentWithBase X):
        """
        AUTHOR: Gonzalo Tornaria (2007-06-20)

        TEST CASES:

        sage: x, y, z, t, r,s = var('x,y,z,t,r,s')
        sage: ZZ.base_extend_canonical(QQ)
        Rational Field
        sage: ZZ[x].base_extend_canonical(QQ)
        Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x].base_extend_canonical(QQ[x])
        Univariate Polynomial Ring in x over Rational Field

        sage: ZZ[x][y].base_extend_canonical(QQ)
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y].base_extend_canonical(QQ[x])
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y].base_extend_canonical(QQ[y])
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y].base_extend_canonical(QQ[x][y])
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

        sage: ZZ[x][y][z].base_extend_canonical(QQ)
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z].base_extend_canonical(QQ[x])
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z].base_extend_canonical(QQ[y])
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z].base_extend_canonical(QQ[z])
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z].base_extend_canonical(QQ[x][y])
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z].base_extend_canonical(QQ[x][z])
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z].base_extend_canonical(QQ[y][z])
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z].base_extend_canonical(QQ[x][y][z])
        Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

        sage: ZZ[x][y][z][t].base_extend_canonical(QQ[x][y][z])
        Univariate Polynomial Ring in t over Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z][t].base_extend_canonical(QQ[y][z])
        Univariate Polynomial Ring in t over Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: ZZ[x][y][z][r][s][t].base_extend_canonical(QQ[r][s][t])
        Univariate Polynomial Ring in t over Univariate Polynomial Ring in s over Univariate Polynomial Ring in r over Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

        # Incorrectly ordered variables will fail.
        # This is consistent with has_coerce_map_from

        sage: ZZ[x][y].base_extend_canonical(QQ[y][x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[z][x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[z][y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[y][x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[y][x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[x][z][y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[y][x][z])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[y][z][x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[z][x][y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x][y][z].base_extend_canonical(QQ[z][y][x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension


        # Some tests for matrix spaces

        sage: MatrixSpace(ZZ,2,2).base_extend_canonical(QQ)
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: MatrixSpace(ZZ[x],2,2).base_extend_canonical(QQ)
        Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
        sage: MatrixSpace(ZZ[x],2,2).base_extend_canonical(QQ[x])
        Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
        sage: MatrixSpace(ZZ[x][y],2,2).base_extend_canonical(QQ)
        Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: MatrixSpace(ZZ[x][y],2,2).base_extend_canonical(QQ[x])
        Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: MatrixSpace(ZZ[x][y],2,2).base_extend_canonical(QQ[y])
        Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
        sage: MatrixSpace(ZZ[x][y],2,2).base_extend_canonical(QQ[x][y])
        Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

        # Some ambiguous extensions

        sage: ZZ[x].base_extend_canonical(ZZ[y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: ZZ[x].base_extend_canonical(QQ[y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: QQ[y].base_extend_canonical(ZZ[x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension
        sage: QQ[y].base_extend_canonical(QQ[x])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension

        sage: ZZ[x][z].base_extend_canonical(ZZ[y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension

        sage: MatrixSpace(ZZ,1,1).base_extend_canonical(ZZ[y])
        Full MatrixSpace of 1 by 1 dense matrices over Univariate Polynomial Ring in y over Integer Ring
        sage: MatrixSpace(ZZ[x],1,1).base_extend_canonical(ZZ[y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension

        sage: FreeModule(ZZ,1).base_extend_canonical(ZZ[y])
        Ambient free module of rank 1 over the integral domain Univariate Polynomial Ring in y over Integer Ring
        sage: FreeModule(ZZ[x],1).base_extend_canonical(ZZ[y])
        Traceback (most recent call last):
        ...
        TypeError: Ambiguous base extension


        """
        #cdef special(t):
        #    return PY_TYPE_CHECK(self, t) and not PY_TYPE_CHECK(X, t)
        #
        #if special(matrix.matrix_space.MatrixSpace_generic) or \
        #   special(sage.modules.free_module.FreeModule_generic):
        #    return self.base_extend(self.base().base_extend_canonical(X))

        # NOTE: if you add more special cases here, make sure you add them down in
        # base_extend_canonical_sym_c()
        from sage.matrix.matrix_space import MatrixSpace_generic
        from sage.modules.free_module import FreeModule_generic
        if (
                PY_TYPE_CHECK(self, MatrixSpace_generic) and not PY_TYPE_CHECK(X, MatrixSpace_generic)
           ) or (
                PY_TYPE_CHECK(self, FreeModule_generic) and not PY_TYPE_CHECK(X, FreeModule_generic)
           ):
            return self.base_extend(self.base().base_extend_canonical(X))

        if X.has_coerce_map_from(self):
            return X
        B = self._base_extend_canonical_rec(X)
        if B is None:
            raise TypeError, "No base extension defined"
        if X._base_extend_canonical_rec(self) is None:
            return B
        raise TypeError, "Ambiguous base extension"

    def base_extend_canonical_sym(self, ParentWithBase X):
        """
        Symmetric version of base_extend_canonical()
        NOTE: this function should not be extended.
        AUTHOR: Gonzalo Tornaria (2007-06-20)
        """
        return  self.base_extend_canonical_sym_c(X)

    cdef base_extend_canonical_sym_c(self, ParentWithBase X):
        from sage.matrix.matrix_space import MatrixSpace_generic
        from sage.modules.free_module import FreeModule_generic

        if (
                PY_TYPE_CHECK(self, MatrixSpace_generic) and not PY_TYPE_CHECK(X, MatrixSpace_generic)
           ) or (
                PY_TYPE_CHECK(self, FreeModule_generic) and not PY_TYPE_CHECK(X, FreeModule_generic)
           ):
                return self.base_extend(self.base().base_extend_canonical(X))

        if (
                not PY_TYPE_CHECK(self, MatrixSpace_generic) and PY_TYPE_CHECK(X, MatrixSpace_generic)
           ) or (
                not PY_TYPE_CHECK(self, FreeModule_generic) and PY_TYPE_CHECK(X, FreeModule_generic)
           ):
                return X.base_extend(X.base().base_extend_canonical(self))

        if self.has_coerce_map_from(X):
            return self
        if X.has_coerce_map_from(self):
            return X

        B1 = self._base_extend_canonical_rec(X)
        B2 = X._base_extend_canonical_rec(self)
        if B1 is None and B2 is None:
            raise TypeError, "No base extension defined"
        if B1 is None:
            return B2
        if B2 is None:
            return B1
        raise TypeError, "Ambiguous base extension"

    def base_extend_recursive(self, X):
        """
        Returns recursive base extension.
        NOTE: this function should not be extended.
        AUTHOR: Gonzalo Tornaria (2007-06-20)
        """
        return  self.base_extend_recursive_c(X)

    # Not sure if this function is ready to use. It should work fine in the non-ambiguous case,
    # but might give some unexpected answers in cases of ambiguity.
    cdef base_extend_recursive_c(self, ParentWithBase X):
        """
        AUTHOR: Gonzalo Tornaria (2007-06-20)
        """
        if X.has_coerce_map_from(self):
            return X
        B = self._base_extend_canonical_rec(X)
        if B is None:
            raise TypeError, "No base extension defined"
        return B

    # Private function for recursion
    cdef _base_extend_canonical_rec(self, ParentWithBase X):
        """
        AUTHOR: Gonzalo Tornaria (2007-06-20)
        """
        cdef ParentWithBase self_b
        self_b = self._base
        if self_b is self:
            return None
        try:
            if X.has_coerce_map_from(self_b):
                return self.base_extend(X)
            else:
                B = X._base_extend_canonical_rec(self_b)
            if B is None:
                B = self_b._base_extend_canonical_rec(X)
                if not B is None:
                    return self.base_extend(B)
            else:
                if not B.has_coerce_map_from(self):
                    return self.base_extend(B)
            B = self_b._base_extend_canonical_rec(X.base())
            if B is None:
                return None
            return self.base_extend(B)

        except TypeError:

            return None


    ############################################################################
    # Homomorphism --
    ############################################################################
    def Hom(self, codomain, cat=None):
        r"""
        self.Hom(codomain, cat=None):

        Return the homspace \code{Hom(self, codomain, cat)} of all
        homomorphisms from self to codomain in the category cat.  The
        default category is \code{self.category()}.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: R.Hom(QQ)
            Set of Homomorphisms from Polynomial Ring in x, y over Rational Field to Rational Field

        Homspaces are defined for very general \sage objects, even elements of familiar rings.
            sage: n = 5; Hom(n,7)
            Set of Morphisms from 5 to 7 in Category of elements of Integer Ring
            sage: z=(2/3); Hom(z,8/1)
            Set of Morphisms from 2/3 to 8 in Category of elements of Rational Field

        This example illustrates the optional third argument:
            sage: QQ.Hom(ZZ, Sets())
            Set of Morphisms from Rational Field to Integer Ring in Category of sets
        """
        try:
            return self._Hom_(codomain, cat)
        except (TypeError, AttributeError):
            pass
        from sage.categories.all import Hom
        return Hom(self, codomain, cat)
