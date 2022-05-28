r"""
Submodules of finite rank free modules

Free modules and submodules of a finite rank free module over a principla ideal
domain have well-defined notion of rank, and they are implemented in
:module:`sage.modules.free_module`. Here submodules with no rank are
implemented. For example, submodules of free modules over multivariate
polynomial rings with more than one variables have no notion of rank.

EXAMPLES::

    sage: S.<x,y,z> = PolynomialRing(QQ)
    sage: M = S**2
    sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
    sage: N
    Submodule of Ambient free module of rank 2 over the integral domain Multivariate Polynomial Ring in x, y, z over Rational Field
    Basis matrix:
    [x - y     z]
    [  y*z   x*z]

AUTHORS:

- Kwankyu Lee (2022-05): initial version

"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.free_module import (basis_seq,
                                      FreeModule_generic_domain,
                                      FreeModule_ambient_domain)

class Submodule_ambient_domain(FreeModule_generic_domain):
    """
    Base class of submodules of ambient free modules over an integral domain.

    INPUT:

    - ``ambient`` -- an ambient free module

    - ``basis`` -- vectors of the ambient free module generating this submodule

    - ``check`` -- boolean; if ``True``, vectors in ``basis`` are checked whether
      they belong to the ambient free module

    A basis of a submodule over an integral domain is just a generating set of
    the submodule. It is not a free basis if there is a nontrivial sygyzy of
    the basis vectors.

    EXAMPLES::

        sage: S.<x,y,z> = PolynomialRing(QQ)
        sage: M = S**2
        sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
        sage: N
        Submodule of Ambient free module of rank 2 over the integral domain Multivariate Polynomial Ring in x, y, z over Rational Field
        Basis matrix:
        [x - y     z]
        [  y*z   x*z]
    """
    def __init__(self, ambient, basis, check=True):
        r"""
        Initialize.

        TESTS::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: N.is_submodule(M)
            True
        """
        if not isinstance(ambient, FreeModule_ambient_domain):
            raise TypeError("ambient (=%s) must be ambient." % ambient)
        self.__ambient_module = ambient
        R = ambient.base_ring()
        R_coord = R

        if check:
            try:
                # convert all basis elements to the ambient module
                basis = [ambient(x) for x in basis]
            except TypeError:
                raise TypeError("each element of basis must be in the ambient free module")

        super().__init__(base_ring=R, rank=None, degree=ambient.degree(), sparse=ambient.is_sparse())

        C = self.element_class
        w = [C(self, x.list(), coerce=False, copy=False) for x in basis]
        self.__basis = basis_seq(self, w)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: M.submodule([vector([x-y,z]), vector([y*z, x*z])])
            Submodule of Ambient free module of rank 2 over the integral domain Multivariate Polynomial Ring in x, y, z over Rational Field
            Basis matrix:
            [x - y     z]
            [  y*z   x*z]
        """
        return "Submodule of %s\n" % self.ambient_module() + \
               "Basis matrix:\n%s" % self.basis_matrix()

    def basis(self):
        """
        Return the basis of this submodule.

        The basis vectors generate this submodule.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: N.basis()
            [
            (x - y, z),
            (y*z, x*z)
            ]
        """
        return self.__basis

    def coordinate_vector(self, v, check=False):
        """
        Return the coordinate vector of ``self``.

        INPUT:

        - ``v`` -- a vector of this submodule

        The coordinate vector of ``v`` is ``v`` itself.

        EXAMPLES::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: M = S**2
            sage: N = M.submodule([vector([x - y, z]), vector([y * z, x * z])])
            sage: v = vector([x - y, z])
            sage: N.coordinate_vector(v)
            (x - y, z)
        """
        if check:
            raise NotImplementedError
        return v

