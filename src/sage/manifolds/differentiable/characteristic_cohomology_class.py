r"""
Characteristic cohomology classes
"""

#******************************************************************************
#       Copyright (C) 2021 Michael Jung <m.jung at vu.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.algebras.finite_gca import FiniteGCAlgebra
from sage.combinat.free_module import IndexedFreeModuleElement
from sage.misc.fast_methods import Singleton
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from .affine_connection import AffineConnection
from .bundle_connection import BundleConnection
from .levi_civita_connection import LeviCivitaConnection


class CharacteristicCohomologyClass(IndexedFreeModuleElement):
    r"""

    """

    def __init__(self, parent, x, name=None, latex_name=None):
        r"""

        """
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._mixed_forms = {}  # dict. of characteristic forms of `self`
                                # (key: bundle connection)
        super().__init__(parent, x)

    def _repr_(self):
        r"""

        """
        if self._name is None:
            name = super()._repr_()
        else:
            name = self._name
        vbundle = self.parent()._vbundle
        name = f'({name})({vbundle._name})'
        return f'Characteristic cohomology class {name} of the {vbundle}'

    def _latex_(self):
        r"""

        """
        if self._latex_name is None:
            latex = super()._latex_()
        else:
            latex = self._latex_name
        vbundle = self.parent()._vbundle
        latex = r'\left(' + latex + r'\right)\right('
        latex += vbundle._latex_name + r'\right)'
        return latex

    def get_form(self, nab):
        r"""

        """
        if nab not in self._mixed_forms:
            dom = nab._domain
            A = dom.mixed_form_algebra()

            # trivial cases
            if self == 1:
                self._mixed_forms[nab] = A(dom._one_scalar_field)
            elif self == 0:
                self._mixed_forms[nab] = A(dom._zero_scalar_field)
            else:  # non-trivial case
                from functools import reduce
                from itertools import chain

                parent = self.parent()

                gen_algorithm = parent._algorithm  # this is a list
                # concatenate generators
                gen_forms = list(chain.from_iterable(a.get(nab)
                                                     for a in gen_algorithm))
                grading = parent.print_options()['sorting_key']
                res = [dom.diff_form_module(i).zero()
                       for i in range(dom._dim + 1)]
                for ind, c in self:
                    deg = grading(ind)
                    gen_pow = [fast_wedge_power(f, i)
                               for f, i in zip(gen_forms, ind)]
                    res[deg] += c * reduce(lambda x, y: x.wedge(y), gen_pow)

                res = A(res)  # convert result into mixed form

                # preparse names
                vbundle = parent._vbundle
                if self._name is None:
                    name = f'({super()._repr_()})'
                else:
                    name = f'({self._name})'
                if self._latex_name is None:
                    latex_name = r'\left(' + super()._latex_() + r'\right)'
                else:
                    latex_name = r'\left(' + self._latex_name + r'\right)'
                # appendix
                append_name = f'({vbundle._name}, {nab._name})'
                append_latex_name = r'\left(' + vbundle._latex_name
                append_latex_name += ', ' + nab._latex_name + r'\right)'

                # set names of components
                from sage.arith.misc import gcd

                step = gcd(parent._degrees)  # step size of (possibly) non-zero
                for i in range(dom._dim // step + 1):
                    # enumerate (possibly) non-zero components
                    comp_name = name + f'_{i}' + append_name
                    comp_latex_name = latex_name + r'_{' + str(i) + '}'
                    comp_latex_name += append_latex_name
                    res[step * i].set_name(name=comp_name,
                                           latex_name=comp_latex_name)

                # set global names
                res._name = name + append_name
                res._latex_name = latex_name + append_latex_name

                res.set_immutable()

                self._mixed_forms[nab] = res  # cache result in dict

        return self._mixed_forms[nab]

    representative = get_form


class CharacteristicCohomologyClassRing(FiniteGCAlgebra):
    r"""

    """
    Element = CharacteristicCohomologyClass

    def __init__(self, base, vbundle):
        r"""

        """
        self._vbundle = vbundle
        self._domain = vbundle._base_space
        dim = self._domain._dim
        rk = vbundle._rank
        if vbundle._field_type == 'complex':
            ran = min(rk, dim // 2)
            names = [f'c_{i}' for i in range(1, ran + 1)]
            degrees = [2 * i for i in range(1, ran + 1)]
            self._algorithm = [ChernAlgorithm()]
        elif vbundle._field_type == 'real':
            ran = min(rk // 2, dim // 4)
            names = [f'p_{i}' for i in range(1, ran + 1)]
            degrees = [4 * i for i in range(1, ran + 1)]
            self._algorithm = [PontryaginAlgorithm()]
            if vbundle.has_orientation():
                # add Euler class generator
                names += ['e']
                degrees += [rk]
                self._algorithm += [EulerAlgorithm()]
                # TODO: add relation e^2=p_k for dim=2*k
        else:
            raise TypeError(f'Characteristic cohomology classes not supported '
                            f'for vector bundles with '
                            f'field type {vbundle._field_type}')

        names = tuple(names)  # hashable
        degrees = tuple(degrees)  # hashable
        super().__init__(base=base, names=names, degrees=degrees,
                         max_degree=dim)

    def _element_constructor_(self, x, name=None, latex_name=None):
        r"""
        Convert ``x`` into ``self``.
        """
        R = self.base_ring()

        if x in R:
            one_basis = self.one_basis()
            d = {one_basis: R(x)}
        elif isinstance(x, CharacteristicCohomologyClass):
            d = x._monomial_coefficients
        # x is an element of the basis enumerated set;
        # This is a very ugly way of testing this
        elif ((hasattr(self._indices, 'element_class') and
               isinstance(self._indices.element_class, type) and
               isinstance(x, self._indices.element_class)) or
              self.parent()(x) == self._indices):
            d = {x: R.one()}
        elif x in self._indices:
            d = {self._indices(x): R.one()}
        else:
            raise TypeError(f"do not know how to make x (= {x}) "
                            f"an element of self (={self})")

        return self.element_class(self, d, name=name, latex_name=latex_name)


# *****************************************************************************
# ALGORITHMS
# *****************************************************************************

def fast_wedge_power(form, n):
    r"""
    Return the wedge product power of `form` using a square-and-wedge
    algorithm.
    """
    if n == 0:
        return form._domain._one_scalar_field
    elif n < 0:
        raise ValueError("'n' must be non-negative")
    val = form
    while not (n & 1):
        print(n)
        val = val.wedge(val)
        n >>= 1

    # Now multiply together the correct factors form^(2^i)
    res = val
    n >>= 1
    while n:
        val = val.wedge(val)
        if n & 1:
            res = val.wedge(res)
        n >>= 1

    return res


class Algorithm_generic(SageObject):
    r"""
    Algorithm class to compute the characteristic forms of the generators.
    """

    @cached_method
    def get(self, nab):
        r"""
        Return the global characteristic forms of the generators w.r.t. a given
        connection.
        """
        if isinstance(nab, AffineConnection):
            vbundle = nab._domain.tangent_bundle()
        elif isinstance(nab, BundleConnection):
            vbundle = nab._vbundle
        else:
            raise TypeError(f'{nab} must be a connection')
        dom = nab._domain
        res = []  # will be specified within first iteration
        for frame in dom._get_min_covering(nab._coefficients):
            cmatrix = [[nab.curvature_form(i, j, frame)
                        for j in vbundle.irange()]
                       for i in vbundle.irange()]
            res_loc = self.get_local(cmatrix)
            if not res:
                # until now, degrees of generators were unknown
                res = [dom.diff_form(loc_form.degree())
                       for loc_form in res_loc]
            for form, loc_form in zip(res, res_loc):
                form.set_restriction(loc_form)
            # TODO: make `res` immutable?
        return res

    @abstract_method
    def get_local(self, cmat):
        r"""
        Abstract method to get the local forms of the generators w.r.t. a given
        curvature matrix `cmat`.
        """
        pass


class ChernAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Chern forms.
    """

    def get_local(self, cmat):
        r"""
        Return the local Chern forms w.r.t. a given curvature matrix.

        .. ALGORITHM::

            The algorithm is based on the Faddeev-LeVerrier algorithm for the
            characteristic polynomial.
        """
        from sage.symbolic.constants import pi
        from sage.libs.pynac.pynac import I

        dom = cmat[0][0]._domain
        rk = len(cmat)
        dim = dom._dim
        ran = min(rk, dim // 2)
        if ran < 1:
            return []  # nothing to compute
        fac = I / (2 * pi)
        res = []
        m = cmat
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / k
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= I / (2 * pi)
            m = [[sum(cmat[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        res.append(-fac * sum(m[i][i] for i in range(rk)) / ran)
        return res


class PontryaginAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Pontryagin forms.
    """

    def get_local(self, cmat):
        r"""
        Return the local Pontryagin forms w.r.t. a given curvature matrix.

        .. ALGORITHM::

            The algorithm is based on the Faddeev-LeVerrier algorithm for the
            characteristic polynomial.
        """
        from sage.symbolic.constants import pi

        dom = cmat[0][0]._domain
        rk = len(cmat)
        dim = dom._dim
        ran = min(rk // 2, dim // 4)
        if ran < 2:
            return []  # nothing to compute
        fac = -1 / (2 * pi) ** 2
        res = []
        m = cmat2 = [[sum(cmat[i][l].wedge(cmat[l][j])
                          for l in range(rk))
                      for j in range(rk)] for i in range(rk)]
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / (2 * k)
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= -1 / (2 * pi) ** 2
            m = [[sum(cmat2[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        res.append(-fac * sum(m[i][i] for i in range(rk)) / (2 * ran))
        return res


class EulerAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Euler forms.
    """

    @cached_method
    def get(self, nab):
        r"""
        Return the global characteristic forms of the generators w.r.t. a given
        connection.
        """
        if not isinstance(nab, LeviCivitaConnection):
            raise TypeError('Euler forms are currently only supported for '
                            'Levi-Civita connections')
        dom = nab._domain
        vbundle = dom.tangent_bundle()
        rk = vbundle._rank
        if not vbundle.has_orientation():
            raise ValueError('Euler forms can only be defined for orientable '
                             'vector bundles')
        if rk % 2 != 0:
            raise ValueError('Euler forms are currently only supported for '
                             'vector bundles with odd rank')
        res = dom.diff_form(rk)
        g = nab._metric
        for frame in dom._get_min_covering(vbundle.orientation()):
            cmatrix = [[nab.curvature_form(i, j, frame)
                        for j in vbundle.irange()]
                       for i in vbundle.irange()]
            res_loc = self.get_local(cmatrix)  # not the local Euler form
            # TODO: multiply cmatrix with metric tensor
            det = g.det(frame)
            if det.is_trivial_zero():
                raise ValueError(f'metric {g} must be non-degenerate')
            sqrt_det = det.abs().sqrt()
            res.set_restriction(res_loc / sqrt_det)  # local Euler form
            # TODO: make `res` immutable?
        return [res]

    def get_local(self, cmat):
        r"""
        Return the normalized Pfaffian w.r.t. a given curvature matrix.

        .. ALGORITHM::

            The algorithm is based on the Bär-Faddeev-LeVerrier algorithm for
            the Pfaffian.
        """
        from sage.symbolic.constants import pi

        rk = len(cmat)
        ran = rk // 2
        m = a = [cmat[i].copy() for i in range(rk)]
        for i in range(0, rk, 2):
            m[i], m[i + 1] = m[i + 1], m[i]  # swap entries
            for k in range(rk):
                m[k][i + 1] = -m[k][i + 1]
        for k in range(1, ran):
            e = -sum(m[i][i] for i in range(rk)) / (2 * k)
            for i in range(rk):
                m[i][i] += e
            m = [[sum(a[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        e = -sum(m[i][i] for i in range(rk)) / (2 * rk)  # Pfaffian mod sign
        e *= (-1 / (2 * pi)) ** ran  # normalize
        return e
