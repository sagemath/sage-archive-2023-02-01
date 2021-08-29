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

                parent = self.parent()
                algorithm = parent._algorithm

                grading = parent.print_options()['sorting_key']
                res = [dom.diff_form_module(i).zero()
                       for i in range(dom._dim + 1)]
                for ind, c in self:
                    deg = grading(ind)
                    gen_pow = [algorithm.get_gen_pow(nab, i, ind[i])
                               for i in range(len(ind))]
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

    def representative(self, nab=None):
        r"""
        Return any representative of ``self``.

        INPUT:

        - ``nab`` -- (default: ``None``) if stated, return the representative
          w.r.t. to the connection ``nab``; otherwise an arbitrary already
          computed representative will be chosen

        """
        if nab is None:
            if not self._mixed_forms:
                raise AttributeError('cannot pick a representative')
            return next(iter(self._mixed_forms.values()))
        return self.get_form(nab)


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
            self._algorithm = ChernAlgorithm()
        elif vbundle._field_type == 'real':
            ran = min(rk // 2, dim // 4)
            names = [f'p_{i}' for i in range(1, ran + 1)]
            degrees = [4 * i for i in range(1, ran + 1)]
            self._algorithm = PontryaginAlgorithm()
            if vbundle.has_orientation():
                # add Euler class generator
                # Euler should be first entry; see `PontryaginEulerAlgorithm`
                names = ['e'] + names
                degrees = [rk] + degrees
                self._algorithm = PontryaginEulerAlgorithm()
                # TODO: add relation e^2=p_k for dim=2*k
        else:
            raise TypeError(f'Characteristic cohomology classes not supported '
                            f'for vector bundles with '
                            f'field type {vbundle._field_type}')

        if not names or not degrees:
            raise ValueError(f'cannot find any generators')

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

    def _repr_(self):
        r"""

        """
        vbundle = self._vbundle
        repr = f'Algebra of characteristic cohomology classes of the {vbundle}'
        return repr

# *****************************************************************************
# ALGORITHMS
# *****************************************************************************

def multiplicative_sequence(q):
    r"""
    Turn the polynomial ``q`` into its multiplicative sequence.

    OUTPUT:

    - A symmetric polynomial representing the multiplicative sequence.
    """
    from sage.combinat.sf.sf import SymmetricFunctions
    from sage.combinat.partition import Partitions
    from sage.misc.misc_c import prod

    R = q.parent().base_ring()
    Sym = SymmetricFunctions(R)
    m = Sym.m()

    # Get the multiplicative sequence in the monomial basis:
    mon_pol = m._from_dict({p: prod(q[i] for i in p)
                            for k in range(len(q.degree()))
                            for p in Partitions(k)})
    return Sym.e()(mon_pol)

def additive_sequence(q, rk):
    r"""
    Turn the polynomial ``q`` into its additive sequence.

    INPUT:

    - ``q`` -- polynomial to turn into its additive sequence.
    - ``rk`` -- rank of the underlying vector bundle

    OUTPUT:

    - A symmetric polynomial representing the additive sequence.
    """
    from sage.combinat.sf.sf import SymmetricFunctions
    from sage.combinat.partition import Partitions

    R = q.parent().base_ring()
    Sym = SymmetricFunctions(R)
    m = Sym.m()

    # Express the additive sequence in the monomial basis, the 0th
    # order term must be treated separately; here comes ``rk`` into play:
    m_dict = {Partitions(0)([]): rk * q[0]}
    m_dict.update({Partitions(k)([k]): q[k] for k in range(1, q.degree())})
    mon_pol = m._from_dict(m_dict)
    return Sym.e()(mon_pol)


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

        OUTPUT:

        - a list containing the generator's global characteristic forms

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
            cmat = [[nab.curvature_form(i, j, frame) for j in vbundle.irange()]
                    for i in vbundle.irange()]
            res_loc = self.get_local(cmat)
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
        curvature matrix ``cmat``.

        OUTPUT:

        - a list containing the generator's local characteristic forms

        """
        pass

    @cached_method
    def get_gen_pow(self, nab, i, n):
        r"""
        Return the `n`-th power of the `i`-th generator's characteristic form
        w.r.t ``nab``.
        """
        if n == 0:
            return nab._domain._one_scalar_field  # no computation necessary
        return fast_wedge_power(self.get(nab)[i], n)


class ChernAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Chern forms.
    """

    def get_local(self, cmat):
        r"""
        Return the local Chern forms w.r.t. a given curvature matrix.

        OUTPUT:

        - a list containing the local characteristic Chern forms

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

        OUTPUT:

        - a list containing the local characteristic Pontryagin forms

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

        OUTPUT:

        - a list containing the global characteristic Euler form

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
            # (G_s * Ω_s)_ij = g(R(.,.)s_i, s_j)
            gcmat = [[sum(g[[frame, i, j]] * nab.curvature_form(j, k, frame)
                          for j in vbundle.irange())
                      for k in vbundle.irange()] for i in vbundle.irange()]
            [res_loc] = self.get_local(gcmat)  # Pf(G_s * Ω_s) mod const.
            # e = 1 / sqrt(|det(G_s)|) * Pf(G_s * Ω_s) mod const.
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

        The normalization is given by the factor
        `\left(\frac{1}{2 \pi}\right)^{\frac{k}{2}}`, where `k` is the
        dimension of the curvature matrix.

        OUTPUT:

        - a list containing the normalized Pfaffian of a given curvature form

        .. NOTE::

            The result is the local Euler form if ``cmat`` is given w.r.t. an
            orthonormal oriented frame.

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
        e = -sum(m[i][i] for i in range(rk)) / (2 * ran)  # Pfaffian mod sign
        e *= (-1 / (2 * pi)) ** ran  # normalize
        return [e]


class PontryaginEulerAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Euler and Pontryagin forms.
    """

    @cached_method
    def get(self, nab):
        r"""
        Return the global characteristic forms of the generators w.r.t. a given
        connection.

        OUTPUT:

        - a list containing the global Euler form in the first entry, and the
          global Pontryagin forms in the remaining entries.

        """
        return EulerAlgorithm().get(nab) + PontryaginAlgorithm().get(nab)

    def get_local(self, cmat):
        r"""
        Return the local Euler and Pontryagin forms w.r.t. a given curvature
        matrix.

        OUTPUT:

        - a list containing the local Euler form in the first entry, and the
          local Pontryagin forms in the remaining entries.

        """
        res = EulerAlgorithm().get_local(cmat)  # first entry is Euler class
        res += PontryaginAlgorithm().get_local(cmat)  # rest Pontryagin
        return res

    @cached_method
    def get_gen_pow(self, nab, i, n):
        r"""
        Return the `n`-th power of the `i`-th generator w.r.t ``nab``.
        """
        if n == 0:
            return nab._domain._one_scalar_field  # no computation necessary
        if i == 0:
            return fast_wedge_power(EulerAlgorithm().get(nab)[0], n)
        return fast_wedge_power(PontryaginAlgorithm().get(nab)[i-1], n)
