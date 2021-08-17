from sage.algebras.finite_gca import FiniteGCAlgebra
from sage.combinat.free_module import IndexedFreeModuleElement
from sage.misc.fast_methods import Singleton
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.manifolds.differentiable.affine_connection import AffineConnection
from sage.manifolds.differentiable.bundle_connection import BundleConnection

class CharacteristicCohomologyClass_Chern(IndexedFreeModuleElement):
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
        self._mixed_forms = {}  # dict. of mixed forms w.r.t. this class
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
        return f'Characteristic cohomology class {name} over the {vbundle}'

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

                c = ChernAlgorithm().get(nab)
                parent = self.parent()
                grading = parent.print_options()['sorting_key']
                res = [dom.diff_form_module(i).zero() for i in range(dom._dim + 1)]
                for ind, coeff in self:
                    deg = grading(ind)
                    gen_pow = [fast_wedge_power(f, i) for f, i in zip(c, ind)]
                    res[deg] += coeff * reduce(lambda x, y: x.wedge(y), gen_pow)

                # prepare result:
                res = A(res)

                # preparse names (put brackets around)
                vbundle = self.parent()._vbundle
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
                for i in range(dom._dim // 2 + 1):
                    comp_name = name + f'_{i}' + append_name
                    comp_latex_name = latex_name + r'_{' + str(i) + '}'
                    comp_latex_name += append_latex_name
                    res[2*i].set_name(name=comp_name,
                                      latex_name=comp_latex_name)

                # set global names
                res._name = name + append_name
                res._latex_name = latex_name + append_latex_name

                res.set_immutable()  # set result immutable

                # add result to dict
                self._mixed_forms[nab] = res

        return self._mixed_forms[nab]

    representative = get_form

class CharacteristicCohomologyClassRing_Chern(FiniteGCAlgebra):
    r"""

    """
    Element = CharacteristicCohomologyClass_Chern

    def __init__(self, base, vbundle):
        r"""

        """
        self._vbundle = vbundle
        self._domain = vbundle._base_space
        dim = self._domain._dim
        ran = min(vbundle._rank, dim // 2)
        names = tuple(f'c_{i}' for i in range(1, ran + 1))
        degrees = tuple(2*i for i in range(1, ran + 1))
        super().__init__(base=base, names=names, degrees=degrees,
                         max_degree=dim)

    def _element_constructor_(self, x, name=None, latex_name=None):
        r"""

        """
        R = self.base_ring()

        if x in R:
            one_basis = self.one_basis()
            d = {one_basis: R(x)}
        elif isinstance(x, CharacteristicCohomologyClass_Chern):
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
            raise TypeError(f"do not know how to make x (= {x}) an element of self (={self})")

        return self.element_class(self, d, name=name, latex_name=latex_name)

#*****************************************************************************
# ALGORITHMS
#*****************************************************************************

def fast_wedge_power(form, n):
    r"""
    Return the wedge product power of `form` using a square-and-wedge algorithm.
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
    Algorithm class to generate characteristic forms.
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
                res = [dom.diff_form(loc_form.degree()) for loc_form in res_loc]
            for form, loc_form in zip(res, res_loc):
                form.set_restriction(loc_form)
            # TODO: make res immutable?
        return res

    @abstract_method
    def get_local(self, cmat):
        pass

class ChernAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Chern forms.
    """
    def get_local(self, cmat):
        r"""
        Return the local Chern forms for a given curvature matrix.

        .. ALGORITHM::

            The algorithm is based on the Faddeev-LeVerrier algorithm for the
            characteristic polynomial.
        """
        from sage.symbolic.constants import pi
        from sage.libs.pynac.pynac import I

        dom = cmat[0][0]._domain
        rk = len(cmat)
        dim = dom._dim
        ran = min(rk, dim//2)
        fac = I / (2*pi)
        res = []
        m = cmat
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / k
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= I / (2*pi)
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
        Return the local Pontryagin forms for a given curvature matrix.

        .. ALGORITHM::

            The algorithm is based on the Faddeev-LeVerrier algorithm for the
            characteristic polynomial.
        """
        from sage.symbolic.constants import pi

        dom = cmat[0][0]._domain
        rk = len(cmat)
        dim = dom._dim
        ran = min(rk//2, dim//4)
        fac = -1 / (2*pi)**2
        res = []
        m = cmat2 = [[sum(cmat[i][l].wedge(cmat[l][j])
                          for l in range(rk))
                      for j in range(rk)] for i in range(rk)]
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / (2*k)
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= -1 / (2*pi)**2
            m = [[sum(cmat2[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        res.append(-fac * sum(m[i][i] for i in range(rk)) / (2*ran))
        return res

class EulerAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Euler forms.
    """
    def get_local(self, cmat):
        r"""
        Return the local Euler form for a given curvature matrix.

        .. ALGORITHM::

            The algorithm is based on the Bär-Faddeev-LeVerrier algorithm for
            the Pfaffian.
        """
        from sage.symbolic.constants import pi

        rk = len(cmat)
        ran = rk // 2
        m = a = [cmat[i].copy() for i in range(rk)]
        for i in range(0, rk, 2):
            m[i], m[i+1] = m[i+1], m[i]  # swap entries
            for k in range(rk):
                m[k][i+1] = -m[k][i+1]
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / (2*k)
            for i in range(rk):
                m[i][i] += c
            m = [[sum(a[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        c = -sum(m[i][i] for i in range(rk)) / (2*rk)  # Pfaffian mod sign
        c *= (-1/(2*pi))**rk  # normalize
        return [c]
