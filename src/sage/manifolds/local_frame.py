r"""
Local Frames

AUTHORS:

- Michael Jung (2019): initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_basis import (FreeModuleBasis,
                                                   FreeModuleCoBasis)
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule

class LocalCoFrame(FreeModuleCoBasis):
    r"""
    Local coframe on a vector bundle.



    INPUT:

    - ``frame`` -- the vector frame dual to the coframe

    EXAMPLES:



    """
    def __init__(self, frame, symbol, latex_symbol=None, indices=None,
                 latex_indices=None):
        r"""
        Construct a local coframe, dual to a given local frame.

        TESTS::



        """
        self._domain = frame._domain
        self._manifold = self._domain.manifold()
        FreeModuleCoBasis.__init__(self, frame, symbol,
                                   latex_symbol=latex_symbol, indices=indices,
                                   latex_indices=latex_indices)
        # The coframe is added to the domain's set of coframes, as well as to
        # all the superdomains' sets of coframes
        for sd in self._domain._supersets:
            sd._coframes.append(self)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

        """
        desc = "Local coframe " + self._name + " "
        desc += "on {}".format(self._basis.module())
        return desc

    def at(self, point):
        r"""

        """
        return self._basis.at(point).dual_basis()

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='up',
                 include_domain=True):
        r"""

        EXAMPLES::

        """
        pass

#******************************************************************************

class LocalFrame(FreeModuleBasis):
    r"""

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = LocalCoFrame

    @staticmethod
    def __classcall_private__(cls, section_module, symbol,
                              latex_symbol=None, from_frame=None, indices=None,
                              latex_indices=None, symbol_dual=None,
                              latex_symbol_dual=None):
        """
        Transform input lists into tuples for the unique representation of
        LocalFrame.

        TESTS::



        """
        if isinstance(symbol, list):
            symbol = tuple(symbol)
        if isinstance(latex_symbol, list):
            latex_symbol = tuple(latex_symbol)
        if isinstance(indices, list):
            indices = tuple(indices)
        if isinstance(latex_indices, list):
            latex_indices = tuple(latex_indices)
        if isinstance(symbol_dual, list):
            symbol_dual = tuple(symbol_dual)
        if isinstance(latex_symbol_dual, list):
            latex_symbol_dual = tuple(latex_symbol_dual)
        return super(LocalFrame, cls).__classcall__(cls, section_module,
                                        symbol, latex_symbol=latex_symbol,
                                        from_frame=from_frame, indices=indices,
                                        latex_indices=latex_indices,
                                        symbol_dual=symbol_dual,
                                        latex_symbol_dual=latex_symbol_dual)

    def __init__(self, section_module, symbol, latex_symbol=None,
                 indices=None, latex_indices=None):
        pass

    def _repr_(self):
        pass

    def _new_instance(self, symbol, latex_symbol=None, indices=None,
                      latex_indices=None, symbol_dual=None,
                      latex_symbol_dual=None):
        pass

    def domain(self):
        pass

    def coframe(self):
        r"""

        """
        return self._coframe

    def new_frame(self, change_of_frame, symbol, latex_symbol=None,
                  indices=None, latex_indices=None, symbol_dual=None,
                  latex_symbol_dual=None):
        pass

    def restrict(self, subdomain):
        pass

    @cached_method
    def structure_coeff(self):
        pass

    def at(self, point):
        pass

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='up',
                 include_domain=True):
        pass



#******************************************************************************

class TrivializationCoFrame(LocalCoFrame):
    def __init__(self, trivialization_frame, symbol, latex_symbol=None,
                 indices=None, latex_indices=None):
        r"""
        Construct a local coframe from a local trivialization.

        TESTS::


        """
        if not isinstance(coord_frame, CoordFrame):
            raise TypeError("the first argument must be a local trivialization "
                            "frame")
        LocalCoFrame.__init__(self, trivialization_frame, symbol,
                              latex_symbol=latex_symbol, indices=indices,
                              latex_indices=latex_indices)
        self._trivialization = trivialization_frame._trivialization

#******************************************************************************

class TrivializationFrame(LocalFrame):
    r"""

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = TrivializationCoFrame

    def __init__(self, trivialization):
        r"""

        """
        from sage.misc.latex import latex
        from sage.manifolds.differentiable.chart import DiffChart
        if not isinstance(chart, DiffChart):
            raise TypeError("the first argument must be a chart")
        self._trivialization = trivialization
        # TODO: Adept LaTeX
        coords = chart[:] # list of all coordinates
        symbol = tuple("d/d" + str(x) for x in coords)
        latex_symbol = tuple(r"\frac{\partial}{\partial" + latex(x) + "}"
                             for x in coords)
        symbol_dual = tuple("d" + str(x) for x in coords)
        latex_symbol_dual = tuple(r"\mathrm{d}" + latex(x) for x in coords)
        VectorFrame.__init__(self,
                             chart._domain.vector_field_module(force_free=True),
                             symbol=symbol, latex_symbol=latex_symbol,
                             symbol_dual=symbol_dual,
                             latex_symbol_dual=latex_symbol_dual)

    def trivialization(self):
        r"""
        Return the underlying trivialization of ``self``.

        OUTPUT:

        -

        EXAMPLES::



        """
        return self._trivialization