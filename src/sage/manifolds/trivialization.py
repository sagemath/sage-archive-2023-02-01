r"""
Trivializations

AUTHORS:

- Michael Jung (2019) : initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

class Trivialization(UniqueRepresentation, SageObject):
    r"""

    """
    def __init__(self, vbundle, domain):
        r"""
        Construct a local trivialization.

        """
        self._base_space = domain.manifold()
        self._vbundle = vbundle
        self._bdl_rank = vbundle.rank()
        self._base_field = vbundle.base_field()
        self._domain = domain
        self._sindex = self._base_space.start_index()
        # Add this trivialization to the atlas of the vector bundle:
        vbundle._atlas.append(self)

    def _repr_(self):
        r"""
        Return the string representation of self.


        """
        description = "Trivialization "
        description += "({}|_{} -> {})".format(self._vbundle._name,
                                               self._domain._name,
                                               self._domain._name)
        return description

    def _latex_(self):
        r"""
        Return the LaTeX representation of self.


        """
        from sage.misc.latex import latex
        latex = r'{} |_{{{}}} \to {} \times {}^{}'.format(self._vbundle._latex_name,
                                        latex(self._domain), latex(self._domain),
                                        latex(self._base_field), self._bdl_rank)
        return latex

    def domain(self):
        r"""
        Return the open subset on which the trivialization is defined.


        """
        return self._domain

    def base_space(self):
        r"""
        Return the manifold on which the trivialization is defined.


        """
        return self._base_space

    def transition_map(self, other, trafo_matrix):
        r"""

        """
        return TransitionMap(self, other, trafo_matrix)

    def vector_bundle(self):
        r"""

        """
        return self._vbundle


# *****************************************************************************

class TransitionMap(SageObject):
    r"""

    """
    def __init__(self, triv1, triv2, matrix):
        r"""

        """
        bs1 = triv1.base_space()
        bs2 = triv2.base_space()
        if bs1 is not bs2:
            raise ValueError("base spaces must coincide")
        self._base_space = bs1

        vb1 = triv1.vector_bundle()
        vb2 = triv2.vector_bundle()
        if vb1 is not vb2:
            raise ValueError("vector bundles must coincide")
        self._vbundle = vb1
        self._bdl_rank = self._vbundle.rank()

        mfd_field = self._base_space.base_field()
        vb_field = self._vbundle.base_field()
        if not mfd_field.is_subring(vb_field):
            raise ValueError("for concrete implementation, manifold's base "
                             "field must be a subfield of the vector bundle's "
                             "base field")

        dom1 = triv1.domain()
        dom2 = triv2.domain()
        dom = dom1.intersection(dom2)
        scal_field_alg = dom.scalar_field_algebra()
        from sage.matrix.matrix_space import MatrixSpace
        matrix_space = MatrixSpace(scal_field_alg, self._bdl_rank)
        self._matrix = matrix_space(matrix)
        self._domain = dom
        self._triv1 = triv1
        self._triv2 = triv2
        self._inverse = None
        self._vbundle._transitions[(triv1, triv2)] = self

    def _repr_(self):
        r"""

        """
        dom1 = self._triv1.domain()
        dom2 = self._triv2.domain()
        dom = self._domain
        description = "({}|_{}, g_{}{})".format(self._vbundle._latex_name,
                                                dom._name, dom1._name,
                                                dom2._name)
        return description

    def _latex_(self):
        r"""

        """
        dom1 = self._triv1.domain()
        dom2 = self._triv2.domain()
        dom = self._domain
        from sage.misc.latex import latex
        latex = r'({}, g_{{{}{}}})'.format(latex(dom), latex(dom1), latex(dom2))
        return latex
    
    def transformation_matrix(self, chart=None):
        r"""

        """
        if chart is None:
            return self._matrix
        else:
            m_list = []
            for row in self._matrix:
                for entry in row:
                    m_list.append(entry.expr(chart))
        from sage.matrix.matrix_space import MatrixSpace
        parent = m_list[0].parent()
        matrix_space = MatrixSpace(parent, self._bdl_rank)
        return matrix_space(m_list)

    def inverse(self):
        r"""

        """
        self._inverse = type(self)(self._triv2, self._triv1, ~self._matrix)
        return self._inverse