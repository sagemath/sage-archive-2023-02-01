r"""
Elements of Hecke modular forms spaces

AUTHORS:

- Jonas Jermann (2013): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from graded_ring_element import FormsRingElement


class FormsElement(FormsRingElement):
    """
    (Hecke) modular forms.
    """

    def __init__(self, parent, rat):
        r"""
        An element of a space of (Hecke) modular forms.

        INPUT:

        - ``parent``     -- a modular form space
        - ``rat``        -- a rational function which corresponds to a
                            modular form in the modular form space

        OUTPUT:

        A (Hecke) modular form element corresponding to the given rational function
        with the given parent space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: MF = ModularForms(n=5, k=20/3, ep=1)
            sage: MF.default_prec(3)
            sage: el = MF(x^5*d-y^2*d)
            sage: el
            q - 9/(200*d)*q^2 + O(q^3)
            sage: el.rat()
            x^5*d - y^2*d
            sage: el.parent()
            ModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: el.rat().parent()
            Fraction Field of Multivariate Polynomial Ring in x, y, z, d over Integer Ring

            sage: subspace = MF.subspace([MF.gen(1)])
            sage: ss_el = subspace(x^5*d-y^2*d)
            sage: ss_el == el
            True
            sage: ss_el.parent()
            Subspace of dimension 1 of ModularForms(n=5, k=20/3, ep=1) over Integer Ring
        """

        super(FormsElement, self).__init__(parent, rat)

        if self.AT(["quasi"])>=self._analytic_type:
            pass
        elif not (\
            self.is_homogeneous() and\
            self._weight  == parent.weight() and\
            self._ep      == parent.ep() ):
                raise ValueError("{} does not correspond to an element of {}.".format(rat, parent))

        from subspace import SubSpaceForms
        if isinstance(parent, SubSpaceForms) and (parent._module is not None):
            try:
                self.coordinate_vector()
            except TypeError:
                raise ValueError("{} does not correspond to an element of {}.".format(rat, parent))

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: QuasiModularForms(n=5, k=10, ep=-1)(x^3*z^3-y^3)
            21/(20*d)*q - 4977/(16000*d^2)*q^2 + 297829/(12800000*d^3)*q^3 + 27209679/(20480000000*d^4)*q^4 + O(q^5)
        """

        return self._qexp_repr()

    # This function is just listed here to emphasize the choice used
    # for the latex representation of ``self``
    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: latex(QuasiModularForms(n=5, k=10, ep=-1)(x^3*z^3-y^3))
            f_{\rho}^{3} E_{2}^{3} -  f_{i}^{3}
        """

        return super(FormsElement, self)._latex_()

    def coordinate_vector(self):
        r"""
        Return the coordinate vector of ``self`` with
        respect to ``self.parent().gens()``.

        Note: This uses the corresponding function of the
        parent. If the parent has not defined a coordinate
        vector function or a module for coordinate vectors
        then an exception is raised by the parent
        (default implementation).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.gen(0).coordinate_vector().parent()
            Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.gen(0).coordinate_vector()
            (1, 0, 0)
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(2)])
            sage: subspace.gen(0).coordinate_vector().parent()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: subspace.gen(0).coordinate_vector()
            (1, 0)
            sage: subspace.gen(0).coordinate_vector() == subspace.coordinate_vector(subspace.gen(0))
            True
        """

        return self.parent().coordinate_vector(self)

    def ambient_coordinate_vector(self):
        r"""
        Return the coordinate vector of ``self`` with
        respect to ``self.parent().ambient_space().gens()``.

        The returned coordinate vector is an element
        of ``self.parent().module()``.        
        
        Mote: This uses the corresponding function of the
        parent. If the parent has not defined a coordinate
        vector function or an ambient module for
        coordinate vectors then an exception is raised
        by the parent (default implementation).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.gen(0).ambient_coordinate_vector().parent()
            Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.gen(0).ambient_coordinate_vector()
            (1, 0, 0)
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(2)])
            sage: subspace.gen(0).ambient_coordinate_vector().parent()
            Vector space of degree 3 and dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            Basis matrix:
            [1 0 0]
            [0 0 1]
            sage: subspace.gen(0).ambient_coordinate_vector()
            (1, 0, 0)
            sage: subspace.gen(0).ambient_coordinate_vector() == subspace.ambient_coordinate_vector(subspace.gen(0))
            True
        """

        return self.parent().ambient_coordinate_vector(self)
