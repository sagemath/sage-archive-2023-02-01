# -*- coding: utf-8 -*-
r"""
KontInfo

This module contains the class :class:`KnotInfoBase` which is derived from :class:`Enum`
and provides all knots and links listed in the databases at https://knotinfo.math.indiana.edu/
and https://linkinfo.sitehost.iu.edu as its items.

Be aware that there are a couple of conventions used differently on KnotInfo as in Sage, especially
concerning the selection of the symmetry version of the link. In our transitions to Sage objects
these are translated in order to avoid confusion about exchanged mirror versions.

Briefly, these differences are:

   - ``pd_notation`` --        KnotInfo: counter clockwise Sage: clockwise, see note in
     :meth:`link`

   - ``homfly_polynomial`` --  KnotInfo: ``v``  Sage: `1/a`, see note in :meth:`homfly_polynomial`.

   - ``braid_notation``    --  This is used accordingly: The crossing of the braid generators are positive
     in both systems. Here it is listed because there could arise confusion from the source where they are
     taken from. There, the braid generators are assumed to have a negative crossing
     (see definition 3  of Gittings, T., "Minimum Braids: A Complete Invariant of Knots and Links
     https://arxiv.org/abs/math/0401051).

REFERENCES:

- https://knotinfo.math.indiana.edu/
- https://linkinfo.sitehost.iu.edu/



AUTHORS:

- Sebastian Oehms August 2020: initial version
"""


##############################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################



from enum import Enum
from sage.misc.cachefunc import cached_method
from sage.misc.sage_eval import sage_eval
from sage.groups.braid import BraidGroup
from sage.knots.knot import Knots
from sage.databases.knotinfo_db import KnotInfoColumnTypes, KnotInfoColumns, db




def is_knotinfo_available():
    r"""
    Return wether the KnotInfo databases are installed or not.

    EXAMPLES::

        sage: from sage.knots.knotinfo import is_knotinfo_available
        sage: is_knotinfo_available()     # optional - database_knotinfo
        True
    """
    return db.is_available()


def eval_knotinfo(string, locals={}, to_tuple=True):
    r"""
    Preparse a string from the KnotInfo database and evaluate it by ``sage_eval``.

    INPUT:

    - ``string``  -- string that gives a value of some database entry
    - ``locals``      -- dictionary of locals passed to ``sage_eval``

    EXAMPLES::

        sage: from sage.knots.knotinfo import KnotInfo, eval_knotinfo
        sage: L = KnotInfo.L4a1_0
        sage: L.braid_notation(original=True)
        '{4, {1, -2, 3, -2, -1, -2, -3, -2}}'
        sage: eval_knotinfo(_)
        (4, (1, -2, 3, -2, -1, -2, -3, -2))
    """
    if to_tuple:
        new_string = string.replace('{', '(')
        new_string = new_string.replace('}', ')')
    else:
        new_string = string.replace('{', '[')
        new_string = new_string.replace('}', ']')
    new_string = new_string.replace(';', ',')
    return sage_eval(new_string, locals=locals)



class KnotInfoBase(Enum):
    r"""
    Enum class to select the knots and links provided by http://www.indiana.edu/~knotinfo

    EXAMPLES::
        sage: from sage.knots.knotinfo import KnotInfo
        sage: [knot.name for knot in KnotInfo if knot.crossing_number() < 5]
        ['K0_1', 'K3_1', 'K4_1', 'L2a1_0', 'L2a1_1', 'L4a1_0', 'L4a1_1']
    """
    @property
    def items(self):
        r"""
        Return an Enum class to select a column item of the KnotInfo database.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: it = L.items
            sage: [i.name for i in it if i.name.endswith('notation')]   # optional - database_knotinfo
            ['dt_notation',
             'conway_notation',
             'two_bridge_notation',
             'gauss_notation',
             'enhanced_gauss_notation',
             'pd_notation',
             'braid_notation',
             'positive_braid_notation',
             'positive_pd_notation',
             'strongly_quasipositive_braid_notation',
             'quasipositive_braid_notation',
             'arc_notation']
        """
        return db.columns()

    @cached_method
    def __getitem__(self, item):
        """
        sage: from sage.knots.knotinfo import KnotInfo
        sage: L = KnotInfo.L4a1_0
        sage: L[L.items.alternating]
        'Y'
        sage: L[L.items.arc_notation]
        '{{6, 4}, {3, 5}, {4, 2}, {1, 3}, {2, 6}, {5, 1}}'
        sage: L[L.items.braid_notation]
        '{4, {1, -2, 3, -2, -1, -2, -3, -2}}'
        sage: L[0]
        Traceback (most recent call last):
        ...
        KeyError: "Item must be an instance of <enum 'KnotInfoColumns'>"
        """
        if not isinstance(item, KnotInfoColumns):
            raise KeyError('Item must be an instance of %s' %(KnotInfoColumns))
        if item.column_type() == KnotInfoColumnTypes.OnlyLinks and self.is_knot():
            raise KeyError('Item not available for knots' %(KnotInfoColumns))
        if item.column_type() == KnotInfoColumnTypes.OnlyKnots and not self.is_knot():
            raise KeyError('Item not available for links' %(KnotInfoColumns))

        l = db.read(item)
        offset = 0
        if item.column_type() == KnotInfoColumnTypes.OnlyLinks:
            offset = self._offset_knots()

        return l[self.value[0]-offset]

    def _offset_knots(self):
        r"""
        Return the list index of the first proper link in a conbined
        list containing knots and proper links together which is the
        case for columns used for KnotInfo and LinkInfo in common.
        This index is exactly the total number of knots recorded
        in KnotInfo.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L._offset_knots()          # optional - database_knotinfo
            2978
        """
        return db.read_num_knots()

    @cached_method
    def _braid_group(self):
        r"""
        Return the braid group corresponding to the braid index
        of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L._braid_group()
            Braid group on 4 strands
        """
        n = self.braid_index()
        if n == 1:
            return BraidGroup(2)
        else:
            return BraidGroup(n)


    @cached_method
    def _homfly_pol_ring(self, var1, var2):
        r"""
        Return the parent Laurent polynomial ring for the Homfly-PT
        polynomial according to Sage's internal one.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_1
            sage: L._homfly_pol_ring('u', 'v')
            Multivariate Laurent Polynomial Ring in u, v over Integer Ring
        """
        K3_1 = Knots().from_table(3,1)
        return K3_1.homfly_polynomial(var1=var1, var2=var2).parent()

    @cached_method
    def pd_notation(self, original=False):
        r"""
        Return the value of column ``pd_notation`` for this
        link as a Python list of Python lists.

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT::

        Python list of python lists each entry of the outer list
        representing a crossing.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.pd_notation()
            [[6, 1, 7, 2], [8, 3, 5, 4], [2, 5, 3, 6], [4, 7, 1, 8]]
            sage: L.pd_notation(original=True)
            '{{6, 1, 7, 2}, {8, 3, 5, 4}, {2, 5, 3, 6}, {4, 7, 1, 8}}'
            sage: K = KnotInfo.K4_1
            sage: K.pd_notation()
            [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]
        """
        if self.is_knot():
            pd_notation = self[self.items.pd_notation]
        else:
            pd_notation = self[self.items.pd_notation_vector]

        if original:
            return pd_notation

        if not pd_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(pd_notation, to_tuple=False)

    @cached_method
    def dt_notation(self, original=False):
        r"""
        Return the value of column ``dt_notation`` for this
        link as a Python list of Python lists.

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT::

        Python list of python lists each entry of the outer list
        representing a crossing.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.dt_notation()
            [[6, 8], [2, 4]]
            sage: L.dt_notation(original=True)
            '[{6, 8}, {2, 4}]'
            sage: L = KnotInfo.L4a1_0
            sage: K = KnotInfo.K4_1
            sage: K.dt_notation()
            [4, 6, 8, 2]
        """
        if self.is_knot():
            dt_notation = self[self.items.dt_notation]
        else:
            dt_notation = self[self.items.dt_code]

        if original:
            return dt_notation

        if not dt_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(dt_notation, to_tuple=False)

    @cached_method
    def gauss_notation(self, original=False):
        r"""
        Return the value of column ``gauss_notation`` for this
        link as a Python list of Python lists.

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT::

        Python list of

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.gauss_notation()
            [[1, -3, 2, -4], [3, -1, 4, -2]]
            sage: L.gauss_notation(original=True)
            '{{1, -3, 2, -4}, {3, -1, 4, -2}}'
        """
        gauss_notation = self[self.items.gauss_notation]
        if original:
            return gauss_notation

        if not gauss_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(gauss_notation, to_tuple=False)

    @cached_method
    def braid_notation(self, original=False):
        r"""
        Return the value of column ``braid_notation`` for this
        link as a Python tuple (Tietze form).

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT::

        Python tuple representing the braid whose closure is ``self``
        in Tietze form.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.braid_notation()
            (1, -2, 3, -2, -1, -2, -3, -2)
            sage: L.braid_notation(original=True)
            '{4, {1, -2, 3, -2, -1, -2, -3, -2}}'
        """
        braid_notation = self[self.items.braid_notation]
        if original:
            return braid_notation

        if not braid_notation:
            # don't forget the unknot
            return (1, -1)

        braid_notation = eval_knotinfo(braid_notation)
        if type(braid_notation) is list:
            # in some cases there are a pair of braid representations
            # in the database. If this is the case we select the
            # corresponding to the braid index.
            if type(braid_notation[0]) is tuple:
                i = self.braid_index()
                for b in braid_notation:
                    if -i < min(b) and max(b) < i:
                        braid_notation = b
                        break

        if not self.is_knot():
            # in linkinfo the braid_notation includes the braid_index as first item of a pair
            braid_notation = braid_notation[1]
        return braid_notation

    @cached_method
    def braid_index(self):
        r"""
        Return the value of column ``braid_index`` for this
        link as a Python int.

        OUTPUT::

        Python int giving the minimum of strands needed to
        represent ``self`` as closure of a braid.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.braid_index()
            4
        """
        if self.is_knot():
            return int(self[self.items.braid_index])
        else:
            braid_notation = self[self.items.braid_notation]
            braid_notation = eval_knotinfo(braid_notation)
            return int(braid_notation[0])

    @cached_method
    def braid_length(self):
        r"""
        Return the value of column ``braid_length`` for this
        link as a Python int.

        OUTPUT::

        Python int giving the minimum length of a braid word
        needed to represent ``self`` as closure of a braid.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.braid_length()
            3
        """
        return int(self[self.items.braid_length])

    @cached_method
    def braid(self):
        r"""
        Return the braid notation of self as an instance of :class:`Braid`.


        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.braid()
            s^3
            sage: K.braid_notation()
            (1, 1, 1)
        """
        return self._braid_group()(self.braid_notation())

    def num_components(self):
        r"""
        Return the number of compoents of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L6a1_0.num_components()
            2
        """
        return self.value[1]

    @cached_method
    def crossing_number(self):
        r"""
        Return the minimal number of crossings.

        .. NOTE::

           In contrast to the number of crossings displayed for instances
           of :class:`Link` this number is the minimum over all possible
           diagrams of the link. The number of crossings displayed in
           the representation string of :class:`Link` referes to the
           special representation which could be larger.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L4a1_0.crossing_number()
            4
            sage: KnotInfo.K3_1.crossing_number()
            3
            sage: Link(KnotInfo.L4a1_0.braid())
            Link with 2 components represented by 8 crossings
        """
        return int(self[self.items.crossing_number])

    def is_knot(self):
        r"""
        Return wether ``self`` is a knot or a proper link.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L7a1_0.is_knot()      # optional - database_knotinfo
            False
            sage: KnotInfo.K6_3.is_knot()        # optional - database_knotinfo
            True
        """
        return self.num_components() == 1

    @cached_method
    def symmetry_type(self):
        r"""
        Return the symmetry type of ``self``.

        From the KnotInfo description page:

            If a knot is viewed as the oriented diffeomorphism
            class of an oriented pair, `K = (S_3, S_1), with `S_i`
            diffeomorphic to `S^i`, there are four oriented knots
            associated to any particular knot `K`. In addition to
            `K` itself, there is the reverse, `K^r = (S_3, -S_1)`,
            the concordance inverse, `-K = (-S_3, -S_1)`, and the
            mirror image, `K^m = (-S3, S1)`. A knot is called
            reversible if `K = K^r`, negative amphicheiral if
            `K = -K`, and positive amphicheiral if `K = K^m`.

            A knot possessing any two of these types of symmetry
            has all three. Thus, in the table, a knot is called
            reversible if that is the only type of symmetry it has,
            and likewise for negative amphicheiral. If it has none
            of these types of symmetry it is called chiral, and if
            it has all three it is called fully amphicheiral.

            For prime knots with fewer than 12 crossings, all
            amphicheiral knots are negative amphicheiral.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: [(L.name, L.symmetry_type()) for L in KnotInfo if L.is_knot() and L.crossing_number() < 6]
            [('K0_1', 'fully amphicheiral'),
            ('K3_1', 'reversible'),
            ('K4_1', 'fully amphicheiral'),
            ('K5_1', 'reversible'),
            ('K5_2', 'reversible')]
        """
        if not self.is_knot():
            raise NotImplementedError('This is only available for knots')
        if not self[self.items.symmetry_type] and self.crossing_number() == 0:
            return 'fully amphicheiral'
        return self[self.items.symmetry_type]


    def is_reversible(self):
        r"""
        Return wether ``self`` is reversible.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K6_3.is_reversible()
            True
        """
        if self.symmetry_type() == 'reversible':
            return True
        if self.symmetry_type() == 'fully amphicheiral':
            return True
        return False

    def is_amphicheiral(self, positive=False):
        r"""
        Return wether ``self`` is amphicheiral.

        INPUT:

        - ``positive`` -- Boolean (default False) wether to check
          if ``self`` is positive or negative amphicheiral (see
          doctest of :meth:`symmetry_type`)

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K12a_427                 # optional - database_knotinfo
            sage: K.is_amphicheiral()                   # optional - database_knotinfo
            False
            sage: K.is_amphicheiral(positive=True)      # optional - database_knotinfo
            True
        """
        if positive:
            if self.symmetry_type() == 'positive amphicheiral':
                return True
        else:
            if self.symmetry_type() == 'negative amphicheiral':
                return True
        if self.symmetry_type() == 'fully amphicheiral':
            return True
        return False


    @cached_method
    def homfly_polynomial(self, var1='L', var2='M', original=False):
        r"""
        Return the value of column ``homfly_polynomial`` for this
        knot or link (in this case the column ``homflypt_polynomial``
        is used) as an instance of the element class according to
        the output of :meth:`homfly_polynomial` of :class:`Link`.

        INPUT:

        - ``var1`` -- (default: ``'L'``) the first variable
        - ``var2`` -- (default: ``'M'``) the second variable
        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT::

        A Laurent polynomial over the integers, more precisely an instance of
        :class:`sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.

        .. NOTE::

            The skein-relation for the Homfly-PT polynomial given on KnotInfo
            differs from the ones used in Sage:

            KnotInfo: P(O) = 1,   ~v P(L+) -  v P(L-) = z P(L0)

            (see: https://knotinfo.math.indiana.edu/descriptions/jones_homfly_kauffman_description/polynomial_defn.html)

            Using Sage's Homfy-PT polynomials with ``normalization='az'``
            the corresponding skein-relation is (see :meth:`homfly_polynomial`
            of :class:`Link`):

            Sage:     P(O) = 1,    a P(L+) - ~a P(L-) = z P(L0)

            Thus, the Homfly-PT polynomial of KnotInfo compares to the one of Sage
            by replacing ``v`` by ``~a``. To keep them comparable this translation is
            performed, as well.


        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K3_1 = KnotInfo.K3_1
            sage: PK3_1 = K3_1.homfly_polynomial(); PK3_1
            L^-2*M^2 + 2*L^-2 - L^-4
            sage: PK3_1 == K3_1.link().homfly_polynomial(normalization='az')
            True
            sage: L4a1_1 = KnotInfo.L4a1_1
            sage: PL4a1_1 = L4a1_1.homfly_polynomial(var1='x', var2='y'); PL4a1_1
            x^-3*y^3 + 3*x^-3*y + x^-3*y^-1 - x^-5*y - x^-5*y^-1
            sage: PL4a1_1 == L4a1_1.link().homfly_polynomial(var1='x', var2='y', normalization='az')
            True

        check the skein-relation given in the doc string of :meth:`homfly_polynomial` of
        :class:`Link` (applied to one of the positive crossings of the right-handed trefoil)::

            sage: R = PK3_1.parent()
            sage: PO = R.one()
            sage: L2a1_1 = KnotInfo.L2a1_1
            sage: PL2a1_1 = L2a1_1.homfly_polynomial()
            sage: a, z = R.gens()
            sage: a*PK3_1 - ~a*PO == z*PL2a1_1
            True

        check the skein-relation from the KnotInfo description page with the original version::

            sage: pK3_1o = K3_1.homfly_polynomial(original=True); pK3_1o
            '(2*v^2-v^4)+ (v^2)*z^2'
            sage: pL2a1_1o = L2a1_1.homfly_polynomial(original=True); pL2a1_1o
            'v/z-v^3/z + v*z'

            sage: R.<v, z> = LaurentPolynomialRing(ZZ)
            sage: PO = R.one()
            sage: PK3_1o   = sage_eval(pK3_1o, locals={'v':v, 'z':z})
            sage: PL2a1_1o = sage_eval(pL2a1_1o, locals={'v':v, 'z':z})
            sage: ~v*PK3_1o - v*PO == z*PL2a1_1o
            True

        TESTS::

            all(L.homfly_polynomial() == L.link().homfly_polynomial(normalization='az') for L in KnotInfo if L.crossing_number() > 0 and L.crossing_number() < 7)
            True
        """
        if self.is_knot():
            homfly_polynomial = self[self.items.homfly_polynomial]
        else:
            homfly_polynomial = self[self.items.homflypt_polynomial]

        if original:
            return homfly_polynomial

        R = self._homfly_pol_ring(var1, var2)
        if not homfly_polynomial and self.crossing_number() == 0:
            return R.one()

        L, M = R.gens()
        lc = {'z':  M, 'v': ~L}  #
        return eval_knotinfo(homfly_polynomial, locals=lc)


    @cached_method
    def link(self, use_item=db.columns().pd_notation):
        r"""
        Return ``self`` as in instance of :class:`Link`.

        INPUT:

        - ``use_item`` -- (optional default ``self.items.pd_notation``)
          instance of :class:`KnotInfoColumns` to choose the column
          that should be used to construct the link. Allowed values
          are:
          -- self.items.pd_notation
          -- self.items.dt_notation    (only for knots)
          -- self.items.gauss_notation (only for knots)
          -- self.items.braid_notation

        .. NOTE::

            We use the PD-notation to construct ``self`` as
            default. This ensures that the number of crossings
            displayed in representation string of the link
            coincides with the crossing number as a topological
            invariant.

            But attention: The convention on how the edges are
            listed are opposite to each other

            KnotInfo: counter clockwise
            Sage:     clockwise

            Therefore, we have to take the mirror_image of the
            link!

            Furthermore, note that the mirror version may depend
            on the used KnotInfo-notation. For example for the
            knot `5_1` the Gauss- and the DT-notation refer to
            the mirror image (see example below).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.link()
            Knot represented by 3 crossings
            sage: _.braid()
            s^3
            sage: _ == K.braid()
            True

            sage: K.link(use_item=K.items.dt_notation)
            Knot represented by 3 crossings
            sage: _.braid()
            s^3

            sage: L = KnotInfo.L4a1_0
            sage: L.link()
            Link with 2 components represented by 4 crossings

            sage: L.link(use_item=L.items.dt_notation)
            Traceback (most recent call last):
            ...
            NotImplementedError: Columns.dt_notation only implemented for knots

        but observe::

            sage: L.link(use_item=L.items.braid_notation)
            Link with 2 components represented by 8 crossings

            sage: K6_1 = KnotInfo.K6_1
            sage: K6_1.link().braid() == K6_1.braid()
            False

            sage: K4_1 = KnotInfo.K4_1
            sage: K4_1.link().pd_code()
            [[4, 1, 5, 2], [8, 5, 1, 6], [6, 4, 7, 3], [2, 8, 3, 7]]
            sage: K4_1.pd_notation()
            [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]

            sage: K5_1 = KnotInfo.K5_1
            sage: K5_1.link().braid()
            s^5
            sage: K5_1.link(K5_1.items.dt_notation).braid()
            s^-5
            sage: K5_1.link(K5_1.items.gauss_notation).braid()
            s^-5
        """
        if not isinstance(use_item, KnotInfoColumns):
            raise TypeError('%s must be an instance of %s' %(use_item, KnotInfoColumns))

        if self.is_knot():
            from sage.knots.knot import Knot as Link
        else:
            from sage.knots.link import Link

        if   use_item == self.items.pd_notation:
            return Link(self.pd_notation()).mirror_image() # for mirror_image see note above
        elif use_item == self.items.braid_notation:
            return Link(self.braid())
        elif use_item == self.items.dt_notation:
            if not self.is_knot():
                raise NotImplementedError('%s only implemented for knots' %use_item)
            from sage.knots.knot import Knots
            return Knots().from_dowker_code(self.dt_notation())
        elif use_item == self.items.gauss_notation:
            if not self.is_knot():
                raise NotImplementedError('%s only implemented for knots' %use_item)
            from sage.knots.knot import Knots
            return Knots().from_gauss_code(self.gauss_notation())
        else:
            raise ValueError('Construction using %s not possible' %use_item)




KnotInfo = KnotInfoBase('KnotInfo', db.read_row_dict())
