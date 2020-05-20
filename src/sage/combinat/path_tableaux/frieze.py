r"""
Frieze Patterns

This is an implementation of the abstract base class
:class:`sage.combinat.pathtableau.pathtableaux`.

This implements the original frieze patterns due to Conway and Coxeter.

In this implementation we have sequences of nonnegative integers.
This follows [CoCo1]_ and [CoCo2]_

REFERENCES:

.. [CoCo1] J.H. Conway and H.S.M. Coxeter
    *Triangulated polygons and frieze patterns*,
    The Mathematical Gazette (1973) 57 p.87-94


.. [CoCo2] J.H. Conway and H.S.M. Coxeter
    *Triangulated polygons and frieze patterns (continued)*,
    The Mathematical Gazette (1973) 57 p.175-183

.. [TJ18] Thorsten Holm and  Peter Jorgensen
    *A p-angulated generalisation of Conway and Coxeter's theorem on frieze patterns*,
    International Mathematics Research Notices (2018)


AUTHORS:

- Bruce Westbury (2019): initial version
"""
#*****************************************************************************
#       Copyright (C) 2019 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from six import add_metaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.parent import Parent
from sage.categories.sets_cat import Sets
#from sage.structure.list_clone import ClonableArray
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux, CylindricalDiagram
#from sage.combinat.combinatorial_map import combinatorial_map
#from sage.combinat.tableau import Tableau, StandardTableau
from sage.rings.integer import Integer
from sage.categories.fields import Fields
from sage.rings.all import ZZ, QQ

###############################################################################

"""

EXAMPLES::

    sage: t = FriezePattern([0,1,2,1,2,3,1,0])
    sage: CylindricalDiagram(t)
     [0, 1, 2, 1, 2, 3, 1, 0]
     ['', 0, 1, 1, 3, 5, 2, 1, 0]
     ['', '', 0, 1, 4, 7, 3, 2, 1, 0]
     ['', '', '', 0, 1, 2, 1, 1, 1, 1, 0]
     ['', '', '', '', 0, 1, 1, 2, 3, 4, 1, 0]
     ['', '', '', '', '', 0, 1, 3, 5, 7, 2, 1, 0]
     ['', '', '', '', '', '', 0, 1, 2, 3, 1, 1, 1, 0]
     ['', '', '', '', '', '', '', 0, 1, 2, 1, 2, 3, 1, 0]

    sage: TestSuite(t).run()

    sage: t = FriezePattern([0,1,2,7,5,3,7,4,1,0])
    sage: CylindricalDiagram(t)
     [0, 1, 2, 7, 5, 3, 7, 4, 1, 0]
     ['', 0, 1, 4, 3, 2, 5, 3, 1, 1, 0]
     ['', '', 0, 1, 1, 1, 3, 2, 1, 2, 1, 0]
     ['', '', '', 0, 1, 2, 7, 5, 3, 7, 4, 1, 0]
     ['', '', '', '', 0, 1, 4, 3, 2, 5, 3, 1, 1, 0]
     ['', '', '', '', '', 0, 1, 1, 1, 3, 2, 1, 2, 1, 0]
     ['', '', '', '', '', '', 0, 1, 2, 7, 5, 3, 7, 4, 1, 0]
     ['', '', '', '', '', '', '', 0, 1, 4, 3, 2, 5, 3, 1, 1, 0]
     ['', '', '', '', '', '', '', '', 0, 1, 1, 1, 3, 2, 1, 2, 1, 0]
     ['', '', '', '', '', '', '', '', '', 0, 1, 2, 7, 5, 3, 7, 4, 1, 0]
    sage: TestSuite(t).run()

    sage: t = FriezePattern([0,1,3,4,5,1,0])
    sage: CylindricalDiagram(t)
     [0, 1, 3, 4, 5, 1, 0]
     ['', 0, 1, 5/3, 7/3, 2/3, 1, 0]
     ['', '', 0, 1, 2, 1, 3, 1, 0]
     ['', '', '', 0, 1, 1, 4, 5/3, 1, 0]
     ['', '', '', '', 0, 1, 5, 7/3, 2, 1, 0]
     ['', '', '', '', '', 0, 1, 2/3, 1, 1, 1, 0]
     ['', '', '', '', '', '', 0, 1, 3, 4, 5, 1, 0]

    sage: TestSuite(t).run()

This constructs the examples from [TJ18]_

    sage: K.<sqrt3> = NumberField(x^2-3)
    sage: t = FriezePattern([0,1,sqrt3,2,sqrt3,1,1,0], field=K)
    sage: CylindricalDiagram(t)
     [0, 1, sqrt3, 2, sqrt3, 1, 1, 0]
     ['', 0, 1, sqrt3, 2, sqrt3, sqrt3 + 1, 1, 0]
     ['', '', 0, 1, sqrt3, 2, sqrt3 + 2, sqrt3, 1, 0]
     ['', '', '', 0, 1, sqrt3, sqrt3 + 2, 2, sqrt3, 1, 0]
     ['', '', '', '', 0, 1, sqrt3 + 1, sqrt3, 2, sqrt3, 1, 0]
     ['', '', '', '', '', 0, 1, 1, sqrt3, 2, sqrt3, 1, 0]
     ['', '', '', '', '', '', 0, 1, sqrt3 + 1, sqrt3 + 2, sqrt3 + 2, sqrt3 + 1, 1, 0]
     ['', '', '', '', '', '', '', 0, 1, sqrt3, 2, sqrt3, 1, 1, 0]

    sage: TestSuite(t).run()

    sage: K.<sqrt2> = NumberField(x^2-2)
    sage: t = FriezePattern([0,1,sqrt2,1,sqrt2,3,2*sqrt2,5,3*sqrt2,1,0], field=K)
    sage: CylindricalDiagram(t)
     [0, 1, sqrt2, 1, sqrt2, 3, 2*sqrt2, 5, 3*sqrt2, 1, 0]
     ['', 0, 1, sqrt2, 3, 5*sqrt2, 7, 9*sqrt2, 11, 2*sqrt2, 1, 0]
     ['', '', 0, 1, 2*sqrt2, 7, 5*sqrt2, 13, 8*sqrt2, 3, sqrt2, 1, 0]
     ['', '', '', 0, 1, 2*sqrt2, 3, 4*sqrt2, 5, sqrt2, 1, sqrt2, 1, 0]
     ['', '', '', '', 0, 1, sqrt2, 3, 2*sqrt2, 1, sqrt2, 3, 2*sqrt2, 1, 0]
     ['', '', '', '', '', 0, 1, 2*sqrt2, 3, sqrt2, 3, 5*sqrt2, 7, 2*sqrt2, 1, 0]
     ['', '', '', '', '', '', 0, 1, sqrt2, 1, 2*sqrt2, 7, 5*sqrt2, 3, sqrt2, 1, 0]
     ['', '', '', '', '', '', '', 0, 1, sqrt2, 5, 9*sqrt2, 13, 4*sqrt2, 3, 2*sqrt2, 1, 0]
     ['', '', '', '', '', '', '', '', 0, 1, 3*sqrt2, 11, 8*sqrt2, 5, 2*sqrt2, 3, sqrt2, 1, 0]
     ['', '', '', '', '', '', '', '', '', 0, 1, 2*sqrt2, 3, sqrt2, 1, sqrt2, 1, sqrt2, 1, 0]
     ['', '', '', '', '', '', '', '', '', '', 0, 1, sqrt2, 1, sqrt2, 3, 2*sqrt2, 5, 3*sqrt2, 1, 0]

    sage: TestSuite(t).run()
"""

@add_metaclass(InheritComparisonClasscallMetaclass)
class FriezePattern(PathTableau):
    """
    An instance is a sequence in the ground field.
    """

    @staticmethod
    def __classcall_private__(cls, fp, field=QQ):
        """This is the preprocessing for creating friezes.

        INPUT:

        - a sequence of elements of the field ''field''

        EXAMPLES::

            sage: FriezePattern([1,2,1,2,3,1])
            [1, 2, 1, 2, 3, 1]

        TESTS::

            sage FriezePattern(2)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce sqrt3 to a rational

            sage: K.<sqrt3> = NumberField(x^2-3)
            sage: t = FriezePattern([0,1,sqrt3,2,sqrt3,1,1,0])
            Traceback (most recent call last):
            ...
            ValueError: [0, 1, sqrt3, 2, sqrt3, 1, 1, 0] is not a sequence in the field Rational Field

            sage: FriezePattern([1,2,1,2,3,1],field=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Integer Ring must be a field
        """
        if not field in Fields:
            raise ValueError("{} must be a field".format(field))

        w = None

        if isinstance(fp, (list,tuple)):
            try:
                w = tuple([field(a) for a in fp])
            except TypeError:
                raise ValueError("{} is not a sequence in the field {}".format(fp, field))

        if w is None:
            raise ValueError("invalid input {}".format(fp))

        return FriezePatterns(field)(w)

    def check(self):
        """
        There is nothing to check.
        """

    def _local_rule(self,i):
        r"""
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = FriezePattern([1,2,1,2,3,1])
            sage: t._local_rule(3)
            [1, 2, 1, 2, 3, 1]

            sage: t = FriezePattern([1,2,1,2,3,1])
            sage: t._local_rule(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not a valid integer
        """

        def _rule(x):
            """
            This is the rule on a sequence of three scalars.
            """
            return (x[0]*x[2]+1)/x[1]

        if not (i > 0 and i < len(self) ):
            raise ValueError("{} is not a valid integer".format(i))

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def is_skew(self):
        r"""
        Return ``True`` if ``self`` is skew and ``False`` if not.

        EXAMPLES::

            sage: FriezePattern([1,2,1,2,3,1]).is_skew()
            False

            sage: FriezePattern([2,2,1,2,3,1]).is_skew()
            True
        """
        return self[0] != 1

    def is_integral(self):
        r"""
        Return ``True`` if all entries of the frieze pattern are positive integers and ``False`` if not.

        EXAMPLES::

            sage: FriezePattern([0,1,2,7,5,3,7,4,1,0]).is_integral()
            True

            sage: FriezePattern([0,1,3,4,5,1,0]).is_integral()
            False

        """
        n = len(self)
        cd = CylindricalDiagram(self).diagram
        for i, a in enumerate(cd):
            v = a[i+1:n+i-2]
            try:
                v = [ Integer(k) for k in v ]
            except TypeError:
                return False
            if any(k <= 0 for k in v):
                return False
        return True

    def plot(self):
        r"""
        If ``self`` is integral then plot the triangulation.

        EXAMPLES::

            sage: FriezePattern([1,2,7,5,3,7,4,1]).plot()
            Graphics object consisting of 25 graphics primitives

        TESTS::

            sage: FriezePattern([1,2,1/7,5,3]).plot()
            Traceback (most recent call last):
            ...
            ValueError: [1, 2, 1/7, 5, 3] must be an integral frieze
        """
        if not self.is_integral():
            raise ValueError("{!s} must be an integral frieze".format(self))
        n = len(self)+1
        cd = CylindricalDiagram(self).diagram
        from sage.plot.plot import Graphics
        from sage.plot.line import line
        from sage.plot.text import text
        from sage.functions.trig import sin, cos
        from sage.all import pi
        G = Graphics()
        G.set_aspect_ratio(1.0)

        vt = [(cos(2*theta*pi/(n)), sin(2*theta*pi/(n))) for theta in range(n+1)]
        for i, p in enumerate(vt):
            G += text(str(i),[1.05*p[0],1.05*p[1]])

        for i, r in enumerate(cd):
            for j, a in enumerate(r[:n-1]):
                if a == 1:
                    G += line([vt[i],vt[j+1]])

        G.axes(False)
        return G

class FriezePatterns(PathTableaux):
    """
    The parent class for FriezePattern.
    """
    def __init__(self, field):
        """
        Initializes the abstract class of all FriezePatterns

        TESTS::

            sage: FriezePattern([1,1]).parent() # indirect test
            <sage.combinat.path_tableaux.frieze.FriezePatterns_with_category object at ...>

        """
        self.field = field

        Parent.__init__(self, category=Sets())

    Element = FriezePattern

