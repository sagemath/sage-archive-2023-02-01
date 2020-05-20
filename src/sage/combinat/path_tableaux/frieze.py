r"""
Unimodular Frieze Patterns

This is an implementation of the abstract base class
:class:`sage.combinat.pathtableau.pathtableaux`.

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
from sage.structure.list_clone import ClonableArray
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
#from sage.combinat.combinatorial_map import combinatorial_map
#from sage.combinat.tableau import Tableau, StandardTableau
from sage.rings.integer import Integer
from sage.rings.all import QQ

###############################################################################

"""

EXAMPLES::

    sage: t = UnimodularFriezePattern([0,1,2,1,2,3,1,0])

    sage: SkewTableau(t.cylindrical_diagram()).pp()
      0  1  2  1  2  3  1  0
      .  0  1  1  3  5  2  1  0
      .  .  0  1  4  7  3  2  1  0
      .  .  .  0  1  2  1  1  1  1  0
      .  .  .  .  0  1  1  2  3  4  1  0
      .  .  .  .  .  0  1  3  5  7  2  1  0
      .  .  .  .  .  .  0  1  2  3  1  1  1  0
      .  .  .  .  .  .  .  0  1  2  1  2  3  1  0

    sage: TestSuite(t).run()

    sage: t = UnimodularFriezePattern([0,1,2,7,5,3,7,4,1,0])

    sage: SkewTableau(t.cylindrical_diagram()).pp()
      0  1  2  7  5  3  7  4  1  0
      .  0  1  4  3  2  5  3  1  1  0
      .  .  0  1  1  1  3  2  1  2  1  0
      .  .  .  0  1  2  7  5  3  7  4  1  0
      .  .  .  .  0  1  4  3  2  5  3  1  1  0
      .  .  .  .  .  0  1  1  1  3  2  1  2  1  0
      .  .  .  .  .  .  0  1  2  7  5  3  7  4  1  0
      .  .  .  .  .  .  .  0  1  4  3  2  5  3  1  1  0
      .  .  .  .  .  .  .  .  0  1  1  1  3  2  1  2  1  0
      .  .  .  .  .  .  .  .  .  0  1  2  7  5  3  7  4  1  0

    sage: TestSuite(t).run()

    sage: t = UnimodularFriezePattern([0,1,3,4,5,1,0])

    sage: SkewTableau(t.cylindrical_diagram()).pp()
      0  1  3  4  5  1  0
      .  0  15/37/32/3  1  0
      .  .  0  1  2  1  3  1  0
      .  .  .  0  1  1  45/3  1  0
      .  .  .  .  0  1  57/3  2  1  0
      .  .  .  .  .  0  12/3  1  1  1  0
      .  .  .  .  .  .  0  1  3  4  5  1  0


    sage: TestSuite(t).run()

This constructs the examples from [TJ18]_

    sage: K.<sqrt3> = NumberField(x^2-3)
    sage: t = UnimodularFriezePattern([0,1,sqrt3,2,sqrt3,1,1,0], field=K)

    sage: SkewTableau(t.cylindrical_diagram()).pp()
      0  1sqrt3  2sqrt3  1  1  0
      .  0  1sqrt3  2sqrt3sqrt3 + 1  1  0
      .  .  0  1sqrt3  2sqrt3 + 2sqrt3  1  0
      .  .  .  0  1sqrt3sqrt3 + 2  2sqrt3  1  0
      .  .  .  .  0  1sqrt3 + 1sqrt3  2sqrt3  1  0
      .  .  .  .  .  0  1  1sqrt3  2sqrt3  1  0
      .  .  .  .  .  .  0  1sqrt3 + 1sqrt3 + 2sqrt3 + 2sqrt3 + 1  1  0
      .  .  .  .  .  .  .  0  1sqrt3  2sqrt3  1  1  0

    sage: TestSuite(t).run()

    sage: K.<sqrt2> = NumberField(x^2-2)
    sage: t = UnimodularFriezePattern([0,1,sqrt2,1,sqrt2,3,2*sqrt2,5,3*sqrt2,1,0], field=K)

    sage: SkewTableau(t.cylindrical_diagram()).pp()
      0  1sqrt2  1sqrt2  32*sqrt2  53*sqrt2  1  0
      .  0  1sqrt2  35*sqrt2  79*sqrt2 112*sqrt2  1  0
      .  .  0  12*sqrt2  75*sqrt2 138*sqrt2  3sqrt2  1  0
      .  .  .  0  12*sqrt2  34*sqrt2  5sqrt2  1sqrt2  1  0
      .  .  .  .  0  1sqrt2  32*sqrt2  1sqrt2  32*sqrt2  1  0
      .  .  .  .  .  0  12*sqrt2  3sqrt2  35*sqrt2  72*sqrt2  1  0
      .  .  .  .  .  .  0  1sqrt2  12*sqrt2  75*sqrt2  3sqrt2  1  0
      .  .  .  .  .  .  .  0  1sqrt2  59*sqrt2 134*sqrt2  32*sqrt2  1  0
      .  .  .  .  .  .  .  .  0  13*sqrt2 118*sqrt2  52*sqrt2  3sqrt2  1  0
      .  .  .  .  .  .  .  .  .  0  12*sqrt2  3sqrt2  1sqrt2  1sqrt2  1  0
      .  .  .  .  .  .  .  .  .  .  0  1sqrt2  1sqrt2  32*sqrt2  53*sqrt2  1  0

    sage: TestSuite(t).run()
"""


@add_metaclass(InheritComparisonClasscallMetaclass)
class UnimodularFriezePattern(PathTableau):
    """
    An instance is the sequence of nonnegative
    integers.
    """

    @staticmethod
    def __classcall_private__(cls, fp, field=QQ):
        """This is the preprocessing for creating paths.

        INPUT:

            - a sequence of nonnegative integers

        EXAMPLES::

            sage: UnimodularFriezePattern([1,2,1,2,3,1])
            [1, 2, 1, 2, 3, 1]

        """
        w = None

        if isinstance(fp, (list,tuple)):
            try:
                w = tuple([field(a) for a in fp])
            except TypeError:
                raise ValueError("%s is not a sequence of integers" % fp)

        if w is None:
            raise ValueError("invalid input %s" % fp)

        return UnimodularFriezePatterns(field)(w)

    def check(self):

        n = len(self)
        if any(a < 0 for a in self):
           raise ValueError( "%s has a negative entry" % (str(self)) )

    def _local_rule(self,i):
        """
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = UnimodularFriezePattern([1,2,1,2,3,1])
            sage: t._local_rule(3)
            [1, 2, 1, 2, 3, 1]
        """

        def _rule(x):
            """
            This is the rule on a sequence of three letters.
            """
            return (x[0]*x[2]+1)/x[1]

        if not (i > 0 and i < len(self) ):
            raise ValueError("%d is not a valid integer" % i)

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def is_skew(self):
        """
        Return ``True`` if ``self`` is skew and ``False`` if not.

        EXAMPLES::

            sage: UnimodularFriezePattern([1,2,1,2,3,1]).is_skew()
            False

            sage: UnimodularFriezePattern([2,2,1,2,3,1]).is_skew()
            True
        """
        return self[0] != 1

    def is_integral(self):
        """
        Return ``True`` if all entries of the frieze pattern are positive integers and ``False`` if not.

        EXAMPLES::

            sage: UnimodularFriezePattern([0,1,2,7,5,3,7,4,1,0]).is_integral()
            True

            sage: UnimodularFriezePattern([0,1,3,4,5,1,0]).is_integral()
            False

        """
        n = len(self)
        cd = self.cylindrical_diagram()
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
        """
        If ``self`` is integral then plot the triangulation.

        EXAMPLES::

            sage: UnimodularFriezePattern([0,1,2,7,5,3,7,4,1,0]).plot()
            Graphics object consisting of 24 graphics primitives

        """
        if not self.is_integral():
            raise ValueError("must be an integral frieze")
        n = len(self)
        cd = self.cylindrical_diagram()
        from sage.plot.plot import Graphics
        from sage.plot.line import line
        from sage.plot.text import text
        from sage.functions.trig import sin, cos
        from sage.all import pi
        G = Graphics()
        G.set_aspect_ratio(1.0)

        vt = [(cos(2*theta*pi/(n-1)), sin(2*theta*pi/(n-1))) for theta in range(n-1)]
        for i, p in enumerate(vt):
            G += text(str(i),[1.05*p[0],1.05*p[1]])

        for i, r in enumerate(cd):
            for j, a in enumerate(r[:n-1]):
                if a == 1:
                    G += line([vt[i],vt[j]])

        G.axes(False)
        return G

class UnimodularFriezePatterns(PathTableaux):

    def __init__(self, field):

        self._field = field

        Parent.__init__(self, category=Sets())

    def field(self):
        return self._field

    Element = UnimodularFriezePattern

