#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This is an impementation of the Category PathTableaux.

In this implementation we have sequences of partitions. These are in
bijection with dual semistandard tableaux. This gives an effective
version of operations on tableaux constructed using jeu-de-taquin.
In the standard constructions of these operations one usually assumes
the tableau is standard.

For rectification and evacuation the operations here
agree with the standard construction. For promotion the construction
here agrees with the standard construction on rectangular standard
tableaux, but, in general, they are different.

The operations here also give the Bender-Knuth involutions and
dual equivalence graphs.

AUTHORS:

- Bruce Westbury (2018): initial version

#*****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


"""

from six import add_metaclass

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent

from sage.categories.pathtableaux import PathTableaux


@add_metaclass(InheritComparisonClasscallMetaclass)
class DualSemistandardTableau(ClonableArray):
    """
       An instance is the sequence of partitions correspond to the
       chain of partitions of a dual semistandard skew tableau.
       
    The acceptable inputs are:
        - a sequence such that each term defines a partition
        - a semistandard skew tableau

    EXAMPLES:
        
    sage: DualSemistandardTableau([[],[1],[2],[2,1]])
    [[], [1], [2], [2, 1]]

    sage: t = SkewTableau([[None,None,None,4,4,5,6,7],[None,2,4,6,7,7,7],[None,4,5,8,8,9],[None,6,7,10],[None,8,8,11],[None],[4]])
    sage: DualSemistandardTableau(t)
    [[6, 1, 1], [6, 1, 1], [6, 2, 1], [6, 2, 1], [7, 3, 2, 1, 1], [7, 3, 3, 1, 1, 1], [7, 4, 3, 2, 1, 1, 1], [7, 4, 4, 2, 2, 2, 2, 1], [7, 5, 5, 3, 3, 2, 2, 1], [7, 5, 5, 3, 3, 3, 2, 1], [7, 5, 5, 4, 3, 3, 2, 1], [7, 5, 5, 5, 3, 3, 2, 1]]

    """
    @staticmethod
    def __classcall_private__(self, ot):
        
        w = None

        if isinstance(ot,(SkewTableau,SemistandardTableau)):
            w = ot.conjugate().to_chain()
            
        if isinstance(ot,(list,tuple)):
            try:
                w = tuple([ Partition(a) for a in ot ])
            except TypeError:
                raise ValueError("%s is not a sequence of partitions." % str(ot) )

        if w == None:
            raise ValueError( "Sorry, not sorry; I don't know what to do with %s." % str(ot) )  
                        
        return DualSemistandardTableaux()(w)

    def _hash_(self):
        return hash(tuple(map(tuple, self)))

    def check(self):
        n = len(self)
        for i in range(n-1):
            h = self[i]
            t = self[i+1]
            if not t.contains(h):
                raise ValueError( "%s must contain %s" % (str(t),str(h)) )
            for r, s in zip(h,t):
                if s > r+1:
                    raise ValueError( "%s / %s is not a vertical strip" % (str(t),str(h)) )
            for a in t[len(h):]:
                if a > 1:
                    raise ValueError( "%s / %s is not a vertical strip" % (str(t),str(h)) )
    @staticmethod
    def _rule(x):
        y = map(list,x)
        m = max([ len(u) for u in y ])
        z = map( lambda u: vector(u + [0]*(m-len(u)) ), y )
        result = list(z[0]-z[1]+z[2])
        result.sort(reverse=true)
        return Partition(result)

    def evaluation(self):
        z = [ p.size() for p in self ]
        return [ z[i+1] - z[i] for i in range(len(self)-1) ]
    
    def to_tableau(self):
        """
        Returns the conjugate skew tableau. This will be semistandard.
        """
        ch = [ p.conjugate() for p in self]
        try:
            return SkewTableau(chain=ch)
        except TypeError:
            return SemistandardTableau(chain=ch)
    
    def is_skew(self):
        """
        Returns True if Tableau is skew and False if not.

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.is_skew()
        False
        """
        return self[0] != Partition([])

    def rectify(self):
        pass
    
    def check_bender_knuth(self,i):
        lhs = self.local_rule(i).to_tableau()
        rhs = self.to_tableau().bender_knuth_involution(i)
        return lhs == rhs
    
    def check_rectify(self):
        lhs = self.rectify().to_tableau()
        rhs = self.to_tableau().rectify()
        return lhs == rhs
        
    def check_evacuation(self):
        lhs = self.evacuation().to_tableau()
        rhs = self.to_tableau().evacuation()
        return lhs == rhs
"""
I wanted to put in checks of the claims I made. However SkewTableaux
does not have the operations of promotion or evacuation

"""
###############################################################################

class DualSemistandardTableaux(UniqueRepresentation,Parent):

    @staticmethod
    def __classcall_private__(cls):
        return super(DualSemistandardTableaux, cls).__classcall__(cls)

    def __init__(self):

        Parent.__init__(self, category=PathTableaux())

    def __contains__(self, ot):

        return isinstance(ot, (list, tuple, DualSemistandardTableau))

    def _element_constructor_(self, ot, check=True):

        if isinstance(ot, DualSemistandardTableaux) and ot.parent() == self:
            return ot

        return self.element_class(self, list(ot))

    Element = DualSemistandardTableau

