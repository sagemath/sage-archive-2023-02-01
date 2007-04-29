"""
Catalogue of Functions and of Parent and Element Constructors
"""

# TODO -- this is not done yet.

class Catalogue:
    """
    Use the new object and tab completion to create new SAGE objects.  It gathers
    together and organizers the various SAGE constructors.
    """
    pass

new = Catalogue()

class Parent:
    """
    Parent structures, i.e., container objects.
    """
    pass
parent = Parent()
new.parent = parent

class Element:
    """
    Elements of various parent structures.
    """
    pass
element = Element()
new.element = element

class Category:
    """
    Categories.
    """
    pass
category = Category()
#new.category = category

class Function:
    """
    Functions
    """
    pass
function = Function()
new.function = function

#######################################################
# Rings
#######################################################

class Ring:
    pass
ring = Ring()

import sage.rings.all
ring.ZZ = sage.rings.all.ZZ
ring.QQ = sage.rings.all.QQ
ring.PolynomialRing = sage.rings.all.PolynomialRing
ring.LaurentSeriesRing = sage.rings.all.LaurentSeriesRing

parent.ring = ring

class RingElts:
    pass
ring = RingElts()
ring.Integer = sage.rings.all.Integer

element.ring = ring

#######################################################
# Groups
#######################################################

import sage.groups.all
class Group:
    pass
group = Group()
group.PermutationGroup = sage.groups.all.PermutationGroup
parent.group = group

#######################################################
# Matrix Groups
#######################################################
class MatrixGroup:
    """
    Matrix groups, i.e., groups that are viewed as a set of matrices.
    """
    pass
matrix_group = MatrixGroup()
matrix_group.GL = sage.groups.all.GL
matrix_group.SL = sage.groups.all.SL
matrix_group.Sp = sage.groups.all.Sp
matrix_group.SU = sage.groups.all.SU
matrix_group.GU = sage.groups.all.GU
matrix_group.SO = sage.groups.all.SO
matrix_group.GO = sage.groups.all.GO
matrix_group.MatrixGroup = sage.groups.all.MatrixGroup
parent.group.matrix_group = matrix_group

#######################################################
# Modules
#######################################################

#######################################################
# Graphs
#######################################################

#######################################################
# Modular forms
#######################################################

#######################################################
# Functions
#######################################################
class Elementary:
    """
    Elementary functions.
    """
elementary = Elementary()
function.elementary = elementary
import sage.calculus.all
elementary.sin = sage.calculus.all.sin
elementary.cos = sage.calculus.all.cos
elementary.log = sage.calculus.all.log
