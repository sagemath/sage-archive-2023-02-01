"""
Finitely presented graded left modules over graded connected algebras
"""


# TODO:
#
# 1. Why do we need to define __bool__ and __eq__ in element.py? They should be taken care of automatically, once we define __nonzero__.
# 2. inhheritance: free modules, elements, etc., should perhaps inherit from fp versions, or maybe both should inherit from generic classes, to consolidate some methods (like degree, _repr_, others?)
# 3. Test with graded modules over other graded rings. (Should be okay, but add some doctests.)
#
# In __classcall__/__init__ for FP_Modules, allow as input a free module or a morphism of free modules?
