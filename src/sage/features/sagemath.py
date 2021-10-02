r"""
Check for SageMath Python modules
"""
from . import PythonModule


class sage__combinat(PythonModule):

    def __init__(self):
        PythonModule.__init__('sage.combinat.combinations')


class sage__graphs(PythonModule):

    def __init__(self):
        PythonModule.__init__('sage.graphs.graph')


class sage__rings__real_double(PythonModule):

    def __init__(self):
        PythonModule.__init__('sage.rings.real_double')


class sage__symbolic(PythonModule):

    def __init__(self):
        PythonModule.__init__('sage.symbolic.expression', spkg="sagemath_symbolics")
