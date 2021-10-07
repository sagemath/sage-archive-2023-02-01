from . import Feature
from join_feature import JoinFeature


class MIPBackend(Feature):
    r"""
    A feature describing whether a :class:`MixedIntegerLinearProgram` backend is available.
    """

    def _is_present(self):
        try:
            from sage.numerical.mip import MixedIntegerLinearProgram
            MixedIntegerLinearProgram(solver=self.name)
            return FeatureTestResult(self, True)
        except Exception:
            return FeatureTestResult(self, False)


class CPLEX(MIPBackend):

    def __init__(self):
        MIPBackend.__init__('cplex',
                            spkg='sage_numerical_backends_cplex')


class Gurobi(MIPBackend):

    def __init__(self):
        MIPBackend.__init__('gurobi',
                            spkg='sage_numerical_backends_gurobi')


class COIN(JoinFeature):

    def __init__(self):
        JoinFeature.__init__('sage_numerical_backends_coin',
                             [MIPBackend('coin')],
                             spkg='sage_numerical_backends_coin')
