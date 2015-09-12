
import sage
from ring import SymbolicRing, SR
from expression import Expression


class GenericSymbolicSubring(SymbolicRing):

    def __init__(self, vars):
        super(GenericSymbolicSubring, self).__init__()
        self._vars_ = set(vars)


    def _repr_variables_(self):
        return ', '.join(str(v) for v in sorted(self._vars_, key=str))


    def is_variable_valid(self, var):
        raise NotImplementedError('Not implemented in this abstract base class')


    def _element_constructor_(self, x):
        expression = super(GenericSymbolicSubring, self)._element_constructor_(x)
        if not all(self.is_variable_valid(var)
                   for var in expression.variables()):
            raise ValueError('%s is not contained in %s' % (x, self))


class SymbolicSubringAcceptingVars(GenericSymbolicSubring):

    def _repr_(self):
        return 'Symbolic Subring accepting variables %s' % \
            (self._repr_variables_())


    def is_variable_valid(self, var):
        return var in self._vars_


class SymbolicSubringRejectingVars(GenericSymbolicSubring):

    def _repr_(self):
        return 'Symbolic Subring rejecting variables %s' % \
            (self._repr_variables_())


    def is_variable_valid(self, var):
        return var not in self._vars_


class SymbolicConstantsSubring(GenericSymbolicSubring):

    def _repr_(self):
        return 'Symbolic Constants Subring'


    def is_variable_valid(self, var):
        return False
