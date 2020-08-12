from . import Polynomial


def get_splitting(f, variables):
    ring = s.ring()
    s = f.set()
    for var_index in variables:
        s1 = Polynomial(s.subset1(var_index))
        s0 = Polynomial(s.subset0(var_index))
        f1 = s0 + s1
        f0 = s0
        var = ring.variable(var_index)
        if f1.constant():
            return dict(b=f1, a=Polynomial(1, ring), x=var, g=f0)
        if f0.constant():
            return dict(b=f0, a=Polynomial(0, ring), x=var, g=f1)
    return None


def nested_canalyzing_function(f, variables=None):
    f = Polynomial(f)
    if variables is None:
        variables = f.vars_as_monomial()
    if not variables.reducible_by(f.vars_as_monomial()):
        raise ValueError
    if variables != f.vars_as_monomial():
        return False
    if len(variables) == 0:
        return True
    splitting = get_splitting(f, variables)
    if splitting is None:
        return False
    rec = nested_canalyzing_function(splitting["g"], variables / splitting["x"
        ])
    if rec:
        if rec is True:
            return (splitting, )
        else:
            return (splitting, rec)
    else:
        return False
