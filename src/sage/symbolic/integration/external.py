from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR

def maxima_integrator(expression, v, a=None, b=None):
    """
        sage: from sage.symbolic.integration.external import maxima_integrator
        sage: maxima_integrator(sin(x), x)
        -cos(x)
        sage: maxima_integrator(cos(x), x)
        sin(x)
        sage: f(x) = function('f', x)
        sage: maxima_integrator(f(x), x)
        integrate(f(x), x)
    """
    if not isinstance(expression, Expression):
        expression = SR(expression)
    if a is None:
        result = expression._maxima_().integrate(v)
    else:
        try:
            result = expression._maxima_().integrate(v, a, b)
        except TypeError, error:
            s = str(error)
            if "divergent" in s or 'Principal Value' in s:
                raise ValueError, "Integral is divergent."
            else:
                raise
    return result.sage()

def sympy_integrator(expression, v, a=None, b=None):
    """
        sage: from sage.symbolic.integration.external import sympy_integrator
        sage: sympy_integrator(sin(x), x)
        -cos(x)
        sage: sympy_integrator(cos(x), x)
        sin(x)
    """
    import sympy
    ex = expression._sympy_()
    v = v._sympy_()
    if a is None:
        result = sympy.integrate(ex, v)
    else:
        result = sympy.integrate(ex, (v, a._sympy_(), b._sympy_()))
    return result._sage_()

def mma_free_integrator(expression, v, a=None, b=None):
    """
        sage: from sage.symbolic.integration.external import mma_free_integrator
        sage: mma_free_integrator(sin(x), x) # optional - requires internet
        -cos(x)
    """
    import urllib, re
    # We need to integrate against x
    vars = [str(x) for x in expression.variables()]
    if any(len(x)>1 for x in vars):
        raise NotImplementedError, "Mathematica online integrator can only handle single letter variables."
    x = SR.var('x')
    if repr(v) != 'x':
        for i in range(ord('a'), ord('z')+1):
            if chr(i) not in vars:
                shadow_x = SR.var(chr(i))
                break
        expression = expression.subs({x:shadow_x}).subs({dvar: x})
    params = urllib.urlencode({'expr': expression._mathematica_init_(), 'random': 'false'})
    page = urllib.urlopen("http://integrals.wolfram.com/index.jsp", params).read()
    page = page[page.index('"inputForm"'):page.index('"outputForm"')]
    page = re.sub("\s", "", page)
    mexpr = re.match(r".*Integrate.*==</em><br/>(.*)</p>", page).groups()[0]
    try:
        ans = SR(mexpr.lower().replace('[', '(').replace(']', ')'))
        if repr(v) != 'x':
            ans = ans.subs({x:v}).subs({shadow_x:x})
        return ans
    except TypeError:
        raise ValueError, "Unable to parse: %s" % mexpr
