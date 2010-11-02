"""
A library of interacts
"""

from sagenb.notebook.interact import interact, slider, range_slider, input_box
from sage.all import sin, plot, point, html, show, latex, SR,exp
x=SR.var('x')

from sage.misc.decorators import sage_wraps
from sage.misc.html import html

def library_interact(f):
    """
    This is a decorator for using interacts in the Sage library.

    EXAMPLES::

        sage: @interacts.library.library_interact
        ... def f(n=5): print n
        ...
        sage: f()  # an interact appears
        <html>...</html>
    """
    @sage_wraps(f)
    def library_wrapper():
       # Maybe program around bug (?) in the notebook:
       html("</pre>")
       # This prints out the relevant html code to make
       # the interact appear:
       interact(f)
    return library_wrapper

@library_interact
def demo(n=tuple(range(10)), m=tuple(range(10))):
    """
    This is a demo interact that sums two numbers.

    INPUT:

        - `n` -- integer slider
        - `m` -- integer slider

    EXAMPLES::

        sage: interacts.library.demo()
        <html>...</html>
    """
    print n+m



@library_interact
def taylor_polynomial(f=input_box(sin(x)*exp(-x)), order=slider(range(1,13))):
    """
    An interact which illustrates the Taylor polynomial approximation
    of various orders around `x=0`.

        - `f` -- function expression
        - ``order`` -- integer slider

    EXAMPLES::

        sage: interacts.calculus.taylor_polynomial()
        <html>...</html>
    """
    x0  = 0
    p   = plot(f,(x,-1,5), thickness=2)
    dot = point((x0,f(x=x0)),pointsize=80,rgbcolor=(1,0,0))
    ft = f.taylor(x,x0,order)
    pt = plot(ft,(-1, 5), color='green', thickness=2)
    html('$f(x)\;=\;%s$'%latex(f))
    html('$\hat{f}(x;%s)\;=\;%s+\mathcal{O}(x^{%s})$'%(x0,latex(ft),order+1))
    show(dot + p + pt, ymin = -.5, ymax = 1)
