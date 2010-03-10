"""
A decorator to enable interacts in the Sage library.
"""

from sagenb.notebook.interact import interact, slider, range_slider
from sage.misc.html import html

def library_interact(f):
    """
    This is a decorator for using interacts in the Sage library.

    EXAMPLES::

        sage: @interacts.library_interact
        ... def f(n=5): print n
        ...
        sage: f()  # an interact appears
        <html>...</html>
    """
    def new_fun():
       # Maybe program around bug (?) in the notebook:
       html("</pre>")
       # This prints out the relevant html code to make
       # the interact appear:
       interact(f)
    # preserve the docstring
    new_fun.__doc__ = f.__doc__
    return new_fun

@library_interact
def demo(n=tuple(range(10)), m=tuple(range(10))):
    """
    This is a demo interact that sums two numbers.

    INPUT:

        - `n` -- integer slider
        - `m` -- integer slider

    EXAMPLES::

        sage: interacts.demo()
        <html>...</html>
    """
    print n+m
