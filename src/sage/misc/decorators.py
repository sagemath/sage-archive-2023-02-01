"""
Decorators

Python decorators for use in Sage.

AUTHORS:

- Tim Dumol (5 Dec 2009) -- initial version.
"""
from functools import partial, wraps

class specialize:
    r"""
    A decorator generator that returns a decorator that in turn
    returns a specialized function for function ``f``. In other words,
    it returns a function that acts like ``f`` with arguments
    ``*args`` and ``**kwargs`` supplied.

    INPUT:

    - ``*args``, ``**kwargs`` -- arguments to specialize the function for.

    OUTPUT:

    - a decorator that accepts a function ``f`` and specializes it
      with ``*args`` and ``**kwargs``

    EXAMPLES::

        sage: f = specialize(5)(lambda x, y: x+y)
        sage: f(10)
        15
        sage: f(5)
        10
        sage: @specialize("Bon Voyage")
        ... def greet(greeting, name):
        ...    print "{0}, {1}!".format(greeting, name)
        sage: greet("Monsieur Jean Valjean")
        Bon Voyage, Monsieur Jean Valjean!
        sage: greet(name = 'Javert')
        Bon Voyage, Javert!
    """
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, f):
        return wraps(f)(partial(f, *self.args, **self.kwargs))
