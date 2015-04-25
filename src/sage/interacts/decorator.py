# -*- encoding: utf-8 -*-
"""
The Interact Decorator
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import inspect
import textwrap
from sage.misc.html import html


def interact(func, **kwds):
    """
    The interact decorator

    This is the global @interact-decorator. It should always return an
    interact implementation appropriate to the UI. Right now only
    SageNB is supported.

    INPUT:

    - ``func`` -- function. The interact function.

    - ``**kwds`` -- optional keyword arguments. To be passed to the
      UI interact decorator.

    EXAMPLES::

        sage: @interact
        ....: def f(x=[1,2,3]):
        ....:     print(x)
        <html>...</html>
    """
    return _sagenb_interact(func, **kwds)


def _sagenb_interact(func, **kwds):
    """
    Wrapper around the SageNB interact decorator

    Reverts to the old ``html`` behavior for SageNB interacts, where
    calling ``html(x)`` prints html to stdout instead of returning a
    :class:`sage.misc.html.HtmlFragment` instance.

    EXAMPLES::

        sage: from sage.interacts.decorator import _sagenb_interact
        sage: g = Graphics()
        sage: @_sagenb_interact
        ....: def f(x=[1,2,3], graph=g):
        ....:     html(x)
        <html>...</html>
        sage: f(2)
        <script type="math/tex">2</script>
    """
    # Evil hack: can only generate function with given argspec via exec
    argspec = inspect.getargspec(func)
    defaults = ['_argspec_defaults[{0}]'.format(i) for i in range(len(argspec.defaults))]
    signature = inspect.formatargspec(argspec[0], argspec[1], argspec[2], defaults)
    func_src = textwrap.dedent("""
    @sagenb_interact
    def sagenb_func{signature}:
        import sage.misc.html
        try:
            sage.misc.html._old_and_deprecated_behavior = True
            return func{args}
        finally:
            sage.misc.html._old_and_deprecated_behavior = False
    """).format(
        signature=signature,
        args=inspect.formatargspec(argspec[0]),
        _argspec_defaults=argspec.defaults,
    )
    from sagenb.notebook.interact import interact as sagenb_interact
    globals = dict(func=func, sagenb_interact=sagenb_interact)
    exec func_src in globals
    return globals['sagenb_func']



