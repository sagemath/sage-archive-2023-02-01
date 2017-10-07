"""
Mixin For Extra Tab Completions

The :class:`ExtraTabCompletion` class helps you to extend the tab
completions for objects. This is used in interfaces with third-party
interpreters where we want, for example, ``gap.[TAB]`` to also list
GAP commands and not only the Python methods of the gap interface. To
use it, just inherit (usually by multiple inheritance) from
:class:`ExtraTabCompletion` and implement a ``_tab_completion`` method
that returns a list of strings. These strings are then included in the
tab completions on instances. It is up to you to make these "fake"
method names work, usually by implementing ``__getattr__``.

EXAMPLES::

    sage: from sage.interfaces.tab_completion import ExtraTabCompletion
    sage: class Foo(ExtraTabCompletion, object):
    ....:      def a(self):
    ....:          return 1
    ....:      def _tab_completion(self):
    ....:          return ['c', 'd']
    sage: f = Foo()
    sage: f.b = 2
    sage: sorted(dir(f))
    [..., '_tab_completion', 'a', 'b', 'c', 'd']
"""

import warnings


class ExtraTabCompletion(object):

    def __dir__(self):
        """
        Add to the dir() output

        This is used by IPython to read off the tab completions.

        EXAMPLES::

            sage: from sage.interfaces.tab_completion import ExtraTabCompletion
            sage: obj = ExtraTabCompletion()
            sage: dir(obj)
            Traceback (most recent call last):
            ...
            NotImplementedError: <class 'sage.interfaces.tab_completion.ExtraTabCompletion'> must implement _tab_completion() method
        """
        try:
            tab_fn = self._tab_completion
        except AttributeError:
            raise NotImplementedError(
                '{0} must implement _tab_completion() method'.format(self.__class__))
        return dir(self.__class__) + self.__dict__.keys() + tab_fn()
