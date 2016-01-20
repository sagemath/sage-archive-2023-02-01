"""
Mixin For Extra Tab Completions

EXAMPLES::

    sage: from sage.interfaces.tab_completion import ExtraTabCompletion
    sage: class Foo(ExtraTabCompletion):
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
