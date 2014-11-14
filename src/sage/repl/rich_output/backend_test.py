# -*- encoding: utf-8 -*-
r"""
Test Backend

This backend is only for doctesting purposes.

EXAMPLES::

    sage: from sage.repl.rich_output import get_display_manager
    sage: dm = get_display_manager()
    sage: from sage.repl.rich_output.backend_test import BackendTest, TestObject
    sage: doctest_backend = dm.switch_backend(BackendTest())
    sage: dm
    The Sage display manager using the test backend
    sage: dm._output_promotions
    {<class 'sage.repl.rich_output.output_basic.OutputPlainText'>: 
     <class 'sage.repl.rich_output.backend_test.TestOutputPlainText'>}
    sage: dm.displayhook(1/2)
    1/2 [TestOutputPlainText]
  
    sage: test = TestObject()
    sage: test
    _repr_ method
    sage: dm.displayhook(test)
    _rich_repr_ method [TestOutputPlainText]
    sage: dm.switch_backend(doctest_backend)   # switch back
    test
"""

from sage.structure.sage_object import SageObject
from sage.repl.rich_output.backend_base import BackendBase
from sage.repl.rich_output.output_catalog import OutputPlainText, OutputImagePng


class TestOutputPlainText(OutputPlainText):

    def __init__(self):
        """
        Backend-specific subclass of the plain text output container.

        Backends must not influence how the display system constucts
        output containers, they can only control how the output
        container is displayed. In particular, we cannot override the
        constructor (only the
        :class:`~sage.repl.rich_output.output_basic.OutputPlainText`
        constructor is used).

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_test import TestOutputPlainText
            sage: TestOutputPlainText()
            Traceback (most recent call last):
            AssertionError
        """
        assert False

    def print_to_stdout(self):
        """
        Write the data to stdout.
        """
        print('{0} [{1}]'.format(self.text.get(), self.__class__.__name__))


class TestObject(SageObject):
    """
    Test object with both :meth:`_repr_` and :meth:`_rich_repr_`
    """

    def _repr_(self):
        """
        Return string representation
        """
        return '_repr_ method'

    def _rich_repr_(self, display_manager):
        tp = display_manager.types
        if tp.OutputPlainText in display_manager.supported_output():
            return tp.OutputPlainText('_rich_repr_ method')

        
class BackendTest(BackendBase):
    
    def _repr_(self):
        return 'test'

    def supported_output(self):
        return set([TestOutputPlainText])

    def display_immediately(self, plain_text, rich_output):
        plain_text.print_to_stdout()
        if plain_text is not rich_output:
            print('rich output type: [{0}]'.format(rich_output.__class__.__name__))
