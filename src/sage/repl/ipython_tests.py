'''
Tests for the IPython integration

First, test the pinfo magic for Python code. This is what IPython
calls when you ask for the single-questionmark help, like `foo?` ::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell(u'from sage.repl.ipython_tests import dummy')
    sage: shell.run_cell(u'%pinfo dummy')
    Type:            function
    String form:     <function dummy at 0x...>
    File:            /.../sage/repl/ipython_tests.py
    Definition:      dummy(argument, optional=None)
    Docstring:
       Dummy Docstring Title
    <BLANKLINE>
       Dummy docstring explanation.
    <BLANKLINE>
       INPUT:
    <BLANKLINE>
       * "argument" -- anything. Dummy argument.
    <BLANKLINE>
       * "optional" -- anything (optional). Dummy optional.
    <BLANKLINE>
       EXAMPLES:
    ...
    Class docstring:
    function(code, globals[, name[, argdefs[, closure]]])
    <BLANKLINE>
    Create a function object from a code object and a dictionary. The
    optional name string overrides the name from the code object. The
    optional argdefs tuple specifies the default argument values. The
    optional closure tuple supplies the bindings for free variables.
    Init docstring:  x.__init__(...) initializes x; see help(type(x)) for signature


Next, test the pinfo magic for Cython code::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell(u'from sage.tests.stl_vector import stl_int_vector')
    sage: shell.run_cell(u'%pinfo stl_int_vector')
    Type:           type
    String form:    <type 'sage.tests.stl_vector.stl_int_vector'>
    File:           /.../sage/tests/stl_vector.pyx
    Docstring:
       Example class wrapping an STL vector
    <BLANKLINE>
       EXAMPLES:
    ...
    <BLANKLINE>
    Init docstring: x.__init__(...) initializes x; see help(type(x)) for signature


Next, test the pinfo2 magic for Python code. This is what IPython
calls when you ask for the double-questionmark help, like `foo??` ::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell(u'from sage.repl.ipython_tests import dummy')
    sage: shell.run_cell(u'%pinfo2 dummy')
    Type:            function
    String form:     <function dummy at 0x...>
    File:            /.../sage/repl/ipython_tests.py
    Definition:      dummy(argument, optional=None)
    Source:
    def dummy(argument, optional=None):
        """
        Dummy Docstring Title
    ...

Next, test the pinfo2 magic for Cython code::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell(u'from sage.tests.stl_vector import stl_int_vector')
    sage: shell.run_cell(u'%pinfo2 stl_int_vector')
    Type:           type
    String form:    <type 'sage.tests.stl_vector.stl_int_vector'>
    File:           /.../sage/tests/stl_vector.pyx
    Source:
    cdef class stl_int_vector(SageObject):
        """
        Example class wrapping an STL vector
    ...
'''

def dummy(argument, optional=None):
    """
    Dummy Docstring Title

    Dummy docstring explanation.

    INPUT:

    - ``argument`` -- anything. Dummy argument.

    - ``optional`` -- anything (optional). Dummy optional.

    EXAMPLES:

        sage: from sage.repl.ipython_tests import dummy
        sage: dummy(1)
        'Source code would be here' 
    """
    return 'Source code would be here'
