"""
Miscellaneous operating system functions
"""

def have_program(program, path=None):
    """
    Return ``True`` if a ``program`` executable is found in the path
    given by ``path``.

    INPUT:

    - ``program`` - a string, the name of the program to check.

    - ``path`` - string or None. Paths to search for ``program``,
      separated by ``os.pathsep``. If ``None``, use the :envvar:`PATH`
      environment variable.

    OUTPUT: bool

    EXAMPLES::

        sage: from sage.misc.sage_ostools import have_program
        sage: have_program('ls')
        True
        sage: have_program('there_is_not_a_program_with_this_name')
        False
        sage: have_program('sage', path=SAGE_ROOT)
        True
        sage: have_program('ls', path=SAGE_ROOT)
        False
    """
    import os
    if path is None:
        path = os.environ.get('PATH', "")
    for p in path.split(os.pathsep):
        try:
            if os.access(os.path.join(p, program), os.X_OK):
                return True
        except OSError:
            pass
    return False
