"""
Miscellaneous operating system functions
"""

def have_program(program, path=None):
    """
    Return ``True`` if ``program`` is found in the path.

    INPUT:

    - ``program`` - a string, the name of the program to check.
    - ``path`` - string or None.  If a nonempty string, use this as
      the setting for the PATH environment variable.

    OUTPUT: bool

    This uses the shell command ``command -v``, which is Posix
    standard and seems to be present and behave consistently on
    various platforms (including Linux, BSD, Mac, Solaris, Cygwin).
    This is in contrast to ``which``, which does not use seem to use
    return codes consistently: if ``program`` does not exist, then on
    some platforms, ``which program`` has a return code of 0, on some
    a return code of 1.

    EXAMPLES::

        sage: from sage.misc.sage_ostools import have_program
        sage: have_program('ls')
        True
        sage: have_program('there_is_not_a_program_with_this_name')
        False
    """
    from subprocess import call, PIPE
    import os
    try:
        if path:
            # env is a copy of the current environment, so modifying
            # it won't affect os.environ.
            env = os.environ.copy()
            env['PATH'] = path
            return not call('command -v %s' % program, shell=True,
                            stdout=PIPE, env=env)
        else:
            return not call('command -v %s' % program, shell=True,
                            stdout=PIPE)
    except OSError:
        return False
