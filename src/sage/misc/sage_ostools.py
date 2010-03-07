"""
Miscellaneous operating system functions
"""

def have_program(program):
    """
    Return ``True`` if ``program`` is found in the path.

    INPUT:

    - ``program`` - a string, the name of the program to check

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
    try:
        return not call('command -v ' + program, shell=True, stdout=PIPE, stderr=PIPE)
    except OSError:
        return False

