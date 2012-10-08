"""
Temporary file handling

AUTHORS:

- Volker Braun, Jeroen Demeyer (2012-10-18): move these functions here
  from sage/misc/misc.py and make them secure, see :trac:`13579`.
"""

#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun@stp.dias.ie>
#       Copyright (C) 2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os
import tempfile
import atexit

def delete_tmpfiles():
    """
    Remove the directory ``SAGE_TMP``.

    TESTS:

    This is automatically run when Sage exits, test this by running a
    separate session of Sage::

        sage: from sage.tests.cmdline import test_executable
        sage: child_SAGE_TMP, err, ret = test_executable(["sage", "-c", "print SAGE_TMP"])
        sage: err, ret
        ('', 0)
        sage: os.path.exists(child_SAGE_TMP)
        False

    The parent directory should exist::

        sage: parent_SAGE_TMP = os.path.normpath(child_SAGE_TMP + '/..')
        sage: os.path.isdir(parent_SAGE_TMP)
        True
    """
    import shutil
    from sage.misc.misc import SAGE_TMP
    shutil.rmtree(str(SAGE_TMP), ignore_errors=True)


# Run when Python shuts down
atexit.register(delete_tmpfiles)


#################################################################
# temporary directory
#################################################################

def tmp_dir(name="dir_", ext=""):
    r"""
    Create and return a temporary directory in
    ``$HOME/.sage/temp/hostname/pid/``

    The temporary directory is deleted automatically when Sage exits.

    INPUT:

    - ``name`` -- (default: ``"dir_"``) A prefix for the directory name.

    - ``ext`` -- (default: ``""``) A suffix for the directory name.

    OUTPUT:

    The absolute path of the temporary directory created, with a
    trailing slash (or whatever the path separator is on your OS).

    EXAMPLES::

        sage: d = tmp_dir('dir_testing_', '.extension')
        sage: d   # random output
        '/home/username/.sage/temp/hostname/7961/dir_testing_XgRu4p.extension/'
        sage: os.chdir(d)
        sage: _ = open('file_inside_d', 'w')

    Temporary directories are unaccessible by other users::

        sage: os.stat(d).st_mode & 0o077
        0
    """
    from sage.misc.misc import SAGE_TMP
    tmp = tempfile.mkdtemp(prefix=name, suffix=ext, dir=str(SAGE_TMP))
    name = os.path.abspath(tmp)
    return name + os.sep


#################################################################
# temporary filename
#################################################################

def tmp_filename(name="tmp_", ext=""):
    r"""
    Create and return a temporary file in
    ``$HOME/.sage/temp/hostname/pid/``

    The temporary file is deleted automatically when Sage exits.

    .. warning::

        If you need a particular file extension always use
        ``tmp_filename(ext=".foo")``, this will ensure that the file
        does not yet exist. If you were to use
        ``tmp_filename()+".foo"``, then you might overwrite an
        existing file!

    INPUT:

    - ``name`` -- (default: ``"tmp_"``) A prefix for the file name.

    - ``ext`` -- (default: ``""``) A suffix for the file name.

    OUTPUT:

    The absolute path of the temporary file created.

    EXAMPLES::

        sage: fn = tmp_filename('just_for_testing_', '.extension')
        sage: fn  # random
        '/home/username/.sage/temp/hostname/8044/just_for_testing_tVVHsn.extension'
        sage: _ = open(fn, 'w')

    Temporary files are unaccessible by other users::

        sage: os.stat(fn).st_mode & 0o077
        0
    """
    from sage.misc.misc import SAGE_TMP
    handle, tmp = tempfile.mkstemp(prefix=name, suffix=ext, dir=str(SAGE_TMP))
    os.close(handle)
    name = os.path.abspath(tmp)
    return name


def graphics_filename(ext='png'):
    """
    Return the next available canonical filename for a plot/graphics
    file.
    """
    i = 0
    while os.path.exists('sage%d.%s'%(i,ext)):
        i += 1
    filename = 'sage%d.%s'%(i,ext)
    return filename
