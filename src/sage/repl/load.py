"""
Load Python, Sage, Cython and Magma files in Sage
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os
import base64

def is_loadable_filename(filename):
    """
    Returns whether a file can be loaded into Sage.  This checks only
    whether its name ends in one of the supported extensions ``.py``,
    ``.pyx``, ``.sage``, ``.spyx``, and ``.m``.  Note: :func:`load`
    assumes the latter signifies a Magma file.

    INPUT:

    - ``filename`` - a string

    OUTPUT:

    - a boolean

    EXAMPLES::

        sage: sage.repl.load.is_loadable_filename('foo.bar')
        False
        sage: sage.repl.load.is_loadable_filename('foo.c')
        False
        sage: sage.repl.load.is_loadable_filename('foo.sage')
        True
        sage: sage.repl.load.is_loadable_filename('foo.m')
        True
    """
    if filename.endswith(('.py', '.pyx', '.sage', '.spyx', '.m')):
        return True
    return False


def load_cython(name):
    """
    Helper function to load a Cython file.

    INPUT:

    - ``name`` -- filename of the Cython file

    OUTPUT:

    - A string with Python code to import the names from the compiled
      module.
    """
    from sage.misc.cython import cython
    mod, dir = cython(name, compile_message=True, use_cache=True)
    import sys
    sys.path.append(dir)
    return 'from {} import *'.format(mod)


def load(filename, globals, attach=False):
    """
    Executes a file in the scope given by ``globals``.  The
    ``filename`` itself is also evaluated in the scope.  If the name
    starts with ``http://``, it is treated as a URL and downloaded.

    .. NOTE::

        For Cython files, the situation is more complicated --
        the module is first compiled to a temporary module ``t`` and
        executed via::

            from t import *

    INPUT:

    - ``filename`` -- a string; a .py, .sage, .pyx, etc., filename,
      URL, or expression that evaluates to one

    - ``globals`` -- a string:object dictionary; the context in which
      to evaluate the ``filename`` and exec its contents

    - ``attach`` -- a boolean (default: False); whether to add the
      file to the list of attached files

    EXAMPLES:

    Note that ``.py`` files are *not* preparsed::

        sage: t = tmp_filename(ext='.py')
        sage: open(t,'w').write("print 'hi', 2/3; z = -2/7")
        sage: z = 1
        sage: sage.repl.load.load(t, globals())
        hi 0
        sage: z
        -1

    A ``.sage`` file *is* preparsed::

        sage: t = tmp_filename(ext='.sage')
        sage: open(t,'w').write("print 'hi', 2/3; z = -2/7")
        sage: z = 1
        sage: sage.repl.load.load(t, globals())
        hi 2/3
        sage: z
        -2/7

    Cython files are *not* preparsed::

        sage: t = tmp_filename(ext='.pyx')
        sage: open(t,'w').write("print 'hi', 2/3; z = -2/7")
        sage: z = 1
        sage: sage.repl.load.load(t, globals())
        Compiling ...
        hi 0
        sage: z
        -1

    If the file isn't a Cython, Python, or a Sage file, a ValueError
    is raised::

        sage: sage.repl.load.load('a.foo',globals())
        Traceback (most recent call last):
        ...
        ValueError: argument (='a.foo') to load or attach must have extension py, pyx, sage, spyx, or m

    A filename given as an expression get evaluated.  This ensures
    that ``load DATA+'foo.sage'`` works in the Notebook, say::

        sage: t=tmp_filename(ext='.py'); open(t,'w').write("print 'hello world'")
        sage: sage.repl.load.load(t, globals())
        hello world

    We load a file given at a remote URL::

        sage: sage.repl.load.load('http://wstein.org/loadtest.py', globals())  # optional - internet
        hi from the net
        5

    We can load files using secure http (https)::

        sage: sage.repl.load.load('https://github.com/jasongrout/minimum_rank/raw/minimum_rank_1_0_0/minrank.py', globals())  # optional - internet

    We attach a file::

        sage: t = tmp_filename(ext='.py')
        sage: open(t,'w').write("print 'hello world'")
        sage: sage.repl.load.load(t, globals(), attach=True)
        hello world
        sage: t in attached_files()
        True

    You can't attach remote URLs (yet)::

        sage: sage.repl.load.load('http://wstein.org/loadtest.py', globals(), attach=True)  # optional - internet
        Traceback (most recent call last):
        ...
        NotImplementedError: you can't attach a URL

    The default search path for loading and attaching files is the
    current working directory, i.e., ``'.'``.  But you can modify the
    path with :func:`load_attach_path`::

        sage: sage.repl.attach.reset(); reset_load_attach_path()
        sage: load_attach_path()
        ['.']
        sage: t_dir = tmp_dir()
        sage: fullpath = os.path.join(t_dir, 'test.py')
        sage: open(fullpath, 'w').write("print 37 * 3")
        sage: load_attach_path(t_dir)
        sage: attach('test.py')
        111
        sage: sage.repl.attach.reset(); reset_load_attach_path() # clean up

    or by setting the environment variable ``SAGE_LOAD_ATTACH_PATH``
    to a colon-separated list before starting Sage::

        $ export SAGE_LOAD_ATTACH_PATH="/path/to/my/library:/path/to/utils"
        $ sage
        sage: load_attach_path()          # not tested
        ['.', '/path/to/my/library', '/path/to/utils']

    Make sure that load handles filenames with spaces in the name or path::

        sage: t = tmp_filename(ext=' b.sage'); open(t,'w').write("print 2")
        sage: sage.repl.load.load(t, globals())
        2
    """
    if attach:
        from sage.repl.attach import add_attached_file

    try:
        filename = eval(filename, globals)
    except Exception:
        # First check if the file exists. The filename may have spaces in
        # its name, but more importantly modified_attached_files calls load
        # with the absolute file path and that may contain spaces in the path
        # As a side effect, this also allows file names with spaces in
        # them, but currently I don't see a way to disallow this case.
        if not os.path.exists(filename) and not os.path.isabs(filename):
            # handle multiple input files separated by spaces, which was
            # maybe a bad idea, but which we have to handle for backwards
            # compatibility.
            v = filename.split()
            if len(v) > 1:
                for file in v:
                    load(file, globals, attach=attach)
                return

    filename = filename.strip()

    if filename.lower().startswith(('http://', 'https://')):
        if attach:
            # But see http://en.wikipedia.org/wiki/HTTP_ETag for how
            # we will do this.
            # http://www.diveintopython.net/http_web_services/etags.html
            raise NotImplementedError("you can't attach a URL")
        from remote_file import get_remote_file
        filename = get_remote_file(filename, verbose=False)

    if not is_loadable_filename(filename):
        raise ValueError('argument (=%r) to load or attach must have extension py, pyx, sage, spyx, or m' % filename)

    fpath = os.path.expanduser(filename)
    if os.path.isabs(fpath):
        if not os.path.exists(fpath):
            raise IOError('did not find file %r to load or attach' % filename)
    else:
        from sage.repl.attach import load_attach_path
        for path in load_attach_path():
            fpath = os.path.join(path, filename)
            fpath = os.path.expanduser(fpath)
            if os.path.exists(fpath):
                break
        else:
            raise IOError('did not find file %r in load / attach search path' \
                % filename)

    if fpath.endswith('.py'):
        if attach:
            add_attached_file(fpath)
        with open(fpath) as f:
            code = compile(f.read(), fpath, 'exec')
            exec(code, globals)
    elif fpath.endswith('.sage'):
        from sage.repl.attach import load_attach_mode
        from sage.repl.preparse import preparse_file_named, preparse_file
        load_debug_mode, attach_debug_mode = load_attach_mode()
        if (attach and attach_debug_mode) or ((not attach) and load_debug_mode):
            # Preparse to a file to enable tracebacks with
            # code snippets. Use preparse_file_named to make
            # the file name appear in the traceback as well.
            # See Trac 11812.
            if attach:
                add_attached_file(fpath)
            with open(preparse_file_named(fpath)) as f:
                code = compile(f.read(), preparse_file_named(fpath), 'exec')
                exec(code, globals)
        else:
            # Preparse in memory only for speed.
            if attach:
                add_attached_file(fpath)
            exec(preparse_file(open(fpath).read()) + "\n", globals)
    elif fpath.endswith('.spyx') or fpath.endswith('.pyx'):
        if attach:
            add_attached_file(fpath)
        exec(load_cython(fpath), globals)
    elif fpath.endswith('.m'):
        # Assume magma for now, though maybe .m is used by maple and
        # mathematica too, and we should really analyze the file
        # further.
        s = globals['magma'].load(fpath)
        i = s.find('\n'); s = s[i+1:]
        print(s)


def load_wrap(filename, attach=False):
    """
    Encodes a load or attach command as valid Python code.

    INPUT:

    - ``filename`` - a string; the argument to the load or attach
      command

    - ``attach`` - a boolean (default: False); whether to attach
      ``filename``, instead of loading it

    OUTPUT:

    - a string

    EXAMPLES::

        sage: sage.repl.load.load_wrap('foo.py', True)
        'sage.repl.load.load(sage.repl.load.base64.b64decode("Zm9vLnB5"),globals(),True)'
        sage: sage.repl.load.load_wrap('foo.sage')
        'sage.repl.load.load(sage.repl.load.base64.b64decode("Zm9vLnNhZ2U="),globals(),False)'
        sage: sage.repl.load.base64.b64decode("Zm9vLnNhZ2U=")
        'foo.sage'
    """
    return 'sage.repl.load.load(sage.repl.load.base64.b64decode("{}"),globals(),{})'.format(
        base64.b64encode(filename), attach)
