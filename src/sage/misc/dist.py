"""
Installing shortcut scripts
"""

from __future__ import print_function
import os

def install_scripts(directory=None, ignore_existing=False):
    r"""
    Running ``install_scripts(directory)`` creates scripts in the
    given directory that run various software components included with
    Sage.  Each of these scripts essentially just runs ``sage --CMD``
    where ``CMD`` is also the name of the script:

    - 'gap' runs GAP
    - 'gp' runs the PARI/GP interpreter
    - 'hg' runs Mercurial
    - 'ipython' runs IPython
    - 'maxima' runs Maxima
    - 'mwrank' runs mwrank
    - 'R' runs R
    - 'singular' runs Singular
    - 'sqlite3' runs SQLite version 3
    - 'kash' runs Kash if it is installed (Kash is an optional Sage
      package)
    - 'M2' runs Macaulay2 if it is installed (Macaulay2 is an
      experimental Sage package)

    This command:

    -  verbosely tells you which scripts it adds, and

    -  will *not* overwrite any scripts you already have in the given
       directory.

    INPUT:

    - ``directory`` - string; the directory into which to put the
      scripts.  This directory must exist and the user must have write
      and execute permissions.

    - ``ignore_existing`` - bool (optional, default False): if True,
      install script even if another version of the program is in your
      path.

    OUTPUT: Verbosely prints what it is doing and creates files in
    ``directory`` that are world executable and readable.

    .. note::

       You may need to run ``sage`` as ``root`` in order to run
       ``install_scripts`` successfully, since the user running
       ``sage`` needs write permissions on ``directory``.  Note
       that one good candidate for ``directory`` is
       ``'/usr/local/bin'``, so from the shell prompt, you could run ::

           sudo sage -c "install_scripts('/usr/local/bin')"

    .. note::

       Running ``install_scripts(directory)`` will be most helpful if
       ``directory`` is in your path.

    AUTHORS:

    - William Stein: code / design

    - Arthur Gaer: design

    - John Palmieri: revision, 2011-07 (:trac:`11602`)

    EXAMPLES::

        sage: install_scripts(str(SAGE_TMP), ignore_existing=True)
        Checking that Sage has the command 'gap' installed
        ...
    """
    if directory is None:
        # We do this since the intended user of install_scripts
        # will likely be pretty clueless about how to use Sage or
        # its help system.
        import sagedoc
        print(sagedoc.format(install_scripts.__doc__))
        print("USAGE: install_scripts('directory')")
        return

    if not os.path.exists(directory):
        print("Error: '{}' does not exist.".format(directory))
        return

    if not os.path.isdir(directory):
        print("Error: '{}' is not a directory.".format(directory))
        return

    if not (os.access(directory, os.W_OK) and os.access(directory, os.X_OK)):
        print("Error: you do not have write permission for '{}'.".format(directory))
        return

    from sage.misc.sage_ostools import have_program
    from sage.env import SAGE_LOCAL
    script_created = False
    SAGE_BIN = os.path.join(SAGE_LOCAL, 'bin')
    # See if 'directory' is already in PATH, and then remove
    # SAGE_LOCAL/bin from PATH so that we can later check whether
    # cmd is available outside of Sage.
    PATH = os.environ['PATH'].split(os.pathsep)
    PATH = [d for d in PATH if os.path.exists(d)]
    dir_in_path = any([os.path.samefile(directory, d) for d in PATH])
    PATH = os.pathsep.join([d for d in PATH if not
                            os.path.samefile(d, SAGE_BIN)])
    for cmd in ['gap', 'gp', 'hg', 'ipython', 'maxima',
              'mwrank', 'R', 'singular', 'sqlite3', 'M2', 'kash']:
        print("Checking that Sage has the command '{}' installed".format(cmd))
        # Check to see if Sage includes cmd.
        cmd_inside_sage = have_program(cmd, path=SAGE_BIN)
        cmd_outside_sage = have_program(cmd, path=PATH)
        if not cmd_inside_sage:
            print("The command '{}' is not available as part ".format(cmd)
                   + "of Sage; not creating script.")
            print()
            continue
        if cmd_outside_sage:
            print("The command '{}' is installed outside of Sage;".format(cmd), end=' ')
            if not ignore_existing:
                print("not creating script.")
                print()
                continue
            print("trying to create script anyway...")
        else:
            print("Creating script for '{}'...".format(cmd))
        # Install shortcut.
        target = os.path.join(directory, cmd)
        if os.path.exists(target):
            print("The file '{}' already exists; not adding script.".format(target))
        else:
            o = open(target,'w')
            o.write('#!/bin/sh\n')
            o.write('exec sage --%s "$@"\n'%cmd)
            o.close()
            print("Created script '{}'".format(target))
            os.system('chmod a+rx {}'.format(target))
            script_created = True
        print()

    if script_created:
        print("Finished creating scripts.")
        print()
        print("You need not do this again even if you upgrade or move Sage.")
        print("The only requirement is that your PATH contains both")
        print("'{}' and the directory containing the command 'sage'.".format(directory))
        if not dir_in_path:
            print()
            print("Warning: '{}' is not currently in your PATH.".format(directory))
            print()
    else:
        print("No scripts created.")
