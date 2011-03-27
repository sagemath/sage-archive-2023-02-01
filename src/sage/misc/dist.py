"""
Installing shortcut scripts
"""

import os

from subprocess import Popen, PIPE

def install_scripts(bin_directory=None):
    r"""
    Run this command as
    ``install_scripts(bin_directory)`` to create scripts
    in the given bin directory that, independently of Sage, run various
    software components included with Sage: ['gap', 'gp', 'singular',
    'maxima', 'M2', 'kash', 'mwrank', 'ipython', 'hg', 'R']

    This command:

    -  verbosely tell you which scripts it adds, and

    -  will *not* overwrite any scripts you already have in the given
       bin directory.

    INPUT:

    -  ``bin_directory`` - string; the directory into
       which to put the scripts

    OUTPUT: Verbosely prints what it is doing and creates files in
    bin_directory that are world executable and readable.

    .. note::

       You may need to run Sage as root in order to run
       ``install_scripts`` successfully, since the user running Sage
       will need write permissions on ``bin_directory``.

    AUTHORS:

    - William Stein: code / design

    - Arthur Gaer: design

    EXAMPLES::

        sage: install_scripts(SAGE_TMP)
        Checking that Sage has the command 'gap' installed
        Created script ...
    """
    if bin_directory is None:
        # We do this since the intended user of install_scripts
        # will likely be pretty clueless about how to use Sage or
        # its help system.
        import sagedoc
        print sagedoc.format(install_scripts.__doc__)
        print "USAGE: install_scripts('bin directory name')"
        return

    if not (os.path.exists(bin_directory) and os.path.isdir(bin_directory)):
        raise RuntimeError, "'%s' must exist and be a directory"%bin_directory

    for c in ['gap', 'gp', 'singular', 'maxima', 'M2', 'kash', \
              'mwrank', 'ipython', 'hg', 'R']:
        print "Checking that Sage has the command '%s' installed"%c
        p = Popen(['which', c], stdout=PIPE, stderr=PIPE)
        path = os.path.realpath(p.communicate()[0].rstrip("\n"))
        error = p.wait()
        if error:
            # the 'which' command came up empty:
            print "The command '%s' is not available; not adding shortcut"%c
        elif not path.startswith(os.path.realpath(os.environ['SAGE_ROOT'])):
            # 'which' returned a path outside of the Sage directory:
            # then the command is already installed, and we shouldn't
            # install the Sage version:
            print "The command '%s' is installed outside of Sage; not adding shortcut"%c
        else:
            # 'which' returned SAGE_ROOT/local/bin/...: create the
            # shortcut if it doesn't exist already:
            target = os.path.join(bin_directory, c)
            if os.path.exists(target):
                print "The file '%s' already exists; not adding shortcut"%(target)
            else:
                o = open(target,'w')
                o.write('#!/bin/sh\n')
                o.write('sage -%s $*\n'%c)
                print "Created script '%s'"%target
                os.system('chmod a+rx %s'%target)
        print

    print "Finished creating scripts."
    print "You need not do this again even if you upgrade or move Sage."
    print "The only requirement is that the command 'sage' is in the PATH."
