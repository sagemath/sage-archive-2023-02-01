"""
Installing shortcut scripts
"""

import os

def install_scripts(bin_directory=None):
    r"""
    Run this command as
    ``install_scripts(bin_directory)`` to create scripts
    in the given bin directory that, independently of Sage, run various
    software components included with Sage: ['gap', 'gp', 'singular',
    'maxima', 'M2', 'kash', 'mwrank', 'ipython', 'hg', 'hgmerge', 'R']

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
    """
    if bin_directory is None:
        # We do this since the intended user of install_scripts
        # will likely be pretty clueless about how to use SAGE or
        # its help system.
        import sagedoc
        print sagedoc.format(install_scripts.__doc__)
        print "USAGE: install_scripts('bin directory name')"
        return

    if not (os.path.exists(bin_directory) and os.path.isdir(bin_directory)):
        raise RuntimeError, "'%s' must exist and be a directory"%bin_directory

    for c in ['gap', 'gp', 'singular', 'maxima', 'M2', 'kash', \
              'mwrank', 'ipython', 'hg', 'hgmerge', 'R']:
        print "\nChecking that SAGE has the command '%s' installed"%c
        if os.system('which %s > /dev/null'%c):
            print "The command '%s' is not available; not adding shortcut"%c
        else:
            target = '%s/%s'%(bin_directory, c)
            if os.path.exists(target):
                print "** Not creating script for '%s' since the file '%s' already exists"%(c, target)
            else:
                o = open(target,'w')
                o.write('#!/bin/sh\n')
                o.write('sage -%s $*\n'%c)
                print "Created script '%s'"%target
                os.system('chmod a+rx %s'%target)

    print "Finished creating scripts."
    print "You need not do this again even if you upgrade or move SAGE."
    print "The only requirement is that the command 'sage' is in the PATH."
