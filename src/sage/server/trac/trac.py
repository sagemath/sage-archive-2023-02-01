import os, sys

import sage.misc.misc as misc

from   sage.misc.viewer     import browser


def trac_create_instance(directory='sage_trac'):
    if os.system('which trac-admin 2>/dev/null 1>/dev/null') != 0:
        print "You must install the optional trac package."
        print "Try something like install_package('trac-2006.09.08'),"
        print "but note that the package name may have a different"
        print "version.  Use optional_packages() to get a list"
        print "of current package names."
        return

    e = os.system('trac-admin "%s" initenv'%directory)
    if e:
        raise RuntimeError, "Error creating trac environment."

def trac(directory='sage_trac',
         port=10000,
         address='localhost',
         open_viewer = False):
    r"""
    Create (if necessary) and start up a trac server.

    The trac server will be served on the given port.

    You must install the optional SAGE trac package. (Use
    optional_packages() for the exact name.)
    """
    if not os.path.exists(directory):
        trac_create_instance(directory)

    if open_viewer:
        cmd = 'sleep 2 && %s http://%s:%s 1>&2 >/dev/null &'%(browser(), address, port)
        os.system(cmd)

    passwd = '%s/conf/passwd'%directory
    if not os.path.exists(passwd) or len(open(passwd).read()) < 2:
        print "*"*80
        print "*** To create new accounts use the htdigest command on the command line ***"
        print "*** Note that htdigest is part of apache. ***"
        print "To add a new user with name username, do this:"
        print "  cd %s"%(os.path.abspath('%s/conf'%directory))
        print "  htdigest passwd %s <username>"%address
        print "Then restart the trac server."
        print "*"*80
        open(passwd,'w').close()

    print "Trac server started."
    print "Open your web browser to http://%s:%s"%(address, port)
    cmd ='tracd --port %s --hostname %s --auth *,%s,%s "%s" '%(port, address,
                                                passwd, address, directory)
    print cmd
    os.system(cmd)
