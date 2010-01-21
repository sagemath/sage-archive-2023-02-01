r"""
Collaborative Wiki Engine

Sage includes the MoinMoin_ Wiki, "an advanced, easy to use and
extensible WikiEngine with a large community of users.  Said in a few
words, it is about collaboration on easily editable web pages."

To start your own math-typesetting-aware wiki server immediately, just
type ``wiki()`` at the command line.  For instructions on how to
upgrade existing MoinMoin wikis, please visit the `MoinMoin How-Tos`_

.. _MoinMoin: http://moinmo.in/
.. _`MoinMoin How-Tos`: http://moinmo.in/HowTo
"""

import os
import socket
import sys
import shutil

import sage.misc.misc as misc
from sage.misc.viewer import browser
from sage.server.misc import print_open_msg

join = os.path.join

def wiki_create_instance(directory='sage_wiki'):
    """
    Creates a new MoinMoin wiki.  This is a modified version of
    MoinMoin with jsMath typesetting preconfigured.

    INPUT:

      - ``directory`` -- a string (default: 'sage_wiki') the directory
        to use for the new wiki.
    """
    share = os.path.join(misc.SAGE_LOCAL, 'share', 'moin')

    directory = os.path.abspath(directory)

    if os.path.exists(directory):
        print "Directory '%s' already exists." % directory
        return

    os.makedirs(directory)

    shutil.copytree(join(share, 'data'), join(directory, 'data'))
    shutil.copytree(join(share, 'underlay'), join(directory, 'underlay'))
    shutil.copytree(join(share, 'server'), join(directory, 'server'))
    shutil.copyfile(join(share, 'config', 'wikiconfig.py'),
                    join(directory, 'wikiconfig.py'))

    wsgi_conf_name = join(directory, 'server', 'moin.wsgi')
    wsgi_conf_fd = open(wsgi_conf_name, 'r+')
    wsgi_conf = wsgi_conf_fd.read()
    wsgi_conf = wsgi_conf.replace("#sys.path.insert(0, '/path/to/wikiconfigdir')",
                                  "sys.path.insert(0, os.path.join(r'{0}'))".format(directory))
#    wsgi_conf = wsgi_conf.replace("#sys.path.insert(0, '/path/to/farmconfigdir')",
#                                  "sys.path.insert(0, join(r'{0}', 'config', 'wikifarm'))".format(directory))
    wsgi_conf_fd.seek(0)
    wsgi_conf_fd.truncate(0)
    wsgi_conf_fd.write(wsgi_conf)
    wsgi_conf_fd.close()

    os.makedirs(join(directory, 'twisted', 'plugins'))
    shutil.copyfile(join(directory, 'server', 'mointwisted.py'),
                   join(directory, 'twisted', 'plugins',
                        'mointwisted_plugin.py'))

    shutil.copyfile(wsgi_conf_name, join(directory, 'moin.py'))

def wiki(directory='sage_wiki',
         port=9000,
         address='localhost',
         fork=False):
    r"""
    Create (if necessary) and start up a MoinMoin wiki.

    The wiki will be served on the given port.

    The moin package contains a modified version of MoinMoin, which
    comes with jsMath's LaTeX typesetting preconfigured.  Use dollar
    signs (\$) to delimit text to be typeset.

    INPUT:

      - ``directory`` -- string (default: 'sage_wiki') directory to
        create/run the instance in

      - ``port`` -- integer (default: 9000) first port to listen to.
        If it is already taken, up to 256 subsequent port numbers will
        be tried.

      - ``address`` -- string (default: 'localhost') address to bind
        to

      - ``fork`` -- boolean (default: False) whether to daemonize the
        process
    """
    if not os.path.exists(directory):
        wiki_create_instance(directory)
    original_path = os.path.abspath(os.curdir)
    os.chdir(directory)

    port = int(port)

    def run(port):
        print_open_msg(address, port)
        daemonize_str = '-n'
        if fork:
            daemonize_str = ''
        e = os.system('twistd %s moin -p %d -a "%s"' % (
                daemonize_str, port, address))
        if e:
            raise socket.error

    for i in range(256):
        try:
            run(port + i)
        except socket.error:
            print "Port %s is already in use.  Trying next port..." % port
        else:
            break

    os.chdir(original_path)
    return True
