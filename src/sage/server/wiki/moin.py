r"""
Wiki Interactive Web Page.

SAGE includes the Moin Moin Wiki interactive web page system standard.
To start your own math-typesetting-aware wiki server immediately, just
type \code{wiki()} at the command line.

The Moin Moin Wiki ``is an advanced, easy to use and extensible
WikiEngine with a large community of users. Said in a few words, it is
about collaboration on easily editable web pages.''
"""

import os, socket, sys

import sage.misc.misc as misc

from   sage.misc.viewer     import browser

# if you change the default sage_wiki, you must also change local/bin/sage-wiki
def wiki_create_instance(directory='sage_wiki'):
    from MoinMoin.server.standalone import StandaloneConfig, run

    share = '%s/share/moin'%misc.SAGE_LOCAL

    if os.path.exists(directory):
        print "Directory '%s' already exists."%directory
        return

    os.makedirs(directory)
    os.system('cp -r %s/data %s/'%(share,directory))
    os.system('cp -r %s/underlay %s'%(share,directory))
    os.system('cp %s/config/wikiconfig.py %s/'%(share,directory))
    os.system('cp %s/server/moin.py %s/'%(share,directory))
    R = open('%s/moin.py'%directory,'r').read()
    R = R.replace('/path/to/wikiconfig',directory)
    open('%s/moin.py'%directory,'w').write(R)

def wiki(directory='sage_wiki',
         port=9000,
         address='localhost',
         open_viewer = False):
    r"""
    Create (if necessary) and start up a Moin Moin wiki.

    The wiki will be served on the given port.

    The moin package contains a modified version of moin moin, which
    comes with jsmath latex typesetting preconfigured; use dollar
    signs to typeset.
    """
    sys.path.insert(0, os.path.abspath(directory))

    from MoinMoin.server.standalone import StandaloneConfig, run

    if not os.path.exists(directory):
        wiki_create_instance(directory)
    os.chdir(directory)

    moin = '%s/share/moin/'%misc.SAGE_LOCAL
    the_port = int(port)

    class Config(StandaloneConfig):
        # Server name
        # Used to create .log, .pid and .prof files
        name = 'moin'

        # Path to moin shared files (default '/usr/share/moin/wiki/htdocs')
        # If you installed with --prefix=PREFIX, use 'PREFIX/share/moin/wiki/htdocs'
        docs = '%s/htdocs'%moin

        # Port
        port = the_port

        # To serve privileged port under 1024 you will have to run as root.
        # Interface (default 'localhost')
        # The default will listen only to localhost.
        # '' will listen to any interface
        interface = address

        # Log (default commented)
        # Log is written to stderr or to a file you specify here.
        ## logPath = name + '.log'

        # Server class (default ThreadPoolServer)
        # 'ThreadPoolServer' - create a constant pool of threads, simplified
        # Apache worker mpm.
        # 'ThreadingServer' - serve each request in a new thread. Much
        # slower for static files.
        # 'ForkingServer' - serve each request on a new child process -
        # experimental, slow.
        # 'SimpleServer' - server one request at a time. Fast, low
        # memory footprint.
        # If you set one of the threading servers and threads are not
        # available, the server will fallback to ForkingServer. If fork is
        # not available, the server will fallback to SimpleServer.
        serverClass = 'ThreadPoolServer'

        # Thread limit (default 10)
        # Limit the number of threads created. Ignored on non threaded servers.
        threadLimit = 10

        # Request queue size (default 50)
        # The size of the socket listen backlog.
        requestQueueSize = 50

        # Properties
        # Allow overriding any request property by the value defined in
        # this dict e.g properties = {'script_name': '/mywiki'}.
        properties = {}

        # Memory profile (default commented)
        # Useful only if you are a developer or interested in moin memory usage
        # A memory profile named 'moin--2004-09-27--01-24.log' is
        # created each time you start the server.
        ## from MoinMoin.util.profile import Profiler
        ## memoryProfile = Profiler(name, requestsPerSample=100, collect=0)

        # Hotshot profile (default commented)
        # Not compatible with threads - use with SimpleServer only.
        ## hotshotProfile = name + '.prof'


    if open_viewer:
        cmd = '%s http://%s:%s 1>&2 >/dev/null &'%(browser(), address, port)
        os.system(cmd)

    for i in range(256):
        try:
            run(Config)
        except socket.error:
            print "Port %s is already in use.  Trying next port..."%(Config.port)
            Config.port += 1
        else:
            break


    return True
