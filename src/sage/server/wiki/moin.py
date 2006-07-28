import os, sys

import sage.misc.misc as misc

def wiki_create_instance(directory='sage_wiki'):
    try:
        from MoinMoin.server.standalone import StandaloneConfig, run
    except ImportError:
        print "You must install the optional moin package."
        print "Try something like install_package('moin-1.5.4'),"
        print "but note that the package name may have a different"
        print "version.  Use optional_packages() to get a list"
        print "of current package names."
        return

    share = '%s/share/moin'%misc.SAGE_LOCAL

    if os.path.exists(directory):
        print "Directory '%s' already exists."%directory
        return

    os.makedirs(directory)
    os.system('cp -r %s/data %s/'%(share,directory))
    os.system('cp -r %s/underlay %s/'%(share,directory))
    os.system('cp %s/config/wikiconfig.py %s/'%(share,directory))
    os.system('cp %s/server/moin.py %s/'%(share,directory))


def wiki(directory='sage_wiki', port=9000, address='localhost'):
    try:
        from MoinMoin.server.standalone import StandaloneConfig, run
    except ImportError:
        print "You must install the optional moin package."
        print "Try something like install_package('moin-1.5.4'),"
        print "but note that the package name may have a different"
        print "version.  Use optional_packages() to get a list"
        print "of current package names."
        return False

    if not os.path.exists(directory):
        wiki_create_instance(directory)
    os.chdir(directory)

    moin = '%s/share/moin/'%misc.SAGE_LOCAL
    sys.path.insert(0, '%s/config/'%moin)
    the_port = port

    class Config(StandaloneConfig):
        # Server name
        # Used to create .log, .pid and .prof files
        name = 'moin'

        # Path to moin shared files (default '/usr/share/moin/wiki/htdocs')
        # If you installed with --prefix=PREFIX, use 'PREFIX/share/moin/wiki/htdocs'
        docs = '%s/htdocs'%moin

        # Port (default 8000)
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

    run(Config)
    return True
