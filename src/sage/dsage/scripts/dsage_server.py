#!/usr/bin/env python
############################################################################
#
#   DSAGE: Distributed SAGE
#
#       Copyright (C) 2006, 2007 Yi Qiang <yqiang@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
############################################################################


import sys
import os
from optparse import OptionParser
import ConfigParser

from twisted.internet import reactor, task
from twisted.spread import pb
from twisted.python import log
from twisted.cred import portal

from sage.dsage.database.jobdb import JobDatabaseZODB, DatabasePruner
from sage.dsage.twisted.pb import Realm
from sage.dsage.twisted.pb import WorkerPBServerFactory
from sage.dsage.twisted.pb import _SSHKeyPortalRoot
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsChecker
from sage.dsage.server.server import DSageServer, DSageWorkerServer

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
# Begin reading configuration
try:
    conf_file = os.path.join(DSAGE_DIR, 'server.conf')
    config = ConfigParser.ConfigParser()
    config.read(conf_file)

    LOG_FILE = config.get('server_log', 'log_file')
    LOG_LEVEL = config.getint('server_log', 'log_level')
    SSL = config.getint('ssl', 'ssl')
    SSL_PRIVKEY = config.get('ssl', 'privkey_file')
    SSL_CERT = config.get('ssl', 'cert_file')
    WORKER_PORT = config.getint('server', 'worker_port')
    CLIENT_PORT = config.getint('server', 'client_port')
    PUBKEY_DATABASE = os.path.expanduser(config.get('auth',
                                                    'pubkey_database'))
except:
    print "Error reading %s, please fix manually run dsage.setup()" % \
    conf_file
    sys.exit(-1)
# End reading configuration

def usage():
    """Prints usage help."""

    # usage options
    usage = ['usage: %progr [options]\n',
              'required options: --jobdir\n\n',
              'Bug reports to <yqiang@gmail.com>']

    parser = OptionParser(usage=''.join(usage))
    parser.add_option("-j", "--jobdir", dest="dir",
                        help="directory containing jobs")

    (options, args) = parser.parse_args()
#    if options.dir == None:
#        parser.print_help()
#        sys.exit(0)
#        parser.error("Please specify a valid job directory with the \
#                    --jobdir flag.")

#     sys.path.append(os.path.abspath(options.dir))
    return options

def startLogging(log_file):
    """This method initializes the logging facilities for the server. """
    if log_file == 'stdout':
        log.startLogging(sys.stdout)
    else:
        print "Logging to file: ", log_file
        server_log = open(log_file, 'a')
        log.startLogging(server_log)

def main():
    """Main execution loop of the server."""

    options = usage()
    jobdb = JobDatabaseZODB()

    # Start to prune out old jobs
    jobdb_pruner = DatabasePruner(jobdb)
    prune_db = task.LoopingCall(jobdb_pruner.prune)
    prune_db.start(60*60*24.0, now=True) # start now, interval is one day

    # Create the main DSage object
    dsage_server = DSageServer(jobdb, log_level=LOG_LEVEL)
    p = _SSHKeyPortalRoot(portal.Portal(Realm(dsage_server)))

    # Get authorized keys
    p.portal.registerChecker(PublicKeyCredentialsChecker(PUBKEY_DATABASE))

    # HACK: unsafeTracebacks should eventually be TURNED off
    client_factory = pb.PBServerFactory(p, unsafeTracebacks=True)

    # Create the PBServerFactory for workers
    dsage_worker = DSageWorkerServer(jobdb, log_level=LOG_LEVEL)
    worker_factory = WorkerPBServerFactory(dsage_worker)

    # We will listen on 2 ports
    # One port that is authenticated so clients can submit new jobs
    # One port for workers to connect to to receive and submit jobs
    if SSL == 1:
        from twisted.internet import ssl
        sslContext = ssl.DefaultOpenSSLContextFactory(
                    SSL_PRIVKEY,
                    SSL_CERT)

        reactor.listenSSL(CLIENT_PORT,
                          client_factory,
                          contextFactory = sslContext)
        reactor.listenSSL(WORKER_PORT,
                          worker_factory,
                          contextFactory = sslContext)

    else:
        reactor.listenTCP(CLIENT_PORT, client_factory)
        reactor.listenTCP(WORKER_PORT, worker_factory)

    dsage_server.client_factory = client_factory
    dsage_server.worker_factory = worker_factory

    # start logging
    startLogging(LOG_FILE)

    reactor.run()

if __name__ == "__main__":
    main()
