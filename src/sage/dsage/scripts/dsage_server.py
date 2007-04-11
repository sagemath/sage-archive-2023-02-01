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

from twisted.internet import reactor, error, ssl, task
from twisted.spread import pb
from twisted.python import log
from twisted.cred import portal

from sage.dsage.database.jobdb import JobDatabaseZODB, JobDatabaseSQLite
from sage.dsage.database.jobdb import DatabasePruner
from sage.dsage.database.clientdb import ClientDatabase
from sage.dsage.database.monitordb import MonitorDatabase
from sage.dsage.twisted.pb import Realm
from sage.dsage.twisted.pb import WorkerPBServerFactory
from sage.dsage.twisted.pb import _SSHKeyPortalRoot
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsChecker
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsCheckerDB
from sage.dsage.server.server import DSageServer, DSageWorkerServer
from sage.dsage.misc.constants import delimiter as DELIMITER
from sage.dsage.__version__ import version

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')

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

def write_stats(dsage_server, stats_file):
    # Put this entire thing in a try block, should not cause the server to die in any way.
    try:
        fname = os.path.join(DSAGE_DIR, stats_file)
        f = open(fname, 'w')
        f.write(dsage_server.generate_xml_stats())
        f.close()
    except Exception, msg:
        print 'Error writing stats: %s' % (msg)
        return

def startLogging(log_file):
    """This method initializes the logging facilities for the server. """
    if log_file == 'stdout':
        log.startLogging(sys.stdout)
    else:
        print "Logging to file: ", log_file
        server_log = open(log_file, 'a')
        log.startLogging(server_log)

def main():
    """
    Main execution loop of the server.

    """

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
        PUBKEY_DATABASE = os.path.expanduser(config.get('auth', 'pubkey_database'))
        STATS_FILE = config.get('general', 'stats_file')
        old_version = config.get('general', 'version')
        if version != old_version:
            raise ValueError, "Incompatible version. You have %s, need %s." % (old_version, version)
    except Exception, msg:
        print msg
        print "Error reading %s, run dsage.setup()" % conf_file
        sys.exit(-1)

    # start logging
    startLogging(LOG_FILE)

    # Job database
    jobdb = JobDatabaseSQLite()

    # Worker database
    monitordb = MonitorDatabase()

    # Client database
    clientdb = ClientDatabase()

    # Create the main DSage object
    dsage_server = DSageServer(jobdb, monitordb, clientdb, log_level=LOG_LEVEL)
    p = _SSHKeyPortalRoot(portal.Portal(Realm(dsage_server)))

    # Credentials checker
    p.portal.registerChecker(PublicKeyCredentialsCheckerDB(clientdb))

    # HACK: unsafeTracebacks should eventually be TURNED off
    client_factory = pb.PBServerFactory(p, unsafeTracebacks=True)

    # Create the looping call that will output the XML file for Dashboard
    tsk1 = task.LoopingCall(write_stats, dsage_server, STATS_FILE)
    tsk1.start(5.0, now=False)

    # Create the PBServerFactory for workers
    # Use this for unauthorized workers
    # dsage_worker = DSageWorkerServer(jobdb, log_level=LOG_LEVEL)
    # worker_factory = WorkerPBServerFactory(dsage_worker)

    dsage_server.client_factory = client_factory

    attempts = 0
    err_msg = """Could not find an open port after 50 attempts."""
    if SSL == 1:
        log.msg('Using SSL...')
        ssl_context = ssl.DefaultOpenSSLContextFactory(SSL_PRIVKEY, SSL_CERT)
        while True:
            if attempts > 50:
                log.err(err_msg)
                log.err('Last attempted port: %s' % (CLIENT_PORT))
            try:
                try:
                    import socket
                    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    s.connect(('', CLIENT_PORT))
                    succeeded = True
                except socket.error, msg:
                    succeeded = False
                if not succeeded:
                    reactor.listenSSL(CLIENT_PORT,
                                      client_factory,
                                      contextFactory = ssl_context)
                    break
                else:
                    raise SystemError
            except (SystemError, error.CannotListenError), msg:
                attempts += 1
                CLIENT_PORT += 1
    else:
        while True:
            if attempts > 50:
                log.err(err_msg)
                log.err('Last attempted port: %s' % (CLIENT_PORT))
            try:
                reactor.listenTCP(CLIENT_PORT, client_factory)
                break
            except error.CannotListenError:
                attempts += 1
                CLIENT_PORT += 1

    log.msg(DELIMITER)
    log.msg('DSAGE Server')
    log.msg('Listening on %s' % (CLIENT_PORT))
    log.msg(DELIMITER)

    reactor.run(installSignalHandlers=1)

if __name__ == "__main__":
    main()
