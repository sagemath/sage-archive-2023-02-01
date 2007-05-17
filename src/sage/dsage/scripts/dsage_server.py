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
#
# import gc
# gc.set_debug(gc.DEBUG_LEAK)

import sys
import os
from optparse import OptionParser
import socket

from twisted.internet import reactor, error, ssl, task
from twisted.spread import pb
from twisted.python import log
from twisted.cred import portal

from sage.dsage.database.jobdb import JobDatabaseSQLite
from sage.dsage.database.clientdb import ClientDatabase
from sage.dsage.database.monitordb import MonitorDatabase
from sage.dsage.twisted.pb import Realm
from sage.dsage.twisted.pb import _SSHKeyPortalRoot
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsCheckerDB
from sage.dsage.server.server import DSageServer
from sage.dsage.misc.constants import DELIMITER
from sage.dsage.misc.constants import DSAGE_DIR
from sage.dsage.misc.config import get_bool
from sage.dsage.misc.config import check_dsage_dir

def usage():
    """
    Prints usage help.

    """

    # usage options
    usage = ['usage: %prog [options]\n',
              'Bug reports to <yqiang@gmail.com>']

    parser = OptionParser(usage=''.join(usage))
    parser.add_option('-p', '--port',
                      dest='port',
                      type='int',
                      default=8081,
                      help='port to listen on')
    parser.add_option('-f', '--logfile',
                      dest='logfile',
                      default=os.path.join(DSAGE_DIR, 'server.log'),
                      help='log file. default=~/.sage/dsage/server.log')
    parser.add_option('-l', '--loglevel',
                      dest='loglevel',
                      type='int',
                      default=0,
                      help='log level, higher means more verbose')
    parser.add_option('--statsfile',
                      dest='statsfile',
                      default=os.path.join(DSAGE_DIR, 'dsage.xml'),
                      help='xml file for dsage statistics. ' +
                           'default=~/.sage/dsage/dsage.xml')
    parser.add_option('--ssl',
                      dest='ssl',
                      default=True,
                      help='enable or disable ssl')
    parser.add_option('-k', '--privkey',
                      dest='privkey',
                      default=os.path.join(DSAGE_DIR, 'cacert.pem'),
                      help='private key for ssl certificate')
    parser.add_option('-c', '--cert',
                      dest='cert',
                      default=os.path.join(DSAGE_DIR, 'pubcert.pem'),
                      help='ssl certificate')
    parser.add_option('-d', '--dbfile',
                      dest='dbfile',
                      default=os.path.join(DSAGE_DIR, 'dsage.db'),
                      help='database file')
    parser.add_option('--job_failures',
                      dest='job_failure_threshold',
                      type='int',
                      default=3,
                      help='sets the threshold for job failures')
    parser.add_option('--noblock',
                      dest='noblock',
                      action='store_true',
                      default=False,
                      help='tells that the server was ' +
                           'started in blocking mode')

    (options, args) = parser.parse_args()

    return options

def write_stats(dsage_server, stats_file):
    try:
        fname = os.path.join(DSAGE_DIR, stats_file)
        f = open(fname, 'w')
        f.write(dsage_server.generate_xml_stats())
        f.close()
    except Exception, msg:
        print 'Error writing stats: %s' % (msg)
        return

def create_manhole():
    """
    This is a manhole backdoor to inspect a running server.
    ONLY turn it on for debugging, otherwise you can get pwn3d.

    """

    from twisted.manhole import telnet
    factory = telnet.ShellFactory()
    factory.username = 'yqiang'
    factory.password = 'foo'
    port = reactor.listenTCP(2000, factory)

    return port

def startLogging(log_file):
    """
    This method initializes the logging facilities for the server.

    """

    check_dsage_dir()

    if log_file == 'stdout':
        log.startLogging(sys.stdout)
        log.msg('WARNING: DSAGE Server ONLY logging to stdout!')
    else:
        server_log = open(log_file, 'a')
        log.startLogging(sys.stdout)
        log.startLogging(server_log)
        log.msg("DSAGE Server: Logging to file: ", log_file)

def main(options):
    """
    Main execution loop of the server.

    """

    # config = get_conf('server')
    # LOG_FILE = config['log_file']
    # LOG_LEVEL = config['log_level']
    # SSL = get_bool(config['ssl'])
    # SSL_PRIVKEY = config['privkey_file']
    # SSL_CERT = config['cert_file']
    # CLIENT_PORT = int(config['client_port'])
    # PUBKEY_DATABASE = os.path.expanduser(config['pubkey_database'])
    # STATS_FILE = config['stats_file']

    LOG_FILE = options.logfile
    LOG_LEVEL = options.loglevel
    SSL = get_bool(options.ssl)
    SSL_PRIVKEY = options.privkey
    SSL_CERT = options.cert
    CLIENT_PORT = options.port
    STATS_FILE = options.statsfile
    DB_FILE = options.dbfile
    FAILURE_THRESHOLD = options.job_failure_threshold

    # start logging
    startLogging(LOG_FILE)

    # Job database
    jobdb = JobDatabaseSQLite(db_file=DB_FILE,
                              job_failure_threshold=FAILURE_THRESHOLD,
                              log_file=LOG_FILE, log_level=LOG_LEVEL)
    # Worker database
    monitordb = MonitorDatabase(db_file=DB_FILE,
                                log_file=LOG_FILE, log_level=LOG_LEVEL)

    # Client database
    clientdb = ClientDatabase()

    # Create the main DSage object
    dsage_server = DSageServer(jobdb, monitordb,
                               clientdb, log_level=LOG_LEVEL)
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
    err_msg = "Could not find an open port after 50 attempts."
    NEW_CLIENT_PORT = CLIENT_PORT
    while True:
        if attempts > 50:
            log.err(err_msg)
            log.err('Last attempted port: %s' % (NEW_CLIENT_PORT))
            sys.exit(-1)
        try:
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.connect(('', NEW_CLIENT_PORT))
                port_used = True
            except socket.error, msg:
                port_used = False
            if not port_used:
                if SSL:
                    ssl_context = ssl.DefaultOpenSSLContextFactory(
                                    SSL_PRIVKEY, SSL_CERT)
                    reactor.listenSSL(NEW_CLIENT_PORT,
                                      client_factory,
                                      contextFactory = ssl_context)
                    break
                else:
                    reactor.listenTCP(NEW_CLIENT_PORT, client_factory)
                    break
            else:
                raise SystemError('Trying to bind to open port: '
                                  + '%s.' % (NEW_CLIENT_PORT))
        except (SystemError, error.CannotListenError):
            attempts += 1
            NEW_CLIENT_PORT += 1
        except Exception, msg:
            print 'Exception: ', msg
            if SSL:
                print 'Error starting server with SSL enabled, please ' + \
                      'check your configuration'
            else:
                print 'Error starting server, please check your configuration'
            sys.exit()

    if CLIENT_PORT != NEW_CLIENT_PORT:
        log.msg(DELIMITER)
        log.msg("***NOTICE***")
        log.msg("Changing listening port in server.conf " +
                "to %s" % (NEW_CLIENT_PORT))
        log.msg(DELIMITER)

    log.msg(DELIMITER)
    log.msg('DSAGE Server')
    log.msg('Started with PID: %s' % (os.getpid()))
    if SSL:
        log.msg('Using SSL: True')
    else:
        log.msg('Using SSL: False')
    log.msg('Listening on port: %s' % (NEW_CLIENT_PORT))
    log.msg(DELIMITER)

    # Code below can be turned on to do countrefs
    # from sage.dsage.misc.countrefs import logInThread
    # logInThread(n=15)
    # reactor.callWhenRunning(create_manhole)
    # import pdb; pdb.set_trace()
    if options.noblock:
        reactor.run(installSignalHandlers=0)
    else:
        reactor.run(installSignalHandlers=1)

if __name__ == "__main__":
    options = usage()
    main(options)
