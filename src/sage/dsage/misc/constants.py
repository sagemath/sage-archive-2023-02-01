import os

DELIMITER = '=' * 50
DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
TMP_WORKER_FILES = os.path.join(DSAGE_DIR, 'tmp_worker_files')
SERVER_LOG = os.path.join(DSAGE_DIR, 'server.log')
WORKER_LOG = os.path.join(DSAGE_DIR, 'worker.log')
DSAGE_LOCAL = os.path.join(os.getenv('SAGE_ROOT'), 'local/dsage')
DSAGE_DB_DIR = os.path.join(DSAGE_DIR, 'db')
DSAGE_DB = os.path.join(DSAGE_DB_DIR, 'dsage.db')

# These are the twisted tac files to be used with twistd
SERVER_TAC = """import sys
import os
from optparse import OptionParser
import socket
import sqlite3

from twisted.application import internet, service
from twisted.internet import reactor, error, task
from twisted.spread import pb
from twisted.python import log
from twisted.cred import portal
from twisted.web2 import server, http, resource, channel, static

from gnutls.constants import *
from gnutls.crypto import *
from gnutls.errors import *
from gnutls.interfaces.twisted import X509Credentials

# from sage.dsage.database.jobdb import JobDatabaseSQLite as JobDatabase
from sage.dsage.database.jobdb import JobDatabaseSA as JobDatabase
from sage.dsage.database.clientdb import ClientDatabaseSA as ClientDatabase
from sage.dsage.database.workerdb import WorkerDatabaseSA as WorkerDatabase
from sage.dsage.database.db_config import init_db_sa as init_db
from sage.dsage.twisted.pb import Realm, _SSHKeyPortalRoot
from sage.dsage.twisted.pubkeyauth import PublicKeyCredentialsCheckerDB
from sage.dsage.server.server import DSageServer
from sage.dsage.misc.constants import DELIMITER, DSAGE_DIR
from sage.dsage.misc.config import check_dsage_dir
from sage.dsage.misc.misc import find_open_port, random_string
from sage.dsage.web.web_server import Toplevel, UserManagement

global open_ports
open_ports = find_open_port()

def start_web_server(dsage_server, dsage_service):
    web_server_port = open_ports.next()
    web_server_port = open_ports.next()
    dsage_server.web_port = web_server_port
    top_level = Toplevel(dsage_server, web_server_port)
    secure_url = random_string(length=20)
    dsage_server.secure_url = secure_url
    # Creates the User Management page on a random url which is
    # printed to the console
    top_level.putChild(secure_url, UserManagement(dsage_server))
    site = server.Site(top_level)

    s = internet.TCPServer(web_server_port, channel.HTTPFactory(site))
    s.setServiceParent(dsage_service)

    dsage_server.web_server_port = web_server_port

    return web_server_port, secure_url

def start_dsage_server(dsage_service):
    DB_FILE = '%s'
    FAILURE_THRESHOLD = %s
    SSL = %s
    LOG_LEVEL = %s
    LOG_FILE = '%s'
    SSL_PRIVKEY = '%s'
    SSL_CERT = '%s'
    SERVER_PORT = %s
    TESTING = %s

    Session = init_db(DB_FILE)
    jobdb = JobDatabase(Session)
    workerdb = WorkerDatabase(Session)
    clientdb = ClientDatabase(Session)

    # Create the main DSage object
    dsage_server = DSageServer(jobdb, workerdb, clientdb, log_level=LOG_LEVEL)

    # Credentials checker
    p = _SSHKeyPortalRoot(portal.Portal(Realm(dsage_server)))
    p.portal.registerChecker(PublicKeyCredentialsCheckerDB(clientdb))

    # HACK: unsafeTracebacks should eventually be TURNED off
    client_factory = pb.PBServerFactory(p, unsafeTracebacks=True)
    dsage_server.client_factory = client_factory

    if TESTING:
        print 'Going into testing mode...'
        dsage_server._testing = True
        clientdb._add_test_client()
    dsage_server.port = SERVER_PORT
    if SSL:
        ## This for OpenSSL, SAGE uses GNUTLS now
        ## ssl_context = ssl.DefaultOpenSSLContextFactory(
        ##                 SSL_PRIVKEY, SSL_CERT)
        ## reactor.listenSSL(NEW_CLIENT_PORT,
        ##                   client_factory,
        ##                   contextFactory = ssl_context)
        cert = X509Certificate(open(SSL_CERT).read())
        key = X509PrivateKey(open(SSL_PRIVKEY).read())
        cred = X509Credentials(cert, key)
        cred.verify_peer = False # Do not verify certs
        cred.session_params.compressions = (COMP_LZO,
                                            COMP_DEFLATE,
                                            COMP_NULL)
        s = internet.TLSServer(SERVER_PORT, client_factory, cred)
        s.setServiceParent(dsage_service)
        dsage_server.ssl = True
    else:
        s = internet.TCPServer(SERVER_PORT, client_factory)
        s.setServiceParent(dsage_service)
        dsage_server.ssl = False

    return dsage_server, SERVER_PORT, SSL

def print_info(dsage_server):
    log.msg(DELIMITER)
    log.msg('DSAGE Server')
    log.msg('Started with PID:', os.getpid())
    if dsage_server.ssl:
        log.msg('Using SSL: True')
    else:
        log.msg('Using SSL: False')
    log.msg('Listening on port:', dsage_server.port)
    log.msg(DELIMITER)
    log.msg('DSAGE Web Server')
    log.msg('http://localhost:' + str(dsage_server.web_port))
    log.msg('')
    log.msg('User Management page')
    log.msg('http://localhost:' + str(dsage_server.web_port) + '/' + dsage_server.secure_url)
    log.msg(DELIMITER)

application = service.Application('DSage Server')
dsage_service = service.MultiService()

dsage_server, dsage_server_port, ssl = start_dsage_server(dsage_service)
web_server_port, secure_url = start_web_server(dsage_server, dsage_service)
dsage_service.setServiceParent(application)

print_info(dsage_server)"""

WORKER_TAC = """"""