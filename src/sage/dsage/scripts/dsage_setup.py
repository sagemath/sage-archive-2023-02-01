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

import os
import ConfigParser
import subprocess
import sys

from sage.dsage.database.clientdb import ClientDatabase
from sage.dsage.twisted.pubkeyauth import get_pubkey_string
from sage.dsage.misc.constants import delimiter as DELIMITER
from sage.dsage.__version__ import version

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
DB_DIR = os.path.join(DSAGE_DIR, 'db/')
SAGE_ROOT = os.getenv('SAGE_ROOT')
DSAGE_VERSION = version

def check_dsage_dir():
    if os.path.exists(DSAGE_DIR):
        return
    else:
        print "Creating " + DSAGE_DIR
        os.mkdir(DSAGE_DIR)

def get_config(type):
    config = ConfigParser.ConfigParser()
    config.add_section('general')
    config.set('general', 'version', DSAGE_VERSION)
    config.add_section('ssl')
    if type == 'client':
        config.add_section('auth')
        config.add_section('log')
    elif type == 'worker':
        config.add_section('uuid')
        config.add_section('log')
    elif type == 'server':
        config.add_section('auth')
        config.add_section('server')
        config.add_section('server_log')
        config.add_section('db')
        config.add_section('db_log')
    return config

def setup_client():
    check_dsage_dir()
    # Get ConfigParser object
    config = get_config('client')

    config.set('auth', 'username', os.getenv('USER'))
    config.set('general', 'server', 'localhost')
    config.set('general', 'port', 8081)
    config.set('ssl', 'ssl', 1)
    config.set('log', 'log_file', 'stdout')
    config.set('log', 'log_level', '0')
    # set public key authentication info
    print DELIMITER
    print "Generating public/private key pair for authentication..."
    print "Your key will be stored in %s/dsage_key"%DSAGE_DIR
    print "Just hit enter when prompted for a passphrase"
    print DELIMITER
    key_file = os.path.join(DSAGE_DIR, 'dsage_key')
    cmd = ["ssh-keygen", "-q", "-trsa", "-f%s" % key_file]
    p = subprocess.call(cmd)
    print "\n"
    conf_file = os.path.join(DSAGE_DIR, 'client.conf')
    config.set('auth', 'privkey_file', key_file)
    config.set('auth', 'pubkey_file', key_file + '.pub')
    config.write(open(conf_file, 'w'))
    print "Client configuration finished.\n"

def setup_worker():
    check_dsage_dir()
     # Get ConfigParser object
    config = get_config('worker')

    config.set('general', 'server', 'localhost')
    config.set('general', 'port', 8081)
    config.set('general', 'nice_level', 20)
    config.set('general', 'workers', 2)
    config.set('uuid', 'id', '')
    config.set('ssl', 'ssl', 1)
    config.set('log', 'log_file', 'stdout')
    config.set('log', 'log_level', '0')
    config.set('general', 'delay', '5')
    config.set('general', 'anonymous', False)
    conf_file = os.path.join(DSAGE_DIR, 'worker.conf')
    config.write(open(conf_file, 'w'))
    print "Worker configuration finished.\n"

def setup_server():
    check_dsage_dir()
    # Get ConfigParser object
    config = get_config('server')
    config.set('server', 'client_port', 8081)
    config.set('ssl', 'ssl', 1)
    config.set('server_log', 'log_file', 'stdout')
    config.set('server_log', 'log_level', '0')
    config.set('db_log', 'log_file', 'stdout')
    config.set('db_log', 'log_level', '0')
    config.set('auth', 'pubkey_database', os.path.join(DB_DIR, 'dsage.db'))
    config.set('db', 'db_file', os.path.join(DB_DIR, 'dsage.db'))
    config.set('db', 'prune_in_days', 7)
    config.set('db', 'stale_in_days', 365)
    config.set('db', 'job_failure_threshold', 2)
    config.set('ssl', 'privkey_file', os.path.join(DSAGE_DIR, 'cacert.pem'))
    config.set('ssl', 'cert_file', os.path.join(DSAGE_DIR, 'privkey.pem'))
    config.set('general', 'stats_file', 'gauge.xml')

    print DELIMITER
    print "Generating SSL certificate for server..."
    cmd = ['openssl genrsa > %s' % (config.get('ssl', 'privkey_file'))]
    subprocess.call(cmd, shell=True)
    cmd = ['openssl req  -config %s -new -x509 -key %s -out %s -days \
           1000' % (os.path.join(SAGE_ROOT,'local/openssl/openssl.cnf'),
                    config.get('ssl', 'privkey_file'),
                    config.get('ssl', 'cert_file'))]
    subprocess.call(cmd, shell=True)
    print DELIMITER

    conf_file = os.path.join(DSAGE_DIR, 'server.conf')
    config.write(open(conf_file, 'w'))

    print "Server configuration finished.\n\n"

    # add default user
    c = ConfigParser.ConfigParser()
    c.read(os.path.join(DSAGE_DIR, 'client.conf'))
    username = c.get('auth', 'username')
    pubkey_file = c.get('auth', 'pubkey_file')
    clientdb = ClientDatabase()
    if clientdb.get_user(username) is None:
        clientdb.add_user(username, pubkey_file)
        print 'Added user %s.\n' % (username)
    else:
        user, key = clientdb.get_user_and_key(username)
        if key != get_pubkey_string(pubkey_file):
            clientdb.del_user(username)
            clientdb.add_user(username, pubkey_file)
            print "User %s's pubkey changed, setting to new one." % (username)
        else:
            print 'User %s already exists.' % (username)

def setup():
    setup_client()
    setup_worker()
    setup_server()
    print "Configuration finished.."

if __name__ == '__main__':
    if len(sys.argv) == 1:
        setup()
    if len(sys.argv) == 2:
        if sys.argv[1] == 'server':
            setup_server()
        elif sys.argv[1] == 'worker':
            setup_worker()
        elif sys.argv[1] == 'client':
            setup_client()

