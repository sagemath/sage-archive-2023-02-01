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

DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
PUBKEY_DATABASE = os.path.abspath(os.path.join(DSAGE_DIR,
                                  'authorized_keys.db'))
DB_DIR = os.path.join(DSAGE_DIR, 'db/')
SAGE_ROOT = os.getenv('SAGE_ROOT')

def check_dsage_dir():
    if os.path.exists(DSAGE_DIR):
        return
    else:
        print "Creating " + DSAGE_DIR
        os.mkdir(DSAGE_DIR)

def get_config(type):
    config = ConfigParser.ConfigParser()
    config.add_section('general')
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
    print "Generating public/private key pair for authentication..."
    print "Your key will be stored in ${DOT_SAGE}/dsage/dsage_key"
    print "Just hit enter when prompted for a passphrase"
    key_file = os.path.join(DSAGE_DIR, 'dsage_key')
    cmd = ["ssh-keygen", "-q", "-trsa", "-f%s" % key_file]
    p = subprocess.call(cmd)
    print "\n"
    conf_file = os.path.join(DSAGE_DIR, 'client.conf')
    config.set('auth', 'privkey_file', key_file)
    config.set('auth', 'pubkey_file', key_file + '.pub')
    config.write(open(conf_file, 'w'))
    print "Client configuration finished."

def setup_worker():
    check_dsage_dir()
     # Get ConfigParser object
    config = get_config('worker')

    config.set('general', 'server', 'localhost')
    config.set('general', 'port', 8082)
    config.set('general', 'nice_level', 20)
    config.set('general', 'workers', 2)
    config.set('uuid', 'id', '')
    config.set('ssl', 'ssl', 1)
    config.set('log', 'log_file', 'stdout')
    config.set('log', 'log_level', '0')
    config.set('general', 'delay', '5')
    conf_file = os.path.join(DSAGE_DIR, 'worker.conf')
    config.write(open(conf_file, 'w'))
    print "Worker configuration finished."

def setup_server():
    check_dsage_dir()
    # Get ConfigParser object
    config = get_config('server')
    config.set('server', 'client_port', 8081)
    config.set('server', 'worker_port', 8082)
    config.set('ssl', 'ssl', 1)
    config.set('server_log', 'log_file', 'stdout')
    config.set('server_log', 'log_level', '0')
    config.set('db_log', 'log_file', 'stdout')
    config.set('db_log', 'log_level', '0')
    config.set('auth', 'pubkey_database', PUBKEY_DATABASE)
    config.set('db', 'db_file', os.path.join(DB_DIR, 'jobdb.fs'))
    config.set('db', 'prune_in_days', 7)
    config.set('db', 'stale_in_days', 365)
    config.set('db', 'failure_threshhold', 5)
    config.set('ssl', 'privkey_file', os.path.join(DSAGE_DIR, 'cacert.pem'))
    config.set('ssl', 'cert_file', os.path.join(DSAGE_DIR, 'privkey.pem'))


    print "Generating SSL certificate for server..."
    # creates a private key
    cmd = ['openssl genrsa > %s' % (config.get('ssl', 'privkey_file'))]
    subprocess.call(cmd, shell=True)
    # creates a self signed SSL certificate
    cmd = ['openssl req  -config %s -new -x509 -key %s -out %s -days \
           1000' % (os.path.join(SAGE_ROOT,'local/openssl/openssl.cnf'),
                    config.get('ssl', 'privkey_file'),
                    config.get('ssl', 'cert_file'))]
    subprocess.call(cmd, shell=True)

    # add default user
    c = ConfigParser.ConfigParser()
    c.read(os.path.join(DSAGE_DIR, 'client.conf'))
    username = c.get('auth', 'username')
    pubkey_file = c.get('auth', 'pubkey_file')
    add_user(username, pubkey_file)

    conf_file = os.path.join(DSAGE_DIR, 'server.conf')
    config.write(open(conf_file, 'w'))

    print "Server configuration finished."
    print "\n"

def add_user(username, pubkey_file):
    f = open(pubkey_file)
    type, key = f.readlines()[0].split()[:2]
    f.close()
    f1 = open(PUBKEY_DATABASE, 'a')
    f1.write(':'.join([username, key]) + '\n')
    f1.close()

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

