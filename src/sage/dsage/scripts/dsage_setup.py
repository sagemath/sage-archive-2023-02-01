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
import random
import socket
import ConfigParser
import subprocess
import sys
import sqlite3

from sage.dsage.database.clientdb import ClientDatabaseSA as ClientDatabase
from sage.dsage.database.db_config import create_schema
from sage.dsage.misc.constants import (DELIMITER, DSAGE_DIR, DSAGE_DB_DIR,
                                       DSAGE_DB)
from sage.dsage.misc.config import check_dsage_dir
from sage.dsage.__version__ import version

from sage.misc.viewer import cmd_exists

DB_DIR = os.path.join(DSAGE_DIR, 'db/')
SAGE_ROOT = os.getenv('SAGE_ROOT')
DSAGE_VERSION = version

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

def add_default_client(Session):
    """
    Adds the default client.

    """

    from twisted.conch.ssh import keys
    from getpass import getuser

    clientdb = ClientDatabase(Session)

    username = getuser()
    pubkey_file = os.path.join(DSAGE_DIR, 'dsage_key.pub')
    pubkey = keys.Key.fromFile(pubkey_file)
    pubkey_str = pubkey.toString(type='openssh')
    if clientdb.get_client(username) is None:
        clientdb.add_client(username, pubkey_str)
        print 'Added user %s.\n' % (username)
    else:
        client = clientdb.get_client(username)
        if client.public_key != pubkey_str:
            clientdb.del_client(username)
            clientdb.add_client(username, pubkey_str)
            print "User %s's pubkey changed, setting to new one." % (username)
        else:
            print 'User %s already exists.' % (username)

def setup_client(testing=False):
    check_dsage_dir()
    key_file = os.path.join(DSAGE_DIR, 'dsage_key')
    if testing:
        cmd = ["ssh-keygen", "-q", "-trsa", "-P ''", "-f%s" % key_file]
        return

    if not cmd_exists('ssh-keygen'):
        print DELIMITER
        print "Could NOT find ssh-keygen."
        print "Aborting."
        return

    print DELIMITER
    print "Generating public/private key pair for authentication..."
    print "Your key will be stored in %s/dsage_key" % DSAGE_DIR
    print "Just hit enter when prompted for a passphrase"
    print DELIMITER

    cmd = ["ssh-keygen", "-q", "-trsa", "-f%s" % key_file]
    ld = os.environ['LD_LIBRARY_PATH']
    try:
        del os.environ['LD_LIBRARY_PATH']
        p = subprocess.call(cmd)
    finally:
        os.environ['LD_LIBRARY_PATH'] = ld

    print "\n"
    print "Client configuration finished.\n"

def setup_worker():
    check_dsage_dir()
    print "Worker configuration finished.\n"

def setup_server(template=None):
    check_dsage_dir()
    print "Choose a domain name for your SAGE notebook server,"
    print "for example, localhost (personal use) or %s (to allow outside connections)." % socket.getfqdn()
    dn = raw_input("Domain name [localhost]: ").strip()
    if dn == '':
        print "Using default localhost"
        dn = 'localhost'

    template_dict = {'organization': 'SAGE (at %s)' % (dn),
                'unit': '389',
                'locality': None,
                'state': 'Washington',
                'country': 'US',
                'cn': dn,
                'uid': 'sage_user',
                'dn_oid': None,
                'serial': str(random.randint(1,2**31)),
                'dns_name': None,
                'crl_dist_points': None,
                'ip_address': None,
                'expiration_days': 10000,
                'email': 'sage@sagemath.org',
                'ca': None,
                'tls_www_client': None,
                'tls_www_server': True,
                'signing_key': True,
                'encryption_key': True,
                }

    if isinstance(template, dict):
        template_dict.update(template)

    s = ""
    for key, val in template_dict.iteritems():
        if val is None:
            continue
        if val == True:
            w = ''
        elif isinstance(val, list):
            w = ' '.join(['"%s"' % x for x in val])
        else:
            w = '"%s"' % val
        s += '%s = %s \n' % (key, w)

    template_file = os.path.join(DSAGE_DIR, 'cert.cfg')
    f = open(template_file, 'w')
    f.write(s)
    f.close()

    # Disable certificate generation -- not used right now anyways
    privkey_file = os.path.join(DSAGE_DIR, 'cacert.pem')
    pubkey_file = os.path.join(DSAGE_DIR, 'pubcert.pem')

    print DELIMITER
    print "Generating SSL certificate for server..."

    if False and os.uname()[0] != 'Darwin' and cmd_exists('openssl'):
        # We use openssl by default if it exists, since it is *vastly*
        # faster on Linux.
        cmd = ['openssl genrsa > %s' % privkey_file]
        print "Using openssl to generate key"
        print cmd[0]
        subprocess.call(cmd, shell=True)
    else:
        cmd = ['certtool --generate-privkey --outfile %s' % privkey_file]
        print "Using certtool to generate key"
        print cmd[0]
        # cmd = ['openssl genrsa > %s' % privkey_file]
        subprocess.call(cmd, shell=True)

    cmd = ['certtool --generate-self-signed --template %s --load-privkey %s \
           --outfile %s' % (template_file, privkey_file, pubkey_file)]
    subprocess.call(cmd, shell=True)
    print DELIMITER

    # Set read only permissions on cert
    os.chmod(os.path.join(DSAGE_DIR, 'cacert.pem'), 0600)

    # create database schemas
    from sage.dsage.database.db_config import init_db_sa as init_db
    Session = init_db(DSAGE_DB)

    # add default user
    add_default_client(Session)

    print "Server configuration finished.\n\n"

def setup(template=None):
    setup_client()
    setup_worker()
    setup_server(template=template)
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

