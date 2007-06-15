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
from sage.dsage.misc.constants import DELIMITER as DELIMITER
from sage.dsage.misc.constants import DSAGE_DIR
from sage.dsage.misc.config import check_dsage_dir
from sage.dsage.__version__ import version

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

def setup_client(testing=False):
    check_dsage_dir()
    key_file = os.path.join(DSAGE_DIR, 'dsage_key')
    if testing:
        cmd = ["ssh-keygen", "-q", "-trsa", "-P ''", "-f%s" % key_file]
        return

    print DELIMITER
    print "Generating public/private key pair for authentication..."
    print "Your key will be stored in %s/dsage_key"%DSAGE_DIR
    print "Just hit enter when prompted for a passphrase"
    print DELIMITER
    cmd = ["ssh-keygen", "-q", "-trsa", "-f%s" % key_file]
    p = subprocess.call(cmd)
    print "\n"
    print "Client configuration finished.\n"

def setup_worker():
    check_dsage_dir()
    print "Worker configuration finished.\n"

def setup_server(template_dict=None):
    check_dsage_dir()
    template_file = 'cert.cfg'
    template = {'organization': 'SAGE',
                'unit': '389',
                'locality': None,
                'state': 'Washington',
                'country': 'US',
                'cn': 'SAGE User',
                'uid': 'sage_user',
                'dn_oid': None,
                'serial': 007,
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
    if isinstance(template_dict, dict):
        template.update(template_dict)

    s = ""
    for key, val in template.iteritems():
        if val is None:
            continue
        if val == True:
            w = ''
        elif isinstance(val, list):
            w = ' '.join(['"%s"'%x for x in val])
        else:
            w = '"%s"'%val
        s += '%s = %s \n'%(key, w)
    f = open(os.path.join(DSAGE_DIR, template_file), 'w')
    f.write(s)
    f.close()
    # Disable certificate generation -- not used right now anyways
    privkey_file = os.path.join(DSAGE_DIR, 'cacert.pem')
    pubkey_file = os.path.join(DSAGE_DIR, 'pubcert.pem')
    print DELIMITER
    print "Generating SSL certificate for server..."
    cmd = ['certtool --generate-privkey --outfile %s' % privkey_file]
    # cmd = ['openssl genrsa > %s' % privkey_file]
    subprocess.call(cmd, shell=True)
    # cmd = ['openssl req  -config %s -new -x509 -key %s -out %s -days \
    #        1000' % (os.path.join(SAGE_ROOT,'local/openssl/openssl.cnf'),
    #                 privkey_file, pubkey_file)]
    cmd = ['certtool --generate-self-signed --template %s --load-privkey %s \
           --outfile %s' % (template_file, privkey_file, pubkey_file)]
    subprocess.call(cmd, shell=True)
    print DELIMITER
    os.chmod(os.path.join(DSAGE_DIR, 'cacert.pem'), 0600)

    # conf_file = os.path.join(DSAGE_DIR, 'server.conf')
    # config.write(open(conf_file, 'w'))

    # add default user
    from twisted.conch.ssh import keys
    import base64

    # c = ConfigParser.ConfigParser()
    # c.read(os.path.join(DSAGE_DIR, 'client.conf'))
    from getpass import getuser
    username = getuser()
    pubkey_file = os.path.join(DSAGE_DIR, 'dsage_key.pub')
    clientdb = ClientDatabase()
    pubkey = base64.encodestring(
                    keys.getPublicKeyString(filename=pubkey_file).strip())
    if clientdb.get_user(username) is None:
        clientdb.add_user(username, pubkey)
        print 'Added user %s.\n' % (username)
    else:
        user, key = clientdb.get_user_and_key(username)
        if key != pubkey:
            clientdb.del_user(username)
            clientdb.add_user(username, pubkey)
            print "User %s's pubkey changed, setting to new one." % (username)
        else:
            print 'User %s already exists.' % (username)

    print "Server configuration finished.\n\n"

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

