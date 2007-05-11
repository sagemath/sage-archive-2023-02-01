##############################################################################
#
#  DSAGE: Distributed SAGE
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
#
##############################################################################

"""
Gets the different configuration options for the various components of DSAGE.

"""

import os
import ConfigParser
import uuid

from sage.dsage.misc.misc import random_str

def check_version(old_version):
    from sage.dsage.__version__ import version
    if version != old_version:
        msg = "Version mismatch:\nHas:\t%s\nNeeds:\t%s" % (old_version,
                                                           version)
        raise ValueError(msg)

def read_conf(config):
    conf = {}
    for sec in config.sections():
        conf.update(dict(config.items(sec)))
    try:
        check_version(conf['version'])
    except ValueError, msg:
        raise ValueError(msg)
    try:
        conf['log_level'] = int(conf['log_level'])
    except Exception, msg:
        pass
    try:
        conf['workers'] = int(conf['workers'])
    except Exception, msg:
        pass
    try:
        conf['delay'] = float(conf['delay'])
    except Exception, msg:
        pass
    try:
        conf['priority'] = int(conf['priority'])
    except Exception, msg:
        pass

    return conf

def get_conf(type, test=False):
    if test:
        DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'tmp')
    else:
        DSAGE_DIR = os.path.join(os.getenv('DOT_SAGE'), 'dsage')
    config = ConfigParser.ConfigParser()
    try:
        if type == 'client':
            conf_file = os.path.join(DSAGE_DIR, 'client.conf')
            DATA =  random_str(length=500)
            config.read(conf_file)
            conf = read_conf(config)
            conf['data'] = DATA
        elif type == 'server':
            conf_file = os.path.join(DSAGE_DIR, 'server.conf')
            config.read(conf_file)
            conf = read_conf(config)
        elif type == 'monitor':
            conf_file = os.path.join(DSAGE_DIR, 'worker.conf')
            config.read(conf_file)
            if len(config.get('uuid', 'id')) != 36:
                config.set('uuid', 'id', str(uuid.uuid1()))
                f = open(conf_file, 'w')
                config.write(f)
                config.read(conf_file)
            conf = read_conf(config)
        elif type == 'jobdb' or type == 'clientdb' or type == 'monitordb':
            conf_file = os.path.join(DSAGE_DIR, 'server.conf')
            config.read(conf_file)
            conf = read_conf(config)
            conf['log_level'] = config.get('db_log', 'log_level')
            conf['log_file'] = config.get('db_log', 'log_file')

        conf['conf_file'] = conf_file

        return conf
    except Exception, msg:
        print msg
        print "Error reading '%s', run dsage.setup()" % conf_file

def get_bool(value):
    boolean_states = {'0': False,
    '1': True,
    'false': False,
    'no': False,
    'off': False,
    'on': True,
    'true': True,
    'yes': True}
    if value.lower() not in boolean_states:
        raise ValueError('Not a boolean: %s' % value)

    return boolean_states[value.lower()]