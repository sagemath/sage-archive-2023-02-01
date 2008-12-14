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
Miscellaneous helper methods.

"""

import sys
import random
import os
import string

from sage.dsage.misc.constants import DSAGE_DIR
from sage.dsage.misc.config import check_dsage_dir

def exec_wrs(script):
    """
    Executes a script which reports a rank as an integer.

    """

    from cStringIO import StringIO
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        exec(script)
        return int(sys.stdout.getvalue())
    finally:
        sys.stdout = old_stdout

def random_string(length):
    """
    Returns a random string

    Parameters:
        length -- the length of the string

    """

    random.seed()
    l = list((length*length) * (string.letters + string.digits))
    random.shuffle(l)
    s = ''.join(random.sample(l, length))

    return s

def get_uuid():
    check_dsage_dir()
    try:
        uuid = open(os.path.join(DSAGE_DIR,'uuid')).read()
    except (OSError, IOError):
        uuid = ""

    return uuid

def set_uuid(uuid):
    from sage.dsage.misc.config import check_dsage_dir
    check_dsage_dir()
    f = open(os.path.join(DSAGE_DIR, 'uuid'), 'w')
    f.write(uuid)
    f.close()

def gen_uuid():
    import uuid
    return str(uuid.uuid1())

def check_uuid(uuid):
    if not isinstance(uuid, str):
        return False
    elif not len(uuid) ==36:
        return False
    else:
        return True

def random_str(length=500):
    """
    Generates a random string.

    INPUT:
    length -- the length of the string

    """

    r_str = [chr(i) for i in [random.randint(65, 123) for n in range(length)]]

    return ''.join(r_str)

def timedelta_to_seconds(time_delta):
    """
    Converts a timedelta object into seconds.

    """

    days, seconds, microseconds = (time_delta.days,
                                   time_delta.seconds,
                                   time_delta.microseconds)

    seconds = float(days*24*60*60 + seconds + (microseconds/10.0**6))

    return seconds

def find_open_port(server='localhost', low=0):
    """
    Tries to find an open port on your machine to use.

    """
    if low == 0:
        low = 8081+(os.getpid()%2609)
    import socket

    port = low
    while(True):
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect((server, port))
            s.close()
            port += 1
        except socket.error, msg:
            if msg[1] == 'Connection refused':
                yield port
                port += 1
            else:
                port += 1
