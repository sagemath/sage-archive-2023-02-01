r"""
This file exposes most of the public facing \code{dsage} functionality.

"""
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

import os
from sage.dsage.dist_functions.all import *
from sage.dsage.misc.constants import DSAGE_DIR
from sage.dsage.dsage import DistributedSage
import sage.plot.plot

def connect(server='localhost', port=8081,
          username=os.getenv('USER'),
          pubkey_file=os.path.join(DSAGE_DIR,'dsage_key.pub'),
          privkey_file=os.path.join(DSAGE_DIR, 'dsage_key'),
          log_level=0, ssl=True, testing=False):
    """
    This object represents a connection to the distributed SAGE server.

    Parameters:
    server -- str (Default: 'localhost')
    port -- int (Default: 8081)
    username -- str
    pubkey_file -- str (Default: ~/.sage/dsage/dsage_key.pub)
    privkey_file -- str (Default: ~/.sage/dsage/dsage_key)
    log_level -- int (Default: 0)
    ssl -- int (Default: 1)

    """

    from sage.dsage.interface.dsage_interface import BlockingDSage

    if sage.plot.plot.DOCTEST_MODE:
        testing = True

    return BlockingDSage(server=server, port=port, username=username,
                         pubkey_file=pubkey_file, privkey_file=privkey_file,
                         log_level=log_level, ssl=ssl, testing=testing)

dsage = DistributedSage()
dsage.connect = connect