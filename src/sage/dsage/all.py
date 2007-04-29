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

from sage.dsage.dsage import dsage
from sage.dsage.dist_functions.all import *

def DSage(server=None, port=None, username=None, pubkey_file=None, privkey_file=None):
    from sage.dsage.interface.dsage_interface import BlockingDSage
    return BlockingDSage(server=server, port=port, username=username,
                         pubkey_file=pubkey_file, privkey_file=privkey_file)

