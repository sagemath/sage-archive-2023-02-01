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

from twisted.spread import pb

class NoJobException(Exception):
    def __init__(self):
        return

    def __str__(self):
        return "No jobs received from server."

class ErrorInResultException(Exception):
    def __init__(self):
        return

    def __str__(self):
        return "Caught exception/error in sage output."

class NotConnectedException(Exception):
    def __init__(self):
        return

    def __str__(self):
        return "Not connected to a remote server."

class BadTypeError(pb.Error):
    pass

class BadJobError(pb.Error):
    pass

class BadUserNameException(pb.Error):
    pass

class BadSignatureException(pb.Error):
    pass

class BadKeyException(pb.Error):
    pass
