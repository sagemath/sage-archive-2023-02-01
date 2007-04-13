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

"""This is a dummy module that keeps a list of connected workers. """

worker_list = []

def add(worker):
    """
    Adds an avatar to worker_list.

    """

    worker_list.append(worker)

def remove(worker):
    """
    Removes the avatar from the worker_list.

    """

    worker_list.remove(worker)
