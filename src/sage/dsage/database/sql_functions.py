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

from twisted.python import log

def table_exists(con, tablename):
    """
    Check if a given table exists.
    If the below query is not None, then the table exists

    """

    query = """SELECT name FROM sqlite_master
    WHERE type = 'table' AND name = ?;
    """

    cur = con.cursor()
    cur.execute(query, (tablename,))
    result = cur.fetchone()
    return result

def create_table(con, tablename, query):
    r"""
    Creates a table given the connection.

    """

    log.msg('Creating table %s...' % tablename)
    con.execute(query)

def add_trigger(con, trigger):
    con.execute(trigger)