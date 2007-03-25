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

from sage.dsage.server.server import DSageServer
from sage.dsage.database.jobdb import JobDatabaseSQLite
from sage.dsage.database.monitordb import MonitorDatabase

def main():
    self.jobdb = JobDatabaseSQLite()
    self.monitordb = MonitorDatabase()
    self.dsage_server = DSageServer(self.jobdb, self.monitordb)

    return self.dsage_server.generate_xml_stats()

if __name__ == '__main__':
    main()