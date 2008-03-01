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
##############################################################################

from xml.etree.ElementTree import (ElementTree as ET,
                                   Element,
                                   SubElement)
from cStringIO import StringIO

class XMLStats(object):
    """
    Generates a XML document containing statistics for the server.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.workerdb = self.dsage_server.workerdb
        self.root = Element('stats')

    def gen_xml(self):
        self.add_workers()
        self.add_mhz()

        xml_stream = StringIO()
        ET(self.root).write(xml_stream)

        return xml_stream.getvalue()

    def add_element(self, element, value):
        SubElement(self.root, element).text = str(value)

    def add_workers(self):
        """
        Adds workers to the root.

        """

        busy_workers = self.workerdb.get_worker_count(connected=True,
                                                       busy=True)
        free_workers = self.workerdb.get_worker_count(connected=True,
                                                       busy=False)

        self.add_element('busy_workers', busy_workers)
        self.add_element('free_workers', free_workers)

    def add_mhz(self):
        """
        Adds the current working MHz to root."""

        working_mhz = self.workerdb.get_cpu_speed(connected=True,
                                                   busy= True)
        total_mhz = self.workerdb.get_cpu_speed(connected=True,
                                                 busy=False)
        online_mhz = self.workerdb.get_cpu_speed(busy=False, connected=True)
        self.add_element('working_mhz', working_mhz)
        self.add_element('total_mhz', total_mhz)