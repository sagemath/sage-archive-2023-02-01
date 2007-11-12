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
                                   SubElement,
                                   dump,
                                   XML)

class XMLStats(self, dsage_server):
    """
    Generates a XML document containing statistics for the server.

    """

    def __init__(self, dsage_server):
        self.dsage_server = dsage_server
        self.monitordb = self.dsage_server.monitordb
        self.root = Element('stats')

    def gen_xml(self):
        self.add_workers()
        self.add_working_mhz()

        return dump(self.root)

    def add_workers(self):
        """
        Adds workers to the root.

        """

        busy_workers = self.monitordb.get_worker_count(connected=True,
                                                       busy=True)
        free_workers = self.monitordb.get_worker_count(connected=True,
                                                       busy=False)
        offline_workers = self.monitordb.get_worker_count(connected=False)
        total_workers = (working_workers +
                         free_workers +
                         disconnected_workers)

        bw = SubElement(self.root, 'busy_workers')
        bw.text = busy_workers
        fw = SubElement(self.root, 'free_workers')
        fw.text = free_workers
        ow = SubElement(self.root, 'offline_workers')
        ow.text = offline_workers
        tw = SubElement(self.root, 'total_workers')
        tw.text = total_workers

    def add_working_mhz(self):
        """
        Adds the current working MHz to root."""

        working_mhz = self.monitordb.get_cpu_speed(connected=True,
                                                   busy=True)
        wm = SubElement(self.root, 'working_mhz')
        wm.text = working_mhz

# def generate_xml_stats(self):
#     """
#     This method returns a an XML document to be consumed by the
#     Mac OS X Dashboard widget and the web server.
#
#     """
#
#     def add_onlineProcessorCount(doc, gauge):
#         onlineProcessorCount = doc.createElement('onlineProcessorCount')
#         gauge.appendChild(onlineProcessorCount)
#         cpu_count = self.monitordb.get_cpu_count(connected=True)
#         c = doc.createTextNode(str(cpu_count))
#         onlineProcessorCount.appendChild(c)
#
#         return doc, onlineProcessorCount
#
#     def add_offlineProcessorCount(doc, gauge):
#         offlineProcessorCount = doc.createElement('offlineProcessorCount')
#         gauge.appendChild(offlineProcessorCount)
#         cpu_count = self.monitordb.get_cpu_count(connected=False)
#         c = doc.createTextNode(str(cpu_count))
#         offlineProcessorCount.appendChild(c)
#
#         return doc, offlineProcessorCount
#
#     def add_workingProcessorCount(doc, gauge):
#         workingProcessorCount = doc.createElement('workingProcessorCount')
#         gauge.appendChild(workingProcessorCount)
#         worker_count = self.monitordb.get_cpu_count(connected=True)
#         pcount = doc.createTextNode(str(worker_count))
#         workingProcessorCount.appendChild(pcount)
#
#         return doc, workingProcessorCount
#
#     def add_workingAgentPercentage(doc, gauge):
#         workingAgentPercentage = doc.createElement(
#                                                 'workingAgentPercentage')
#         gauge.appendChild(workingAgentPercentage)
#         working_workers = self.monitordb.get_worker_count(connected=True,
#                                                           busy=True)
#         free_workers = self.monitordb.get_worker_count(connected=True,
#                                                        busy=False)
#         disconnected_workers = self.monitordb.get_worker_count(
#                                connected=False,
#                                busy=False)
#         total_workers = (working_workers +
#                          free_workers +
#                          disconnected_workers)
#
#         if total_workers != 0:
#             worker_percentage = (float(working_workers / total_workers) *
#             100)
#         else:
#             worker_percentage = 0.0
#         percentage = doc.createTextNode(str(worker_percentage))
#         workingAgentPercentage.appendChild(percentage)
#
#         return doc, workingAgentPercentage
#
#     def add_date(doc, gauge):
#         date = datetime.datetime.now()
#
#         year = doc.createElement('Year')
#         gauge.appendChild(year)
#         year.appendChild(doc.createTextNode(str(date.year)))
#
#         seconds = doc.createElement('Seconds')
#         gauge.appendChild(seconds)
#         seconds.appendChild(doc.createTextNode(str(date.second)))
#
#         minutes = doc.createElement('Minutes')
#         gauge.appendChild(minutes)
#         minutes.appendChild(doc.createTextNode(str(date.minute)))
#
#         return doc, year, seconds, minutes
#
#     doc = xml.dom.minidom.Document()
#     doc, gauge = create_gauge(doc)
#
#     add_onlineAgentCount(doc, gauge)
#     add_offlineAgentCount(doc, gauge)
#     add_availableAgentCount(doc, gauge)
#     add_unavailableAgentCount(doc, gauge)
#     add_totalAgentCount(doc, gauge)
#     add_workingAgentCount(doc, gauge)
#     add_workingAgentPercentage(doc, gauge)
#     add_onlineProcessorCount(doc, gauge)
#     add_availableProcessorCount(doc, gauge)
#     add_unavailableProcessorCount(doc, gauge)
#     add_workingProcessorCount(doc, gauge)
#     add_workingMegaHertz(doc, gauge)
#
#     add_date(doc, gauge)
#     s = cStringIO.StringIO()
#     doc.writexml(s, newl='\n')
#
#     return s.getvalue()