#!/usr/bin/env python
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

import sys
import os
from twisted.spread import pb
from twisted.internet import     utils, defer
from twisted.python import log

class HostInfo(pb.Copyable, pb.RemoteCopy):
    """
    Class to gather computer specifications on the running host.

    """

    def __str__(self):
        return str(self.host_info)

    def __repr__(self):
        return str(self.host_info)

    def _catch_failure(self, failure):
        log.msg("Error: ", failure.getErrorMessage())
        log.msg("Traceback: ", failure.printTraceback())

    def _gotsysctl(self, output, host_info):
        d = defer.Deferred()
        output = output.split('\n')
        for line in output:
            if line == '':
                continue
            l = line.strip()
            if '=' in l:
                l = l.split('=')
            elif ':' in line:
                l = l.split(':')

            l = [li.strip() for li in l]
            if l[0] == 'hw.cpufrequency':
                host_info['cpu Mhz'] = (int(l[1]) / 1000000.0)
            elif l[0] == 'hw.availcpu':
                host_info['cpus'] = int(l[1])
            elif l[0] == 'hw.physmem':
                host_info['MemTotal'] = int(l[1])
            elif l[0] == 'hw.usermem':
                host_info['MemFree'] = (host_info['MemTotal'] -
                                           int(l[1])) / 1024 / 1024
            elif l[0] == 'machdep.cpu.brand_string':
                host_info['cpu model'] = l[1]

        self.canonical_info(host_info)
        d.callback(self.host_info)
        return d

    def get_host_info(self):
        platform = sys.platform
        host_info = {}
        if platform == 'linux2' or platform == 'linux':
            d = defer.Deferred()
            host_info['os'] = platform
            cpuinfo = open('/proc/cpuinfo','r').readlines()

            cpus = 0
            for line in cpuinfo:
                line = line.replace('\n','')
                line = line.replace('\t','')
                if line == '':
                    continue
                if 'processor' in line:
                    cpus += 1
                s = line.split(':')
                if s != ['\n']:
                    host_info[s[0].strip()] = s[1].strip()

            host_info['cpus'] = cpus

            # uptime
            uptime = open('/proc/uptime', 'r').readline().split(' ')
            host_info['uptime'] = uptime[0]

            # memory info
            meminfo = open('/proc/meminfo', 'r').readlines()
            for line in meminfo:
                s = line.split(':')
                if s != ['\n']:
                    host_info[s[0].strip()] = s[1].strip()

            # hostname
            hostname = os.popen('hostname').readline().strip()
            host_info['hostname'] = hostname

            # kernel version
            kernel_version = os.popen('uname -r').readline().strip()
            host_info['kernel_version'] = kernel_version

            host_info['os'] = platform
            d.callback(self.canonical_info(host_info))
            return d

        if platform == 'darwin':
            cmd = '/usr/sbin/sysctl'
            args = ('-a', 'hw')
            d = utils.getProcessOutput(cmd, args, errortoo=1)
            d.addCallback(self._gotsysctl, host_info)
            d.addErrback(self._catch_failure)
            return d

    def canonical_info(self, platform_host_info):
        """Standarize host info so we can parse it easily"""

        unify_info = {'model name': 'cpu_model',
                          'cpu Mhz': 'cpu_speed',
                          'MemTotal': 'mem_total',
                          'MemFree': 'mem_free',
                          'kernel_version': 'os',
                          'cpu family': 'cpu_model',
                          'cpu model': 'cpu_model',
                          'cache size': 'cpu_cache_size',
                          'fpu': 'fpu',
                          'hostname': 'hostname',
                          'cpus': 'cpus',
                          'ip': 'ip',
                          'os': 'os'
                     }

        canonical_info = {}
        for k,v in platform_host_info.iteritems():
            try:
                canonical_info[unify_info[k]] = v
            except KeyError:
                pass

        try:
            import sage.version
            canonical_info['sage_version'] = sage.version.version
        except ImportError:
            canonical_info['sage_version'] = 'unknown'

        self.host_info = canonical_info

        return self.host_info

class ClassicHostInfo(object):
    """
    Class to gather computer specifications on the running host.

    """

    def __init__(self):
        self.host_info = self.get_host_info()

    def __str__(self):
        return str(self.host_info)

    def __repr__(self):
        return str(self.host_info)

    def get_host_info(self):
        import platform
        host_info = {}
        #  uname = (system,node,release,version,machine,processor)
        uname = platform.uname()
        host_info['os'] = uname[0]

        if host_info['os'] == 'Linux':
            # Get CPU related data
            cpuinfo = open('/proc/cpuinfo','r').readlines()
            cpus = 0
            for line in cpuinfo:
                if 'processor' in line:
                    cpus += 1
                s = line.split(':')
                if s != ['\n']:
                    key = s[0].strip()
                    value = s[1].strip()
                    if key == 'cpu MHz':
                        host_info[key] = int(float(value))
                    else:
                        host_info[key] = value
            host_info['cpus'] = cpus

            # perform architecture specific modifications of host_info
            if uname[5] == 'ppc':
                host_info['clock'] = host_info['clock'].replace('MHz', '')
                host_info['cpu MHz'] = int(float(host_info['clock']))
                host_info['model name'] = host_info['cpu']
            elif uname[5] == 'ia64':
                host_info['model name'] = host_info['family']

            # Get memory related date
            meminfo = open('/proc/meminfo', 'r').readlines()
            for line in meminfo:
                s = line.split(':')
                if s != ['\n']:
                    key = s[0].strip()
                    value = s[1]
                    if key == 'MemTotal':
                        mem_total = int(int(value.split()[0].strip()) / 1024)
                        host_info[key] = mem_total
                    elif key == 'MemFree':
                        mem_free = int(int(value.split()[0].strip()) / 1024)
                        host_info[key] = mem_free
                    else:
                        host_info[key] = value.strip()
        elif host_info['os'] == 'Darwin':
            for line in os.popen('sysctl -a hw machdep').readlines():
                l = line.strip()
                if '=' in l:
                    l = l.split('=')
                if ':' in line:
                    l = l.split(':')

                l = [li.strip() for li in l]
                if l[0] == 'hw.cpufrequency':
                    host_info['cpu MHz'] = str(int(l[1]) / 1000000)
                elif l[0] == 'hw.availcpu':
                    host_info['cpus'] = int(l[1])
                elif l[0] == 'hw.physmem':
                    host_info['MemTotal'] = int(int(l[1]) / (1024*1024))
                elif l[0] == 'hw.usermem':
                    mem_total = int(host_info['MemTotal'])
                    user_mem = int(l[1]) / (1024*1024)
                    mem_free = int(mem_total - user_mem)
                    host_info['MemFree'] = mem_free
                elif l[0] == 'machdep.cpu.brand_string':
                    host_info['model name'] = l[1]
                elif l[0] == 'hw.model': # OS X PPC
                    host_info['model name'] = l[1]

        host_info['hostname'] = os.uname()[1]
        host_info['kernel_version'] = os.uname()[2]
        host_info['uptime'] = uptime = os.popen('uptime').readline().strip()

        return self.canonical_info(host_info)

    def canonical_info(self, platform_host_info):
        """
        Standarize host info so we can parse it easily.

        """

        unify_info = {'model name': 'cpu_model',
                      'cpu MHz': 'cpu_speed',
                      'MemTotal': 'mem_total',
                      'MemFree': 'mem_free',
                      'kernel_version': 'kernel_version',
                      'cache size': 'cpu_cache_size',
                      'fpu': 'fpu',
                      'hostname': 'hostname',
                      'cpus': 'cpus',
                      'ip': 'ip',
                      'os': 'os',
                      'uptime': 'uptime'}
        canonical_info = {}
        for k,v in platform_host_info.iteritems():
            try:
                canonical_info[unify_info[k]] = v
            except KeyError:
                pass
        try:
            import sage.version
            canonical_info['sage_version'] = sage.version.version
        except ImportError:
            canonical_info['sage_version'] = 'unknown'

        return canonical_info

if __name__ == '__main__':
    h = ClassicHostInfo()
    print h