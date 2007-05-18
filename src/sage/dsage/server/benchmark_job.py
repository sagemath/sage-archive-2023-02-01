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

from sage.dsage.database.job import Job

class BenchmarkJob(object):
    """
    This job is sent to workers as a way to benchmark their performance.

    """

    def __init__(self):
        pass

    def get_job(self):
        job = Job()
        job.code = """import time
start = time.time()
os.system("openssl speed rsa1024")
end = time.time()
DSAGE_RESULT = end - start
"""

        return job