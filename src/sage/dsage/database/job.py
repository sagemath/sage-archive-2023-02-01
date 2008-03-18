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
#
############################################################################

import datetime
import cPickle
import zlib
import bz2
import os
import copy
from getpass import getuser


class Job(object):
    """
    Defines a Job that gets distributed to clients.

    """

    def __init__(self, job_id=None, name='Unamed', username=getuser(),
                 code='', timeout=0, kind='sage', priority=5):
        """
        Represents a job.

        :type job_id: string
        :param job_id: unique identifier for a job

        :type name: string
        :param name: name given to a job, not unique

        :type code: string
        :param code: the code that needs to be executed

        :type parent: string
        :param parent: the job_id of another job

        :type username: string
        :param username: username of person who created job

        :type timeout: integer
        :param timeout: upper bound for number of seconds this job takes

        :type priority: integer
        :param priority: a jobs priority from 1-5, 1 being the highest

        :type kind: string
        :param kind: kind of the job (file, string, generator)

        """

        self.job_id = job_id
        self.name = name
        self.username = username
        self.code = code
        self.timeout = timeout
        self.priority = 5
        self.kind = kind
        self.status = 'new'
        self.data = []
        self.output = ''
        self.result = None
        self.creation_time = datetime.datetime.now()
        self.start_time = None
        self.wall_time = None
        self.update_time = None
        self.finish_time = None
        self.killed = False
        self.failures = 0
        self.uuid = ''

    def __str__(self):
        return "<Job('%s', %s)>" % (self.job_id, self.username)

    def __repr__(self):
        return self.__str__()

    def attach(self, var, obj, file_name=None):
        """
        Attaches an object to a job.

        Parameters:
        var -- the variable name you'd like the worker to use
        obj -- the object you want to attach
        filename -- optional, if your object is a saved sobj

        """

        if file_name is not None:
            try:
                s = open(file_name, 'rb').read()
                s = zlib.compress(s)
            except:
                print 'Unable to load %s. ' % file_name
                return
        else:
            try:
                s = cPickle.dumps(obj, 2)
                s = zlib.compress(s)
            except cPickle.PicklingError:
                print 'Unable to attach your object.'
        self.data.append((var, s, 'object'))

    def attach_file(self, file_name):
        """
        Attach a file to a job.

        Parameters:
        file_name -- obvious

        """

        f = open(file_name, 'rb').read()
        f = zlib.compress(f)

        # Strip out any hard coded path in the file name
        file_name = os.path.split(file_name)[1]
        self.data.append((file_name, f, 'file'))

    def _reduce(self):
        """
        Returns a _reduced form of Job.jdict to be sent over the network.

        """

        # dump and compress the data of the job
        jdict = copy.deepcopy(self.__dict__)
        jdict['data'] = cPickle.dumps(self.data, 2)
        jdict['data'] = zlib.compress(jdict['data'])
        # We do not compress jdict['result'] since it's already compressed
        jdict['result'] = self.result

        # Remove attributes that sqlalchemy put there for us
        for key in jdict.keys():
            if key.startswith('_'):
                del jdict[key]

        return jdict

def expand_job(jdict):
    """
    This method recreates a Job object given a jdict.

    :type jdict: dictionary
    :param jdict: the job dictionary

    """

    if jdict is None:
        return None

    job_id = jdict['job_id']
    job = Job(job_id=job_id)

    # decompress and load data
    try:
        jdict['data'] = zlib.decompress(jdict['data'])
        jdict['data'] = cPickle.loads(jdict['data'])
    except (KeyError, TypeError):
        jdict['data'] = None

    try:
        jdict['result'] = zlib.decompress(jdict['result'])
        jdict['result'] = cPickle.loads(jdict['result'])
    except (KeyError, TypeError):
        jdict['result'] = None

    for k, v in jdict.iteritems():
        setattr(job, k, v)

    return job