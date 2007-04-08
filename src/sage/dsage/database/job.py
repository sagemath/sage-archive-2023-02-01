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

from persistent import Persistent

class Job(Persistent):
    r"""
    Defines a Job that gets distributed to clients.

    """

    def __init__(self, id_=None, name=None, code=None, parent=None,
                 username=None, type_='sage'):
        r"""
        Creates a new job.

        Parameters:
        name -- job name
        id -- job id (must be unique)
        code -- job code
        parent -- sets if the job is a parent job
        username -- author of the job (not neccesarily unique)
        type -- the type of the job (file, string, generator)
                defaults to string

        """

        self.jdict = {}

        # Job keywords
        self.jdict['job_id'] = id_
        self.jdict['name'] = name
        self.jdict['code'] = code
        self.jdict['username'] = username
        self.jdict['data'] = []
        self.jdict['output'] = ''
        self.jdict['worker_info'] = None
        # Valid status keywords:
        # new, completed, incomplete, processing
        self.jdict['status'] = 'new'
        self.jdict['creation_time'] = datetime.datetime.now()
        self.jdict['update_time'] = None
        self.jdict['finish_time'] = None
        self.jdict['killed'] = False
        self.jdict['type'] = type_
        self.jdict['result'] = None # result should be a pickled object
        self.jdict['failures'] = 0
        self.jdict['verifiable'] = False # is this job easily verified?

        # These might become deprecated
        self.jdict['parent'] = parent
        self.jdict['children'] = []

    def __str__(self):
        return str(self.jdict)

    def __setattr__(self, name, value):
        if name == 'jdict':
            if not self.__dict__.has_key('jdict'):
                self.__dict__[name] = value
            else:
                raise ValueError, 'Do not reassign Job.jdict.'
        else:
            Persistent.__setattr__(self, name, value)

    def num_of_children(self):
        return len(self.jdict['children'])

    def get_name(self):
        return self.jdict['name']
    def set_name(self, value):
        if not isinstance(value, str):
            raise TypeError
        self.jdict['name'] = value
    name = property(fget=get_name, fset=set_name, fdel=None, doc='Job name')

    def get_code(self):
        return self.jdict['code']
    def set_code(self, value):
        if not isinstance(value, str):
            raise TypeError
        self.jdict['code'] = value
    code = property(fget=get_code, fset=set_code, fdel=None,
                    doc='Job code')

    def get_id(self):
        return self.jdict['job_id']
    def set_id(self, value):
        if not isinstance(value, str):
            raise TypeError
        self.jdict['job_id'] = value
    job_id = property(fget=get_id, fset=set_id, fdel=None, doc='Job ID')

    def get_status(self):
        return self.jdict['status']
    def set_status(self, value):
        # statuses = ['new', 'completed', 'incomplete', 'processing']
        # if not value in statuses:
        #    raise TypeError
        #if value == 'completed':
        #    self.finish_time = datetime.datetime.now()
        self.jdict['status'] = value
    status = property(fget=get_status, fset=set_status, fdel=None,
                      doc='Job status')

    def get_result(self):
        # loads the result
        try:
            result = cPickle.loads(self.jdict['result'])
        except Exception, msg1:
            try:
                result = cPickle.loads(zlib.decompress(self.jdict['result']))
            except Exception, msg2:
                try:
                    result = cPickle.loads(bz2.decompress(self.jdict['result']))
                except:
                    result = self.jdict['result']
        return result
    def set_result(self, value):
        self.jdict['result'] = value
    result = property(fget=get_result, fset=set_result, fdel=None,
                      doc='Job result')

    def get_data(self):
       return self.jdict['data']
    def set_data(self, value):
       self.jdict['data'] = value
    data = property(fget=get_data, fset=set_data, fdel=None,
                  doc='Job data')

    def get_output(self):
        return self.jdict['output']
    def set_output(self, value):
        self.jdict['output'] = value
    output = property(fget=get_output, fset=set_output, fdel=None,
                      doc='Job output')

    def get_username(self):
        return self.jdict['username']
    def set_username(self, value):
        self.jdict['username'] = value
    username = property(fget=get_username, fset=set_username, fdel=None,
                      doc='Job author')

    def get_finish_time(self):
        return self.jdict['finish_time']
    def set_finish_time(self, value):
        if not isinstance(value, datetime.datetime):
            raise TypeError
        self.jdict['finish_time'] = value
    finish_time = property(fget=get_finish_time, fset=set_finish_time,
                           fdel=None, doc='Job finish time')

    def get_update_time(self):
        return self.jdict['update_time']
    def set_update_time(self, value):
        if not isinstance(value, datetime.datetime):
            raise TypeError
        self.jdict['update_time'] = value
    update_time = property(fget=get_update_time, fset=set_update_time,
                            fdel=None, doc='Job updated time')

    def get_creation_time(self):
        return self.jdict['creation_time']

    def set_creation_time(self, value):
        if not isinstance(value, datetime.datetime):
            raise TypeError
        self.jdict['creation_time'] = value
    creation_time = property(fget=get_creation_time, fset=set_creation_time,
                             fdel=None, doc='Job creation time')

    def get_type(self):
        return self.jdict['type']
    def set_type(self, value):
        types = ['sage', 'spyx', 'py']
        if value not in types:
            raise TypeError

        self.jdict['type'] = value
    type = property(fget=get_type, fset=set_type,
                    fdel=None, doc='Job type')

    def get_failures(self):
        return self.jdict['failures']
    def set_failures(self, value):
        if not isinstance(value, int):
            raise TypeError
        self.jdict['failures'] = value
    failures = property(fget=get_failures, fset=set_failures,
                        fdel=None, doc='Number of failures')

    def get_killed(self):
        return self.jdict['killed']
    def set_killed(self, value):
        if not isinstance(value, bool):
            raise TypeError
        self.jdict['killed'] = value
    killed = property(fget=get_killed, fset=set_killed,
                      fdel=None, doc='Job killed status')

    def get_worker_info(self):
        return self.jdict['worker_info']
    def set_worker_info(self, value):
        self.jdict['worker_info'] = value
    worker_info = property(fget=get_worker_info, fset=set_worker_info,
                           fdel=None, doc='Worker info')

    def get_verifiable(self):
        return self.jdict['verifiable']
    def set_verifiable(self, value):
        if not isinstance(value, bool):
            raise TypeError
        self.jdict['verifiable'] = value

    def attach(self, var, obj, file_name=None):
        r"""
        Attaches an object to a job.

        Parameters:
        var -- the variable name you'd like the worker to use
        obj -- the object you want to attach
        filename -- optional, if your object is a saved sobj

        """

        if file_name:
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
                return
        self.jdict['data'].append((var, s, 'object'))

    def attach_file(self, file_name):
        r"""
        Attach a file to a job.

        Parameters:
        file_name -- obvious

        """

        f = open(file_name, 'rb').read()
        f = zlib.compress(f)

        # Strip out any hard coded path in the file name
        file_name = os.path.split(file_name)[1]
        self.jdict['data'].append((file_name, f, 'file'))

    def pickle(self):
        r"""
        Returns a pickled representation of self.

        """

        s = cPickle.dumps(self, 2)
        s = zlib.compress(s)

        return s

    def unpickle(self, pickled_job):
        r"""
        Returns the unpickled version of myself.

        """

        return cPickle.loads(zlib.decompress(pickled_job))

    def reduce(self):
        r"""
        Returns a reduced form of Job.jdict to be sent over the network.

        """

        # TODO: Figure out what attributes are safe to delete

        # dump and compress the data of the job
        jdict = copy.deepcopy(self.jdict)
        jdict['data'] = cPickle.dumps(self.jdict['data'], 2)
        jdict['data'] = zlib.compress(jdict['data'])

        return jdict

def expand_job(jdict):
    r"""
    This method recreates a Job object given a jdict.

    """

    if jdict is None:
        return None

    job = Job()

    # decompress and load data
    try:
        jdict['data'] = zlib.decompress(jdict['data'])
        jdict['data'] = cPickle.loads(jdict['data'])
    except:
        jdict['data'] = None
    job.jdict.update(jdict)

    return job