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

import datetime
import copy
import cPickle
import cStringIO
import zlib
import tarfile

from twisted.internet import task
from dsage.database.job import Job

from sage.all import *

class DistributedFunction(object):
    r"""
    Parent class for all classes that wish to use Distributed SAGE.

    Parameters:
        dsage -- a DSage() object

    """

    def __init__(self, dsage):
        self.DSage = dsage
        self.result = None
        self.outstanding_jobs = [] # unsubmitted jobs
        self.waiting_jobs = [] # unprocessed jobs
        self.processed_jobs = [] # processed jobs
        self.start_time = None
        self.end_time= None
        self.name = None
        self.job_files = []
        self.checker_task = None
        self._done = False

    def __getstate__(self):
        d = copy(self.__dict__)
        d['DSage'] = None
        d['checker_task'] = None
        return d

    def get_done(self):
        return self._done
    def set_done(self, value):
        if value:
            self.end_time = datetime.datetime.now()
        self._done = value
    done = property(fget=get_done, fset=set_done, fdel=None,
                    doc="Set the jobs status")

    def get_time(self):
        if self.end_time == None:
            return "End time not available"
        return str(self.end_time - self.start_time)
    time = property(fget=get_time, fset=None, fdel=None,
                     doc="Time it took for job to complete")

    def save(self, filename=None, compress=True):
        r"""
        Saves your distributed job to a sobj.

        """

        if filename is None:
            filename = str(self.name)
        if not filename.endswith('.sobj'):
            filename += '.sobj'
        s = cPickle.dumps(self, 2)
        if compress:
            s = zlib.compress(s)
        f = open(filename, 'wb')
        f.write(s)
        f.close()
        return filename

    def restore(self, dsage):
        if dsage.remoteobj is None:
            # XXX This is a hack because dsage.remoteobj is not set yet
            from twisted.internet import reactor
            reactor.callLater(1.0, self.restore, dsage)
            return
        self.DSage = dsage
        for job in self.waiting_jobs:
            job.remoteobj = self.DSage.remoteobj
        if not len(self.outstanding_jobs) == 0:
            self.submit_jobs(self.name)
        self.checker_task = task.LoopingCall(self.check_results)
        self.checker_task.start(2.0, now=True)

    def submit_job(self, job, job_name='job'):
        if isinstance(job, Job):
            self.waiting_jobs.append(self.DSage.send_job(job))
        else:
            self.waiting_jobs.append(self.DSage.eval(job, job_name=job_name))

    def submit_jobs(self, job_name='job'):
        for job in self.outstanding_jobs:
           self.submit_job(job, job_name)
        self.outstanding_jobs = []

    def start(self):
        self.start_time = datetime.datetime.now()
        self.submit_jobs(self.name)
        self.checker_task = task.LoopingCall(self.check_results)
        self.checker_task.start(1.0, now=False)

    def check_results(self):
        for wrapped_job in self.waiting_jobs:
            wrapped_job.syncJob()
            # job.job refers to the remove job
            if wrapped_job._job.status == 'completed':
                self.waiting_jobs.remove(wrapped_job)
                self.process_result(wrapped_job._job)
                self.processed_jobs.append(wrapped_job)
        if self.done:
            # kill the jobs in the waiting queue
            for wrapped_job in self.waiting_jobs:
                wrapped_job.kill()
            self.checker_task.stop()
            return
        self.done = len(self.waiting_jobs) == 0

class DistributedTestFunction(DistributedFunction):
    def __init__(self, DSage, n, name='DistributedTestFunction'):
        DistributedFunction.__init__(self, DSage)
        self.n = n
        self.name = name
        self.result = 0
        self.results = []
        self.outstanding_jobs = ["print %s"%i for i in range(1,n+1)]

    def process_result(self, job):
        self.result += int(job.output)


class DistributedFactor(DistributedFunction):
    r"""
    DistributedFactor uses ECM and QSIEVE to find factors of numbers.

       DistributedFactor will first perform trial division on the number and
       then use ECM.

       AUTHORS:
           Robert Bradshaw
           Yi Qiang
    """
    def __init__(self, DSage, n, concurrent=10, verbosity=0,
                 trial_division_limit=10000, name='DistributedFactor'):
        r"""
        Parameters:
            DSage -- an instance of a dsage connection
            n -- the square-free number to be factored
            concurrent -- number of parallel jobs to run
            trial_division_limit -- perform trial division up to this number
                                    before attempting ecm. Defaults to 10000
                                    which finishes quite quickly Set to -1
                                    to skip (if the number is known to
                                    contain no small factors)
            name, verbosity  -- obvious

        """
        DistributedFunction.__init__(self, DSage)
        self.id = "ecm_factor(%s)" % (n)
        self.n = n
        self.prime_factors = []
        self.cur_B1 = 2000
        self.curve_count = 50
        self.concurrent = concurrent
        self.verbosity = verbosity
        self.name = name

        # Trial division first to peel off some factors
        for d in prime_range(2, trial_division_limit):
            while d.divides(n):
                self.prime_factors.append(d)
                n = n // d
        if n == 1:
            self.done = True
        else:
            self.composite_factors = [n]
            self.outstanding_jobs = [self.qsieve_job()]
            for i in range(concurrent-1):
                self.outstanding_jobs.append(self.ecm_job())

    def next_job(self):
        return self.ecm_job()

    def qsieve_job(self):
        n = max(self.composite_factors)
        if n < 10**40:
            return self.ecm_job()
        job = Job(file="""
n = %s
if is_prime(n):
    result = [[n], [True], {}, 'primality']
else:
    v, t = qsieve(n)
    result = [v, [0]*len(v), {}, 'qsieve']
    save(result, 'result')
    DSAGE_RESULT = 'result.sobj'
""" % n, name='qsieve')
        job.n = str(n) # otherwise get some weird twisted class
        job.algorithm = 'qsieve'
        return job

    def ecm_job(self):
        try:
            self.i += 1
        except AttributeError:
            self.i = 0
        n = self.composite_factors[self.i % len(self.composite_factors)]
        rate_multiplier = float(self.concurrent / len(self.composite_factors))
        job = Job(file="""
n = %s
if is_prime(n):
    result = [[n], [True], {}, 'primality']
else:
    e = ECM()
    result = [e.find_factor(n, B1=%s, c=%s, I=%s), e.primality, e.last_params, 'ecm']
save(result, 'result')
DSAGE_RESULT = 'result.sobj'

""" % (n, self.cur_B1, self.curve_count, rate_multiplier), name='ecm' )
        job.n = n
        job.algorithm = 'ecm'
        return job



    def process_result(self, job):
        r"""
        For each factor m of n found by the worker, record them by
            1) Dividing each element x of composite_factors by gcd(x,m)
            2) Storing (non-trivial) gcd(x,m) to the composite factor list.
            3) Storing m in prime_factors or composite_factors according to
               its classification given by the worker.
        If the factorization is not yet complete, spawn another job.

        """
        if prod(self.prime_factors) == self.n:
            print 'Found all prime factors. \r'
            self.done = True
            return
        if self.verbosity > 2:
            print "process_result()", job, job.output
            print job.result
        result = job.result
        if self.verbosity > 1:
            print "factors:", self.prime_factors, self.composite_factors

        # If result is unexpected...
        try:
            factors, primality, params, algorithm = result
        except:
            print 'Error in processing result.'
            return
        try:
            self.cur_B1 = max(self.cur_B1, int(params['B1']))
        except KeyError:
            pass

        found_new_factors = False

        if len(factors) > 1 or primality[0]:
            if self.n % prod(factors) != 0:
                self.prime_factors.append("BAD FACTORS")
            # switch the order of indices
            result = [(result[0][i],
                       result[1][i]) for i in range(len(result[0]))]
            # Only works for square-free numbers
            for r, is_prime_factor in result:
                for p in self.prime_factors:
                    if p.divides(r):
                        r = r // p
                if r == 1:
                    continue
                if r in self.composite_factors:
                    if is_prime_factor:
                        self.composite_factors.remove(r)
                        self.prime_factors.append(r)
                else:
                    for m in self.composite_factors:
                        if is_prime_factor:
                            if r.divides(m):
                                self.composite_factors.remove(m)
                                # we know m != p from above
                                self.composite_factors.append(m//r)
                                self.prime_factors.append(r)
                                break
                        else:
                            g = gcd(r, m)
                            if g > 1:
                                self.composite_factors.remove(m)
                                self.composite_factors.append(g)
                                m = m // g
                                if m != 1:
                                    self.composite_factors.append(m)
                                r = r // g
                                if r == 1:
                                    break
            if self.verbosity > 1:
                print "factors:", self.prime_factors, self.composite_factors

        if len(self.composite_factors) == 0:
            self.result = self.prime_factors
            self.done = True
        else:
            qsieve_count = 0
            for wrapped_job in self.waiting_jobs:
                if wrapped_job._job.algorithm =='qsieve':
                    if ZZ(wrapped_job._job.n) not in self.composite_factors:
                        if self.verbosity > 2:
                            print "killing qsieve(%s)" % wrapped_job._job.n
                        wrapped_job.kill()
                        self.waiting_jobs.remove(wrapped_job)
                    else:
                        qsieve_count += 1

            if qsieve_count == 0:
                self.submit_job(self.qsieve_job(), self.name)
            self.submit_job(self.ecm_job(), self.name)


class DistributedPOVRay(DistributedFunction):
    r"""
    DistributedPOVRay distributes rendering of a .pov file.

    Parameters:
        DSage -- a DSage object
        name -- the name of your job
        files -- a list of files required for your pov job, including
                the pov file itself
        splits -- number of jobs you want to split into
        **kwargs -- parameters you wish to pass to povray
            Note: You must pass in at least a width and height.

    OUTPUT:
        A number of .ppm files (depending on split) and a final .ppm file
        which is the combination of all the rendered parts.

    AUTHOR:
        Yi Qiang

    """

    def __init__(self, DSage, name, files, splits, **kwargs):
        DistributedFunction.__init__(self, DSage)
        self.name = name
        self.files = files
        for f in files:
            if f.endswith('.pov'):
                self.pov_fname = f

        self.splits = splits
        self.width = kwargs['W']
        self.height = kwargs['H']
        self.kwargs = kwargs

        # Figure out how to split up the image into jobs
        self.remainder = self.height % self.splits
        self.step = self.height // self.splits
        self.sr = 1
        self.er = self.sr + self.step
        self.n = 0
        self.n1 = 0
        for i in xrange(splits):
            job = self.next_job()
            if not job == None:
                self.outstanding_jobs.append(job)

    def next_job(self):
        if self.n + 1 == self.splits:
            self.er = self.height
        if self.n + 1 > self.splits:
            return

        job_file = "povray('%s', outfile='%s_%04d.ppm', " % (self.pov_fname,
                                                             self.name,
                                                             self.n)

        for k, v in self.kwargs.iteritems():
            job_file += '%s=%s, ' % (k, v)

        job_file += 'SR=%s, ER=%s, ' % (self.sr, self.er)
        job_file += ')'
        job_file += '\n'
        job_file += "tmp = open('%s_%04d.ppm', 'rb').read()\n" % (self.name,
                                                                  self.n)
        job_file += "save(tmp, '%s_%04d')\n" % (self.name, self.n)
        job_file += "DSAGE_RESULT = '%s_%04d' + '.sobj'\n" % (self.name,
                                                              self.n)
        job_file += '\n'

        self.job_files.append(job_file)

        job = Job(file=job_file, name='%s_%04d.ppm' % (self.name, self.n))
        for file in self.files:
            job.attach_file(file)
        self.n += 1

        self.sr = self.er + 1
        self.er += self.step

        return job

    def process_result(self, job):
        if self.done:
            return
        f = open(job.name, 'wb')
        f.write(job.result)
        f.close()
        self.n1 += 1
        if len(self.waiting_jobs) == 0:
            print "Got all the images, stiching them now.\n"
            cmd = "combineppm %s > %s.ppm" % (self.name, self.name)
            os.system(cmd)
            self.done = True


class DistributedPointCount(DistributedFunction):
      r"""
      DistributedPointCount

      """

      def __init__(self, DSage, a_invariants, primeLimit):
          DistributedFunction.__init__(self, DSage)
          self.a_invariants = a_invariants
          self.primeLimit = primeLimit
          self.points = {}
          for i in primes(primeLimit):
              job = self.next_job(i)
              if not job == None:
                  self.outstanding_jobs.append(job)

      def next_job(self, i):
	  code = """
E = EllipticCurve(GF(%s), a_invariants)
k = len(E.points())
result = [%s, k]
save(result, 'result.sobj')
DSAGE_RESULT = 'result.sobj'
""" % (i, i)
	  job = Job()
	  job.file = code
	  job.attach('a_invariants', self.a_invariants)
	  job.name = 'point counting'
	  return job

      def process_result(self, job):
          self.points[job.result[0]]=job.result[1]
