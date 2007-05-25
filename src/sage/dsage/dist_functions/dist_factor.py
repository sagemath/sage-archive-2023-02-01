from sage.dsage.database.job import Job
from sage.dsage.dist_functions.dist_function import DistributedFunction
from sage.dsage.interface.dsage_interface import JobWrapper

from sage.all import *

class DistributedFactor(DistributedFunction):
    """
    DistributedFactor uses ECM and QSIEVE to find factors of numbers.

       DistributedFactor will first perform trial division on the number and
       then use ECM.

       AUTHORS:
           Robert Bradshaw
           Yi Qiang
    """

    def __init__(self, DSage, n, concurrent=10, verbosity=0,
                 trial_division_limit=1000000, name='DistributedFactor',
                 use_qsieve=False):
        """
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
        self.n = n
        self.prime_factors = []
        self.composite_factors = []
        self.cur_B1 = 2000
        self.curve_count = 50
        self.concurrent = concurrent
        self.verbosity = verbosity
        self.name = name
        self.use_qsieve = use_qsieve
        # Trial division first to peel off some factors
        for d in prime_range(2, trial_division_limit):
            while d.divides(n):
                self.prime_factors.append(d)
                n = n // d
        if n == 1:
            self.done = True
        elif is_prime(n): # The last value might be prime
            self.done = True
            self.prime_factors.append(n)
        else:
            self.composite_factors.append(n)
            if self.use_qsieve:
                self.outstanding_jobs = [self.qsieve_job()]
            for i in range(concurrent-1):
                self.outstanding_jobs.append(self.ecm_job())

    def next_job(self):
        return self.ecm_job()

    def qsieve_job(self):
        n = max(self.composite_factors)
        if n < 10**40:
            return self.ecm_job()
        job = Job(code="""
n = %s
if is_prime(n):
    result = [[n], [True], {}, 'primality']
else:
    v, t = qsieve(n)
    DSAGE_RESULT = [v, [0]*len(v), {}, 'qsieve']
""" % n, name='qsieve')
        job.n = int(n) # otherwise cPickle will crash
        job.algorithm = 'qsieve'
        job.verifiable = True
        job.type = 'qsieve'
        job.timeout = 60*60*24

        return job

    def ecm_job(self):
        try:
            self.i += 1
        except AttributeError:
            self.i = 0
        n = self.composite_factors[self.i % len(self.composite_factors)]
        rate_multiplier = float(self.concurrent / len(self.composite_factors))
        job = Job(code="""
n = %s
if is_prime(n):
    DSAGE_RESULT = [[n], [True], {}, 'primality']
else:
    e = ECM()
    DSAGE_RESULT = [e.find_factor(n, B1=%s, c=%s, I=%s), e.primality, e.last_params, 'ecm']

""" % (n, self.cur_B1, self.curve_count, rate_multiplier), name='ecm' )
        job.n = int(n)
        job.algorithm = 'ecm'
        job.verifiable = True
        job.type = 'ecm'
        job.timeout = 60*60*24
        return job

    def process_result(self, job):
        """
        For each factor m of n found by the worker, record them by
            1) Dividing each element x of composite_factors by gcd(x,m)
            2) Storing (non-trivial) gcd(x,m) to the composite factor list.
            3) Storing m in prime_factors or composite_factors according to
               its classification given by the worker.
        If the factorization is not yet complete, spawn another job.

        """

        if prod(self.prime_factors) == self.n:
            self.done = True
            return

        for factor in self.composite_factors:
            if is_prime(factor):
                self.prime_factors.append(factor)
                self.composite_factors.remove(factor)

        if self.verbosity > 2:
            print "process_result(): ", job, job.output, job.result
        result = job.result
        if self.verbosity > 1:
            print "factors:", self.prime_factors, self.composite_factors

        # If result is unexpected...
        try:
            factors, primality, params, algorithm = result
        except Exception, msg:
            print 'Error in processing result.'
            print result
            print msg
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
            result = [(result[0][i], result[1][i]) for i in range(len(result[0]))]
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
            self.qsieve_count = 0
            to_be_removed_jobs = []
            for wrapped_job in self.waiting_jobs:
                if wrapped_job.algorithm == 'qsieve':
                    if ZZ(wrapped_job.n) not in self.composite_factors:
                        if self.verbosity > 2:
                            print "killing qsieve(%s)" % wrapped_job.n
                        if isinstance(wrapped_job, JobWrapper):
                            wrapped_job.kill()
                        else:
                            wrapped_job.async_kill()
                        to_be_removed_jobs.append(wrapped_job)
                    else:
                        self.qsieve_count += 1
            for job in to_be_removed_jobs:
                self.waiting_jobs.remove(job)
            if self.use_qsieve:
                if self.qsieve_count == 0:
                    self.submit_job(self.qsieve_job(), self.name, async=True)
            self.submit_job(self.ecm_job(), self.name, async=True)

        self.prime_factors.sort()
        self.composite_factors.sort()