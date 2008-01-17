"""nodoctest"""

"""
Enumeration of Totally Real Fields: DSage components

AUTHORS:
    -- John Voight (2007-11-04):
        * Changes made after DSage bug fixes.
    -- John Voight (2007-10-27):
        * Initial version.
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein and John Voight
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
#*****************************************************************************

import time
from sage.rings.number_field.totallyreal_data import tr_data
from sage.rings.number_field.totallyreal import timestr, weed_fields
from sage.libs.pari.gen import pari

class totallyreal_dsage:
    counts = [0,0,0,0]
    S = []
    cputime = 0
    walltime = 0

    def __init__(self, D=None, timeout = 30*60):
        r"""
        Initializes the totallyreal DSage instance.

        INPUT:
        D -- (default: None) a DSage instance
        timeout -- (default: 30*60) the default timeout

        OUTPUT:
        the totallyreal DSage instance

        EXAMPLES:
        First, in the server window:

        sage: dsage.setup_all()
        [...]
        sage: dsage.server()
        [...]

        Now, setup workers.  On each machine that wants to be a worker:

        sage: dsage.setup_worker()
        [...]
        sage: dsage.worker()
        [...]

        Finally, at a client:

        sage: Dtr = totallyreal_dsage()
        """

        if D == None:
            self.D = DSage()
        else:
            self.D = D
        self.timeout = timeout

    def enumerate(self, n, B, k, load=True):
        """
        Runs the enumeration algorithm in enumerate_totallyreal_fields in
        a distributed computing environment.
        Once coefficients the coefficients a[k+1...n] are obtained,
        the rest are performed by workers.

        INPUT:
        n -- integer, the degree
        B -- integer, the discriminant bound
        k -- integer, the sequence length

        OUTPUT:
        the updated totallyreal dsage instance with jobs running.

        EXAMPLES:
        We enumerate totally real fields of degree 6 and root discriminant <= 15.
        The workers are given coefficient ranges up to k = 3, i.e.
        [?,?,?,?,!,!,1].

        sage: Dtr.enumerate(6, 15^6, 4)
        [...]
        sage: Dtr.compile_fields()

        [[300125, x^6 - x^5 - 7*x^4 + 2*x^3 + 7*x^2 - 2*x - 1],
         [371293, x^6 - x^5 - 5*x^4 + 4*x^3 + 6*x^2 - 3*x - 1],
         [434581, x^6 - 2*x^5 - 4*x^4 + 5*x^3 + 4*x^2 - 2*x - 1],
         [453789, x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1],
         [485125, x^6 - 2*x^5 - 4*x^4 + 8*x^3 + 2*x^2 - 5*x + 1],
         [592661, x^6 - x^5 - 5*x^4 + 4*x^3 + 5*x^2 - 2*x - 1],
         [703493, x^6 - 2*x^5 - 5*x^4 + 11*x^3 + 2*x^2 - 9*x + 1],
         [722000, x^6 - x^5 - 6*x^4 + 7*x^3 + 4*x^2 - 5*x + 1]]

        """

        # We first precompute the a[k+1..n] possibilities.
        # Initialize.
        self.n = n
        self.B = B
        self.k = k
        T = tr_data(n,B,[])
        f_out = [0]*n + [1]

        # Increment, and halt at level k.
        # Feed tasks to the workers, one for each a.
        self.prestr = 'DSAGE_RESULT=enumerate_totallyreal_fields(' + str(n) + ',' + str(B) + ','
        self.jobs = []
        # Small tweak.
        if k <> n-1:
            T.incr(f_out, haltk=k)

        while f_out[n] <> 0:
            a = f_out[k:n+1]
            # Dsage handles it all here.
            if load:
                job = self.D.eval(self.prestr + str(a) + ',return_seqs=True)', timeout = self.timeout)
                self.jobs.append([a,job])
            else:
                self.jobs.append([a,''])
            T.incr(f_out, haltk=k)

    def init_jobs(self, n, B, A, split=False):
        self.n = n
        self.B = B
        f_out = [0]*n + [1]
        self.jobs = []
        self.prestr = 'DSAGE_RESULT=enumerate_totallyreal_fields(' + str(n) + ',' + str(B) + ','

        if split:
            for a0 in A:
                T = tr_data(n, B, a0)
                f_out = [0]*n + [1]
                k_new = n-len(a0)
                T.incr(f_out, 0, haltk=k_new)
                while f_out[self.n] <> 0:
                    a = f_out[k_new:self.n+1]
                    # Dsage handles it all here.
                    job = self.D.eval(self.prestr + str(a) + ',return_seqs=True)', timeout = self.timeout)
                    self.jobs.append([a,job])
                    T.incr(f_out, haltk=k_new)
            self.k = k_new
        else:
            for a in A:
                job = self.D.eval(self.prestr + str(a) + ',return_seqs=True)', timeout = self.timeout)
                self.jobs.append([a,job])

    def recover(self):
        jobs = self.D.get_my_jobs(True)
        find_str = '(' + str(self.n) + ',' + str(self.B)
        jobs_v = [j[0] for j in self.jobs]
        for i in range(len(jobs)):
            jc = jobs[i].code
            try:
                jc.index(find_str)
                v = eval(jc[jc.find('['):jc.find(']')+1])
                self.jobs[jobs_v.index(v)][1] = jobs[i]
            except ValueError:
                v = []

    def load_jobs(self, A, split=False):
        self.init_jobs(self.n, self.B, A, split)

    def restart_job(self, i, force=False):
        if force or self.jobs[i][1].status == 'completed':
            job = self.D.eval(self.prestr + str(self.jobs[i][0]) + ',return_seqs=True)', timeout = self.timeout)
            self.jobs[i][1] = job

    def restart_all_jobs(self, force=False):
        for i in range(len(self.jobs)):
            self.restart_job(i, force)

    def split_job(self, i, force=False):
        if force or self.jobs[i][1].output.find('computation timed out') <> -1:
            job = self.jobs.pop(i)
            T = tr_data(self.n, self.B, job[0])
            f_out = [0]*self.n + [1]
            k_new = self.n-len(job[0])
            T.incr(f_out, 0, haltk=k_new)
            while f_out[self.n] <> 0:
                a = f_out[k_new:self.n+1]
                # Dsage handles it all here.
                job = self.D.eval(self.prestr + str(a) + ',return_seqs=True)', timeout = self.timeout)
                self.jobs.append([a,job])
                T.incr(f_out, haltk=k_new)

    def split_all_jobs(self, force=False):
        if force:
            for i in range(len(self.jobs)):
                self.split_job(0,force)
        else:
            i = 0
            while i < len(self.jobs):
                if self.jobs[i][1].output.find('computation timed out') <> -1:
                    self.split_job(i,True)
                else:
                    i += 1

    def status(self):
        r"""
        Prints the current status of all jobs.
        """

        for i in range(len(self.jobs)):
            if type(self.jobs[i][1]) <> str:
                print self.jobs[i][0], self.jobs[i][1].status

    def compile_fields(self, write_result=True):
        r"""
        Pop any finished jobs and compiles into the master list.
        If save_result == True, then save the output of the result.
        """

        filename_start = str(self.n) + '/' + str(self.n) + '-' + str(self.B) + '-'

        i = 0
        # For each completed job, add the fields to the list.
        while i < len(self.jobs):
            if type(self.jobs[i][1]) <> str:
                # self.jobs[i][1].get_job()
                if self.jobs[i][1].status == 'completed' and not self.jobs[i][1].result in ['None', '']:
                    print "Compiling job", i, "with", self.jobs[i][0]
                    job = self.jobs[i]

                    # Add the timings.
                    self.cputime += job[1].cpu_time
                    self.walltime += job[1].wall_time

                    if write_result:
                        fsock = open(filename_start + str(job[0]).replace(' ','') + '.out', 'w')
                        fsock.write("Cpu time = " + timestr(job[1].cpu_time) + "\n")
                        fsock.write("Wall time = " + timestr(job[1].wall_time) + "\n")

                    # Output from dsage comes as a string, so convert.
                    S = job[1].result
                    # Add the counts of polynomials checked.
                    for j in range(4):
                        self.counts[j] += S[0][j]
                    if write_result:
                        fsock.write("Counts: " + str(S[0]) + "\n")
                        fsock.write("Total number of fields: " + str(len(S[1])) + "\n\n")

                    # Convert the sequences to pari objects, and compile.
                    S = S[1]
                    for j in range(len(S)):
                        S[j][1] = pari(str(S[j][1])).Polrev()
                    self.S += S
                    if write_result:
                        for j in range(len(S)):
                            fsock.write(str(S[j]) + "\n")

                    fsock.close()
                    self.jobs.pop(i)

                # Otherwise, continue.
                # Note we do not add 1 to i if we just popped jobs!
                else:
                    i += 1
            else:
                i += 1

        # Now check for copies of fields and isomorphic fields.
        self.S.sort()
        i = 0
        while i < len(self.S)-1:
            if self.S[i] == self.S[i+1]:
                self.S.pop(i)
            else:
                i += 1
        weed_fields(self.S)

        if write_result:
            fsock = open(filename_start + 'all.out', 'w')
            fsock.write("Cpu time: " + timestr(self.cputime) + "\n")
            fsock.write("Wall time: " + timestr(self.walltime) + "\n")
            fsock.write("Counts: " + str(self.counts) + "\n")
            fsock.write("Total number of fields: " + str(len(self.S)) + "\n\n")

            if len(self.jobs) > 0:
                fsock.write("Incomplete jobs:\n")
                for i in range(len(self.jobs)):
                    fsock.write(str(self.jobs[i][0]) + "\n")
                fsock.write("\n")

            for j in range(len(self.S)):
                fsock.write(str([self.S[j][0],self.S[j][1].reverse().Vec()]) + "\n")
            fsock.close()

            fsock = open(filename_start + 'reload.dat', 'w')
            fsock.write(str([[self.jobs[i][0] for i in range(len(self.jobs))],
                             [[s[0], s[1].reverse().Vec()] for s in self.S],
                             self.counts, self.cputime, self.walltime]))
            fsock.close()
            print self.num_jobs(), "remaining..."

    def reload(self, n, B, load=True):
        self.n = n
        self.B = B
        self.prestr = 'DSAGE_RESULT=enumerate_totallyreal_fields(' + str(n) + ',' + str(B) + ','
        filename_start = str(n) + '/' + str(n) + '-' + str(B) + '-'

        fsock = open(filename_start + 'reload.dat', 'r')
        data = eval(fsock.read())
        A = data[0]
        self.S = data[1]
        for i in range(len(self.S)):
            self.S[i][1] = pari(self.S[i][1]).Polrev()
        self.counts = data[2]
        self.cputime = data[3]
        self.walltime = data[4]
        if load:
            self.load_jobs(A)
        else:
            self.jobs = [[a,''] for a in A]

    def wait_split_save(self):
        while True:
            time.sleep(self.timeout)
            self.compile_fields()
            self.split_all_jobs()

    def num_jobs(self):
        return len(self.jobs)
