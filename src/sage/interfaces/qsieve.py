"""
Interface to Bill Hart's Quadratic Sieve
"""

import os

import sage.rings.integer

def qsieve(n, block=True, time=False, verbose=False):
    """
    Run Hart's quadratic sieve and return the distinct proper factors
    of the integer n that it finds.

    INPUT:
        n -- an integer with at least 40 digits
        block -- (default: True) if True, you must wait until the
            sieve computation is complete before using SAGE further.
            If False, SAGE will run while the sieve computation
            runs in parallel.
        time -- (default: False) if True, time the command using
            the UNIX "time" command (which you might have to install).
        verbose -- (default: False) if True, print out verbose
            logging information about what happened during
            the Sieve run (for non-blocking Sieve, verbose information
            is always available via the log() method.)

    OUTPUT:
        list -- a list of the distinct proper factors of n found
        str -- the time in cpu seconds that the computation took, as given
               by the command line time command.  (If time is False,
               this is always an empty string.)

    EXAMPLES:
        sage: k = 19; n = next_prime(10^k)*next_prime(10^(k+1))
        sage: factor(n)  # (currently) uses PARI
        10000000000000000051 * 100000000000000000039
        sage: v, t = qsieve(n, time=True)   # uses the sieve
        sage: v
        [10000000000000000051, 100000000000000000039]
        sage: t   # random output
    """
    Z = sage.rings.integer.Integer
    n = Z(n)
    if len(str(n)) < 40:
        raise ValueError, "n must have at least 40 digits"
    if block:
        return qsieve_block(n, time, verbose)
    else:
        return qsieve_nonblock(n, time)

def qsieve_block(n, time, verbose=False):
    """
    Compute the factorization of n using Hart's quadratic Sieve
    blocking until complete.
    """
    if time:
        t = 'time '
    else:
        t = ''
    out = os.popen('echo "%s" | %s QuadraticSieve 2>&1'%(n,t)).read()
    z = data_to_list(out, n, time=time)
    if verbose:
        print z[-1]
    return z[:2]

def data_to_list(out, n, time):
    """
    Convert output of Hart's sieve and n to a list and time.

    INPUT:
        out -- snapshot of text output of Hart's QuadraticSieve program
        n -- the integer being factored

    OUTPUT:
        list -- proper factors found so far
        str -- time information
    """
    i = out.find('FACTORS:')
    if i == -1:
        verbose = ''
    else:
        verbose = out[:i]
        out = out[i+len('FACTORS:')+1:].strip()
    if time:
        w = out.split('\n')
        for i in range(len(w)):
            if 'user' in w[i]:
                break
        if i < len(w):
            t = w[i].strip()
            out = '\n'.join([w[j] for j in range(i)])
        else:
            t = ''
    else:
        t = ''
    Z = sage.rings.integer.Integer
    v = out.split()
    v = list(set([Z(m) for m in v if Z(m) != n]))
    v.sort()
    return v, t, verbose


import pexpect
import monitor
class qsieve_nonblock:
    """
    A non-blocking version of Hart's quadratic sieve.

    The sieve starts running when you create the object, but you can
    still use SAGE in parallel:

        sage: k = 19; n = next_prime(10^k)*next_prime(10^(k+1))
        sage: q = qsieve(n, block=False)
        sage: q     # random output
        Proper factors so far: []
        sage: q     # random output
        ([10000000000000000051, 100000000000000000039], '0.21')
        sage: q.list()
        [10000000000000000051, 100000000000000000039]
        sage: q.time()    # random output
        '0.21'
    """
    def __init__(self, n, time):
        self._n = n
        if time:
            cmd = 'time QuadraticSieve'
        else:
            cmd = 'QuadraticSieve'
        self._p = pexpect.spawn(cmd)
        monitor.monitor(self._p.pid)
        self._p.sendline(str(self._n)+'\n\n\n')
        self._done = False
        self._out = ''
        self._time = ''
        self._do_time = time

    def n(self):
        """
        Return the integer that is being factored.
        """
        return self._n

    def pid(self):
        """
        Return the PIN id of the QuadraticSieve process (actually
        of the time process that spawns the sieve process).
        """
        return self._p.pid

    def done(self):
        """
        Return True if the sieve process has completed.
        """
        return self._done

    def __repr__(self):
        """
        Return a text representation of self.
        """
        if self._done:
            v = data_to_list(self._get(), self._n, self._do_time)
            if self._do_time:
                return str(v[:2])
            else:
                return str(v[0])
        else:
            return 'Proper factors so far: %s'%self.list()

    def cputime(self):
        """
        Return the time in seconds (as a string) that it took to
        factor n, or return '?' if the factorization has not
        completed or the time is unknown.
        """
        if not self._do_time:
            raise ValueError, "you have to start the seive with the option time=True in order to get timing information"
        try:
            return data_to_list(self._get(), self._n, self._do_time)[1]
        except IndexError:
            return '?'
    time = cputime

    def log(self):
        """
        Return all output of running the sieve so far.
        """
        return self._get()

    def __getitem__(self, i):
        """
        Return the i-th factor (in sorted order) found so far.
        """
        return self.list()[i]

    def __len__(self):
        """
        Return the number of factors found so far.  If q is the
        Sieve object, type len(q) to see the number of factors.
        """
        return len(self.list())

    def list(self):
        """
        Return the a list of the factors found so far, as SAGE
        integers.
        """
        try:
            return data_to_list(self._get(), self._n, self._do_time)[0]
        except IndexError:
            return []

    def quit(self):
        """
        Terminate the QuadraticSieve process, in case you want
        to give up on computing this factorization.
        """
        self._p.close()
        self._done = True

    def _get(self, timeout=0.1):
        """
        Used internally to get information about what has been
        computed so far.
        """
        if self._done:
            return self._out
        e = self._p
        try:
            e.expect('xxx', timeout=timeout)
        except pexpect.TIMEOUT, msg:
            pass
        except pexpect.EOF, msg:
            pass
            self._done = True
            self._p.close()
        self._out += e.before
        return self._out



