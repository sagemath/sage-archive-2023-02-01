r"""
Interface to SAGE

This is an expect interface to \emph{another} copy of the \sage
interpreter.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import cPickle, os, time

from expect import Expect, ExpectElement, FunctionElement, tmp
import sage.misc.preparser

from sage.structure.sage_object import dumps, loads, load


class Sage(Expect):
    r"""
    Expect interface to the \sage interpreter itself.

    INPUT:
        server -- (optional); if specified runs SAGE on a remote machine with
                  address.  You must have ssh keys setup so you can login to
                  the remote machine by typing "ssh remote_machine" and no password, e.g.:
                     cd; ssh-keygen -t rsa; scp .ssh/id_rsa.pub remote_machine:.ssh/authorized_keys2

                  The version of SAGE should be the same as on the
                  local machine, since pickling is used to move data
                  between the two SAGE process.

    EXAMPLES:
    We create an interface to a copy of \sage.  This copy of \sage runs
    as an external process with its own memory space, etc.

        sage: s = Sage()

    Create the element 2 in our new copy of \sage, and cubeit.
        sage: a = s(2)
        sage: a^3
        8

    Create a vector space of dimension $4$, and compute its generators:
        sage: V = s('QQ^4')
        sage: V.gens()
        ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))

    Note that V is a not a vector space, it's a wrapper around an object
    (which happens to be a vector space), in another running instance
    of \sage.
        sage: type(V)
        <class 'sage.interfaces.sage0.SageElement'>
        sage: V.parent()
        Sage
        sage: g = V.0;  g
        (1, 0, 0, 0)
        sage: g.parent()
        Sage

    We can still get the actual parent by using the name attribute of g,
    which is the variable name of the object in the child process.
        sage: s('%s.parent()'%g.name())
        Vector space of dimension 4 over Rational Field

    Note that the memory space is completely different.
        sage: x = 10
        sage: s('x = 5')
        5
        sage: x
        10
        sage: s('x')
        5

    We can have the child interpreter itself make another child
    \sage process, so now three copies of \sage are running:
        sage: s3 = s('Sage()')
        sage: a = s3(10)
        sage: a
        10

    This $a=10$ is in a subprocess of a subprocesses of your original \sage.

        sage: _ = s.eval('%s.eval("x=8")'%s3.name())
        sage: s3('"x"')
        8
        sage: s('x')
        5
        sage: x
        10

    The double quotes are needed because the call to s3 first evaluates
    its arguments using the s interpeter, so the call to s3 is passed
    \code{s('"x"')}, which is the string \code{"x"} in the s interpreter.
    """
    def __init__(self, logfile   = None,
                       preparse  = True,
                       python    = False,
                       init_code = None,
                       server    = None,
                       **kwds):
        if python:
            if server:
                # It is important that the cleaner stay running no matter what.
                command = "nohup sage -cleaner & sage -python -u"
            else:
                command = "sage -python -u"
            prompt = ">>>"
            if init_code is None:
                init_code = ['from sage.all import *', 'import cPickle']
        else:
            command = "sage"
            prompt = "sage: "
            if init_code is None:
                init_code = ['import cPickle']

        Expect.__init__(self,
                        name = 'sage',
                        prompt = prompt,
                        command = command,
                        restart_on_ctrlc = False,
                        logfile = logfile,
                        init_code = init_code,
                        server = server,
                        **kwds
                        )
        self._preparse = preparse
        self._is_local = (server is None)

    def is_local(self):
        return self._is_local

    def trait_names(self):
        return eval(self.eval('globals().keys()'))

    def quit(self, verbose=False):
        import signal
        if not self._expect is None:
            pid = self._expect.pid
            if verbose:
                print "Exiting spawned %s process (pid=%s)."%(self, pid)
            try:
                for i in range(10):   # multiple times, since clears out junk injected with ._get, etc.
                    self._expect.sendline(chr(3))  # send ctrl-c
                    self._expect.sendline('quit_sage(verbose=%s)'%verbose)
                    self._so_far(wait=0.2)

                os.killpg(pid, 9)
                os.kill(pid, 9)

            except (RuntimeError, OSError), msg:
                pass

            try:
                os.killpg(pid, 9)
                os.kill(pid, 9)
            except OSError:
                pass

            try:
                self._expect.close(signal.SIGQUIT)
            except Exception:
                pass
            self._expect = None

            F = '%s/tmp/%s'%(os.environ['SAGE_ROOT'], pid)
            if os.path.exists(F):
                O = open(F).read()

    def _remote_tmpfile(self):
        try:
            return self.__remote_tmpfile
        except AttributeError:
            self.__remote_tmpfile = eval(self.eval('import sage.interfaces.expect as e; e.tmp'))
            return self.__remote_tmpfile

    def _send_tmpfile_to_server(self):
        cmd = 'scp "%s" %s:"%s" 1>&2 2>/dev/null'%(tmp, self._server, self._remote_tmpfile())
        #print cmd
        os.system(cmd)

    def _get_object_from_server_tmpfile(self):
        cmd = 'scp %s:"%s_get.sobj" "%s_get.sobj" 1>&2 2>/dev/null'%( self._server, self._remote_tmpfile(), tmp)
        #print cmd
        os.system(cmd)
        return load(tmp + "_get")

    def __call__(self, x):
        if isinstance(x, SageElement) and x.parent() is self:
            return x
        if isinstance(x, str):
            return SageElement(self, x)

        if self.is_local():
            open(tmp,'w').write(cPickle.dumps(x,2))
            return SageElement(self, 'cPickle.load(open("%s"))'%tmp)
        else:
            open(tmp,'w').write(dumps(x))   # my dumps is compressed by default
            self._send_tmpfile_to_server()
            return SageElement(self, 'loads(open("%s").read())'%self._remote_tmpfile())

    def __reduce__(self):
        return reduce_load_Sage, tuple([])

    def _quit_string(self):
        return 'from sage.misc.misc import delete_tmpfiles; delete_tmpfiles()'

    def preparse(self, x):
        return sage.misc.preparser.preparse(x)

    def eval(self, line, strip=True):
        """
        Send the code x to a second instance of the \sage interpreter and
        return the output as a string.

        This allows you to run two completely independent copies of \sage
        at the same time in a unified way.

        INPUT:
            line -- input line of code
            strip -- ignored
        """
        if self._preparse:
            line = self.preparse(line)
        return Expect.eval(self, line)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s'%(var,value)
        out = self.eval(cmd)
        if 'Traceback' in out:
            raise TypeError, "Error executing code in SAGE\nCODE:\n\t%s\nSAGE ERROR:\n\t%s"%(cmd, out)

    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval('print %s'%var)

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
    #    self.eval('del %s'%var)

    def _contains(self, v1, v2):
        return self.eval('%s in %s'%(v1,v2))

    def _is_true_string(self, t):
        return t == "True"

    def console(self):
        sage0_console()

    def version(self):
        return eval(sage0_version())

    def _object_class(self):
        return SageElement

    def new(self, x):
        return SageElement(self, x)


class SageElement(ExpectElement):
    def __getattr__(self, attrname):
        self._check_valid()
        return SageFunction(self, attrname)

    def _sage_(self):
        P = self.parent()
        if P.is_remote():
            P.eval('save(%s, "%s_get")'%(self.name(), P._remote_tmpfile()))
            return P._get_object_from_server_tmpfile()
        else:
            P.eval('dumps(%s, "%s")'%(self.name(), tmp))
            return loads(open(tmp).read())

class SageFunction(FunctionElement):
    def __call__(self, *args):
        P = self._obj.parent()
        a = [P(x) for x in args]
        b = ','.join([x.name() for x in a])
        z = SageElement(P, '%s.%s(%s)'%(self._obj._name, self._name, b))
        return z

    def __repr__(self):
        return str(self.eval('%s.%s'%(self._obj._name, self._name)))



sage0 = Sage()

def reduce_load_Sage():
    return sage0

def reduce_load_element(s):
    return sage0('loads(%s)'%s)


import os
def sage0_console():
    os.system('sage')

def sage0_version():
    return sage0('version()')

#def irun(filename):
#    """
#    Run the script in filename step-by-step displaying each input line
#
#    This does not work right with for loops, which span multiple lines.
#    """
#    print 'Interactive runing "%s"'%filename
#    for L in open(filename).xreadlines():
#        raw_input("sage: "+L[:-1])
#        print sage0(L)[:-1]

