"""
Common Interface Functionality through Pexpect

See the examples in the other sections for how to use specific
interfaces. The interface classes all derive from the generic
interface that is described in this section.

AUTHORS:

- William Stein (2005): initial version

- William Stein (2006-03-01): got rid of infinite loop on startup if
  client system missing

- Felix Lawrence (2009-08-21): edited ._sage_() to support lists and float exponents in foreign notation.

- Simon King (2010-09-25): Expect._local_tmpfile() depends on
  Expect.pid() and is cached; Expect.quit() clears that cache,
  which is important for forking.

- Jean-Pierre Flori (2010,2011): Split non Pexpect stuff into a parent class.

- Simon King (2010-11-23): Ensure that the interface is started again
  after a crash, when a command is executed in _eval_line. Allow
  synchronisation of the GAP interface.

"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import os
import sys
import weakref
import time
import gc
import operator
import quit
import cleaner
from random import randrange

########################################################
# Important note: We use Pexpect 2.0 *not* Pexpect 2.1.
# For reasons I don't understand, pexpect2.1 is much
# worse than pexpect 2.0 for everything Sage does.
########################################################
import pexpect
from pexpect import ExceptionPexpect
from sage.interfaces.interface import Interface, InterfaceElement, InterfaceFunction, InterfaceFunctionElement, AsciiArtString

from sage.structure.sage_object import SageObject
from sage.structure.parent_base import ParentWithBase
from sage.structure.element import RingElement

import sage.misc.sage_eval
from sage.misc.misc import SAGE_EXTCODE, verbose, SAGE_TMP_INTERFACE, LOCAL_IDENTIFIER
from sage.misc.object_multiplexer import Multiplex

BAD_SESSION = -2

failed_to_start = []

# The subprocess is a shared resource.  In a multi-threaded
# environment, there would have to be a lock to control access to the
# subprocess.  Fortunately, Sage does not use Python threads.
# Unfortunately, the combination of the garbage collector and __del__
# methods gives rise to the same issues.  So, in places where we need
# to do a sequence of operations on the subprocess and make sure
# nothing else intervenes (for example, when we write a command and
# then read back the result) we need to disable the garbage collector.
# See TRAC #955 for a more detailed description of this problem.

# To turn off the garbage collector for a particular region of code,
# do:
#   with gc_disabled():
#       ... your code goes here ...
# The garbage collector will be returned to its original state
# whenever the code exits by any means (falling off the end, executing
# "return", "break", or "continue", raising an exception, ...)

class gc_disabled(object):
    """
    This is a "with" statement context manager. Garbage collection is
    disabled within its scope. Nested usage is properly handled.

    EXAMPLES::

        sage: import gc
        sage: from sage.interfaces.expect import gc_disabled
        sage: gc.isenabled()
        True
        sage: with gc_disabled():
        ...       print gc.isenabled()
        ...       with gc_disabled():
        ...           print gc.isenabled()
        ...       print gc.isenabled()
        False
        False
        False
        sage: gc.isenabled()
        True
    """
    def __enter__(self):
        self._enabled = gc.isenabled()
        gc.disable()

    def __exit__(self, ty, val, tb):
        if self._enabled:
            gc.enable()
        return False

class Expect(Interface):
    """
    Expect interface object.
    """
    def __init__(self, name, prompt, command=None, server=None, server_tmpdir=None,
                 ulimit = None, maxread=100000,
                 script_subdirectory="", restart_on_ctrlc=False,
                 verbose_start=False, init_code=[], max_startup_time=None,
                 logfile = None, eval_using_file_cutoff=0,
                 do_cleaner=True, remote_cleaner=False, path=None,
                 terminal_echo=True):

        Interface.__init__(self, name)
        self.__is_remote = False
        self.__remote_cleaner = remote_cleaner
        if command is None:
            command = name
        if server is not None:
            if ulimit:
                command = 'sage-native-execute ssh -t %s "ulimit %s; %s"'%(server, ulimit, command)
            else:
                command = "sage-native-execute ssh -t %s %s"%(server, command)
            self.__is_remote = True
#            eval_using_file_cutoff = 0  # don't allow this!
            if verbose_start:
                print "Using remote server"
                print command
            self._server = server
            if server_tmpdir is None:
                # TO DO: Why default to /tmp/? Might be better to use the expect process itself to get a tmp folder
                print "No remote temporary directory (option server_tmpdir) specified, using /tmp/ on "+server
                self.__remote_tmpdir = "/tmp/"
            else:
                self.__remote_tmpdir = server_tmpdir
        else:
            self._server = None
        self.__do_cleaner = do_cleaner
        self.__maxread = maxread
        self._eval_using_file_cutoff = eval_using_file_cutoff
        self.__script_subdirectory = script_subdirectory
        self.__command = command
        self._prompt = prompt
        self._restart_on_ctrlc = restart_on_ctrlc
        self.__verbose_start = verbose_start
        if path is not None:
            self.__path = path
        elif script_subdirectory is None:
            self.__path = '.'
        else:
            self.__path = os.path.join(SAGE_EXTCODE,name,self.__script_subdirectory)
        self.__initialized = False
        self.__seq = -1
        self._expect = None
        self._session_number = 0
        self.__init_code = init_code
        self.__max_startup_time = max_startup_time

        #Handle the log file
        if isinstance(logfile, basestring):
            logfile = open(logfile,'w')
        self.__logfile = logfile

        quit.expect_objects.append(weakref.ref(self))
        self._available_vars = []
        self._terminal_echo = terminal_echo

    def _get(self, wait=0.1, alternate_prompt=None):
        if self._expect is None:
            self._start()
        E = self._expect
        wait=float(wait)
        try:
            if alternate_prompt is None:
                E.expect(self._prompt, timeout=wait)
            else:
                E.expect(alternate_prompt, timeout=wait)
        except pexpect.TIMEOUT as msg:
            return False, E.before
        except pexpect.EOF as msg:
            return True, E.before
        except Exception as msg:   # weird major problem!
            return True, E.before
        return True, E.before

    def _send(self, cmd):
        if self._expect is None:
            self._start()
        E = self._expect
        self.__so_far = ''
        E.sendline(cmd)

    def is_running(self):
        """
        Return True if self is currently running.
        """
        if self._expect is None:
            return False
        try:
            os.kill(self._expect.pid,0)
        except OSError:
            # This means the process is not running
            return False
        return True

    def _so_far(self, wait=0.1, alternate_prompt=None):
        """
        Return whether done and output so far and new output since last
        time called.
        """
        done, new = self._get(wait=wait, alternate_prompt=alternate_prompt)
        try:
            if done:
                #if new is not None:
                X = self.__so_far + new
                del self.__so_far
                return True, X, new
            #new = self._expect.before
            try:
                self.__so_far += new
            except (AttributeError, TypeError):
                self.__so_far = new
            return False, self.__so_far, new
        except AttributeError as msg:   # no __so_far
            raise RuntimeError(msg)

    def is_remote(self):
        return self.__is_remote

    def is_local(self):
        return not self.__is_remote

    def user_dir(self):
        return self.__path

    def _change_prompt(self, prompt):
        self._prompt = prompt

#    (pdehaye 20070819: this was used by some interfaces but does not work well remotely)
#    def _temp_file(self, x):
#        T = self.__path + "/tmp/"
#        if not os.path.exists(T):
#            os.makedirs(T)
#        return T + str(x)

    def path(self):
        return self.__path

    def expect(self):
        if self._expect is None:
            self._start()
        return self._expect

    def pid(self):
        """
        Return the PID of the underlying sub-process.

        REMARK:

        If the interface terminates unexpectedly, the original
        PID will still be used. But if it was terminated using
        :meth:`quit`, a new sub-process with a new PID is
        automatically started.

        EXAMPLE::

            sage: pid = gap.pid()
            sage: gap.eval('quit;')
            ''
            sage: pid == gap.pid()
            True
            sage: gap.quit()
            sage: pid == gap.pid()
            False

        """
        if self._expect is None:
            self._start()
        return self._expect.pid

    def _install_hints(self):
        r"""
        Hints for installing needed slave program on your computer.

        There are no hints by default.
        """
        return ''

    def _install_hints_ssh(self):
        r"""
        Hints for installing passwordless authentication on your
        computer...
        """
        # Written by Paul-Olivier Dehaye 2007/08/23
        return """
In order for Sage (on "local") to launch a "slave" process on "remote", the following command needs to work from local's console, without the need to enter any password:

       "ssh -t remote slave",

where "slave" could be "math" (for text-mode Mathematica), "gap", "magma", "sage", "maple", etc.

This thus requires passwordless authentication to be setup, which can be done with commands like these:
        cd; ssh-keygen -t rsa; scp .ssh/id_rsa.pub remote:.ssh/authorized_keys2\n
(WARNING: this would overwrite your current list of authorized keys on "remote")

In many cases, the server that can actually run "slave" is not accessible from the internet directly, but you have to hop through an intermediate trusted server, say "gate".
If that is your case, get help with _install_hints_ssh_through_gate().

"""

    def _install_hints_ssh_through_gate(self):
        r"""
        Hints for installing passwordless authentication through a gate
        """
        # Written by Paul-Olivier Dehaye 2007/08/23
        return """

 We assume you would like to run a "slave" process  on a machine called "remote" from a machine running Sage called "local". We also assume "remote" can only be accessed from "local" by ssh'ing first to "gate" (this is a fairly common setup). Sometimes, "gate" and "remote" have a shared filesystem, and this helps a bit.

  Note: You cannot just create shell scripts on "local" and "gate" that would use two successive SSH connections to "remote" in order to simulate running "slave" locally. This is because Sage will sometimes use files (and scp)  to communicate with "remote", which shell scripts would not take care of.

You need to setup:
 * passwordless authentication to "gate" from "local"
 * add passwordless authentication to "remote" from "local",
   for instance by appending the file local:~/.ssh/id_rsa.pub to remote:~/.ssh/authorized_keys2 and logging in once
      (this is only needed if "remote" and "gate" don\'t share filesystem)
 * add a few lines to your local:~/.ssh/ssh_config. Mine look like

       Host remote_for_sage
            ProxyCommand ssh gate nc -w 1 remote 22

That's it, normally.

The last step tells ssh that whenever an ssh connection is required to
the host "remote_for_sage", it should tunnel it through "gate". Any
attempt to scp-connect to "remote_for_sage" will use ssh and thus
this configuration file, and properly channel those file transfers
through the tunnel.

A good test is to attempt an scp connection from the command-line
of "local" to "remote_for_sage" as if no tunnel through "gate" was
required. No password should be asked for the second time around.

Finally, we created the new name "remote_for_sage" for "remote",
but this name only exists locally. this is to avoid interfering
with any other program that might already ssh to "remote" in
their own way.

If this all works, you can then make calls like:
         math = Mathematica(server="remote_for_sage")

"""


    def _do_cleaner(self):
        try:
            return self.__do_cleaner
        except AttributeError:
            return False

    def _start(self, alt_message=None, block_during_init=True):
        from sage.misc.misc import sage_makedirs
        self.quit()  # in case one is already running
        global failed_to_start

        self._session_number += 1
        current_path = os.path.abspath('.')
        dir = self.__path
        sage_makedirs(dir)
        os.chdir(dir)

        #If the 'SAGE_PEXPECT_LOG' environment variable is set and
        #the current logfile is None, then set the logfile to be one
        #in .sage/pexpect_logs/
        if self.__logfile is None and 'SAGE_PEXPECT_LOG' in os.environ:
            from sage.env import DOT_SAGE
            logs = '%s/pexpect_logs'%DOT_SAGE
            sage_makedirs(logs)

            filename = '%s/%s-%s-%s-%s.log'%(logs, self.name(), os.getpid(), id(self), self._session_number)
            self.__logfile = open(filename, 'w')

        cmd = self.__command

        if self.__verbose_start:
            print cmd
            print "Starting %s"%cmd.split()[0]

        try:
            if self.__remote_cleaner and self._server:
                c = 'sage-native-execute  ssh %s "nohup sage -cleaner"  &'%self._server
                os.system(c)

            # Unset some environment variables for the children to
            # reduce the chances they do something complicated breaking
            # the terminal interface.
            # See Trac #12221 and #13859.
            pexpect_env = dict(os.environ)
            pexpect_del_vars = ['TERM', 'COLUMNS']
            for i in pexpect_del_vars:
                try:
                    del pexpect_env[i]
                except KeyError:
                    pass
            self._expect = pexpect.spawn(cmd, logfile=self.__logfile, env=pexpect_env)
            if self._do_cleaner():
                cleaner.cleaner(self._expect.pid, cmd)

        except (ExceptionPexpect, pexpect.EOF, IndexError):
            self._expect = None
            self._session_number = BAD_SESSION
            failed_to_start.append(self.name())
            raise RuntimeError("Unable to start %s because the command '%s' failed.\n%s"%(
                self.name(), cmd, self._install_hints()))

        os.chdir(current_path)
        self._expect.timeout = self.__max_startup_time
        if not self._terminal_echo:
            self._expect.setecho(0)

        #self._expect.setmaxread(self.__maxread)
        self._expect.maxread = self.__maxread
        self._expect.delaybeforesend = 0
        try:
            self._expect.expect(self._prompt)
        except (pexpect.TIMEOUT, pexpect.EOF) as msg:
            self._expect = None
            self._session_number = BAD_SESSION
            failed_to_start.append(self.name())
            raise RuntimeError("Unable to start %s"%self.name())
        self._expect.timeout = None
        with gc_disabled():
            if block_during_init:
                for X in self.__init_code:
                    self.eval(X)
            else:
                for X in self.__init_code:
                    self._send(X)


    def clear_prompts(self):
        while True:
            try:
                self._expect.expect(self._prompt, timeout=0.1)
            except pexpect.TIMEOUT:
                return

    def __del__(self):
        try:
            if self._expect is None:
                return
            try:
                self.quit()
            except (TypeError, AttributeError):
                pass

            # The following programs around a bug in pexpect.
            def dummy(): pass
            try:
                self._expect.close = dummy
            except Exception as msg:
                pass
        except Exception as msg:
            pass

    def quit(self, verbose=False, timeout=0.25):
        """
        EXAMPLES::

            sage: a = maxima('y')
            sage: maxima.quit()
            sage: a._check_valid()
            Traceback (most recent call last):
            ...
            ValueError: The maxima session in which this object was defined is no longer running.
        """
        self._session_number += 1
        try:
            del self.__local_tmpfile
        except AttributeError:
            pass
        if self._expect is None:
            return
        # Send a kill -9 to the process *group*.
        # this is *very useful* when external binaries are started up
        # by shell scripts, and killing the shell script doesn't
        # kill the binary.
        E = self._expect
        if verbose:
            if self.is_remote():
                    print "Exiting spawned %s process (local pid=%s, running on %s)"%(self,E.pid,self._server)
            else:
                print "Exiting spawned %s process."%self
        try:
#            TO DO: This should be implemented or the remote tmp will get crowded
#            self._remove_remote_tmpfile()
            E.sendline(self._quit_string())
            self._so_far(wait=timeout)
            # In case of is_remote(), killing the local "ssh -t" also kills the remote process it initiated
            os.killpg(E.pid, 9)
            os.kill(E.pid, 9)
        except (RuntimeError, OSError) as msg:
            pass
        self._expect = None
        return

    def _quit_string(self):
        return 'quit'

    def _local_tmpfile(self):
        """
        Return a filename that is used to buffer long command lines for this interface

        INPUT:

        ``e`` -- an expect interface instance

        OUTPUT:

        A string that provides a temporary filename and is unique for the
        given interface.

        TEST:

        The filename is cached::

            sage: gap._local_tmpfile() is gap._local_tmpfile()
            True

        The following two problems were fixed in #10004.

        1. Different interfaces have different temp-files::

            sage: gap._local_tmpfile() != singular._local_tmpfile()
            True

        2. Interface instances in different branches of a parallelised
           function have different temp-files::

            sage: @parallel
            ... def f(n):
            ...     return gap._local_tmpfile()
            ...
            sage: L = [t[1] for t in f(range(5))]
            sage: len(set(L))
            5

        The following used to fail::

            sage: s = gap._local_tmpfile()
            sage: L = [t[1] for t in f(range(5))]
            sage: len(set(L))
            5

        AUTHOR:

        - Simon King (2010-09): Making the tmp-file unique for the interface instance

        """
        try:
            return self.__local_tmpfile
        except AttributeError:
            self.__local_tmpfile = os.path.join(SAGE_TMP_INTERFACE, 'tmp' + str(self.pid()))
            return self.__local_tmpfile

    def _remote_tmpdir(self):
        return self.__remote_tmpdir

    def _remote_tmpfile(self):
        try:
            return self.__remote_tmpfile
        except AttributeError:
            self.__remote_tmpfile = self._remote_tmpdir()+"/interface_%s:%s"%(LOCAL_IDENTIFIER,self.pid())
            return self.__remote_tmpfile

    def _send_tmpfile_to_server(self, local_file=None, remote_file=None):
        if local_file is None:
            local_file = self._local_tmpfile()
        if remote_file is None:
            remote_file = self._remote_tmpfile()
        cmd = 'scp "%s" %s:"%s" 1>&2 2>/dev/null'%(local_file, self._server, remote_file)
#        print cmd
        os.system(cmd)

    def _get_tmpfile_from_server(self, local_file=None,remote_file=None):
        if local_file is None:
            local_file = self._local_tmpfile()
        if remote_file is None:
            remote_file = self._remote_tmpfile()
        cmd = 'scp %s:"%s" "%s" 1>&2 2>/dev/null'%( self._server, remote_file, local_file)
#        print cmd
        os.system(cmd)

    def _remove_tmpfile_from_server(self):
        if not (self.__remote_tmpfile is None):
            raise NotImplementedError

    def read(self, filename):
        r"""
        EXAMPLES::

            sage: filename = tmp_filename()
            sage: f = open(filename, 'w')
            sage: f.write('x = 2\n')
            sage: f.close()
            sage: octave.read(filename)  # optional - octave
            sage: octave.get('x')        #optional
            ' 2'
            sage: import os
            sage: os.unlink(filename)
        """
        self.eval(self._read_in_file_command(filename))

    def _read_in_file_command(self, filename):
        raise NotImplementedError

    def _eval_line_using_file(self, line, restart_if_needed=True):
        """
        Evaluate a line of commands, using a temporary file.

        REMARK:

        By default, this is called when a long command is
        evaluated in :meth:`eval`.

        If the command can not be evaluated since the interface
        has crashed, it is automatically restarted and tried
        again *once*.

        INPUT:

        - ``line`` -- (string) a command.
        - ``restart_if_needed`` - (optional bool, default ``True``) --
          If it is ``True``, the command evaluation is evaluated
          a second time after restarting the interface, if an
          ``EOFError`` occured.

        TESTS::

            sage: singular._eval_line_using_file('def a=3;')
            '< "...";'
            sage: singular('a')
            3
            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line_using_file('def a=3;')
            Singular crashed -- automatically restarting.
            '< "...";'
            sage: singular('a')
            3
            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line_using_file('def a=3;', restart_if_needed=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Singular terminated unexpectedly while reading in a large line...

        We end by triggering a re-start of Singular, since otherwise
        the doc test of another method would fail by a side effect.
        ::

            sage: singular(3)
            Singular crashed -- automatically restarting.
            3

        """
        F = open(self._local_tmpfile(), 'w')
        F.write(line+'\n')
        F.close()
        tmp_to_use = self._local_tmpfile()
        if self.is_remote():
            self._send_tmpfile_to_server()
            tmp_to_use = self._remote_tmpfile()
        try:
            s = self._eval_line(self._read_in_file_command(tmp_to_use), allow_use_file=False, restart_if_needed=False)
        except pexpect.EOF as msg:
            if self._quit_string() in line:
                # we expect to get an EOF if we're quitting.
                return ''
            elif restart_if_needed==True: # the subprocess might have crashed
                try:
                    self._synchronize()
                    return self._post_process_from_file(self._eval_line_using_file(line, restart_if_needed=False))
                except RuntimeError as msg:
                    raise RuntimeError('%s terminated unexpectedly while reading in a large line:\n%s'%(self,msg[0]))
                except TypeError:
                    pass
            raise RuntimeError('%s terminated unexpectedly while reading in a large line'%self)
        except RuntimeError as msg:
            if self._quit_string() in line:
                if self._expect is None or not self._expect.isalive():
                    return ''
                raise
            if restart_if_needed==True and (self._expect is None or not self._expect.isalive()):
                try:
                    self._synchronize()
                    return self._post_process_from_file(self._eval_line_using_file(line, restart_if_needed=False))
                except TypeError:
                    pass
                except RuntimeError as msg:
                    raise RuntimeError('%s terminated unexpectedly while reading in a large line'%self)
            if "Input/output error" in msg[0]: # This occurs on non-linux machines
                raise RuntimeError('%s terminated unexpectedly while reading in a large line'%self)
            raise RuntimeError('%s terminated unexpectedly while reading in a large line:\n%s'%(self,msg[0]))
        return self._post_process_from_file(s)

    def _post_process_from_file(self, s):
        return s

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=True):
        """
        Evaluate a line of commands.

        REMARK:

        By default, a long command (length exceeding ``self._eval_using_file_cutoff``)
        is evaluated using :meth:`_eval_line_using_file`.

        If the command can not be evaluated since the interface
        has crashed, it is automatically restarted and tried
        again *once*.

        If the optional ``wait_for_prompt`` is ``False`` then even a very
        long line will not be evaluated by :meth:`_eval_line_using_file`,
        since this does not support the ``wait_for_prompt`` option.

        INPUT:

        - ``line`` -- (string) a command.
        - ``allow_use_file`` (optional bool, default ``True``) --
          allow to evaluate long commands using :meth:`_eval_line_using_file`.
        - ``wait_for_prompt`` (optional bool, default ``True``) --
          wait until the prompt appears in the sub-process' output.
        - ``restart_if_needed`` (optional bool, default ``True``) --
          If it is ``True``, the command evaluation is evaluated
          a second time after restarting the interface, if an
          ``EOFError`` occured.

        TESTS::

            sage: singular._eval_line('def a=3;')
            'def a=3;'
            sage: singular('a')
            3
            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line('def a=3;')
            Singular crashed -- automatically restarting.
            'def a=3;'
            sage: singular('a')
            3
            sage: singular.eval('kill a')
            'kill a;'

        We are now sending a command that would run forever. But since
        we declare that we are not waiting for a prompt, we can interrupt
        it without a KeyboardInterrupt. At the same time, we test that
        the line is not forwarded to :meth:`_eval_line_using_file`, since
        that method would not support the ``wait_for_prompt`` option.
        For reasons which are currently not understood, the ``interrupt``
        test usually returns immediately, but sometimes it takes a very
        long time on the same system. ::

            sage: cutoff = singular._eval_using_file_cutoff
            sage: singular._eval_using_file_cutoff = 4
            sage: singular._eval_line('for(int i=1;i<=3;i++){i=1;};', wait_for_prompt=False)
            ''
            sage: singular.interrupt(timeout=3)  # sometimes very slow (up to 60s on sage.math, 2012)
            False
            sage: singular._eval_using_file_cutoff = cutoff

        Last, we demonstrate that by default the execution of a command
        is tried twice if it fails the first time due to a crashed
        interface::

            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line_using_file('def a=3;', restart_if_needed=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Singular terminated unexpectedly while reading in a large line...

        Since the test of the next method would fail, we re-start
        Singular now. ::

            sage: singular(2+3)
            Singular crashed -- automatically restarting.
            5

        """
        if allow_use_file and wait_for_prompt and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line)
        try:
            if self._expect is None:
                self._start()
            E = self._expect
            try:
                if len(line) >= 4096:
                    raise RuntimeError("Sending more than 4096 characters with %s on a line may cause a hang and you're sending %s characters"%(self, len(line)))
                E.sendline(line)
                if wait_for_prompt == False:
                    return ''

            except OSError as msg:
                if restart_if_needed:
                    # The subprocess most likely crashed.
                    # If it's really still alive, we fall through
                    # and raise RuntimeError.
                    if sys.platform.startswith('sunos'):
                        # On (Open)Solaris, we might need to wait a
                        # while because the process might not die
                        # immediately. See Trac #14371.
                        for t in [0.5, 1.0, 2.0]:
                            if E.isalive():
                                time.sleep(t)
                            else:
                                break
                    if not E.isalive():
                        try:
                            self._synchronize()
                        except (TypeError, RuntimeError):
                            pass
                        return self._eval_line(line,allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt, restart_if_needed=False)
                raise RuntimeError, "%s\nError evaluating %s in %s"%(msg, line, self), sys.exc_info()[2]

            if len(line)>0:
                try:
                    if isinstance(wait_for_prompt, basestring):
                        E.expect(wait_for_prompt)
                    else:
                        E.expect(self._prompt)
                except pexpect.EOF as msg:
                    try:
                        if self.is_local():
                            tmp_to_use = self._local_tmpfile()
                        else:
                            tmp_to_use = self._remote_tmpfile()
                        if self._read_in_file_command(tmp_to_use) in line:
                            raise pexpect.EOF(msg)
                    except NotImplementedError:
                        pass
                    if self._quit_string() in line:
                        # we expect to get an EOF if we're quitting.
                        return ''
                    elif restart_if_needed==True: # the subprocess might have crashed
                        try:
                            self._synchronize()
                            return self._eval_line(line,allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt, restart_if_needed=False)
                        except (TypeError, RuntimeError):
                            pass
                    raise RuntimeError("%s\n%s crashed executing %s"%(msg,self, line))
                if self._terminal_echo:
                    out = E.before
                else:
                    out = E.before.rstrip('\n\r')
                    if out == '':   # match bug with echo
                        out = line
            else:
                if self._terminal_echo:
                    out = '\n\r'
                else:
                    out = ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt("Ctrl-c pressed while running %s"%self)
        if self._terminal_echo:
            i = out.find("\n")
            j = out.rfind("\r")
            return out[i+1:j].replace('\r\n','\n')
        else:
            return out.replace('\r\n','\n')

    def _keyboard_interrupt(self):
        print "Interrupting %s..."%self
        if self._restart_on_ctrlc:
            try:
                self._expect.close(force=1)
            except pexpect.ExceptionPexpect as msg:
                raise pexpect.ExceptionPexpect( "THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
            self._start()
            raise KeyboardInterrupt("Restarting %s (WARNING: all variables defined in previous session are now invalid)"%self)
        else:
            self._expect.sendline(chr(3))  # send ctrl-c
            self._expect.expect(self._prompt)
            self._expect.expect(self._prompt)
            raise KeyboardInterrupt("Ctrl-c pressed while running %s"%self)

    def interrupt(self, tries=20, timeout=0.3, quit_on_fail=True):
        E = self._expect
        if E is None:
            return True
        success = False
        try:
            for i in range(tries):
                E.sendline(chr(3))
                E.sendline(self._quit_string())
                try:
                    E.expect(self._prompt, timeout=timeout)
                    success= True
                    break
                except (pexpect.TIMEOUT, pexpect.EOF) as msg:
                    #print msg
                    pass
        except Exception as msg:
            pass
        if success:
            pass
        elif quit_on_fail:
            self.quit()
        return success

    ###########################################################################
    # BEGIN Synchronization code.
    ###########################################################################

    def _before(self):
        r"""
        Return the previous string that was sent through the interface.

        EXAMPLES::

            sage: singular(2+3)
            5
            sage: singular._before()
            '5\r\n'
        """
        return self._expect.before

    def _interrupt(self):
        for i in range(15):
            try:
                self._sendstr('quit;\n'+chr(3))
                self._expect_expr(timeout=2)
            except pexpect.TIMEOUT:
                pass
            except pexpect.EOF:
                self._crash_msg()
                self.quit()
            else:
                return

    def _expect_expr(self, expr=None, timeout=None):
        r"""
        Wait for a given expression expr (which could be a regular
        expression or list of regular expressions) to appear in the output
        for at most timeout seconds.

        Use ``r._expect.before`` to see what was put in the
        output stream before the expression.

        INPUT:


        -  ``expr`` - None or a string or list of strings
           (default: None)

        -  ``timeout`` - None or a number (default: None)


        EXAMPLES:

        We test all of this using the R interface. First we put
        10 + 15 in the input stream::

            sage: r._sendstr('abc <- 10 +15;\n')

        Here an exception is raised because 25 hasn't appeared yet in the
        output stream. The key thing is that this doesn't lock, but instead
        quickly raises an exception.

        ::

            sage: t = walltime()
            sage: try:
            ....:    r._expect_expr('25', timeout=0.5)
            ....: except Exception:
            ....:    print 'Did not get expression'
            Did not get expression

        A quick consistency check on the time that the above took::

            sage: w = walltime(t); w > 0.4 and w < 10
            True

        We tell R to print abc, which equals 25.

        ::

            sage: r._sendstr('abc;\n')

        Now 25 is in the output stream, so we can wait for it.

        ::

            sage: r._expect_expr('25')

        This gives us everything before the 25.

        ::

            sage: r._expect.before
            'abc;\r\n[1] '

        We test interrupting ``_expect_expr`` using the GP interface,
        see #6661.  Unfortunately, this test doesn't work reliably using
        Singular, see #9163 and the follow-up #10476.
        The ``gp.eval('0')`` in this test makes sure that ``gp`` is
        running, so a timeout of 1 second should be sufficient. ::

            sage: print sage0.eval("dummy=gp.eval('0'); alarm(1); gp._expect_expr('1')")  # long time
            Control-C pressed.  Interrupting PARI/GP interpreter. Please wait a few seconds...
            ...
            AlarmInterrupt:
        """
        if expr is None:
            # the following works around gap._prompt_wait not being defined
            expr = (hasattr(self,'_prompt_wait') and self._prompt_wait) or self._prompt
        if self._expect is None:
            self._start()
        try:
            if timeout:
                i = self._expect.expect(expr,timeout=timeout)
            else:
                i = self._expect.expect(expr)
            if i > 0:
                v = self._expect.before
                self.quit()
                raise ValueError("%s\nComputation failed due to a bug in %s -- NOTE: Had to restart."%(v, self))
        except KeyboardInterrupt as err:
            i = 0
            while True:
                try:
                    print "Control-C pressed.  Interrupting %s. Please wait a few seconds..."%self
                    self._sendstr('quit;\n'+chr(3))
                    self._sendstr('quit;\n'+chr(3))
                    self.interrupt()
                    self.interrupt()
                except KeyboardInterrupt:
                    i += 1
                    if i > 10:
                        break
                    pass
                else:
                    break
            raise err

    def _sendstr(self, str):
        r"""
        Send a string to the pexpect interface, autorestarting the expect
        interface if anything goes wrong.

        INPUT:


        -  ``str`` - a string


        EXAMPLES: We illustrate this function using the R interface::

            sage: r._sendstr('a <- 10;\n')
            sage: r.eval('a')
            '[1] 10'

        We illustrate using the singular interface::

            sage: singular._sendstr('int i = 5;')
            sage: singular('i')
            5
        """
        if self._expect is None:
            self._start()
        try:
            os.write(self._expect.child_fd, str)
        except OSError:
            self._crash_msg()
            self.quit()
            self._sendstr(str)

    def _crash_msg(self):
        r"""
        Show a message if the interface crashed.

        EXAMPLE::

            sage: singular._crash_msg()
            Singular crashed -- automatically restarting.

        ::

            sage: singular('2+3')
            5
            sage: singular._sendstr('quit;\n')   # make it so that singular appears to die.
            sage: singular('2+3')
            Singular crashed -- automatically restarting.
            5
        """
        print "%s crashed -- automatically restarting."%self

    def _synchronize(self, cmd='1+%s;\n'):
        """
        Synchronize pexpect interface.

        This put a random integer (plus one!) into the output stream, then
        waits for it, thus resynchronizing the stream. If the random
        integer doesn't appear within 1 second, the interface is sent
        interrupt signals.

        This way, even if you somehow left the interface in a busy state
        computing, calling _synchronize gets everything fixed.

        EXAMPLES: We observe nothing, just as it should be::

            sage: r._synchronize()

        TESTS: This illustrates a synchronization bug being fixed (thanks
        to Simon King and David Joyner for tracking this down)::

            sage: R.<x> = QQ[]; f = x^3 + x + 1;  g = x^3 - x - 1; r = f.resultant(g); gap(ZZ); singular(R)
            Integers
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
        """
        if self._expect is None:
            return
        rnd = randrange(2147483647)
        s = str(rnd+1)
        cmd = cmd%rnd
        self._sendstr(cmd)
        try:
            self._expect_expr(timeout=0.5)
            if not s in self._expect.before:
                self._expect_expr(s,timeout=0.5)
                self._expect_expr(timeout=0.5)
        except pexpect.TIMEOUT:
            self._interrupt()
        except pexpect.EOF:
            self._crash_msg()
            self.quit()

    ###########################################################################
    # END Synchronization code.
    ###########################################################################

    def eval(self, code, strip=True, synchronize=False, locals=None, allow_use_file=True,
             split_lines="nofile", **kwds):
        """
        INPUT:


        -  ``code``       -- text to evaluate

        -  ``strip``      -- bool; whether to strip output prompts,
                             etc. (ignored in the base class).

        - ``locals``      -- None (ignored); this is used for compatibility
                             with the Sage notebook's generic system interface.

        - ``allow_use_file`` -- bool (default: True); if True and ``code`` exceeds an
                                interface-specific threshold then ``code`` will be communicated
                                via a temporary file rather that the character-based interface.
                                If False then the code will be communicated via the character interface.

        - ``split_lines`` -- Tri-state (default: "nofile"); if "nofile" then ``code`` is sent line by line
                             unless it gets communicated via a temporary file.
                             If True then ``code`` is sent line by line, but some lines individually
                             might be sent via temporary file. Depending on the interface, this may transform
                             grammatical ``code`` into ungrammatical input.
                             If False, then the whole block of code is evaluated all at once.

        -  ``**kwds``     -- All other arguments are passed onto the _eval_line
                             method. An often useful example is reformat=False.
        """
        if synchronize:
            try:
                self._synchronize()
            except AttributeError:
                pass

        if strip:
            try:
                code = self._strip_prompt(code)
            except AttributeError:
                pass

        if not isinstance(code, basestring):
            raise TypeError('input code must be a string.')

        #Remove extra whitespace
        code = code.strip()

        try:
            with gc_disabled():
                if (split_lines is "nofile" and allow_use_file and
                        self._eval_using_file_cutoff and len(code) > self._eval_using_file_cutoff):
                    return self._eval_line_using_file(code)
                elif split_lines:
                    return '\n'.join([self._eval_line(L, allow_use_file=allow_use_file, **kwds)
                                        for L in code.split('\n') if L != ''])
                else:
                    return self._eval_line(code, allow_use_file=allow_use_file, **kwds)
        # DO NOT CATCH KeyboardInterrupt, as it is being caught
        # by _eval_line
        # In particular, do NOT call self._keyboard_interrupt()
        except TypeError as s:
            raise TypeError('error evaluating "%s":\n%s'%(code,s))

    ############################################################
    #         Functions for working with variables.
    #  The first three must be overloaded by derived classes,
    #  and the definition depends a lot on the class.  But
    #  the functionality one gets from this is very nice.
    ############################################################

    def _object_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect._object_class(maxima)
            <class 'sage.interfaces.expect.ExpectElement'>
        """
        return ExpectElement

    def _function_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect._function_class(maxima)
            <class 'sage.interfaces.expect.ExpectFunction'>
        """
        return ExpectFunction

    def _function_element_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect._function_element_class(maxima)
            <class 'sage.interfaces.expect.FunctionElement'>
        """
        return FunctionElement


class ExpectFunction(InterfaceFunction):
    """
    Expect function.
    """
    pass


class FunctionElement(InterfaceFunctionElement):
    """
    Expect function element.
    """
    pass


def is_ExpectElement(x):
    return isinstance(x, ExpectElement)

class ExpectElement(InterfaceElement):
    """
    Expect element.
    """
    def __init__(self, parent, value, is_name=False, name=None):
        RingElement.__init__(self, parent)
        self._create = value
        if parent is None: return     # means "invalid element"
        # idea: Joe Wetherell -- try to find out if the output
        # is too long and if so get it using file, otherwise
        # don't.
        if isinstance(value, basestring) and parent._eval_using_file_cutoff and \
           parent._eval_using_file_cutoff < len(value):
            self._get_using_file = True

        if is_name:
            self._name = value
        else:
            try:
                self._name = parent._create(value, name=name)
            # Convert ValueError and RuntimeError to TypeError for
            # coercion to work properly.
            except (RuntimeError, ValueError) as x:
                self._session_number = -1
                raise TypeError, x, sys.exc_info()[2]
            except BaseException:
                self._session_number = -1
                raise
        self._session_number = parent._session_number

    def __hash__(self):
        """
        Returns the hash of self. This is a default implementation of hash
        which just takes the hash of the string of self.
        """
        return hash('%s%s'%(self, self._session_number))


    def _check_valid(self):
        """
        Check that this object is valid, i.e., the session in which this
        object is defined is still running. This is relevant for
        interpreters that can't be interrupted via ctrl-C, hence get
        restarted.
        """
        try:
            P = self.parent()
            if P is None or P._session_number == BAD_SESSION or self._session_number == -1 or \
                          P._session_number != self._session_number:
                raise ValueError("The %s session in which this object was defined is no longer running."%P.name())
        except AttributeError:
            raise ValueError("The session in which this object was defined is no longer running.")
        return P

    def __del__(self):
        try:
            self._check_valid()
        except ValueError:
            return
        try:
            if hasattr(self,'_name'):
                P = self.parent()
                if P is not None:
                    P.clear(self._name)

        except (RuntimeError, ExceptionPexpect) as msg:    # needed to avoid infinite loops in some rare cases
            #print msg
            pass

#    def _sage_repr(self):
#TO DO: this could use file transfers when self.is_remote()


class StdOutContext:
    """
    A context in which all communation between Sage and a subprocess
    interfaced via pexpect is printed to stdout.
    """
    def __init__(self, interface, silent=False, stdout=None):
        """
        Construct a new context in which all communation between Sage
        and a subprocess interfaced via pexpect is printed to stdout.

        INPUT:

        - ``interface`` - the interface whose communcation shall be dumped.

        - ``silent`` - if ``True`` this context does nothing

        - ``stdout`` - optional parameter for alternative stdout device (default: ``None``)

        EXAMPLE::

            sage: from sage.interfaces.expect import StdOutContext
            sage: with StdOutContext(gp):
            ...       gp('1+1')
            ...
            sage=...
        """
        self.interface = interface
        self.silent = silent
        self.stdout = stdout if stdout else sys.stdout

    def __enter__(self):
        """
        EXAMPLE::

            sage: from sage.interfaces.expect import StdOutContext
            sage: with StdOutContext(singular):
            ...       singular.eval('1+1')
            ...
            1+1;
            ...
        """
        if self.silent:
            return
        if self.interface._expect is None:
            self.interface._start()
        self._logfile_backup = self.interface._expect.logfile
        if self.interface._expect.logfile:
            self.interface._expect.logfile = Multiplex(self.interface._expect.logfile, self.stdout)
        else:
            self.interface._expect.logfile = Multiplex(self.stdout)

    def __exit__(self, typ, value, tb):
        """
        EXAMPLE::

            sage: from sage.interfaces.expect import StdOutContext
            sage: with StdOutContext(gap):
            ...       gap('1+1')
            ...
            $sage...
        """
        if self.silent:
            return
        self.interface._expect.logfile.flush()
        self.stdout.write("\n")
        self.interface._expect.logfile = self._logfile_backup

import os
def console(cmd):
    os.system(cmd)



