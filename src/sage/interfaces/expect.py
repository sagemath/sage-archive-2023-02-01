"""
Common Interface Functionality

See the examples in the other sections for how to use specific
interfaces. The interface classes all derive from the generic
interface that is described in this section.

AUTHORS:

- William Stein (2005): initial version

- William Stein (2006-03-01): got rid of infinite loop on startup if
  client system missing
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

from __future__ import with_statement

import os
import weakref
import time
import gc
from random import randrange

########################################################
# Important note: We use Pexpect 2.0 *not* Pexpect 2.1.
# For reasons I don't understand, pexpect2.1 is much
# worse than pexpect 2.0 for everything SAGE does.
########################################################
import pexpect
from pexpect import ExceptionPexpect


from sage.structure.sage_object import SageObject
from sage.structure.parent_base import ParentWithBase
import  sage.structure.element

import sage.misc.sage_eval

import quit

import cleaner

import os

from sage.misc.misc import SAGE_ROOT, verbose, SAGE_TMP_INTERFACE, LOCAL_IDENTIFIER
from sage.structure.element import RingElement
BAD_SESSION = -2

failed_to_start = []

#tmp_expect_interface_local='%s/tmp'%SAGE_TMP_INTERFACE

def tmp_expect_interface_local():
    return '%s/tmp'%SAGE_TMP_INTERFACE + str(os.getpid())

## On some platforms, e.g., windows, this can easily take 10 seconds!?!  Terrible.  And
## it should not be necessary or used anyways.
## def _absolute(cmd):
##     c = cmd.split()
##     s  = c[0]
##     t = os.popen('which %s'%s).read().strip()
##     if len(t) == 0:
##         raise RuntimeError
##     return ' '.join([t] + c[1:])

# The subprocess is a shared resource.  In a multi-threaded
# environment, there would have to be a lock to control access to the
# subprocess.  Fortunately, SAGE does not use Python threads.
# Unfortunately, the combination of the garbage collector and __del__
# methods gives rise to the same issues.  So, in places where we need
# to do a sequence of operations on the subprocess and make sure
# nothing else intervenes (for example, when we write a command and
# then read back the result) we need to disable the garbage collector.
# See TRAC #955 for a more detailed description of this problem.

# This class is intended to be used with the "with" statement found
# in Python 2.5 and above.  To use it, add the following line at the top
# of your file:
#   from __future__ import with_statement
# Then to turn off the garbage collector for a particular region of code,
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

class AsciiArtString(str):
    def __init__(self, x):
        str.__init__(self, x)

    def __repr__(self):
        return str(self)

class PropTypeError(Exception):
    pass


class Expect(ParentWithBase):
    """
    Expect interface object.
    """
    def __init__(self, name, prompt, command=None, server=None, server_tmpdir=None,
                 ulimit = None, maxread=100000,
                 script_subdirectory="", restart_on_ctrlc=False,
                 verbose_start=False, init_code=[], max_startup_time=30,
                 logfile = None, eval_using_file_cutoff=0,
                 do_cleaner = True, remote_cleaner = False, path=None):

        self.__is_remote = False
        self.__remote_cleaner = remote_cleaner
        if command == None:
            command = name
        if not server is None:
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
        self.__name = name
        self.__coerce_name = '_' + name.lower() + '_'
        self.__command = command
        self._prompt = prompt
        self._restart_on_ctrlc = restart_on_ctrlc
        self.__verbose_start = verbose_start
        if not path is None:
            self.__path = path
        elif script_subdirectory is None:
            self.__path = '.'
        else:
            self.__path = '%s/data/extcode/%s/%s'%(SAGE_ROOT,self.__name,
                                                   self.__script_subdirectory)
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
        ParentWithBase.__init__(self, self)

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
        except pexpect.TIMEOUT, msg:
            return False, E.before
        except pexpect.EOF, msg:
            return True, E.before
        except Exception, msg:   # weird major problem!
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
                #if not new is None:
                X = self.__so_far + new
                del self.__so_far
                return True, X, new
            #new = self._expect.before
            try:
                self.__so_far += new
            except (AttributeError, TypeError):
                self.__so_far = new
            return False, self.__so_far, new
        except AttributeError, msg:   # no __so_far
            raise RuntimeError(msg)

    def is_remote(self):
        return self.__is_remote

    def is_local(self):
        return not self.__is_remote

    def user_dir(self):
        return self.__path

    def _repr_(self):
        return self.__name.capitalize()

    def _change_prompt(self, prompt):
        self._prompt = prompt

#    (pdehaye 20070819: this was used by some interfaces but does not work well remotely)
#    def _temp_file(self, x):
#        T = self.__path + "/tmp/"
#        if not os.path.exists(T):
#            os.makedirs(T)
#        return T + str(x)

    def name(self, new_name=None):
        return self.__name

    def path(self):
        return self.__path

    def expect(self):
        if self._expect is None:
            self._start()
        return self._expect

    def interact(self):
        r"""
        This allows you to interactively interact with the child
        interpreter. Press Ctrl-D or type 'quit' or 'exit' to exit and
        return to Sage.

        .. note::

           This is completely different than the console() member
           function. The console function opens a new copy of the
           child interepreter, whereas the interact function gives you
           interactive access to the interpreter that is being used by
           Sage. Use sage(xxx) or interpretername(xxx) to pull objects
           in from sage to the interpreter.
        """
        if self._expect is None:
            self._start()
        import sage.misc.preparser_ipython
        sage.misc.preparser_ipython.switch_interface_general(self)

    def _pre_interact(self):
        pass

    def _post_interact(self):
        pass

    def pid(self):
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
(WARNING: this would overwrite your current list of auhorized keys on "remote")

In many cases, the server that can actually run "slave" is not accessible from the internet directly, but you have to hop through an intermediate trusted server, say "gate".
If that is your case, get help with _install_hints_ssh_through_gate().

"""

    def _install_hints_ssh_through_gate(self):
        r"""
        Hints for installing passwordless authentication through a gate
        """
        # Written by Paul-Olivier Dehaye 2007/08/23
        return """

 We assume you would like to run a "slave" process  on a machine called "remote" from a machine running SAGE called "local". We also assume "remote" can only be accessed from "local" by ssh'ing first to "gate" (this is a fairly common setup). Sometimes, "gate" and "remote" haved a shared filesystem, and this helps a bit.

  Note: You cannot just create shell scripts on "local" and "gate" that would use two successive SSH connections to "remote" in order to simulate running "slave" locally. This is because SAGE will sometimes use files (and scp)  to communicate with "remote", which shell scripts would not take care of.

You need to setup:
 * passwordless authentication to "gate" from "local"
 * add passwordless authentication to "remote" from "local",
   for instance by appending the file local:~/.ssh/id_rsa.pub to remote:~/.ssh/authorized_keys2 and logging in once
      (this is only needed if "remote" and "gate" don\'t share fylesystems)
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
        self.quit()  # in case one is already running
        global failed_to_start

        self._session_number += 1
        current_path = os.path.abspath('.')
        dir = self.__path
        if not os.path.exists(dir):
            os.makedirs(dir)
        os.chdir(dir)

        #If the 'SAGE_PEXPECT_LOG' environment variable is set and
        #the current logfile is None, then set the logfile to be one
        #in .sage/pexpect_logs/
        if self.__logfile is None and 'SAGE_PEXPECT_LOG' in os.environ:
            from sage.misc.all import DOT_SAGE
            logs = '%s/pexpect_logs'%DOT_SAGE
            if not os.path.exists(logs):
                os.mkdir(logs)

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
            self._expect = pexpect.spawn(cmd, logfile=self.__logfile)
            if self._do_cleaner():
                cleaner.cleaner(self._expect.pid, cmd)

        except (ExceptionPexpect, pexpect.EOF, IndexError):
            self._expect = None
            self._session_number = BAD_SESSION
            failed_to_start.append(self.__name)
            raise RuntimeError, "Unable to start %s because the command '%s' failed.\n%s"%(
                self.__name, cmd, self._install_hints())

        os.chdir(current_path)
        self._expect.timeout = self.__max_startup_time
        #self._expect.setmaxread(self.__maxread)
        self._expect.maxread = self.__maxread
        self._expect.delaybeforesend = 0
        try:
            self._expect.expect(self._prompt)
        except (pexpect.TIMEOUT, pexpect.EOF), msg:
            self._expect = None
            self._session_number = BAD_SESSION
            failed_to_start.append(self.__name)
            raise RuntimeError, "Unable to start %s"%self.__name
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
            except Exception, msg:
                pass
        except Exception, msg:
            pass

    def cputime(self):
        """
        CPU time since this process started running.
        """
        raise NotImplementedError

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
        except (RuntimeError, OSError), msg:
            pass
        self._expect = None
        return

    def _quit_string(self):
        return 'quit'

    def _local_tmpfile(self):
        try:
            return self.__local_tmpfile
        except AttributeError:
            self.__local_tmpfile = tmp_expect_interface_local()
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
            sage: octave.read(filename)  #optional -- requires Octave
            sage: octave.get('x')        #optional
            ' 2'
            sage: import os
            sage: os.unlink(filename)
        """
        self.eval(self._read_in_file_command(filename))

    def _read_in_file_command(self, filename):
        raise NotImplementedError

    def _eval_line_using_file(self, line):
        F = open(self._local_tmpfile(), 'w')
        F.write(line+'\n')
        F.close()
        tmp_to_use = self._local_tmpfile()
        if self.is_remote():
            self._send_tmpfile_to_server()
            tmp_to_use = self._remote_tmpfile()
        try:
            s = self._eval_line(self._read_in_file_command(tmp_to_use), allow_use_file=False)
        except pexpect.EOF, msg:
            if self._quit_string() in line:
                # we expect to get an EOF if we're quitting.
                return ''
            raise RuntimeError, '%s terminated unexpectedly while reading in a large line'%self
        return self._post_process_from_file(s)

    def _post_process_from_file(self, s):
        return s

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True):
        if allow_use_file and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line)
        try:
            if self._expect is None:
                self._start()
            E = self._expect
            try:
                if len(line) >= 4096:
                    raise RuntimeError, "Sending more than 4096 characters with %s on a line may cause a hang and you're sending %s characters"%(self, len(line))
                E.sendline(line)
                if wait_for_prompt == False:
                    return ''

            except OSError, msg:
                raise RuntimeError, "%s\nError evaluating %s in %s"%(msg, line, self)

            if len(line)>0:
                try:
                    if isinstance(wait_for_prompt, basestring):
                        E.expect(wait_for_prompt)
                    else:
                        E.expect(self._prompt)
                except pexpect.EOF, msg:
                    try:
                        if self.is_local():
                            tmp_to_use = self._local_tmpfile()
                        else:
                            tmp_to_use = self._remote_tmpfile()
                        if self._read_in_file_command(tmp_to_use) in line:
                            raise pexpect.EOF, msg
                    except NotImplementedError:
                        pass
                    if self._quit_string() in line:
                        # we expect to get an EOF if we're quitting.
                        return ''
                    raise RuntimeError, "%s\n%s crashed executing %s"%(msg,self, line)
                out = E.before
            else:
                out = '\n\r'
        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt, "Ctrl-c pressed while running %s"%self
        i = out.find("\n")
        j = out.rfind("\r")
        return out[i+1:j].replace('\r\n','\n')

    def _keyboard_interrupt(self):
        print "Interrupting %s..."%self
        if self._restart_on_ctrlc:
            try:
                self._expect.close(force=1)
            except pexpect.ExceptionPexpect, msg:
                raise pexpect.ExceptionPexpect( "THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
            self._start()
            raise KeyboardInterrupt, "Restarting %s (WARNING: all variables defined in previous session are now invalid)"%self
        else:
            self._expect.sendline(chr(3))  # send ctrl-c
            self._expect.expect(self._prompt)
            self._expect.expect(self._prompt)
            raise KeyboardInterrupt, "Ctrl-c pressed while running %s"%self

    def interrupt(self, tries=20, timeout=0.3, quit_on_fail=True):
        E = self._expect
        if E is None:
            return True
        success = False
        try:
            for i in range(tries):
                E.sendline(self._quit_string())
                E.sendline(chr(3))
                try:
                    E.expect(self._prompt, timeout=timeout)
                    success= True
                    break
                except (pexpect.TIMEOUT, pexpect.EOF), msg:
                    #print msg
                    pass
        except Exception, msg:
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
        Return the previous string that was send through the interface.

        EXAMPLES::

            sage: singular(2+3)
            5
            sage: singular._before()
            'print(sage...);\r\n5\r\n'
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


        EXAMPLES: We test all of this using the R interface. First we put
        10 + 15 in the input stream::

            sage: r._sendstr('abc <- 10 +15;\n')

        Here an exception is raised because 25 hasn't appeared yet in the
        output stream. The key thing is that this doesn't lock, but instead
        quickly raises an exception.

        ::

            sage: t = walltime()
            sage: try: r._expect_expr('25', timeout=0.5)
            ... except: print 'Did not get expression'
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
        """
        if expr is None:
            expr = self._prompt_wait
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
                raise ValueError, "%s\nComputation failed due to a bug in %s -- NOTE: Had to restart."%(v, self)
        except KeyboardInterrupt, msg:
            i = 0
            while True:
                try:
                    print "Control-C pressed.  Interrupting R. Please wait a few seconds..."
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
            raise KeyboardInterrupt, msg

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
            sage: import os
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

    def eval(self, code, strip=True, synchronize=False, locals=None, **kwds):
        """
        INPUT:


        -  ``code`` - text to evaluate

        -  ``strip`` - bool; whether to strip output prompts,
           etc. (ignored in the base class).

        - ``locals`` - None (ignored); this is used for compatibility with the
          Sage notebook's generic system interface.

        -  ``**kwds`` - All other arguments are passed onto
           the _eval_line method. An often useful example is
           reformat=False.
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
            raise TypeError, 'input code must be a string.'

        #Remove extra whitespace
        code = code.strip()

        try:
            with gc_disabled():
                return '\n'.join([self._eval_line(L, **kwds) for L in code.split('\n') if L != ''])
        except KeyboardInterrupt:
            # DO NOT CATCH KeyboardInterrupt, as it is being caught
            # by _eval_line
            # In particular, do NOT call self._keyboard_interrupt()
            raise
        except TypeError, s:
            raise TypeError, 'error evaluating "%s":\n%s'%(code,s)

    def execute(self, *args, **kwds):
        return self.eval(*args, **kwds)

    def __call__(self, x, name=None):

        r"""
        Create a new object in self from x.

        The object X returned can be used like any SAGE object, and
        wraps an object in self.  The standard arithmetic operators
        work.  Morever if foo is a function then
                      X.foo(y,z,...)
        calls foo(X, y, z, ...) and returns the corresponding object.

        EXAMPLES::

            sage: gp(2)
            2
            sage: gp('2')
            2
            sage: a = gp(2); gp(a) is a
            True

        """
        cls = self._object_class()

        #Handle the case when x is an object
        #in some interface.
        if isinstance(x, ExpectElement):
            if x.parent() is self:
                return x

            #We convert x into an object in this
            #interface by first going through Sage.
            try:
                return self(x._sage_())
            except (NotImplementedError, TypeError):
                pass

        if isinstance(x, basestring):
            return cls(self, x, name=name)
        try:
            return self._coerce_from_special_method(x)
        except TypeError:
            raise
        except AttributeError, msg:
            pass
        try:
            return self._coerce_impl(x, use_special=False)
        except TypeError, msg:
            try:
                return cls(self, str(x), name=name)
            except TypeError, msg2:
                raise TypeError, msg

    def _coerce_from_special_method(self, x):
        """
        Tries to coerce to self by calling a special underscore method.

        If no such method is defined, raises an AttributeError instead of a
        TypeError.
        """
        s = '_%s_'%self.name()
        if s == '_pari_':
            s = '_gp_'
        try:
            return (x.__getattribute__(s))(self)
        except AttributeError:
            return self(x._interface_init_())


    def _coerce_impl(self, x, use_special=True):
        if isinstance(x, (int, long)):
            import sage.rings.all
            return self(sage.rings.all.Integer(x))
        elif isinstance(x, float):
            import sage.rings.all
            return self(sage.rings.all.RDF(x))
        if use_special:
            try:
                return self._coerce_from_special_method(x)
            except AttributeError, msg:
                pass

        if isinstance(x, (list, tuple)):
            A = []
            z = []
            cls = self._object_class()
            for v in x:
                if isinstance(v, cls):
                    A.append(v.name())
                    z.append(v)
                else:
                    w = self(v)
                    A.append(w.name())
                    z.append(w)
            X = ','.join(A)
            r = self.new('%s%s%s'%(self._left_list_delim(), X, self._right_list_delim()))
            r.__sage_list = z   # do this to avoid having the entries of the list be garbage collected
            return r


        raise TypeError, "unable to coerce element into %s"%self.name()

    def new(self, code):
        return self(code)

    ###################################################################
    # these should all be appropriately overloaded by the derived class
    ###################################################################

    def _left_list_delim(self):
        return "["

    def _right_list_delim(self):
        return "]"

    def _assign_symbol(self):
        return "="

    def _equality_symbol(self):
        raise NotImplementedError

    # For efficiency purposes, you should definitely override these
    # in your derived class.
    def _true_symbol(self):
        try:
            return self.__true_symbol
        except AttributeError:
            self.__true_symbol = self.eval('1 %s 1'%self._equality_symbol())

    def _false_symbol(self):
        try:
            return self.__false_symbol
        except AttributeError:
            self.__false_symbol = self.eval('1 %s 2'%self._equality_symbol())

    def _lessthan_symbol(self):
        return '<'

    def _greaterthan_symbol(self):
        return '>'


    ############################################################
    #         Functions for working with variables.
    #  The first three must be overloaded by derived classes,
    #  and the definition depends a lot on the class.  But
    #  the functionality one gets from this is very nice.
    ############################################################

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s%s%s;'%(var,self._assign_symbol(), value)
        self.eval(cmd)

    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval(var)

    def get_using_file(self, var):
        r"""
        Return the string representation of the variable var in self,
        possibly using a file. Use this if var has a huge string
        representation, since it may be way faster.

        .. warning::

           In fact unless a special derived class implements this, it
           will *not* be any faster. This is the case for this class
           if you're reading it through introspection and seeing this.
        """
        return self.get(var)

    def clear(self, var):
        """
        Clear the variable named var.
        """
        self._available_vars.append(var)

    def _next_var_name(self):
        if len(self._available_vars) != 0:
            v = self._available_vars[0]
            del self._available_vars[0]
            return v
        self.__seq += 1
        return "sage%s"%self.__seq

    def _create(self, value, name=None):
        name = self._next_var_name() if name is None else name
        self.set(name, value)
        return name

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

    def _convert_args_kwds(self, args=None, kwds=None):
        """
        Converts all of the args and kwds to be elements of this
        interface.

        EXAMPLES::

            sage: args = [5]
            sage: kwds = {'x': 6}
            sage: args, kwds = gap._convert_args_kwds(args, kwds)
            sage: args
            [5]
            sage: map(type, args)
            [<class 'sage.interfaces.gap.GapElement'>]
            sage: type(kwds['x'])
            <class 'sage.interfaces.gap.GapElement'>
        """
        args = [] if args is None else args
        kwds = {} if kwds is None else kwds
        if not isinstance(args, list):
            args = [args]
        for i, arg in enumerate(args):
            if not isinstance(arg, ExpectElement) or arg.parent() is not self:
                args[i] = self(arg)
        for key, value in kwds.iteritems():
            if not isinstance(value, ExpectElement) or value.parent() is not self:
                kwds[key] = self(value)

        return args, kwds

    def _check_valid_function_name(self, function):
        """
        Checks to see if function is a valid function name in this
        interface. If it is not, an exception is raised. Otherwise, nothing
        is done.

        EXAMPLES::

            sage: gap._check_valid_function_name('SymmetricGroup')
            sage: gap._check_valid_function_name('')
            Traceback (most recent call last):
            ...
            ValueError: function name must be nonempty
            sage: gap._check_valid_function_name('__foo')
            Traceback (most recent call last):
            ...
            AttributeError
        """
        if function == '':
            raise ValueError, "function name must be nonempty"
        if function[:2] == "__":
            raise AttributeError

    def function_call(self, function, args=None, kwds=None):
        """
        EXAMPLES::

            sage: maxima.quad_qags(x, x, 0, 1, epsrel=1e-4)
            [0.5,5.5511151231257...E-15,21,0]
            sage: maxima.function_call('quad_qags', [x, x, 0, 1], {'epsrel':'1e-4'})
            [0.5,5.5511151231257...E-15,21,0]
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        return self.new("%s(%s)"%(function, ",".join([s.name() for s in args]+
                                                     ['%s=%s'%(key,value.name()) for key, value in kwds.items()])))

    def call(self, function_name, *args, **kwds):
        return self.function_call(function_name, args, kwds)

    def _contains(self, v1, v2):
        raise NotImplementedError

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return self._function_class()(self, attrname)

    def __cmp__(self, other):
        """
        Compare to pseudo-tty interfaces. To interfaces compare equal if
        and only if they are identical objects (this is a critical
        constrait so that caching of representations of objects in
        interfaces works correctly). Otherwise they are never equal.

        EXAMPLES::

            sage: sage.calculus.calculus.maxima == maxima
            False
            sage: maxima == maxima
            True
        """
        if self is other:
            return 0
        c = cmp(type(self), type(other))
        if c:
            return c
        return -1  # sucky, but I can't think of anything better; it is important that different interfaces to the same system still compare differently; unfortunately there is nothing to distinguish them.

    def console(self):
        raise NotImplementedError

    def help(self, s):
        return AsciiArtString('No help on %s available'%s)


class ExpectFunction(SageObject):
    """
    Expect function.
    """
    def __init__(self, parent, name):
        self._parent = parent
        self._name = name

    def __repr__(self):
        return "%s"%self._name

    def __call__(self, *args, **kwds):
        return self._parent.function_call(self._name, list(args), kwds)

    def _sage_doc_(self):
        """
        EXAMPLES::

            sage: gp.gcd._sage_doc_()
            'gcd(x,{y}): greatest common divisor of x and y.'
        """
        M = self._parent
        return M.help(self._name)




class FunctionElement(SageObject):
    """
    Expect function element.
    """
    def __init__(self, obj, name):
        self._obj = obj
        self._name = name

    def __repr__(self):
        return "%s"%self._name

    def __call__(self, *args, **kwds):
        return self._obj.parent().function_call(self._name, [self._obj] + list(args), kwds)

    def help(self):
        print self._sage_doc_()

    def _sage_doc_(self):
        """
        EXAMPLES::

            sage: gp(2).gcd._sage_doc_()
            'gcd(x,{y}): greatest common divisor of x and y.'
        """
        M = self._obj.parent()
        return M.help(self._name)


def is_ExpectElement(x):
    return isinstance(x, ExpectElement)



class ExpectElement(RingElement):
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
            except (TypeError, KeyboardInterrupt, RuntimeError, ValueError), x:
                self._session_number = -1
                raise TypeError, x
        self._session_number = parent._session_number

    def _latex_(self):
#        return "\\begin{verbatim}%s\\end{verbatim}"%self
        string = str(self)
        if not '|' in string:
            delim = '|'
        elif not '#' in string:
            delim = '#'
        elif not '@' in string:
            delim = '@'
        elif not '~' in string:
            delim = '~'
        return "\\verb%s%s%s"%(delim, string, delim)

    def __iter__(self):
        for i in range(1, len(self)+1):
            yield self[i]

    def __len__(self):
        raise NotImplementedError

    def __reduce__(self):
        return reduce_load, (self.parent(), self._reduce())

    def _reduce(self):
        return repr(self)

    def _r_action(self, x):   # used for coercion
        raise AttributeError

    def __call__(self, *args):
        self._check_valid()
        P = self.parent()
        return getattr(P, self.name())(*args)

    def __contains__(self, x):
        P = self._check_valid()
        if not isinstance(x, ExpectElement) or x.parent() is not self.parent():
            x = P.new(x)
        return P._contains(x.name(), self.name())


    def _sage_doc_(self):
        """
        EXAMPLES::

            sage: gp(2)._sage_doc_()
            '2'
        """
        return str(self)

    def __hash__(self):
        """
        Returns the hash of self. This is a defualt implementation of hash
        which just takes the hash of the string of self.
        """
        return hash('%s%s'%(self, self._session_number))

    def __cmp__(self, other):
        P = self.parent()
        if P.eval("%s %s %s"%(self.name(), P._equality_symbol(),
                                 other.name())) == P._true_symbol():
            return 0
        elif P.eval("%s %s %s"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
            return -1
        elif P.eval("%s %s %s"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
            return 1

        # everything is supposed to be comparable in Python, so we define
        # the comparison thus when no comparable in interfaced system.
        if (hash(self) < hash(other)):
            return -1
        else:
            return 1

    def _matrix_(self, R):
        raise NotImplementedError

    def _vector_(self, R):
        raise NotImplementedError

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
                raise ValueError, "The %s session in which this object was defined is no longer running."%P.name()
        except AttributeError:
            raise ValueError, "The session in which this object was defined is no longer running."
        return P

    def __del__(self):
        try:
            self._check_valid()
        except ValueError:
            return
        try:
            if hasattr(self,'_name'):
                P = self.parent()
                if not (P is None):
                    P.clear(self._name)

        except (RuntimeError, ExceptionPexpect), msg:    # needed to avoid infinite loops in some rare cases
            #print msg
            pass

    def _sage_(self):
        """
        Attempt to return a Sage version of this object.
        This is a generic routine that just tries to evaluate
        the repr(self).

        EXAMPLES::

            sage: gp(1/2)._sage_()
            1/2
            sage: _.parent()
            Rational Field
        """
        #TO DO: this could use file transfers when self.is_remote()
        try:
            return sage.misc.sage_eval.sage_eval(repr(self).replace('\n',''))
        except:
            raise NotImplementedError


    def sage(self):
        """
        Attempt to return a Sage version of this object.

        EXAMPLES::

            sage: gp(1/2).sage()
            1/2
            sage: _.parent()
            Rational Field
        """
        return self._sage_()

    def __repr__(self):
        try:
            self._check_valid()
        except ValueError:
            return '(invalid object -- defined in terms of closed session)'
        try:
            if self._get_using_file:
                s = self.parent().get_using_file(self._name)
        except AttributeError:
            s = self.parent().get(self._name)
        if s.__contains__(self._name):
            if hasattr(self, '__custom_name'):
                s =  s.replace(self._name, self.__dict__['__custom_name'])
        return s

    def __getattr__(self, attrname):
        P = self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        return P._function_element_class()(self, attrname)

    def get_using_file(self):
        """
        Return this element's string representation using a file. Use this
        if self has a huge string representation. It'll be way faster.

        EXAMPLES::

            sage: a = maxima(str(2^1000))
            sage: a.get_using_file()
            '10715086071862673209484250490600018105614048117055336074437503883703510511249361224931983788156958581275946729175531468251871452856923140435984577574698574803934567774824230985421074605062371141877954182153046474983581941267398767559165543946077062914571196477686542167660429831652624386837205668069376'
        """
        try:
            self._check_valid()
        except ValueError:
            return '(invalid object -- defined in terms of closed session)'
        return self.parent().get_using_file(self._name)

    def hasattr(self, attrname):
        """
        Returns whether the given attribute is already defined by this
        object, and in particular is not dynamically generated.

        EXAMPLES::

            sage: m = maxima('2')
            sage: m.hasattr('integral')
            True
            sage: m.hasattr('gcd')
            False
        """
        return not isinstance(getattr(self, attrname), FunctionElement)

    def attribute(self, attrname):
        """
        If this wraps the object x in the system, this returns the object
        x.attrname. This is useful for some systems that have object
        oriented attribute access notation.

        EXAMPLES::

            sage: g = gap('SO(1,4,7)')
            sage: k = g.InvariantQuadraticForm()
            sage: k.attribute('matrix')
            [ [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ],
              [ 0*Z(7), 0*Z(7), Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ] ]

        ::

            sage: e = gp('ellinit([0,-1,1,-10,-20])')
            sage: e.attribute('j')
            -122023936/161051
        """
        P = self._check_valid()
        return P('%s.%s'%(self.name(), attrname))

    def __getitem__(self, n):
        P = self._check_valid()
        if not isinstance(n, tuple):
            return P.new('%s[%s]'%(self._name, n))
        else:
            return P.new('%s[%s]'%(self._name, str(n)[1:-1]))

    def __int__(self):
        """
        EXAMPLES::

            sage: int(maxima('1'))
            1
            sage: type(_)
            <type 'int'>
        """
        return int(repr(self))

    def bool(self):
        P = self.parent()
        t = P._true_symbol()
        cmd = '%s %s %s'%(self._name, P._equality_symbol(), t)
        return P.eval(cmd) == t

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: bool(maxima(0))
            False
            sage: bool(maxima(1))
            True
        """
        return self.bool()

    def __long__(self):
        """
        EXAMPLES::

            sage: m = maxima('1')
            sage: long(m)
            1L
        """
        return long(repr(self))

    def __float__(self):
        """
        EXAMPLES::

            sage: m = maxima('1/2')
            sage: m.__float__()
            0.5
            sage: float(m)
            0.5
        """
        return float(repr(self))

    def _integer_(self, ZZ=None):
        """
        EXAMPLES::

            sage: m = maxima('1')
            sage: m._integer_()
            1
            sage: _.parent()
            Integer Ring
            sage: QQ(m)
            1
        """
        import sage.rings.all
        return sage.rings.all.Integer(repr(self))

    def _rational_(self):
        """
        EXAMPLES::

            sage: m = maxima('1/2')
            sage: m._rational_()
            1/2
            sage: _.parent()
            Rational Field
            sage: QQ(m)
            1/2
        """
        import sage.rings.all
        return sage.rings.all.Rational(repr(self))

    def name(self, new_name=None):
        """
        Returns the name of self. If new_name is passed in, then this
        function returns a new object identitical to self whose name is
        new_name.

        Note that this can overwrite existing variables in the system.

        EXAMPLES::

            sage: x = r([1,2,3]); x
            [1] 1 2 3
            sage: x.name()
            'sage3'
            sage: x = r([1,2,3]).name('x'); x
            [1] 1 2 3
            sage: x.name()
            'x'

        ::

            sage: s5 = gap.SymmetricGroup(5).name('s5')
            sage: s5
            SymmetricGroup( [ 1 .. 5 ] )
            sage: s5.name()
            's5'
        """
        if new_name is not None:
            if not isinstance(new_name, str):
                raise TypeError, "new_name must be a string"
            p = self.parent()
            p.set(new_name, self._name)
            return p._object_class()(p, new_name, is_name=True)

        return self._name

    def gen(self, n):
        P = self._check_valid()
        return P.new('%s.%s'%(self._name, int(n)))

    def _operation(self, operation, right):
        P = self._check_valid()
        try:
            return P.new('%s %s %s'%(self._name, operation, right._name))
        except Exception, msg:
            raise TypeError, msg

    def _add_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f + g
            sin(x)+cos(x)
            sage: f + 2
            cos(x)+2
            sage: 2 + f
            cos(x)+2
        """
        return self._operation("+", right)

    def _sub_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f - g
            cos(x)-sin(x)
            sage: f - 2
            cos(x)-2
            sage: 2 - f
            2-cos(x)
        """
        return self._operation('-', right)

    def _mul_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f*g
            cos(x)*sin(x)
            sage: 2*f
            2*cos(x)
        """
        return self._operation('*', right)

    def _div_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f/g
            cos(x)/sin(x)
            sage: f/2
            cos(x)/2
        """
        return self._operation("/", right)

    def __pow__(self, n):
        """
        EXAMPLES::

            sage: a = maxima('2')
            sage: a^(3/4)
            2^(3/4)
        """
        P = self._check_valid()
        if not hasattr(n, 'parent') or P is not n.parent():
            n = P(n)
        return self._operation("^", n)


def reduce_load(parent, x):
    return parent(x)

import os
def console(cmd):
    os.system(cmd)



