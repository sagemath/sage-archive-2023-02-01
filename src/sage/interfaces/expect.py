"""
Common Interface Functionality

See the examples in the other sections for how to use specific
interfaces.  The interface classes all derive from the generic
interface that is described in this section.

AUTHORS:
    -- William Stein (2005): initial version
    -- William Stein (2006-03-01): got rid of infinite loop on startup
                                   if client system missing
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
import weakref
import time

########################################################
# Important note: We use Pexpect 2.0 *not* Pexpect 2.1.
# For reasons I don't understand, pexpect2.1 is much
# worse than pexpect 2.0 for everything SAGE does.
########################################################
import pexpect

from sage.structure.sage_object import SageObject
from sage.structure.parent_base import ParentWithBase
import  sage.structure.element

import sage.misc.sage_eval

import quit

import monitor
EXPECT_MONITOR_INTERVAL=5  # kill any slave processes at most 5 seconds after parent dies.

from sage.misc.misc import SAGE_ROOT, verbose, SAGE_TMP_INTERFACE
from sage.structure.element import RingElement
BAD_SESSION = -2

failed_to_start = []

tmp='%s/tmp'%SAGE_TMP_INTERFACE

## On some platforms, e.g., windows, this can easily take 10 seconds!?!  Terrible.  And
## it should not be necessary or used anyways.
## def _absolute(cmd):
##     c = cmd.split()
##     s  = c[0]
##     t = os.popen('which %s'%s).read().strip()
##     if len(t) == 0:
##         raise RuntimeError
##     return ' '.join([t] + c[1:])


class Expect(ParentWithBase):
    """
    Expect interface object.
    """
    def __init__(self, name, prompt, command=None, server=None, maxread=100000,
                 script_subdirectory="", restart_on_ctrlc=False,
                 verbose_start=False, init_code=[], max_startup_time=30,
                 logfile = None, eval_using_file_cutoff=0,
                 do_monitor = True):

        self.__is_remote = False
        if command == None:
            command = name
        if server != None:
            command = "ssh -t %s %s"%(server, command)
            self.__is_remote = True
            eval_using_file_cutoff = 0  # don't allow this!
            #print command
            self._server = server
        self.__do_monitor = do_monitor
        self.__maxread = maxread
        self._eval_using_file_cutoff = eval_using_file_cutoff
        self.__script_subdirectory = script_subdirectory
        self.__name = name
        self.__coerce_name = '_' + name.lower() + '_'
        self.__command = command
        self._prompt = prompt
        self._restart_on_ctrlc = restart_on_ctrlc
        self.__verbose_start = verbose_start
        if script_subdirectory == None:
            self.__path = '.'
        else:
            self.__path = '%s/data/extcode/%s/%s'%(SAGE_ROOT,self.__name, self.__script_subdirectory)
        self.__initialized = False
        self.__seq = -1
        self._expect = None
        self._session_number = 0
        self.__init_code = init_code
        self.__max_startup_time = max_startup_time
        if isinstance(logfile, str):
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
        Return whether done and output so far and new output since last time called.
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
        except AttributeError:   # no __so_far
            raise RuntimeError

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

    def _temp_file(self, x):
        T = self.__path + "/tmp/"
        if not os.path.exists(T):
            os.makedirs(T)
        return T + str(x)

    def name(self):
        return self.__name

    def path(self):
        return self.__path

    def _continuation_prompt(self):
        return False

    def expect(self):
        if self._expect is None:
            self._start()
        return self._expect

    def interact(self):
        r"""
        This allows you to interactively interact with the child
        interpreter.  Press Ctrl-D or type 'quit' or 'exit'
        to exit and return to SAGE.

        \note{This is completely different than the console() member
        function.  The console function opens a new copy of the child
        interepreter, whereas the interact function gives you
        interactive access to the interpreter that is being used by
        SAGE.   Use sage(xxx) or interpretername(xxx) to pull objects
        in from sage to the interpreter.}
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
        """
        Hints for installing needed slave program on your computer.

        There are no hints by default.
        """
        return ''

    def _do_monitor(self):
        try:
            return self.__do_monitor
        except AttributeError:
            return False

    def _start(self, alt_message=None, block_during_init=True):
        self.quit()  # in case one is already running
        global failed_to_start
        if self.__name in failed_to_start:
            if alt_message:
                raise RuntimeError, alt_message
            else:
                raise RuntimeError, 'Unable to start %s (%s failed to start during this SAGE session; not attempting to start again)\n%s'%(self.__name, self.__name, self._install_hints())

        self._session_number += 1
        current_path = os.path.abspath('.')
        dir = self.__path
        if not os.path.exists(dir):
            os.makedirs(dir)
        os.chdir(dir)

        cmd = self.__command
##         try:
##             cmd = _absolute(self.__command)
##         except RuntimeError:
##             failed_to_start.append(self.__name)
##             raise RuntimeError, "%s\nCommand %s not available."%(
##                  self._install_hints(), self.__name)

        if self.__verbose_start:
            print cmd
            print "Starting %s"%cmd.split()[0]

        try:
            self._expect = pexpect.spawn(cmd, logfile=self.__logfile)
            if self._do_monitor() and not self.__is_remote:
                monitor.monitor(self._expect.pid, EXPECT_MONITOR_INTERVAL)

        except (pexpect.ExceptionPexpect, pexpect.EOF, IndexError):
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
            print 1
            self._expect.expect(self._prompt)
            print 2
        except (pexpect.TIMEOUT, pexpect.EOF), msg:
            print 3
            self._expect = None
            self._session_number = BAD_SESSION
            failed_to_start.append(self.__name)
            print msg
            raise RuntimeError, "Unable to start %s"%self.__name
        self._expect.timeout = None
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
                print msg
        except Exception, msg:
            print msg

    def cputime(self):
        """
        CPU time since this process started running.
        """
        raise NotImplementedError

    def quit(self, verbose=False):
        if self._expect is None:
            return
        # Send a kill -9 to the process *group*.
        # this is *very useful* when external binaries are started up
        # by shell scripts, and killing the shell script doesn't
        # kill the binary.
        if verbose:
            print "Exiting spawned %s process."%self
        try:
            self._expect.sendline(self._quit_string())
            self._so_far(wait=0.25)
            os.killpg(self._expect.pid, 9)
            os.kill(self._expect.pid, 9)
        except OSError, msg:
            pass
        try:
            self._expect.close(9)
        except Exception:
            pass
        self._expect = None

    def _quit_string(self):
        return 'quit'

    def _read_in_file_command(self, filename):
        raise NotImplementedError

    def _eval_line_using_file(self, line, tmp):
        F = open(tmp, 'w')
        F.write(line+'\n')
        F.close()
        try:
            s = self._eval_line(self._read_in_file_command(tmp), allow_use_file=False)
        except pexpect.EOF, msg:
            if self._quit_string() in line:
                # we expect to get an EOF if we're quitting.
                return ''
            raise RuntimeError, '%s terminated unexpectedly while reading in a large line'%self
        return self._post_process_from_file(s)

    def _post_process_from_file(self, s):
        return s

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True):
        #if line.find('\n') != -1:
        #    raise ValueError, "line must not contain any newlines"
        if allow_use_file and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line, tmp)
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
                    if isinstance(wait_for_prompt, str):
                        E.expect(wait_for_prompt)
                    else:
                        E.expect(self._prompt)
                except pexpect.EOF, msg:
                    try:
                        if self._read_in_file_command(tmp) in line:
                            raise pexpect.EOF, msg
                    except NotImplementedError:
                        pass
                    if self._quit_string() in line:
                        # we expect to get an EOF if we're quitting.
                        return ''
                    raise RuntimeError, "%s\n%s crashed executing %s"%(msg,
                                                   self, line)
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
            except pexpect.ExceptionPexpect:
                print "WARNING: -- unable to kill %s. You may have to do so manually."%self
                pass
            self._start()
            raise KeyboardInterrupt, "Restarting %s (WARNING: all variables defined in previous session are now invalid)"%self
        else:
            self._expect.sendline(chr(3))  # send ctrl-c
            self._expect.expect(self._prompt)
            self._expect.expect(self._prompt)
            raise KeyboardInterrupt, "Ctrl-c pressed while running %s"%self

    def eval(self, code, strip=True):
        """
        INPUT:
            code -- text to evaluate
            strip -- bool; whether to strip output prompts, etc.
                     (ignored in the base class).
        """
        code = str(code)
        code = code.strip()
        try:
            return '\n'.join([self._eval_line(L) for L in code.split('\n') if L != ''])
        except KeyboardInterrupt:
            self._keyboard_interrupt()
        except TypeError, s:
            return 'error evaluating "%s":\n%s'%(code,s)

    def execute(self, *args, **kwds):
        return self.eval(*args, **kwds)

    def __call__(self, x):
        r"""
        Create a new object in self from x.

        The object X returned can be used like any SAGE object, and
        wraps an object in self.  The standard arithmetic operators
        work.  Morever if foo is a function then
                      X.foo(y,z,...)
        calls foo(X, y, z, ...) and returns the corresponding object.
        """
        cls = self._object_class()

        if isinstance(x, cls) and x.parent() is self:
            return x
        if isinstance(x, str):
            return cls(self, x)
        try:
            return self._coerce_impl(x)
        except TypeError, msg:
            try:
                return cls(self, str(x))
            except TypeError, msg2:
                raise TypeError, msg


    def _coerce_impl(self, x):
        s = '_%s_'%self.name()
        if s == '_pari_':
            s = '_gp_'
        try:
            return (x.__getattribute__(s))(self)
        except AttributeError, msg:
            pass
        try:
            return x._interface_(self)
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

    def _create(self, value):
        name = self._next_var_name()
        #self._last_name = name
        self.set(name, value)
        return name

    #def _last_created_varname(self):
    #    return self._last_name

    def _object_class(self):
        return ExpectElement

    def function_call(self, function, args=[]):
        if function == '':
            raise ValueError, "function name must be nonempty"
        if function[:2] == "__":
            raise AttributeError
        if not isinstance(args, list):
            args = [args]
        for i in range(len(args)):
            if not isinstance(args[i], ExpectElement):
                args[i] = self.new(args[i])
        return self.new("%s(%s)"%(function, ",".join([s.name() for s in args])))

    def call(self, function_name, *args):
        return self.function_call(function_name, args)

    def _contains(self, v1, v2):
        raise NotImplementedError

    def _is_true_string(self, s):
        raise NotImplementedError

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return ExpectFunction(self, attrname)

    def __cmp__(self, other):
        return cmp(type(self), type(other))

    def console(self):
        raise NotImplementedError

    def help(self, s):
        print 'No help on %s available'%s


class ExpectFunction(SageObject):
    """
    Expect function.
    """
    def __init__(self, parent, name):
        self._parent = parent
        self._name = name

    def __repr__(self):
        return "%s"%self._name

    def __call__(self, *args):
        return self._parent.function_call(self._name, list(args))


class FunctionElement(SageObject):
    """
    Expect function element.
    """
    def __init__(self, obj, name):
        self._obj = obj
        self._name = name

    def __repr__(self):
        return "%s"%self._name

    def __call__(self, *args):
        return self._obj.parent().function_call(self._name, [self._obj] + list(args))

    def help(self):
        print self._sage_doc_()

    def _sage_doc_(self):
        return ''


def is_ExpectElement(x):
    return isinstance(x, ExpectElement)

class ExpectElement(RingElement):
    """
    Expect element.
    """
    def __init__(self, parent, value, is_name=False):
        RingElement.__init__(self, parent)
        self._create = value
        if parent is None: return     # means "invalid element"
        # idea: Joe Wetherell -- try to find out if the output
        # is too long aqnd if so get it using file, otherwise
        # don't.
        if isinstance(value, str) and parent._eval_using_file_cutoff and \
           parent._eval_using_file_cutoff < len(value):
            self._get_using_file = True

        if is_name:
            self._name = value
        else:
            try:
                self._name = parent._create(value)
            except (TypeError, KeyboardInterrupt, RuntimeError, ValueError), x:
                self._session_number = -1
                raise TypeError, x
        self._session_number = parent._session_number

    def _latex_(self):
        return "\\begin{verbatim}%s\\end{verbatim}"%self


    def __iter__(self):
        for i in range(1, len(self)+1):
            yield self[i]

    def __len__(self):
        raise NotImplementedError

    def __reduce__(self):
        return reduce_load, (self.parent(), self._reduce())

    def _reduce(self):
        return str(self)

    def _r_action(self, x):   # used for coercion
        raise AttributeError

    def __call__(self, *args):
        self._check_valid()
        P = self.parent()
        return getattr(P, self.name())(*args)

    def _sage_doc_(self):
        return str(self)

    def __cmp__(self, other):
        P = self.parent()
        if P.eval("%s %s %s"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
            return -1
        elif P.eval("%s %s %s"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
            return 1
        elif P.eval("%s %s %s"%(self.name(), P._equality_symbol(),
                                 other.name())) == P._true_symbol():
            return 0
        else:
            return -1  # everything is supposed to be comparable in Python, so we define
                       # the comparison thus when no comparable in interfaced system.

    def _matrix_(self, R):
        raise NotImplementedError

    def _vector_(self, R):
        raise NotImplementedError

    def _check_valid(self):
        """
        Check that this object is valid, i.e., the session in which
        this object is defined is still running.  This is relevant for
        interpreters that can't be interrupted via ctrl-C, hence get
        restarted.
        """
        try:
            P = self.parent()
            if P is None is None or P._session_number == BAD_SESSION or self._session_number == -1 or \
                          (P._restart_on_ctrlc and P._session_number != self._session_number):
                raise ValueError, "The %s session in which this object was defined is no longer running."%P.name()
        except AttributeError:
            raise ValueError, "The session in which this object was defined is no longer running."
        return P

    def __contains__(self, x):
        P = self._check_valid()
        if not isinstance(x, ExpectElement) or x.parent() != self.parent():
            x = P.new(x)
        t = P._contains(x.name(), self.name())
        return P._is_true_string(t)

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
        except RuntimeError, msg:    # needed to avoid infinite loops in some rare cases
            #print msg
            pass

    def _sage_(self):
        """
        Attempt to return a SAGE version of this object.
        """
        return self.sage()

    def sage(self):
        return sage.misc.sage_eval.sage_eval(str(self))

    def __repr__(self):
        try:
            self._check_valid()
        except ValueError:
            return '(invalid object -- defined in terms of closed session)'
        try:
            if self._get_using_file:
                return self.parent().get_using_file(self._name)
        except AttributeError:
            return self.parent().get(self._name)

    def __getattr__(self, attrname):
        self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        return FunctionElement(self, attrname)

    def hasattr(self, attrname):
        """
        Returns whether the given attribute is already defined by this object,
        and in particular is not dynamically generated.
        """
        return not isinstance(getattr(self, attrname), FunctionElement)

    def attribute(self, attrname):
        """
        If this wraps the object x in the system, this returns the object
        x.attrname.  This is useful for some systems that have object
        oriented attribute access notation.

        EXAMPLES:
            sage: g = gap('SO(1,4,7)')
            sage: k = g.InvariantQuadraticForm()
            sage: k.attribute('matrix')
            [ [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ],
              [ 0*Z(7), 0*Z(7), Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ] ]

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
        return int(str(self))

    def bool(self):
        P = self.parent()
        t = P._true_symbol()
        cmd = '%s %s %s'%(self._name, P._equality_symbol(), t)
        return P.eval(cmd) == t


    def __long__(self):
        return long(str(self))

    def __float____(self):
        return float(str(self))

    def _integer_(self):
        import sage.rings.all
        return sage.rings.all.Integer(str(self))

    def _rational_(self):
        import sage.rings.all
        return sage.rings.all.Rational(str(self))

    def name(self):
        return self._name

    def gen(self, n):
        P = self._check_valid()
        return P.new('%s.%s'%(self._name, int(n)))

    def _add_(self, right):
        P = self._check_valid()
        return P.new('%s + %s'%(self._name, right._name))

    def _sub_(self, right):
        P = self._check_valid()
        return P.new('%s - %s'%(self._name, right._name))

    def _mul_(self, right):
        P = self._check_valid()
        return P.new('%s * %s'%(self._name, right._name))

    def _div_(self, right):
        P = self._check_valid()
        return P.new('%s / %s'%(self._name, right._name))

    def __pow__(self, n):
        P = self._check_valid()
        if isinstance(n, ExpectElement):
            return P.new('%s ^ %s'%(self._name,n._name))
        else:
           return P.new('%s ^ %s'%(self._name,n))


def reduce_load(parent, x):
    return parent(x)

import os
def console(cmd):
    os.system(cmd)


