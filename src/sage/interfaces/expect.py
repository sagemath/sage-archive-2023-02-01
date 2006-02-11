"""
Common Interface Functionality

See the examples in the other sections for how to use specific
interfaces.  The interface classes all derive from the generic
interface that is described in this section.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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
import pexpect
import weakref
from sage.ext.sage_object import SageObject
from sage.structure.element import Element_cmp_
from sage.misc.sage_eval import sage_eval

import quit

import sage.rings.coerce as coerce
from sage.misc.misc import SAGE_ROOT, verbose, SAGE_TMP_INTERFACE
from sage.structure.element import RingElement
BAD_SESSION = -2

failed_to_start = []

tmp='%s/tmp'%SAGE_TMP_INTERFACE

class Expect(SageObject):
    def __init__(self, name, prompt, command=None, server=None, maxread=100000,
                 script_subdirectory="", restart_on_ctrlc=False,
                 verbose_start=False, init_code=[], max_startup_time=30,
                 logfile = None, eval_using_file_cutoff=0):

        self.__is_remote = False
        if command == None:
            command = name
        if server != None:
            command = "ssh -t %s %s"%(server, command)
            self.__is_remote = True
            eval_using_file_cutoff = 0  # don't allow this!
            #print command
            self._server = server
        self.__maxread = maxread
        self._eval_using_file_cutoff = eval_using_file_cutoff
        self.__script_subdirectory = script_subdirectory
        self.__name = name
        self.__coerce_name = '_' + name + '_'
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

    def is_remote(self):
        return self.__is_remote

    def is_local(self):
        return not self.__is_remote

    def user_dir(self):
        return self.__path

    def __repr__(self):
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

    def expect(self):
        if self._expect is None:
            self._start()
        return self._expect

    def interact(self):
        r"""
        This allows you to interactively interact with the child
        interpreter.  Press Ctrl-] to exit and return to SAGE.

        {\em This function does not work very well. }

        \note{This is completely different than the console() member
        function.  The console function opens a new copy of the child
        interepreter, whereas the interact function gives you
        interactive access to the interpreter that is being used by
        SAGE.}
        """
        if not hasattr(self, '__expect'):
            self._start()
        self._eval_line('')
        self._eval_line('')
        self._expect.interact()
        self._eval_line('')
        self._eval_line('')

    def _start(self, alt_message=None):
        global failed_to_start
        if self.__name in failed_to_start:
            if alt_message:
                raise RuntimeError, alt_message
            else:
                raise RuntimeError, 'Unable to start %s (%s failed to start during this SAGE session; not attempting to start again)'%(self.__name, self.__name)

        self._session_number += 1
        current_path = os.path.abspath('.')
        dir = self.__path
        if not os.path.exists(dir):
            os.makedirs(dir)
        os.chdir(dir)
        if self.__verbose_start:
            print "Starting %s"%self.__command.split()[0]
        try:
            self._expect = pexpect.spawn(self.__command, logfile=self.__logfile)
        except pexpect.ExceptionPexpect:
            self._expect = None
            self._session_number = BAD_SESSION
            failed_to_start.append(self.__name)
            raise RuntimeError, "Unable to start %s because the command '%s' failed."%(self.__name, self.__command)
        os.chdir(current_path)
        self._expect.timeout = self.__max_startup_time
        #self._expect.setmaxread(self.__maxread)
        self._expect.maxread = self.__maxread
        self._expect.delaybeforesend = 0
        try:
            self._expect.expect(self._prompt)
        except pexpect.TIMEOUT:
            self._expect = None
            self._session_number = BAD_SESSION
            failed_to_start.append(self.__name)
            raise RuntimeError, "Unable to start %s"%self.__name
        self._expect.timeout = 9999999  # don't make this bigger, or it breaks on OS X
        for X in self.__init_code:
            self.eval(X)

    def clear_prompts(self):
        while True:
            try:
                self._expect.expect(self._prompt, timeout=0.1)
            except pexpect.TIMEOUT:
                return

    def __del__(self):
        if self._expect is None:
            return
        try:
            self.quit()
        except:
            pass
        # The following programs around a bug in pexpect.
        def dummy(): pass
        try:
            self._expect.close = dummy
        except AttributeError:
            pass

    def quit(self):
        if self._expect is None:
            return
        try:
            self._expect.timeout = 1  # give it at most 1 second to quit
            #self._expect.sendline(self._quit_string()+';;')
            self._eval_line(self._quit_string(), allow_use_file=False)
            self._expect.timeout = 9999999
        except (OSError, pexpect.TIMEOUT, pexpect.EOF, RuntimeError): # intentional EOF or RuntimeError if already quit
            pass

    def _quit_string(self):
        return 'quit'

    def _read_in_file_command(self, filename):
        raise NotImplementedError

    def _eval_line_using_file(self, line, tmp):
        F = open(tmp, 'w')
        F.write(line)
        F.close()
        s = self._eval_line(self._read_in_file_command(tmp), allow_use_file=False)
        return s

    def _eval_line(self, line, allow_use_file=True):
        if line.find('\n') != -1:
            raise ValueError, "line must not contain any newlines"
        if allow_use_file and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line, tmp)
        try:
            if self._expect is None:
                self._start()
            E = self._expect
            try:
                E.sendline(line)
            except OSError:
                return RuntimeError, "Error evaluating %s in %s"%(line, self)
            if len(line)>0:
                try:
                    E.expect(self._prompt)
                except pexpect.EOF, msg:
                    if self._quit_string() in line:
                        # we expect to get an EOF if we're quitting.
                        return ''
                    print "** %s crashed or quit executing '%s' **"%(self, line)
                    print "Restarting %s and trying again"%self
                    self._expect.close(force=1)
                    self._start()
                    if line != '':
                        return self._eval_line(line, allow_use_file=allow_use_file)
                    else:
                        return ''
                    #raise RuntimeError, "%s crashed executing %s"%(self, line)
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

    def eval(self, code):
        try:
            return '\n'.join([self._eval_line(L) for L in str(code).split('\n')])
        except KeyboardInterrupt:
            self._keyboard_interrupt()
        except TypeError, s:
            return 'error evaluating "%s":\n%s'%(code,s)

    def _coerce_(self, x):
        return self(x)

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

        elif (not isinstance(x, ExpectElement) and hasattr(x, self.__coerce_name)) or  \
                (isinstance(x, ExpectElement) and x.hasattr(self.__coerce_name)):
            return getattr(x, self.__coerce_name)(self)

        elif isinstance(x, (list, tuple)):
            A = []
            for v in x:
                if isinstance(v, cls):
                    A.append(v.name())
                else:
                    A.append(str(v))
            X = ','.join(A)
            return self.new('%s%s%s'%(self._left_list_delim(), X, self._right_list_delim()))

        return cls(self, x)

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
        return "=="

    # These below could easily be wrong and it not be obvious, so we require
    # they be overloaded.  The above would obviously bomb if not correct.
    def _true_symbol(self):
        raise NotImplementedError

    def _false_symbol(self):
        raise NotImplementedError

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
        self.set(name, value)
        return name

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
        if self is other:
            return 0
        return -1

    def console(self):
        raise NotImplementedError


class ExpectFunction(SageObject):
    def __init__(self, parent, name):
        self._parent = parent
        self._name = name

    def __repr__(self):
        return "Function %s"%self._name

    def __call__(self, *args):
        return self._parent.function_call(self._name, list(args))


class FunctionElement(SageObject):
    def __init__(self, obj, name):
        self._obj = obj
        self._name = name

    def __repr__(self):
        return "member function %s"%self._name

    def __call__(self, *args):
        return self._obj.parent().function_call(self._name, [self._obj] + list(args))

def is_ExpectElement(x):
    return isinstance(x, ExpectElement)

class ExpectElement(Element_cmp_, RingElement):
    def __init__(self, parent, value, is_name=False):
        RingElement.__init__(self, parent)
        self._create = value
        if parent is None: return     # means "invalid element"
        if isinstance(value, str) and parent._eval_using_file_cutoff and \
           parent._eval_using_file_cutoff > len(value):
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

    def _cmp_(self, other):
        #if not (isinstance(other, ExpectElement) and other.parent() is self.parent()):
        #    return coerce.cmp(self, other)
        P = self.parent()
        if P.eval("%s < %s"%(self.name(), other.name())) == P._true_symbol():
            return -1
        elif P.eval("%s > %s"%(self.name(), other.name())) == P._true_symbol():
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

    def __contains__(self, x):
        self._check_valid()
        P = self.parent()
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
        return sage_eval(str(self))

    def __repr__(self):
        self._check_valid()
        try:
            if self._get_using_file:
                return self.parent().get_using_file(self._name)
        except AttributeError:
            return self.parent().get(self._name)

    def __getattr__(self, attrname):
        self._check_valid()
        return FunctionElement(self, attrname)

    def hasattr(self, attrname):
        """
        Returns whether the given attribute is already defined by this object,
        and in particular is not dynamically generated.
        """
        return not isinstance(getattr(self, attrname), FunctionElement)


    def __getitem__(self, n):
        self._check_valid()
        if not isinstance(n, tuple):
            return self.parent().new('%s[%s]'%(self._name, n))
        else:
            return self.parent().new('%s[%s]'%(self._name, str(n)[1:-1]))

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
        self._check_valid()
        return self.parent().new('%s.%s'%(self._name, int(n)))

    def _add_(self, right):
        self._check_valid()
        return self.parent().new('%s + %s'%(self._name, right._name))

    def _sub_(self, right):
        self._check_valid()
        return self.parent().new('%s - %s'%(self._name, right._name))

    def _mul_(self, right):
        self._check_valid()
        return self.parent().new('%s * %s'%(self._name, right._name))

    def _div_(self, right):
        self._check_valid()
        return self.parent().new('%s / %s'%(self._name, right._name))

    def __pow__(self, n):
        self._check_valid()
        if isinstance(n, ExpectElement):
            return self.parent().new('%s ^ %s'%(self._name,n._name))
        else:
            return self.parent().new('%s ^ %s'%(self._name,n))


def reduce_load(parent, x):
    return parent(x)

import os
def console(cmd):
    os.system(cmd)


