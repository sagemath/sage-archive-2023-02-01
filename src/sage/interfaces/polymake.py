r"""
Interface to Polymake

"""


#*****************************************************************************
#       Copyright (C) 2017 Simon King <simon.king@uni-jena.de>
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
from __future__ import print_function
from __future__ import absolute_import

import os
import re
import sys
import six

from sage.structure.parent import Parent
from .expect import console, Expect, ExpectElement, ExpectFunction, FunctionElement

from sage.env import SAGE_EXTCODE, DOT_SAGE
from sage.misc.misc import get_verbose
from sage.misc.cachefunc import cached_method
from sage.interfaces.tab_completion import ExtraTabCompletion

import pexpect
from random import randrange

from time import sleep
from six.moves import range
from six import reraise as raise_
import warnings

_name_pattern = re.compile('SAGE[0-9]+')

_available_polymake_answers = {
    0: "returns prompt",
    1: "returns continuation prompt",
    2: "requests interactive input",
    3: "kills computation",
    4: "raises error",
    5: "issues warning",
    6: "shows additional information",
    7: "lost connection",
    8: "fails to respond timely"
        }

class PolymakeError(RuntimeError):
    """
    Raised if Polymake yields an error message.

    TESTS::

        sage: polymake.eval('print foo;')    # optional polymake
        Traceback (most recent call last):
        ...
        PolymakeError: Unquoted string "foo" may clash with future reserved word at input line 1.

    """
    pass

def polymake_console():
    """
    Spawn a new Polymake command-line session.

    EXAMPLES::

        sage: from sage.interfaces.polymake import polymake_console
        sage: polymake_console()        # not tested
        Welcome to polymake version ...
        ...
        Ewgenij Gawrilow, Michael Joswig (TU Berlin)
        http://www.polymake.org

        This is free software licensed under GPL; see the source for copying conditions.
        There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

        Press F1 or enter 'help;' for basic instructions.

        Application polytope currently uses following third-party software packages:
        4ti2, bliss, cdd, latte, libnormaliz, lrs, permlib, ppl, sketch, sympol, threejs, tikz, topcom, tosimplex
        For more details:  show_credits;
        polytope >

    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%polymake magics instead.')
    os.system(os.getenv('SAGE_POLYMAKE_COMMAND') or 'polymake')

class Polymake(ExtraTabCompletion, Expect):
    r"""
    Interface to the Polymake interpreter.

    In order to use this interface, you need to either install the
    optional Polymake package for Sage, or install Polymake system-wide
    on your computer.

    Type ``polymake.[tab]`` for a list of most functions
    available from your Polymake install. Type
    ``polymake.Function?`` for Polymake's help about a given ``Function``
    Type ``polymake(...)`` to create a new Magma
    object, and ``polymake.eval(...)`` to run a string using
    Polymake and get the result back as a string.

    EXAMPLES::

    """
    def __init__(self, script_subdirectory=None,
                 logfile=None, server=None,server_tmpdir=None,
                 seed=None, command=None):
        if command is None:
            command = "env TERM=dumb {}".format(os.getenv('SAGE_POLYMAKE_COMMAND') or 'polymake')
        prompt = 'polytope > '
        Expect.__init__(self,
                        name="polymake",
                        prompt=prompt,
                        command=command,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=False,
                        logfile=logfile,
                        eval_using_file_cutoff=100)

        self._seed = seed
        self.__tab_completion = {}

    @cached_method
    def version(self):
        import subprocess
        return subprocess.check_output(["polymake", "--version"], stderr=subprocess.STDOUT).split()[2]
    # Pickling etc

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.polymake import polymake
            sage: polymake.__reduce__()
            (<function reduce_load_Polymake at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            Polymake
        """
        return reduce_load_Polymake, tuple([])

    def _object_class(self):
        return PolymakeElement

    def _function_element_class(self):
        return PolymakeFunctionElement

    #~ def _function_class(self):
        #~ return PolymakeFunction

    def function_call(self, function, args=None, kwds=None):
        """
        EXAMPLES::

            sage: maxima.quad_qags(x, x, 0, 1, epsrel=1e-4)
            [0.5,0.55511151231257...e-14,21,0]
            sage: maxima.function_call('quad_qags', [x, x, 0, 1], {'epsrel':'1e-4'})
            [0.5,0.55511151231257...e-14,21,0]
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        s = self._function_call_string(function,
                                       [s.name() for s in args],
                                       ['%s=>%s'%(key,value.name()) for key, value in kwds.items()])
        return self(s)

    def _function_call_string(self, function, args, kwds):
        """
        Returns the string used to make function calls.

        EXAMPLES::

            sage: maxima._function_call_string('diff', ['f(x)', 'x'], [])
            'diff(f(x),x)'
        """
        if kwds:
            if args:
                return "%s(%s, %s);"%(function, ",".join(list(args)), ",".join(list(kwds)))
            return "%s(%s);"%(function, ",".join(list(kwds)))
        return "%s(%s);"%(function, ",".join(list(args)))

    def console(self):
        """
        Spawn a new Polymake command-line session.

        EXAMPLES::

            sage: polymake.console()  # not tested
        """
        polymake_console()

    # Methods concerning interface communication

    def _install_hints(self):
        return "Please install the optional Polymake package for sage or install Polymake system-wide"

    def _start(self, alt_message=None):
        Expect._start(self, alt_message)
        self.application("polytope")
        self.eval('use Scalar::Util qw(reftype);')
        self.eval('use Scalar::Util qw(blessed);')

    def _quit_string(self):
        return "exit;"

    def _assign_symbol(self):
        return "="

    def _equality_symbol(self):
        return "=="

    #~ def _true_symbol(self):
    #~ def _false_symbol(self):
    #~ def _lessthan_symbol(self):
    #~ def _greaterthan_symbol(self):
    #~ def _inequality_symbol(self):
    #~ def _exponent_symbol(self):  # scientific notation for floats, the default "e" is fine


    #~ def _continuation_prompt(self):
        #~ # Problem: it is of the form "application (number)> ", where number is the number of the current line
        #~ return NotImplemented

    def _read_in_file_command(self, filename):
        return 'script "{}";\n'.format(filename)

    def _keyboard_interrupt(self):
        if not self.is_running():
            raise KeyboardInterrupt
        print("Interrupting %s..." % self)
        while True:
            try:
                self._expect.send(chr(3))
            except pexpect.ExceptionPexpect as msg:
                raise pexpect.ExceptionPexpect("THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
            i = self._expect.expect_list(self._prompt, timeout=1)
            if i==0:
                break
            elif i==7:  # EOF
                warnings.warn("Polymake {} during keyboard interrupt".format(_available_polymake_answers[i]), RuntimeWarning)
                self._crash_msg()
                self.quit()
            elif i==8:  # Timeout
                self.quit()
                raise RuntimeError("{} interface is not responding. We closed it".format(self))
            else:
                warnings.warn("Polymake {} during keyboard interrupt".format(_available_polymake_answers[i]), RuntimeWarning)
        raise KeyboardInterrupt("Ctrl-c pressed while running %s"%self)

    def _synchronize(self):
        """
        TEST::

            sage: from sage.interfaces.polymake import polymake
            sage: Q = polymake.cube(4)
            sage: polymake('"ok"')
            ok
            sage: polymake._expect.sendline()
            1

        Now the interface is badly out of sync::

            sage: polymake('"foobar"')
            <BLANKLINE>
            sage: Q.typeof()
            ('', 'Polymake::polytope::Polytope__Rational', '$SAGE...')
            sage: Q.typeof.clear_cache()

        After synchronisation, things work again as expected::

            sage: polymake._synchronize()
            sage: polymake('"back to normal"')
            back to normal
            sage: Q.typeof()
            ('Polymake::polytope::Polytope__Rational', 'ARRAY', '$SAGE...')

        """
        if not self.is_running():
            return
        rnd = randrange(2147483647)
        res = str(rnd+1)
        cmd='print 1+{};'+self._expect.linesep
        self._sendstr(cmd.format(rnd))
        pat = self._expect.expect(self._prompt,timeout=0.5)
        # 0: normal prompt
        # 1: continuation prompt
        # 2: user input expected when requestion "help"
        # 3: what we are looking for when interrupting a computation
        # 4: error
        # 5: warning
        # 6: anything but an error or warning, thus, an information
        # 7: unexpected end of the stream
        # 8: (expected) timeout
        if pat == 8: # timeout
            warnings.warn("{} unexpectedly {} during synchronisation.".format(self, _available_polymake_answers[pat]), RuntimeWarning)
            self.interrupt()
            # ... but we continue, as that probably means we currently are at the end of the buffer
        elif pat == 7: # EOF
            self._crash_msg()
            self.quit()
        elif pat == 0:
            # We got the right prompt, but perhaps in a wrong position in the stream
            # The result of the addition should appear *before* our prompt
            if not res in self._expect.before:
                try:
                    warnings.warn("{} seems out of sync: The expected output did not appear before reaching the next prompt.".format(self))
                    while True:
                        i = self._expect.expect_list(self._prompt, timeout=0.1)
                        if i==8: # This time, we do expect a timeout
                            return
                        elif i>0:
                            raise RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[i]))
                except pexpect.TIMEOUT:
                    warnings.warn("A timeout has occured when synchronising {}.".format(self), RuntimeWarning)
                    self._interrupt()
                except pexpect.EOF:
                    self._crash_msg()
                    self.quit()
            else:
                return
        else:
            raise RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[pat]))

    def _next_var_name(self):
        r"""
        Returns the next unused variable name.

        This is for "scalar" objects only!

        """
        if len(self._available_vars) != 0:
            return self._available_vars.pop(0)
        try:
            self.__seq += 1
        except AttributeError:
            self.__seq = 0
        return r'SAGE%s'%self.__seq

    def clear(self, var):
        """
        Clear the variable named var.
        """
        self._available_vars.append(_name_pattern.search(var).group())

    def _create(self, value, name=None):
        name = self._next_var_name() if name is None else name
        self.set(name, value)
        # If value is a list, then @name is now equal to that list.
        # Otherwise, value is obtained by $name[0]. So, we modify
        # the name returned by _create so that it can be used to
        # access the wrapped value.
        if self.eval('print scalar @{};'.format(name)).strip() == '1':
            return '$'+name+'[0]'
        return '@'+name

    def set(self, var, value):
        """
        Set the variable var to the given value.

        Eventually, ``var`` is a reference to ``value``.
        """
        if isinstance(value, six.string_types):
            value = value.strip().rstrip(';').strip()
        cmd = '@%s%s(%s);'%(var,self._assign_symbol(), value)
        self.eval(cmd)

    def get(self, cmd):
        return self.eval("print {};".format(cmd)).strip()

    def help(self, topic, pager=True):
        H = self.eval('help("{}");\n'.format(topic))
        if pager:
            from IPython.core.page import page
            page(H, start = 0)
        else:
            return H

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=True):
        line = line.strip()
        if allow_use_file and wait_for_prompt and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line)
        try:
            if not self.is_running():
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
                raise_(RuntimeError, "%s\nError evaluating %s in %s"%(msg, line, self), sys.exc_info()[2])

            p_warnings = []
            p_errors = []
            have_warning = False
            have_error = False
            have_log = False
            if len(line)>0:
                first = True
                while True:
                    try:
                        if isinstance(wait_for_prompt, six.string_types):
                            pat = E.expect(wait_for_prompt)
                        else:
                            pat = E.expect_list(self._prompt)
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
                    if self._terminal_echo and first:
                        i = out.find("\n")
                        j = out.rfind("\r")
                        out = out[i+1:j].replace('\r\n','\n')
                    else:
                        out = out.strip().replace('\r\n','\n')
                    first = False
                    if have_error:
                        p_errors.append(out)
                        have_error = False
                        out = ""
                    elif have_warning:
                        p_warnings.append(out)
                        have_warning = False
                        out = ""
                    elif have_log:
                        if get_verbose() > 0:
                            print(out)
                        have_log = False
                        out = ""
                    # 0: normal prompt
                    # 1: continuation prompt
                    # 2: user input expected when requestion "help"
                    # 3: what we are looking for when interrupting a computation
                    # 4: error
                    # 5: warning
                    # 6: anything but an error or warning, thus, an information
                    # 7: unexpected end of the stream
                    # 8: (expected) timeout
                    if pat == 0:
                        have_log = False
                        have_error = False
                        have_warning = False
                        if E.buffer:
                            if not E.buffer.strip():
                                E.send(chr(3))
                                pat = E.expect_list(self._prompt)
                                if E.buffer or pat:
                                    raise RuntimeError("Couldn't return to prompt after command '{}'".format(line))
                        break
                    elif pat == 1: # unexpected continuation prompt
                        # Return to normal prompt
                        i = pat
                        E.send(chr(3))
                        i = E.expect_list(self._prompt)
                        assert i==0, "Command '{}': Couldn't return to normal prompt after Polymake {}. Instead, Polymake {}".format(line,_available_polymake_answers[pat],_available_polymake_answers[i])
                        raise SyntaxError("Incomplete Polymake command '{}'".format(line))
                    elif pat == 2: # request for user interaction
                        # Return to normal prompt
                        warnings.warn("{} expects user interaction. We abort and return the options that {} provides.".format(self,self))
                        i = pat
                        while i:
                            self._expect.send(chr(3))
                            i = self._expect.expect(self._prompt, timeout=0.1)
                        # User interaction is expected to happen when requesting help
                        if line.startswith('help'):
                            out = os.linesep.join(out.split(os.linesep)[:-1])
                            break
                        else:
                            RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[pat]))
                    elif pat == 3: # killed by signal
                        i = pat
                        while pat != 0:
                            E.send(chr(3))
                            i = E.expect_list(self._prompt)
                        RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[pat]))
                    elif pat == 4: # Polymake error
                        have_error = True
                    elif pat == 5: # Polymake warning
                        have_warning = True
                    elif pat == 6: # apparently Polymake prints a comment
                        have_log = True
                    elif pat == 7: # we have reached the end of the buffer
                        warnings.warn("Polymake unexpectedly {}".format(_available_polymake_answers[pat]), RuntimeWarning)
                        E.buffer = E.before + E.after + E.buffer
                        break
                    else: # timeout or some other problem
                        raise RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[pat]))
            else:
                out = ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt("Ctrl-c pressed while running %s"%self)
        for w in p_warnings:
            warnings.warn(w, RuntimeWarning)
        for e in p_errors:
            raise PolymakeError(e)
        return out

    def cputime(self, t=None):
        return NotImplemented

    def _tab_completion(self):
        """
        Returns a list of Polymake function names.

        NOTE:

        - The list of functions depends on the current application. The
          result is cached, of course separately for each application.
        - It is generally not the case that all the returned function names
          can actually successfully be called.

        """
        if not self.is_running():
            self._start()
        try:
            return self.__tab_completion[self._application]
        except KeyError:
            pass
        s = self.eval("apropos '';").split(self._expect.linesep)
        out = []
        for name in s:
            if name.startswith("/function"):
                out.append(name.split("/")[-1])
        self.__tab_completion[self._application] = sorted(out)
        return self.__tab_completion[self._application]

    # Polymake specific methods

    def application(self, app):
        assert app in ["common", "fulton", "group", "matroid", "topaz", "fan", "graph", "ideal", "polytope", "tropical"], "Unknown Polymake application '{}'".format(app)
        if not self.is_running():
            self._start()
        self._application = app
        patterns = ["{} > ".format(app),            # 0: normal prompt
                    "{} \([0-9]+\)> ".format(app),  # 1: continuation prompt
                    "Please choose ".format(app),   # 2: user input expected when requesting "help"
                    "killed by signal",             # 3: what we are looking for when interrupting a computation
                    "polymake: +ERROR: +",          # 4: error
                    "polymake: +WARNING: +",        # 5: warning
                    "polymake: +",                  # 6: anything but an error or warning, thus, an information
                    pexpect.EOF,                    # 7: unexpected end of the stream
                    pexpect.TIMEOUT]                # 8: timeout

        self._change_prompt(self._expect.compile_pattern_list(patterns))
        self._sendstr('application "{}";{}'.format(app, self._expect.linesep))
        pat = self._expect.expect_list(self._prompt)
        if pat:
            raise RuntimeError("When changing the application, Polymake unexpectedly {}".format(_available_polymake_answers[pat]))

    def new_object(self, name, *args, **kwds):
        try:
            f = self.__new[name]
        except AttributeError:
            self.__new = {}
            f = self.__new[name] = self._function_class()(self, "new %s"%name)
        except KeyError:
            f = self.__new[name] = self._function_class()(self, "new %s"%name)
        return f(*args, **kwds)

polymake = Polymake()

def reduce_load_Polymake():
    """
    Returns the Polymake interface object defined in sage.interfaces.polymake.

    EXAMPLES::

        sage: from sage.interfaces.polymake import reduce_load_Polymake
        sage: reduce_load_Polymake()
        Polymake
    """
    return polymake

########################################
## Elements

from warnings import warn

class PolymakeElement(ExtraTabCompletion, ExpectElement):
    def _repr_(self):
        T1, T2 = self.typeof()
        P = self._check_valid()
        name = self._name
        if T1:
            Temp = self.typename()
            if Temp:
                T1 = Temp
        if T1 in ['Matrix', 'Vector']:
            out = P.get(name).strip()
        elif 'RuleChain' in T1:
            out = os.linesep.join(P.get('join("##",{}->list)'.format(name)).split('##'))
        else:
            try:
                out = P.get('{}->description'.format(name)).strip()
            except PolymakeError:
                out = ''
            if os.linesep in out:
                out = ''
        if not out:
            if "Polytope" == T1:
                out = "{}[{}]".format(P.get("{}->type->full_name".format(name)) or "PolymakeElement", _name_pattern.search(name).group())
            elif T1=='' and T2=='ARRAY':
                out = P.get('@{}'.format(name)).strip()
            elif T1=='' and T2=='HASH':
                out = P.get('%{}'.format(name)).strip()
            else:
                out = P.get(name).strip()
        return out

    def __cmp__(self, other):
        """
        Comparison of Polymake elements.

        TESTS:


        """
        P = self._check_valid()
        if P.eval("print $%s %s $%s;"%(self.name(), P._equality_symbol(), other.name())).strip() == P._true_symbol():
            return 0
        if P.eval("print $%s %s $%s;"%(self.name(), P._lessthan_symbol(), other.name())).strip() == P._true_symbol():
            return -1
        if P.eval("print $%s %s $%s;"%(self.name(), P._greaterthan_symbol(), other.name())).strip() == P._true_symbol():
            return 1
        return -2 # that's supposed to be an error value.

    def bool(self):
        """
        Return whether this Polymake element is equal to ``True``.

        EXAMPLES::

            sage: from sage.interfaces.polymake import polymake
            sage: polymake(0).bool()    # optional polymake
            False
            sage: polymake(1).bool()    # optional polymake
            True

        """
        P = self._check_valid()
        t = P._true_symbol()
        cmd = '%s %s %s;'%(self._name, P._equality_symbol(), t)
        return P.get(cmd) == t

    def known_properties(self):
        P = self._check_valid()
        try:
            return sorted(P.get('join(", ", {}->list_properties)'.format(self._name)).split(', '))
        except PolymakeError:
            return []

    @cached_method
    def _member_list(self):
        ### return the members of a "big" object.
        P = self._check_valid()
        try:
            P.eval('$SAGETMP = typeof {+'+self._name+'};')
        except (TypeError, PolymakeError):  # this happens for a perl type that isn't a Polymake type
            return []
        cmd = 'print join(", ", sorted_uniq(sort { $a cmp $b } map { keys %{$_->properties} }$SAGETMP, @{$SAGETMP->super}));'
        try:
            out = P.eval(cmd).split(', ')
        except PolymakeError, msg:
            return []
        return sorted(out)

    def typename(self):
        P = self._check_valid()
        try:
            return P.eval('print {}->type->name;'.format(self._name))
        except PolymakeError as msg:
            return ''

    def full_typename(self):
        P = self._check_valid()
        try:
            return P.eval('print {}->type->full_name;'.format(self._name))
        except PolymakeError as msg:
            return ''

    def qualified_typename(self):
        P = self._check_valid()
        try:
            return P.eval('print {}->type->qualified_name;'.format(self._name))
        except PolymakeError as msg:
            return ''

    def _tab_completion(self):
        return sorted(self._member_list()+self.parent()._tab_completion())

    def __getattr__(self, attrname):
        P = self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        if attrname not in P._tab_completion():
            if attrname in self._member_list():
                try:
                    return P('{}->{}'.format(self._name, attrname))
                except (TypeError, PolymakeError):
                    raise AttributeError
            else:
                try:
                    return P._function_element_class()(self, '{}->{}'.format(self._name, attrname), memberfunction=True)
                except (TypeError, PolymakeError):
                    raise AttributeError
        return P._function_element_class()(self, attrname, memberfunction=False)

    def get_member_function(self, attrname):
        P = self._check_valid()
        return P._function_element_class()(self, '{}->{}'.format(self._name, attrname), memberfunction=True)

    def get_member(self, attrname):
        P = self._check_valid()
        return P('%s->%s'%(self.name(), attrname))

    def __getitem__(self, key):
        P = self._check_valid()
        _, T = self.typeof()
        if self._name.startswith('@'):
            return P('${}[{}]'.format(self._name[1:], key))
        if T=='ARRAY':
            return P('{}[{}]'.format(self._name, key))
        if T=='HASH':
            try:
                if key.parent() is self.parent():
                    key = key._name
                else:
                    key = str(key)
            except AttributeError:
                key = str(key)
            return P(name+"{"+key+"}")
        raise NotImplementedError("Cannot get items from Perl type {}".format(T))

    def __len__(self):
        P = self._check_valid()
        T1,T2 = self.typeof()
        name = self._name
        if T2=='ARRAY':
            return int(P.eval('print scalar @{+%s};'%name))
        if T2=='HASH':
            return int(P.eval('print scalar keys %{+%s};'%name))
        if T1:
            raise TypeError("Don't know how to compute the length of {} object".format(T1))
        raise TypeError("Don't know how to compute the length of {} object".format(T2))

    @cached_method
    def typeof(self):
        P = self._check_valid()
        name = self._name
        return P.eval('print ref({});'.format(name)), P.eval('print reftype({});'.format(name))

    def _sage_doc_(self):
        """
        EXAMPLES::

        """
        P = self._check_valid()
        # according to Julian Pfeifle, the only case in which the fully qualified
        # typename would not provide the doc.
        Tname = self.typename()
        Tqname = ""
        Tfname = ""
        if Tname == 'Polytope':
            doc = P.eval('help "Polytope";')
        else:
            Tqname = self.qualified_typename()
            try:
                doc = P.eval('help "{}";'.format(Tqname))
            except PolymakeError:
                doc = ''
        #~ T1, T2 = self.typeof()
        #~ if not T1:
            #~ if T2:
                #~ try:
                    #~ P.eval('help "{}";'.format(T2))
                #~ except PolymakeError:
                    #~ return "{} element of unknown type".format(P)
        #~ try:
            #~ return P.eval('help "{}";'.format(T1))
        #~ except PolymakeError:
            #~ pass
        #~ names = T1.split('::')
        #~ # Search backwards (from more specific to more general topics)
        #~ # in what is separated by "::", and in each topic search forward
        #~ # in the names separated by "_".
        #~ for name in reversed(names):
            #~ try:
                #~ return P.eval('help "{}";'.format(name))
            #~ except PolymakeError:
                #~ pass
            #~ if "__" in name:
                #~ for part in name.split("_"):
                    #~ if not part:
                        #~ continue
                    #~ try:
                        #~ return P.eval('help "{}";'.format(part))
                    #~ except PolymakeError:
                        #~ pass
        if doc:
            return doc
        return "Undocumented polymake type '{}'".format(self.full_typename())

class PolymakeFunctionElement(FunctionElement):
    """
    A callable member of a Polymake element.
    """
    def __init__(self, obj, name, memberfunction=False):
        self._obj = obj
        self._name = name
        self._is_memberfunc = memberfunction
    def _repr_(self):
        if self._is_memberfunc:
            return "Member function '{}' of {} object".format(self._name.split("->")[-1], self._obj.typeof()[0])
        return "{} (bound to {} object)".format(self._name, self._obj.typeof()[0])

    def __call__(self, *args, **kwds):
        if self._is_memberfunc:
            return self._obj._check_valid().function_call(self._name, list(args), kwds)
        return self._obj._check_valid().function_call(self._name, [self._obj] + list(args), kwds)

    def _sage_doc_(self):
        P = self._obj._check_valid()
        return P.help(self._name.split("->")[-1])

