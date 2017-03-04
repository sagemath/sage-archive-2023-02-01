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

class PolymakeError(RuntimeError):
    """
    Raised if Polymake yields an error message
    """
    pass

def polymake_console():
    """
    Spawn a new Polymake command-line session.

    EXAMPLES::

    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%polymake magics instead.')
    os.system(os.getenv('SAGE_POLYMAKE_COMMAND') or 'polymake')

class Polymake(ExtraTabCompletion, Expect):
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

    #~ def _function_element_class(self):
        #~ return PolymakeFunctionElement

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
            return "%s(%s, {%s})"%(function, ",".join(list(args)), ",".join(list(kwds)))
        return "%s(%s)"%(function, ",".join(list(args)))

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
        return 'load_commands "{}";\n'.format(filename)

    def _keyboard_interrupt(self):
        if self._expect is None:
            raise KeyboardInterrupt
        print("Interrupting %s..." % self)
        cl = self._expect.compile_pattern_list(["killed by signal", self._prompt])
        while True:
            try:
                self._expect.send(chr(3))
            except pexpect.ExceptionPexpect as msg:
                raise pexpect.ExceptionPexpect("THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
            try:
                if self._expect.expect_list(cl, timeout=0.5): # 1, if self._prompt was found
                    break
            except pexpect.TIMEOUT:
                raise RuntimeError("{} interface is not responding".format(self))
        raise KeyboardInterrupt("Ctrl-c pressed while running %s"%self)

    def _synchronize(self, cmd='print 1+{};\n'):
        """
        TEST::

            sage: from sage.interfaces.polymake import polymake
            sage: Q = polymake.cube(4)
            sage: polymake('"ok"')
            ok
            sage: polymake._expect.sendline()
            1

        Now the interface is badly of sync::

            sage: polymake('"foobar"')
            <BLANKLINE>
            sage: Q.typeof()
            ('', 'Polymake::polytope::Polytope__Rational\nprint reftype($SAGE...);;')
            sage: Q.typeof.clear_cache()

        After synchronisation, things work again as expected::

            sage: polymake._synchronize()
            sage: polymake('"back to normal"')
            back to normal
            sage: Q.typeof()
            ('Polymake::polytope::Polytope__Rational', 'ARRAY')

        """
        if self._expect is None:
            return
        rnd = randrange(2147483647)
        s = str(rnd+1)
        self._sendstr(cmd.format(rnd))
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
        return r'$SAGE%s'%self.__seq

    def help(self, topic, pager=True):
        H = self.eval('help("{}");\n'.format(topic))
        if pager:
            from IPython.core.page import page
            page(H, start = 0)
        else:
            return H

    def eval(self, x, **kwds):
        s = Expect.eval(self, x, **kwds).lstrip().split(os.linesep)
        if s[0].startswith('polymake:  ERROR'):
            raise PolymakeError("\n".join(s))
        elif s[0].startswith('polymake: '):
            # presumably there is a comment to be stripped
            i = 0
            c = []
            while s:
                next_comment = s.pop(0)
                if next_comment:
                    c.append(next_comment)
                else:
                    break
            if get_verbose() > 0:
                for next_comment in c:
                    print(next_comment)
                print()
        return os.linesep.join(s)

    def get(self, var):
        return self.eval("print {};".format(var))

    def cputime(self, t=None):
        return NotImplemented

    @cached_method
    def _tab_completion(self):
        """
        Returns a list of Polymake function names.

        NOTE:

        It is not always the case that the returned functions are
        actually available in the current application.
        """
        s = self.eval("apropos '';").split(os.linesep)
        out = []
        for name in s:
            if name.startswith("/function"):
                out.append(name.split("/")[-1])
        return sorted(out)

    # Polymake specific methods

    def application(self, app):
        assert app in ["common", "fulton", "group", "matroid", "topaz", "fan", "graph", "ideal", "polytope", "tropical"], "Unknown Polymake application '{}'".format(app)
        if not self.is_running():
            self._start()
        self._application = app
        self._change_prompt("{} > ".format(app))
        self._sendstr('application "{}";{}'.format(app, os.linesep))
        self._expect.expect(self._prompt)

    def set_seed(self, seed=None):
        """
        Sets the seed for Polymake interface.
        The seed should be an integer.

        EXAMPLES::

        """
        if seed is None:
            seed = self.rand_seed()
        self.eval("srand(%d);" % seed)
        self._seed = seed
        return seed


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
        if 'Matrix' in T1 or 'Vector' in T1:
            return P.eval('print '+self._name+";").strip()
        if 'Polytope' in T1:
            return P.eval('print @{}[1];'.format(self._name)).strip()
        if 'RuleChain' in T1:
            return os.linesep.join(P.eval('print join("##",{}->list);'.format(self._name)).split('##'))
        if T1=='' and T2=='ARRAY':
            return P.eval('print @{};'.format(self._name)).strip()
        if T1=='' and T2=='HASH':
            return P.eval('print %@{};'.format(self._name)).strip()
        return P.eval('print '+self._name+";").strip()

    def __cmp__(self, other):
        """
        Comparison of Polymake elements.

        TESTS:


        """
        P = self._check_valid()
        if P.get("%s %s %s"%(self.name(), P._equality_symbol(), other.name())) == P._true_symbol():
            return 0
        if P.get("%s %s %s"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
            return -1
        if P.get("%s %s %s"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
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
        cmd = '%s %s %s'%(self._name, P._equality_symbol(), t)
        return P.get(cmd) == t

    @cached_method
    def _member_list(self):
        ### return the members of a "big" object.
        P = self._check_valid()
        try:
            X = P('typeof {};'.format(self._name))
        except (TypeError, PolymakeError):  # this happens for a perl type that isn't a Polymake type
            return []
        cmd = 'print join(", ", sorted_uniq(sort { $a cmp $b } map { keys %{$_->properties} }'+X._name + ', @{'+X._name+'->super}));'
        try:
            out = P.eval(cmd).split(', ')
        except PolymakeError, msg:
            return []
        return sorted(out)

    def known_properties(self):
        P = self._check_valid()
        try:
            return sorted(P.eval('print join(", ", {}->list_properties);'.format(self._name)).split(', '))
        except PolymakeError:
            return []

    def get_schedule(self, name):
        P = self._check_valid()
        return P('{}->get_schedule("{}")'.format(self._name, name))

    def _tab_completion(self):
        return sorted(self._member_list()+self.parent()._tab_completion())

    def __getattr__(self, attrname):
        P = self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        if attrname not in P._tab_completion():
            try:
                return P('{}->{}'.format(self._name, attrname))
            except (TypeError, PolymakeError):
                pass
        return P._function_element_class()(self, attrname)

    def attribute(self, attrname):
        P = self._check_valid()
        return P('%s->%s'%(self.name(), attrname))

    def __getitem__(self, key):
        P = self._check_valid()
        _,T = self.typeof()
        if T=='ARRAY':
            return P('@{}[{}]'.format(self._name, key))
        if T=='HASH':
            try:
                if key.parent() is self.parent():
                    key = key._name
                else:
                    key = str(key)
            except AttributeError:
                key = str(key)
            return P(self._name+"{"+key+"}")
        raise NotImplementedError("Cannot get items from Perl type {}".format(T))

    def __len__(self):
        P = self._check_valid()
        T1,T2 = self.typeof()
        if T2=='ARRAY':
            if self._name.startswith('@'):
                name = self._name[1:]
            else:
                name = self._name
            return int(P.eval('print scalar @{};'.format(name)))
        if T2=='HASH':
            if self._name.startswith('%'):
                name = self._name[1:]
            else:
                name = self._name
            return int(P.eval('print scalar keys %{};'.format(name)))
        if T1:
            raise TypeError("Don't know how to compute the length of {} object".format(T1))
        raise TypeError("Don't know how to compute the length of {} object".format(T2))

    @cached_method
    def typeof(self):
        P = self._check_valid()
        P.eval('use Scalar::Util qw(reftype);')
        return P.get('ref({});'.format(self._name)), P.get('reftype({});'.format(self._name))
