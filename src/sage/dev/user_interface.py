"""
User Interface

Provides small class that is used for displaying messages as well
as prompting for user input.
"""
from __future__ import print_function

import os

from getpass import getpass

try:
    import struct
    import fcntl
    import termios
    def _ioctl_GWINSZ(fd):
        """
        gets window size of the terminal at the
        file descriptor ``fd``

        TESTS::

            sage: import os
            sage: from sage.dev.user_interface import _ioctl_GWINSZ
            sage: try:                                    # not tested
            ....:     fd = os.open(os.ctermid(), os.O_RDONLY)
            ....:     _ioctl_GWINSZ(fd)
            ....: finally:
            ....:     os.close(fd)
            (48, 194)
        """
        try:
            return struct.unpack('hh',
                    fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
        except IOError:
            return

except ImportError:
    _ioctl_GWINSZ = lambda fd: None

class UserInterface(object):
    def select(self, prompt, options, default=None):
        """
        should ask user to select from list of options and return
        selected option

        INPUT:

        - ``prompt`` -- prompt to display

        - ``options`` -- iterable that contains list of options

        - ``default`` -- default option

        TESTS::

            sage: UI = sage.dev.user_interface.UserInterface()
            sage: s = UI.select("Should I delete your home directory?",
            ....:         ("yes","no","maybe"), default=1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def confirm(self, question, default_no=False):
        """
        should ask yes/no question from the developer and return the
        boolean of the response

        INPUT:

        - ``question`` -- a string

        - ``default_no`` -- boolean whether to default to no

        TESTS::

            sage: UI = sage.dev.user_interface.UserInterface()
            sage: UI.confirm("Should I delete your home directory?")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.select(question, ("yes", "no"), int(default_no)) == "yes"

    def get_input(self, prompt):
        """
        gets input from developer after displaying prompt

        TESTS::

            sage: UI = sage.dev.user_interface.UserInterface()
            sage: s = UI.get_input("What do you want for dinner?")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def get_password(self, prompt):
        """
        gets password from developer after displaying prompt

        TESTS::

            sage: UI = sage.dev.user_interface.UserInterface()
            sage: p = UI.get_password("What is the passphrase for your safe?")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def show(self, message):
        """
        should display message to user

        TESTS::

            sage: UI = sage.dev.user_interface.UserInterface()
            sage: UI.show("I ate fillet mignon for dinner.")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def edit(self, filename):
        """
        should drop user into editor with filename open

        TESTS::

            sage: UI = sage.dev.user_interface.UserInterface()
            sage: UI.edit("filename")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class CmdLineInterface(UserInterface):
    def _std_values(self, prompt, options, default):
        """
        returns standard version for prompt, options, and default

        TESTS::

            sage: UI = sage.dev.user_interface.CmdLineInterface()
            sage: UI._std_values("Should I delete your home directory?",
            ....:         ("yes","no","maybe"), default=1)
            ('Should I delete your home directory? [yes/No/maybe] ', ('yes', 'no', 'maybe'), 'no')
        """
        if options is not None:
            options = tuple(opt.lower() for opt in options)
            if options:
                prompt += " ["
                prompt += "/".join(opt if i != default else opt.capitalize()
                                    for i, opt in enumerate(options))
                prompt += "]"
                if default is not None:
                    default = options[default]
            else:
                options = None
        prompt += " "
        return prompt, options, default

    def _get_input(self,
            prompt,
            options=None,
            default=None,
            input_func=raw_input):
        """
        general implementation for :meth:`switch`, :meth:`get_input`,
        and :meth:`get_password`

        INPUT:

        - ``prompt`` -- a string

        - ``options`` -- a list of strings or None

        - ``default`` -- a string or None

        - ``input_func`` -- function to get input from user

        TESTS::

            sage: UI = sage.dev.user_interface.CmdLineInterface()
            sage: UI._get_input("Should I delete your home directory?", # not tested
            ....:         ("yes","no","maybe"), default=0)
            Should I delete your home directory? [Yes/no/maybe]
        """
        prompt, options, default = self._std_values(prompt, options, default)

        while True:
            s = input_func(prompt)

            if options is None:
                return s

            if len(s.strip()) == 0:
                if default is None:
                    self.show("Please enter an option.")
                    continue
                else:
                    return default

            itr = (opt for opt in options if opt.startswith(s))
            try:
                ret = itr.next()
            except StopIteration:
                self.show("Please specify an allowable option.")
                continue

            try:
                ret = itr.next()
                self.show("Please disambiguate between options.")
            except StopIteration:
                return ret

    def select(self, prompt, options, default=None):
        """
        asks user to select from list of options and returns selected
        option

        INPUT:

        - ``prompt`` -- prompt to display

        - ``options`` -- iterable that contains list of options

        - ``default`` -- default option

        TESTS::

            sage: UI = sage.dev.user_interface.DoctestInterface()
            sage: UI.append("n")
            sage: s = UI.select("Should I delete your home directory?",
            ....:         ("yes","no","maybe"), default=2)
            Should I delete your home directory? [yes/no/Maybe] n
            sage: print s
            no
        """
        return self._get_input(prompt, options, default)

    def get_input(self, prompt):
        """
        get input from developer after displaying prompt

        TESTS::

            sage: UI = sage.dev.user_interface.DoctestInterface()
            sage: UI.append("fillet mignon")
            sage: s = UI.get_input("What do you want for dinner?")
            What do you want for dinner? fillet mignon
            sage: print s
            fillet mignon
        """
        return self._get_input(prompt)

    def get_password(self, prompt):
        """
        gets password from developer after displaying prompt

        TESTS::

            sage: UI = sage.dev.user_interface.DoctestInterface()
            sage: UI.append("a number")
            sage: p = UI.get_password("What is the key combo for your safe?")
            What is the key combo for your safe?
            sage: print p
            a number
        """
        return self._get_input(prompt, input_func=getpass)

    def _get_dimensions(self):
        """
        trys to return the dimensions of the terminal

        TESTS::

            sage: UI = sage.dev.user_interface.CmdLineInterface()
            sage: UI._get_dimensions()        # not tested
            (48, 194)
        """
        dim = _ioctl_GWINSZ(0) or _ioctl_GWINSZ(1) or _ioctl_GWINSZ(2)

        if dim is None:
            try:
                fd = os.open(os.ctermid(), os.O_RDONLY)
                dim = _ioctl_GWINSZ(fd)
            finally:
                os.close(fd)

        if dim is None:
            raise EnvironmentError("cannot determine dimensions of terminal")

        return tuple(int(x) for x in dim)

    def show(self, message):
        """
        displays message to user

        TESTS::

            sage: UI = sage.dev.user_interface.CmdLineInterface()
            sage: UI.show("I ate fillet mignon for dinner.")
            I ate fillet mignon for dinner.
        """
        try:
            height, width = self._get_dimensions()
        except EnvironmentError:
            height, width = float('inf'), float('inf')

        message = message.strip().splitlines()
        message = [line.rstrip() for line in message]
        if (len(message)+2 <= height and
                max(len(line) for line in message) <= width):
            print(*message, sep='\n')
        else:
            message = '\n'.join(message)+'\n'
            try:
                self._pager(message)
            except AttributeError:
                import pydoc
                self._pager = pydoc.getpager()
                self._pager(message)

    def edit(self, filename):
        """
        drops user into editor with filename open

        TESTS::

            sage: import os, tempfile
            sage: tmp = tempfile.mkstemp()[1]
            sage: UI = sage.dev.user_interface.CmdLineInterface()
            sage: os.environ['EDITOR'] = 'echo "stuff to put in file" >'
            sage: UI.edit(tmp)
            sage: print open(tmp,'r').read()
            stuff to put in file
            sage: os.unlink(tmp)
        """
        r = os.system(" ".join((os.environ.get('EDITOR', 'nano'), filename)))
        if r: return r

class DoctestInterface(CmdLineInterface, list):
    def _get_input(self, *args, **kwds):
        """
        overwrites `input_func` with custom function for doctesting

        TESTS::

            sage: UI = sage.dev.user_interface.DoctestInterface()
            sage: r = UI._get_input("Should I delete your home directory?",
            ....:         ("yes","no","maybe"), default=0)
            Should I delete your home directory? [Yes/no/maybe]
            sage: print r
            yes
        """
        old_input_func = kwds.get('input_func')
        def input_func(prompt):
            ret = "" if not self else self.pop()
            if old_input_func is getpass:
                self.show(prompt)
            else:
                self.show(prompt+ret)
            return ret
        kwds['input_func'] = input_func

        return super(DoctestInterface, self)._get_input(*args, **kwds)
