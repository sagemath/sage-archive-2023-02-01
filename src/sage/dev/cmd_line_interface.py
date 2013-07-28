from __future__ import print_function
from subprocess import check_call, CalledProcessError
from getpass import getpass
import os

class CmdLineInterface(UserInterface):
    def _ioctl_GWINSZ(self, fd):
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
            import struct
            import fcntl
            import termios
        except ImportError:
            return None

        try:
            return struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
        except IOError:
            return None

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
        if sum(len(l) // width + 1 for l in message) + 2 <= height:
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
        try:
            check_call([os.environ.get('EDITOR', 'nano'), filename])
        except CalledProcessError:
            raise RuntimeError("Editor returned non-zero exit value")
