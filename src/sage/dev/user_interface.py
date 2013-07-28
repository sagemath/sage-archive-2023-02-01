r"""
User Interface

This file provides an abstract base class that is used for displaying messages
and prompting for user input.

.. SEEALSO::

    :class:`CmdLineInterface` for an implementation of this class

AUTHORS:

- TODO: add authors from github's history and trac's history

"""
#*****************************************************************************
#       Copyright (C) 2013 TODO
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

class UserInterface(object):
    r"""
    An abstract base class for displaying messages and prompting the user for
    input.

    TESTS::

        sage: from sage.dev.user_interface import UserInterface
        sage: UI = UserInterface()
        sage: type(UI)

    """
    def select(self, prompt, options, default=None):
        r"""
        Ask the user to select from list of options and return selected option.

        INPUT:

        - ``prompt`` -- a string, prompt to display

        - ``options`` -- iterable of strings, the options

        - ``default`` -- a string or ``None`` (default: ``None``), the default
          option

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: UserInterface().select("Should I delete your home directory?", ("yes","no","maybe"), default=1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    def confirm(self, question, default_no=False):
        r"""
        Ask a yes/no question and return the response as a boolean.

        INPUT:

        - ``question`` -- a string

        - ``default_no`` -- boolean (default: ``False``), whether to default to
          no

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: UserInterface().confirm("Should I delete your home directory?")
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        return self.select(question, ("yes", "no"), int(default_no)) == "yes"

    def get_input(self, prompt):
        r"""
        Read input after displaying ``prompt``.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: UserInterface().get_input("What do you want for dinner?")
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    def get_password(self, prompt):
        r"""
        Ask for a password after displaying ``prompt``.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: UserInterface().get_password("What is the passphrase for your safe?")
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    def show(self, message):
        r"""
        Display ``message``.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: UserInterface().show("I ate filet mignon for dinner.")
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    def edit(self, filename):
        r"""
        Drop user into an editor with ``filename`` open.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: UserInterface().edit("filename")
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError
