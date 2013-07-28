"""
User Interface

Abstract base class that is used for displaying messages and prompting for user
input.
"""

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
