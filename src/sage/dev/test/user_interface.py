class DoctestInterface(CmdLineInterface, list):
    def __init__(self, *args, **kwds):
        CmdLineInterface.__init__(self, *args, **kwds)
        self._edits = 0

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

    def edit(self, filename):
        lines = open(filename).read()
        self._edits += 1
        lines = lines.replace("Summary: \n", "Summary: ticket number %s\n" % self._edits)
        lines = lines.replace("ADD DESCRIPTION", "A description after edit number %s" % self._edits)
        open(filename, 'w').write(lines)
