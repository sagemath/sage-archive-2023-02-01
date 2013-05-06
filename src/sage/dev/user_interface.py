"""
User Interface
"""
from getpass import getpass

class CmdLineInterface(object):
    def __init__(self):
        self._answer_stack = []

    def get_input(self, prompt, options=None, default=None, dryrun=False, strip=None):
        """
        Get input from the developer.

        INPUT:

        - ``prompt`` -- a string
        - ``options`` -- a list of strings or None
        - ``default`` -- a string or None
        - ``strip`` -- a string to strip off the beginning of options
        - ``dryrun`` -- boolean

        EXAMPLES::

            sage: UI = sage.dev.user_interface.CmdLineInterface()
            sage: UI.get_input("Should I delete your home directory?",
            ... ["yes","no","maybe"], default="m", dryrun=True)
            Should I delete your home directory? [yes/no/Maybe]
            'maybe'
        """
        if default is not None:
            for i, opt in enumerate(options):
                if opt.startswith(default):
                    default = i
                    break
            else:
                default = None
        if options is not None:
            options = list(options)
            if len(options) > 0:
                for i, option in enumerate(options):
                    if i == default:
                        options[i] = str(option).capitalize()
                    else:
                        options[i] = str(option).lower()
                prompt += " [" + "/".join(options) + "] "
                if default is not None:
                    options[default] = options[default].lower()
                if strip is not None:
                    for i, option in enumerate(options):
                        if option.startswith(strip):
                            options[i] = option[len(strip):]
            else:
                prompt += " "
                options = None
        if dryrun:
            print prompt
            return options[default]
        while True:
            if not self._answer_stack:
                s = raw_input(prompt)
            else:
                s = self._answer_stack.pop()
                print "%s%s"%(prompt,s)
            if strip is not None and s.startswith(strip):
                s = s[len(strip):]
            if options is None:
                return s
            if len(s.strip()) == 0:
                if default is None:
                    print "Please enter an option"
                    continue
                else:
                    return options[default]
            found = -1
            for i, opt in enumerate(options):
                if opt.startswith(s):
                    if found == -1:
                        found = i
                    else:
                        break
            else:
                if found != -1:
                    return options[found]
            if found == -1:
                print "Please specify an allowable option"
            else:
                print "Please disambiguate between options"

    def confirm(self, question, default_yes=True, dryrun=False):
        """
        Ask a yes/no question from the developer.

        INPUT:

        - ``question`` -- a string
        - ``default_yes`` -- boolean whether to default to yes
        - ``dryrun`` -- boolean

        EXAMPLES::

            sage: UI = sage.dev.user_interface.CmdLineInterface()
            sage: UI.confirm("Should I delete your home directory?",
            ... dryrun=True)
            Should I delete your home directory? [Yes/no]
            True
        """
        ok = self.get_input(question, ["yes","no"],
                            "yes" if default_yes else "no", dryrun)
        return ok == "yes"

    def get_password(self, prompt):
        if not self._answer_stack:
            return getpass(prompt)
        else:
            ret = self._answer_stack.pop()
            print prompt
            return ret

    def show(self, message):
        print message
