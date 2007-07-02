"""nodoctest
Configuration
"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

class Configuration:

    def __init__(self):
        self.confs = {}

    def __repr__(self):
        return 'Configuration: %s'%self.confs

    def defaults(self):
        raise NotImplementedError

    def __getitem__(self, key):
        try:
            return self.confs[key]
        except KeyError:
            if self.defaults().has_key(key):
                A = self.defaults()[key]
                self.confs[key] = A
                return A
            else:
                raise KeyError, "No key '%s' and no default for this key"%key

    def __setitem__(self, key, value):
        if not isinstance(value, (str, int)):
            # Make sure the input is a string or Python int.  Nothing else is allowed.
            # This will make having a text file or database format, etc., later
            # much easier.
            try:
                value = int(value)
            except (TypeError, ValueError):
                value = str(value)
        self.confs[key] = value

    def html_conf_form(self, action):
        D = self.defaults()
        C = self.confs
        K = list(set(self.confs.keys() + D.keys()))
        K.sort()
        options = ''
        for key in K:
            options += '<tr><td>%s</td><td><input type="text" name="%s" value="%s"></td></tr>\n'%(key, key, self[key])
        s = """
        <form method="post" action="%s" enctype="multipart/form-data">
        <input type="submit" value="Submit">
        <table border=0 cellpadding=5 cellspacing=2>
%s
        </table>
        </form>
        """%(action, options)
        return s
