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
