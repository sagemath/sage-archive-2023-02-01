"""
Automated database installation.

This file defines a command "database_install":
           -- takes name of db as argument
           -- if db sitting in SAGE_ROOT, installs it
           -- if not, downloads it with wget from modular
               and installs it.
"""


from sage.misc.misc import SAGE_SHARE

class GenericDatabaseInstaller:
    def __init__(self, name):
        self._name = name

    def __repr__(self):
        return "Database installer for database %s"%self._name

    def name(self):
        return self._name

    def directory(self):
        """
        Returns the directory that contains this database.
        """
        return "%s/%s"%(SAGE_SHARE, self._name)

    def archive_filename(self):
        """
        Returns the filename of the database archive.
        """
        return 'db-%s.tar'%self._name

    def get_archive_file(self):
        """
        Makes sure that the archive file is in the SAGE_ROOT
        directory.
        """
        filename = self.archive_filename()

    def install(self):
        F = self.archive_filename()
        raise NotImplementedError



def database_install(name):
    """
    Install the database name.

    INPUT:
        name -- string
    OUTPUT:
        installs the database so it is available to Sage.
        (You may have to restart Sage.)
    """
    i = name.find('.')
    if i != -1:
        name = name[:i]
        print("Truncating database name to '%s'"%name)
    D = GenericDatabaseInstaller(name)
    D.install()

