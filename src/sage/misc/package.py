r"""
Sage package management commands

A Sage package has the extension .spkg. It is a tarball that is
(usually) bzip2 compressed that contains arbitrary data and an
spkg-install file. An Sage package typically has the following
components:


-  spkg-install - shell script that is run to install the package

-  Sage.txt - file that describes how the package was made, who
   maintains it, etc.

-  sage - directory with extra patched version of files that needed
   during the install


Use the :func:``optional_packages`` command to list all
optional packages available on the central Sage server.

Actually installing the packages should be done via the Sage command
line, using the following commands:

- ``sage -i PACKAGE_NAME`` -- install the given package

- ``sage -f PACKAGE_NAME`` -- re-install the given package, even if it
  was already installed
"""

import os

__installed_packages = None


def install_package(package=None, force=None):
    """
    Return a list of all packages that have been installed into this
    Sage install.

    EXAMPLES::

        sage: install_package()
        [...'atlas...'python...]

    .. seealso::

        :func:`standard_packages`, :func:`optional_packages`, :func:`experimental_packages`
    """
    if package is not None:
        from sage.misc.stopgap import stopgap
        stopgap("Installing Sage packages using 'install_package()' is obsolete.\nRun 'sage -i {}' from a shell prompt instead.".format(package), 16759)
        return

    global __installed_packages
    if __installed_packages is None:
        import sage.env
        __installed_packages = sorted(os.listdir(sage.env.SAGE_SPKG_INST))
    return __installed_packages


def is_package_installed(package):
    """
    Return true if ``package`` is installed.

    EXAMPLES::

        sage: is_package_installed('pari')
        True

    Giving just the beginning of the package name is not good enough::

        sage: is_package_installed('matplotli')
        False

    Otherwise, installing "pillow" will cause this function to think
    that "pil" is installed, for example.
    """
    return any(p.split('-')[0] == package for p in install_package())

def package_versions(package_type, local=False):
    r"""
    Return version information for each Sage package.

    INPUT:

    - ``package_type`` (string) -- one of `"standard"`, `"optional"` or
      `"experimental"`

    - ``local`` (boolean) -- only query local data (no internet needed)

    For packages of the given type, return a dictionary whose entries
    are of the form ``'package': (installed, latest)``, where
    ``installed`` is the installed version (or ``None`` if not
    installed) and ``latest`` is the latest available version. If the
    package has a directory in ``SAGE_ROOT/build/pkgs/``, then
    ``latest`` is determined by the file ``package-version.txt`` in
    that directory.  If ``local`` is False, then Sage's servers are
    queried for package information.

    EXAMPLES::

        sage: std = package_versions('standard', local=True)
        sage: 'gap' in std
        True
        sage: std['zn_poly']
        ('0.9.p11', '0.9.p11')
    """
    if package_type not in ['standard','optional','experimental']:
        raise ValueError("'package_type' must be one of 'standard','optional','experimental'.")

    cmd = 'sage-list-packages {} --dump'.format(package_type)
    if local:
        cmd += " --local"
    X = os.popen(cmd).read().split('\n')[:-1]

    versions = {}
    for line in X:
        line = line.split(' ')
        installed = line[2]
        if installed == 'not_installed': installed = None
        versions[line[0]] = (installed, line[1])
    return versions

def _package_lists_from_sage_output(package_type, local=False):
    r"""
    Helper function for :func:`standard_packages`, :func:`optional_packages` and
    :func:`experimental_packages`.

    INPUT:

    - ``package_type`` (string) -- one of `"standard"`, `"optional"` or
      `"experimental"`

    - ``local`` (boolean) -- only query local data (no internet needed)

    OUTPUT:

    The function returns a pair of lists ``(installed,not_installed)``
    with the corresponding packages' name, sorted alphabetically. If
    ``local`` is ``True``, then the list of all packages is downloaded
    from the server; otherwise, the list is extracted from the
    packages in ``SAGE_ROOT/build/pkgs/`

    EXAMPLE::

        sage: from sage.misc.package import standard_packages
        sage: installed, not_installed = standard_packages()  # optional internet

    Local check that all standard packages are installed::

        sage: from sage.misc.package import _package_lists_from_sage_output
        sage: installed, not_installed = _package_lists_from_sage_output('standard',local=True)
        sage: bool(not_installed)
        False
        sage: 'glpk' in installed
        True
    """
    installed     = []
    not_installed = []
    versions = package_versions(package_type, local)
    for p in versions:
        if versions[p][0] is None:
            not_installed.append(p)
        else:
            installed.append(p)

    return sorted(installed), sorted(not_installed)

def standard_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed standard packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed standard packages (as a list)

    -  NOT installed standard packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    EXAMPLE::

        sage: from sage.misc.package import standard_packages
        sage: installed, not_installed = standard_packages() # optional internet
        sage: installed[0], installed[-1]                    # optional internet
        ('atlas', 'zn_poly')
        sage: 'mercurial' in not_installed                   # optional internet
        True

    """
    return _package_lists_from_sage_output('standard')

def optional_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed optional packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed optional packages (as a list)

    -  NOT installed optional packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    EXAMPLE::

        sage: from sage.misc.package import optional_packages
        sage: installed, not_installed = optional_packages() # optional internet
        sage: min(installed+not_installed)                   # optional internet
        '4ti2'
        sage: max(installed+not_installed)                   # optional internet
        'zeromq'
    """
    return _package_lists_from_sage_output('optional')

def experimental_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed experimental packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed experimental packages (as a list)

    -  NOT installed experimental packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    EXAMPLE::

        sage: from sage.misc.package import experimental_packages
        sage: installed, not_installed = experimental_packages() # optional internet
        sage: min(installed+not_installed)                   # optional internet
        'PyQt4'
        sage: max(installed+not_installed)                   # optional internet
        'yassl'
    """
    return _package_lists_from_sage_output('experimental')


def upgrade():
    """
    Obsolete function, run 'sage --upgrade' from a shell prompt instead.
    """
    from sage.misc.stopgap import stopgap
    stopgap("Upgrading Sage using 'upgrade()' is obsolete.\nRun 'sage --upgrade' from a shell prompt instead.", 16759)
    return


class PackageNotFoundError(RuntimeError):
    """
    This class defines the exception that should be raised when a
    function, method, or class cannot detect a Sage package that it
    depends on.

    This exception should be raised with a single argument, namely
    the name of the package.

    When an ``PackageNotFoundError`` is raised, this means one of the
    following:

    - The required optional package is not installed.

    - The required optional package is installed, but the relevant
      interface to that package is unable to detect the package.

    EXAMPLES::

        sage: from sage.misc.package import PackageNotFoundError
        sage: raise PackageNotFoundError("my_package")
        Traceback (most recent call last):
        ...
        PackageNotFoundError: the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell
    """
    def __str__(self):
        """
        Return the actual error message.

        EXAMPLES::

            sage: from sage.misc.package import PackageNotFoundError
            sage: str(PackageNotFoundError("my_package"))
            "the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell"
        """
        return ("the package {0!r} was not found. "
            "You can install it by running 'sage -i {0}' in a shell"
            .format(self.args[0]))
