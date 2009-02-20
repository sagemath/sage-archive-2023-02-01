r"""
SAGE package management commands

A SAGE package has the extension .spkg.  It is a tarball that is
(usually) bzip2 compressed that contains arbitrary data and an
spkg-install file.  An SAGE package typically has the following
components:
\begin{itemize}
  \item spkg-install -- shell script that is run to install the package
  \item SAGE.txt -- file that describes how the package was made, who maintains it, etc.
  \item sage -- directory with extra patched version of files that needed during the install
\end{itemize}

Use the \code{install_package} command to install a new package, and use
\code{optional_packages} to list all optional packages available on the
central SAGE server.   The \code{upgrade} command upgrades all \emph{standard}
packages -- there is no auto-upgrade command for optional packages.

All package management can also be done via the SAGE command line.
"""

import os

__installed_packages = None


def install_all_optional_packages(force=True, dry_run=False):
    """
    Install all available optional spkg's in the official Sage spkg
    repository.  Returns a list of all spkg's that *fail* to install.

    INPUT:
        force -- bool (default: True); whether to force reinstall of
                 spkg's that are already installed.
        dry_run -- bool (default: False); if True, just list the
                   packages that would be installed in order, but
                   don't actually install them.
    OUTPUT:
        list of strings

    NOTE: This is designed mainly for testing purposes.  This also
    doesn't do anything with respect to dependencies -- the packages
    are installed in alphabetical order.  Dependency issues will be
    dealt with in a future version.

    AUTHOR:
        -- William Stein (2008-12)

    EXAMPLES:
        sage: sage.misc.package.install_all_optional_packages(dry_run=True)  # optional - internet
        Installing ...
        []
    """
    # Get list of all packages from the package server
    installed, not_installed = optional_packages()
    if force:
        v = installed + not_installed
    else:
        v = not_installed
    failed = []

    # sort the list of packages in alphabetical order
    v.sort()

    # install each one
    for pkg in v:
        try:
            print "Installing %s..."%pkg
            if not dry_run:
                # only actually do the install of the package if dry_run is False
                install_package(pkg, force=force)
        except ValueError, msg:
            # An error occured -- catch exception and record this in
            # the list of failed installed.
            print "*"*70
            print "FAILED to install '%s'"%pkg
            print "*"*70
            failed.append(pkg)
    return failed

def install_package(package=None, force=False):
    """
    Install a package or return a list of all packages
    that have been installed into this SAGE install.

    You must have an internet connection.  Also, you will have to
    restart SAGE for the changes to take affect.

    It is not needed to provide the version number.

    INPUT:
        package -- optional; if specified, install the
                   given package.  If not, list all
                   installed packages.

    IMPLEMENTATION: calls 'sage -f'.

    RELATED COMMANDS:
        optional_packages -- list of all optional packages
        upgrade -- upgrade to latest version of core packages
                   (optional packages are not automatically upgraded).
    """
    global __installed_packages
    if os.uname()[0][:6] == 'CYGWIN':
        print "install_package may not work correctly under Microsoft Windows"
        print "since you can't change an opened file.  Quit all"
        print "instances of sage and use 'sage -i' instead or"
        print "use the force option to install_package."
        return
    if package is None:
        if __installed_packages is None:
            X = os.popen('sage -f').read().split('\n')
            i = X.index('Currently installed packages:')
            X = [Y for Y in X[i+1:] if Y != '']
            X.sort()
            __installed_packages = X
        return __installed_packages
    # Get full package name
    if force:
        S = [P for P in standard_packages()[0] if P.startswith(package)]
        O = [P for P in optional_packages()[0] if P.startswith(package)]
        E = [P for P in experimental_packages()[0] if P.startswith(package)]
    else:
        S,O,E = [], [], []
    S.extend([P for P in standard_packages()[1] if P.startswith(package)])
    O.extend([P for P in optional_packages()[1] if P.startswith(package)])
    E.extend([P for P in experimental_packages()[1] if P.startswith(package)])
    L = S+O+E
    if len(L)>1:
        if force:
            print "Possible package names starting with '%s' are:"%(package)
        else:
            print "Possible names of non-installed packages starting with '%s':"%(package)
        for P in L:
            print " ", P
        raise ValueError, "There is more than one package name starting with '%s'. Please specify!"%(package)
    if len(L)==0:
        if not force:
            if is_package_installed(package):
                raise ValueError, "Package is already installed. Try install_package('%s',force=True)"%(package)
        raise ValueError, "There is no package name starting with '%s'."%(package)
    os.system('sage -f "%s"'%(L[0]))
    __installed_packages = None


def is_package_installed(package):
    """
    Return true if a package starting with the given string is installed.

    EXAMPLES:
        sage: is_package_installed('sage')
        True
    """
    return any(p.startswith(package) for p in install_package())

def standard_packages():
    """
    Return two lists.  The first contains the installed and the second
    contains the not-installed standard packages that are available
    from the SAGE repository.      You must have an internet connection.

    OUTPUT:
        -- installed standard packages (as a list)
        -- NOT installed standard packages (as a list)

    Use \code{install_package(package_name)} to install or re-install
    a given package.

    RELATED COMMANDS:
        install_package -- list of all standard packages
        upgrade -- upgrade to latest version of core packages
                   (standard packages are not automatically upgraded).
    """
    R = os.popen('sage -standard').read()
    X = R.split('\n')
    try:
        i = X.index('INSTALLED:')
        j = X.index('NOT INSTALLED:')
    except ValueError, msg:
        print R
        print "standard package list (shown above) appears to be currently not available or corrupted (network error?)."
        return [], []

    installed = []
    for k in X[i+1:]:
        if k == '':
            break
        installed.append(k)

    not_installed = []
    for k in X[j+1:]:
        if k == '':
            break
        not_installed.append(k)
    return installed, not_installed

def optional_packages():
    """
    Return two lists.  The first contains the installed and the second
    contains the not-installed optional packages that are available
    from the SAGE repository.      You must have an internet connection.

    OUTPUT:
        -- installed optional packages (as a list)
        -- NOT installed optional packages (as a list)

    Use \code{install_package(package_name)} to install or re-install
    a given package.

    RELATED COMMANDS:
        install_package -- list of all optional packages
        upgrade -- upgrade to latest version of core packages
                   (optional packages are not automatically upgraded).
    """
    R = os.popen('sage -optional').read()
    X = R.split('\n')
    try:
        i = X.index('INSTALLED:')
        j = X.index('NOT INSTALLED:')
    except ValueError, msg:
        print R
        print "Optional package list (shown above) appears to be currently not available or corrupted (network error?)."
        return [], []

    installed = []
    for k in X[i+1:]:
        if k == '':
            break
        installed.append(k)

    not_installed = []
    for k in X[j+1:]:
        if k == '':
            break
        not_installed.append(k)
    return installed, not_installed

def experimental_packages():
    """
    Return two lists.  The first contains the installed and the second
    contains the not-installed experimental packages that are available
    from the SAGE repository.      You must have an internet connection.

    OUTPUT:
        -- installed experimental packages (as a list)
        -- NOT installed experimental packages (as a list)

    Use \code{install_package(package_name)} to install or re-install
    a given package.

    RELATED COMMANDS:
        install_package -- list of all experimental packages
        upgrade -- upgrade to latest version of core packages
                   (experimental packages are not automatically upgraded).
    """
    R = os.popen('sage -experimental').read()
    X = R.split('\n')
    try:
        i = X.index('INSTALLED:')
        j = X.index('NOT INSTALLED:')
    except ValueError, msg:
        print R
        print "experimental package list (shown above) appears to be currently not available or corrupted (network error?)."
        return [], []

    installed = []
    for k in X[i+1:]:
        if k == '':
            break
        installed.append(k)

    not_installed = []
    for k in X[j+1:]:
        if k == '':
            break
        not_installed.append(k)
    return installed, not_installed

#################################################################
# Upgrade to latest version of SAGE.
#################################################################


def upgrade():
    """
    Download and build the latest version of SAGE.

    You must have an internet connection.  Also, you will have to
    restart SAGE for the changes to take affect.

    This upgrades to the latest version of core packages (optional
    packages are not automatically upgraded).

    This will not work on systems that don't have a C compiler.

    RELATED COMMANDS:
        install_package -- list of all optional packages
        optional_packages -- list of all optional packages
    """
    global __installed_packages
    if os.uname()[0][:6] == 'CYGWIN':
        print "Upgrade may not work correctly under Microsoft Windows"
        print "since you can't change an opened file.  Quit all"
        print "instances of sage and use 'sage -upgrade' instead."
        return
    os.system('sage -upgrade')
    __installed_packages = None
    print "You should quit and restart SAGE now."


def package_mesg(package_name):
    mesg  = 'To install the package %s type install_package("%s")\n'%(package_name, package_name)
    mesg += 'at the sage prompt.  Note, the version number might\n'
    mesg += 'change; if so, type optional_packages() to see a list \n'
    mesg += 'of possibilities.   All this requires an internet connection.'
    mesg += 'For more help, type optional_packages?'
    return mesg

