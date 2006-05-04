import os

def install_package(package=None):
    """
    Install a package or return a list of all packages
    that have been installed into this SAGE install.

    You must have an internet connection.  Also, you will have to
    restart SAGE for the changes to take affect.

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
    if package is None:
        X = os.popen('sage -f').read().split('\n')
        i = X.index('Currently installed packages:')
        X = [Y for Y in X[i+1:] if Y != '']
        X.sort()
        return X
    os.system('sage -i "%s"'%package)

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
    X = os.popen('sage -optional').read().split('\n')
    i = X.index('INSTALLED:')
    j = X.index('NOT INSTALLED:')
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
    os.system('sage -upgrade')
    print "You should quit and restart SAGE now."


def package_mesg(package_name):
    mesg  = 'To install the package %s type install_package("%s")\n'%(package_name, package_name)
    mesg += 'at the sage prompt.  Note, the version number might\n'
    mesg += 'change; if so, type optional_packages() to see a list \n'
    mesg += 'of possibilities.   All this requires an internet connection.'
    mesg += 'For more help, type optional_packages?'
    return mesg

