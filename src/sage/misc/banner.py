r"""
Sage version and banner info
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.version

def version(clone = False):
    """
    Return the version of Sage.

    INPUT:
       nothing
    OUTPUT:
       str

    EXAMPLES:
       sage: version()
       'Sage Version ..., Release Date: ...'
       sage: version(clone=True)
       ('Sage Version ..., Release Date: ...',
        'Mercurial clone branch: ...')
    """
    import os
    branch = os.popen("ls -l $SAGE_ROOT/devel/sage").read().split()[-1][5:]
    v = 'Sage Version %s, Release Date: %s'%(sage.version.version, sage.version.date)
    if clone:
        return v,"Mercurial clone branch: %s"%branch
    return v


def banner_text():
    """
    Text for the Sage banner.

    INPUT:
        None

    OUTPUT:
        string containing the banner message

    EXAMPLES:
        sage: print sage.misc.banner.banner_text()
        ----------------------------------------------------------------------
        | Sage Version ..., Release Date: ...                  |
        | Type notebook() for the GUI, and license() for information.        |
        ----------------------------------------------------------------------
    """
    bars = "-"*70
    s = bars
    s += "\n| %-66s |\n"%version()
    s += "| %-66s |\n"%'Type notebook() for the GUI, and license() for information.'
    #s += "| %-66s |\n"%'Distributed under the GNU General Public License V2.'
    s += bars
    return s


def banner():
    """
    Print the Sage banner.

    INPUT:
        None

    OUTPUT:
        None

    EXAMPLES:
        sage: banner()
        ----------------------------------------------------------------------
        | Sage Version ..., Release Date: ...                  |
        | Type notebook() for the GUI, and license() for information.        |
        ----------------------------------------------------------------------
    """
    print banner_text()


def version_dict():
    """
    A dictionary describing the version of Sage.

    INPUT:
        nothing
    OUTPUT:
        dictionary with keys 'major', 'minor', 'tiny', 'prerelease'

    This process the Sage version string and produces a dictionary.
    It expects the Sage version to be in one of these forms:
       N.N
       N.N.N
       N.N.N.N
       N.N.str
       N.N.N.str
       N.N.N.N.str
    where 'N' stands for an integer and 'str' stands for a string.
    The first integer is stored under the 'major' key and the second
    integer under 'minor'.  If there is one more integer, it is stored
    under 'tiny'; if there are two more integers, then they are stored
    together as a float N.N under 'tiny'.  If there is a string, then
    the key 'prerelease' returns True.

    For example, if the Sage version is '3.2.1', then the dictionary
    is {'major': 3, 'minor': 2, 'tiny': 1, 'prerelease': False}.
    If the Sage version is '3.2.1.2', then the dictionary is
    {'major': 3, 'minor': 2, 'tiny': 1.2, 'prerelease': False}.
    If the Sage version is '3.2.alpha0', then the dictionary is
    {'major': 3, 'minor': 2, 'tiny': 0, 'prerelease': True}.

    EXAMPLES:
        sage: from sage.misc.banner import version_dict
        sage: print "Sage major version is %s" % version_dict()['major']
        Sage major version is ...
        sage: version_dict()['major'] == int(sage.version.version.split('.')[0])
        True
    """
    from sage.rings.integer import Integer
    v = sage.version.version.split('.')
    dict = {}
    dict['major'] = int(v[0])
    dict['minor'] = int(v[1])
    dict['tiny'] = 0
    dict['prerelease'] = (isinstance(v[-1], str))
    if len(v) >= 3:
        if isinstance(v[2], (int, Integer, long)):
            dict['tiny'] = int(v[2])
            if len(v) >= 4:
                if isinstance(v[3], (int, Integer, long)):
                    dict['tiny'] += 0.1 * v[3]
                else:
                    dict['alpha'] = v[3]
    return dict


def require_version(major, minor=0, tiny=0, prerelease=False,
                    print_message=False):
    """
    True if Sage version is at least major.minor.tiny.

    INPUT:
        major -- integer
        minor -- integer (optional, default = 0)
        tiny -- float (optional, default = 0)
        prerelease -- boolean (optional, default = False)
        print_message -- boolean (optional, default = False)

    OUTPUT:
        True if major.minor.tiny is <= version of Sage, False otherwise

    For example, if the Sage version number is 3.1.2, then
    require_version(3, 1, 3) will return False, while
    require_version(3, 1, 2) will return True.
    If the Sage version is 3.1.2.alpha0, then
    require_version(3, 1, 1) will return True, while, by default,
    require_version(3, 1, 2) will return False.  Note, though, that
    require_version(3, 1, 2, prerelease=True) will return True:
    if the optional argument prerelease is True, then a prerelease
    version of Sage counts as if it were the released version.

    If optional argument print_message is True and this function
    is returning False, print a warning message.

    EXAMPLES:
        sage: from sage.misc.banner import require_version
        sage: require_version(2, 1, 3)
        True
        sage: require_version(821, 4)
        False
        sage: require_version(821, 4, print_message=True)
        This code requires at least version 821.4 of Sage to run correctly.
        You are running version ...
        False
    """
    vers = version_dict()
    prerelease_checked = (prerelease if vers['prerelease'] else True)
    if (vers['major'] > major
        or (vers['major'] == major and vers['minor'] > minor)
        or (vers['major'] == major and vers['minor'] == minor
            and vers['tiny'] > tiny)
        or (vers['major'] == major and vers['minor'] == minor
            and vers['tiny'] == tiny and prerelease_checked)):
        return True
    else:
        if print_message:
            print "This code requires at least version",
            print "%g" % (major + 0.1 * minor + 0.01 * tiny,),
            print "of Sage to run correctly."
            print "You are running version %s." % sage.version.version
        return False
