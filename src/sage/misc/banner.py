r"""nodoctest
SAGE version and banner info
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


def version(clone = False):
    """
    Return the version of SAGE.

    INPUT:
       nothing
    OUTPUT:
       str

    EXAMPLES:
       sage: version()
       'SAGE Version ..., Release Date: ...'
       sage: version(clone=True)
       ('SAGE Version ..., Release Date: ...',
        'Mercurial clone branch: ...')
    """
    import os
    branch = os.popen("ls -l devel/sage").read().split()[-1][5:]
    import sage.version
    v = 'SAGE Version %s, Release Date: %s'%(sage.version.version, sage.version.date)
    if clone:
        return v,"Mercurial clone branch: %s"%branch
    return v

def banner_text():
    bars = "-"*70
    s = bars
    s += "\n| %-66s |\n"%version()
    s += "| %-66s |\n"%'Type notebook() for the GUI, and license() for information.'
    #s += "| %-66s |\n"%'Distributed under the GNU General Public License V2.'
    s += bars
    return s


def banner():
    print banner_text()


