"""
Base class for Dedekind domains
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

from sage.misc.superseded import deprecation
deprecation(20011, "the module sage.rings.dedekind_domain is deprecated and will be removed")

from sage.rings.ring import DedekindDomain

def is_DedekindDomain(R):
    """
    Check to see if ``R`` is a :class:`DedekindDomain`.

    EXAMPLES::

        sage: import sage.rings.dedekind_domain
        doctest:...: DeprecationWarning: the module sage.rings.dedekind_domain is deprecated and will be removed
        See http://trac.sagemath.org/20011 for details.
        sage: sage.rings.dedekind_domain.is_DedekindDomain(DedekindDomain(QQ))
        True
    """
    return isinstance(R, DedekindDomain)
