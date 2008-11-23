#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
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

class MagmaExpr(str):
    def __repr__(self):
        return str(self)

def magma_free_eval(code, strip=True, columns=0):
    """
    Use the free online MAGMA calculator to evaluate the given
    input code and return the answer as a string.

    LIMITATIONS: The code must evaluate in at most 20 seconds
    and there is a limitation on the amount of RAM.

    EXAMPLES:
        sage: magma_free("Factorization(9290348092384)")  # optional - internet
        [ <2, 5>, <290323377887, 1> ]
    """
    import urllib
    url = "http://magma.maths.usyd.edu.au/calc/"
    code = "SetColumns(%s);\n"%columns + code
    urldata = urllib.urlencode({'input':code})
    results = urllib.urlopen(url, urldata).read()
    if strip:
        i = results.find('-----\n\n') + 7
        j = results.rfind('\n\nTotal time')
    else:
        i = results.find('Magma V')
        j = results.rfind('</textarea>')
    class MagmaExpr(str):
       def __repr__(self):
          return str(self)
    return MagmaExpr(results[i:j])

class MagmaFree:
    """
    Evaluate MAGMA code without requiring that MAGMA be installed
    on your computer by using the free online MAGMA calculator.

    EXAMPLES:
        sage: magma_free("Factorization(9290348092384)")  # optional - internet
        [ <2, 5>, <290323377887, 1> ]
    """
    def eval(self, x):
        return magma_free_eval(x)
    def __call__(self, code, strip=True, columns=0):
        return magma_free_eval(code, strip=strip, columns=columns)

magma_free = MagmaFree()
