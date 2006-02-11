r"""
Interface to mwrank
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

import os, weakref
from expect import Expect

instances={}
def Mwrank(options="", server=None):
    """
    Create and return an mwrank interpreter, with given options.

    INPUT:
       options -- string; passed when starting mwrank.  The format is
       q p<precision> v<verbosity> b<hlim_q> x<naux>  c<hlim_c> l t o s d>]
    """
    global instances
    try:
        X = instances[options]()
        if X:
            return X
    except KeyError:
        pass
    X = Mwrank_class(options, server=server)
    instances[options] = weakref.ref(X)
    return X


class Mwrank_class(Expect):
    """
    Interface to the Mwrank interpreter.
    """
    def __init__(self, options="", server=None):
        """
        INPUT:
           options -- string; passed when starting mwrank.  The format is
           q p<precision> v<verbosity> b<hlim_q> x<naux>  c<hlim_c> l t o s d>]
        """
        Expect.__init__(self,
                        name = 'mwrank',
                        prompt = 'Enter curve: ',
                        command = "mwrank %s"%options,
                        server = server,
                        maxread = 10000,
                        restart_on_ctrlc = True,
                        verbose_start = False)

    def __reduce__(self):
        return reduce_load_Mwrank, tuple([])

    def __call__(self, cmd):
        return self.eval(cmd)

    def console(self):
        mwrank_console()


# An instance
mwrank = Mwrank()

def reduce_load_mwrank():
    return mwrank

import os
def mwrank_console():
    os.system('mwrank')

