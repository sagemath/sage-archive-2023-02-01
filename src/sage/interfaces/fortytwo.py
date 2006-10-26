"""
Interface to 4ti2.
"""

import os

import sage.misc.misc

seq = 0
tmp = None

def tmp_dir():
    global tmp
    if tmp is None:
        tmp = sage.misc.misc.tmp_dir('4ti2')
    return tmp

bin='%s/local/lib/4ti2/'%os.environ['SAGE_ROOT']

class FortyTwo:
    def __init__(self, x):
        global seq
        self._seq = seq
        seq = seq + 1
        self._dir = tmp_dir()
        self._x = x

    def _input_text(self):
        try:
            return self.__input_text
        except AttributeError:
            pass
        x = self._x
        nr = x.nrows()
        nc = x.ncols()
        v  = x.list()
        s = '%s %s\n%s'%(nr,nc,v)
        s = s.replace(',','').replace('[','').replace(']','')
        s = self.__input_text
        return s

    def groebner(self):
        try:
            return self.__groebner
        except AttributeError:
            pass
        s = self._input_text()
        open('%s/%s'%(self._dir, self._seq), 'w').write(s)
        e = os.system('cd "%s"; "%s"/groebner %s'%(self._dir, bin, self._seq))
        if e:
            raise RuntimeError
        self.__groebner = open('%s/%s.gro'%(self._dir, self._seq)).read()
        return self.__groebner
