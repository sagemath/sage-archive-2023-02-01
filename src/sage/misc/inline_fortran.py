import __builtin__
import tempfile
import os

count=0


class InlineFortran:
    def __init__(self,globals):
        self.globals=globals
        self.library_paths=[]
        self.libraries=[]
        self._verbose = False

    def __repr__(self):
        return "Interface to Fortran compiler"

    def eval(self,x):
        import os
        global count
        # On linux g77_shared should be a script that runs gfortran -shared
        # On OSX it should be a script that runs gfortran -bundle -undefined dynamic_lookup
        path = os.environ['SAGE_LOCAL']+'/bin/sage-g77_shared'
        from numpy import f2py
        name='fortran_module_%d'%count
        s_lib_path=""
        s_lib=""
        for s in self.library_paths:
            s_lib_path=s_lib_path+"-L"+s+" "

        for s in self.libraries:
            s_lib=s_lib +"-l"+s + " "

        # if the first line has !f90 as a commment gfortran will treat it as
        # fortran 90 code
        if x[0:4]=='!f90':
            fname = os.path.join(tempfile.mktemp()+'.f90')
        else:
            fname = os.path.join(tempfile.mktemp()+'.f')

        f2py.compile(x, name, '--quiet --f77exec=%s --f90exec=%s'%(path,path)+ \
                     ' '+s_lib_path+' '+s_lib,source_fn=fname)
        if self._verbose:
            print "OUTPUT from Fortran"
            print f.read()

        count += 1
        m=__builtin__.__import__(name)
        for k, x in m.__dict__.iteritems():
            if k[0] != '_':
                self.globals[k] = x

    def verbose(self, option):
        """
        EXAMPLES:
            sage: fortran.verbose(False)
            sage: fortran.verbose(True)
        """
        self._verbose = bool(option)

    def add_library(self,s):
       self.libraries.append(s)

    def add_library_path(self,s):
       self.library_paths.append(s)

