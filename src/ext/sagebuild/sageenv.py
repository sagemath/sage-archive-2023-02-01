import os, signal, sys, time, thread, threading, tempfile
from build.all import *

class sage_env(enviroment):
    def __init__(self):
        enviroment.__init__(self)
        self.choose_blas()
        self.setup_options()
        self.setup_defaults()

    def choose_blas(self):
        ## Choose cblas library -- note -- make sure to update sage/misc/cython.py
        ## if you change this!!
        if os.environ.has_key('SAGE_BLAS'):
            self.options['BLAS']=os.environ['SAGE_BLAS']
            self.options['BLAS2']=os.environ['SAGE_BLAS']
        elif os.path.exists('%s/lib/libatlas.so'%os.environ['SAGE_LOCAL']):
            self.options['BLAS']='cblas'
            self.options['BLAS2']='atlas'
        elif os.path.exists('/usr/lib/libcblas.dylib') or \
             os.path.exists('/usr/lib/libcblas.so'):
            self.options['BLAS']='cblas'
            self.options['BLAS2']='atlas'
        elif os.path.exists('/usr/lib/libblas.dll.a'):
            self.options['BLAS']='gslcblas'
            self.options['BLAS2']='gslcblas'
        else:
            # This is very slow  (?), but *guaranteed* to be available.
            self.options['BLAS']='gslcblas'
            self.options['BLAS2']='gslcblas'

    def setup_options(self):
        if not os.environ.has_key('SAGE_ROOT'):
            print "    ERROR: The environment variable SAGE_ROOT must be defined."
            sys.exit(1)
        else:
            SAGE_ROOT  = os.environ['SAGE_ROOT']
            self.options['SAGE_ROOT'] = SAGE_ROOT
            SAGE_LOCAL = self.options['SAGE_ROOT'] + '/local/'
            self.options['SAGE_LOCAL'] = SAGE_LOCAL
            SAGE_DEVEL = self.options['SAGE_ROOT'] + '/devel/'
            self.options['SAGE_DEVEL'] = SAGE_DEVEL
        if not os.environ.has_key('SAGE_VERSION'):
            self.options['SAGE_VERSION']=0
        else:
            self.options['SAGE_VERSION'] = os.environ['SAGE_VERSION']
        self.options['SITE_PACKAGES'] = '%slib/python/site-packages/'%self.options['SAGE_LOCAL']
        if not os.path.exists(self.options['SITE_PACKAGES']):
            self.options['SITE_PACKAGES'] = '%slib/python2.5/site-packages/'%self.options['SAGE_LOCAL']
            if not os.path.exists(self.options['SITE_PACKAGES']):
                self.options['SITE_PACKAGES'] = '%slib/python2.4/site-packages/'%self.options['SAGE_LOCAL']
                if not os.path.exists(self.options['SITE_PACKAGES']) and os.environ.has_key('SAGE_DEBIAN'):
                    self.options['SITE_PACKAGES'] = '/usr/lib/python2.5/site-packages/'
                if not os.path.exists(self.options['SITE_PACKAGES']):
                    raise RuntimeError, "Unable to find site-packages directory (see setup.py file in sage python code)."
        if not self.options.has_key('SAGE_DEBIAN'):
            self.options['SAGE_DEBIAN'] = False
        self.options['SITE_PACKAGES_REL'] = self.options['SITE_PACKAGES'][len(self.options['SAGE_LOCAL']):]
        try:
            self.options['UNAME'] = os.environ['UNAME']
        except:
            pass
    def setup_defaults(self):
        self.defaults['include_dirs'] = ['%sinclude'%self.options['SAGE_LOCAL'], \
        '%sinclude/csage'%self.options['SAGE_LOCAL'], \
        '%ssage/sage/ext'%self.options['SAGE_DEVEL'], \
        '%sinclude/python2.5'%self.options['SAGE_LOCAL']]

        if self.options['SAGE_DEBIAN']:
            self.debian_include_dirs=["/usr/include","/usr/include/numpy","/usr/include/FLINT","/usr/include/givaro", "/usr/include/gsl","/usr/include/fplll","/usr/include/eclib","/usr/include/gmp++","/usr/include/linbox","/usr/include/NTL","/usr/include/pari","/usr/include/qd","/usr/include/singular","/usr/include/singular/singular","/usr/include/symmetrica","/usr/include/polybori"]
            self.defaults['include_dirs'].expand(self.debian_include_dirs)
        else:
            self.debian_include_dirs=[]
        self.defaults['prelibraries'] = ['csage']
        self.defaults['postlibraries'] = ['stdc++', 'ntl']
        if self.options['SAGE_DEBIAN']:
            self.defaults['library_dirs'] =  ['/usr/lib','/usr/lib/eclib','/usr/lib/singular','/usr/lib/R/lib','%s/lib' % self.options['SAGE_LOCAL'] ]
        else:
            self.defaults['library_dirs'] = ['%slib' % self.options['SAGE_LOCAL'] ]

class debian_env(sage_env):
    def __init__(self):
        self.options = { }
        self.options['SAGE_DEBIAN'] = True
        sage_env.__init__(self)

def create_env():
    if 'SAGE_DEBIAN' in os.environ:
        env = debian_env()
    else:
        env = sage_env()
    return env
