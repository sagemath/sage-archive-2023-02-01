#!/usr/bin/env python

import os, sys
from distutils.core import setup
from distutils.extension import Extension


#########################################################
### Configure SAGE_ROOT
#########################################################

if not os.environ.has_key('SAGE_ROOT'):
    print "    ERROR: The environment variable SAGE_ROOT must be defined."
    sys.exit(1)
else:
    SAGE_ROOT  = os.environ['SAGE_ROOT']
    SAGE_LOCAL = SAGE_ROOT + '/local/'
    SAGE_DEVEL = SAGE_ROOT + '/devel/'

#########################################################
### BLAS setup
#########################################################

## Choose cblas library -- note -- make sure to update sage/misc/cython.py
## if you change this!!
if os.environ.has_key('SAGE_BLAS'):
    BLAS=os.environ['SAGE_BLAS']
    BLAS2=os.environ['SAGE_BLAS']
elif os.path.exists('%s/lib/libatlas.so'%os.environ['SAGE_LOCAL']):
    BLAS='cblas'
    BLAS2='atlas'
elif os.path.exists('/usr/lib/libcblas.dylib') or \
     os.path.exists('/usr/lib/libcblas.so'):
    BLAS='cblas'
    BLAS2='atlas'
elif os.path.exists('/usr/lib/libblas.dll.a'):
    BLAS='gslcblas'
    BLAS2='gslcblas'
else:
    # This is very slow  (?), but *guaranteed* to be available.
    BLAS='gslcblas'
    BLAS2='gslcblas'


#########################################################
### Debian-related stuff
#########################################################

if os.environ.has_key('SAGE_DEBIAN'):
    debian_include_dirs=["/usr/include",
                         "/usr/include/cudd",
                         "/usr/include/eclib",
                         "/usr/include/FLINT",
                         "/usr/include/fplll",
                         "/usr/include/givaro",
                         "/usr/include/gmp++",
                         "/usr/include/gsl",
                         "/usr/include/linbox",
                         "/usr/include/NTL",
                         "/usr/include/numpy",
                         "/usr/include/pari",
                         "/usr/include/polybori",
                         "/usr/include/polybori/groebner",
                         "/usr/include/qd",
                         "/usr/include/singular",
                         "/usr/include/singular/singular",
                         "/usr/include/symmetrica",
                         "/usr/include/zn_poly"]
    include_dirs = include_dirs + debian_include_dirs
else:
    debian_include_dirs=[]

#########################################################
### Commonly used include directories
#########################################################

numpy_include_dirs = [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include']

#############################################################
### List of modules
###
### Note that the list of modules is sorted alphabetically
### by extension name. Please keep this list sorted when
### adding new modules!
###
#############################################################

ext_modules = [

    ################################
    ##
    ## sage.calculus
    ##
    ################################

    Extension('sage.calculus.var',
              sources = ['sage/calculus/var.pyx']),

    ################################
    ##
    ## sage.categories
    ##
    ################################

    Extension('sage.categories.action',
              sources = ['sage/categories/action.pyx']),

    Extension('sage.categories.functor',
              sources = ['sage/categories/functor.pyx']),

    Extension('sage.categories.map',
              sources = ['sage/categories/map.pyx']),

    Extension('sage.categories.morphism',
              sources = ['sage/categories/morphism.pyx']),

    ################################
    ##
    ## sage.coding
    ##
    ################################

    Extension('sage.coding.binary_code',
              sources = ['sage/coding/binary_code.pyx']),

    ################################
    ##
    ## sage.combinat
    ##
    ################################

    Extension('sage.combinat.expnums',
              sources = ['sage/combinat/expnums.pyx'],
              libraries = ['gmp']),

    Extension('sage.combinat.matrices.dancing_links',
              sources = ['sage/combinat/matrices/dancing_links.pyx'],
              libraries = ["stdc++"],
              language='c++'),

    Extension('sage.combinat.partitions',
              sources = ['sage/combinat/partitions.pyx',
                         'sage/combinat/partitions_c.cc'],
              libraries = ['qd', 'gmp', 'mpfr'],
              depends = ['sage/combinat/partitions_c.h'],
              language='c++'),

    ################################
    ##
    ## sage.ext
    ##
    ################################

    Extension('sage.ext.fast_eval',
              sources = ['sage/ext/fast_eval.pyx']),

    Extension('sage.ext.interactive_constructors_c',
              sources = ['sage/ext/interactive_constructors_c.pyx']),

    Extension('sage.ext.multi_modular',
              sources = ['sage/ext/multi_modular.pyx'],
              libraries=['gmp']),

    Extension('sage.ext.sig',
              sources = ['sage/ext/sig.pyx']),

    ################################
    ##
    ## sage.finance
    ##
    ################################

    Extension('sage.finance.fractal',
              sources = ['sage/finance/fractal.pyx']),

    Extension('sage.finance.markov_multifractal_cython',
              sources = ['sage/finance/markov_multifractal_cython.pyx']),

    Extension('sage.finance.time_series',
              sources = ['sage/finance/time_series.pyx'],
              include_dirs = numpy_include_dirs),


    ################################
    ##
    ## sage.graphs
    ##
    ################################
    # sage.graphs

    Extension('sage.graphs.chrompoly',
              sources = ['sage/graphs/chrompoly.pyx']),

    Extension('sage.graphs.graph_fast',
              sources = ['sage/graphs/graph_fast.pyx'],
              libraries = ['gmp']),

    Extension('sage.graphs.graph_isom',
              sources = ['sage/graphs/graph_isom.pyx']),

    Extension('sage.graphs.planarity',
              sources = ['sage/graphs/planarity.pyx',
                         'sage/graphs/planarity/graphEmbed.c',
                         'sage/graphs/planarity/graphIO.c',
                         'sage/graphs/planarity/graphIsolator.c',
                         'sage/graphs/planarity/graphNonplanar.c',
                         'sage/graphs/planarity/graphPreprocess.c',
                         'sage/graphs/planarity/graphStructure.c',
                         'sage/graphs/planarity/graphTests.c',
                         'sage/graphs/planarity/listcoll.c',
                         'sage/graphs/planarity/planarity.c',
                         'sage/graphs/planarity/stack.c'],
              depends = ['sage/graphs/planarity/appconst.h',
                         'sage/graphs/planarity/graph.h',
                         'sage/graphs/planarity/listcoll.h',
                         'sage/graphs/planarity/platformTime.h',
                         'sage/graphs/planarity/stack.h']),

    Extension('sage.graphs.base.sparse_graph',
              sources = ['sage/graphs/base/sparse_graph.pyx']),

        ################################
        ##
        ## sage.graphs.base
        ##
        ################################

    Extension('sage.graphs.base.c_graph',
              sources = ['sage/graphs/base/c_graph.pyx']),

    Extension('sage.graphs.base.dense_graph',
              sources = ['sage/graphs/base/dense_graph.pyx']),

    ################################
    ##
    ## sage.groups
    ##
    ################################

    Extension('sage.groups.group',
              sources = ['sage/groups/group.pyx']),

    Extension('sage.groups.perm_gps.permgroup_element',
              sources = ['sage/groups/perm_gps/permgroup_element.pyx']),

        ###################################
        ##
        ## sage.groups.perm_gps.partn_ref
        ##
        ###################################

    Extension('sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label',
              sources = ['sage/groups/perm_gps/partn_ref/automorphism_group_canonical_label.pyx']),

    Extension('sage.groups.perm_gps.partn_ref.double_coset',
              sources = ['sage/groups/perm_gps/partn_ref/double_coset.pyx']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_binary',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_binary.pyx']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_graphs',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_graphs.pyx']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_matrices',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_matrices.pyx']),

    ################################
    ##
    ## sage.gsl
    ##
    ################################

    Extension('sage.gsl.callback',
              sources = ['sage/gsl/callback.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.dwt',
              sources = ['sage/gsl/dwt.pyx'],
              libraries=['gsl',BLAS],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.fft',
              sources = ['sage/gsl/fft.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.gsl_array',
              sources = ['sage/gsl/gsl_array.pyx'],
              libraries=['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.integration',
              sources = ['sage/gsl/integration.pyx'],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')],
              libraries=['gsl',BLAS, BLAS2]),

    Extension('sage.gsl.interpolation',
              sources = ['sage/gsl/interpolation.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.ode',
              sources = ['sage/gsl/ode.pyx'],
              libraries=['gsl',BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.probability_distribution',
              sources = ['sage/gsl/probability_distribution.pyx'],
              libraries=['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    ################################
    ##
    ## sage.libs
    ##
    ################################

    Extension('sage.libs.flint.flint',
              sources = ["sage/libs/flint/flint.pyx"],
              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],
              include_dirs = [SAGE_ROOT+'/local/include/FLINT/'],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              depends = [SAGE_ROOT + "/local/include/FLINT/flint.h"]),

    Extension('sage.libs.flint.fmpz_poly',
              sources = ["sage/libs/flint/fmpz_poly.pyx"],
              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],
              include_dirs = [SAGE_ROOT+'/local/include/FLINT/'],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              depends = [SAGE_ROOT + "/local/include/FLINT/flint.h"]),

    Extension('sage.libs.fplll.fplll',
              sources = ['sage/libs/fplll/fplll.pyx'],
              libraries = ['gmp', 'mpfr', 'stdc++', 'fplll'],
              language="c++",
              include_dirs = [SAGE_ROOT +'/local/include/fplll']),

    Extension('sage.libs.linbox.linbox',
              sources = ['sage/libs/linbox/linbox.pyx'],
              # For this to work on cygwin, linboxwrap *must* be
              # before ntl.
              libraries = ['linboxsage', 'ntl', 'linbox', 'gmp', 'gmpxx',
                           'stdc++', 'givaro', BLAS, BLAS2],
              language = 'c++'),

    Extension('sage.libs.libecm',
              sources = ['sage/libs/libecm.pyx'],
              libraries = ['ecm', 'gmp'],
              depends = [SAGE_ROOT + "/local/include/ecm.h"]),

    Extension('sage.libs.mwrank.mwrank',
              sources = ["sage/libs/mwrank/mwrank.pyx",
                         "sage/libs/mwrank/wrap.cc"],
              define_macros = [("NTL_ALL",None)],
              depends = ["sage/libs/mwrank/wrap.h"],
              libraries = ["curvesntl", "g0nntl", "jcntl", "rankntl",
                           "ntl", "gmp", "gmpxx", "stdc++", "m", "pari"]),

    Extension('sage.libs.pari.gen',
              sources = ["sage/libs/pari/gen.pyx"],
              libraries = ['pari', 'gmp']),

    Extension('sage.libs.singular.singular',
              sources = ['sage/libs/singular/singular.pyx'],
              libraries = ['m', 'readline', 'singular', 'singfac',
                           'singcf', 'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs = [SAGE_ROOT +'/local/include/singular'],
              depends = [SAGE_ROOT + "/local/include/libsingular.h"]),

    Extension('sage.libs.symmetrica.symmetrica',
              sources = ["sage/libs/symmetrica/%s"%s for s in ["symmetrica.pyx"]],
              include_dirs = ['/usr/include/malloc/'],
              libraries = ["symmetrica"]),

        ###################################
        ##
        ## sage.libs.cremona
        ##
        ###################################

    Extension('sage.libs.cremona.homspace',
              sources = ["sage/libs/cremona/homspace.pyx"],
              libraries = ['g0nntl', 'jcntl', 'gmpxx', 'ntl', 'gmp',
                           'm', 'stdc++', 'pari', 'curvesntl'],
              language='c++',
              define_macros = [("NTL_ALL",None)]),

    Extension('sage.libs.cremona.mat',
              sources = ["sage/libs/cremona/mat.pyx"],
              libraries = ['g0nntl', 'jcntl', 'gmpxx', 'ntl',
                           'gmp', 'm', 'stdc++', ],
              language='c++',
              define_macros = [("NTL_ALL",None)]),

    Extension('sage.libs.cremona.newforms',
              sources = ["sage/libs/cremona/newforms.pyx"],
              libraries = ['g0nntl', 'jcntl', 'gmpxx', 'ntl', 'gmp',
                           'm', 'stdc++', 'pari', 'curvesntl'],
              language='c++',
              define_macros = [("NTL_ALL",None)]),

        ###################################
        ##
        ## sage.libs.ntl
        ##
        ###################################

    # NOTE: It is *very* important (for cygwin) that csage be the
    # first library listed for all ntl extensions below.

    Extension('sage.libs.ntl.ntl_GF2',
              sources = ["sage/libs/ntl/ntl_GF2.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2E',
              sources = ["sage/libs/ntl/ntl_GF2E.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2EContext',
              sources = ["sage/libs/ntl/ntl_GF2EContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2EX',
              sources = ["sage/libs/ntl/ntl_GF2EX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2X',
              sources = ["sage/libs/ntl/ntl_GF2X.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_p',
              sources = ["sage/libs/ntl/ntl_lzz_p.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_pContext',
              sources = ["sage/libs/ntl/ntl_lzz_pContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_pX',
              sources = ["sage/libs/ntl/ntl_lzz_pX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_mat_GF2',
              sources = ["sage/libs/ntl/ntl_mat_GF2.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_mat_GF2E',
              sources = ["sage/libs/ntl/ntl_mat_GF2E.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_mat_ZZ',
              sources = ["sage/libs/ntl/ntl_mat_ZZ.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ',
              sources = ["sage/libs/ntl/ntl_ZZ.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZX',
              sources = ["sage/libs/ntl/ntl_ZZX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_p',
              sources = ["sage/libs/ntl/ntl_ZZ_p.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pContext',
              sources = ["sage/libs/ntl/ntl_ZZ_pContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pE',
              sources = ["sage/libs/ntl/ntl_ZZ_pE.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pEContext',
              sources = ["sage/libs/ntl/ntl_ZZ_pEContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pEX',
              sources = ["sage/libs/ntl/ntl_ZZ_pEX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pX',
              sources = ["sage/libs/ntl/ntl_ZZ_pX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),



    ################################
    ##
    ## sage.matrix
    ##
    ################################

    Extension('sage.matrix.action',
              sources = ['sage/matrix/action.pyx']),

    Extension('sage.matrix.change_ring',
              sources = ['sage/matrix/change_ring.pyx'],
              libraries=[BLAS, BLAS2, 'gmp'],
              include_dirs = numpy_include_dirs),

    Extension('sage.matrix.matrix',
              sources = ['sage/matrix/matrix.pyx']),

    Extension('sage.matrix.matrix0',
              sources = ['sage/matrix/matrix0.pyx']),

    Extension('sage.matrix.matrix1',
              sources = ['sage/matrix/matrix1.pyx']),

    Extension('sage.matrix.matrix2',
              sources = ['sage/matrix/matrix2.pyx']),

    Extension('sage.matrix.matrix_complex_double_dense',
              sources = ['sage/matrix/matrix_complex_double_dense.pyx'],
              libraries=[BLAS, BLAS2],
              include_dirs = numpy_include_dirs),

    Extension('sage.matrix.matrix_cyclo_dense',
              sources = ['sage/matrix/matrix_cyclo_dense.pyx'],
              language = "c++",
              libraries=['ntl', 'gmp']),

    #Extension('sage.matrix.matrix_cyclo_sparse',
    #          sources = ['sage/matrix/matrix_cyclo_sparse.pyx']),

    Extension('sage.matrix.matrix_dense',
              sources = ['sage/matrix/matrix_dense.pyx']),

    #Extension('sage.matrix.matrix_domain_dense',
    #          sources = ['sage/matrix/matrix_domain_dense.pyx']),

    #Extension('sage.matrix.matrix_domain_sparse',
    #          sources = ['sage/matrix/matrix_domain_sparse.pyx']),

    Extension('sage.matrix.matrix_double_dense',
              sources = ['sage/matrix/matrix_double_dense.pyx'],
              libraries=[BLAS, BLAS2],
              include_dirs = numpy_include_dirs),

    Extension('sage.matrix.matrix_generic_dense',
              sources = ['sage/matrix/matrix_generic_dense.pyx']),

    Extension('sage.matrix.matrix_generic_sparse',
              sources = ['sage/matrix/matrix_generic_sparse.pyx']),

    Extension('sage.matrix.matrix_integer_2x2',
              sources = ['sage/matrix/matrix_integer_2x2.pyx'],
              libraries = ['gmp']),

    # TODO -- change to use BLAS at some point.
    Extension('sage.matrix.matrix_integer_dense',
              sources = ['sage/matrix/matrix_integer_dense.pyx'],
              # order matters for cygwin!!
              libraries = ['iml', 'gmp', 'm', BLAS, BLAS2]),

    Extension('sage.matrix.matrix_integer_sparse',
              sources = ['sage/matrix/matrix_integer_sparse.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_mod2_dense',
              sources = ['sage/matrix/matrix_mod2_dense.pyx'],
              libraries = ['gmp','m4ri', 'png12', 'gd'],
              depends = [SAGE_ROOT + "/local/include/png.h"]),

    Extension('sage.matrix.matrix_modn_dense',
              sources = ['sage/matrix/matrix_modn_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_modn_sparse',
              sources = ['sage/matrix/matrix_modn_sparse.pyx']),

    Extension('sage.matrix.matrix_mpolynomial_dense',
              sources = ['sage/matrix/matrix_mpolynomial_dense.pyx'],
              libraries = ['m', 'readline', 'singular', 'singcf', 'singfac',
                           'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs = [SAGE_ROOT +'/local/include/singular'],
              depends = [SAGE_ROOT + "/local/include/libsingular.h"]),

    #Extension('sage.matrix.matrix_pid_dense',
    #          sources = ['sage/matrix/matrix_pid_dense.pyx']),

    #Extension('sage.matrix.matrix_pid_sparse',
    #          sources = ['sage/matrix/matrix_pid_sparse.pyx']),

    Extension('sage.matrix.matrix_rational_dense',
              sources = ['sage/matrix/matrix_rational_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_rational_sparse',
              sources = ['sage/matrix/matrix_rational_sparse.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_real_double_dense',
              sources = ['sage/matrix/matrix_real_double_dense.pyx'],
              libraries=[BLAS, BLAS2],
              include_dirs = numpy_include_dirs),

    Extension('sage.matrix.matrix_sparse',
              sources = ['sage/matrix/matrix_sparse.pyx']),

    Extension('sage.matrix.matrix_symbolic_dense',
              sources = ['sage/matrix/matrix_symbolic_dense.pyx']),

    Extension('sage.matrix.matrix_window',
              sources = ['sage/matrix/matrix_window.pyx']),

    Extension('sage.matrix.matrix_window_modn_dense',
              sources = ['sage/matrix/matrix_window_modn_dense.pyx']),

    Extension('sage.matrix.misc',
              sources = ['sage/matrix/misc.pyx'],
              libraries=['mpfr','gmp']),

    Extension('sage.matrix.strassen',
              sources = ['sage/matrix/strassen.pyx']),

    #Extension('sage.matrix.padics.matrix_padic_capped_relative_dense',
    #          sources = ['sage/matrix/padics/matrix_padic_capped_relative_dense.pyx']),

    ################################
    ##
    ## sage.media
    ##
    ################################

    Extension('sage.media.channels',
              sources = ['sage/media/channels.pyx']),

    ################################
    ##
    ## sage.misc
    ##
    ################################

    Extension('sage.misc.allocator',
              sources = ['sage/misc/allocator.pyx']),

    Extension('sage.misc.citation',
              sources = ['sage/misc/citation.pyx']),

    Extension('sage.misc.cython_c',
              sources = ['sage/misc/cython_c.pyx']),

    Extension('sage.misc.derivative',
              sources = ['sage/misc/derivative.pyx']),

    Extension('sage.misc.fpickle',
              sources = ['sage/misc/fpickle.pyx']),

    Extension('sage.misc.misc_c',
              sources = ['sage/misc/misc_c.pyx']),

    Extension('sage.misc.parser',
              sources = ['sage/misc/parser.pyx']),

    Extension('sage.misc.pickle_old',
              sources = ['sage/misc/pickle_old.pyx']),

    Extension('sage.misc.randstate',
              sources = ['sage/misc/randstate.pyx']),

    Extension('sage.misc.refcount',
              sources = ['sage/misc/refcount.pyx']),

    Extension('sage.misc.reset',
              sources = ['sage/misc/reset.pyx']),

    Extension('sage.misc.sage_timeit_class',
              sources = ['sage/misc/sage_timeit_class.pyx']),

    Extension('sage.misc.sagex_ds',
              sources = ['sage/misc/sagex_ds.pyx']),

    Extension('sage.misc.search',
              sources = ['sage/misc/search.pyx']),

    Extension('sage.misc.session',
              sources = ['sage/misc/session.pyx']),


    ################################
    ##
    ## sage.modular
    ##
    ################################

    Extension('sage.modular.congroup_pyx',
              sources = ['sage/modular/congroup_pyx.pyx']),

    Extension('sage.modular.modsym.apply',
              sources = ['sage/modular/modsym/apply.pyx'],
              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],
              include_dirs = [SAGE_ROOT+'/local/include/FLINT/'],
              extra_compile_args=["-std=c99",  "-D_XPG6"],
              depends = [SAGE_ROOT + "/local/include/FLINT/flint.h"]),

    Extension('sage.modular.modsym.heilbronn',
              sources = ['sage/modular/modsym/heilbronn.pyx'],
              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],
              include_dirs = [SAGE_ROOT+'/local/include/FLINT/'],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              depends = [SAGE_ROOT + "/local/include/FLINT/flint.h"]),

    Extension('sage.modular.modsym.p1list',
              sources = ['sage/modular/modsym/p1list.pyx'],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.modules
    ##
    ################################

    Extension('sage.modules.free_module_element',
              sources = ['sage/modules/free_module_element.pyx']),

    Extension('sage.modules.module',
              sources = ['sage/modules/module.pyx']),

    Extension('sage.modules.vector_complex_double_dense',
              ['sage/modules/vector_complex_double_dense.pyx'],
              libraries = [BLAS, BLAS2],
              include_dirs = numpy_include_dirs),

    Extension('sage.modules.vector_double_dense',
              ['sage/modules/vector_double_dense.pyx'],
              libraries = [BLAS, BLAS2],
              include_dirs = numpy_include_dirs),

    Extension('sage.modules.vector_integer_dense',
              sources = ['sage/modules/vector_integer_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.modules.vector_modn_dense',
              sources = ['sage/modules/vector_modn_dense.pyx']),

    Extension('sage.modules.vector_rational_dense',
              sources = ['sage/modules/vector_rational_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.modules.vector_real_double_dense',
              ['sage/modules/vector_real_double_dense.pyx'],
              libraries = [BLAS, BLAS2],
              include_dirs = numpy_include_dirs),

    # Extension('sage.modules.vector_rational_sparse',
    #           sources = ['sage/modules/vector_rational_sparse.pyx'],
    #           libraries = ['gmp']),


    ################################
    ##
    ## sage.plot
    ##
    ################################

    Extension('sage.plot.plot3d.base',
              sources = ['sage/plot/plot3d/base.pyx']),

    Extension('sage.plot.plot3d.index_face_set',
              sources = ['sage/plot/plot3d/index_face_set.pyx']),

    Extension('sage.plot.plot3d.parametric_surface',
              sources = ['sage/plot/plot3d/parametric_surface.pyx']),

    Extension('sage.plot.plot3d.shapes',
              sources = ['sage/plot/plot3d/shapes.pyx']),

    Extension('sage.plot.plot3d.transform',
              sources = ['sage/plot/plot3d/transform.pyx']),

    ################################
    ##
    ## sage.rings
    ##
    ################################

    Extension('sage.rings.bernmm',
              sources = ['sage/rings/bernmm.pyx',
                         'sage/rings/bernmm/bern_modp.cpp',
                         'sage/rings/bernmm/bern_modp_util.cpp',
                         'sage/rings/bernmm/bern_rat.cpp'],
              libraries = ['gmp', 'ntl', 'stdc++', 'pthread'],
              depends = ['sage/rings/bernmm/bern_modp.h',
                         'sage/rings/bernmm/bern_modp_util.h',
                         'sage/rings/bernmm/bern_rat.h'],
              language = 'c++',
              define_macros=[('USE_THREADS', '1')]),

    Extension('sage.rings.bernoulli_mod_p',
              sources = ['sage/rings/bernoulli_mod_p.pyx'],
              libraries=['ntl','stdc++'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.complex_double',
              sources = ['sage/rings/complex_double.pyx'],
              libraries = ['gsl', BLAS, BLAS2, 'pari', 'gmp']),

    Extension('sage.rings.complex_interval',
              sources = ['sage/rings/complex_interval.pyx'],
              libraries = ['mpfi', 'mpfr', 'gmp']),

    Extension('sage.rings.complex_number',
              sources = ['sage/rings/complex_number.pyx'],
              libraries = ['mpfr', 'gmp']),

    Extension('sage.rings.integer',
              sources = ['sage/rings/integer.pyx'],
              libraries=['ntl', 'gmp', 'pari']),

    Extension('sage.rings.integer_ring',
              sources = ['sage/rings/integer_ring.pyx'],
              libraries=['ntl', 'gmp']),

    Extension('sage.rings.fast_arith',
              sources = ['sage/rings/fast_arith.pyx'],
              libraries=['gmp','pari','csage']),

    Extension('sage.rings.finite_field_givaro',
              sources = ["sage/rings/finite_field_givaro.pyx"],
              # this order is needed to compile under windows.
              libraries = ['givaro', 'gmpxx', 'gmp', 'm', 'stdc++', ],
              language='c++'),

    Extension('sage.rings.finite_field_ntl_gf2e',
              sources = ['sage/rings/finite_field_ntl_gf2e.pyx'],
              libraries = ['ntl', 'gmp'],
              language = 'c++'),

    Extension('sage.rings.fraction_field_element',
              sources = ['sage/rings/fraction_field_element.pyx']),

    Extension('sage.rings.integer_mod',
              sources = ['sage/rings/integer_mod.pyx'],
              libraries = ['gmp']),

    Extension('sage.rings.laurent_series_ring_element',
              sources = ['sage/rings/laurent_series_ring_element.pyx']),

    Extension('sage.rings.memory',
              sources = ['sage/rings/memory.pyx'],
              libraries=['gmp','stdc++']),

    Extension('sage.rings.morphism',
              sources = ['sage/rings/morphism.pyx']),

    Extension('sage.rings.power_series_mpoly',
              sources = ['sage/rings/power_series_mpoly.pyx']),

    Extension('sage.rings.power_series_poly',
              sources = ['sage/rings/power_series_poly.pyx']),

    Extension('sage.rings.power_series_ring_element',
              sources = ['sage/rings/power_series_ring_element.pyx']),

    Extension('sage.rings.rational',
              sources = ['sage/rings/rational.pyx'],
              libraries=['ntl', 'gmp']),

    Extension('sage.rings.real_double',
              sources = ['sage/rings/real_double.pyx'],
              libraries = ['gsl', 'gmp', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.rings.real_lazy',
              sources = ['sage/rings/real_lazy.pyx']),

    Extension('sage.rings.real_mpfi',
              sources = ['sage/rings/real_mpfi.pyx'],
              libraries = ['mpfi', 'mpfr', 'gmp']),

    Extension('sage.rings.real_mpfr',
              sources = ['sage/rings/real_mpfr.pyx'],
              libraries = ['mpfr', 'pari', 'gmp']),

    Extension('sage.rings.real_rqdf',
              sources = ["sage/rings/real_rqdf.pyx"],
              libraries = ['qd', 'm', 'stdc++','gmp','mpfr' ],
              language='c++'),

    Extension('sage.rings.residue_field',
              sources = ['sage/rings/residue_field.pyx']),

    Extension('sage.rings.ring',
              sources = ['sage/rings/ring.pyx']),

        ################################
        ##
        ## sage.rings.number_field
        ##
        ################################

    Extension('sage.rings.number_field.number_field_base',
              sources = ['sage/rings/number_field/number_field_base.pyx']),

    Extension('sage.rings.number_field.number_field_element',
              sources = ['sage/rings/number_field/number_field_element.pyx'],
              libraries=['ntl','gmp'],
              language = 'c++'),

    Extension('sage.rings.number_field.number_field_element_quadratic',
              sources = ['sage/rings/number_field/number_field_element_quadratic.pyx'],
              libraries=['gmp'],
              language = 'c++'),

    Extension('sage.rings.number_field.number_field_morphisms',
              sources = ['sage/rings/number_field/number_field_morphisms.pyx']),

    Extension('sage.rings.number_field.totallyreal',
              sources = ['sage/rings/number_field/totallyreal.pyx'],
              libraries = ['gmp', 'pari']),

    Extension('sage.rings.number_field.totallyreal_data',
              sources = ['sage/rings/number_field/totallyreal_data.pyx'],
              libraries = ['gmp']),

        ################################
        ##
        ## sage.rings.padics
        ##
        ################################

    Extension('sage.rings.padics.local_generic_element',
              sources = ['sage/rings/padics/local_generic_element.pyx']),

    #Extension('sage.rings.padics.morphism',
    #          sources = ['sage/rings/padics/morphism.pyx'],
    #          libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
    #          language='c++'),

    Extension('sage.rings.padics.padic_base_generic_element',
              sources = ['sage/rings/padics/padic_base_generic_element.pyx']),

    Extension('sage.rings.padics.padic_capped_absolute_element',
              sources = ['sage/rings/padics/padic_capped_absolute_element.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.padics.padic_capped_relative_element',
              sources = ['sage/rings/padics/padic_capped_relative_element.pyx'],
              libraries=['gmp', 'csage']),

    Extension('sage.rings.padics.padic_ext_element',
              sources = ['sage/rings/padics/padic_ext_element.pyx'],
              libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_fixed_mod_element',
              sources = ['sage/rings/padics/padic_fixed_mod_element.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.padics.padic_generic_element',
              sources = ['sage/rings/padics/padic_generic_element.pyx']),

    Extension('sage.rings.padics.padic_printing',
              sources = ['sage/rings/padics/padic_printing.pyx'],
              libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_CA_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CA_element.pyx'],
              libraries = ['gmp','ntl','csage','gmpxx','m','stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_CR_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CR_element.pyx'],
              libraries=['gmp','ntl','csage','gmpxx','m','stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_element.pyx'],
              libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_FM_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_FM_element.pyx'],
              libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
              language='c++'),

    Extension('sage.rings.padics.pow_computer',
              sources = ['sage/rings/padics/pow_computer.pyx'],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.rings.padics.pow_computer_ext',
              sources = ['sage/rings/padics/pow_computer_ext.pyx'],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.rings.padics.rigid_functions',
              sources = ['sage/rings/padics/rigid_functions.pyx']),

        ################################
        ##
        ## sage.rings.polynomial
        ##
        ################################

    Extension('sage.rings.polynomial.cyclotomic',
              sources = ['sage/rings/polynomial/cyclotomic.pyx']),

    Extension('sage.rings.polynomial.laurent_polynomial',
              sources = ['sage/rings/polynomial/laurent_polynomial.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial',
              sources = ['sage/rings/polynomial/multi_polynomial.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial_ideal_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_ideal_libsingular.pyx'],
              libraries = ['m', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs = [SAGE_ROOT +'/local/include/singular'],
              depends = [SAGE_ROOT + "/local/include/libsingular.h"]),

    Extension('sage.rings.polynomial.multi_polynomial_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_libsingular.pyx'],
              libraries = ['m', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs = [SAGE_ROOT +'/local/include/singular'],
              depends = [SAGE_ROOT + "/local/include/libsingular.h"]),

    Extension('sage.rings.polynomial.multi_polynomial_ring_generic',
              sources = ['sage/rings/polynomial/multi_polynomial_ring_generic.pyx']),

    Extension('sage.rings.polynomial.polydict',
              sources = ['sage/rings/polynomial/polydict.pyx']),

    Extension('sage.rings.polynomial.polynomial_compiled',
               sources = ['sage/rings/polynomial/polynomial_compiled.pyx']),

    Extension('sage.rings.polynomial.polynomial_element',
              sources = ['sage/rings/polynomial/polynomial_element.pyx']),

    Extension('sage.rings.polynomial.polynomial_gf2x',
              sources = ['sage/rings/polynomial/polynomial_gf2x.pyx'],
              libraries = ['ntl', 'stdc++', 'gmp'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.polynomial.polynomial_zmod_flint',
              sources = ['sage/rings/polynomial/polynomial_zmod_flint.pyx'],
              libraries = ["csage", "flint", "gmp", "gmpxx", "ntl", "zn_poly"],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              include_dirs = [SAGE_ROOT+'/local/include/FLINT/'],
              depends = [SAGE_ROOT + "/local/include/FLINT/flint.h"]),

    Extension('sage.rings.polynomial.polynomial_integer_dense_flint',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_flint.pyx'],
              language = 'c++',
              libraries = ["csage", "flint", "gmp", "gmpxx", "ntl"],
              include_dirs = [SAGE_ROOT+'/local/include/FLINT/'],
              depends = [SAGE_ROOT + "/local/include/FLINT/flint.h"]),

    Extension('sage.rings.polynomial.polynomial_integer_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_ntl.pyx'],
              libraries = ['ntl', 'stdc++', 'gmp'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.polynomial.polynomial_modn_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_modn_dense_ntl.pyx'],
              libraries = ['ntl', 'stdc++', 'gmp'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.polynomial.pbori',
              sources = ['sage/rings/polynomial/pbori.pyx'],
              libraries=['polybori','pboriCudd','groebner'],
              include_dirs = [SAGE_ROOT+'/local/include/cudd',
                              SAGE_ROOT+'/local/include/polybori',
                              SAGE_ROOT+'/local/include/polybori/groebner'],
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_real_mpfr_dense',
              sources = ['sage/rings/polynomial/polynomial_real_mpfr_dense.pyx'],
              libraries = ['mpfr', 'gmp']),

    Extension('sage.rings.polynomial.real_roots',
              sources = ['sage/rings/polynomial/real_roots.pyx'],
              libraries=['mpfr', 'qd'],
              include_dirs = numpy_include_dirs),

    ################################
    ##
    ## sage.schemes
    ##
    ################################

    Extension('sage.schemes.hyperelliptic_curves.hypellfrob',
              sources = ['sage/schemes/hyperelliptic_curves/hypellfrob.pyx',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.cpp',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.cpp',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.cpp'],
              libraries = ['ntl', 'stdc++', 'gmp', 'zn_poly'],
              depends = ['sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.h',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.h',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.h'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/',
                              'sage/schemes/hyperelliptic_curves/hypellfrob/']),


    ################################
    ##
    ## sage.stats
    ##
    ################################

    Extension('sage.stats.hmm.chmm',
              sources = ['sage/stats/hmm/chmm.pyx'],
              libraries = ['ghmm'],
              include_dirs=numpy_include_dirs),

    Extension('sage.stats.hmm.hmm',
              sources = ['sage/stats/hmm/hmm.pyx'],
              libraries = ['ghmm'],
              include_dirs=numpy_include_dirs),

    ################################
    ##
    ## sage.structure
    ##
    ################################

    Extension('sage.structure.category_object',
              sources = ['sage/structure/category_object.pyx']),

    Extension('sage.structure.coerce',
              sources = ['sage/structure/coerce.pyx']),

    Extension('sage.structure.coerce_actions',
              sources = ['sage/structure/coerce_actions.pyx']),

    Extension('sage.structure.coerce_dict',
              sources = ['sage/structure/coerce_dict.pyx']),

    Extension('sage.structure.coerce_maps',
              sources = ['sage/structure/coerce_maps.pyx']),

    Extension('sage.structure.element',
              sources = ['sage/structure/element.pyx']),

    Extension('sage.structure.factory',
              sources = ['sage/structure/factory.pyx']),

    Extension('sage.structure.generators',
              sources = ['sage/structure/generators.pyx']),

    Extension('sage.structure.mutability',
              sources = ['sage/structure/mutability.pyx']),

    Extension('sage.structure.parent',
              sources = ['sage/structure/parent.pyx']),

    Extension('sage.structure.parent_base',
              sources = ['sage/structure/parent_base.pyx']),

    Extension('sage.structure.parent_gens',
              sources = ['sage/structure/parent_gens.pyx']),

    Extension('sage.structure.parent_old',
              sources = ['sage/structure/parent_old.pyx']),

    Extension('sage.structure.sage_object',
              sources = ['sage/structure/sage_object.pyx']),

    Extension('sage.structure.wrapper_parent',
              sources = ['sage/structure/wrapper_parent.pyx']),


    ################################
    ##
    ## sage.symbolic
    ##
    ################################

    Extension('sage.symbolic.constants',
              sources = ['sage/symbolic/constants.pyx'],
              language = 'c++',
              libraries = ["pynac"]),

    Extension('sage.symbolic.expression',
              sources = ['sage/symbolic/expression.pyx'],
              language = 'c++',
              libraries = ["pynac"]),

    Extension('sage.symbolic.function',
              sources = ['sage/symbolic/function.pyx'],
              language = 'c++',
              libraries = ["pynac"]),

    Extension('sage.symbolic.pynac',
              sources = ['sage/symbolic/pynac.pyx'],
              language = 'c++',
              libraries = ["pynac"]),

    Extension('sage.symbolic.ring',
              sources = ['sage/symbolic/ring.pyx'],
              language = 'c++',
              libraries = ["pynac"]),




    ]

