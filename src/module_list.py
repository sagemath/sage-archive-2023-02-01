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
    debian_include_dirs=["/usr/include","/usr/include/numpy","/usr/include/FLINT","/usr/include/givaro", "/usr/include/gsl","/usr/include/fplll","/usr/include/eclib","/usr/include/gmp++","/usr/include/linbox","/usr/include/NTL","/usr/include/pari","/usr/include/qd","/usr/include/singular","/usr/include/singular/singular","/usr/include/symmetrica","/usr/include/polybori","/usr/include/cudd","/usr/include/polybori/groebner","/usr/include/zn_poly"]
    include_dirs = include_dirs + debian_include_dirs
else:
    debian_include_dirs=[]

#########################################################
### List of modules
#########################################################

ext_modules = [ \

    Extension('sage.structure.sage_object',
              sources = ['sage/structure/sage_object.pyx']),

    Extension('sage.structure.category_object',
              sources = ['sage/structure/category_object.pyx']),

    Extension('sage.structure.parent',
              sources = ['sage/structure/parent.pyx']),

    Extension('sage.structure.parent_old',
              sources = ['sage/structure/parent_old.pyx']),

    Extension('sage.structure.parent_base',
              sources = ['sage/structure/parent_base.pyx']),

    Extension('sage.structure.parent_gens',
              sources = ['sage/structure/parent_gens.pyx']),

    Extension('sage.structure.generators',
              sources = ['sage/structure/generators.pyx']),

    Extension('sage.structure.coerce',
              sources = ['sage/structure/coerce.pyx']),

    Extension('sage.structure.coerce_actions',
              sources = ['sage/structure/coerce_actions.pyx']),

    Extension('sage.structure.coerce_maps',
              sources = ['sage/structure/coerce_maps.pyx']),

    Extension('sage.structure.coerce_dict',
              sources = ['sage/structure/coerce_dict.pyx']),

    Extension('sage.structure.element',
              sources = ['sage/structure/element.pyx']),

    Extension('sage.categories.map',
              sources = ['sage/categories/map.pyx']),

    Extension('sage.categories.morphism',
              sources = ['sage/categories/morphism.pyx']),

    Extension('sage.modules.free_module_element',
              ['sage/modules/free_module_element.pyx']),

    Extension('sage.modules.complex_double_vector',
              ['sage/modules/complex_double_vector.pyx'],
              libraries = ['gsl', BLAS, BLAS2, 'pari', 'gmp'],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.modules.real_double_vector',
              ['sage/modules/real_double_vector.pyx'],
              libraries = ['gsl', BLAS, BLAS2, 'pari','gmp'],
              define_macros = [('GSL_DISABLE_DEPRECATED','1')],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.modules.vector_integer_dense',
              ['sage/modules/vector_integer_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.modules.vector_modn_dense',
              ['sage/modules/vector_modn_dense.pyx']),

    Extension('sage.modules.vector_rational_dense',
              ['sage/modules/vector_rational_dense.pyx'],
              libraries = ['gmp']),

#     Extension('sage.modules.vector_rational_sparse',
#               ['sage/modules/vector_rational_sparse.pyx'],
#               libraries = ['gmp']),

    Extension('sage.libs.pari.gen',
              sources = ["sage/libs/pari/gen.pyx"],
              libraries = ['pari', 'gmp']),

    Extension("sage.libs.mwrank.mwrank",
              sources = ["sage/libs/mwrank/mwrank.pyx",
                         "sage/libs/mwrank/wrap.cc"],
              define_macros = [("NTL_ALL",None)],
              depends = ["sage/libs/mwrank/wrap.h"],
              libraries = ["curvesntl", "g0nntl", "jcntl", "rankntl",
                           "ntl", "gmp", "gmpxx", "stdc++", "m", "pari"]),

    Extension('sage.libs.flint.fmpz_poly',
              sources = ["sage/libs/flint/fmpz_poly.pyx"],
              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/include/FLINT/'],
              extra_compile_args=["-std=c99", "-D_XPG6"]),

    # NOTE: It is *very* important (for cygwin) that csage be the first library
    # listed below for ntl.
    Extension('sage.libs.ntl.ntl_ZZ',
              sources = ["sage/libs/ntl/ntl_ZZ.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZX',
              sources = ["sage/libs/ntl/ntl_ZZX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pContext',
              sources = ["sage/libs/ntl/ntl_ZZ_pContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_p',
              sources = ["sage/libs/ntl/ntl_ZZ_p.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pX',
              sources = ["sage/libs/ntl/ntl_ZZ_pX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pEContext',
              sources = ["sage/libs/ntl/ntl_ZZ_pEContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pE',
              sources = ["sage/libs/ntl/ntl_ZZ_pE.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pEX',
              sources = ["sage/libs/ntl/ntl_ZZ_pEX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_pContext',
              sources = ["sage/libs/ntl/ntl_lzz_pContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_p',
              sources = ["sage/libs/ntl/ntl_lzz_p.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_pX',
              sources = ["sage/libs/ntl/ntl_lzz_pX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2',
              sources = ["sage/libs/ntl/ntl_GF2.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2X',
              sources = ["sage/libs/ntl/ntl_GF2X.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2EContext',
              sources = ["sage/libs/ntl/ntl_GF2EContext.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2E',
              sources = ["sage/libs/ntl/ntl_GF2E.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2EX',
              sources = ["sage/libs/ntl/ntl_GF2EX.pyx"],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_mat_ZZ',
              sources = ["sage/libs/ntl/ntl_mat_ZZ.pyx"],
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



    Extension('sage.matrix.matrix',
              ['sage/matrix/matrix.pyx']),

    Extension('sage.matrix.action',
              ['sage/matrix/action.pyx']),

    Extension('sage.matrix.misc',
              ['sage/matrix/misc.pyx'],
              libraries=['mpfr','gmp']),

    Extension('sage.matrix.matrix_dense',
              ['sage/matrix/matrix_dense.pyx']),

    Extension('sage.matrix.matrix_generic_dense',
              ['sage/matrix/matrix_generic_dense.pyx']),

    Extension('sage.matrix.matrix_sparse',
              ['sage/matrix/matrix_sparse.pyx']),

    Extension('sage.matrix.matrix_generic_sparse',
              ['sage/matrix/matrix_generic_sparse.pyx']),

    #matrix_padic_capped_relative_dense =
    #Extension('sage.matrix.padics.matrix_padic_capped_relative_dense',
    #          ['sage/matrix/padics/matrix_padic_capped_relative_dense.pyx']),

    # TODO -- change to use BLAS at some point.
    Extension('sage.matrix.matrix_integer_dense',
              ['sage/matrix/matrix_integer_dense.pyx'],
              # order matters for cygwin!!
              libraries = ['iml', 'gmp', 'm', BLAS, BLAS2]),

    Extension('sage.matrix.matrix_rational_dense',
              ['sage/matrix/matrix_rational_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_rational_sparse',
              ['sage/matrix/matrix_rational_sparse.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_cyclo_dense',
              ['sage/matrix/matrix_cyclo_dense.pyx'],
              language = "c++",
              libraries=['ntl', 'gmp']),

    #Extension('sage.matrix.matrix_cyclo_sparse',
    #          ['sage/matrix/matrix_cyclo_sparse.pyx']),

    #Extension('sage.matrix.matrix_domain_dense',
    #          ['sage/matrix/matrix_domain_dense.pyx']),

    #Extension('sage.matrix.matrix_domain_sparse',
    #          ['sage/matrix/matrix_domain_sparse.pyx']),

    #Extension('sage.matrix.matrix_pid_dense',
    #          ['sage/matrix/matrix_pid_dense.pyx']),

    #Extension('sage.matrix.matrix_pid_sparse',
    #          ['sage/matrix/matrix_pid_sparse.pyx']),

    Extension('sage.matrix.matrix_integer_2x2',
              ['sage/matrix/matrix_integer_2x2.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_double_dense',
              ['sage/matrix/matrix_double_dense.pyx'],
              libraries=[BLAS, BLAS2],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.matrix.matrix_integer_sparse',
              ['sage/matrix/matrix_integer_sparse.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_real_double_dense',
              ['sage/matrix/matrix_real_double_dense.pyx'],
              libraries=[BLAS, BLAS2],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.matrix.change_ring',
              ['sage/matrix/change_ring.pyx'],
              libraries=[BLAS, BLAS2, 'gmp'],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.matrix.matrix_complex_double_dense',
              ['sage/matrix/matrix_complex_double_dense.pyx'],
              libraries=[BLAS, BLAS2],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.libs.linbox.linbox',
              ['sage/libs/linbox/linbox.pyx'],
              # For this to work on cygwin, linboxwrap
              # *must* be before ntl.
              libraries = ['linboxsage', 'ntl', 'linbox', 'gmp', 'gmpxx',
                           'stdc++', 'givaro', BLAS, BLAS2],
              language = 'c++'),

    Extension('sage.matrix.matrix_modn_dense',
              ['sage/matrix/matrix_modn_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_modn_sparse',
              ['sage/matrix/matrix_modn_sparse.pyx']),

    Extension('sage.matrix.matrix_mod2_dense',
              ['sage/matrix/matrix_mod2_dense.pyx'],
              libraries = ['gmp','m4ri', 'png', 'gd']),

    Extension('sage.matrix.matrix_mpolynomial_dense',
              ['sage/matrix/matrix_mpolynomial_dense.pyx'],
              libraries = ['m', 'readline', 'singular', 'singcf', 'singfac',
                           'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs=debian_include_dirs + [SAGE_ROOT +'/local/include/singular']),

    Extension('sage.matrix.matrix_symbolic_dense',
              ['sage/matrix/matrix_symbolic_dense.pyx']),

    Extension('sage.libs.cremona.mat',
              sources = ["sage/libs/cremona/mat.pyx"],
              libraries = ['g0nntl', 'jcntl', 'gmpxx', 'ntl',
                           'gmp', 'm', 'stdc++', ],
              language='c++',
              define_macros = [("NTL_ALL",None)]),

    Extension('sage.libs.cremona.homspace',
              sources = ["sage/libs/cremona/homspace.pyx"],
              libraries = ['g0nntl', 'jcntl', 'gmpxx', 'ntl', 'gmp',
                           'm', 'stdc++', 'pari', 'curvesntl'],
              language='c++',
              define_macros = [("NTL_ALL",None)]),

    Extension('sage.libs.cremona.newforms',
              sources = ["sage/libs/cremona/newforms.pyx"],
              libraries = ['g0nntl', 'jcntl', 'gmpxx', 'ntl', 'gmp',
                           'm', 'stdc++', 'pari', 'curvesntl'],
              language='c++',
              define_macros = [("NTL_ALL",None)]),

    Extension('sage.rings.finite_field_givaro',
              sources = ["sage/rings/finite_field_givaro.pyx"],
              # this order is needed to compile under windows.
              libraries = ['givaro', 'gmpxx', 'gmp', 'm', 'stdc++', ],
              language='c++'),

    Extension('sage.rings.finite_field_ntl_gf2e',
              sources = ['sage/rings/finite_field_ntl_gf2e.pyx'],
              libraries = ['ntl', 'gmp'],
              language = 'c++'),

    Extension('sage.libs.singular.singular',
              sources = ['sage/libs/singular/singular.pyx'],
              libraries = ['m', 'readline', 'singular', 'singfac',
                           'singcf', 'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs=debian_include_dirs + [SAGE_ROOT +'/local/include/singular']),

    Extension('sage.libs.fplll.fplll',
              sources = ['sage/libs/fplll/fplll.pyx'],
              libraries = ['gmp', 'mpfr', 'stdc++', 'fplll'],
              language="c++",
              include_dirs=debian_include_dirs + [SAGE_ROOT +'/local/include/fplll']),

    Extension('sage.gsl.dwt',
              ['sage/gsl/dwt.pyx'],
              libraries=['gsl',BLAS],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    ################ GSL wrapping ######################
    Extension('sage.gsl.gsl_array',
              ['sage/gsl/gsl_array.pyx'],
              libraries=['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.ode',
              ['sage/gsl/ode.pyx'],
              libraries=['gsl',BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.fft',
              ['sage/gsl/fft.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.interpolation',
              ['sage/gsl/interpolation.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.callback',
              ['sage/gsl/callback.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.probability_distribution',
              ['sage/gsl/probability_distribution.pyx'],
              libraries=['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.integration',
              ['sage/gsl/integration.pyx'],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')],
              libraries=['gsl',BLAS, BLAS2]),

    Extension('sage.rings.real_double',
              ['sage/rings/real_double.pyx'],
              libraries = ['gsl', 'gmp', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.rings.complex_double',
              ['sage/rings/complex_double.pyx'],
              libraries = ['gsl', BLAS, BLAS2, 'pari', 'gmp']),

    Extension('sage.rings.real_rqdf',
              sources = ["sage/rings/real_rqdf.pyx"],
              libraries = ['qd', 'm', 'stdc++','gmp','mpfr' ],
              language='c++'),

    Extension('sage.rings.complex_number',
              ['sage/rings/complex_number.pyx'],
              libraries = ['mpfr', 'gmp']),

    Extension('sage.misc.sagex_ds',
              ['sage/misc/sagex_ds.pyx']),

    Extension('sage.misc.citation',
              ['sage/misc/citation.pyx']),

    Extension('sage.libs.symmetrica.symmetrica',
              sources = ["sage/libs/symmetrica/%s"%s for s in ["symmetrica.pyx"]],
              include_dirs=debian_include_dirs + ['/usr/include/malloc/'],
              libraries = ["symmetrica"]),

    Extension('sage.finance.time_series',
              ['sage/finance/time_series.pyx'],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.finance.markov_multifractal_cython',
              ['sage/finance/markov_multifractal_cython.pyx']),

    Extension('sage.finance.fractal',
              ['sage/finance/fractal.pyx']),


    Extension('sage.stats.hmm.hmm',
              ['sage/stats/hmm/hmm.pyx'],
              libraries = ['ghmm'],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.stats.hmm.chmm',
              ['sage/stats/hmm/chmm.pyx'],
              libraries = ['ghmm'],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']),

    Extension('sage.media.channels',
              sources = ['sage/media/channels.pyx']),

    Extension('sage.ext.sig',
              sources = ['sage/ext/sig.pyx']),

    Extension('sage.ext.fast_eval',
              sources = ['sage/ext/fast_eval.pyx']),

    Extension('sage.rings.fast_arith',
              sources = ['sage/rings/fast_arith.pyx'],
              libraries=['gmp','pari','csage']),

    Extension('sage.ext.multi_modular',
              sources = ['sage/ext/multi_modular.pyx'],
              libraries=['gmp']),

    Extension('sage.modular.congroup_pyx',
              sources = ['sage/modular/congroup_pyx.pyx']),

    Extension('sage.categories.functor',
              sources = ['sage/categories/functor.pyx']),

    Extension('sage.categories.action',
              sources = ['sage/categories/action.pyx']),

    Extension('sage.modules.module',
              sources = ['sage/modules/module.pyx']),

    Extension('sage.rings.ring',
              sources = ['sage/rings/ring.pyx']),

    Extension('sage.rings.polynomial.cyclotomic',
              sources = ['sage/rings/polynomial/cyclotomic.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial',
              sources = ['sage/rings/polynomial/multi_polynomial.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial_ring_generic',
              sources = ['sage/rings/polynomial/multi_polynomial_ring_generic.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_libsingular.pyx'],
              libraries = ['m', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs=debian_include_dirs + [SAGE_ROOT +'/local/include/singular']),

    Extension('sage.rings.polynomial.multi_polynomial_ideal_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_ideal_libsingular.pyx'],
              libraries = ['m', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs=debian_include_dirs + [SAGE_ROOT +'/local/include/singular']),

    Extension('sage.groups.group',
              sources = ['sage/groups/group.pyx']), \

    Extension('sage.groups.perm_gps.permgroup_element',
              sources = ['sage/groups/perm_gps/permgroup_element.pyx']), \

    Extension('sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label',
              sources = ['sage/groups/perm_gps/partn_ref/automorphism_group_canonical_label.pyx']), \

    Extension('sage.groups.perm_gps.partn_ref.double_coset',
              sources = ['sage/groups/perm_gps/partn_ref/double_coset.pyx']), \

    Extension('sage.groups.perm_gps.partn_ref.refinement_graphs',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_graphs.pyx']), \

    Extension('sage.groups.perm_gps.partn_ref.refinement_binary',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_binary.pyx']), \

    Extension('sage.groups.perm_gps.partn_ref.refinement_matrices',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_matrices.pyx']), \

    Extension('sage.structure.sage_object',
              sources = ['sage/structure/sage_object.pyx']), \

    Extension('sage.structure.parent',
              sources = ['sage/structure/parent.pyx']), \

    Extension('sage.structure.parent_base',
              sources = ['sage/structure/parent_base.pyx']), \

    Extension('sage.structure.parent_gens',
              sources = ['sage/structure/parent_gens.pyx']), \

    Extension('sage.ext.interactive_constructors_c',
              sources = ['sage/ext/interactive_constructors_c.pyx']), \

    Extension('sage.misc.derivative',
              sources = ['sage/misc/derivative.pyx']), \

    Extension('sage.misc.cython_c',
              sources = ['sage/misc/cython_c.pyx']), \

    Extension('sage.misc.misc_c',
              sources = ['sage/misc/misc_c.pyx']), \

    Extension('sage.misc.session',
              sources = ['sage/misc/session.pyx']), \

    Extension('sage.misc.parser',
              ['sage/misc/parser.pyx']), \

    Extension('sage.misc.refcount',
              sources = ['sage/misc/refcount.pyx']), \

    Extension('sage.misc.randstate',
              sources = ['sage/misc/randstate.pyx']), \

    Extension('sage.rings.real_mpfr',
              sources = ['sage/rings/real_mpfr.pyx'],
              libraries = ['mpfr', 'pari', 'gmp']), \

    Extension('sage.rings.real_lazy',
              sources = ['sage/rings/real_lazy.pyx']), \

    Extension('sage.rings.real_mpfi',
              sources = ['sage/rings/real_mpfi.pyx'],
              libraries = ['mpfi', 'mpfr', 'gmp']), \

    Extension('sage.rings.complex_interval',
              sources = ['sage/rings/complex_interval.pyx'],
              libraries = ['mpfi', 'mpfr', 'gmp']), \

    Extension('sage.rings.residue_field',
              sources = ['sage/rings/residue_field.pyx']), \

    Extension('sage.rings.integer',
              sources = ['sage/rings/integer.pyx'],
              libraries=['ntl', 'gmp', 'pari']), \

    Extension('sage.misc.allocator',
              sources = ['sage/misc/allocator.pyx']), \

    Extension('sage.rings.integer_ring',
              sources = ['sage/rings/integer_ring.pyx'],
              libraries=['ntl', 'gmp']), \

    Extension('sage.libs.libecm',
              sources = ['sage/libs/libecm.pyx'],
              libraries = ['ecm', 'gmp']), \

    Extension('sage.rings.padics.pow_computer',
              sources = ['sage/rings/padics/pow_computer.pyx'],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'), \

    Extension('sage.rings.padics.pow_computer_ext',
              sources = ['sage/rings/padics/pow_computer_ext.pyx'],
              libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"],
              language='c++'), \

    Extension('sage.rings.padics.local_generic_element',
              sources = ['sage/rings/padics/local_generic_element.pyx']),

    Extension('sage.rings.padics.padic_generic_element',
              sources = ['sage/rings/padics/padic_generic_element.pyx']),

    Extension('sage.rings.padics.padic_base_generic_element',
              sources = ['sage/rings/padics/padic_base_generic_element.pyx']),

    Extension('sage.rings.padics.padic_fixed_mod_element',
              sources = ['sage/rings/padics/padic_fixed_mod_element.pyx'],
              libraries=['gmp']),

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

    Extension('sage.rings.padics.padic_ZZ_pX_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_element.pyx'],
              libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_FM_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_FM_element.pyx'],
              libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_CR_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CR_element.pyx'],
              libraries=['gmp','ntl','csage','gmpxx','m','stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_CA_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CA_element.pyx'],
              libraries = ['gmp','ntl','csage','gmpxx','m','stdc++'],
              language='c++'),

    Extension('sage.rings.padics.padic_printing',
              sources = ['sage/rings/padics/padic_printing.pyx'],
              libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
              language='c++'),

    Extension('sage.rings.padics.rigid_functions',
              sources = ['sage/rings/padics/rigid_functions.pyx']),

    #Extension('sage.rings.padics.morphism',
    #          sources = ['sage/rings/padics/morphism.pyx'],
    #          libraries=['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'],
    #          language='c++'),

    Extension('sage.rings.memory', \
              sources = ['sage/rings/memory.pyx'], \
              libraries=['gmp','stdc++']), \

    Extension('sage.rings.bernoulli_mod_p',
              sources = ['sage/rings/bernoulli_mod_p.pyx'],
              libraries=['ntl','stdc++'],
              language = 'c++',
              include_dirs=debian_include_dirs + ['sage/libs/ntl/']), \

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
              define_macros=[('USE_THREADS', '1')]), \

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
              include_dirs=debian_include_dirs + ['sage/libs/ntl/', 'sage/schemes/hyperelliptic_curves/hypellfrob/']), \

    Extension('sage.rings.polynomial.polynomial_compiled',
               sources = ['sage/rings/polynomial/polynomial_compiled.pyx']), \

    Extension('sage.rings.polynomial.polynomial_element',
              sources = ['sage/rings/polynomial/polynomial_element.pyx']), \

    Extension('sage.rings.polynomial.polynomial_integer_dense_flint',
                 sources = ['sage/rings/polynomial/polynomial_integer_dense_flint.pyx'],
                 language = 'c++',
                 libraries = ["csage", "flint", "gmp", "gmpxx", "ntl"],
                 include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/include/FLINT/']),

    Extension('sage.rings.polynomial.polynomial_integer_dense_ntl',
                 sources = ['sage/rings/polynomial/polynomial_integer_dense_ntl.pyx'],
                 libraries = ['ntl', 'stdc++', 'gmp'],
                 language = 'c++',
                 include_dirs=debian_include_dirs + ['sage/libs/ntl/']), \

    Extension('sage.rings.polynomial.polynomial_modn_dense_ntl',
                 sources = ['sage/rings/polynomial/polynomial_modn_dense_ntl.pyx'],
                 libraries = ['ntl', 'stdc++', 'gmp'],
                 language = 'c++',
                 include_dirs=debian_include_dirs + ['sage/libs/ntl/']), \

    Extension('sage.rings.polynomial.polynomial_gf2x',
                 sources = ['sage/rings/polynomial/polynomial_gf2x.pyx'],
                 libraries = ['ntl', 'stdc++', 'gmp'],
                 language = 'c++',
                 include_dirs=debian_include_dirs + ['sage/libs/ntl/']), \

    Extension('sage.rings.polynomial.polynomial_real_mpfr_dense',
                 sources = ['sage/rings/polynomial/polynomial_real_mpfr_dense.pyx'],
                 libraries = ['mpfr', 'gmp']), \

    Extension('sage.rings.power_series_ring_element',
              sources = ['sage/rings/power_series_ring_element.pyx']), \

    Extension('sage.rings.power_series_poly',
              sources = ['sage/rings/power_series_poly.pyx']), \

    Extension('sage.rings.power_series_mpoly',
              sources = ['sage/rings/power_series_mpoly.pyx']), \

    Extension('sage.rings.laurent_series_ring_element',
              sources = ['sage/rings/laurent_series_ring_element.pyx']), \

    Extension('sage.rings.rational',
              sources = ['sage/rings/rational.pyx'],
              libraries=['ntl', 'gmp']), \

    Extension('sage.rings.sparse_poly',
              sources = ['sage/rings/sparse_poly.pyx'],
              libraries=['gmp']), \

    Extension('sage.rings.polynomial.polydict',
              sources = ['sage/rings/polynomial/polydict.pyx']), \

    Extension('sage.rings.polynomial.real_roots',
              sources = ['sage/rings/polynomial/real_roots.pyx'],
              libraries=['mpfr', 'qd']), \

    Extension('sage.rings.polynomial.laurent_polynomial',
              sources = ['sage/rings/polynomial/laurent_polynomial.pyx']), \


    Extension('sage.rings.number_field.number_field_element',
              sources = ['sage/rings/number_field/number_field_element.pyx'],
              libraries=['ntl','gmp'],
              language = 'c++'), \

    Extension('sage.rings.number_field.number_field_element_quadratic',
              sources = ['sage/rings/number_field/number_field_element_quadratic.pyx'],
              libraries=['gmp'],
              language = 'c++'), \
                    # Needs c++ because it has c++ members in class struct

    Extension('sage.rings.number_field.number_field_base',
              sources = ['sage/rings/number_field/number_field_base.pyx']), \

    Extension('sage.rings.number_field.totallyreal',
              ['sage/rings/number_field/totallyreal.pyx'],
              libraries = ['gmp', 'pari']), \

    Extension('sage.rings.number_field.totallyreal_data',
              ['sage/rings/number_field/totallyreal_data.pyx'],
              libraries = ['gmp'],
              include_dirs = [SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy']), \

    Extension('sage.rings.morphism',
              sources = ['sage/rings/morphism.pyx']), \

    Extension('sage.structure.wrapper_parent',
              sources = ['sage/structure/wrapper_parent.pyx']), \

    Extension('sage.misc.search',
              ['sage/misc/search.pyx']), \

    Extension('sage.misc.reset',
              ['sage/misc/reset.pyx']), \

    Extension('sage.calculus.var',
              ['sage/calculus/var.pyx']), \

    Extension('sage.symbolic.ring',
                 sources = ['sage/symbolic/ring.pyx'],
                 language = 'c++',
                 libraries = ["pynac"]), \

    Extension('sage.symbolic.expression',
                 sources = ['sage/symbolic/expression.pyx'],
                 language = 'c++',
                 libraries = ["pynac"]), \

    Extension('sage.symbolic.pynac',
                 sources = ['sage/symbolic/pynac.pyx'],
                 language = 'c++',
                 libraries = ["pynac"]), \

    Extension('sage.symbolic.constants',
                 sources = ['sage/symbolic/constants.pyx'],
                  language = 'c++',
                  libraries = ["pynac"]), \

    Extension('sage.symbolic.function',
                 sources = ['sage/symbolic/function.pyx'],
                 language = 'c++',
                 libraries = ["pynac"]), \

    Extension('sage.modular.modsym.heilbronn',
              ['sage/modular/modsym/heilbronn.pyx'],
              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/include/FLINT/'],
              extra_compile_args=["-std=c99"]), \

    Extension('sage.modular.modsym.apply',
              ['sage/modular/modsym/apply.pyx'],
              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/include/FLINT/'],
              extra_compile_args=["-std=c99"]), \

    Extension('sage.modular.modsym.p1list',
              ['sage/modular/modsym/p1list.pyx'],
              libraries = ['gmp']), \

    Extension('sage.structure.mutability',
              ['sage/structure/mutability.pyx']
              ), \

    Extension('sage.matrix.matrix0',
              ['sage/matrix/matrix0.pyx']
              ), \

    Extension('sage.matrix.matrix1',
              ['sage/matrix/matrix1.pyx']
              ), \

    Extension('sage.matrix.matrix2',
              ['sage/matrix/matrix2.pyx']
              ), \

    Extension('sage.matrix.matrix',
              ['sage/matrix/matrix.pyx']
              ), \

    Extension('sage.matrix.matrix_window',
              ['sage/matrix/matrix_window.pyx']
              ), \

    Extension('sage.matrix.matrix_window_modn_dense',
              ['sage/matrix/matrix_window_modn_dense.pyx']
              ), \

    Extension('sage.matrix.strassen',
              ['sage/matrix/strassen.pyx']
              ), \

    Extension('sage.rings.integer_mod',
              ['sage/rings/integer_mod.pyx'],
              libraries = ['gmp']
              ), \

    Extension('sage.combinat.expnums',
              ['sage/combinat/expnums.pyx'],
              libraries = ['gmp']
              ), \

    Extension('sage.combinat.partitions',
              ['sage/combinat/partitions.pyx',
               'sage/combinat/partitions_c.cc'],
              libraries = ['qd', 'gmp', 'mpfr'],
              depends = ['sage/combinat/partitions_c.h'],
              language='c++'
              ), \

    Extension('sage.combinat.matrices.dancing_links',
              ['sage/combinat/matrices/dancing_links.pyx'],
              libraries = ["stdc++"],
              language='c++'
              ), \

    Extension('sage.graphs.base.c_graph',
              ['sage/graphs/base/c_graph.pyx']
              ), \

    Extension('sage.graphs.base.sparse_graph',
              ['sage/graphs/base/sparse_graph.pyx']
              ), \

    Extension('sage.graphs.base.dense_graph',
              ['sage/graphs/base/dense_graph.pyx']
              ), \

    Extension('sage.graphs.graph_fast',
              ['sage/graphs/graph_fast.pyx'],
              libraries = ['gmp']
              ), \

    Extension('sage.graphs.graph_isom',
              ['sage/graphs/graph_isom.pyx']
              ), \

    Extension('sage.graphs.planarity',
              ['sage/graphs/planarity.pyx',
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
                         'sage/graphs/planarity/stack.h']
              ), \

    Extension('sage.graphs.chrompoly',
              ['sage/graphs/chrompoly.pyx']
              ), \

    Extension('sage.coding.binary_code',
              ['sage/coding/binary_code.pyx']
              ), \


    Extension('sage.plot.plot3d.base',
              ['sage/plot/plot3d/base.pyx']
              ), \
    Extension('sage.plot.plot3d.transform',
              ['sage/plot/plot3d/transform.pyx']
              ), \
    Extension('sage.plot.plot3d.index_face_set',
              ['sage/plot/plot3d/index_face_set.pyx']
              ), \
    Extension('sage.plot.plot3d.parametric_surface',
              ['sage/plot/plot3d/parametric_surface.pyx']
              ), \

    Extension('sage.plot.plot3d.shapes',
              ['sage/plot/plot3d/shapes.pyx']
              ), \

    Extension('sage.rings.polynomial.pbori',
              sources = ['sage/rings/polynomial/pbori.pyx'],
              libraries=['polybori','pboriCudd','groebner'],
              include_dirs=debian_include_dirs + [SAGE_ROOT+'/local/include/cudd',
                                                  SAGE_ROOT+'/local/include/polybori',
                                                  SAGE_ROOT+'/local/include/polybori/groebner'],
              language = 'c++'),

    Extension('sage.misc.sage_timeit_class',
              ['sage/misc/sage_timeit_class.pyx']),

    Extension('sage.misc.fpickle',
              ['sage/misc/fpickle.pyx']),
    ]

