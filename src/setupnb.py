#!/usr/bin/env python

import distutils.sysconfig, os, sys
from distutils.core import setup, Extension


if len(sys.argv) > 1 and sys.argv[1] == "sdist":
    sdist = True
else:
    sdist = False

if not os.environ.has_key('SAGE_VERSION'):
    SAGE_VERSION=0
else:
    SAGE_VERSION = os.environ['SAGE_VERSION']

SITE_PACKAGES = '%s/lib/python/site-packages/'%SAGE_LOCAL
if not os.path.exists(SITE_PACKAGES):
    SITE_PACKAGES = '%s/lib/python2.5/site-packages/'%SAGE_LOCAL
    if not os.path.exists(SITE_PACKAGES):
        SITE_PACKAGES = '%s/lib/python2.4/site-packages/'%SAGE_LOCAL
        if not os.path.exists(SITE_PACKAGES):
            raise RuntimeError, "Unable to find site-packages directory (see setup.py file in sage python code)."

SITE_PACKAGES_REL=SITE_PACKAGES[len(SAGE_LOCAL)+5:]

if not os.path.exists('build/sage'):
    os.makedirs('build/sage')

sage_link = SITE_PACKAGES + '/sage'
if not os.path.islink(sage_link) or not os.path.exists(sage_link):
    os.system('rm -rf "%s"'%sage_link)
    os.system('cd %s; ln -sf ../../../../devel/sage/build/sage .'%SITE_PACKAGES)


include_dirs = ['%s/include'%SAGE_LOCAL, '%s/include/python'%SAGE_LOCAL, \
                '%s/sage/sage/ext'%SAGE_DEVEL]

#####################################################

hanke = Extension(name = "sage.libs.hanke.hanke",
              sources = ["sage/libs/hanke/hanke.pyx",
                         "sage/libs/hanke/wrap.cc",
                         "sage/libs/hanke/Matrix_mpz/Matrix_mpz.cc",
                         "sage/libs/hanke/Matrix_mpz/CountLocal2.cc",
                         "sage/libs/hanke/Matrix_mpz/CountLocal.cc",
                         "sage/libs/hanke/Matrix_mpz/Local_Constants.cc",
                         "sage/libs/hanke/Matrix_mpz/Local_Density_Front.cc",
                         "sage/libs/hanke/Matrix_mpz/Local_Density_Congruence.cc",
                         "sage/libs/hanke/Matrix_mpz/Local_Normal.cc",
                         "sage/libs/hanke/Matrix_mpz/Local_Invariants.cc",
                         "sage/libs/hanke/Utilities/string_utils.cc",
                         "sage/libs/hanke/GMP_class_extras/mpz_class_extras.cc",
                         "sage/libs/hanke/GMP_class_extras/vectors.cc" ],
                   libraries = ["gmp", "gmpxx", "stdc++"])


# NOTE: It is *very* important (for cygwin) that csage be the first library
# listed below for ntl.
ntl = Extension('sage.libs.ntl.ntl',
                 sources = ["sage/libs/ntl/ntl.pyx"],
                 libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"]
                 )

mwrank =  Extension("sage.libs.mwrank.mwrank",
                    sources = ["sage/libs/mwrank/mwrank.pyx",
                         "sage/libs/mwrank/wrap.cc"],
                    define_macros = [("NTL_ALL",None)],
                    libraries = ["mwrank", "ntl", "gmp", "gmpxx", "stdc++", "m", "pari"])

pari = Extension('sage.libs.pari.gen',
                 sources = ["sage/libs/pari/gen.pyx"],
                 libraries = ['pari', 'gmp'])


givaro_gfq = Extension('sage.rings.finite_field_givaro',
                       sources = ["sage/rings/finite_field_givaro.pyx"],
                       libraries = ['givaro', 'gmpxx', 'gmp', 'm', 'stdc++', ],   # this order is needed to compile under windows.
                       language='c++'
                       )


qd = Extension('sage.rings.real_rqdf',
                       sources = ["sage/rings/real_rqdf.pyx"],
                       libraries = ['qd', 'm', 'stdc++','gmp','mpfr' ],
                       language='c++'
                       )

matrix = Extension('sage.matrix.matrix', ['sage/matrix/matrix.pyx'])

matrix_misc = Extension('sage.matrix.misc', ['sage/matrix/misc.pyx'],
                        libraries=['gmp'])

matrix_dense = Extension('sage.matrix.matrix_dense',
                         ['sage/matrix/matrix_dense.pyx'])

matrix_sparse = Extension('sage.matrix.matrix_sparse',
                          ['sage/matrix/matrix_sparse.pyx'])

matrix_generic_dense = Extension('sage.matrix.matrix_generic_dense',
                                 ['sage/matrix/matrix_generic_dense.pyx'])

matrix_generic_sparse = Extension('sage.matrix.matrix_generic_sparse',
                                  ['sage/matrix/matrix_generic_sparse.pyx'])

matrix_domain_dense = Extension('sage.matrix.matrix_domain_dense',
                                ['sage/matrix/matrix_domain_dense.pyx'])

matrix_domain_sparse = Extension('sage.matrix.matrix_domain_sparse',
              ['sage/matrix/matrix_domain_sparse.pyx'])

matrix_pid_dense = Extension('sage.matrix.matrix_pid_dense',
                       ['sage/matrix/matrix_pid_dense.pyx'])

matrix_pid_sparse = Extension('sage.matrix.matrix_pid_sparse',
                       ['sage/matrix/matrix_pid_sparse.pyx'])


matrix_integer_2x2 = Extension('sage.matrix.matrix_integer_2x2',
                                 ['sage/matrix/matrix_integer_2x2.pyx'],
                                 libraries = ['gmp'])

linbox = Extension('sage.libs.linbox.linbox',
                   ['sage/libs/linbox/linbox.pyx'],
                   # For this to work on cygwin, linboxwrap *must* be before ntl.
                   libraries = ['linboxsage', 'ntl', 'linbox', 'gmp', 'gmpxx', 'stdc++', 'givaro', BLAS],
                   language = 'c++')

libsingular = Extension('sage.libs.singular.singular',
                        sources = ['sage/libs/singular/singular.pyx'],
                        libraries = ['gmp', 'm', 'readline', 'singular', 'singfac', 'singcf', 'omalloc', 'givaro', 'gmpxx'],
                        language="c++",
                        include_dirs=[SAGE_ROOT +'/local/include/singular']
                        )


matrix_modn_dense = Extension('sage.matrix.matrix_modn_dense',
                              ['sage/matrix/matrix_modn_dense.pyx'],
                              libraries = ['gmp'])

matrix_mod2_dense = Extension('sage.matrix.matrix_mod2_dense',
                              ['sage/matrix/matrix_mod2_dense.pyx',
                               'sage/libs/m4ri/packedmatrix.c',
                               'sage/libs/m4ri/matrix.c',
                               'sage/libs/m4ri/brilliantrussian.c',
                               'sage/libs/m4ri/grayflex.c',],
                              libraries = ['gmp'])

matrix_modn_sparse = Extension('sage.matrix.matrix_modn_sparse',
                               ['sage/matrix/matrix_modn_sparse.pyx'])

matrix_field_dense = Extension('sage.matrix.matrix_field_dense',
                       ['sage/matrix/matrix_field_dense.pyx'])

matrix_field_sparse = Extension('sage.matrix.matrix_field_sparse',
                       ['sage/matrix/matrix_field_sparse.pyx'])

matrix_rational_dense = Extension('sage.matrix.matrix_rational_dense',
                                  ['sage/matrix/matrix_rational_dense.pyx'],
                                 libraries = ['gmp'])

matrix_integer_sparse = Extension('sage.matrix.matrix_integer_sparse',
                                  ['sage/matrix/matrix_integer_sparse.pyx'],
                                  libraries = ['gmp'])

matrix_rational_sparse = Extension('sage.matrix.matrix_rational_sparse',
                                  ['sage/matrix/matrix_rational_sparse.pyx'],
                                 libraries = ['gmp'])

# TODO -- change to use BLAS at some point.
matrix_integer_dense = Extension('sage.matrix.matrix_integer_dense',
                                 ['sage/matrix/matrix_integer_dense.pyx'],
                                  libraries = ['iml', 'gmp', 'm', BLAS])  # order matters for cygwin!!

matrix_real_double_dense=Extension('sage.matrix.matrix_real_double_dense',
   ['sage/matrix/matrix_real_double_dense.pyx'],libraries=['gsl',BLAS],
   define_macros=[('GSL_DISABLE_DEPRECATED','1')],include_dirs=[SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy'])

matrix_complex_double_dense=Extension('sage.matrix.matrix_complex_double_dense',
   ['sage/matrix/matrix_complex_double_dense.pyx'],libraries=['gsl',BLAS],
   define_macros=[('GSL_DISABLE_DEPRECATED','1')],include_dirs=[SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy'])


solve = Extension('sage.matrix.solve',['sage/matrix/solve.pyx'],libraries = ['gsl',BLAS],define_macros =
   [('GSL_DISABLE_DEPRECATED','1')])

matrix_cyclo_dense = Extension('sage.matrix.matrix_cyclo_dense',
                               ['sage/matrix/matrix_cyclo_dense.pyx'])

matrix_rational_sparse = Extension('sage.matrix.matrix_rational_sparse',
                                   ['sage/matrix/matrix_rational_sparse.pyx'],
                                   libraries = ['gmp'])

matrix_cyclo_sparse = Extension('sage.matrix.matrix_cyclo_sparse',
                                   ['sage/matrix/matrix_cyclo_sparse.pyx'])

#matrix_padic_capped_relative_dense = Extension('sage.matrix.padics.matrix_padic_capped_relative_dense',
#                                               ['sage/matrix/padics/matrix_padic_capped_relative_dense.pyx'])

complex_number = Extension('sage.rings.complex_number',
			    ['sage/rings/complex_number.pyx'],
			    libraries = ['mpfr', 'gmp'])

free_module_element = Extension('sage.modules.free_module_element',
                                ['sage/modules/free_module_element.pyx'])

################ GSL wrapping ######################
gsl_probability=Extension('sage.gsl.probability_distribution',['sage/gsl/probability_distribution.pyx'],libraries=['gsl',BLAS],define_macros=[('GSL_DISABLE_DEPRECATED','1')])
gsl_integration=Extension('sage.gsl.integration',['sage/gsl/integration.pyx'],define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl',BLAS])

gsl_ode = Extension('sage.gsl.ode',['sage/gsl/ode.pyx'],libraries=['gsl',BLAS],define_macros=[('GSL_DISABLE_DEPRECATED','1')])

gsl_fft = Extension('sage.gsl.fft',
                ['sage/gsl/fft.pyx'],
                libraries = ['gsl', BLAS],define_macros=[('GSL_DISABLE_DEPRECATED','1')])

gsl_interpolation = Extension('sage.gsl.interpolation',
                ['sage/gsl/interpolation.pyx'],
                libraries = ['gsl', BLAS],
define_macros=[('GSL_DISABLE_DEPRECATED','1')])

gsl_callback = Extension('sage.gsl.callback',
                ['sage/gsl/callback.pyx'],
                libraries = ['gsl', BLAS]
,define_macros=[('GSL_DISABLE_DEPRECATED','1')])

real_double = Extension('sage.rings.real_double',
                ['sage/rings/real_double.pyx'],
                libraries = ['gsl', 'gmp', BLAS],define_macros=[('GSL_DISABLE_DEPRECATED','1')])

complex_double = Extension('sage.rings.complex_double',
                           ['sage/rings/complex_double.pyx'],
                           libraries = ['gsl', BLAS, 'pari', 'gmp'])

real_double_vector = Extension('sage.modules.real_double_vector',['sage/modules/real_double_vector.pyx'],
                              libraries = ['gsl',BLAS,'pari','gmp'],define_macros = [('GSL_DISABLE_DEPRECAED','1')],include_dirs=[SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy'])

complex_double_vector = Extension('sage.modules.complex_double_vector',['sage/modules/complex_double_vector.pyx'],
                           libraries = ['gsl', BLAS, 'pari', 'gmp'],define_macros=[('GSL_DISABLE_DEPRECATED','1')],include_dirs=[SAGE_ROOT+'/local/lib/python2.5/site-packages/numpy/core/include/numpy'])


vector_integer_dense = Extension('sage.modules.vector_integer_dense',
                                 ['sage/modules/vector_integer_dense.pyx'],
                                 libraries = ['gmp'])

vector_modn_dense = Extension('sage.modules.vector_modn_dense',
                                 ['sage/modules/vector_modn_dense.pyx'])

vector_rational_sparse = Extension('sage.modules.vector_rational_sparse',
                                 ['sage/modules/vector_rational_sparse.pyx'],
                                 libraries = ['gmp'])

vector_rational_dense = Extension('sage.modules.vector_rational_dense',
                                 ['sage/modules/vector_rational_dense.pyx'],
                                 libraries = ['gmp'])

gsl_array = Extension('sage.gsl.gsl_array',['sage/gsl/gsl_array.pyx'],
                libraries=['gsl',BLAS],define_macros=[('GSL_DISABLE_DEPRECATED','1')])

gsl_ode = Extension('sage.gsl.ode',['sage/gsl/ode.pyx'],libraries=['gsl',BLAS],
                define_macros=[('GSL_DISABLE_DEPRECATED','1')])


dwt = Extension('sage.gsl.dwt',['sage/gsl/dwt.pyx'],
                 libraries=['gsl',BLAS],
                 define_macros=[('GSL_DISABLE_DEPRECATED','1')])


sagex_ds = Extension('sage.misc.sagex_ds', ['sage/misc/sagex_ds.pyx'])


#####################################################

ext_modules = [ \
    free_module_element,

    complex_double_vector,
    real_double_vector,

    vector_integer_dense,
    vector_modn_dense,
    vector_rational_dense,

    #vector_rational_sparse,

    pari,

    mwrank,

    ntl,

    matrix,

    matrix_misc,

    matrix_dense,
    matrix_generic_dense,

    matrix_sparse,
    matrix_generic_sparse,

##     matrix_domain_dense,
##     matrix_domain_sparse,

##     matrix_pid_dense,
##     matrix_pid_sparse,

##     matrix_field_dense,
##     matrix_field_sparse,

     matrix_integer_dense,
     matrix_rational_dense,
     matrix_rational_sparse,
     matrix_integer_2x2,
     matrix_integer_sparse,
     matrix_real_double_dense,
     matrix_complex_double_dense,
#     matrix_padic_capped_relative_dense,
     solve,
     linbox,
     matrix_modn_dense,
     matrix_modn_sparse,
     matrix_mod2_dense,
     givaro_gfq, \

     libsingular, \

##     matrix_rational_sparse,

##     matrix_cyclo_dense,
##     matrix_cyclo_sparse,

    dwt,

    gsl_array,
    gsl_ode,
    gsl_fft,
    gsl_interpolation,
    gsl_callback,
    gsl_probability,
    gsl_integration,
    real_double,
    complex_double,
    qd,

    complex_number,

    sagex_ds,

    Extension('sage.media.channels',
              sources = ['sage/media/channels.pyx']), \

    Extension('sage.ext.sig',
              sources = ['sage/ext/sig.pyx']), \

    Extension('sage.ext.arith',
              sources = ['sage/ext/arith.pyx']), \

    Extension('sage.ext.arith_gmp',
              sources = ['sage/ext/arith_gmp.pyx'],
              libraries=['gmp']), \

    Extension('sage.ext.multi_modular',
              sources = ['sage/ext/multi_modular.pyx'],
              libraries=['gmp']), \

    Extension('sage.structure.coerce',
              sources = ['sage/structure/coerce.pyx']), \

    Extension('sage.modular.congroup_pyx',
              sources = ['sage/modular/congroup_pyx.pyx', \
                         'sage/ext/arith.pyx']), \

    Extension('sage.structure.element',
              sources = ['sage/structure/element.pyx']), \

    Extension('sage.categories.morphism',
              sources = ['sage/categories/morphism.pyx']), \

    Extension('sage.categories.functor',
              sources = ['sage/categories/functor.pyx']), \

    Extension('sage.categories.action',
              sources = ['sage/categories/action.pyx']), \

    Extension('sage.modules.module',
              sources = ['sage/modules/module.pyx']), \

    Extension('sage.rings.ring',
              sources = ['sage/rings/ring.pyx']), \

    Extension('sage.rings.polynomial.multi_polynomial',
              sources = ['sage/rings/polynomial/multi_polynomial.pyx']
              ), \

    Extension('sage.rings.polynomial.multi_polynomial_ring_generic',
              sources = ['sage/rings/polynomial/multi_polynomial_ring_generic.pyx']
              ), \

    Extension('sage.rings.polynomial.multi_polynomial_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_libsingular.pyx'],
              libraries = ['gmp', 'm', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx'],
              language="c++",
              include_dirs=[SAGE_ROOT +'/local/include/singular']), \

    Extension('sage.rings.polynomial.multi_polynomial_ideal_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_ideal_libsingular.pyx'],
              libraries = ['gmp', 'm', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx'],
              language="c++",
              include_dirs=[SAGE_ROOT +'/local/include/singular']), \

    Extension('sage.groups.group',
              sources = ['sage/groups/group.pyx']), \

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

    Extension('sage.misc.cython_c',
              sources = ['sage/misc/cython_c.pyx']), \

    Extension('sage.rings.real_mpfr',
              sources = ['sage/rings/real_mpfr.pyx', 'sage/rings/ring.pyx'],
              libraries = ['mpfr', 'pari', 'gmp']), \

    Extension('sage.rings.real_mpfi',
              sources = ['sage/rings/real_mpfi.pyx'],
              libraries = ['mpfi', 'mpfr', 'gmp']), \

    Extension('sage.rings.integer',
              sources = ['sage/ext/arith.pyx', 'sage/rings/integer.pyx', \
                         'sage/ext/mpn_pylong.c', 'sage/ext/mpz_pylong.c'],
              libraries=['ntl', 'gmp']), \

    Extension('sage.rings.integer_ring',
              sources = ['sage/rings/integer_ring.pyx'],
              libraries=['ntl', 'gmp']), \

    Extension('sage.rings.padics.pow_computer',
              sources = ['sage/rings/padics/pow_computer.pyx'],
              libraries=['gmp']),
    Extension('sage.rings.padics.local_generic_element',
              sources = ['sage/rings/padics/local_generic_element.pyx']),
    Extension('sage.rings.padics.padic_generic_element',
              sources = ['sage/rings/padics/padic_generic_element.pyx']),
    Extension('sage.rings.padics.padic_base_generic_element',
              sources = ['sage/rings/padics/padic_base_generic_element.pyx']),
    Extension('sage.rings.padics.padic_fixed_mod_element',
              sources = ['sage/rings/padics/padic_fixed_mod_element.pyx', \
                         'sage/rings/padics/padic_generic_element.c'],
              libraries=['gmp']),
    Extension('sage.rings.padics.padic_capped_absolute_element',
              sources = ['sage/rings/padics/padic_capped_absolute_element.pyx', \
                         'sage/rings/padics/padic_generic_element.c'],
              libraries=['gmp']),
    Extension('sage.rings.padics.padic_capped_relative_element',
              sources = ['sage/rings/padics/padic_capped_relative_element.pyx', \
                         'sage/rings/padics/padic_generic_element.c'],
              libraries=['gmp']),


    Extension('sage.rings.memory', \
              sources = ['sage/rings/memory.pyx'], \
              libraries=['gmp']), \

    Extension('sage.rings.bernoulli_mod_p',
              sources = ['sage/rings/bernoulli_mod_p.pyx', 'sage/ext/arith.pyx'],
              libraries=['ntl'],
              include_dirs=['sage/libs/ntl/']), \

    Extension('sage.schemes.hyperelliptic_curves.frobenius',
                 sources = ['sage/schemes/hyperelliptic_curves/frobenius.pyx',
                            'sage/schemes/hyperelliptic_curves/frobenius_cpp.cpp'],
                 libraries = ['ntl', 'stdc++', 'gmp'],
                 language = 'c++',
                 include_dirs=['sage/libs/ntl/']), \

    Extension('sage.rings.polynomial.polynomial_compiled',
               sources = ['sage/rings/polynomial/polynomial_compiled.pyx']), \

    Extension('sage.rings.polynomial.polynomial_element',
              sources = ['sage/rings/polynomial/polynomial_element.pyx']), \

    Extension('sage.rings.power_series_ring_element',
              sources = ['sage/rings/power_series_ring_element.pyx']), \

    Extension('sage.rings.power_series_poly',
              sources = ['sage/rings/power_series_poly.pyx']), \

    Extension('sage.rings.power_series_mpoly',
              sources = ['sage/rings/power_series_mpoly.pyx']), \

    Extension('sage.rings.laurent_series_ring_element',
              sources = ['sage/rings/laurent_series_ring_element.pyx']), \

    Extension('sage.rings.rational',
              sources = ['sage/rings/rational.pyx',
                         'sage/ext/arith.pyx', \
                         'sage/rings/integer.pyx', \
                         'sage/ext/mpn_pylong.c', 'sage/ext/mpz_pylong.c'],
              libraries=['ntl', 'gmp']), \

    Extension('sage.rings.sparse_poly',
              sources = ['sage/rings/sparse_poly.pyx'],
              libraries=['gmp']), \

    Extension('sage.rings.polynomial.polydict',
              sources = ['sage/rings/polynomial/polydict.pyx']), \

    Extension('sage.rings.number_field.number_field_element',
              sources = ['sage/rings/number_field/number_field_element.pyx'],
              libraries=['ntl','gmp'],
              language = 'c++'), \

    Extension('sage.misc.search',
              ['sage/misc/search.pyx']), \

    Extension('sage.misc.reset',
              ['sage/misc/reset.pyx']), \

    Extension('sage.calculus.var',
              ['sage/calculus/var.pyx']), \

    Extension('sage.modular.modsym.heilbronn',
              ['sage/modular/modsym/heilbronn.pyx',
               'sage/modular/modsym/p1list.pyx',
               'sage/ext/arith.pyx'],
              libraries = ['gmp', 'm']), \

    Extension('sage.modular.modsym.p1list',
              ['sage/modular/modsym/p1list.pyx',
               'sage/ext/arith.pyx'],
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
              libraries = ['gmp', 'mpfr'],
              language='c++'
              ), \

    Extension('sage.graphs.graph_fast',
              ['sage/graphs/graph_fast.pyx'],
              libraries = ['gmp']
              ), \

    Extension('sage.graphs.graph_isom',
              ['sage/graphs/graph_isom.pyx']
              ), \

    Extension('sage.graphs.bruhat_sn',
              ['sage/graphs/bruhat_sn.pyx']
              ), \

    ]


#mpc = Extension('sage.rings.mpc',
#              sources = ['sage/rings/mpc.pyx', 'sage/rings/ring.pyx'],
#              libraries = ['mpc', 'mpfr', 'gmp'])


extra_compile_args = [ ]

if NO_WARN and distutils.sysconfig.get_config_var('CC').startswith("gcc"):
    extra_compile_args = ['-w']

if DEVEL:
    extra_compile_args.append('-ggdb')
    ext_modules.append(hanke)
    #ext_modules.append(mpc)

for m in ext_modules:
    m.libraries = ['csage'] + m.libraries + ['stdc++', 'ntl']
    m.library_dirs += ['%s/lib' % SAGE_LOCAL]


######################################################################
# CODE for generating C/C++ code from Cython and doing dependency
# checking, etc.  In theory distutils would run Cython, but I don't
# trust it at all, and it won't have the more sophisticated dependency
# checking that we need.
######################################################################

def is_older(file1, file2):
    """
    Return True if either file2 does not exist or is older than file1.

    If file1 does not exist, always return False.
    """
    if not os.path.exists(file1):
        return False
    if not os.path.exists(file2):
        return True
    if os.path.getmtime(file2) < os.path.getmtime(file1):
        return True
    return False


def need_to_cython(filename, outfile):
    """
    INPUT:
        filename -- The name of a cython file in the SAGE source tree.
        outfile -- The name of the corresponding c or cpp file in the build directory.

    OUTPUT:
        bool -- whether or not outfile must be regenerated.
    """
    if is_older(filename, outfile):   # outfile is older than filename
        return True
    base =  os.path.splitext(filename)[0]
    pxd = base+'.pxd'
    if is_older(pxd, outfile):   # outfile is older than pxd file (if it exists)
        return True

    ## comment this out to turn on dependency checking!
    # return False

    # Now we look inside the file to see what it cimports or include.
    # If any of these files are newer than outfile, we rebuild
    # outfile.
    S = open(filename).readlines()
    if os.path.exists(pxd):
        S += open(pxd).readlines()
    # Take the lines that begin with cimport (it won't hurt to
    # have extra lines)
    C = [x.strip() for x in S if 'cimport' in x]
    for A in C:
        # Deduce the full module name.
        # The only allowed forms of cimport are:
        #        cimport a.b.c.d
        #        from a.b.c.d cimport e
        # Also, the cimport can be both have no dots (current directory) or absolute.
        # The multiple imports with parens, e.g.,
        #        import (a.b.c.d, e.f.g)
        # of Python are not allowed in Cython.
        # In both cases, the module name is the second word if we split on whitespace.
        try:
            A = A.strip().split()[1]
        except IndexError:
            # illegal statement or not really a cimport (e.g., cimport could
            # be in a comment or something)
            continue

        # Now convert A to a filename, e.g., a/b/c/d
        #
        if '.' in A:
            # It is an absolute cimport.
            A = A.replace('.','/') + '.pxd'
        else:
            # It is a relative cimport.
            A =  os.path.split(base)[0] + '/' + A + '.pxd'
        # Check to see if a/b/c/d.pxd exists and is newer than filename.
        # If so, we have to regenerate outfile.  If not, we're safe.
        if os.path.exists(A) and is_older(A, outfile):
            print "\nRegenerating %s because it depends on %s."%(outfile, A)
            return True # yep we must rebuild

    # OK, next we move on to include pxi files.
    # If they change, we likewise must rebuild the pyx file.
    I = [x for x in S if 'include' in x]
    # The syntax for include is like this:
    #       include "../a/b/c.pxi"
    # I.e., it's a quoted *relative* path to a pxi file.
    for A in I:
        try:
            A = A.strip().split()[1]
        except IndexError:
            # Illegal include statement or not really an include
            # (e.g., include could easily be in a comment or
            # something).  No reason to crash setup.py!
            continue
        # Strip the quotes from either side of the pxi filename.
        A = A.strip('"').strip("'")
        # Now take filename (the input argument to this function)
        # and strip off the last part of the path and stick
        # A onto it to get the filename of the pxi file relative
        # where setup.py is being run.
        A = os.path.split(filename)[0] + '/' + A
        # Finally, check to see if filename is older than A
        if os.path.exists(A) and is_older(A, outfile):
            print "\nBuilding %s because it depends on %s."%(outfile, A)
            return True

    # We do not have to rebuild.
    return False

def process_cython_file(f, m):
    """
    INPUT:
        f -- file name
        m -- Extension module description (i.e., object of type Extension).
    """
    # This is a cython file, so process accordingly.
    pyx_inst_file = '%s/%s'%(SITE_PACKAGES, f)
    if is_older(f, pyx_inst_file):
        print "%s --> %s"%(f, pyx_inst_file)
        os.system('cp %s %s 2>/dev/null'%(f, pyx_inst_file))
    outfile = f[:-4] + ".c"
    if m.language == 'c++':
        outfile += 'pp'

    if need_to_cython(f, outfile):
        cmd = "cython --embed-positions -I%s %s"%(os.getcwd(), f)
        print cmd
        ret = os.system(cmd)
        if ret != 0:
            print "sage: Error running cython."
            sys.exit(1)
        # If the language for the extension is C++,
        # then move the resulting output file to have the correct extension.
        # (I don't know how to tell Cython to do this automatically.)
        if m.language == 'c++':
            os.system('mv %s.c %s'%(f[:-4], outfile))
    return [outfile]


def cython(ext_modules):
    for m in ext_modules:
        m.extra_compile_args += extra_compile_args
#        m.extra_link_args += extra_link_args
        new_sources = []
        for i in range(len(m.sources)):
            f = m.sources[i]
            s = open(f).read()
            if f[-4:] == ".pyx":
                new_sources += process_cython_file(f, m)
            else:
                new_sources.append(f)
        m.sources = new_sources



if not sdist:
    cython(ext_modules)

setup(name        = 'sage',

      version     =  SAGE_VERSION,

      description = 'SAGE: System for Algebra and Geometry Experimentation',

      license     = 'GNU Public License (GPL)',

      author      = 'William Stein',

      author_email= 'wstein@gmail.com',

      url         = 'http://www.sagemath.org',

      packages    = ['sage',

                     'sage.algebras',

                     'sage.catalogue',

                     'sage.categories',

                     'sage.coding',

                     'sage.combinat',

                     'sage.crypto',

                     'sage.databases',

                     'sage.ext',

                     'sage.calculus',

                     'sage.functions',

                     'sage.geometry',

                     'sage.games',

                     'sage.gsl',

                     'sage.graphs',

                     'sage.groups',
                     'sage.groups.abelian_gps',
                     'sage.groups.matrix_gps',
                     'sage.groups.perm_gps',

                     'sage.interfaces',

                     'sage.lfunctions',

                     'sage.libs',
                     'sage.libs.hanke',
                     'sage.libs.linbox',
                     'sage.libs.mwrank',
                     'sage.libs.ntl',
                     'sage.libs.pari',
                     'sage.libs.singular',

                     'sage.matrix',
#                     'sage.matrix.padics',
                     'sage.media',
                     'sage.misc',

                     'sage.modules',

                     'sage.modular',
                     'sage.modular.abvar',
                     'sage.modular.hecke',
                     'sage.modular.modform',
                     'sage.modular.modsym',
                     'sage.modular.ssmod',

                     'sage.monoids',

                     'sage.plot',
                     'sage.plot.mpl3d',

                     'sage.probability',

                     'sage.quadratic_forms',
                     'sage.quadratic_forms.genera',

                     'sage.rings',
                     'sage.rings.number_field',
                     'sage.rings.padics',
                     'sage.rings.polynomial',
                     'sage.rings.polynomial.padics',

                     'sage.tests',

                     'sage.sets',

                     'sage.schemes',
                     'sage.schemes.generic',
                     'sage.schemes.jacobians',
                     'sage.schemes.plane_curves',
                     'sage.schemes.plane_quartics',
                     'sage.schemes.elliptic_curves',
                     'sage.schemes.hyperelliptic_curves',

                     'sage.server',
                     'sage.server.server1',
                     'sage.server.notebook',
                     'sage.server.notebook.compress',
                     'sage.server.wiki',
                     'sage.server.trac',

                     'sage.structure',

                     'sage.dsage',
                     'sage.dsage.tests',
                     'sage.dsage.database',
                     'sage.dsage.database.tests',
                     'sage.dsage.server',
                     'sage.dsage.server.tests',
                     'sage.dsage.interface',
                     'sage.dsage.interface.tests',
                     'sage.dsage.errors',
                     'sage.dsage.twisted',
                     'sage.dsage.twisted.tests',
                     'sage.dsage.dist_functions',
                     'sage.dsage.dist_functions.tests',
                     'sage.dsage.misc',
                     'sage.dsage.misc.tests',
                     'sage.dsage.scripts'
                     ],

      scripts = ['sage/dsage/scripts/dsage_server.py',
                 'sage/dsage/scripts/dsage_worker.py',
                 'sage/dsage/scripts/dsage_setup.py'
                ],

      ext_modules = ext_modules,
      include_dirs = include_dirs)



