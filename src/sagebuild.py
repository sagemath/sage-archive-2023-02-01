##########################################################################
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#
#    * Neither the name of Gary Furnish nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##########################################################################
# -- HOWTO ADD NEW CUSTOM CYTHON EXTENSIONS -- ##
# If the extension does not require custom compiler or linker options,
# PBuild automatically handles everything.  Otherwise, a config_gcc_file
# is needed in the extensions section of this file.
# If the compiled file is C++, set languagge='C++'
# If define macros are needed, set define_macros=deflist (where deflist
# is some list of define macros such as ['def1','def2'])
# If custom include directories are needed, set include_dirs to a list
# of include directories such as ['include','/usr/include/foo']
# If custom library directories are needed, set library_dirs to a list
# of library directories such as ['/lib','/lib/foo']
# If custom libraries are needed, set libraries to a list of libraries such
# libraries=['csage',ntl'].  It is important to note that the preceeding lib
# and the ending .so or .a are not included.
# The first two arguments of config_gcc_file should always be env, pyx_pextn_dict
# The third argument is the path to the file
# Thus a sample call might be
# def config_gcc_file(env, pyx_pextn_dict, path, language = 'C++', define_macros = defmacros, \
# include_dirs = incdirlist, library_dirs = libdirlist, libraries = liblist)
##########################################################################

import os, signal, sys, time, thread, threading, tempfile, shutil

from build.all import *
verbose=1

cacheinfo = None

def abs_sage_path(env, path):
    """
        A convinent function for making a function relative to devel/sage relative to $SAGE_ROOT
    """
    return 'devel/sage/' + path
def real_sage_root(env, path):
    """
        Return the real path of $SAGE_ROOT
    """
    return os.path.realpath(env.options['SAGE_ROOT'] + path)
def real_sage_local(env, path):
    """
        Return the real path of $SAGE_LOCAL
    """
    return os.path.realpath(env.options['SAGE_LOCAL'] + path)
def dict_insert_src(dict, dictlookup, src):
    """
        A utility function that adds a source file to the sources for a given extension if it exists in a dictionary
    """
    obj = dict.get(dictlookup, 0)
    if obj == 0:
        return
    obj.sources.insert(0,src)

def build_sage_clean(env):
    """
        This function gets called when we want sage to build all
    """
    try:
        os.remove(env.options['SAGE_ROOT'] + "/sagebuild.chc")
    except:
        pass
    shutil.rmtree('devel/sage/build', True)

def buildsage(env, gccc):
    """
        This function actually builds sage
    """
    TM = taskmanager.TM
    #acquire a list of all pyx files in the tree
    efw = extfilewalker()
    efw.addcallback('.pyx',lambda x: True)
    pyx_list = efw.walk('devel/sage/sage')
    for entry in pyx_list:
        if entry.find('#')!=-1:
            pyx_list.remove(entry)
    #set the absolute devel directory if necessary
    develdir = os.path.realpath('devel/sage/')
    if verbose>100:
        print 'Devel Dir: ' + develdir
    #Create the cython compiler
    cyc = Cython_compiler()
    #We want -I. so that it correctly finds include files by default
    cyc.set_option_val('-I'+'.',None)
    #Create the build extension dict for all cython extensions
    pyx_pextn_dict = cyc.get_build_extensions(env, pyx_list, cwd = 'devel/sage')

    ###################################################################
    #Configure all extensions
    ###################################################################

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/flint/fmpz_poly.pyx",include_dirs = [real_sage_local(env, 'include/FLINT/')], libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++"],prop_options={ str(GCC_extension_object): {'-std':'c99' } } )

    ntllibs = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"]

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZ.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZX.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZ_pContext.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZ_p.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZ_pX.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZ_pEContext.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZ_pE.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_ZZ_pEX.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_lzz_pContext.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_lzz_p.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_lzz_pX.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_GF2.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_GF2X.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_GF2EContext.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_GF2E.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_GF2EX.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_mat_ZZ.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_mat_GF2E.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/ntl/ntl_GF2EContext.pyx",language="C++", libraries=ntllibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/mwrank/mwrank.pyx", language='C++', define_macros=[("NTL_ALL",None)], libraries=["curvesntl", "g0nntl", "jcntl", "rankntl", "ntl", "gmp", "gmpxx", "stdc++", "m", "pari"] )

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/pari/gen.pyx",libraries=['pari', 'gmp'])

    cremonalibs = ['g0nntl', 'jcntl', 'gmpxx', 'ntl', 'gmp', 'm', 'stdc++', 'pari', 'curvesntl']

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/cremona/mat.pyx",language="C++", define_macros=[("NTL_ALL",None)], libraries=cremonalibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/cremona/homspace.pyx",language="C++", define_macros=[("NTL_ALL",None)], libraries=cremonalibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/cremona/newforms.pyx",language="C++", define_macros=[("NTL_ALL",None)], libraries=cremonalibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/rings/finite_field_givaro.pyx",language="C++", libraries=['givaro', 'gmpxx', 'gmp', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/finite_field_ntl_gf2e.pyx',language="C++", libraries=['ntl', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/rings/real_rqdf.pyx",language="C++", libraries=['qd', 'm', 'stdc++','gmp','mpfr','qd', 'csage' ])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/misc.pyx',libraries=['mpfr','gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_integer_2x2.pyx',libraries=['gmp'])

    BLAS = env.options['BLAS']
    BLAS2 = env.options['BLAS2']

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/libs/linbox/linbox.pyx',language="C++", libraries=['linboxsage', 'ntl', 'linbox', 'gmp', 'gmpxx', 'stdc++', 'givaro', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/libs/singular/singular.pyx',language='C++', include_dirs=[real_sage_local(env,'include/singular')], libraries=['m', 'readline', 'singular', 'singfac', 'singcf', 'omalloc', 'givaro', 'gmpxx', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/libs/fplll/fplll.pyx',language='C++', include_dirs = [real_sage_local(env, 'include/fplll')], libraries=['gmp', 'mpfr', 'stdc++', 'fplll'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_modn_dense.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_mod2_dense.pyx',libraries = ['gmp','m4ri'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_rational_dense.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_integer_sparse.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_rational_sparse.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_integer_dense.pyx',libraries = ['iml', 'gmp', 'm', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_real_double_dense.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], include_dirs = [real_sage_local(env,'lib/python2.5/site-packages/numpy/core/include/numpy')], libraries=[BLAS, BLAS2, 'gsl'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/change_ring.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=[BLAS, BLAS2, 'gsl', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_complex_double_dense.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], include_dirs=[real_sage_local(env,'lib/python2.5/site-packages/numpy/core/include/numpy')], libraries=['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/solve.pyx',define_macros = [('GSL_DISABLE_DEPRECATED','1')], libraries = ['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_mpolynomial_dense.pyx',language="C++", include_dirs = [real_sage_local(env,'include/singular')], libraries = ['m', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/complex_number.pyx',libraries = ['mpfr', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/probability_distribution.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/integration.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/ode.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', BLAS])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/fft.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/interpolation.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/callback.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/real_double.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', 'gmp', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/complex_double.pyx',libraries = ['gsl', BLAS, BLAS2, 'pari', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/modules/real_double_vector.pyx',define_macros = [('GSL_DISABLE_DEPRECAED','1')], include_dirs=[real_sage_local(env,'lib/python2.5/site-packages/numpy/core/include/numpy')], libraries = ['gsl', BLAS, BLAS2, 'pari','gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/modules/complex_double_vector.pyx',define_macros = [('GSL_DISABLE_DEPRECAED','1')], include_dirs=[real_sage_local(env,'lib/python2.5/site-packages/numpy/core/include/numpy')], libraries = ['gsl', BLAS, BLAS2, 'pari','gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/modules/vector_integer_dense.pyx',libraries=['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/modules/vector_rational_dense.pyx',libraries=['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/gsl_array.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl', BLAS, BLAS2])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/ode.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl',BLAS])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/gsl/dwt.pyx',define_macros=[('GSL_DISABLE_DEPRECATED','1')], libraries=['gsl',BLAS])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/libs/symmetrica/symmetrica.pyx',libraries=['symmetrica'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/ext/arith_gmp.pyx',libraries=['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/ext/multi_modular.pyx',libraries=['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/polynomial/multi_polynomial_libsingular.pyx',language='C++', include_dirs = [real_sage_local(env,'include/singular')], libraries = ['m', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/polynomial/multi_polynomial_ideal_libsingular.pyx',language='C++', include_dirs = [real_sage_local(env,'include/singular')], libraries = ['m', 'readline', 'singular', 'singcf', 'singfac', 'omalloc', 'givaro', 'gmpxx', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/structure/sage_object.pyx',libraries=['ntl'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/real_mpfr.pyx',libraries = ['mpfi', 'mpfr', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/real_mpfi.pyx',libraries = ['mpfi', 'mpfr', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/complex_interval.pyx',libraries = ['mpfi', 'mpfr', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/integer.pyx',libraries = ['ntl', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/integer_ring.pyx',libraries = ['ntl', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/libs/libecm.pyx',libraries = ['ecm', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/pow_computer.pyx',language='C++', libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/pow_computer_ext.pyx',language='C++', libraries = ["csage", "ntl", "gmp", "gmpxx", "m", "stdc++"])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_fixed_mod_element.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_capped_absolute_element.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_capped_relative_element.pyx',libraries = ['gmp', 'csage'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_ext_element.pyx',language = 'C++', libraries = ['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_ZZ_pX_element.pyx',language = 'C++', libraries = ['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_ZZ_pX_FM_element.pyx',language = 'C++', libraries = ['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_ZZ_pX_CR_element.pyx',language = 'C++', libraries = ['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_ZZ_pX_CA_element.pyx',language = 'C++', libraries = ['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/padics/padic_printing.pyx',language = 'C++', libraries = ['gmp', 'ntl', 'csage', 'gmpxx', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/memory.pyx',libraries = ['gmp','stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/bernoulli_mod_p.pyx',language = 'C++', include_dirs = [abs_sage_path(env,'devel/sage/libs/ntl/')], libraries=['ntl','stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/schemes/hyperelliptic_curves/hypellfrob.pyx',language = 'C++', include_dirs = [abs_sage_path(env,'devel/sage/libs/ntl/'), abs_sage_path(env,'sage/schemes/hyperelliptic_curves/hypellfrob/')], libraries = ['ntl','stdc++', 'gmp', 'zn_poly'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/polynomial/polynomial_integer_dense_flint.pyx', language = 'C++', include_dirs = [real_sage_local(env, 'include/FLINT')], libraries=["csage", "flint", "gmp", "gmpxx", "ntl"])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/polynomial/polynomial_integer_dense_ntl.pyx',language = 'C++', include_dirs = [abs_sage_path(env,'devel/sage/libs/ntl/')], libraries = ['ntl','stdc++', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/polynomial/polynomial_modn_dense_ntl.pyx',language = 'C++', include_dirs = [abs_sage_path(env,'devel/sage/libs/ntl/')], libraries = ['ntl','stdc++', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/rational.pyx',libraries = ['ntl','gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/sparse_poly.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/polynomial/real_roots.pyx',libraries = ['mpfr', 'qd'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/number_field/number_field_element.pyx',language = 'C++', libraries=['ntl','gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/number_field/number_field_element_quadratic.pyx',language = 'C++', libraries=['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/number_field/totallyreal_data.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/modular/modsym/heilbronn.pyx',libraries = ['gmp','m'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/modular/modsym/p1list.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/integer_mod.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/combinat/expnums.pyx',libraries = ['gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/combinat/partitions.pyx',language='C++', libraries = ['qd', 'gmp', 'mpfr'])

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/combinat/matrices/dancing_links.pyx",language="C++")

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/graphs/graph_fast.pyx',libraries = ['gmp'])

    #config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/symbolics/symbolicarithmetic.pyx',libraries = ['glite'])
    #config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/symbolics/symbolicmularithmetic.pyx',libraries = ['glite'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/polynomial/pbori.pyx',language='C++', include_dirs = [real_sage_local(env,'include/cudd'),real_sage_local(env,'include/polybori'), real_sage_local(env,'include/polybori/groebner')], libraries=['polybori','pboriCudd','groebner'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/finance/time_series.pyx',include_dirs=[real_sage_local(env,'lib/python2.5/site-packages/numpy/core/include/numpy')])

    ###################################################################
    #End standard configurations
    ###################################################################

    tempdir = os.path.realpath(env.options['SAGE_ROOT'])

    #A dictionary of gcc objects extensions
    gcceo_dict = { }
    #A dictionary of gcc shared object extensions
    gcceso_dict = { }

    #Manually create the mwrank C extension
    mwrankcc = GCC_extension_object(gccc, env, ["devel/sage/sage/libs/mwrank/wrap.cc"], tempdir+"/devel/sage/build/temp/sage/libs/mwrank", language='C++', define_macros=[("NTL_ALL",None)], options = { '-fPIC':None } )

    #Create the Hypellfrob extensions
    hypellfrob_cpp = GCC_extension_object(gccc, env,["devel/sage/sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.cpp"], tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob/", language='C++', include_dirs=[abs_sage_path(env, 'sage/libs/ntl'), abs_sage_path(env,'sage/schemes/hyperelliptic_curves/hypellfrob/') ], options = { '-fPIC':None })
    recurrences_zn_poly = GCC_extension_object(gccc, env,[abs_sage_path(env,"sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.cpp")], tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob/", language='C++', include_dirs=[abs_sage_path(env, 'sage/libs/ntl'), abs_sage_path(env,'sage/schemes/hyperelliptic_curves/hypellfrob/') ], options = { '-fPIC':None } )
    recurrences_ntl = GCC_extension_object(gccc, env,["devel/sage/sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.cpp"], tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob/", language='C++', include_dirs=[abs_sage_path(env, 'sage/libs/ntl'), abs_sage_path(env,'sage/schemes/hyperelliptic_curves/hypellfrob/') ], options = { '-fPIC':None } )

    #Manually create the C paritions extension
    partitions_c = GCC_extension_object(gccc, env,["devel/sage/sage/combinat/partitions_c.cc"], tempdir+"/devel/sage/build/temp/sage/combinat", language='C++', options = { '-fPIC':None } )

    #If we don't already have a temporary build directory for planarity, make it
    safemkdirs(tempdir+"/devel/sage/build/temp/sage/graphs/planarity")

    #Locate all files that can possibley have
    depfw = extfilewalker()
    depfw.addcallback('.pyx',lambda x: True)
    depfw.addcallback('.pxd',lambda x: True)
    depfw.addcallback('.pxi',lambda x: True)
    depfw.addcallback('.py',lambda x: x.find(".doctest") ==-1 )
    dep_list = depfw.walk('devel/sage/sage')
    for entry in dep_list:
        if entry.find('#')!=-1:
            dep_list.remove(entry)
    #Set the dependency locater function
    funcdict = { '.pyx':get_cython_file_deps, '.pxd':get_cython_file_deps, 'pxi':get_cython_file_deps }
    #Set possible include directories for various file extensions
    includedict = { '.pyx':[abs_sage_path(env,'')], '.pxd':[abs_sage_path(env,'')], '.pxi':[abs_sage_path(env,'')] }
    #Set cacheinfo (used for saving the state at the end of the build)
    global cacheinfo
    cacheinfo = get_compile_list(dep_list, funcdict, includedict, "sagebuild.chc", env.options['SAGE_ROOT'])
    #Add the task that runs on exit to write out the build cache
    TM.addfinishtask(sagebuild_exit)
    pyx_ext_dict  = { }
    for x in cacheinfo[0]:
        if os.path.splitext(x)[1]==".pyx":
            pyx_ext_dict[x]=pyx_pextn_dict[x]
    if verbose>100:
        print '---------------PYX_EXT_DICT------------'
        print pyx_ext_dict
        print '----------------------------------------'

    #If we are building planarity
    if pyx_ext_dict.has_key('devel/sage/sage/graphs/planarity.pyx'):
        #Locate all planarity C files
        efwplan = extfilewalker()
        efwplan.addcallback('.c',lambda x: True)
        planarity_list = efwplan.walk('devel/sage/sage/graphs/planarity')
        planarity_ext = list()
        for x in planarity_list:
            #Create their extensions and enable them
            ext =  GCC_extension_object(gccc, env, [x], tempdir+"/devel/sage/build/temp/sage/graphs/planarity", options = { '-fPIC':None } )
            ext.generate_action(env).enable()
            planarity_ext.append(ext)

    #If we are building mwrank.pyx, build the mwrank c extensions
    if pyx_ext_dict.has_key("devel/sage/sage/libs/mwrank/mwrank.pyx"):
        mwrankcc.generate_action(env).enable()
    else:
        pass

    #If we are building hypellfrob, make sure the temporary directories exist and enable the C extensions
    if pyx_ext_dict.has_key('devel/sage/sage/schemes/hyperelliptic_curves/hypellfrob.pyx'):
        try:
            os.makedirs(hypellfrob_cpp.get_out_dir())
        except:
            pass
        hypellfrob_cpp.generate_action(env).enable()
        recurrences_ntl.generate_action(env).enable()
        recurrences_zn_poly.generate_action(env).enable()

    #If we are building partions.pyx, we should build the partitions C extension
    if pyx_ext_dict.has_key('devel/sage/sage/combinat/partitions.pyx'):
        partitions_c.generate_action(env).enable()

    #For all Cython extensions, create gcc object file extensions
    for x in pyx_ext_dict.values():
        q = x.generate_action(env)
        gcc_temp_dir = tempdir + '/devel/sage/build/temp/sage/' + (os.path.split(x.sources[0])[0])[len('devel/sage/sage/'):]
        gcceo = GCC_extension_object(gccc, env, [x], gcc_temp_dir, options = { '-fPIC':None } )
        gcceo_dict[x.sources[0]] = gcceo
        gccaction = gcceo.generate_action(env)
        gccaction.dep_register(q)
        q.enable()

    #Create temporary directories for all object file extensions and enable them
    for gcceo in gcceo_dict.values():
        try:
            os.makedirs(gcceo.get_out_dir())
        except:
            pass
        gccaction = gcceo.generate_action(env)
        gccaction.enable()

    #For all object file extensions, create shared object file extensions and enable them
    for gcceo in gcceo_dict.values():
        gcc_build_dir = develdir + '/build' + (os.path.split(gcceo.sources[0])[0])[len('devel/sage'):]
        gcceso = GCC_extension_shared_object(gccc,env, [gcceo], gcc_build_dir)
        gcceso_dict[gcceso.sources[0] ] = gcceso
        gcceso.generate_action(env).dep_register(gcceo.generate_action(env))

    #Link mwrank C files if needed
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/libs/mwrank/mwrank.o", mwrankcc.outfile)
    #Setup dependencies for mwrank if needed
    try:
        gcceso_dict[tempdir + "/devel/sage/build/temp/sage/libs/mwrank/mwrank.o"].generate_action(env).dep_register(mwrankcc.generate_action(env))
    except KeyError:
        pass

    #Link hypellfrob C files if needed
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o", hypellfrob_cpp.outfile)
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o", recurrences_ntl.outfile)
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o", recurrences_zn_poly.outfile)
    #Setup dependencies for hypellfrob shared object if needed
    try:
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o"].generate_action(env).dep_register(hypellfrob_cpp.generate_action(env))
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o"].generate_action(env).dep_register(recurrences_ntl.generate_action(env))
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o"].generate_action(env).dep_register(recurrences_zn_poly.generate_action(env))
    except KeyError:
        pass

    #Link parition C file if needed
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/combinat/partitions.o", partitions_c.outfile)
    #Setup dependencies for parition if needed
    try:
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/combinat/partitions.o"].generate_action(env).dep_register(mwrankcc.generate_action(env))
    except KeyError:
        pass

    #If we are building planarity, then we should link the C extensions into the planarity shared object
    if pyx_ext_dict.has_key('devel/sage/sage/graphs/planarity.pyx'):
        for x in planarity_ext:
            dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/graphs/planarity.o",x.outfile)

    #Make output directorys in site-packages if they don't already exist
    for gcceso in gcceso_dict.values():
        try:
            os.makedirs(gcceso.get_out_dir())
        except:
            pass
        gccesoaction = gcceso.generate_action(env)
        gccesoaction.enable()


    if verbose>10:
        print 'Copying *.py files'
    #Copy all python files to site-packages
    for filenm in cacheinfo[0]:
        if os.path.splitext(filenm)[1]==".py":
            filepart = (os.path.split(filenm)[0])
            newdir = env.options['SAGE_ROOT'] + '/devel/sage/build/sage/' + filepart.replace('devel/sage/sage','')
            filename = os.path.split(filenm)[1]
            if filename[0]=='.':
                continue
            safemkdirs(newdir)
            cmd = 'cp %s %s/%s' % (filenm, newdir, filename)
            if verbose>30:
                print cmd
            os.system( cmd )
    #Setup the site-packages symlink if it doesn't already exist
    safesymlink('../../../../devel/sage/build/sage','local/lib/python/site-packages/sage')
    #Handle DSage
    safemkdirs('local/dsage')
    safesymlink('../../devel/sage/sage/dsage/web', 'local/dsage/web')

def sagebuild_exit():
    """
        This gets called once the entire build system is finished.
        Thus we stick the calls to update the cached values in here
        so they only happen on success
    """
    os.chdir(os.environ['SAGE_ROOT'])
    global cacheinfo
    if verbose>10:
        print 'Writing Sage Cache'
    write_cache(cacheinfo[0], cacheinfo[1], cacheinfo[2], cacheinfo[3], cacheinfo[4], cacheinfo[5], "sagebuild.chc")
