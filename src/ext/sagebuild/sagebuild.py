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
import os, signal, sys, time, thread, threading, tempfile, pprint

from build.all import *
verbose=1

cacheinfo = None

def abs_sage_path(env, path):
    return 'devel/sage/' + path
def real_sage_root(env, path):
    return os.path.realpath(env.options['SAGE_ROOT'] + path)
def real_sage_local(env, path):
    return os.path.realpath(env.options['SAGE_LOCAL'] + path)
def dict_insert_src(dict, dictlookup, src):
    obj = dict.get(dictlookup, 0)
    if obj == 0:
        return
    obj.sources.insert(0,src)

def build_sage_clean(env):
    try:
        os.remove(env.options['SAGE_ROOT'] + "/sagebuild.chc")
    except:
        pass

def buildsage(env, gccc):
    TM = taskmanager.TM
    efw = extfilewalker()
    efw.addcallback('.pyx',lambda x: True)
    pyx_list = efw.walk('devel/sage/sage')
    develdir = os.path.realpath('devel/sage/')
    if verbose>100:
        print 'Devel Dir: ' + develdir
    cyc = Cython_compiler()
    cyc.set_option_val('-I'+'.',None)
    pyx_pextn_dict = cyc.get_build_extensions(env, pyx_list, cwd = 'devel/sage')

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

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/mwrank/mwrank.pyx",define_macros=[("NTL_ALL",None)], libraries=["curvesntl", "g0nntl", "jcntl", "rankntl", "ntl", "gmp", "gmpxx", "stdc++", "m", "pari"] )

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/pari/gen.pyx",libraries=['pari', 'gmp'])

    cremonalibs = ['g0nntl', 'jcntl', 'gmpxx', 'ntl', 'gmp', 'm', 'stdc++', 'pari', 'curvesntl']

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/cremona/mat.pyx",language="C++", define_macros=[("NTL_ALL",None)], libraries=cremonalibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/libs/cremona/homspace.pyx",language="C++", define_macros=[("NTL_ALL",None)], libraries=cremonalibs)

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/rings/finite_field_givaro.pyx",language="C++", libraries=['givaro', 'gmpxx', 'gmp', 'm', 'stdc++'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/rings/finite_field_ntl_gf2e.pyx',language="C++", libraries=['ntl', 'gmp'])

    config_gcc_file(env,pyx_pextn_dict,"devel/sage/sage/rings/real_rqdf.pyx",language="C++", libraries=['qd', 'm', 'stdc++','gmp','mpfr','qd', 'csage' ])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/misc.pyx',libraries=['mpfr','gmp'])

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/matrix/matrix_integer_2x2.pyx',libraries=['gmp'])

    BLAS = env.options['BLAS']
    BLAS2 = env.options['BLAS2']

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/libs/linbox/linbox.pyx',language="C++", libraries=['linboxwrap', 'ntl', 'linbox', 'gmp', 'gmpxx', 'stdc++', 'givaro', BLAS, BLAS2])

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

    config_gcc_file(env,pyx_pextn_dict,'devel/sage/sage/interfaces/libecm.pyx',libraries = ['ecm', 'gmp'])

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

    tempdir = os.path.realpath(env.options['SAGE_ROOT'])

    gcceo_dict = { }
    gcceso_dict = { }

    mwrankcc = GCC_extension_object(gccc, env, ["devel/sage/sage/libs/mwrank/wrap.cc"], tempdir+"/devel/sage/build/temp/sage/libs/mwrank", define_macros=[("NTL_ALL",None)], options = { '-fPIC':None } )


    hypellfrob_cpp = GCC_extension_object(gccc, env,["devel/sage/sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.cpp"], tempdir+"sage/schemes/hyperelliptic_curves/hypellfrob/", language='C++', include_dirs=[abs_sage_path(env, 'sage/libs/ntl'), abs_sage_path(env,'sage/schemes/hyperelliptic_curves/hypellfrob/') ], options = { '-fPIC':None })
    recurrences_zn_poly = GCC_extension_object(gccc, env,[abs_sage_path(env,"sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.cpp")], tempdir+"sage/schemes/hyperelliptic_curves/hypellfrob/", language='C++', include_dirs=[abs_sage_path(env, 'sage/libs/ntl'), abs_sage_path(env,'sage/schemes/hyperelliptic_curves/hypellfrob/') ], options = { '-fPIC':None } )
    recurrences_ntl = GCC_extension_object(gccc, env,["devel/sage/sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.cpp"], tempdir+"sage/schemes/hyperelliptic_curves/hypellfrob/", language='C++', include_dirs=[abs_sage_path(env, 'sage/libs/ntl'), abs_sage_path(env,'sage/schemes/hyperelliptic_curves/hypellfrob/') ], options = { '-fPIC':None } )

    partitions_c = GCC_extension_object(gccc, env,["devel/sage/sage/combinat/partitions_c.cc"], tempdir+"/devel/sage/build/temp/sage/combinat", language='C++', options = { '-fPIC':None } )
    try:
        os.makedirs(tempdir+"/devel/sage/build/temp/sage/graphs/planarity")
    except:
        pass
    depfw = extfilewalker()
    depfw.addcallback('.pyx',lambda x: True)
    depfw.addcallback('.pxd',lambda x: True)
    depfw.addcallback('.pxi',lambda x: True)
    depfw.addcallback('.py',lambda x: x.find(".doctest") ==-1 )
    dep_list = depfw.walk('devel/sage/sage')
    funcdict = { '.pyx':get_cython_file_deps, '.pxd':get_cython_file_deps, 'pxi':get_cython_file_deps }
    includedict = { '.pyx':[abs_sage_path(env,'')], '.pxd':[abs_sage_path(env,'')], '.pxi':[abs_sage_path(env,'')] }
    global cacheinfo
    cacheinfo = get_compile_list(dep_list, funcdict, includedict, "sagebuild.chc", env.options['SAGE_ROOT'])
    TM.addfinishtask(sagebuild_exit)
    pyx_ext_dict  = { }
    for x in cacheinfo[0]:
        if os.path.splitext(x)[1]==".pyx":
            pyx_ext_dict[x]=pyx_pextn_dict[x]
    if verbose>100:
        print '---------------PYX_EXT_DICT------------'
        print pyx_ext_dict
        print '----------------------------------------'
    if pyx_ext_dict.has_key('devel/sage/sage/graphs/planarity.pyx'):
        efwplan = extfilewalker()
        efwplan.addcallback('.c',lambda x: True)
        planarity_list = efwplan.walk('devel/sage/sage/graphs/planarity')
        planarity_ext = list()
        for x in planarity_list:
            ext =  GCC_extension_object(gccc, env, [x], tempdir+"/devel/sage/build/temp/sage/graphs/planarity", options = { '-fPIC':None } )
            ext.generate_action(env).enable()
            planarity_ext.append(ext)
    if pyx_ext_dict.has_key("devel/sage/sage/libs/mwrank/mwrank.pyx"):
        mwrankcc.generate_action(env).enable()
    else:
        pass
    if pyx_ext_dict.has_key('devel/sage/sage/schemes/hyperelliptic_curves/hypellfrob.pyx'):
        try:
            os.makedirs(hypellfrob_cpp.get_out_dir())
        except:
            pass
        hypellfrob_cpp.generate_action(env).enable()
        recurrences_ntl.generate_action(env).enable()
        recurrences_zn_poly.generate_action(env).enable()
    if pyx_ext_dict.has_key('devel/sage/sage/combinat/partitions.pyx'):
        partitions_c.generate_action(env).enable()

    for x in pyx_ext_dict.values():
        q = x.generate_action(env)
        gcc_temp_dir = tempdir + '/devel/sage/build/temp/sage/' + (os.path.split(x.sources[0])[0])[len('devel/sage/sage/'):]
        gcceo = GCC_extension_object(gccc, env, [x], gcc_temp_dir, options = { '-fPIC':None } )
        gcceo_dict[x.sources[0]] = gcceo
        gccaction = gcceo.generate_action(env)
        gccaction.dep_register(q)
        q.enable()
    for gcceo in gcceo_dict.values():
        try:
            os.makedirs(gcceo.get_out_dir())
        except:
            pass
        gccaction = gcceo.generate_action(env)
        gccaction.enable()

    for gcceo in gcceo_dict.values():
        gcc_build_dir = develdir + '/build' + (os.path.split(gcceo.sources[0])[0])[len('devel/sage'):]
        gcceso = GCC_extension_shared_object(gccc,env, [gcceo], gcc_build_dir)
        gcceso_dict[gcceso.sources[0] ] = gcceso
        gcceso.generate_action(env).dep_register(gcceo.generate_action(env))
    dict_insert_src(gcceso_dict, env.options['SAGE_ROOT']+"/devel/sage/build/temp/sage/libs/mwrank/mwrank.o", mwrankcc.outfile)
    try:
        gcceso_dict[env.options['SAGE_ROOT'] + "/devel/sage/build/temp/sage/libs/mwrank/mwrank.o"].generate_action(env).dep_register(mwrankcc.generate_action(env))
    except KeyError:
        pass
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o", hypellfrob_cpp.outfile)
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o", recurrences_ntl.outfile)
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o", recurrences_zn_poly.outfile)
    try:
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o"].generate_action(env).dep_register(hypellfrob_cpp.generate_action(env))
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o"].generate_action(env).dep_register(recurrences_ntl.generate_action(env))
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/schemes/hyperelliptic_curves/hypellfrob.o"].generate_action(env).dep_register(recurrences_zn_poly.generate_action(env))
    except KeyError:
        pass
    dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/combinat/partitions.o", partitions_c.outfile)
    try:
        gcceso_dict[tempdir+"/devel/sage/build/temp/sage/combinat/partitions.o"].generate_action(env).dep_register(mwrankcc.generate_action(env))
    except KeyError:
        pass
    if pyx_ext_dict.has_key(abs_sage_path(env,'/devel/sage/graphs/planarity.pyx')):
        for x in planarity_ext:
            dict_insert_src(gcceso_dict, tempdir+"/devel/sage/build/temp/sage/graphs/planarity.o",x.outfile)


    for gcceso in gcceso_dict.values():
        try:
            os.makedirs(gcceso.get_out_dir())
        except:
            pass
        gccesoaction = gcceso.generate_action(env)
        gccesoaction.enable()


    if verbose>10:
        print 'Copying *.py files'
    for filenm in cacheinfo[0]:
        if os.path.splitext(filenm)[1]==".py":
            filepart = (os.path.split(filenm)[0])
            newdir = env.options['SAGE_ROOT'] + '/devel/sage/build/sage/' + filepart.replace('devel/sage/sage','')
            filename = os.path.split(filenm)[1]
            if filename[0]=='.':
                continue
            try:
                os.makedirs(newdir)
            except:
                pass
            cmd = 'cp %s %s/%s' % (filenm, newdir, filename)
            if verbose>30:
                print cmd
            os.system( cmd )

def sagebuild_exit():
    os.chdir(os.environ['SAGE_ROOT'])
    global cacheinfo
    if verbose>10:
        print 'Writing Sage Cache'
    write_cache(cacheinfo[0], cacheinfo[1], cacheinfo[2], cacheinfo[3], cacheinfo[4], cacheinfo[5], "sagebuild.chc")
