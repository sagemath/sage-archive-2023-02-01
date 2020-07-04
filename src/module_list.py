import os
from distutils.extension import Extension
from sage.env import SAGE_LOCAL

SAGE_INC = os.path.join(SAGE_LOCAL, 'include')

#########################################################
### pkg-config setup
#########################################################

import pkgconfig

# CBLAS can be one of multiple implementations
cblas_pc = pkgconfig.parse('cblas')
cblas_libs = cblas_pc['libraries']
cblas_library_dirs = cblas_pc['library_dirs']
cblas_include_dirs = cblas_pc['include_dirs']

# TODO: Remove Cygwin hack by installing a suitable cblas.pc
if os.path.exists('/usr/lib/libblas.dll.a'):
    cblas_libs = ['gslcblas']

# LAPACK can be one of multiple implementations
lapack_pc = pkgconfig.parse('lapack')
lapack_libs = lapack_pc['libraries']
lapack_library_dirs = lapack_pc['library_dirs']
lapack_include_dirs = lapack_pc['include_dirs']

# GD image library
gd_pc = pkgconfig.parse('gdlib')
gd_libs = gd_pc['libraries']
gd_library_dirs = gd_pc['library_dirs']
gd_include_dirs = gd_pc['include_dirs']

# PNG image library
png_pc = pkgconfig.parse('libpng')
png_libs = png_pc['libraries']
png_library_dirs = png_pc['library_dirs']
png_include_dirs = png_pc['include_dirs']

# zlib
try:
    zlib_pc = pkgconfig.parse('zlib')
except pkgconfig.PackageNotFoundError:
    from collections import defaultdict
    zlib_pc = defaultdict(list, {'libraries': ['z']})
zlib_libs = zlib_pc['libraries']
zlib_library_dirs = zlib_pc['library_dirs']
zlib_include_dirs = zlib_pc['include_dirs']

#########################################################
### M4RI flags
#########################################################

m4ri_pc = pkgconfig.parse('m4ri')
m4ri_libs = m4ri_pc['libraries']
m4ri_library_dirs = m4ri_pc['library_dirs']
m4ri_include_dirs = m4ri_pc['include_dirs']

m4ri_extra_compile_args = pkgconfig.cflags('m4ri').split()
try:
    m4ri_extra_compile_args.remove("-pedantic")
except ValueError:
    pass

from sage.env import cython_aliases
aliases = cython_aliases()
arb_dylib_name = aliases["ARB_LIBRARY"]

#############################################################
### List of modules
###
### Note that the list of modules is sorted alphabetically
### by extension name. Please keep this list sorted when
### adding new modules!
###
#############################################################

from sage_setup.optional_extension import OptionalExtension

ext_modules = [

    ################################
    ##
    ## sage.algebras
    ##
    ################################

    Extension('*', ['sage/algebras/**/*.pyx']),

    ################################
    ##
    ## sage.arith
    ##
    ################################

    Extension('*', ['sage/arith/*.pyx']),

    ################################
    ##
    ## sage.calculus
    ##
    ################################

    Extension('*', ['sage/calculus/**/*.pyx']),

    ################################
    ##
    ## sage.categories
    ##
    ################################

    Extension('*', ['sage/categories/**/*.pyx']),

    ################################
    ##
    ## sage.coding
    ##
    ################################

    Extension('sage.coding.codecan.codecan',
              sources = ['sage/coding/codecan/codecan.pyx']),

    Extension('*', ['sage/coding/**/*.pyx']),

    ################################
    ##
    ## sage.combinat
    ##
    ################################

    Extension('*', ['sage/combinat/**/*.pyx']),

    Extension('sage.combinat.subword_complex_c',
              sources=['sage/combinat/subword_complex_c.pyx']),

    ################################
    ##
    ## sage.cpython
    ##
    ################################

    Extension('*', ['sage/cpython/*.pyx']),

    ################################
    ##
    ## sage.crypto
    ##
    ################################

    Extension('*', ['sage/crypto/*.pyx']),

    ################################
    ##
    ## sage.data_structures
    ##
    ################################

    Extension('*', ['sage/data_structures/*.pyx']),

    ################################
    ##
    ## sage.docs
    ##
    ################################

    Extension('*', ['sage/docs/*.pyx']),

    ################################
    ##
    ## sage.dynamics
    ##
    ################################

    Extension('sage.dynamics.arithmetic_dynamics.projective_ds_helper',
              sources = ['sage/dynamics/arithmetic_dynamics/projective_ds_helper.pyx']),

    Extension('sage.dynamics.complex_dynamics.mandel_julia_helper',
                sources = ['sage/dynamics/complex_dynamics/mandel_julia_helper.pyx']),

    ################################
    ##
    ## sage.ext
    ##
    ################################

    Extension('*', ['sage/ext/**/*.pyx']),

    ################################
    ##
    ## sage.finance
    ##
    ################################

    Extension('*', ['sage/finance/*.pyx']),

    ################################
    ##
    ## sage.functions
    ##
    ################################

    Extension('sage.functions.prime_pi',
        sources = ['sage/functions/prime_pi.pyx']),

    ################################
    ##
    ## sage.games
    ##
    ################################

    Extension('*', ['sage/games/*.pyx']),

    ################################
    ##
    ## sage.geometry
    ##
    ################################

    Extension('*', ['sage/geometry/**/*.pyx']),

    ################################
    ##
    ## sage.graphs
    ##
    ################################

    Extension('sage.graphs.asteroidal_triples',
              sources = ['sage/graphs/asteroidal_triples.pyx']),

    Extension('sage.graphs.chrompoly',
              sources = ['sage/graphs/chrompoly.pyx']),

    Extension('sage.graphs.cliquer',
              sources = ['sage/graphs/cliquer.pyx']),

    Extension('sage.graphs.centrality',
              sources = ['sage/graphs/centrality.pyx']),

    Extension('sage.graphs.independent_sets',
              sources = ['sage/graphs/independent_sets.pyx']),

    Extension('sage.graphs.graph_decompositions.fast_digraph',
              sources = ['sage/graphs/graph_decompositions/fast_digraph.pyx']),

    Extension('sage.graphs.graph_decompositions.vertex_separation',
              sources = ['sage/graphs/graph_decompositions/vertex_separation.pyx']),

    Extension('sage.graphs.graph_decompositions.graph_products',
              sources = ['sage/graphs/graph_decompositions/graph_products.pyx']),

    Extension('sage.graphs.convexity_properties',
              sources = ['sage/graphs/convexity_properties.pyx']),

    Extension('sage.graphs.comparability',
              sources = ['sage/graphs/comparability.pyx']),

    Extension('sage.graphs.generic_graph_pyx',
              sources = ['sage/graphs/generic_graph_pyx.pyx']),

    Extension('sage.graphs.traversals',
              sources = ['sage/graphs/traversals.pyx']),

    Extension('sage.graphs.graph_generators_pyx',
              sources = ['sage/graphs/graph_generators_pyx.pyx']),

    Extension('sage.graphs.distances_all_pairs',
              sources = ['sage/graphs/distances_all_pairs.pyx']),

    Extension('sage.graphs.base.graph_backends',
              sources = ['sage/graphs/base/graph_backends.pyx']),

    Extension('sage.graphs.base.static_dense_graph',
              sources = ['sage/graphs/base/static_dense_graph.pyx']),

    Extension('sage.graphs.base.static_sparse_graph',
              sources = ['sage/graphs/base/static_sparse_graph.pyx']),

    Extension('sage.graphs.base.static_sparse_backend',
              sources = ['sage/graphs/base/static_sparse_backend.pyx']),

    Extension('sage.graphs.graph_coloring',
              sources = ['sage/graphs/graph_coloring.pyx']),

    Extension('sage.graphs.line_graph',
              sources = ['sage/graphs/line_graph.pyx']),

    Extension('sage.graphs.weakly_chordal',
              sources = ['sage/graphs/weakly_chordal.pyx']),

    Extension('sage.graphs.matchpoly',
              sources = ['sage/graphs/matchpoly.pyx']),

    OptionalExtension("sage.graphs.mcqd",
              ["sage/graphs/mcqd.pyx"],
              package = 'mcqd'),

    OptionalExtension("sage.graphs.bliss",
              ["sage/graphs/bliss.pyx"],
              package = 'bliss'),

    Extension('sage.graphs.planarity',
              sources = ['sage/graphs/planarity.pyx']),

    Extension('sage.graphs.strongly_regular_db',
              sources = ['sage/graphs/strongly_regular_db.pyx']),

    Extension('sage.graphs.graph_decompositions.rankwidth',
              sources = ['sage/graphs/graph_decompositions/rankwidth.pyx'],
              libraries=['rw']),

    Extension('sage.graphs.graph_decompositions.bandwidth',
              sources = ['sage/graphs/graph_decompositions/bandwidth.pyx']),

    Extension('sage.graphs.graph_decompositions.cutwidth',
              sources = ['sage/graphs/graph_decompositions/cutwidth.pyx']),

    OptionalExtension('sage.graphs.graph_decompositions.tdlib',
              sources = ['sage/graphs/graph_decompositions/tdlib.pyx'],
              language="c++",
              package = 'tdlib'),

    Extension('sage.graphs.graph_decompositions.clique_separators',
              sources = ['sage/graphs/graph_decompositions/clique_separators.pyx']),

    Extension('sage.graphs.spanning_tree',
              sources = ['sage/graphs/spanning_tree.pyx']),

    Extension('sage.graphs.path_enumeration',
              sources = ['sage/graphs/path_enumeration.pyx']),

    Extension('sage.graphs.connectivity',
          sources = ['sage/graphs/connectivity.pyx']),

    Extension('sage.graphs.trees',
              sources = ['sage/graphs/trees.pyx']),

    Extension('sage.graphs.genus',
              sources = ['sage/graphs/genus.pyx']),

    Extension('sage.graphs.hyperbolicity',
              sources = ['sage/graphs/hyperbolicity.pyx']),

    Extension('sage.graphs.base.c_graph',
              sources = ['sage/graphs/base/c_graph.pyx']),

    Extension('sage.graphs.base.sparse_graph',
              sources = ['sage/graphs/base/sparse_graph.pyx']),

    Extension('sage.graphs.base.dense_graph',
              sources = ['sage/graphs/base/dense_graph.pyx']),

    Extension('sage.graphs.base.boost_graph',
              sources = ['sage/graphs/base/boost_graph.pyx']),

    Extension('sage.graphs.views',
              sources = ['sage/graphs/views.pyx']),

    ################################
    ##
    ## sage.groups
    ##
    ################################

    Extension('*', ['sage/groups/**/*.pyx']),

    ################################
    ##
    ## sage.interacts
    ##
    ################################

    Extension('*', ['sage/interacts/*.pyx']),

    ################################
    ##
    ## sage.interfaces
    ##
    ################################

    OptionalExtension("sage.interfaces.primecount",
              ["sage/interfaces/primecount.pyx"],
              package = "primecount"),

    Extension('*', ['sage/interfaces/*.pyx']),

    ################################
    ##
    ## sage.lfunctions
    ##
    ################################

    Extension('sage.lfunctions.zero_sums',
              sources = ['sage/lfunctions/zero_sums.pyx']),

    ################################
    ##
    ## sage.libs
    ##
    ################################

    OptionalExtension('sage.libs.coxeter3.coxeter',
              sources = ['sage/libs/coxeter3/coxeter.pyx'],
              package = 'coxeter3'),

    Extension('sage.libs.ecl',
              sources = ["sage/libs/ecl.pyx"]),

    OptionalExtension("sage.libs.fes",
             ["sage/libs/fes.pyx"],
             package = 'fes'),

    Extension('sage.libs.flint.flint',
              sources = ["sage/libs/flint/flint.pyx"]),

    Extension('sage.libs.flint.fmpz_poly',
              sources = ["sage/libs/flint/fmpz_poly.pyx"]),

    Extension('sage.libs.flint.arith',
              sources = ["sage/libs/flint/arith.pyx"]),

    Extension("sage.libs.glpk.error",
             ["sage/libs/glpk/error.pyx"]),

    Extension('sage.libs.gmp.pylong',
              sources = ['sage/libs/gmp/pylong.pyx']),

    Extension('sage.libs.braiding',
                      sources = ["sage/libs/braiding.pyx"]),

    Extension('sage.libs.homfly',
                      sources = ["sage/libs/homfly.pyx"]),

    OptionalExtension('sage.libs.sirocco',
                      sources = ["sage/libs/sirocco.pyx"],
                      package="sirocco"),

    Extension('*', ['sage/libs/linbox/*.pyx']),

    Extension('sage.libs.lcalc.lcalc_Lfunction',
              sources = ['sage/libs/lcalc/lcalc_Lfunction.pyx']),

    Extension('sage.libs.libecm',
              sources = ['sage/libs/libecm.pyx']),

    Extension('sage.libs.lrcalc.lrcalc',
              sources = ["sage/libs/lrcalc/lrcalc.pyx"]),

    OptionalExtension("sage.libs.meataxe",
              sources = ['sage/libs/meataxe.pyx'],
              package = 'meataxe'),

    Extension('*', ['sage/libs/pari/*.pyx']),

    Extension('sage.libs.ppl',
              sources = ['sage/libs/ppl.pyx']),

    Extension('*', ['sage/libs/pynac/*.pyx']),

    Extension('sage.libs.ratpoints',
              sources = ["sage/libs/ratpoints.pyx"]),

    Extension('sage.libs.readline',
              sources = ['sage/libs/readline.pyx']),

    Extension('*', sources = ['sage/libs/singular/*.pyx']),

    Extension('sage.libs.symmetrica.symmetrica',
              sources = ["sage/libs/symmetrica/symmetrica.pyx"]),

    Extension('sage.libs.mpmath.utils',
              sources = ["sage/libs/mpmath/utils.pyx"]),

    Extension('sage.libs.mpmath.ext_impl',
              sources = ["sage/libs/mpmath/ext_impl.pyx"]),

    Extension('sage.libs.mpmath.ext_main',
              sources = ["sage/libs/mpmath/ext_main.pyx"]),

    Extension('sage.libs.mpmath.ext_libmp',
              sources = ["sage/libs/mpmath/ext_libmp.pyx"]),

    ###################################
    ##
    ## sage.libs.arb
    ##
    ###################################

    Extension('*', ["sage/libs/arb/*.pyx"]),

    ###################################
    ##
    ## sage.libs.eclib
    ##
    ###################################

    Extension('*', ["sage/libs/eclib/*.pyx"]),

    ################################
    ##
    ## sage.libs.gap
    ##
    ################################

    Extension('*', ["sage/libs/gap/*.pyx"]),

    ###################################
    ##
    ## sage.libs.gsl
    ##
    ###################################

    Extension('*', ["sage/libs/gsl/*.pyx"]),

    ###################################
    ##
    ## sage.libs.ntl
    ##
    ###################################

    Extension('*', ["sage/libs/ntl/*.pyx"]),

    ################################
    ##
    ## sage.matrix
    ##
    ################################

    Extension('sage.matrix.action',
              sources = ['sage/matrix/action.pyx']),

    Extension('sage.matrix.args',
              sources = ['sage/matrix/args.pyx']),

    Extension('sage.matrix.echelon_matrix',
              sources = ['sage/matrix/echelon_matrix.pyx']),

    Extension('sage.matrix.change_ring',
              sources = ['sage/matrix/change_ring.pyx']),

    Extension('sage.matrix.constructor',
              sources = ['sage/matrix/constructor.pyx']),

    Extension('sage.matrix.matrix',
              sources = ['sage/matrix/matrix.pyx']),

    Extension('sage.matrix.matrix0',
              sources = ['sage/matrix/matrix0.pyx']),

    Extension('sage.matrix.matrix1',
              sources = ['sage/matrix/matrix1.pyx']),

    Extension('sage.matrix.matrix2',
              sources = ['sage/matrix/matrix2.pyx']),

    Extension("sage.matrix.matrix_complex_ball_dense",
              ["sage/matrix/matrix_complex_ball_dense.pyx"]),

    Extension('sage.matrix.matrix_complex_double_dense',
              sources = ['sage/matrix/matrix_complex_double_dense.pyx']),

    Extension('sage.matrix.matrix_cyclo_dense',
              sources = ['sage/matrix/matrix_cyclo_dense.pyx']),

    Extension('sage.matrix.matrix_gap',
              sources = ['sage/matrix/matrix_gap.pyx']),

    Extension('sage.matrix.matrix_dense',
              sources = ['sage/matrix/matrix_dense.pyx']),

    Extension('sage.matrix.matrix_double_dense',
              sources = ['sage/matrix/matrix_double_dense.pyx']),

    Extension('sage.matrix.matrix_generic_dense',
              sources = ['sage/matrix/matrix_generic_dense.pyx']),

    Extension('sage.matrix.matrix_generic_sparse',
              sources = ['sage/matrix/matrix_generic_sparse.pyx']),

    Extension('sage.matrix.matrix_integer_dense',
              sources = ['sage/matrix/matrix_integer_dense.pyx']),

    Extension('sage.matrix.matrix_integer_sparse',
              sources = ['sage/matrix/matrix_integer_sparse.pyx']),

    Extension('sage.matrix.matrix_mod2_dense',
              sources = ['sage/matrix/matrix_mod2_dense.pyx']),

    Extension('sage.matrix.matrix_gf2e_dense',
              sources = ['sage/matrix/matrix_gf2e_dense.pyx']),

    Extension('sage.matrix.matrix_modn_dense_float',
              sources = ['sage/matrix/matrix_modn_dense_float.pyx']),

    Extension('sage.matrix.matrix_modn_dense_double',
              sources = ['sage/matrix/matrix_modn_dense_double.pyx']),

    Extension('sage.matrix.matrix_modn_sparse',
              sources = ['sage/matrix/matrix_modn_sparse.pyx']),

    Extension('sage.matrix.matrix_mpolynomial_dense',
              sources = ['sage/matrix/matrix_mpolynomial_dense.pyx']),

    Extension('sage.matrix.matrix_polynomial_dense',
              sources = ['sage/matrix/matrix_polynomial_dense.pyx']),

    Extension('sage.matrix.matrix_rational_dense',
              sources = ['sage/matrix/matrix_rational_dense.pyx']),

    Extension('sage.matrix.matrix_rational_sparse',
              sources = ['sage/matrix/matrix_rational_sparse.pyx']),

    Extension('sage.matrix.matrix_real_double_dense',
              sources = ['sage/matrix/matrix_real_double_dense.pyx']),

    Extension('sage.matrix.matrix_sparse',
              sources = ['sage/matrix/matrix_sparse.pyx']),

    Extension('sage.matrix.matrix_symbolic_dense',
              sources = ['sage/matrix/matrix_symbolic_dense.pyx']),

    Extension('sage.matrix.matrix_window',
              sources = ['sage/matrix/matrix_window.pyx']),

    OptionalExtension("sage.matrix.matrix_gfpn_dense",
              sources = ['sage/matrix/matrix_gfpn_dense.pyx'],
              package = 'meataxe'),

    Extension('sage.matrix.misc',
              sources = ['sage/matrix/misc.pyx']),

    Extension('sage.matrix.strassen',
              sources = ['sage/matrix/strassen.pyx']),

    ################################
    ##
    ## sage.matroids
    ##
    ################################

    Extension('*', ['sage/matroids/*.pyx']),

    ################################
    ##
    ## sage.media
    ##
    ################################

    Extension('*', ['sage/media/*.pyx']),

    ################################
    ##
    ## sage.misc
    ##
    ################################

    Extension('*', ['sage/misc/*.pyx']),

    ################################
    ##
    ## sage.modular
    ##
    ################################

    Extension('*', ['sage/modular/**/*.pyx']),

    ################################
    ##
    ## sage.modules
    ##
    ################################

    Extension('*', ['sage/modules/**/*.pyx']),

    ################################
    ##
    ## sage.numerical
    ##
    ################################


    Extension("sage.numerical.mip",
              ["sage/numerical/mip.pyx"]),

    Extension("sage.numerical.linear_functions",
              ["sage/numerical/linear_functions.pyx"]),

    Extension("sage.numerical.linear_tensor_element",
              ["sage/numerical/linear_tensor_element.pyx"]),

    Extension("sage.numerical.gauss_legendre",
              ["sage/numerical/gauss_legendre.pyx"]),

    Extension("sage.numerical.sdp",
              ["sage/numerical/sdp.pyx"]),

    Extension("sage.numerical.backends.generic_backend",
              ["sage/numerical/backends/generic_backend.pyx"]),

    Extension("sage.numerical.backends.generic_sdp_backend",
              ["sage/numerical/backends/generic_sdp_backend.pyx"]),

    Extension("sage.numerical.backends.glpk_backend",
              ["sage/numerical/backends/glpk_backend.pyx"]),

    Extension("sage.numerical.backends.glpk_exact_backend",
              ["sage/numerical/backends/glpk_exact_backend.pyx"]),

    Extension("sage.numerical.backends.ppl_backend",
              ["sage/numerical/backends/ppl_backend.pyx"]),

    Extension("sage.numerical.backends.cvxopt_backend",
              ["sage/numerical/backends/cvxopt_backend.pyx"]),

    Extension("sage.numerical.backends.cvxopt_sdp_backend",
              ["sage/numerical/backends/cvxopt_sdp_backend.pyx"]),

    Extension("sage.numerical.backends.glpk_graph_backend",
              ["sage/numerical/backends/glpk_graph_backend.pyx"]),

    Extension("sage.numerical.backends.interactivelp_backend",
              ["sage/numerical/backends/interactivelp_backend.pyx"]),

    ################################
    ##
    ## sage.plot
    ##
    ################################

    Extension('*', ['sage/plot/**/*.pyx']),

    ################################
    ##
    ## sage.probability
    ##
    ################################

    Extension('*', ['sage/probability/*.pyx']),

    ################################
    ##
    ## sage.quadratic_forms
    ##
    ################################

    Extension('*', ['sage/quadratic_forms/*.pyx']),

    ###############################
    ##
    ## sage.quivers
    ##
    ###############################

    Extension('*', ['sage/quivers/*.pyx']),

    ################################
    ##
    ## sage.rings
    ##
    ################################

    Extension('sage.rings.sum_of_squares',
              sources = ['sage/rings/sum_of_squares.pyx']),

    Extension('sage.rings.bernmm',
              sources = ['sage/rings/bernmm.pyx']),

    Extension('sage.rings.bernoulli_mod_p',
              sources = ['sage/rings/bernoulli_mod_p.pyx']),

    Extension("sage.rings.complex_arb",
              ["sage/rings/complex_arb.pyx"]),

    Extension('sage.rings.complex_double',
              sources = ['sage/rings/complex_double.pyx']),

    Extension('sage.rings.complex_interval',
              sources = ['sage/rings/complex_interval.pyx']),

    Extension('sage.rings.complex_number',
              sources = ['sage/rings/complex_number.pyx']),

    Extension('sage.rings.integer',
              sources = ['sage/rings/integer.pyx']),

    Extension('sage.rings.integer_ring',
              sources = ['sage/rings/integer_ring.pyx']),

    Extension('sage.rings.factorint',
              sources = ['sage/rings/factorint.pyx']),

    Extension('sage.rings.fast_arith',
              sources = ['sage/rings/fast_arith.pyx']),

    Extension('sage.rings.fraction_field_element',
              sources = ['sage/rings/fraction_field_element.pyx']),

    Extension('sage.rings.fraction_field_FpT',
              sources = ['sage/rings/fraction_field_FpT.pyx']),

    Extension('sage.rings.laurent_series_ring_element',
              sources = ['sage/rings/laurent_series_ring_element.pyx']),

    Extension('sage.rings.morphism',
              sources = ['sage/rings/morphism.pyx']),

    Extension('sage.rings.complex_mpc',
              sources = ['sage/rings/complex_mpc.pyx']),

    Extension('sage.rings.noncommutative_ideals',
              sources = ['sage/rings/noncommutative_ideals.pyx']),

    Extension('sage.rings.power_series_mpoly',
              sources = ['sage/rings/power_series_mpoly.pyx']),

    Extension('sage.rings.power_series_poly',
              sources = ['sage/rings/power_series_poly.pyx']),

    Extension('sage.rings.power_series_pari',
              sources = ['sage/rings/power_series_pari.pyx']),

    Extension('sage.rings.power_series_ring_element',
              sources = ['sage/rings/power_series_ring_element.pyx']),

    Extension('sage.rings.tate_algebra_element',
              sources = ['sage/rings/tate_algebra_element.pyx']),

    Extension('sage.rings.tate_algebra_ideal',
              sources = ['sage/rings/tate_algebra_ideal.pyx']),

    Extension('sage.rings.puiseux_series_ring_element',
              sources = ['sage/rings/puiseux_series_ring_element.pyx']),

    Extension('sage.rings.rational',
              sources = ['sage/rings/rational.pyx']),

    Extension('sage.rings.real_double',
              sources = ['sage/rings/real_double.pyx']),

    Extension('sage.rings.real_interval_absolute',
              sources = ['sage/rings/real_interval_absolute.pyx']),

    Extension("sage.rings.real_arb",
              ["sage/rings/real_arb.pyx"]),

    Extension('sage.rings.real_lazy',
              sources = ['sage/rings/real_lazy.pyx']),

    Extension('sage.rings.real_mpfi',
              sources = ['sage/rings/real_mpfi.pyx']),

    Extension('sage.rings.real_mpfr',
              sources = ['sage/rings/real_mpfr.pyx']),

    Extension('sage.rings.finite_rings.residue_field',
              sources = ['sage/rings/finite_rings/residue_field.pyx']),

    Extension('sage.rings.ring',
              sources = ['sage/rings/ring.pyx']),

    Extension('sage.rings.ring_extension',
              sources = ['sage/rings/ring_extension.pyx']),

    Extension('sage.rings.ring_extension_element',
              sources = ['sage/rings/ring_extension_element.pyx']),

    Extension('sage.rings.ring_extension_morphism',
              sources = ['sage/rings/ring_extension_morphism.pyx']),

    Extension('sage.rings.ring_extension_conversion',
              sources = ['sage/rings/ring_extension_conversion.pyx']),

    Extension('*', ['sage/rings/convert/*.pyx']),

    ################################
    ##
    ## sage.rings.finite_rings
    ##
    ################################

    Extension('sage.rings.finite_rings.finite_field_base',
              sources = ['sage/rings/finite_rings/finite_field_base.pyx']),

    Extension('sage.rings.finite_rings.element_base',
              sources = ['sage/rings/finite_rings/element_base.pyx']),

    Extension('sage.rings.finite_rings.integer_mod',
              sources = ['sage/rings/finite_rings/integer_mod.pyx']),

    Extension('sage.rings.finite_rings.element_givaro',
              sources = ["sage/rings/finite_rings/element_givaro.pyx"]),

    Extension('sage.rings.finite_rings.element_ntl_gf2e',
              sources = ['sage/rings/finite_rings/element_ntl_gf2e.pyx']),

    Extension('sage.rings.finite_rings.element_pari_ffelt',
              sources = ['sage/rings/finite_rings/element_pari_ffelt.pyx']),

    Extension('sage.rings.finite_rings.hom_finite_field',
              sources = ["sage/rings/finite_rings/hom_finite_field.pyx"]),

    Extension('sage.rings.finite_rings.hom_prime_finite_field',
              sources = ["sage/rings/finite_rings/hom_prime_finite_field.pyx"]),

    Extension('sage.rings.finite_rings.hom_finite_field_givaro',
              sources = ["sage/rings/finite_rings/hom_finite_field_givaro.pyx"]),

    ################################
    ##
    ## sage.rings.function_field
    ##
    ################################

    Extension('sage.rings.function_field.element',
              sources = ['sage/rings/function_field/element.pyx']),

    ################################
    ##
    ## sage.rings.number_field
    ##
    ################################

    Extension('sage.rings.number_field.number_field_base',
              sources = ['sage/rings/number_field/number_field_base.pyx']),

    Extension('sage.rings.number_field.number_field_element',
              sources = ['sage/rings/number_field/number_field_element.pyx']),

    Extension('sage.rings.number_field.number_field_element_quadratic',
              sources = ['sage/rings/number_field/number_field_element_quadratic.pyx']),

    Extension('sage.rings.number_field.number_field_morphisms',
              sources = ['sage/rings/number_field/number_field_morphisms.pyx']),

    Extension('sage.rings.number_field.totallyreal',
              sources = ['sage/rings/number_field/totallyreal.pyx']),

    Extension('sage.rings.number_field.totallyreal_data',
              sources = ['sage/rings/number_field/totallyreal_data.pyx']),

    ################################
    ##
    ## sage.rings.padics
    ##
    ################################

    Extension('sage.rings.padics.morphism',
              sources = ['sage/rings/padics/morphism.pyx']),

    Extension('sage.rings.padics.common_conversion',
              sources = ['sage/rings/padics/common_conversion.pyx']),

    Extension('sage.rings.padics.local_generic_element',
              sources = ['sage/rings/padics/local_generic_element.pyx']),

    Extension('sage.rings.padics.padic_capped_absolute_element',
              sources = ['sage/rings/padics/padic_capped_absolute_element.pyx']),

    Extension('sage.rings.padics.padic_capped_relative_element',
              sources = ['sage/rings/padics/padic_capped_relative_element.pyx']),

    Extension('sage.rings.padics.padic_floating_point_element',
              sources = ['sage/rings/padics/padic_floating_point_element.pyx']),

    Extension('sage.rings.padics.padic_ext_element',
              sources = ['sage/rings/padics/padic_ext_element.pyx']),

    Extension('sage.rings.padics.padic_fixed_mod_element',
              sources = ['sage/rings/padics/padic_fixed_mod_element.pyx']),

    Extension('sage.rings.padics.padic_generic_element',
              sources = ['sage/rings/padics/padic_generic_element.pyx']),

    Extension('sage.rings.padics.padic_printing',
              sources = ['sage/rings/padics/padic_printing.pyx']),

    Extension('sage.rings.padics.padic_ZZ_pX_CA_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CA_element.pyx']),

    Extension('sage.rings.padics.padic_ZZ_pX_CR_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CR_element.pyx']),

    Extension('sage.rings.padics.padic_ZZ_pX_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_element.pyx']),

    Extension('sage.rings.padics.padic_ZZ_pX_FM_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_FM_element.pyx']),

    Extension('sage.rings.padics.pow_computer',
              sources = ['sage/rings/padics/pow_computer.pyx']),

    Extension('sage.rings.padics.pow_computer_ext',
              sources = ['sage/rings/padics/pow_computer_ext.pyx']),

    Extension('sage.rings.padics.pow_computer_flint',
              sources = ['sage/rings/padics/pow_computer_flint.pyx']),

    Extension('sage.rings.padics.pow_computer_relative',
              sources = ['sage/rings/padics/pow_computer_relative.pyx']),

    Extension('sage.rings.padics.qadic_flint_CR',
              sources = ['sage/rings/padics/qadic_flint_CR.pyx']),

    Extension('sage.rings.padics.qadic_flint_CA',
              sources = ['sage/rings/padics/qadic_flint_CA.pyx']),

    Extension('sage.rings.padics.qadic_flint_FM',
              sources = ['sage/rings/padics/qadic_flint_FM.pyx']),

    Extension('sage.rings.padics.qadic_flint_FP',
              sources = ['sage/rings/padics/qadic_flint_FP.pyx']),

    Extension('sage.rings.padics.relative_ramified_FM',
              sources = ['sage/rings/padics/relative_ramified_FM.pyx']),
    Extension('sage.rings.padics.relative_ramified_CA',
              sources = ['sage/rings/padics/relative_ramified_CA.pyx']),
    Extension('sage.rings.padics.relative_ramified_CR',
              sources = ['sage/rings/padics/relative_ramified_CR.pyx']),
    Extension('sage.rings.padics.relative_ramified_FP',
              sources = ['sage/rings/padics/relative_ramified_FP.pyx']),

    ################################
    ##
    ## sage.rings.polynomial
    ##
    ################################

    Extension('sage.rings.polynomial.cyclotomic',
              sources = ['sage/rings/polynomial/cyclotomic.pyx']),

    Extension('sage.rings.polynomial.evaluation',
              sources = ['sage/rings/polynomial/evaluation.pyx']),

    Extension('sage.rings.polynomial.laurent_polynomial',
              sources = ['sage/rings/polynomial/laurent_polynomial.pyx']),

    Extension('sage.rings.polynomial.hilbert',
              sources = ['sage/rings/polynomial/hilbert.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial',
              sources = ['sage/rings/polynomial/multi_polynomial.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial_ideal_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_ideal_libsingular.pyx']),

    Extension('sage.rings.polynomial.plural',
              sources = ['sage/rings/polynomial/plural.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_libsingular.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial_ring_base',
              sources = ['sage/rings/polynomial/multi_polynomial_ring_base.pyx']),

    Extension('sage.rings.polynomial.polynomial_number_field',
              sources = ['sage/rings/polynomial/polynomial_number_field.pyx']),

    Extension('sage.rings.polynomial.polydict',
              sources = ['sage/rings/polynomial/polydict.pyx']),

    Extension('sage.rings.polynomial.polynomial_complex_arb',
               sources = ['sage/rings/polynomial/polynomial_complex_arb.pyx']),

    Extension('sage.rings.polynomial.polynomial_compiled',
               sources = ['sage/rings/polynomial/polynomial_compiled.pyx']),

    Extension('sage.rings.polynomial.polynomial_element',
              sources = ['sage/rings/polynomial/polynomial_element.pyx']),

    Extension('sage.rings.polynomial.polynomial_gf2x',
              sources = ['sage/rings/polynomial/polynomial_gf2x.pyx']),

    Extension('sage.rings.polynomial.polynomial_zz_pex',
              sources = ['sage/rings/polynomial/polynomial_zz_pex.pyx']),

    Extension('sage.rings.polynomial.polynomial_zmod_flint',
              sources = ['sage/rings/polynomial/polynomial_zmod_flint.pyx']),

    Extension('sage.rings.polynomial.polynomial_integer_dense_flint',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_flint.pyx']),

    Extension('sage.rings.polynomial.polynomial_integer_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_ntl.pyx']),

    Extension('sage.rings.polynomial.polynomial_rational_flint',
              sources = ['sage/rings/polynomial/polynomial_rational_flint.pyx']),

    Extension('sage.rings.polynomial.polynomial_modn_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_modn_dense_ntl.pyx']),

    Extension('sage.rings.polynomial.polynomial_ring_homomorphism',
              sources = ['sage/rings/polynomial/polynomial_ring_homomorphism.pyx']),

    Extension('sage.rings.polynomial.pbori',
              sources = ['sage/rings/polynomial/pbori.pyx'],
              depends = [SAGE_INC + "/polybori/" + hd + ".h" for hd in ["polybori", "config"]]),

    Extension('sage.rings.polynomial.polynomial_real_mpfr_dense',
              sources = ['sage/rings/polynomial/polynomial_real_mpfr_dense.pyx']),

    Extension('sage.rings.polynomial.real_roots',
              sources = ['sage/rings/polynomial/real_roots.pyx']),

    Extension('sage.rings.polynomial.refine_root',
              sources = ['sage/rings/polynomial/refine_root.pyx']),

    Extension('sage.rings.polynomial.symmetric_reduction',
              sources = ['sage/rings/polynomial/symmetric_reduction.pyx']),

    Extension('sage.rings.polynomial.skew_polynomial_element',
              sources = ['sage/rings/polynomial/skew_polynomial_element.pyx']),
    
    Extension('sage.rings.polynomial.skew_polynomial_finite_order',
              sources = ['sage/rings/polynomial/skew_polynomial_finite_order.pyx']),

    Extension('sage.rings.polynomial.skew_polynomial_finite_field',
              sources = ['sage/rings/polynomial/skew_polynomial_finite_field.pyx']),

    # Note that weil_polynomials includes distutils directives in order to support
    # conditional OpenMP compilation (by uncommenting lines)
    Extension('sage.rings.polynomial.weil.weil_polynomials',
              sources = ['sage/rings/polynomial/weil/weil_polynomials.pyx']),


    ################################
    ##
    ## sage.rings.semirings
    ##
    ################################

    Extension('sage.rings.semirings.tropical_semiring',
              sources = ['sage/rings/semirings/tropical_semiring.pyx']),

    ################################
    ##
    ## sage.sat
    ##
    ################################

    Extension('sage.sat.solvers.satsolver',
              sources = ['sage/sat/solvers/satsolver.pyx']),

    ################################
    ##
    ## sage.schemes
    ##
    ################################

    Extension('*', ['sage/schemes/**/*.pyx']),

    ################################
    ##
    ## sage.sets
    ##
    ################################

    Extension('*', ['sage/sets/*.pyx']),

    ################################
    ##
    ## sage.stats
    ##
    ################################

    Extension('*', ['sage/stats/**/*.pyx']),

    ################################
    ##
    ## sage.structure
    ##
    ################################

    Extension('*', ['sage/structure/*.pyx']),

    ################################
    ##
    ## sage.symbolic
    ##
    ################################

    Extension('*', ['sage/symbolic/*.pyx']),

    ################################
    ##
    ## sage.tests
    ##
    ################################

    Extension('*', ['sage/tests/**/*.pyx'])

]
