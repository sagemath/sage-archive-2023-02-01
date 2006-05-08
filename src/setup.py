#!/usr/bin/env python
DEVEL = False

import distutils.sysconfig, os, sys
from distutils.core import setup, Extension

if len(sys.argv) > 1 and sys.argv[1] == "sdist":
    sdist = True
else:
    sdist = False

NO_WARN = True

if not os.environ.has_key('SAGE_ROOT'):
    print "    ERROR: The environment variable SAGE_ROOT must be defined."
    sys.exit(1)
else:
    SAGE_LOCAL = os.environ['SAGE_ROOT'] + '/local/'

if not os.environ.has_key('SAGE_VERSION'):
    SAGE_VERSION=0
else:
    SAGE_VERSION = os.environ['SAGE_VERSION']

SITE_PACKAGES = '%s/lib/python2.4/site-packages/'%SAGE_LOCAL

if not os.path.exists(SITE_PACKAGES):
    SITE_PACKAGES=None
    #raise RuntimeError, "Unable to find site-packages directory (see setup.py file in sage python code)."

ec =    Extension('sage.libs.ec.ec',
              sources = ["sage/libs/ec/ec.pyx"] +  \
                        ["sage/libs/ec/%s"%s for s in \
                         ["analrank.c", "apcompute.c", "arderivs.c",
                          "arintern.c", "arith.c", "artwists.c",
                          "arutil.c", "checkit.c", "degphi.c",
                          "diskio.c", "docurve.c", "dodisk.c",
                          "equation.c", "exotic.c", "fixit.c",
                          "iisog2.c", "iisog3.c", "isog.c", "isog2.c",
                          "isog23.c", "isog24.c", "isog3.c", "isog5.c",
                          "isog52.c", "isog713.c", "isogNprime.c",
                          "isoggen.c", "isogsort.c", "isogx0.c",
                          "isogx0branch.c", "isogx0branch1.c",
                          "isogx0branch2.c", "isogx0branch3.c",
                          "isogx0branch4.c", "isogx0branch5.c",
                          "isogx0branch6.c", "isogx0getd.c",
                          "isogx0period.c", "readit.c",
                          "special.c", "util.c"]],
              libraries = ["pari", "m"])

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

ntl = Extension('sage.libs.ntl.ntl',
                 sources = ["sage/libs/ntl/ntl.pyx", "sage/libs/ntl/wrap.cc"],
                 libraries = ["ntl", "gmp", "gmpxx", "m", "stdc++"]
                 )

mwrank =  Extension("sage.libs.mwrank.mwrank",
                    sources = ["sage/libs/mwrank/mwrank.pyx",
                         "sage/libs/mwrank/wrap.cc"],
                    define_macros = [("NTL_ALL",None)],
                    libraries = ["mwrank", "ntl", "gmp", "gmpxx", "stdc++", "m", "pari"])

pari = Extension('sage.libs.pari.gen',
                 sources = ["sage/libs/pari/gen.pyx"],
                 libraries = ['pari', 'gmp'])

cf = Extension('sage.libs.cf.cf',
               sources = ["sage/libs/cf/cf.pyxe", "sage/libs/cf/ftmpl_inst.cc"],
               libraries = ['cf', 'cfmem', 'gmp', 'stdc++', 'm']
               )

ext_modules = [ \
    ec, \

    pari, \

    mwrank, \

    ntl, \

    cf, \

    Extension('sage.ext.arith',
              sources = ['sage/ext/arith.pyx']), \

    Extension('sage.ext.arith_gmp',
              sources = ['sage/ext/arith_gmp.pyx'],
              libraries=['gmp']), \

    Extension('sage.ext.coerce',
              sources = ['sage/ext/coerce.pyx']), \

    Extension('sage.ext.congroup_pyx',
              sources = ['sage/ext/congroup_pyx.pyx', \
                         'sage/ext/arith.pyx']), \

    Extension('sage.ext.element',
              sources = ['sage/ext/element.pyx']), \

    Extension('sage.ext.module',
              sources = ['sage/ext/module.pyx']), \

    Extension('sage.ext.ring',
              sources = ['sage/ext/ring.pyx']), \

    Extension('sage.ext.group',
              sources = ['sage/ext/group.pyx']), \

    Extension('sage.ext.sage_object',
              sources = ['sage/ext/sage_object.pyx']), \

    Extension('sage.ext.gens',
              sources = ['sage/ext/gens.pyx']), \

    # this is twice, that correct?
    Extension('sage.ext.mpfr',
              sources = ['sage/ext/mpfr.pyx', 'sage/ext/ring.pyx'],
              libraries = ['mpfr', 'gmp']), \

    Extension('sage.ext.mpfr',
              sources = ['sage/ext/mpfr.pyx', 'sage/ext/ring.pyx'],
              libraries = ['mpfr', 'gmp']), \

    Extension('sage.ext.integer',
              sources = ['sage/ext/arith.pyx', 'sage/ext/integer.pyx', \
                         'sage/ext/mpn_pylong.c', 'sage/ext/mpz_pylong.c'],
              libraries=['gmp']), \

    Extension('sage.ext.intmod_pyx',
              sources = ['sage/ext/intmod_pyx.pyx']), \

    Extension('sage.ext.polynomial_pyx',
              sources = ['sage/ext/polynomial_pyx.pyx',
                         'sage/ext/arith_gmp.pyx'],
              libraries=['gmp']), \

    Extension('sage.ext.rational',
              sources = ['sage/ext/rational.pyx',
                         'sage/ext/arith.pyx', \
                         'sage/ext/integer.pyx', \
                         'sage/ext/mpn_pylong.c', 'sage/ext/mpz_pylong.c'],
              libraries=['gmp']), \

    Extension('sage.ext.sparse_poly',
              sources = ['sage/ext/sparse_poly.pyx'],
              libraries=['gmp']), \

    Extension('sage.ext.sparse_matrix_pyx',
              ['sage/ext/sparse_matrix_pyx.pyx',
               'sage/ext/integer.pyx',
               'sage/ext/rational.pyx',
               'sage/ext/arith.pyx',
               'sage/ext/mpn_pylong.c', 'sage/ext/mpz_pylong.c'],
              libraries=['gmp']), \

    Extension('sage.ext.dense_matrix_pyx',
              ['sage/ext/dense_matrix_pyx.pyx',
               'sage/ext/integer.pyx',
               'sage/ext/rational.pyx',
               'sage/ext/arith.pyx',
               'sage/ext/mpn_pylong.c', 'sage/ext/mpz_pylong.c'],
              libraries=['gmp']), \

    Extension('sage.ext.search',
              ['sage/ext/search.pyx']), \

    Extension('sage.ext.heilbronn',
              ['sage/ext/heilbronn.pyx',
               'sage/ext/p1list.pyx',
               'sage/ext/arith.pyx'],
              libraries = ['gmp', 'm']), \

    Extension('sage.ext.p1list',
              ['sage/ext/p1list.pyx',
               'sage/ext/arith.pyx'],
              libraries = ['gmp']), \
    ]


mpc = Extension('sage.ext.mpc',
              sources = ['sage/ext/mpc.pyx', 'sage/ext/ring.pyx'],
              libraries = ['mpc', 'mpfr', 'gmp'])


extra_compile_args = [ ]

if NO_WARN and distutils.sysconfig.get_config_var('CC').startswith("gcc"):
    extra_compile_args = ['-w']

if DEVEL:
    extra_compile_args.append('-ggdb')
    ext_modules.append(hanke)
    #ext_modules.append(mpc)

for m in ext_modules:
    m.sources += ['sage/ext/interrupt.c']

include_dirs = ['%s/include'%SAGE_LOCAL,
                '%s/include/python2.4'%SAGE_LOCAL]

extra_link_args =  ['-L%s/lib'%SAGE_LOCAL]

def need_to_create(file1, file2):
    """
    Return True if either file2 does not exist or is older than file1.
    """
    if not os.path.exists(file2):
        return True
    if os.path.getctime(file2) <= os.path.getctime(file1):
        return True
    return False


def process_pyrexembed_file(f):
    # This is a pyrexembed file, so process accordingly.
    dir, base = os.path.split(f[:-5])
    tmp = '%s/.tmp_pyrexembed'%dir
    if os.path.exists(tmp) and not os.path.isdir(tmp):
        print "Please delete file '%s' in %s"%(tmp, dir)
        sys.exit(1)
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    pyxe_file = "%s/%s.pyxe"%(tmp, base)

    # The following files will be produced by pyrexembed.
    cpp_file = "%s/%s_embed.cpp"%(dir, base)
    c_file = "%s/%s.c"%(dir, base)
    pyx_file = "%s/%s.pyx"%(dir,base)
    pyx_embed_file = "%s/%s.pyx"%(tmp, base)
    h_file = "%s/%s_embed.h"%(tmp, base)
    if need_to_create(f, pyx_file) or need_to_create(f, cpp_file) or need_to_create(f, h_file):
        os.system('cp -p %s %s'%(f, pyxe_file))
        os.system('cp -p %s/*.pxi %s'%(dir, tmp))
        os.system('cp -p %s/*.pxd %s'%(dir, tmp))
        os.system('cp -p %s/*.h %s'%(dir, tmp))
        cmd = "pyrexembed %s"%pyxe_file
        print cmd
        ret = os.system(cmd)
        if ret != 0:
            print "Error running pyrexembed."
            sys.exit(ret)
        process_pyrex_file(pyx_embed_file)
        cmd = 'cp -p %s/*.pyx %s/; cp -p %s/*.c %s/; cp -p %s/*.h %s/; cp -p %s/*.cpp %s/'%(tmp, dir, tmp, dir, tmp, dir, tmp, dir)
        print cmd
        os.system(cmd)
    return [cpp_file, c_file]

def process_pyrex_file(f):
    # This is a pyrex file, so process accordingly.
    if SITE_PACKAGES:
        pyx_inst_file = '%s/%s'%(SITE_PACKAGES, f)
        if need_to_create(f, pyx_inst_file):
            print "%s --> %s"%(f, pyx_inst_file)
            os.system('cp %s %s 2>/dev/null'%(f, pyx_inst_file))
    c_file = f[:-4] + ".c"
    if need_to_create(f, c_file):
        cmd = "pyrexc %s"%f
        print cmd
        ret = os.system(cmd)
        if ret != 0:
            print "Error running pyrexc."
            sys.exit(ret)
    return [c_file]


def pyrex(ext_modules):
    for m in ext_modules:
        m.extra_compile_args += extra_compile_args
        m.extra_link_args += extra_link_args
        new_sources = []
        for i in range(len(m.sources)):
            f = m.sources[i]
            s = open(f).read()
            if f[-5:] == '.pyxe':# and s.find("#embed") != -1 and s.find('#}embed') != -1:
                new_sources = process_pyrexembed_file(f)
            elif f[-4:] == ".pyx":
                new_sources += process_pyrex_file(f)
            else:
                new_sources.append(f)
        m.sources = new_sources


#############################################
# Update interrupt.h and interrupt.pxi files
for D in os.listdir("sage/libs/"):
    if os.path.isdir('sage/libs/%s'%D):
        os.system("cp sage/ext/interrupt.h sage/libs/%s/"%D)
        os.system("cp sage/ext/interrupt.pxi sage/libs/%s/"%D)


##########################################


if not sdist:
    pyrex(ext_modules)


setup(name        = 'sage',

      version     =  SAGE_VERSION,

      description = 'SAGE: System for Algebra and Geometry Experimentation',

      license     = 'GNU Public License (GPL)',

      author      = 'William Stein',

      author_email= 'wstein@math.ucsd.edu',

      url         = 'http://modular.ucsd.edu/sage',

      packages    = ['sage',

                     'sage.algebras',

                     'sage.categories',

                     'sage.coding',

                     'sage.crypto',

                     'sage.databases',

                     'sage.edu',

                     'sage.ext',

                     'sage.functions',

                     'sage.geometry',

                     'sage.groups',
                     'sage.groups.abelian_gps',
                     'sage.groups.matrix_gps',
                     'sage.groups.perm_gps',

                     'sage.interfaces',

                     'sage.lfunctions',

                     'sage.libs',
                     'sage.libs.hanke',
                     'sage.libs.mwrank',
                     'sage.libs.ntl',
                     'sage.libs.ec',
                     'sage.libs.pari',
                     'sage.libs.cf',

                     'sage.matrix',

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

                     'sage.rings',
                     'sage.rings.number_field',

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

                     'sage.structure',
                     ],

      ext_modules = ext_modules,
      include_dirs = include_dirs)



