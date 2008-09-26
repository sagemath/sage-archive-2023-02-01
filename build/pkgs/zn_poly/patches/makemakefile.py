#
# This python script is called by configure to generate the makefile
# for zn_poly.
#

# --------------------------------------------------------------------------
# various lists of modules


# These are the modules that go into the actual zn_poly library. They get
# compiled in both optimised and debug modes.
lib_modules = ["array", "invert", "ks_support", "mulmid", "mulmid_ks", "misc",
               "mpn_mulmid", "mul", "mul_fft", "mul_fft_dft", "mul_ks",
               "nuss", "pack", "pmf", "pmfvec_fft", "tuning", "zn_mod"]
lib_modules = ["src/" + x for x in lib_modules]

# These are modules containing various test routines. They get compiled
# in debug mode only.
test_modules = ["test", "ref_mul", "invert-test", "pmfvec_fft-test",
                "mulmid_ks-test", "mpn_mulmid-test", "mul_fft-test",
                "mul_ks-test", "nuss-test", "pack-test"]
test_modules = ["test/" + x for x in test_modules]

# These are modules containing various profiling routines. They get compiled
# in optimised mode only.
prof_modules = ["prof_main", "profiler", "array-profile", "invert-profile",
                "mulmid-profile", "mpn_mulmid-profile", "mul-profile",
                "negamul-profile"]
prof_modules = ["profile/" + x for x in prof_modules]

# These are modules containing profiling routines for NTL. They get compiled
# in optimised mode only, with the C++ compiler.
cpp_prof_modules = ["profile/ntl-profile"]

# These are modules containing dummy routines replacing the NTL ones, if
# we're not compiling with NTL support.
noncpp_prof_modules = ["profile/ntl-profile-dummy"]

# These are modules shared by the test and profiling code. They get compiled
# in both debug and optimised mode.
testprof_modules = ["test/support"]

# Profiling targets. Each X has a main file "X-main.c" which is linked
# against prof_main.c. They are compiled once with PROFILE_NTL defined
# and once without.
prof_progs = ["array-profile", "invert-profile", "mulmid-profile",
              "mpn_mulmid-profile", "mul-profile", "negamul-profile"]
prof_progs = ["profile/" + x for x in prof_progs]

# These are modules used in the tuning program; they get compiled only in
# optimised mode.
tune_modules = ["tune", "mulmid-tune", "mpn_mulmid-tune", "mul-tune",
                "mul_ks-tune", "mulmid_ks-tune", "nuss-tune"]
tune_modules = ["tune/" + x for x in tune_modules]

# Demo programs.
demo_progs = ["demo/bernoulli/bernoulli"]

# These are the headers that need to be copied to the install location.
install_headers = ["include/zn_poly.h", "include/wide_arith.h"]

# These are the other headers.
other_headers = ["include/support.h", "include/profiler.h",
                 "include/zn_poly_internal.h"]


# --------------------------------------------------------------------------
# read command line options

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--prefix", dest="prefix", default="/usr/local")
parser.add_option("--cflags", dest="cflags", default="-O2")
parser.add_option("--ldflags", dest="ldflags", default="")
parser.add_option("--gmp-prefix", dest="gmp_prefix", default="/usr/local")
parser.add_option("--ntl-prefix", dest="ntl_prefix", default="/usr/local")
parser.add_option("--use-flint", dest="use_flint", action="store_true",
                  default=False)
parser.add_option("--flint-prefix", dest="flint_prefix", default="/usr/local")

options, args = parser.parse_args()


gmp_include_dir = options.gmp_prefix + "/include"
gmp_lib_dir = options.gmp_prefix + "/lib"

ntl_include_dir = options.ntl_prefix + "/include"
ntl_lib_dir = options.ntl_prefix + "/lib"

if options.use_flint:
   flint_include_dir = options.flint_prefix + "/include"
   flint_lib_dir = options.flint_prefix + "/lib"

cflags = options.cflags
ldflags = options.ldflags
prefix = options.prefix

includes = "-I" + gmp_include_dir + " -I./include"
libs = "-L" + gmp_lib_dir + " -lgmp -lm"

if options.use_flint:
   includes = includes + " -I" + flint_include_dir
   libs = libs + " -L" + flint_lib_dir + " -lflint"
   cflags = cflags + " -std=c99 -DZNP_USE_FLINT"

cpp_includes = includes + " -I" + ntl_include_dir
cpp_libs = libs + " -L" + ntl_lib_dir + " -lntl"


# --------------------------------------------------------------------------
# generate the makefile

import time

print "#"
print "# Do not edit directly -- this file was auto-generated"
print "# by makemakefile.py on " + time.strftime("%a, %d %b %Y %H:%M:%S +0000",
                                                 time.gmtime())
print "#"
print

print
print "CC = gcc"
print "CFLAGS = " + cflags
print "LDFLAGS = " + ldflags
print "INCLUDES = " + includes
print "LIBS = " + libs

print
print "CPP = g++"
print "CPPFLAGS = " + cflags
print "CPP_INCLUDES = " + cpp_includes
print "CPP_LIBS = " + cpp_libs

print
print "HEADERS = " + " ".join(install_headers + other_headers)
print "LIBOBJS = " + " ".join([x + ".o" for x in lib_modules])
print "TESTOBJS = " + " ".join([x + "-DEBUG.o" for x in \
   lib_modules + test_modules + testprof_modules])
print "PROFOBJS = " + " ".join([x + ".o" for x in \
   lib_modules + prof_modules + noncpp_prof_modules + testprof_modules])
print "CPP_PROFOBJS = " + " ".join([x + ".o" for x in \
   lib_modules + prof_modules + cpp_prof_modules + testprof_modules])
print "TUNEOBJS = " + " ".join([x + ".o" for x in \
   lib_modules + tune_modules + testprof_modules + prof_modules + \
   noncpp_prof_modules if x != "profile/prof_main"])

print "all: libzn_poly.a"
print
print "test: test/test"
print "tune: tune/tune"
print
print "check: test"
print "\ttest/test -quick all"
print
print "install:"
print "\tmkdir -p %s/include/zn_poly" % prefix
print "\tmkdir -p %s/lib" % prefix
print "\tcp libzn_poly.a %s/lib" % prefix
print "\tcp include/zn_poly.h %s/include/zn_poly" % prefix
print "\tcp include/wide_arith.h %s/include/zn_poly" % prefix
print
print "clean:"
print "\trm -f *.o"
print "\trm -f test/*.o"
print "\trm -f profile/*.o"
print "\trm -f tune/*.o"
print "\trm -f src/*.o"
print "\trm -f demo/bernoulli/*.o"
print "\trm -f libzn_poly.a"
print "\trm -f libzn_poly.dylib"
print "\trm -f libzn_poly*.so*"
print "\trm -f test/test"
print "\trm -f tune/tune"
for x in prof_progs:
   print "\trm -f " + x
   print "\trm -f " + x + "-ntl"
for x in demo_progs:
   print "\trm -f " + x


print
print
print "##### library targets"
print
print "libzn_poly.a: $(LIBOBJS)"
print "\tar -r libzn_poly.a $(LIBOBJS)"
print "\tranlib libzn_poly.a"
print
print "libzn_poly.dylib: $(LIBOBJS)"
print "\t$(CC) -single_module -fPIC -dynamiclib -o libzn_poly.dylib " \
      "$(LIBOBJS) $(LIBS)"
print
print "libzn_poly.dylib64: $(LIBOBJS)"
print "\t$(CC) -m64 -single_module -fPIC -dynamiclib -o libzn_poly.dylib $(LIBOBJS) $(LIBS)"
print
print "libzn_poly.so: $(LIBOBJS)"
print "\t$(CC) -shared -Wl,-soname,libzn_poly-`cat VERSION`.so " \
      "-o libzn_poly-`cat VERSION`.so $(LIBOBJS) $(LIBS)"
print "\t ln -sf libzn_poly-`cat VERSION`.so libzn_poly.so"

print
print
print "##### test program"
print
print "test/test: $(TESTOBJS) $(HEADERS)"
print "\t$(CC) -g $(LDFLAGS) -o test/test $(TESTOBJS) $(LIBS)"

print
print
print "##### profiling programs"
print
for x in prof_progs:
   print "%s-main.o: %s-main.c $(HEADERS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(INCLUDES) -DNDEBUG -o %s-main.o -c %s-main.c" \
         % (x, x)
   print
   print "%s: %s-main.o $(PROFOBJS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(LDFLAGS) -o %s %s-main.o $(PROFOBJS) $(LIBS)" \
         % (x, x)
   print
   print "%s-main-ntl.o: %s-main.c $(HEADERS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(INCLUDES) -DPROFILE_NTL -DNDEBUG " \
         "-o %s-main-ntl.o -c %s-main.c" % (x, x)
   print
   print "%s-ntl: %s-main-ntl.o $(CPP_PROFOBJS)" % (x, x)
   print "\t$(CPP) $(CPPFLAGS) $(LDFLAGS) -o %s-ntl %s-main-ntl.o " \
         "$(CPP_PROFOBJS) $(CPP_LIBS)" % (x, x)
   print

print
print
print "##### tuning utility"
print
print "tune/tune: $(TUNEOBJS)"
print "\t$(CC) $(CFLAGS) $(LDFLAGS) -o tune/tune $(TUNEOBJS) $(LIBS)"


print
print
print "##### demo programs"
for x in demo_progs:
   print
   print "%s: %s.o $(LIBOBJS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(LDFLAGS) -o %s %s.o $(LIBOBJS) $(LIBS)" % (x, x)


print
print
print "##### object files (with debug code)"
for x in lib_modules + test_modules + testprof_modules + demo_progs:
   print
   print "%s-DEBUG.o: %s.c $(HEADERS)" % (x, x)
   print "\t$(CC) -g $(CFLAGS) $(INCLUDES) -DDEBUG -o %s-DEBUG.o -c %s.c" \
         % (x, x)

print
print
print "##### object files (no debug code)"
for x in lib_modules + prof_modules + testprof_modules + \
                       tune_modules + demo_progs:
   print
   print "%s.o: %s.c $(HEADERS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(INCLUDES) -DNDEBUG -o %s.o -c %s.c" % (x, x)

print
print
print "##### object files (C++, no debug code)"
for x in cpp_prof_modules:
   print
   print "%s.o: %s.c $(HEADERS)" % (x, x)
   print "\t$(CPP) $(CPPFLAGS) $(CPP_INCLUDES) -DNDEBUG -o %s.o -c %s.c" \
         % (x, x)


### end of file
