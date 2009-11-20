#!/bin/sh
# Determine the type of C compiler, which can later
# be used to determine the flags the compiler
# will want. Do this by defining pre-processing a bit
# of C++ code, and checking what are defined.

# The Sun and GNU compilers have been tested.

# The HP and GNU compilers have not been tested, but
# use information gained from the documentation.

# HP-UX C and C++ compiler.
# http://docs.hp.com/en/7730/newhelp0610/preprocess.htm

# IBM Compiler Reference - XL C/C++ for AIX, V10.1
# http://www-01.ibm.com/support/docview.wss?uid=swg27012860&aid=1

# Using HP C++ for Tru64 UNIX and Linux Alpha
# http://h30097.www3.hp.com/cplus/ugu_impl.html#implem_chap


# First, make sure the enviroment variable CC is defined.

if [ -z "$CC" ]; then
   echo "Sorry, you should define the enivronment variable CC"
   exit 1
fi

# Create a test file. It does not need to be a complete
# C++ file, as it is only pre-processed. So there is no
# need for a 'main'

TESTFILE=/tmp/test.$$.c

# The flags for the GNU compilers do not change with
# operating system, so there is no need to worry too
# much about what system this is on.

# However, many compilers will make themselves appears to be
# GNU compilers, such as those from Intel and Sun, when
# they are really improved versions of the basic GCC.
# So care is needed to identify exactly what compiler it is.

echo "#ifdef __GNUC__"                   > $TESTFILE
echo "#ifdef __INTEL_COMPILER"           >> $TESTFILE
echo "it_is_the_Intel_improved_GCC"      >> $TESTFILE
echo "#elseif __SUNPRO_C        "        >> $TESTFILE
echo "it_is_the_Sun_improved_GCC"        >> $TESTFILE
echo "#else                  "           >> $TESTFILE
echo "it_is_the_GNU_compiler"            >> $TESTFILE
echo "#endif"                            >> $TESTFILE
echo "#endif"                            >> $TESTFILE

echo "#ifndef __GNUC__"                  >> $TESTFILE
echo "#ifdef __SUNPRO_C"                 >> $TESTFILE
echo "it_is_the_Sun_compiler"            >> $TESTFILE
echo "#endif"                            >> $TESTFILE
echo "#endif"                            >> $TESTFILE

# Not found any offical documention to suggest __DECC
# would be defined, but __DECCC is defined for C++,
# and a Google on __DECC shows many hits.
echo "#ifdef __digital__"                >> $TESTFILE
echo "#ifdef __DECC"                     >> $TESTFILE
echo "it_is_the_HP_Tru64_compiler"       >> $TESTFILE
echo "#endif"                            >> $TESTFILE
echo "#endif"                            >> $TESTFILE

echo "#ifdef __linux__"                  >> $TESTFILE
echo "#ifdef __DECCC"                   >> $TESTFILE
echo "it_is_the_HP_Alpha_Linux_compiler" >> $TESTFILE
echo "#endif"                            >> $TESTFILE
echo "#endif"                            >> $TESTFILE

echo "#ifdef __HP_cc"                    >> $TESTFILE
echo "it_is_the_HP_HPUX_compiler"        >> $TESTFILE
echo "#endif"                            >> $TESTFILE

echo "#ifdef __xlC__"                    >> $TESTFILE
echo "it_is_the_IBM_AIX_compiler"        >> $TESTFILE
echo "#endif"                            >> $TESTFILE

${CC} -E $TESTFILE | grep it_is_the_Intel_improved_GCC >/dev/null 2>&1
if [ $? = 0 ]; then
   echo Intel_improved_GCC
   rm $TESTFILE
   exit 0
fi

${CC} -E $TESTFILE | grep it_is_the_Sun_improved_GCC >/dev/null 2>&1
if [ $? = 0 ]; then
   echo Sun_improved_GCC
   rm $TESTFILE
   exit 0
fi

${CC} -E $TESTFILE | grep it_is_the_GNU_compiler  >/dev/null 2>&1
if [ $? = 0 ]; then
   echo GNU
   rm $TESTFILE
   exit 0
fi

${CC} -E $TESTFILE | grep it_is_the_Sun_compiler  >/dev/null 2>&1
if [ $? = 0 ]; then
   echo Sun_on_Solaris
   rm $TESTFILE
   exit 0
fi

# HP make compilers for linux, HP-UX and tru64, which complicates matters.

${CC} -E $TESTFILE | grep it_is_the_HP_Tru64_compiler >/dev/null 2>&1
if [ $? = 0 ]; then
   echo HP_on_Tru64
   rm $TESTFILE
   exit 0
fi

${CC} -E $TESTFILE | grep it_is_the_HP_Alpha_Linux_compiler >/dev/null 2>&1
if [ $? = 0 ]; then
   echo HP_on_Alpha_Linux
   rm $TESTFILE
   exit 0
fi

${CC} -E $TESTFILE | grep it_is_the_HP_HPUX_compiler >/dev/null 2>&1
if [ $? = 0 ]; then
   echo HP_on_HPUX
   rm $TESTFILE
   exit 0
fi

${CC} -E $TESTFILE | grep it_is_the_IBM_AIX_compiler >/dev/null 2>&1
if [ $? = 0 ]; then
   echo IBM_on_AIX
   rm $TESTFILE
   exit 0
fi

# Exit saying 'Unknown' if the type of the compiler could not be found.
echo Unknown
rm $TESTFILE
exit 0


