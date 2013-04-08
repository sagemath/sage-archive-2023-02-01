# This is a makefile for the entire library.
# Copyright (C) 2004-2008 by Pawel Pilarczyk.
# This is free software. No warranty. Consult 'license.txt' for details.
# Created on October 26, 2004 by PwP. Last revision: September 27, 2005.


# ===========================================================================
# Settings.
# ===========================================================================

# Set up the home directory.
HOME = ../

# Determine the target configuration (see the file 'makecfg' for details).
include makecfg

# Include the compiler-dependent definitions.
include config/${TARGET}

# Determine the paths to files.
include makedirs


# ===========================================================================
# A list of modules that are part of the library.
# ===========================================================================

MODULES = ${wildcard ${SRC}*}
MODULES := ${MODULES:${SRC}%=%}
MODULES := ${filter-out CVS,${MODULES}}


# ===========================================================================
# A list of source .cpp files in all the modules of the library.
# ===========================================================================

FILES = ${foreach MODULE,${MODULES},${wildcard ${SRC}${MODULE}/*.cpp}}
#FILES := ${notdir ${FILES}}
FILES := ${basename ${FILES}}


# ===========================================================================
# A list of all the programs compiled with the library.
# ===========================================================================

PROGDEMO = ${wildcard ${EXAMPLES}*/makefile}
PROGDEMO := ${notdir ${PROGDEMO:%/makefile=%}}
PROGDEMO := ${filter-out CVS,${PROGDEMO}}

PROGTEST = ${wildcard ${TESTS}*/makefile}
PROGTEST := ${notdir ${PROGTEST:%/makefile=%}}
PROGTEST := ${filter-out CVS,${PROGTEST}}

PROGS = ${wildcard ${PROGRAMS}*/makefile}
PROGS := ${notdir ${PROGS:%/makefile=%}}
PROGS := ${filter-out CVS,${PROGS}}

PROGPRIV = ${wildcard ${PRIVATE}*/makefile}
PROGPRIV := ${notdir ${PROGPRIV:%/makefile=%}}
PROGPRIV := ${filter-out CVS,${PROGPRIV}}


# ===========================================================================
# Targets.
# ===========================================================================

# The main target: everything.
.PHONY: all
all: objdir libdir objects library programs

# Set the library file name.
LIBRARY = ${LIB}lib${LIBNAME}${LIBEXT}


# ===========================================================================
# Make the object files.
# ===========================================================================

# The object directory.
.PHONY: objdir
objdirfile = $(wildcard ${OBJ})
objdir:
ifeq (${objdirfile},)
	-${MKDIR} ${OBJ}
else
endif

# Make the object files.
.PHONY: objects
objects: ${MODULES}

.PHONY: ${MODULES}
${MODULES}: %:
	${MAKE} -f makemain target=${TARGET} module=$@


# ===========================================================================
# Make the library.
# ===========================================================================

# The library directory.
.PHONY: libdir
libdirfile = $(wildcard ${LIB})
libdir:
ifeq (${libdirfile},)
	-${MKDIR} ${LIB}
else
endif

.PHONY: library
library: objects ${LIBRARY}

# Determine the list of all the object files in the project.
OBJFILES = ${FILES:${SRC}%=${OBJ}%${OBJEXT}}

# Make the library.
${LIBRARY}: ${OBJFILES}
	${MAKELIB} ${LIBRARY} ${OBJFILES}
	${IDXLIB} ${LIBRARY}


# ===========================================================================
# Make all the example, test and other programs.
# ===========================================================================

.PHONY: examples
examples: objects library ${PROGDEMO}

.PHONY: tests
tests: objects library ${PROGTEST}

.PHONY: programs
programs: objects library ${PROGS}

.PHONY: private
private: objects library ${PROGPRIV}

.PHONY: ${PROGDEMO}
${PROGDEMO}: %:
	${MAKE} -C ${EXAMPLES}$@

.PHONY: ${PROGTEST}
${PROGTEST}: %:
	${MAKE} -C ${TESTS}$@

.PHONY: ${PROGS}
${PROGS}: %:
	${MAKE} -C ${PROGRAMS}$@

.PHONY: ${PROGPRIV}
${PROGPRIV}: %:
	${MAKE} -C ${PRIVATE}$@


# ===========================================================================
# Clean all the object files, libraries, executables and dependencies.
# ===========================================================================

# delete all the object files
.PHONY: cleanobj
cleanobj: ${PROGDEMO:%=CLEANOBJ_%} ${PROGTEST:%=CLEANOBJ_%} \
	${PROGS:%=CLEANOBJ_%} ${PROGPRIV:%=CLEANOBJ_%}
	-${REMOVE} ${OBJFILES}

.PHONY: ${PROGDEMO:%=CLEANOBJ_%}
${PROGDEMO:%=CLEANOBJ_%}: %:
	${MAKE} -C ${EXAMPLES}${@:CLEANOBJ_%=%} cleanobj

.PHONY: ${PROGTEST:%=CLEANOBJ_%}
${PROGTEST:%=CLEANOBJ_%}: %:
	${MAKE} -C ${TESTS}${@:CLEANOBJ_%=%} cleanobj

.PHONY: ${PROGS:%=CLEANOBJ_%}
${PROGS:%=CLEANOBJ_%}: %:
	${MAKE} -C ${PROGRAMS}${@:CLEANOBJ_%=%} cleanobj

.PHONY: ${PROGPRIV:%=CLEANOBJ_%}
${PROGPRIV:%=CLEANOBJ_%}: %:
	${MAKE} -C ${PRIVATE}${@:CLEANOBJ_%=%} cleanobj

# delete all the executables
.PHONY: cleanexe
cleanexe: ${PROGDEMO:%=CLEANEXE_%} ${PROGTEST:%=CLEANEXE_%} \
	${PROGS:%=CLEANEXE_%} ${PROGPRIV:%=CLEANEXE_%}

.PHONY: ${PROGDEMO:%=CLEANEXE_%}
${PROGDEMO:%=CLEANEXE_%}: %:
	${MAKE} -C ${EXAMPLES}${@:CLEANEXE_%=%} cleanexe

.PHONY: ${PROGTEST:%=CLEANEXE_%}
${PROGTEST:%=CLEANEXE_%}: %:
	${MAKE} -C ${TESTS}${@:CLEANEXE_%=%} cleanexe

.PHONY: ${PROGS:%=CLEANEXE_%}
${PROGS:%=CLEANEXE_%}: %:
	${MAKE} -C ${PROGRAMS}${@:CLEANEXE_%=%} cleanexe

.PHONY: ${PROGPRIV:%=CLEANEXE_%}
${PROGPRIV:%=CLEANEXE_%}: %:
	${MAKE} -C ${PRIVATE}${@:CLEANEXE_%=%} cleanexe

# delete all the library files
.PHONY: cleanlib
cleanlib:
	-${REMOVE} ${MODULES:%=${LIB}lib${LIBNAME}%${LIBEXT}} ${LIBRARY}
ifeq (${target},bcc)
	-${REMOVE} ${LIB}lib${LIBNAME}*.BAK
endif

# delete all the dependency files (generated with "filedeps")
.PHONY: cleandeps
cleandeps:
	-${REMOVE} ${MODULES:%=auto_dep/%} \
	${PROGDEMO:%=${EXAMPLES}%/makedeps} \
	${PROGTEST:%=${TESTS}%/makedeps} \
	${PROGS:%=${PROGRAMS}%/makedeps} \
	${PROGPRIV:%=${PRIVATE}%/makedeps}

.PHONY: purge
purge: cleanobj cleanlib cleanexe

.PHONY: clean
clean: cleanobj


# That's all, folks!

