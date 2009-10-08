###################################################################
#
# Copyright (C) 2009,  Ralf Hemmecke <ralf@hemmecke.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
###################################################################
#
# This Makefile is intended to build two Sage packages which (when
# installed into Sage) will make FriCAS and its interface to Aldor
# available.
#
# This Makefile assumes that Aldor is installed both
# a) for its own execution, i.e. creating the spkg and
# b) for the installation of the aldor interface at the
#    user side.
#
# In both cases Aldor must be installed in such a way that the FriCAS
# configure script can find it.
# Furthermore, it assumes that 'sage' calls a recent version of sage.
# It should be called like:
#
#   make FRICAS_SRC=/path/to/fricas-sources FRICAS_LISP=ecl
#
# If FRICAS_SRC is not given, the Makefile checks whether "./fricas"
# contains the FriCAS sources and uses this if successful.
#
# If FRICAS_LISP is not given, then ecl is assumed.
#
# This Makefile then builds FriCAS and the Aldor interface inside
# a "./build" directory.
#
# After a successful build the script collects all data and creates
# two Sage packages:
#
#   fricas-X.Y.Z.spkg
#   fricasaldor-X.Y.Z.spkg
#
# inside a new directory "./spkg".
# You should then upload these two files to
#
#   http://sagemath.org/packages/optional/
#
# in order to make the files available for automatic download in Sage.
#
# It is intended that a user (who has sage and aldor already
# installed) then simply has to say 'sage -i fricasaldor-X.Y.Z' in
# order to install fricas and the fricas-aldor interface.
#
# We need two .spkg here since fricasaldor-X.Y.Z.spkg will contain the
# file libaxiom.al which must be distributed under the Aldor Public
# License 2.0 (http://www.aldor.org/AldorPublicLicense2_0.html)
# because it contains code from aldor.org.
#
###################################################################

# The makefile must be called via "make FRICAS_SRC=/path/to/fricas-sources".
# If the FRICAS_SRC variable is not set, the local directory "fricas" will
# be assumed.
ifndef FRICAS_SRC
FRICAS_SRC=$(shell if test -d fricas; then echo fricas; fi)
endif

# Figure out the absolute path of $(FRICAS_SRC).
abs_srcdir = $(shell cd $(FRICAS_SRC); pwd)

# The absolute path to the configure script.
CONFIGURE=$(abs_srcdir)/configure

# Which lisp variant should be used for creating the cached lisp files.
# Although ecl will be used in Sage, there should be no problem to
# chose another lisp for the creation of the the spkg, i.e.
# FRICAS_LISP=sbcl should also work.
# This variable should give the full path to the executable. Otherwise
# you must make sure that the respective lisp is in your PATH.
FRICAS_LISP=ecl

# If the NOWEB variable has not been given we assume a noweb tarball
# in the zips subdirectory.
ifndef NOWEB
NOWEB=zips/noweb-2.10a.tgz
endif
NOWEBPATH=$(shell cd `dirname $(NOWEB)`; pwd)
NOWEBFILE=$(shell basename $(NOWEB))
NW = $(NOWEBPATH)/$(NOWEBFILE)

# If the SPADHELP variable has not been given we assume that the helpfiles
# are in the share/spadhelp subdirectory.
ifndef SPADHELP
SPADHELP=share/spadhelp
endif
# Make the directory an absolute path.
HELP_DIR=$(shell cd $(SPADHELP); pwd)

# We need Sage to actually create the spkg from the directories.
SAGE=sage

# We need Mercurial to prepare the spkg.
HG=hg

# This is going to be the version of the SPKG.
V:=$(shell grep "^PACKAGE_VERSION=" ${CONFIGURE} 2>/dev/null | sed "s/PACKAGE_VERSION=//;s/'//g")

# Abbreviate the spkg directories.
F  = spkg/fricas-$V
FA = spkg/fricasaldor-$V

all: $F.spkg $(FA).spkg

prerequisites: check-fricas-src check-configure check-noweb check-helpfiles check-sage check-hg

check-fricas-src:
	@(if test -z "$(FRICAS_SRC)"; then \
	  echo "Could not find FriCAS source dir."; \
	  echo "Call via 'make FRICAS_SRC=/path/to/fricas-sources'."; \
	  exit 1; \
	fi)

check-configure:
	@(if ! test -f "$(CONFIGURE)"; then \
	  echo "Could not find configure file in FriCAS source dir."; \
	  echo "Call via 'make FRICAS_SRC=/path/to/fricas-sources'."; \
	  exit 1; \
	fi)

check-noweb:
	@(if ! test -f "$(NW)"; then \
	  echo "Could not find the noweb tarball."; \
	  echo "Call via"; \
	  echo "  make FRICAS_SRC=/path/to/fricas-sources NOWEB=/path/to/noweb-tarball"; \
	  echo "or execute 'make noweb' to get the noweb tarball from the Internet."; \
	  exit 1; \
	fi)

check-helpfiles:
	@(if ! test -f "$(HELP_DIR)/help.help"; then \
	  echo "Could not find the helpfiles."; \
	  echo "Call via"; \
	  echo "  make FRICAS_SRC=/path/to/fricas-sources SPADHELP=/path/to/spadhelp-dir"; \
	  echo "or execute 'make spadhelp' to get the helpfiles tarball from the Internet."; \
	  exit 1; \
	fi)

check-sage:
	@(if ! $(SAGE) --version >/dev/null; then \
	  echo "Could not find the sage executable."; \
	  echo "Call via 'make SAGE=/path/to/sage'."; \
	  exit 1; \
	fi)

check-hg:
	@(if ! hg --version >/dev/null; then \
	  echo "Could not find the hg executable."; \
	  echo "Mercurial must be installed on your computer."; \
	  echo "Call via 'make HG=/path/to/hg'."; \
	  exit 1; \
	fi)

# Clone the templates
$(FA).clone $F.clone: spkg/%.clone: prerequisites spkg/.dir
	-rm -rf spkg/$*
	hg clone . spkg/$*
	touch $@

build/src/aldor/al/libaxiom.al: build/.dir $(CONFIGURE)
	(cd build; $(CONFIGURE) --enable-aldor --with-lisp=$(FRICAS_LISP))
	(cd build; $(MAKE))

# Make graphics files and update ht.db.
# Note that this step requires an X server and will temporarily pop up
# lots of windows with graphics on your screen.
$F.gphts: build/src/aldor/al/libaxiom.al $F.clone
	(cd build/src/paste; $(MAKE) gphts)
	(cd build; $(MAKE))
	touch $@

# Copy the fricas src tarball with cached lisp and graphics.
$F.copy: $F.gphts
	(cd build; $(abs_srcdir)/src/scripts/mkdist.sh --copy_lisp --copy_gphts --copy_phts --copy_noweb="$(NW)" --copy_help=$(HELP_DIR))
	mv build/dist $F/src

# Copy the information into the appropriate spkg/src dir.
$(FA)/src/axiom.as: build/src/aldor/axiom.as $(FA).clone
	-mkdir $(@D)
	cp $< $@

$(FA)/src/libaxiom.al: build/src/aldor/al/libaxiom.al $(FA)/src/axiom.as
	cp $< $@

$(FA).copy: $(FA)/src/libaxiom.al $(FA)/src/axiom.as

# Adjust the spkg data.
$(FA).adjust $F.adjust: spkg/%-$V.adjust: spkg/%-$V.copy
	(cd spkg/$*-$V && cp spkg-install.$* spkg-install)
	(cd spkg/$*-$V && cp SPKG.txt.$* SPKG.txt)
	touch $@

$(FA).spkg $F.spkg: spkg/%.spkg: prerequisites $(CONFIGURE) spkg/%.adjust
	(cd spkg; sage -pkg $*)


###################################################################
# Creation of build subdirectories.
%/.dir:
	mkdir -p $(patsubst %/.dir,%,$@)
	touch -t 199901010000 $@


###################################################################
# Cleanup
clean:
	-rm -rf build src spkg

distclean: clean
	-rm -rf zips share

###################################################################
# This 'noweb' target is just to get noweb over the Internet if no
# noweb variable has been given.
noweb: zips/$(NOWEB)
zips/$(NOWEB):
	-mkdir zips
	(cd zips; wget http://axiom.svn.sourceforge.net/svnroot/axiom/trunk/axiom/$(NOWEB))

###################################################################
# FriCAS helpfiles are needed for the )help command to work.
spadhelp: share/spadhelp/help.help
share/spadhelp/help.help: share/fricas-spadhelp.tgz
	cd share; tar zxf fricas-spadhelp.tgz
share/fricas-spadhelp.tgz: share/.dir
	(cd share; wget http://sage.math.washington.edu/home/hemmecke/pub/fricas-spadhelp.tgz)
