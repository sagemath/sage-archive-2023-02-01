# -*- mode: makefile; coding: utf-8 -*-
# Copyright © 2005 Matthew A. Nicholson <matt@matt-land.com>
# Description: Builds and cleans packages which have a SConstruct file
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# 02111-1307 USA.


ifndef _cdbs_bootstrap
_cdbs_scripts_path ?= /usr/lib/cdbs
_cdbs_rules_path ?= /usr/share/cdbs/1/rules
_cdbs_class_path ?= /usr/share/cdbs/1/class
endif

ifndef _cdbs_class_scons
_cdbs_class_scons := 1

include $(_cdbs_rules_path)/buildcore.mk$(_cdbs_makefile_suffix)
#include $(_cdbs_class_path)/scons-vars.mk$(_cdbs_makefile_suffix)
include debian/cdbs/scons-vars.mk$(_cdbs_makefile_suffix)

DEB_PHONY_RULES += scons-clean

common-build-arch common-build-indep:: debian/stamp-scons-build
debian/stamp-scons-build:
	$(DEB_SCONS_INVOKE) $(DEB_SCONS_BUILD_TARGET) $(DEB_SCONS_OPTIONS) $(DEB_SCONS_BUILD_OPTIONS)
	touch debian/stamp-scons-build

clean:: scons-clean
scons-clean::
	$(DEB_SCONS_INVOKE) $(DEB_SCONS_CLEAN_TARGET) $(DEB_SCONS_OPTIONS) --keep-going --clean || true
	rm -f debian/stamp-scons-build
	rm -rf .sconf_temp/
	rm -f .sconsign.dblite config.log

common-install-arch common-install-indep:: common-install-impl
common-install-impl::
	@if test -n "$(DEB_SCONS_INSTALL_TARGET)"; then \
	  echo $(DEB_SCONS_ENVVARS) scons --directory $(DEB_BUILDDIR) $(DEB_SCONS_INSTALL_TARGET); \
	  $(DEB_SCONS_INVOKE) $(DEB_SCONS_INSTALL_TARGET) $(DEB_SCONS_OPTIONS) $(DEB_SCONS_INSTALL_OPTIONS); \
	 else \
	   echo "DEB_SCONS_INSTALL_TARGET unset, skipping default scons.mk common-install target"; \
	 fi

ifeq (,$(findstring nocheck,$(DEB_BUILD_OPTIONS)))
common-post-build-arch common-post-build-indep:: common-post-build-impl
common-post-build-impl::
	@if test -n "$(DEB_SCONS_CHECK_TARGET)"; then \
	  echo $(DEB_SCONS_INVOKE) $(DEB_SCONS_CHECK_TARGET); \
	  $(DEB_SCONS_INVOKE) $(DEB_SCONS_CHECK_TARGET) $(DEB_SCONS_OPTIONS); \
	else \
	   echo "DEB_SCONS_CHECK_TARGET unset, not running checks"; \
	fi
endif

endif
