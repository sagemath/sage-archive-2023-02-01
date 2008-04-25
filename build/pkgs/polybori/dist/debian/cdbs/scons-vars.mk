# -*- mode: makefile; coding: utf-8 -*-
# Copyright © 2005 Matthew A. Nicholson <matt@matt-land.com>
# Description: Defines useful variables for packages which have a SConstruct
#              file
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

ifndef _cdbs_class_scons_vars
_cdbs_class_scons_vars := 1

include $(_cdbs_class_path)/langcore.mk$(_cdbs_makefile_suffix)

DEB_SCONS_ENVVARS =
DEB_SCONS_INVOKE = $(DEB_SCONS_ENVVARS) scons --directory $(DEB_BUILDDIR) CFLAGS=$(if $(CFLAGS_$(cdbs_curpkg)),"$(CFLAGS_$(cdbs_curpkg))","$(CFLAGS)") CXXFLAGS=$(if $(CXXFLAGS_$(cdbs_curpkg)),"$(CXXFLAGS_$(cdbs_curpkg))","$(CXXFLAGS)")

# general options (passed on all scons commands)
DEB_SCONS_OPTIONS =

# build target and options (only passed on build)
DEB_SCONS_BUILD_TARGET =
DEB_SCONS_BUILD_OPTIONS =

# install target and options (only passed on install)
DEB_SCONS_INSTALL_TARGET = install
DEB_SCONS_INSTALL_OPTIONS =

# clean target
DEB_SCONS_CLEAN_TARGET = .

DEB_SCONS_CHECK_TARGET =

endif
