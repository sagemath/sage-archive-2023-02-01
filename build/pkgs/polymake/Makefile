#  Copyright (c) 1997-2006
#  Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
#  http://www.math.tu-berlin.de/polymake,  mailto:polymake@math.tu-berlin.de
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation; either version 2, or (at your option) any
#  later version: http://www.gnu.org/licenses/gpl.txt.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#-----------------------------------------------------------------------------
#  $Project: polymake $$Id: Makefile,v 1.126.2.5 2006/03/14 16:42:44 gawrilow Exp $

###########################################################
#
#  Building and installation
#
.PHONY: all all-java configure clean clean-bin clean-java distclean TAGS RTAGS all-perlx perlx \
        install install-devel install-bin install-shared install-perlx install-doc install-devel-bin install-devel-shared

ProjectTop := $(CURDIR)
_allow_lacking_BuildDir := y
include support/utils.make

ifeq ($(wildcard $(BuildDir)/conf.make),)
  BuildDir :=
endif

Configure := $(if $(Arch), Arch="$(Arch)") TOP=$(CURDIR) MAKE=$(MAKE) $(PERL) support/configure.pl

ifndef BuildDir

all install% config-% :
	@$(Configure) -quick
	@$(MAKE) --no-print-directory $@

configure :
	@$(Configure)

else  # BuildDir exists

### load configuration
include $(BuildDir)/conf.make

PerlxDir := $(BuildDir)/perlx-$(PerlVersion)-$(PerlArch)

override _has_export := $(wildcard .export)

ifdef _has_export
  include .export
endif

Apps := $(filter-out $(NoInstall) CVS, $(notdir $(wildcard apps/*)))
Modules := $(filter-out $(NoInstall) CVS, $(notdir $(wildcard modules/*)))
DefaultAppName := $(firstword $(filter polytope, $(Apps)) $(Apps))

ifneq ($(Apps),)

define _MakeGMP
	$(if $(filter-out -,$(GMPbuildDir)), $(MAKE) -C $(GMPbuildDir) $(1))

endef
define _MakeAppBin
	$(MAKE) -C $(BuildDir)/apps/$(1) AppName=$(1) $(2)

endef

ifdef JavaBuildDir
  JavaIgnore := JavaUnconfigured

  define _MakeAppJava
	$(if $(wildcard apps/$(1)/java/*.java),$(MAKE) -j1 -C $(JavaBuildDir)/$(PackageBase)/$(1) JavaIgnore="$(JavaIgnore)" $(2),:)

  endef
else
  _MakeAppJava :=
endif

define _MakeAppsBin
	$(foreach a,$(Apps),$(call _MakeAppBin,$a,BeFirst=$(filter $a,$(firstword $(Apps))) $(1)))
endef

Jars := $(wildcard $(foreach a, $(Apps) $(Modules), java_build/$a.jar))

define _MakeAppsJava
	$(if $(JavaBuildDir), \
	  $(foreach a,$(Apps),$(call _MakeAppJava,$a,BeFirst=$(filter $a,$(firstword $(Apps))) $(1))), \
	  $(if $(filter install,$(1)), \
	    $(if $(Jars), \
	      $(call _CreateDir,-U $(InstallTop)/jars) \
	      $(PERL) support/install.pl -m 444 $(Jars) $(InstallTop)/jars)))

endef

all : $(if $(_has_export),perlx,all-perlx)
	@+$(call _MakeGMP,all)
	@+$(call _MakeAppsBin,all)
	@+$(call _MakeAppsJava,all)

all-perlx:
	[ -d $(PerlxDir) ] || $(PERL) $(INSTALL_PL) -d -m 755 $(PerlxDir)
	cd $(PerlxDir); TOP=$(CURDIR) $(PERL) $(CURDIR)/perl/ext/Makefile.PL
	@$(MAKE) -C $(PerlxDir) all

all-java clean-java install-java:
	@+$(call _MakeAppsJava,$(patsubst %-java,%,$@))

else  # !Apps

all:
	@$(MAKE) -C $(BuildDir)/lib all-lib AppName=lib PTLonly=yes
endif

### installation

define _InstallSubdir
$(PERL) support/install.pl -m $(DirMask) -U -X CVS $(2) $(1) $(InstallTop)/$(1)

endef
define _CreateDir
$(PERL) support/install.pl -d -m $(DirMask) $(1)

endef

ifneq ($(Apps),)

define _InstallDirLinks
$(PERL) support/install_links.pl -m 111 $(1) $(InstallLinks) $(InstallTop) $(DESTDIR)

endef
define _InstallLinks
	$(call _InstallDirLinks, $(InstallArch)/$(1)/bin)
	$(if $(wildcard $(1)/scripts), $(call _InstallDirLinks, $(InstallTop)/$(1)/scripts))
endef
ifdef _has_export
define _FilterRules
$(PERL) support/filter_rulefiles $(1) $(InstallTop) NoInstall

endef
else
_FilterRules :=
endif

define _InstallShared
	$(if $(wildcard $(1)/scripts), $(call _InstallSubdir,$(1)/scripts,-P $(PerlExe)))
	$(foreach d, perllib rules, $(if $(wildcard $(1)/$d), $(call _InstallSubdir,$(1)/$d)))
	$(if $(wildcard $(1)/rules), $(call _FilterRules,$(1)/rules))
endef

install: install-shared install-bin install-doc
install-devel: install-devel-shared install-devel-bin

install-bin: all install-perlx
	@+$(call _MakeGMP, install-libLTLIBRARIES)
	@+$(call _MakeAppsBin, install)
ifdef InstallLinks
	$(call _CreateDir,$(InstallLinks))
	$(call _InstallDirLinks, $(InstallTop)/perl)
	$(foreach m, $(Modules), $(call _InstallLinks,modules/$m))
	$(foreach a, $(Apps), $(call _InstallLinks,apps/$a))
	$(PERL) support/convert_main_script $(InstallLinks) \
		PerlExe=$(PerlExe) $(if $(PREFIX),PREFIX=$(PREFIX)) \
		InstallTop=$(FinalInstallTop) InstallArch=$(FinalInstallArch) $(if $(Arch),Arch="$(Arch)") \
		DefaultAppName=$(DefaultAppName)
endif

install-perlx : all-perlx
	$(call _CreateDir,$(InstallArch))
ifndef RPM_BUILD_ROOT
	@$(MAKE) -C $(PerlxDir) pure_install InstallDir=$(FinalInstallArch)
endif

install-devel-bin:
	@+$(call _MakeAppsBin, install-lib)
	@rm -f $(InstallArch)/Makefile
	@{ sed -e '/Install.*=/ s/^/override /' \
	       -e '/INSTALL_PL=/ s:=.*:=$$(InstallTop)/support/install.pl:' \
	       -e '/GMPbuildDir=/ s:=.*:=-:' \
	       -e '/DESTDIR=/ { s/^/override /; q; }' \
	       $(BuildDir)/conf.make; \
	   echo 'override ProjectTop=$$(InstallTop)'; \
	   echo 'ExtendApp ?= $(DefaultAppName)'; \
	   echo 'include $$(ProjectTop)/support/build.make'; \
	 } >$(InstallArch)/Makefile
	@chmod 444 $(InstallArch)/Makefile

install-shared:
	$(call _CreateDir,$(InstallTop))
	[ -z "$(AutoconfFile)" -o ! -f "$(AutoconfFile)" ] || $(PERL) support/install.pl -m 644 $(AutoconfFile) $(InstallTop)
	$(call _InstallSubdir,perl,-X polymake -X testsuite -X ext)
	$(foreach m, $(Modules), $(call _InstallShared,modules/$m))
	$(foreach a, $(Apps), $(call _InstallShared,apps/$a))
	$(call _MakeAppsJava,install)
	$(call _InstallSubdir,povray)
	$(call _InstallSubdir,scripts)

install-devel-shared:
	$(call _InstallSubdir,support)
	@+$(call _MakeAppsBin, install-include)

ifneq ($(wildcard doc),)
install-doc : perlx
	$(MAKE) -C doc install PerlxDir=$(BuildDir)/perlx
else
install-doc:
	$(PERL) support/install.pl -m $(DirMask) -U doc_build $(InstallDoc)
endif

else  # !Apps

install: install-shared install-bin

install-shared:
	$(call _CreateDir,$(InstallTop))
	@$(MAKE) -C $(BuildDir)/lib install-include AppName=lib PTLonly=yes

install-bin:
	$(call _CreateDir,$(InstallArch))
	@$(MAKE) -C $(BuildDir)/lib install-lib AppName=lib PTLonly=yes

endif

### manual re-configuration
configure :
	@$(Configure) $(BuildDir)/conf.make

config-% :
	@echo $($*)

ifdef _has_export

perlx : all-perlx
	@$(MAKE) -C $(PerlxDir) pure_install InstallDir=$(BuildDir)

### maintenance
tagsFLAGS = -R -e -f $@ --exclude=CVS --exclude='.\#*' --exclude='\#*' --exclude='*~'

TAGS:
	ctags $(tagsFLAGS) --language-force=c++ lib {apps,modules}/*/{src,include}
RTAGS:
	ctags $(tagsFLAGS) --language-force=perl --exclude=testsuite --exclude=ext perl \
			   {apps,modules}/*/{perllib,rules} scripts apps/*/scripts support/*.pl
	ctags -e -f $@ --append --language-force=c perl/ext/*.{xs,h}

include support/export.make

endif  # _has_export

endif  # BuildDir exists

### cleanup

clean-bin:
	@+$(call _MakeGMP, clean)
ifdef BuildDir
	@+$(call _MakeAppsBin, clean)
	rm -rf $(PerlxDir) $(BuildDir)/perlx/$(PerlVersion)/$(PerlArch)
endif

clean : clean-bin
ifdef JavaBuildDir
	@+$(call _MakeAppsJava, clean)
endif
ifneq ($(wildcard doc),)
	@$(MAKE) -C doc clean PerlxDir=none
endif

distclean : clean
	[ ! -d lib/PTL ] || rm -rf build*


# Local Variables:
# mode: Makefile
# End:
