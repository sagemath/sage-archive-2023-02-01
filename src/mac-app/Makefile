# Makefile for building the Mac App
#
# This assumes that you have already built Sage. This Makefile then
# compiles the Mac App UI and combines it with the existing Sage
# install into a DMG file.
#
# The resulting DMG is as relocatable as the original Sage install,
# that is, not at all by default. It will contain paths pointing to
# the Sage tree that you started with (OK if you are only debugging
# the Mac App, just don't delete the original Sage directory tree).
#
# To build a relocatable Sage binary you have to use
# https://github.com/sagemath/binary-pkg


ifdef SAGE_ROOT
    $(error This makefile must be run outside of a sage shell)
endif

SAGE_ROOT:=$(shell cd ../.. && pwd)

ifndef SAGE_VERSION
    SAGE_VERSION:=$(shell source ../bin/sage-version.sh && echo $$SAGE_VERSION)
endif

# Meaningful values, on Intel, are "i386" and "x86_64". To prepare a
# 64-bit binary distribution on an older 64-bit OSX machine that boots
# into a 32-bit system, one would set SAGE_APP_TARGET_ARCH=x86_64
ifndef SAGE_APP_TARGET_ARCH
    SAGE_APP_TARGET_ARCH:=$(shell uname -m)
endif

OSX_VERSION:=$(shell sw_vers -productVersion)

TARGET:=sage-$(SAGE_VERSION)-OSX_$(OSX_VERSION)_$(SAGE_APP_TARGET_ARCH)
APP_NAME:=SageMath-$(SAGE_VERSION).app
APP:=$(TARGET)/$(APP_NAME)
DMG:=$(TARGET).app.dmg

CONFIGURATION=Debug
SDKROOT=

SOURCE=*.h *.m

default: $(DMG)
app: $(APP)/Contents/Resources/sage

# Display the variables
debug:
	@echo SAGE_ROOT=$(SAGE_ROOT)
	@echo SAGE_VERSION=$(SAGE_VERSION)
	@echo OSX_VERSION=$(OSX_VERSION)
	@echo SAGE_APP_TARGET_ARCH=$(SAGE_APP_TARGET_ARCH)
	@echo CONFIGURATION=$(CONFIGURATION)
	@echo APP=$(APP)

# Build the Mac App UI binary (not: Sage itself)
$(APP): $(SOURCE)
	mkdir -p $(TARGET)
	xcodebuild -target Sage -configuration $(CONFIGURATION) \
	    ARCHS=$(SAGE_APP_TARGET_ARCH) \
	    ONLY_ACTIVE_ARCH=NO \
	    SDKROOT=$(SDKROOT)
	mkdir -p $(APP)
	cp -pRL build/$(CONFIGURATION)/Sage.app/ $(APP)
	@# Info.plist is a binary plist, so convert it for processing with sed.
	@# I would just change it to be an xml plist, but xcode changes it back.
	plutil -convert xml1 $(APP)/Contents/Info.plist
	sed -i '' "s/SAGE_VERSION/$(SAGE_VERSION)/" $(APP)/Contents/Info.plist

# Copy Sage into the App directory structure
$(APP)/Contents/Resources/sage: $(APP)
	mkdir -p $(APP)/Contents/Resources/sage
	rsync -av --progress \
	    $(SAGE_ROOT)/ $@ \
	    --exclude .git/ \
	    --exclude upstream/ \
	    --exclude src/mac-app/build \
	    --exclude src/mac-app/$(TARGET)
	@# We touch the director here so that we don't keep rsyncing just
	@# because $(APP) is newer than (the rsync preserves time)
	test -d $@ && touch $@

# .DS_Store tells finder the icon position and background image to use
# when opening the dmg.	
$(TARGET)/.DS_Store: tools/createDSStore/createDSStore.py
	@# On some old Mac OS X, createDSStore fails when there already
	@# is an existing file
	rm -f $(TARGET)/.DS_Store
	tools/createDSStore/createDSStore.py $(TARGET) $(TARGET) $(APP_NAME)

$(TARGET)/.background.png: sage-dmg-background.png
	cp sage-dmg-background.png $(TARGET)/.background.png

$(TARGET)/Applications:
	rm -f $(TARGET)/Applications
	ln -s /Applications/ $(TARGET)/Applications

# Add the readme file
$(TARGET)/README.txt: $(SAGE_ROOT)/src/bin/sage-README-osx.txt
	cp $^ $@

non_app_files: $(TARGET)/README.txt $(TARGET)/.DS_Store $(TARGET)/.background.png $(TARGET)/Applications

sage_symlink: $(APP)/Contents/Resources/sage/sage
	rm -f $(APP)/sage
	ln -s Contents/Resources/sage/sage $(APP)/sage

# Compress the built directory into a DMG
# Dependencies need to be in that order because non_app_files
# needs the directory created by $(APP)/Contents/Resources/sage
$(DMG): $(APP)/Contents/Resources/sage sage_symlink non_app_files
	hdiutil create -srcfolder $(TARGET) -volname $(TARGET) -format UDBZ $(DMG)

.PHONY: clean default debug app non_app_files

clean:
	rm -rf $(APP) build $(DMG)
