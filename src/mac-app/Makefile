
ifdef SAGE_ROOT
    $(error This makefile must be run outside of a sage shell)
endif

SAGE_ROOT=$(shell cd ../.. && pwd)

ifndef SAGE_VERSION
    SAGE_VERSION=$(shell source ../bin/sage-version.sh && echo $$SAGE_VERSION)
endif

ifndef SAGE_APP_TARGET_ARCH
    SAGE_APP_TARGET_ARCH=$(shell uname -m)
endif

OSX_VERSION=$(shell sw_vers -productVersion)

TARGET=sage-$(SAGE_VERSION)-OSX_$(OSX_VERSION)_$(SAGE_APP_TARGET_ARCH)
APP=$(TARGET)/SageMath-$(SAGE_VERSION).app
DMG=sage-$(SAGE_VERSION)-OSX_$(OSX_VERSION)_$(SAGE_APP_TARGET_ARCH).app.dmg

CONFIGURATION=Debug
SDKROOT=

SOURCE=*.h *.m


default: $(DMG)

debug:
	@echo SAGE_ROOT=$(SAGE_ROOT)
	@echo SAGE_VERSION=$(SAGE_VERSION)
	@echo OSX_VERSION=$(OSX_VERSION)
	@echo SAGE_APP_TARGET_ARCH=$(SAGE_APP_TARGET_ARCH)
	@echo CONFIGURATION=$(CONFIGURATION)
	@echo APP=$(APP)

$(APP): $(SOURCE)
	mkdir -p $(TARGET)
	xcodebuild -target Sage -configuration $(CONFIGURATION) \
	    ARCHS=$(SAGE_APP_TARGET_ARCH) \
	    ONLY_ACTIVE_ARCH=NO \
	    SDKROOT=$(SDKROOT)
	cp -pRL build/$(CONFIGURATION)/Sage.app $(APP)
	# Info.plist is a binary plist, so convert it for processing with sed.
	# I would just change it to be an xml plist, but xcode changes it back.
	plutil -convert xml1 $(APP)/Contents/Info.plist
	sed -i '' "s/SAGE_VERSION/$(SAGE_VERSION)/" $(APP)/Contents/Info.plist

$(APP)/Contents/Resources/sage: $(APP)
	mkdir -p $(APP)/Contents/Resources/sage
	rsync -av --progress \
	    $(SAGE_ROOT)/ $@ \
	    --exclude src/mac-app/build \
	    --exclude src/mac-app/$(TARGET)

$(TARGET)/README.txt: $(SAGE_ROOT)/src/bin/sage-README-osx.txt
	cp $^ $@

$(DMG): $(APP)/Contents/Resources/sage $(TARGET)/README.txt
	hdiutil create -srcfolder $(TARGET) -format UDBZ $(DMG)

.PHONY: clean default debug

clean:
	rm -rf $(APP) build $(DMG)
