# Main Makefile for Sage.

# The default target ("all") builds Sage and the whole (HTML) documentation.
#
# Target "build" just builds Sage.
#
# See below for targets to build the documentation in other formats,
# to run various types of test suites, and to remove parts of the build etc.

PIPE = spkg/pipestatus


all: start doc  # indirectly depends on build

build:
	cd spkg && \
	"../$(PIPE)" \
		"env SAGE_PARALLEL_SPKG_BUILD='$(SAGE_PARALLEL_SPKG_BUILD)' ./install all 2>&1" \
		"tee -a ../install.log"
	./sage -b

build-serial: SAGE_PARALLEL_SPKG_BUILD = no
build-serial: build

# Start Sage if the file local/lib/sage-started.txt does not exist
# (i.e. when we just installed Sage for the first time)
start: build
	[ -f local/lib/sage-started.txt ] || local/bin/sage-starts

# You can choose to have the built HTML version of the documentation link to
# the PDF version. To do so, you need to build both the HTML and PDF versions.
# To have the HTML version link to the PDF version, do
#
# $ ./sage --docbuild all html
# $ ./sage --docbuild all pdf
#
# For more information on the docbuild utility, do
#
# $ ./sage --docbuild -H
doc: doc-html

doc-html: build
	$(PIPE) "./sage --docbuild --no-pdf-links all html $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a dochtml.log"

doc-html-jsmath: build
	$(PIPE) "./sage --docbuild --no-pdf-links all html -j $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a dochtml.log"

doc-pdf: build
	$(PIPE) "./sage --docbuild all pdf $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a docpdf.log"

doc-clean:
	@echo "Deleting devel/sage/doc/output..."
	rm -rf devel/sage/doc/output

clean:
	@echo "Deleting spkg/build..."
	rm -rf spkg/build
	mkdir -p spkg/build
	@echo "Deleting spkg/archive..."
	rm -rf spkg/archive
	mkdir -p spkg/archive

distclean: clean
	@echo "Deleting all remaining traces of builds, tests etc. ..."
	rm -rf local
	rm -rf spkg/installed
	rm -rf spkg/logs
	rm -rf spkg/optional
	rm -f install.log
	rm -f dochtml.log docpdf.log
	rm -f test.log testall.log testlong.log ptest.log ptestlong.log
	rm -f start.log
	rm -f spkg/parallel_make.cfg
	rm -rf data
	rm -rf dist
	rm -rf devel
	rm -rf doc
	rm -rf examples
	rm -rf sage-python
	rm -rf spkg/build
	rm -rf spkg/archive
	rm -rf matplotlibrc
	rm -rf tmp
	rm -f .BUILDSTART

micro_release:
	. local/bin/sage-env && local/bin/sage-micro_release

text-expand:
	./spkg/base/text-expand

text-collapse:
	./spkg/base/text-collapse

TESTPRELIMS = local/bin/sage-starts
TESTDIRS = devel/sage/doc/common devel/sage/doc/de devel/sage/doc/en devel/sage/doc/fr devel/sage/doc/ru devel/sage/sage

test: all # i.e. build and (HTML) doc
	$(TESTPRELIMS)
	. local/bin/sage-env && sage-maketest

check: test

testall: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -t --sagenb --optional $(TESTDIRS) 2>&1" "tee -a testall.log"

testlong: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -t --sagenb --long $(TESTDIRS) 2>&1" "tee -a testlong.log"

testalllong: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -t --sagenb --optional --long $(TESTDIRS) 2>&1" "tee -a testalllong.log"

ptest: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -tp --sagenb $(TESTDIRS) 2>&1" "tee -a ptest.log"

ptestall: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -tp --sagenb --optional $(TESTDIRS) 2>&1" "tee -a ptestall.log"

ptestlong: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -tp --sagenb --long $(TESTDIRS) 2>&1" "tee -a ptestlong.log"

ptestalllong: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -tp --sagenb --optional --long $(TESTDIRS) 2>&1" "tee -a ptestalllong.log"


testoptional: testall # just an alias

testoptionallong: testalllong # just an alias

ptestoptional: ptestall # just an alias

ptestoptionallong: ptestalllong # just an alias


install:
	echo "Experimental use only!"
	if [ "$(DESTDIR)" = "" ]; then \
		echo "Set DESTDIR"; \
		exit 1; \
	fi
	mkdir -p $(DESTDIR)
	mkdir -p $(DESTDIR)/sage
	mkdir -p $(DESTDIR)/bin/
	cp -rpv * $(DESTDIR)/sage/
	python local/bin/sage-hardcode_sage_root $(DESTDIR)/sage/sage "$(DESTDIR)"/sage
	cp $(DESTDIR)/sage/sage $(DESTDIR)/bin/
	cd $(DESTDIR)/bin/; ./sage -c


.PHONY: all build build-serial start install \
	doc doc-html doc-html-jsmath doc-pdf \
	doc-clean clean	distclean \
	test check testoptional testall testlong testoptionallong testallong \
	ptest ptestoptional ptestall ptestlong ptestoptionallong ptestallong
