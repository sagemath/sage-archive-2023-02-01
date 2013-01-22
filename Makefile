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

# ssl: build Sage, and also install pyOpenSSL. This is necessary for
# running the secure notebook. This make target requires internet
# access. Note that this requires that your system have OpenSSL
# libraries and headers installed. See README.txt for more
# information.
ssl: all
	./sage -i pyopenssl

build-serial: SAGE_PARALLEL_SPKG_BUILD = no
build-serial: build

# Start Sage if the file local/etc/sage-started.txt does not exist
# (i.e. when we just installed Sage for the first time).
start: build
	[ -f local/etc/sage-started.txt ] || local/bin/sage-starts

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

doc-html-mathjax: build
	$(PIPE) "./sage --docbuild --no-pdf-links all html -j $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a dochtml.log"

# Keep target 'doc-html-jsmath' for backwards compatibility.
doc-html-jsmath: doc-html-mathjax

doc-pdf: build
	$(PIPE) "./sage --docbuild all pdf $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a docpdf.log"

doc-clean:
	@echo "Deleting devel/sage/doc/output..."
	rm -rf devel/sage/doc/output

clean:
	@echo "Deleting spkg/build..."
	rm -rf spkg/build
	@echo "Deleting spkg/archive..."
	rm -rf spkg/archive

distclean: clean
	@echo "Deleting all remaining traces of builds, tests etc. ..."
	rm -rf local
	rm -f spkg/Makefile
	rm -rf spkg/installed
	rm -rf spkg/logs
	rm -rf spkg/optional
	rm -f install.log
	rm -f dochtml.log docpdf.log
	rm -f test.log testall.log testlong.log ptest.log ptestlong.log
	rm -f start.log
	rm -f spkg/parallel_make.cfg
	rm -rf dist
	rm -rf devel
	rm -rf doc
	rm -rf examples
	rm -rf sage-python
	rm -rf matplotlibrc
	rm -rf tmp
	rm -f .BUILDSTART

micro_release:
	bash -c ". spkg/bin/sage-env && local/bin/sage-micro_release"

text-expand:
	./spkg/bin/text-expand

text-collapse:
	./spkg/bin/text-collapse

TESTPRELIMS = local/bin/sage-starts
# The [a-z][a-z] matches all directories whose names consist of two
# lower-case letters, to match the language directories.
TESTDIRS = devel/sage/doc/common devel/sage/doc/[a-z][a-z] devel/sage/sage

test: all # i.e. build and doc
	$(TESTPRELIMS)
	$(PIPE) "./sage -t --sagenb $(TESTDIRS) 2>&1" "tee -a test.log"

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
		echo >&2 "Set the environment variable DESTDIR to the install path."; \
		exit 1; \
	fi
	# Make sure we remove only an existing directory. If $(DESTDIR)/sage is
	# a file instead of a directory then the mkdir statement later will fail
	if [ -d "$(DESTDIR)"/sage ]; then \
		rm -rf "$(DESTDIR)"/sage; \
	fi
	mkdir -p "$(DESTDIR)"/sage
	mkdir -p "$(DESTDIR)"/bin
	cp -Rp * "$(DESTDIR)"/sage
	rm -f "$(DESTDIR)"/bin/sage
	ln -s ../sage/sage "$(DESTDIR)"/bin/sage
	"$(DESTDIR)"/bin/sage -c # Run sage-location


.PHONY: all build build-serial start install \
	doc doc-html doc-html-jsmath doc-html-mathjax doc-pdf \
	doc-clean clean	distclean \
	test check testoptional testall testlong testoptionallong testallong \
	ptest ptestoptional ptestall ptestlong ptestoptionallong ptestallong
