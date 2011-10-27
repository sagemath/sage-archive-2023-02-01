# Main Makefile for Sage.

# The default target ("all") builds Sage and the whole (HTML) documentation.
#
# Target "build" just builds Sage.
#
# See below for targets to build the documentation in other formats,
# to run various types of test suites, and to remove parts of the build etc.

# TODO:
#   * Consistently do "./sage -b" before running tests, and *before*
#     building/updating the documentation.
#   * Shorten description of NUM_THREADS below?

# NUM_THREADS is the number of threads to use for parallel testing (and
# sometime in the future, parallel building).  If this is 0, then it
# will be set to the number of processors, with a default maximum of 8
# -- see sage-ptest.
#
# The detection of number of processors might not be reliable on some
# platforms. On a Sun SPARC T5240 (t2.math), the number of processors
# reported by multiprocessing.cpu_count() might not correspond to the
# actual number of processors. See ticket #6283.
# Python's multiprocessing.cpu_count() actually returns the number of
# *hardware threads*, which is >= number of cores.
#
# WARNING: if your machine has <= 8 cpus (according to cpu_count() and
# you *don't* want to use that many threads for parallel doctesting,
# change the value of NUM_THREADS to a (sensible) positive integer. If
# cpu_count() reports > 8, then if NUM_THREADS is 0, only 8 threads will
# be used. The default value is zero.
NUM_THREADS = 0 # 0 interpreted as min(8, multiprocessing.cpu_count())

PIPE = spkg/pipestatus


all: start doc  # indirectly depends on build

# $(PIPE):
#	# We could generate it here if it doesn't exist, or a specific version
#	# depending on the shell (version), or even redefine $(PIPE) here.
#	# The following would require $(PIPE) to be a phony target:
#	test -x $@ # or make it executable if it exists; sanity check only anyway

build: $(PIPE)
	cd spkg && \
	"../$(PIPE)" \
		"env SAGE_PARALLEL_SPKG_BUILD='$(SAGE_PARALLEL_SPKG_BUILD)' ./install all 2>&1" \
		"tee -a ../install.log"

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
# $ ./sage -docbuild all html
# $ ./sage -docbuild all pdf
#
# For more information on the docbuild utility, do
#
# $ ./sage -docbuild -H
doc: doc-html

doc-html: build # (already) indirectly depends on $(PIPE)
	$(PIPE) "./sage -docbuild --no-pdf-links all html $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a dochtml.log"

doc-html-jsmath: build # (already) indirectly depends on $(PIPE)
	$(PIPE) "./sage -docbuild --no-pdf-links all html -j $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a dochtml.log"

doc-pdf: build # (already) indirectly depends on $(PIPE)
	$(PIPE) "./sage -docbuild all pdf $(SAGE_DOCBUILD_OPTS) 2>&1" "tee -a docpdf.log"

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
	rm -rf data
	rm -rf dist
	rm -rf devel
	rm -rf doc
	rm -rf examples
	rm -rf sage-python
	rm -rf spkg/build
	rm -rf spkg/archive
	rm -rf ipython
	rm -rf matplotlibrc
	rm -rf tmp
	rm -f .BUILDSTART

micro_release:
	. local/bin/sage-env && local/bin/sage-micro_release

text-expand:
	./spkg/base/text-expand

text-collapse:
	./spkg/base/text-collapse

TESTPRELIMS = . local/bin/sage-env && sage-starts &&
TESTDIRS = devel/sage/doc/common devel/sage/doc/de devel/sage/doc/en devel/sage/doc/fr devel/sage/doc/ru devel/sage/sage

test: all # i.e. build and (HTML) doc
	@# $(TESTPRELIMS) (also) puts sage-maketest into the path
	$(TESTPRELIMS) sage-maketest

check: test

testall: all # i.e. build and (HTML) doc, and also indirectly $(PIPE)
	./sage -b
	$(PIPE) "$(TESTPRELIMS) ./sage -t -sagenb -optional $(TESTDIRS) 2>&1" "tee -a testall.log"

testlong: all # i.e. build and (HTML) doc, and also indirectly $(PIPE)
	./sage -b
	$(PIPE) "$(TESTPRELIMS) ./sage -t -sagenb -long $(TESTDIRS) 2>&1" "tee -a testlong.log"

ptest: all # i.e. build and (HTML) doc, and also indirectly $(PIPE)
	$(PIPE) "$(TESTPRELIMS) ./sage -tp $(NUM_THREADS) -sagenb $(TESTDIRS) 2>&1" "tee -a ptest.log"

ptestall: all # i.e. build and (HTML) doc, and also indirectly $(PIPE)
	$(PIPE) "$(TESTPRELIMS) ./sage -tp $(NUM_THREADS) -sagenb -optional $(TESTDIRS) 2>&1" "tee -a ptestall.log"

ptestlong: all # i.e. build and (HTML) doc, and also indirectly $(PIPE)
	$(PIPE) "$(TESTPRELIMS) ./sage -tp $(NUM_THREADS) -sagenb -long $(TESTDIRS) 2>&1" "tee -a ptestlong.log"

testoptional: testall # just an alias

ptestoptional: ptestall # just an alias


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


.PHONY: all build build-serial start \
	doc doc-html doc-html-jsmath doc-pdf \
	doc-clean clean	distclean \
	test check testoptional testlong ptest ptestall ptestlong \
	install testall ptestoptional
